# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

base_dir <- analysis_path()
index_dir <- file.path(base_dir, "Translation_indexes_fixed")
limma_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05")
out_dir <- file.path(limma_dir, "Scanning_score_start_codon_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

psite_file <- file.path(index_dir, "transcript_psite_matrix_long_ALL_samples.csv")
metric_file <- file.path(index_dir, "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")
gtf_file <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
limma_scan_file <- file.path(limma_dir, "Results", "scanning_score", "Interaction_limma_all_genes.csv")

p_cut <- 0.05
lfc_cut <- 0.7
window_min <- -150L
window_max <- 150L
summary_window <- c(-100L, 30L)

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

condition_from_sample <- function(sample) {
  fifelse(grepl("SU8R-DMSO|SU8-R-DMSO", sample), "Resistant DMSO",
  fifelse(grepl("SU8R-Vin|SU8-R-Vin", sample), "Resistant VCR",
  fifelse(grepl("SU8-DMSO", sample), "Sensitive DMSO",
  fifelse(grepl("SU8-Vin", sample), "Sensitive VCR", NA_character_))))
}

message("Reading limma scanning interaction results...")
scan <- fread(limma_scan_file)
scan <- scan[!is.na(gene_name) & gene_name != ""]
scan[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
scan[, scanning_group := fifelse(
  P.Value < p_cut & logFC >= lfc_cut, "Scanning Up",
  fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Scanning Down", "Scanning NS")
)]
scan <- scan[, .(
  gene_id_clean,
  gene_name,
  scanning_logFC = logFC,
  scanning_p = P.Value,
  scanning_fdr = adj.P.Val,
  scanning_group
)]

message("Finding top covered transcript per gene...")
tx_map <- unique(fread(metric_file, select = c("transcript", "gene_id_clean", "gene_name")))
tx_map <- tx_map[!is.na(gene_name) & gene_name != ""]
psite_top <- fread(psite_file, select = c("transcript", "codon_pos", "psite_count"))
psite_top <- merge(psite_top, tx_map, by = "transcript")
top_tx <- psite_top[, .(
  total_psites = sum(psite_count, na.rm = TRUE),
  n_positions = uniqueN(codon_pos)
), by = .(gene_id_clean, gene_name, transcript)][
  order(gene_id_clean, -total_psites, -n_positions)
][, .SD[1], by = .(gene_id_clean, gene_name)]
rm(psite_top)

message("Parsing GTF transcript CDS start coordinates...")
gtf <- fread(
  gtf_file,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
)
gtf <- gtf[feature %in% c("exon", "CDS")]
gtf[, transcript_id := grab_attr(attributes, "transcript_id")]
gtf[, transcript_version := grab_attr(attributes, "transcript_version")]
gtf[, transcript := paste0(transcript_id, ".", transcript_version)]

transcript_cds_start <- function(dt) {
  strand <- dt$strand[1]
  exons <- unique(dt[feature == "exon", .(start, end)])
  cds <- dt[feature == "CDS"]
  if (!nrow(exons) || !nrow(cds)) return(NA_integer_)

  if (strand == "+") {
    exons <- exons[order(start, end)]
    cds_site <- min(cds$start)
    hit <- which(exons$start <= cds_site & exons$end >= cds_site)[1]
    if (is.na(hit)) return(NA_integer_)
    before <- if (hit > 1) sum(exons$end[seq_len(hit - 1)] - exons$start[seq_len(hit - 1)] + 1L) else 0L
    before + (cds_site - exons$start[hit] + 1L)
  } else if (strand == "-") {
    exons <- exons[order(-end, -start)]
    cds_site <- max(cds$end)
    hit <- which(exons$start <= cds_site & exons$end >= cds_site)[1]
    if (is.na(hit)) return(NA_integer_)
    before <- if (hit > 1) sum(exons$end[seq_len(hit - 1)] - exons$start[seq_len(hit - 1)] + 1L) else 0L
    before + (exons$end[hit] - cds_site + 1L)
  } else {
    NA_integer_
  }
}

cds_start <- gtf[
  ,
  .(
    cds_start_tx = transcript_cds_start(.SD),
    strand = strand[1]
  ),
  by = transcript
]
cds_start <- cds_start[!is.na(cds_start_tx)]

analysis_tx <- merge(scan, top_tx, by = c("gene_id_clean", "gene_name"))
analysis_tx <- merge(analysis_tx, cds_start, by = "transcript")
analysis_tx <- analysis_tx[total_psites > 0 & n_positions >= 20]
analysis_tx[, scanning_group := factor(scanning_group, levels = c("Scanning Down", "Scanning NS", "Scanning Up"))]

message("Reading SSU P-site positions around CDS start...")
psite <- fread(psite_file, select = c("sample", "fraction", "transcript", "codon_pos", "psite_cpm", "psite_count"))
psite <- psite[fraction == "SSU" & transcript %in% analysis_tx$transcript]
psite <- merge(psite, analysis_tx, by = "transcript")
psite[, rel_pos := as.integer(codon_pos - cds_start_tx)]
psite <- psite[rel_pos >= window_min & rel_pos <= window_max]
psite[, condition := condition_from_sample(sample)]
psite <- psite[!is.na(condition)]
psite[, condition := factor(condition, levels = c("Sensitive DMSO", "Sensitive VCR", "Resistant DMSO", "Resistant VCR"))]

message("Completing zero-filled start-codon window...")
samples <- unique(psite[, .(sample, condition)])
grid <- CJ(
  transcript = unique(analysis_tx$transcript),
  sample = samples$sample,
  rel_pos = seq(window_min, window_max)
)
grid <- merge(grid, samples, by = "sample")
grid <- merge(grid, analysis_tx, by = "transcript")
psite_small <- psite[, .(transcript, sample, rel_pos, psite_cpm, psite_count)]
psite_full <- merge(grid, psite_small, by = c("transcript", "sample", "rel_pos"), all.x = TRUE)
psite_full[is.na(psite_cpm), psite_cpm := 0]
psite_full[is.na(psite_count), psite_count := 0]

cond_gene_window <- psite_full[
  rel_pos >= summary_window[1] & rel_pos <= summary_window[2],
  .(
    mean_ssu_start_cpm = mean(psite_cpm, na.rm = TRUE),
    sum_ssu_start_cpm = sum(psite_cpm, na.rm = TRUE),
    sum_ssu_start_count = sum(psite_count, na.rm = TRUE)
  ),
  by = .(gene_id_clean, gene_name, transcript, condition, scanning_group, scanning_logFC, scanning_p, scanning_fdr)
]

wide <- dcast(
  cond_gene_window,
  gene_id_clean + gene_name + transcript + scanning_group + scanning_logFC + scanning_p + scanning_fdr ~ condition,
  value.var = "mean_ssu_start_cpm"
)
for (nm in c("Sensitive DMSO", "Sensitive VCR", "Resistant DMSO", "Resistant VCR")) {
  if (!nm %in% names(wide)) wide[, (nm) := 0]
  wide[is.na(get(nm)), (nm) := 0]
}
wide[, resistant_vcr_response := `Resistant VCR` - `Resistant DMSO`]
wide[, sensitive_vcr_response := `Sensitive VCR` - `Sensitive DMSO`]
wide[, ssu_start_interaction_delta := resistant_vcr_response - sensitive_vcr_response]

cor_all <- suppressWarnings(cor.test(wide$scanning_logFC, wide$ssu_start_interaction_delta, method = "spearman", exact = FALSE))
group_tests <- rbindlist(list(
  data.table(
    comparison = "Scanning Up vs Scanning NS",
    p_value = tryCatch(wilcox.test(
      wide[scanning_group == "Scanning Up", ssu_start_interaction_delta],
      wide[scanning_group == "Scanning NS", ssu_start_interaction_delta],
      exact = FALSE
    )$p.value, error = function(e) NA_real_)
  ),
  data.table(
    comparison = "Scanning Up vs Scanning Down",
    p_value = tryCatch(wilcox.test(
      wide[scanning_group == "Scanning Up", ssu_start_interaction_delta],
      wide[scanning_group == "Scanning Down", ssu_start_interaction_delta],
      exact = FALSE
    )$p.value, error = function(e) NA_real_)
  ),
  data.table(
    comparison = "Scanning Down vs Scanning NS",
    p_value = tryCatch(wilcox.test(
      wide[scanning_group == "Scanning Down", ssu_start_interaction_delta],
      wide[scanning_group == "Scanning NS", ssu_start_interaction_delta],
      exact = FALSE
    )$p.value, error = function(e) NA_real_)
  )
))
group_tests[, BH_p_value := p.adjust(p_value, method = "BH")]

group_summary <- wide[, .(
  n_genes = .N,
  median_ssu_start_interaction_delta = median(ssu_start_interaction_delta, na.rm = TRUE),
  mean_ssu_start_interaction_delta = mean(ssu_start_interaction_delta, na.rm = TRUE),
  median_scanning_logFC = median(scanning_logFC, na.rm = TRUE)
), by = scanning_group]

metagene <- psite_full[
  ,
  .(mean_ssu_cpm = mean(psite_cpm, na.rm = TRUE)),
  by = .(scanning_group, condition, gene_name, rel_pos)
][
  ,
  .(
    mean_ssu_cpm = mean(mean_ssu_cpm, na.rm = TRUE),
    se_ssu_cpm = sd(mean_ssu_cpm, na.rm = TRUE) / sqrt(.N),
    n_genes = .N
  ),
  by = .(scanning_group, condition, rel_pos)
]

meta_wide <- dcast(
  metagene,
  scanning_group + rel_pos ~ condition,
  value.var = "mean_ssu_cpm"
)
for (nm in c("Sensitive DMSO", "Sensitive VCR", "Resistant DMSO", "Resistant VCR")) {
  if (!nm %in% names(meta_wide)) meta_wide[, (nm) := 0]
  meta_wide[is.na(get(nm)), (nm) := 0]
}
meta_wide[, ssu_interaction_delta := (`Resistant VCR` - `Resistant DMSO`) - (`Sensitive VCR` - `Sensitive DMSO`)]

fwrite(analysis_tx, file.path(out_dir, "scanning_validation_top_transcripts_with_cds_start.csv"))
fwrite(wide, file.path(out_dir, "gene_level_ssu_start_window_interaction_delta.csv"))
fwrite(group_summary, file.path(out_dir, "ssu_start_window_group_summary.csv"))
fwrite(group_tests, file.path(out_dir, "ssu_start_window_group_wilcoxon_tests.csv"))
fwrite(data.table(
  test = "Spearman scanning logFC vs SSU start-window interaction delta",
  n_genes = nrow(wide),
  rho = unname(cor_all$estimate),
  p_value = cor_all$p.value,
  window = paste(summary_window, collapse = " to ")
), file.path(out_dir, "scanning_logFC_vs_ssu_start_window_correlation.csv"))
fwrite(meta_wide, file.path(out_dir, "metagene_ssu_start_region_interaction_delta_by_scanning_group.csv"))

group_cols <- c("Scanning Down" = "#2C7BB6", "Scanning NS" = "#9CA3AF", "Scanning Up" = "#D7191C")

p_meta <- ggplot(meta_wide, aes(rel_pos, ssu_interaction_delta, color = scanning_group)) +
  geom_hline(yintercept = 0, color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey35") +
  geom_smooth(se = FALSE, method = "loess", span = 0.18, linewidth = 0.9) +
  scale_color_manual(values = group_cols, name = NULL) +
  labs(
    title = "SSU P-site signal around the start codon by scanning-score class",
    subtitle = "Interaction-style SSU delta: (Resistant VCR - Resistant DMSO) - (Sensitive VCR - Sensitive DMSO)",
    x = "Position relative to annotated CDS start (nt)",
    y = "SSU P-site CPM interaction delta"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "A_scanning_group_start_codon_SSU_interaction_metagene.png"), p_meta, width = 8.2, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "A_scanning_group_start_codon_SSU_interaction_metagene.pdf"), p_meta, width = 8.2, height = 5.2, bg = "white")

p_box <- ggplot(wide, aes(scanning_group, ssu_start_interaction_delta, fill = scanning_group)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.78, width = 0.62) +
  geom_jitter(width = 0.13, alpha = 0.28, size = 0.8) +
  scale_fill_manual(values = group_cols, guide = "none") +
  labs(
    title = "Start-region SSU density does not automatically follow the scanning score",
    subtitle = "Window: -100 to +30 nt around annotated CDS start",
    x = NULL,
    y = "SSU start-region interaction delta"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "B_ssu_start_window_delta_by_scanning_group.png"), p_box, width = 6.6, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "B_ssu_start_window_delta_by_scanning_group.pdf"), p_box, width = 6.6, height = 5.2, bg = "white")

p_scatter <- ggplot(wide, aes(ssu_start_interaction_delta, scanning_logFC, color = scanning_group)) +
  geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", color = "grey45") +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(alpha = 0.55, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.65) +
  scale_color_manual(values = group_cols, name = NULL) +
  labs(
    title = "Scanning score versus positional SSU start-region signal",
    subtitle = sprintf("Spearman rho = %.3f, p = %.3g", unname(cor_all$estimate), cor_all$p.value),
    x = "SSU start-region interaction delta (-100 to +30 nt)",
    y = "limma scanning-score interaction logFC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "C_scanning_logFC_vs_SSU_start_window_delta.png"), p_scatter, width = 7.2, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "C_scanning_logFC_vs_SSU_start_window_delta.pdf"), p_scatter, width = 7.2, height = 5.5, bg = "white")

cat("\nOutput directory:\n", out_dir, "\n", sep = "")
cat("\nGenes with usable top transcript and CDS start:", nrow(wide), "\n")
cat("\nScanning group summary:\n")
print(group_summary)
cat("\nGroup Wilcoxon tests:\n")
print(group_tests)
cat("\nCorrelation scanning logFC vs positional SSU start-window interaction delta:\n")
print(data.table(
  n_genes = nrow(wide),
  rho = unname(cor_all$estimate),
  p_value = cor_all$p.value,
  window = paste(summary_window, collapse = " to ")
))
