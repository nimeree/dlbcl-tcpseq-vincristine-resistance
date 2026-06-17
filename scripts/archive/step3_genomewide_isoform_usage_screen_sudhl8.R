# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)
library(ggplot2)
library(scales)

rdata <- input_path("SUDHL.RData")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
base_out <- analysis_path("LongRead_TRA2A")
out_dir <- file.path(base_out, "Step3_GenomeWide_isoform_usage_screen")
target_file <- file.path(base_out, "Step2_TRA2A_targets_isoform_usage", "TRA2A_eCLIP_target_genes_from_peak_gene_overlap.csv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

message("Loading long-read matrices...")
e <- new.env()
load(rdata, envir = e)
annotation <- as.data.table(e$Annotation)[Condition == "SUDHL8"]
samples <- annotation$Sample

message("Parsing transcript annotation...")
gtf_dt <- fread(
  gtf,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
)
tx_dt <- gtf_dt[feature == "transcript"]
attrs <- tx_dt$attributes
tx_map <- data.table(
  gene_id = grab_attr(attrs, "gene_id"),
  gene_name = grab_attr(attrs, "gene_name"),
  gene_biotype = grab_attr(attrs, "gene_biotype"),
  transcript_id = grab_attr(attrs, "transcript_id"),
  transcript_version = grab_attr(attrs, "transcript_version"),
  transcript_name = grab_attr(attrs, "transcript_name"),
  transcript_biotype = grab_attr(attrs, "transcript_biotype")
)
tx_map[, gene_id := sub("[.][0-9]+$", "", gene_id)]
tx_map[, row_id := paste0(transcript_id, ".", transcript_version)]
tx_map <- tx_map[row_id %in% rownames(e$count_matrix)]

message("Building SUDHL8 transcript proportions...")
counts <- e$count_matrix[tx_map$row_id, samples, drop = FALSE]
counts[is.na(counts)] <- 0
count_dt <- as.data.table(as.table(counts))
setnames(count_dt, c("row_id", "Sample", "count"))
count_dt[, count := as.numeric(count)]
count_dt <- merge(count_dt, tx_map, by = "row_id", all.x = TRUE)
count_dt <- merge(count_dt, annotation, by = "Sample", all.x = TRUE)
count_dt <- count_dt[!is.na(gene_id)]

gene_sample_totals <- count_dt[, .(gene_total_count = sum(count)), by = .(gene_id, gene_name, Sample, Type, Replicate)]
count_dt <- merge(count_dt, gene_sample_totals, by = c("gene_id", "gene_name", "Sample", "Type", "Replicate"))
count_dt[, proportion := fifelse(gene_total_count > 0, count / gene_total_count, NA_real_)]

message("Selecting eligible multi-isoform genes...")
gene_summary <- gene_sample_totals[
  ,
  .(
    total_count = sum(gene_total_count),
    mean_original_count = mean(gene_total_count[Type == "Original"]),
    mean_resistant_count = mean(gene_total_count[Type == "Resistant"]),
    n_samples_detected = sum(gene_total_count > 0)
  ),
  by = .(gene_id, gene_name)
]
gene_tx_n <- unique(count_dt[count > 0, .(gene_id, row_id)])[, .(n_detected_transcripts = .N), by = gene_id]
gene_summary <- merge(gene_summary, gene_tx_n, by = "gene_id", all.x = TRUE)
gene_summary[is.na(n_detected_transcripts), n_detected_transcripts := 0L]

eligible_genes <- gene_summary[
  n_detected_transcripts >= 2 &
    total_count >= 30 &
    n_samples_detected >= 3,
  gene_id
]

message("Testing genome-wide isoform usage for ", length(eligible_genes), " genes...")
test_gene <- function(gid) {
  dt <- count_dt[gene_id == gid]
  agg <- dcast(dt, Type ~ row_id, value.var = "count", fun.aggregate = sum)
  mat <- as.matrix(agg[, -"Type"])
  rownames(mat) <- agg$Type
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (ncol(mat) < 2 || nrow(mat) < 2) return(NULL)

  p_value <- tryCatch({
    suppressWarnings(chisq.test(mat)$p.value)
  }, error = function(e) NA_real_)

  prop_summary <- dt[
    ,
    .(mean_prop = mean(proportion, na.rm = TRUE), mean_count = mean(count)),
    by = .(Type, row_id, transcript_name, transcript_biotype)
  ]
  wide <- dcast(prop_summary, row_id + transcript_name + transcript_biotype ~ Type, value.var = "mean_prop")
  for (nm in c("Original", "Resistant")) {
    if (!nm %in% names(wide)) wide[, (nm) := 0]
    wide[is.na(get(nm)), (nm) := 0]
  }
  wide[, resistant_minus_original := Resistant - Original]

  top_original <- wide[which.max(Original), transcript_name]
  top_resistant <- wide[which.max(Resistant), transcript_name]
  top_shift <- wide[which.max(abs(resistant_minus_original))]

  totals <- gene_sample_totals[gene_id == gid, .(mean_total = mean(gene_total_count)), by = Type]
  original_total <- totals[Type == "Original", mean_total]
  resistant_total <- totals[Type == "Resistant", mean_total]

  data.table(
    gene_id = gid,
    gene_name = dt$gene_name[1],
    gene_biotype = dt$gene_biotype[1],
    n_detected_transcripts = uniqueN(dt$row_id[dt$count > 0]),
    total_count = sum(dt$count),
    original_mean_total_count = ifelse(length(original_total), original_total, NA_real_),
    resistant_mean_total_count = ifelse(length(resistant_total), resistant_total, NA_real_),
    total_count_log2FC_resistant_vs_original = log2((ifelse(length(resistant_total), resistant_total, 0) + 1) / (ifelse(length(original_total), original_total, 0) + 1)),
    top_original_isoform = top_original,
    top_resistant_isoform = top_resistant,
    dominant_isoform_switch = top_original != top_resistant,
    max_abs_proportion_shift = max(abs(wide$resistant_minus_original), na.rm = TRUE),
    max_shift_isoform = top_shift$transcript_name,
    max_shift_transcript_id = top_shift$row_id,
    max_shift_biotype = top_shift$transcript_biotype,
    max_shift_original_prop = top_shift$Original,
    max_shift_resistant_prop = top_shift$Resistant,
    max_shift_resistant_minus_original = top_shift$resistant_minus_original,
    chi_square_pvalue = p_value
  )
}

genome_tests <- rbindlist(lapply(eligible_genes, test_gene), fill = TRUE)
genome_tests[, padj := p.adjust(chi_square_pvalue, method = "BH")]

target_genes <- fread(target_file)
target_genes <- unique(target_genes[, .(gene_id, TRA2A_eCLIP_target = TRUE, source_cell_lines, n_source_cell_lines)])
genome_tests <- merge(genome_tests, target_genes, by = "gene_id", all.x = TRUE)
genome_tests[is.na(TRA2A_eCLIP_target), TRA2A_eCLIP_target := FALSE]
genome_tests[is.na(source_cell_lines), source_cell_lines := ""]
genome_tests[is.na(n_source_cell_lines), n_source_cell_lines := 0L]
genome_tests[, switched_stringent := padj < 0.05 & max_abs_proportion_shift >= 0.20 & total_count >= 50]
genome_tests[, switched_relaxed := padj < 0.05 & max_abs_proportion_shift >= 0.10 & total_count >= 50]
setorder(genome_tests, padj, -max_abs_proportion_shift)

target_enrichment <- function(flag_col) {
  tab <- table(
    switched = genome_tests[[flag_col]],
    tra2a_target = genome_tests$TRA2A_eCLIP_target
  )
  ft <- fisher.test(tab)
  data.table(
    definition = flag_col,
    switched_genes = sum(genome_tests[[flag_col]]),
    target_in_switched = sum(genome_tests[[flag_col]] & genome_tests$TRA2A_eCLIP_target),
    target_not_switched = sum(!genome_tests[[flag_col]] & genome_tests$TRA2A_eCLIP_target),
    non_target_switched = sum(genome_tests[[flag_col]] & !genome_tests$TRA2A_eCLIP_target),
    non_target_not_switched = sum(!genome_tests[[flag_col]] & !genome_tests$TRA2A_eCLIP_target),
    odds_ratio = unname(ft$estimate),
    pvalue = ft$p.value
  )
}
enrichment <- rbindlist(list(target_enrichment("switched_stringent"), target_enrichment("switched_relaxed")))

message("Writing genome-wide outputs...")
write.csv(genome_tests, file.path(out_dir, "SUDHL8_genomewide_isoform_usage_screen.csv"), row.names = FALSE)
write.csv(genome_tests[switched_stringent == TRUE], file.path(out_dir, "SUDHL8_stringent_isoform_switch_genes.csv"), row.names = FALSE)
write.csv(genome_tests[switched_relaxed == TRUE], file.path(out_dir, "SUDHL8_relaxed_isoform_switch_genes.csv"), row.names = FALSE)
write.csv(enrichment, file.path(out_dir, "TRA2A_target_enrichment_among_isoform_switches.csv"), row.names = FALSE)

volcano_dt <- copy(genome_tests)
volcano_dt[, neg_log10_padj := -log10(pmax(padj, .Machine$double.xmin))]
volcano_dt[, class := fifelse(switched_stringent & TRA2A_eCLIP_target, "Switched + TRA2A target",
                              fifelse(switched_stringent, "Switched",
                                      fifelse(TRA2A_eCLIP_target, "TRA2A target", "Other")))]
volcano_dt[, class := factor(class, levels = c("Other", "TRA2A target", "Switched", "Switched + TRA2A target"))]

p1 <- ggplot(volcano_dt, aes(x = max_abs_proportion_shift, y = neg_log10_padj, color = class)) +
  geom_point(alpha = 0.65, size = 1.4) +
  scale_color_manual(values = c("Other" = "grey72", "TRA2A target" = "#4C78A8", "Switched" = "#E45756", "Switched + TRA2A target" = "#7B3294")) +
  geom_vline(xintercept = 0.20, linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35) +
  labs(
    title = "Genome-wide SUDHL8 isoform-usage screen",
    subtitle = "Original vs Resistant; points are expressed multi-isoform genes",
    x = "Maximum isoform proportion shift",
    y = "-log10(BH-adjusted P)",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(out_dir, "SUDHL8_genomewide_isoform_usage_screen_volcano.png"), p1, width = 9, height = 6, dpi = 300)

top_genes <- head(genome_tests[switched_stringent == TRUE][order(padj, -max_abs_proportion_shift), gene_id], 20)
if (length(top_genes) > 0) {
  top_dt <- count_dt[gene_id %in% top_genes]
  top_dt[, Type := factor(Type, levels = c("Original", "Resistant"))]
  top_sum <- top_dt[
    ,
    .(prop = sum(count) / unique(gene_total_count)),
    by = .(Sample, Type, gene_id, gene_name, transcript_name)
  ]
  gene_levels <- genome_tests[gene_id %in% top_genes][order(padj, -max_abs_proportion_shift), gene_name]
  top_sum[, gene_name := factor(gene_name, levels = gene_levels)]

  p2 <- ggplot(top_sum, aes(x = Type, y = prop, fill = transcript_name)) +
    geom_col(width = 0.78, color = "white", linewidth = 0.1, position = "fill") +
    facet_wrap(~ gene_name, scales = "free_y") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Top genome-wide SUDHL8 isoform switches",
      subtitle = "Stacked transcript proportions, Original vs Resistant",
      x = NULL,
      y = "Isoform proportion",
      fill = "Transcript"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"), panel.grid.major.x = element_blank())
  ggsave(file.path(out_dir, "top_genomewide_isoform_switches_SUDHL8.png"), p2, width = 11, height = 9, dpi = 300)
}

readme <- c(
  "Step 3 first-pass genome-wide isoform-usage screen.",
  "Input: <THESIS_INPUT_DIR>/SUDHL.RData transcript count matrix.",
  "Samples: SUDHL8 Original vs SUDHL8 Resistant only.",
  "Eligible genes: expressed genes with at least 2 detected transcripts, total count >= 30, detected in at least 3 samples.",
  "Statistical screen: chi-square test on aggregated transcript counts by condition, BH adjusted genome-wide.",
  "Effect size: maximum absolute mean within-gene isoform proportion shift across transcripts.",
  "Stringent switched definition: padj < 0.05, max_abs_proportion_shift >= 0.20, total_count >= 50.",
  "Relaxed switched definition: padj < 0.05, max_abs_proportion_shift >= 0.10, total_count >= 50.",
  "TRA2A target enrichment uses Step 2 ENCODE TRA2A eCLIP peak-overlap target genes.",
  "This is a first-pass screen, not a replacement for DRIMSeq/IsoformSwitchAnalyzeR."
)
writeLines(readme, file.path(out_dir, "README_step3_first_pass_method.txt"))

cat("\nGenome-wide eligible genes:", nrow(genome_tests), "\n")
cat("Stringent switched genes:", genome_tests[switched_stringent == TRUE, .N], "\n")
cat("Relaxed switched genes:", genome_tests[switched_relaxed == TRUE, .N], "\n")
cat("TRA2A targets among eligible genes:", genome_tests[TRA2A_eCLIP_target == TRUE, .N], "\n")
cat("\nTRA2A target enrichment:\n")
print(enrichment)
cat("\nTop stringent switched genes:\n")
print(head(genome_tests[switched_stringent == TRUE], 25))
cat("\nOutput directory:\n")
cat(out_dir, "\n")
