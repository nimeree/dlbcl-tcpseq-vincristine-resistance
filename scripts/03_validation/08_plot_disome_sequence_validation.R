# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
  library(ggplot2)
})

FASTA <- input_path("Homo_sapiens.GRCh38.cds.all.fa.gz")
GENE_METRICS <- analysis_path("Translation_indexes_fixed", "Gene_Level_Clean", "gene_level_clean_collision_complete_8_samples.csv")
OUT_DIR <- analysis_path("Translation_indexes_fixed", "Validation_Plots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, width, height) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = width, height = height)
}

count_pattern <- function(x, pat) {
  n <- nchar(x)
  k <- nchar(pat)
  if (is.na(x) || n < k) return(0L)
  sum(substring(x, seq_len(n - k + 1L), seq_len(n - k + 1L) + k - 1L) == pat)
}

count_rxk <- function(x) {
  n <- nchar(x)
  if (is.na(x) || n < 3L) return(0L)
  tri <- substring(x, seq_len(n - 2L), seq_len(n - 2L) + 2L)
  sum(substr(tri, 1L, 1L) == "R" & substr(tri, 3L, 3L) == "K")
}

codon_vector <- function(dna) {
  s <- toupper(as.character(dna))
  n_codons <- floor(nchar(s) / 3)
  if (n_codons < 1L) return(character())
  substring(s, seq(1L, by = 3L, length.out = n_codons), seq(3L, by = 3L, length.out = n_codons))
}

format_p <- function(p) {
  ifelse(p < 2.2e-16, "p < 2.2e-16", paste0("p = ", signif(p, 3)))
}

message("Reading CDS FASTA")
cds <- readDNAStringSet(FASTA)
headers <- names(cds)
gene_id_clean <- sub("\\.\\d+$", "", sub(".* gene:([^ ]+).*", "\\1", headers))
gene_symbol <- sub(".* gene_symbol:([^ ]+).*", "\\1", headers)
gene_symbol[gene_symbol == headers] <- NA_character_

message("Computing codon-frequency rare-codon proxy")
codons_list <- lapply(cds, codon_vector)
all_codons <- unlist(codons_list, use.names = FALSE)
all_codons <- all_codons[grepl("^[ACGT]{3}$", all_codons) & !(all_codons %in% c("TAA", "TAG", "TGA"))]
codon_freq <- sort(table(all_codons) / length(all_codons))
rare_codons <- names(codon_freq)[seq_len(ceiling(length(codon_freq) * 0.20))]

message("Computing transcript-level motif densities")
aa <- suppressWarnings(translate(cds, if.fuzzy.codon = "X", no.init.codon = TRUE))
aa_chr <- sub("\\*$", "", as.character(aa))

tx_seq <- data.table(
  gene_id_clean = gene_id_clean,
  gene_symbol = gene_symbol,
  aa = aa_chr,
  codons = codons_list
)
tx_seq[, aa_len := nchar(aa)]
tx_seq[, cds_kb := pmax(aa_len * 3 / 1000, 1e-9)]
tx_seq[, `:=`(
  PPP_count = vapply(aa, count_pattern, integer(1), pat = "PPP"),
  PPG_count = vapply(aa, count_pattern, integer(1), pat = "PPG"),
  PPD_count = vapply(aa, count_pattern, integer(1), pat = "PPD"),
  RXK_count = vapply(aa, count_rxk, integer(1)),
  rare_codon_fraction_proxy = vapply(codons, function(x) {
    x <- x[grepl("^[ACGT]{3}$", x) & !(x %in% c("TAA", "TAG", "TGA"))]
    if (!length(x)) return(NA_real_)
    mean(x %in% rare_codons)
  }, numeric(1))
)]
tx_seq[, polyproline_count := PPP_count + PPG_count + PPD_count]
tx_seq[, `:=`(
  PPP_density_per_kb = PPP_count / cds_kb,
  polyproline_density_per_kb = polyproline_count / cds_kb,
  RXK_density_per_kb = RXK_count / cds_kb
)]

gene_seq <- tx_seq[!is.na(gene_id_clean) & gene_id_clean != "", .(
  gene_symbol = names(sort(table(gene_symbol), decreasing = TRUE))[1],
  aa_len = as.numeric(median(aa_len, na.rm = TRUE)),
  PPP_density_per_kb = as.numeric(median(PPP_density_per_kb, na.rm = TRUE)),
  polyproline_density_per_kb = as.numeric(median(polyproline_density_per_kb, na.rm = TRUE)),
  RXK_density_per_kb = as.numeric(median(RXK_density_per_kb, na.rm = TRUE)),
  rare_codon_fraction_proxy = as.numeric(median(rare_codon_fraction_proxy, na.rm = TRUE))
), by = gene_id_clean]

message("Reading clean collision matrix")
metrics <- fread(GENE_METRICS)
gene_collision_all <- metrics[, .(
  gene_name = names(sort(table(gene_name), decreasing = TRUE))[1],
  collision_score = median(collision_score, na.rm = TRUE),
  collision_score_DMSO = median(collision_score[treatment == "DMSO"], na.rm = TRUE),
  collision_score_VCR = median(collision_score[treatment == "VCR"], na.rm = TRUE),
  rs_core_cpm = median(rs_core_cpm, na.rm = TRUE)
), by = gene_id_clean]

joined <- merge(gene_collision_all, gene_seq, by = "gene_id_clean", all.x = FALSE)
fwrite(joined, file.path(OUT_DIR, "disome_gene_level_sequence_validation_joined_data.csv"))

predictors <- c("PPP_density_per_kb", "polyproline_density_per_kb", "RXK_density_per_kb", "rare_codon_fraction_proxy")
predictor_labels <- c(
  PPP_density_per_kb = "PPP density per kb CDS",
  polyproline_density_per_kb = "PPG/PPD/PPP density per kb CDS",
  RXK_density_per_kb = "Arg-X-Lys density per kb CDS",
  rare_codon_fraction_proxy = "Rare codon fraction proxy"
)

run_cor <- function(ycol, ylab) {
  rbindlist(lapply(predictors, function(pred) {
    x <- joined[[pred]]
    y <- joined[[ycol]]
    ok <- is.finite(x) & is.finite(y)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    data.table(
      collision_metric = ylab,
      predictor = pred,
      predictor_label = predictor_labels[pred],
      n = sum(ok),
      spearman_rho = unname(ct$estimate),
      spearman_p = ct$p.value
    )
  }))
}

cor_results <- rbindlist(list(
  run_cor("collision_score", "All samples median collision_score"),
  run_cor("collision_score_DMSO", "DMSO median collision_score"),
  run_cor("collision_score_VCR", "VCR median collision_score")
))
fwrite(cor_results, file.path(OUT_DIR, "disome_gene_level_sequence_validation_correlations.csv"))

plot_long <- melt(
  joined,
  id.vars = c("gene_id_clean", "gene_name", "collision_score", "collision_score_DMSO"),
  measure.vars = predictors,
  variable.name = "predictor",
  value.name = "predictor_value"
)
plot_long <- plot_long[is.finite(predictor_value) & is.finite(collision_score)]
plot_long[, predictor_label := factor(predictor_labels[predictor], levels = predictor_labels)]

ann <- cor_results[collision_metric == "All samples median collision_score"]
ann[, predictor_label := factor(predictor_label, levels = predictor_labels)]
ann_pos <- plot_long[, .(
  x = quantile(predictor_value, 0.04, na.rm = TRUE),
  y = quantile(collision_score, 0.96, na.rm = TRUE)
), by = predictor_label]
ann <- merge(ann, ann_pos, by = "predictor_label", all.x = TRUE)
ann[, label := paste0("Spearman rho = ", sprintf("%.3f", spearman_rho), "\n", format_p(spearman_p), "\nn = ", n)]

p <- ggplot(plot_long, aes(x = predictor_value, y = collision_score)) +
  geom_point(alpha = 0.22, size = 0.65, color = "#2C7A7B") +
  geom_smooth(method = "lm", se = TRUE, color = "#B8323B", fill = "#F1C9CE", linewidth = 0.65) +
  geom_label(
    data = ann,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.5,
    fill = "white",
    linewidth = 0.2
  ) +
  facet_wrap(~ predictor_label, scales = "free_x", ncol = 2) +
  labs(
    title = "Gene-Level Sequence Motif Density vs Collision Score",
    subtitle = "Clean collision-complete gene set; gene-level medians across all samples",
    x = "Sequence feature density",
    y = "Collision score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 15)
  )
save_plot(p, "disome_gene_level_motif_density_collision_validation", 10.5, 7.8)

cat("\nGene-level sequence validation for collision_score\n")
cat("==================================================\n")
cat("Clean collision genes joined to CDS features:", nrow(joined), "\n")
cat("Rare codon proxy uses bottom 20% codons by frequency in the Ensembl CDS FASTA:\n")
cat(paste(rare_codons, collapse = ", "), "\n\n")
print(cor_results[, .(
  collision_metric,
  predictor = predictor_label,
  n,
  spearman_rho = round(spearman_rho, 4),
  spearman_p = signif(spearman_p, 3)
)])
cat("\nOutputs written to:", OUT_DIR, "\n")
