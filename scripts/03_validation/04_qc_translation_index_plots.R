# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

INFILE <- analysis_path("Translation_indexes_fixed", "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")
OUT_DIR <- analysis_path("Translation_indexes_fixed", "QC_Metric_Confidence")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, width, height) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = width, height = height)
}

parse_sample_meta <- function(sample) {
  x <- data.table(sample = sample)
  x[, cell_line := fifelse(grepl("^SU8R", sample, ignore.case = TRUE), "Resistant", "Sensitive")]
  x[, treatment := fifelse(grepl("-Vin-", sample, ignore.case = TRUE), "Vin", "DMSO")]
  x[, replicate := fifelse(grepl("Rep1", sample, ignore.case = TRUE), "Rep1",
    fifelse(grepl("Rep2", sample, ignore.case = TRUE), "Rep2", NA_character_))]
  x[, condition := paste(cell_line, treatment, sep = "_")]
  x
}

metric_cols <- c(
  "ribosome_efficiency_score",
  "protein_output_score",
  "rs_core_cpm",
  "rs_rate",
  "collision_index",
  "ssu_scanning_index",
  "initiation_rate_index",
  "elongation_rate_index",
  "total_translation_rate_proxy"
)

message("Loading metrics")
cols <- unique(c(
  "sample", "pair_key", "transcript", "gene_id_clean", "gene_name",
  "n_cds", "n_core", "initiation_rate_index_stable",
  metric_cols
))
dt <- fread(INFILE, select = cols)
for (cc in intersect(metric_cols, names(dt))) dt[, (cc) := as.numeric(get(cc))]
dt[, stable_mask := !is.na(initiation_rate_index_stable)]

meta <- parse_sample_meta(unique(dt$sample))
dt <- merge(dt, meta, by = "sample", all.x = TRUE)
dt[, sample_label := paste(condition, replicate, sep = "_")]
dt[, sample_label := factor(sample_label, levels = meta[order(condition, replicate), paste(condition, replicate, sep = "_")])]

message("Writing stable-mask coverage")
stable_summary <- dt[, .(
  total_rows = .N,
  stable_rows = sum(stable_mask, na.rm = TRUE),
  stable_fraction = mean(stable_mask, na.rm = TRUE),
  genes_total = uniqueN(gene_id_clean),
  genes_stable = uniqueN(gene_id_clean[stable_mask == TRUE])
), by = .(sample, sample_label, cell_line, treatment, replicate, condition)]
stable_summary[, stable_gene_fraction := genes_stable / genes_total]
stable_summary <- stable_summary[order(condition, replicate)]
fwrite(stable_summary, file.path(OUT_DIR, "stable_mask_coverage_by_sample.csv"))

condition_summary <- stable_summary[, .(
  mean_stable_fraction = mean(stable_fraction, na.rm = TRUE),
  mean_stable_gene_fraction = mean(stable_gene_fraction, na.rm = TRUE),
  min_stable_gene_fraction = min(stable_gene_fraction, na.rm = TRUE),
  max_stable_gene_fraction = max(stable_gene_fraction, na.rm = TRUE),
  mean_genes_total = mean(genes_total, na.rm = TRUE),
  mean_genes_stable = mean(genes_stable, na.rm = TRUE)
), by = condition][order(condition)]
fwrite(condition_summary, file.path(OUT_DIR, "stable_mask_coverage_by_condition.csv"))

p_stable <- ggplot(stable_summary, aes(x = sample_label, y = stable_gene_fraction, fill = condition)) +
  geom_col(width = 0.72, alpha = 0.9) +
  geom_text(aes(label = scales::percent(stable_gene_fraction, accuracy = 0.1)), vjust = -0.3, size = 3.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Stable Metric Coverage by Sample",
    subtitle = "Fraction of genes passing n_cds and n_core count thresholds",
    x = NULL,
    y = "Stable gene fraction"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none", panel.grid.minor = element_blank())
save_plot(p_stable, "stable_mask_gene_fraction_by_sample", 9.5, 5)

message("Plotting per-sample metric distributions")
dist_metrics <- c("ribosome_efficiency_score", "protein_output_score", "rs_core_cpm", "rs_rate", "collision_index")
dist_dt <- melt(
  dt,
  id.vars = c("sample_label", "condition"),
  measure.vars = dist_metrics,
  variable.name = "metric",
  value.name = "value"
)
dist_dt <- dist_dt[is.finite(value)]
dist_dt[, plot_value := value]
dist_dt[metric %in% c("rs_core_cpm", "rs_rate", "collision_index"), plot_value := log10(pmax(value, 0) + 1)]
dist_dt[, metric_label := fcase(
  metric == "ribosome_efficiency_score", "Ribosome efficiency score",
  metric == "protein_output_score", "Protein output score",
  metric == "rs_core_cpm", "log10(RS core CPM + 1)",
  metric == "rs_rate", "log10(RS rate + 1)",
  metric == "collision_index", "log10(Collision index + 1)",
  default = metric
)]
dist_dt[, metric_label := factor(metric_label, levels = c(
  "Ribosome efficiency score",
  "Protein output score",
  "log10(RS core CPM + 1)",
  "log10(RS rate + 1)",
  "log10(Collision index + 1)"
))]

p_violin <- ggplot(dist_dt, aes(x = sample_label, y = plot_value, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.76, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9, linewidth = 0.2) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 1) +
  labs(
    title = "Per-Sample Metric Distributions",
    subtitle = "Transcript/gene rows from fixed translation-index matrix",
    x = NULL,
    y = "Metric value"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none", panel.grid.minor = element_blank())
save_plot(p_violin, "metric_distributions_violin_by_sample", 10, 13)

p_density <- ggplot(dist_dt, aes(x = plot_value, color = sample_label, group = sample_label)) +
  geom_density(linewidth = 0.45, alpha = 0.85) +
  facet_wrap(~ metric_label, scales = "free", ncol = 2) +
  labs(
    title = "Per-Sample Metric Density Profiles",
    x = "Metric value",
    y = "Density",
    color = "Sample"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7), panel.grid.minor = element_blank())
save_plot(p_density, "metric_distributions_density_by_sample", 12, 9)

make_gene_matrix <- function(metric) {
  x <- dt[is.finite(get(metric)), .(
    value = median(get(metric), na.rm = TRUE)
  ), by = .(gene_id_clean, sample_label)]
  wide <- dcast(x, gene_id_clean ~ sample_label, value.var = "value")
  sample_cols <- setdiff(names(wide), "gene_id_clean")
  mat <- as.matrix(wide[, ..sample_cols])
  rownames(mat) <- wide$gene_id_clean
  storage.mode(mat) <- "numeric"
  list(wide = wide, mat = mat)
}

run_pca <- function(metric, transform = identity, min_complete = 8) {
  gm <- make_gene_matrix(metric)
  mat <- transform(gm$mat)
  keep <- rowSums(is.finite(mat)) >= min_complete
  mat <- mat[keep, , drop = FALSE]
  complete <- complete.cases(mat)
  mat <- mat[complete, , drop = FALSE]
  mat <- mat[apply(mat, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  pct <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  scores <- as.data.table(pca$x[, 1:3, drop = FALSE], keep.rownames = "sample_label")
  scores <- merge(scores, unique(dt[, .(sample_label, sample, cell_line, treatment, replicate, condition)]), by = "sample_label", all.x = TRUE)
  list(scores = scores, pct = pct, n_genes = nrow(mat))
}

message("Running PCA")
pca_specs <- list(
  ribosome_efficiency_score = list(metric = "ribosome_efficiency_score", label = "Ribosome efficiency score", transform = identity),
  rs_core_cpm = list(metric = "rs_core_cpm", label = "log10(RS core CPM + 1)", transform = function(x) log10(x + 1)),
  protein_output_score = list(metric = "protein_output_score", label = "Protein output score", transform = identity)
)
pca_summary <- list()
for (nm in names(pca_specs)) {
  spec <- pca_specs[[nm]]
  res <- run_pca(spec$metric, spec$transform)
  scores <- res$scores
  scores[, metric := spec$label]
  scores[, n_genes_used := res$n_genes]
  fwrite(scores, file.path(OUT_DIR, paste0("PCA_scores_", nm, ".csv")))
  pca_summary[[nm]] <- data.table(
    metric = spec$label,
    n_genes_used = res$n_genes,
    PC1_percent = res$pct[1],
    PC2_percent = res$pct[2],
    PC3_percent = res$pct[3]
  )
  p_pca <- ggplot(scores, aes(x = PC1, y = PC2, color = condition, shape = replicate, label = replicate)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(vjust = -1.05, size = 3.2, show.legend = FALSE) +
    labs(
      title = paste0("PCA of ", spec$label),
      subtitle = sprintf("Genes used: %s; PC1 %.1f%%, PC2 %.1f%%", format(res$n_genes, big.mark = ","), res$pct[1], res$pct[2]),
      x = sprintf("PC1 (%.1f%%)", res$pct[1]),
      y = sprintf("PC2 (%.1f%%)", res$pct[2])
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  save_plot(p_pca, paste0("PCA_", nm), 7.2, 5.8)
}
fwrite(rbindlist(pca_summary), file.path(OUT_DIR, "PCA_variance_summary.csv"))

message("Computing sample correlation matrices")
cor_metrics <- c("ribosome_efficiency_score", "protein_output_score", "rs_core_cpm", "rs_rate")
cor_summary <- list()
for (metric in cor_metrics) {
  gm <- make_gene_matrix(metric)
  mat <- gm$mat
  if (metric %in% c("rs_core_cpm", "rs_rate")) mat <- log10(mat + 1)
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
  fwrite(
    as.data.table(cor_mat, keep.rownames = "sample_label"),
    file.path(OUT_DIR, paste0("sample_spearman_correlation_matrix_", metric, ".csv"))
  )
  cor_long <- as.data.table(as.table(cor_mat))
  setnames(cor_long, c("sample_1", "sample_2", "spearman_rho"))
  cor_long[, metric := metric]
  cor_long <- merge(cor_long, unique(dt[, .(sample_1 = sample_label, condition_1 = condition, replicate_1 = replicate)]), by = "sample_1")
  cor_long <- merge(cor_long, unique(dt[, .(sample_2 = sample_label, condition_2 = condition, replicate_2 = replicate)]), by = "sample_2")
  cor_long[, pair_type := fifelse(sample_1 == sample_2, "self",
    fifelse(condition_1 == condition_2, "same_condition", "different_condition"))]
  cor_summary[[metric]] <- cor_long[sample_1 < sample_2]

  p_cor <- ggplot(cor_long, aes(x = sample_1, y = sample_2, fill = spearman_rho)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2.6) +
    scale_fill_gradient2(low = "#B8323B", mid = "#F7F7F7", high = "#2C7A7B", midpoint = 0.5, limits = c(0, 1), name = "Spearman rho") +
    labs(
      title = paste0("Sample Correlation Matrix: ", metric),
      subtitle = "Gene-level medians; Spearman correlation",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  save_plot(p_cor, paste0("sample_correlation_heatmap_", metric), 8.5, 7.5)
}

cor_pairs <- rbindlist(cor_summary, fill = TRUE)
fwrite(cor_pairs, file.path(OUT_DIR, "sample_spearman_correlation_pairs_long.csv"))
cor_pair_summary <- cor_pairs[pair_type != "self", .(
  n_pairs = .N,
  median_spearman = median(spearman_rho, na.rm = TRUE),
  min_spearman = min(spearman_rho, na.rm = TRUE),
  max_spearman = max(spearman_rho, na.rm = TRUE)
), by = .(metric, pair_type)]
fwrite(cor_pair_summary, file.path(OUT_DIR, "sample_spearman_correlation_pair_summary.csv"))

message("QC complete. Outputs written to: ", OUT_DIR)
