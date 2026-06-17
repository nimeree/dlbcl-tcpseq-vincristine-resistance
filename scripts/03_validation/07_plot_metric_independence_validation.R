# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

INFILE <- analysis_path("Translation_indexes_fixed", "Gene_Level_Clean", "gene_level_clean_translation_metrics_all_samples.csv")
OUT_DIR <- analysis_path("Translation_indexes_fixed", "Validation_Plots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, width, height) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = width, height = height)
}

format_p <- function(p) {
  fifelse(
    is.na(p),
    "p = NA",
    fifelse(p < 1e-4, "p < 0.0001", paste0("p = ", signif(p, 3)))
  )
}

mode_first <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

dt <- fread(INFILE)
dt <- dt[baseline_cpm_line > 0]

gene_level <- dt[, .(
  gene_name = mode_first(gene_name),
  ribosome_efficiency_score = median(ribosome_efficiency_score, na.rm = TRUE),
  protein_output_score = median(protein_output_score, na.rm = TRUE),
  collision_score = median(collision_score, na.rm = TRUE),
  scanning_score = median(scanning_score, na.rm = TRUE),
  complete_te_protein_set = all(complete_te_protein_set)
), by = gene_id_clean]

plot_dt <- gene_level[complete_te_protein_set == TRUE]

top_like_genes <- unique(c(
  grep("^RPL|^RPS", plot_dt$gene_name, value = TRUE),
  "EEF1A1", "EEF1B2", "EEF1D", "EEF1G", "EEF2",
  "EIF3A", "EIF3B", "EIF3C", "EIF3D", "EIF3E", "EIF3F", "EIF3G", "EIF3H", "EIF3I", "EIF3J", "EIF3K", "EIF3L", "EIF3M",
  "PABPC1", "NPM1", "TPT1"
))
top_like_genes <- setdiff(top_like_genes, c("TOP1", "TOP2A", "TOP2B", "TOP3A", "TOP3B", "TOPBP1", "TOPORS"))

plot_dt[, top_group := fifelse(gene_name %chin% top_like_genes, "TOP-like mRNA genes", "Background")]
plot_dt[, top_group := factor(top_group, levels = c("Background", "TOP-like mRNA genes"))]

long <- melt(
  plot_dt,
  id.vars = c("gene_id_clean", "gene_name", "top_group"),
  measure.vars = c("ribosome_efficiency_score", "protein_output_score"),
  variable.name = "metric",
  value.name = "value"
)
long <- long[is.finite(value)]
long[, metric_label := factor(
  fifelse(metric == "ribosome_efficiency_score", "Ribosome efficiency score", "Protein output score"),
  levels = c("Ribosome efficiency score", "Protein output score")
)]

stats <- rbindlist(lapply(unique(long$metric), function(mm) {
  x <- long[metric == mm & top_group == "TOP-like mRNA genes", value]
  y <- long[metric == mm & top_group == "Background", value]
  alt <- if (mm == "ribosome_efficiency_score") "less" else "greater"
  data.table(
    metric = mm,
    metric_label = if (mm == "ribosome_efficiency_score") "Ribosome efficiency score" else "Protein output score",
    n_top = length(x),
    n_background = length(y),
    median_top = median(x, na.rm = TRUE),
    median_background = median(y, na.rm = TRUE),
    p = suppressWarnings(wilcox.test(x, y, alternative = alt, exact = FALSE)$p.value)
  )
}))
stats[, metric_label := factor(metric_label, levels = levels(long$metric_label))]
stats[, label := paste0(
  "TOP median = ", sprintf("%.2f", median_top),
  "\nBackground median = ", sprintf("%.2f", median_background),
  "\n", format_p(p),
  "\nn = ", n_top, " vs ", n_background
)]
ann_pos <- long[, .(
  x = 1.5,
  y = quantile(value, 0.985, na.rm = TRUE)
), by = metric_label]
stats <- merge(stats, ann_pos, by = "metric_label", all.x = TRUE)

p_top <- ggplot(long, aes(x = top_group, y = value, fill = top_group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.72, linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.92, linewidth = 0.25) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 8, color = "black") +
  geom_label(
    data = stats,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.4,
    fill = "white",
    linewidth = 0.2
  ) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Background" = "#A7B0BA", "TOP-like mRNA genes" = "#2C7A7B")) +
  labs(
    title = "TOP-like mRNA Genes Validate Distinct Translation Metrics",
    subtitle = "Gene-level medians; TOP-like set excludes topoisomerase genes",
    x = NULL,
    y = "Metric value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 18, hjust = 1),
    panel.grid.minor = element_blank()
  )
save_plot(p_top, "top_like_mrna_metric_validation", 10, 5.8)

cor_dt <- gene_level[complete_te_protein_set == TRUE, .(
  ribosome_efficiency_score,
  protein_output_score,
  collision_score,
  scanning_score
)]
cor_mat <- cor(cor_dt, use = "pairwise.complete.obs", method = "spearman")
metric_names <- c(
  ribosome_efficiency_score = "Ribosome\nefficiency",
  protein_output_score = "Protein\noutput",
  collision_score = "Collision",
  scanning_score = "Scanning"
)
cor_long <- as.data.table(as.table(cor_mat))
setnames(cor_long, c("metric_x", "metric_y", "spearman_rho"))
cor_long[, `:=`(
  metric_x_label = factor(metric_names[as.character(metric_x)], levels = metric_names),
  metric_y_label = factor(metric_names[as.character(metric_y)], levels = rev(metric_names)),
  label = sprintf("%.2f", spearman_rho)
)]

p_cor <- ggplot(cor_long, aes(x = metric_x_label, y = metric_y_label, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#B8323B",
    mid = "#F7F7F7",
    high = "#2C7A7B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  coord_fixed() +
  labs(
    title = "Translation Metrics Capture Distinct Biological Axes",
    subtitle = "Gene-level Spearman correlations across complete TE/protein-output genes",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right"
  )
save_plot(p_cor, "translation_metric_independence_correlation_heatmap", 7.2, 5.8)

fwrite(stats[, .(
  metric = as.character(metric_label),
  n_top,
  n_background,
  median_top,
  median_background,
  wilcox_p = p
)], file.path(OUT_DIR, "top_like_mrna_metric_validation_stats.csv"))
fwrite(as.data.table(cor_mat, keep.rownames = "metric"), file.path(OUT_DIR, "translation_metric_independence_spearman_matrix.csv"))

cat("Wrote TOP validation and metric independence plots to:", OUT_DIR, "\n")
print(stats[, .(metric = as.character(metric_label), n_top, n_background, median_top, median_background, p)])
print(round(cor_mat, 3))
