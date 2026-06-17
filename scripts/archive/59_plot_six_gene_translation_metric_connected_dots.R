# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

infile <- analysis_path("Translation_indexes_fixed", "Gene_Level_Clean", "gene_level_clean_translation_metrics_all_samples.csv")
out_dir <- analysis_path("Translation_indexes_fixed", "Six_gene_metric_connected_dot_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes <- c("TAX1BP1", "TRA2A", "DEPP1", "HIVEP1", "ZNF266", "H2AC18")

metric_map <- data.table(
  metric = c("scanning_score", "ribosome_efficiency_score", "collision_score"),
  metric_label = c("Scanning", "Ribosome\nengagement", "Collision")
)

cell_cols <- c("Sensitive" = "#2C7FB8", "Resistant" = "#D7301F")
treatment_levels <- c("DMSO", "VCR")
condition_levels <- c("Sensitive_DMSO", "Sensitive_VCR", "Resistant_DMSO", "Resistant_VCR")

dt <- fread(infile)
if (!"sample_id" %in% names(dt)) dt[, sample_id := sample]
if (!"ribosome_engagement_score" %in% names(dt)) {
  dt[, ribosome_engagement_score := ribosome_efficiency_score]
}

dt[, treatment := fifelse(treatment %chin% c("Vin", "VIN", "VCR"), "VCR", treatment)]
dt[, cell_line := fifelse(cell_line %chin% c("Sensitive", "Resistant"), cell_line, cell_line_from_sample)]
dt[, condition := paste(cell_line, treatment, sep = "_")]
dt[, condition := factor(condition, levels = condition_levels)]
dt[, treatment := factor(treatment, levels = treatment_levels)]
dt[, cell_line := factor(cell_line, levels = c("Sensitive", "Resistant"))]

plot_dt <- dt[gene_name %chin% genes]
missing_genes <- setdiff(genes, unique(plot_dt$gene_name))
if (length(missing_genes)) {
  warning("Missing requested genes from input table: ", paste(missing_genes, collapse = ", "))
}

plot_long <- melt(
  plot_dt,
  id.vars = c("gene_name", "sample_id", "condition", "cell_line", "treatment", "replicate"),
  measure.vars = metric_map$metric,
  variable.name = "metric",
  value.name = "metric_value"
)
plot_long <- plot_long[is.finite(metric_value)]
plot_long <- merge(plot_long, metric_map, by = "metric", all.x = TRUE)
plot_long[, metric_label := factor(metric_label, levels = metric_map$metric_label)]
plot_long[, gene_name := factor(gene_name, levels = genes)]

mean_dt <- plot_long[, .(
  mean_value = mean(metric_value, na.rm = TRUE)
), by = .(gene_name, metric, metric_label, cell_line, treatment)]
mean_dt <- mean_dt[is.finite(mean_value)]

dodge <- position_dodge(width = 0.34)

make_gene_plot <- function(gene) {
  gene_points <- plot_long[gene_name == gene]
  gene_means <- mean_dt[gene_name == gene]

  ggplot() +
    geom_line(
      data = gene_means,
      aes(x = treatment, y = mean_value, group = cell_line, color = cell_line),
      position = dodge,
      linewidth = 0.55,
      alpha = 0.95,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    geom_point(
      data = gene_points,
      aes(x = treatment, y = metric_value, color = cell_line),
      position = position_jitterdodge(jitter.width = 0.055, jitter.height = 0, dodge.width = 0.34, seed = 11),
      size = 1.15,
      alpha = 0.62,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    geom_point(
      data = gene_means,
      aes(x = treatment, y = mean_value, color = cell_line),
      position = dodge,
      size = 2.55,
      alpha = 1,
      na.rm = TRUE,
      show.legend = TRUE
    ) +
    facet_wrap(~ metric_label, nrow = 1, scales = "free_y", drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_color_manual(values = cell_cols, limits = names(cell_cols), breaks = names(cell_cols), drop = FALSE) +
    labs(
      title = gene,
      x = "Treatment",
      y = "Metric value",
      color = "Cell line"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(size = 9.5),
      axis.text = element_text(size = 8.5),
      legend.position = "bottom",
      plot.margin = margin(4, 4, 4, 4)
    )
}

plots <- lapply(genes, make_gene_plot)
names(plots) <- genes

combined <- wrap_plots(plots, nrow = 2, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Translation metric scores across treatment and cell line",
    subtitle = "Small dots are biological replicates; large dots are group means; lines connect DMSO to VCR means within each cell line"
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 11)
  )

out_png <- file.path(out_dir, "six_gene_translation_metric_connected_dot_grid.png")
out_pdf <- file.path(out_dir, "six_gene_translation_metric_connected_dot_grid.pdf")
ggsave(out_png, combined, width = 16, height = 12, dpi = 300, bg = "white")
ggsave(out_pdf, combined, width = 16, height = 12, bg = "white")

fwrite(plot_long, file.path(out_dir, "six_gene_translation_metric_plot_long_data.csv"))
fwrite(mean_dt, file.path(out_dir, "six_gene_translation_metric_group_means.csv"))

cat("\nSix-gene translation metric connected dot plot complete.\n")
cat("Input:\n", infile, "\n", sep = "")
cat("PNG:\n", out_png, "\n", sep = "")
cat("PDF:\n", out_pdf, "\n", sep = "")
cat("\nRows by gene and condition:\n")
print(plot_long[, .N, by = .(gene_name, condition)][order(gene_name, condition)])
