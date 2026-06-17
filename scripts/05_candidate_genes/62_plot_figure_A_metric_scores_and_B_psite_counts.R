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

base_dir <- analysis_path()
metric_file <- file.path(base_dir, "Translation_indexes_fixed", "Gene_Level_Clean", "gene_level_clean_translation_metrics_all_samples.csv")
psite_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
psite_file <- file.path(psite_dir, "psite_gene_counts_long_by_fraction_sample.csv")
out_dir <- file.path(base_dir, "Figure_A_B_gene_metric_and_psite_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_cols <- c("Sensitive" = "#2C7FB8", "Resistant" = "#D7301F")
treatment_levels <- c("DMSO", "VCR")
cell_levels <- c("Sensitive", "Resistant")
dodge_width <- 0.34

standardize_sample_fields <- function(dt) {
  if (!"sample_id" %in% names(dt)) dt[, sample_id := sample]
  dt[, treatment := fifelse(grepl("Vin|VCR", condition, ignore.case = TRUE), "VCR", "DMSO")]
  dt[, cell_line := fifelse(grepl("^Resistant", condition), "Resistant", "Sensitive")]
  dt[, treatment := factor(treatment, levels = treatment_levels)]
  dt[, cell_line := factor(cell_line, levels = cell_levels)]
  dt
}

make_connected_plot <- function(point_dt, mean_dt, facet_col, title, y_lab) {
  dodge <- position_dodge(width = dodge_width)
  ggplot() +
    geom_line(
      data = mean_dt,
      aes(x = treatment, y = mean_value, group = cell_line, color = cell_line),
      position = dodge,
      linewidth = 0.6,
      alpha = 0.95,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    geom_point(
      data = point_dt,
      aes(x = treatment, y = value, color = cell_line),
      position = position_jitterdodge(
        jitter.width = 0.055,
        jitter.height = 0,
        dodge.width = dodge_width,
        seed = 11
      ),
      size = 1.15,
      alpha = 0.55,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    geom_point(
      data = mean_dt,
      aes(x = treatment, y = mean_value, color = cell_line),
      position = dodge,
      size = 2.65,
      alpha = 1,
      na.rm = TRUE,
      show.legend = TRUE
    ) +
    facet_wrap(as.formula(paste("~", facet_col)), nrow = 1, scales = "free_y", drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_color_manual(values = cell_cols, limits = names(cell_cols), breaks = names(cell_cols), drop = FALSE) +
    labs(title = title, x = "Treatment", y = y_lab, color = "Cell line") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      strip.text = element_text(face = "bold", size = 8.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(size = 9.5),
      axis.text = element_text(size = 8.5),
      legend.position = "bottom",
      plot.margin = margin(4, 4, 4, 4)
    )
}

make_row_figure <- function(plots, title, subtitle = NULL) {
  wrap_plots(plots, nrow = 1, ncol = length(plots), guides = "collect") +
    plot_annotation(title = title, subtitle = subtitle) &
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 10.5)
    )
}

# Figure A: composite translation metric scores.
figure_a_genes <- c("TAX1BP1", "TRA2A", "DEPP1")
metric_map <- data.table(
  metric = c("scanning_score", "ribosome_efficiency_score", "collision_score"),
  metric_label = c("Scanning", "Ribosome\nengagement", "Collision")
)

metric_dt <- fread(metric_file)
metric_dt <- standardize_sample_fields(metric_dt)
metric_dt <- metric_dt[gene_name %chin% figure_a_genes]

metric_long <- melt(
  metric_dt,
  id.vars = c("gene_name", "sample_id", "cell_line", "treatment"),
  measure.vars = metric_map$metric,
  variable.name = "metric",
  value.name = "value"
)
metric_long <- merge(metric_long, metric_map, by = "metric", all.x = TRUE)
metric_long[, metric_label := factor(metric_label, levels = metric_map$metric_label)]
metric_long[, gene_name := factor(gene_name, levels = figure_a_genes)]
metric_mean <- metric_long[is.finite(value), .(
  mean_value = mean(value, na.rm = TRUE)
), by = .(gene_name, metric_label, cell_line, treatment)]

figure_a_plots <- lapply(figure_a_genes, function(g) {
  make_connected_plot(
    metric_long[gene_name == g],
    metric_mean[gene_name == g],
    "metric_label",
    g,
    "Score"
  )
})

figure_a <- make_row_figure(
  figure_a_plots,
  "Figure A. Translation metric scores",
  "Small dots are biological replicates; large dots are group means; lines connect DMSO to VCR means within each cell line"
)

figure_a_png <- file.path(out_dir, "Figure_A_TAX1BP1_TRA2A_DEPP1_metric_scores.png")
figure_a_pdf <- file.path(out_dir, "Figure_A_TAX1BP1_TRA2A_DEPP1_metric_scores.pdf")
ggsave(figure_a_png, figure_a, width = 16, height = 5, dpi = 300, bg = "white")
ggsave(figure_a_pdf, figure_a, width = 16, height = 5, bg = "white")

fwrite(metric_long, file.path(out_dir, "Figure_A_metric_scores_long_data.csv"))
fwrite(metric_mean, file.path(out_dir, "Figure_A_metric_scores_group_means.csv"))

# Figure B: fraction-level P-site counts, normalized as log2 CPM + 1.
figure_b_genes <- c("HIVEP1", "ZNF266", "H2AC18")
fraction_levels <- c("SSU", "RS", "DS")

psite_dt <- fread(psite_file)
psite_dt <- standardize_sample_fields(psite_dt)
psite_dt[, fraction := factor(fraction, levels = fraction_levels)]

lib_dt <- rbindlist(lapply(fraction_levels, function(frac) {
  f <- file.path(psite_dir, paste0("Fraction_", frac), paste0(frac, "_library_sizes_after_gene_filter.csv"))
  if (!file.exists(f)) stop("Missing library-size file: ", f)
  d <- fread(f)
  d[, fraction := frac]
  d[, .(sample, fraction, lib_psites)]
}), use.names = TRUE)
lib_dt[, fraction := factor(fraction, levels = fraction_levels)]

psite_dt <- merge(psite_dt, lib_dt, by = c("sample", "fraction"), all.x = TRUE)
psite_dt[, normalised_count := log2((psite_count + 0.5) / pmax(lib_psites, 1) * 1e6 + 1)]
psite_dt <- psite_dt[gene_name %chin% figure_b_genes]
psite_dt[, gene_name := factor(gene_name, levels = figure_b_genes)]
setnames(psite_dt, "normalised_count", "value")

psite_mean <- psite_dt[is.finite(value), .(
  mean_value = mean(value, na.rm = TRUE)
), by = .(gene_name, fraction, cell_line, treatment)]

figure_b_plots <- lapply(figure_b_genes, function(g) {
  make_connected_plot(
    psite_dt[gene_name == g],
    psite_mean[gene_name == g],
    "fraction",
    g,
    "Normalised P-site count (log2)"
  )
})

figure_b <- make_row_figure(
  figure_b_plots,
  "Figure B. Fraction-resolved P-site counts",
  "Normalised count is log2(CPM + 1), using fraction/sample P-site library sizes"
)

figure_b_png <- file.path(out_dir, "Figure_B_HIVEP1_ZNF266_H2AC18_fraction_psite_counts.png")
figure_b_pdf <- file.path(out_dir, "Figure_B_HIVEP1_ZNF266_H2AC18_fraction_psite_counts.pdf")
ggsave(figure_b_png, figure_b, width = 16, height = 5, dpi = 300, bg = "white")
ggsave(figure_b_pdf, figure_b, width = 16, height = 5, bg = "white")

fwrite(psite_dt, file.path(out_dir, "Figure_B_fraction_psite_log2CPM_long_data.csv"))
fwrite(psite_mean, file.path(out_dir, "Figure_B_fraction_psite_log2CPM_group_means.csv"))

cat("\nFigure A and B exports complete.\n")
cat("Figure A PNG:\n", figure_a_png, "\n", sep = "")
cat("Figure B PNG:\n", figure_b_png, "\n", sep = "")
cat("\nFigure A finite replicate counts:\n")
print(metric_long[is.finite(value), .N, by = .(gene_name, metric_label, cell_line, treatment)][order(gene_name, metric_label, cell_line, treatment)])
cat("\nFigure B finite replicate counts:\n")
print(psite_dt[is.finite(value), .N, by = .(gene_name, fraction, cell_line, treatment)][order(gene_name, fraction, cell_line, treatment)])
