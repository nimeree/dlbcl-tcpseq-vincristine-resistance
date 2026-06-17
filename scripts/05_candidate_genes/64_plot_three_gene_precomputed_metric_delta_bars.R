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

out_dir <- analysis_path("Figure_A_B_gene_metric_and_psite_panels")
limma_dir <- analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05", "Results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_dt <- data.table(
  gene_name = rep(c("TRA2A", "SEC24C", "MAPKBP1"), each = 12),
  contrast = rep(rep(c("Baseline", "Sensitive VCR", "Resistant VCR", "Interaction"), each = 3), times = 3),
  metric = rep(c("Scanning", "Ribosome engagement", "Collision"), times = 12),
  median_delta = c(
    -1.11,  5.29,  2.16,
    -0.32,  0.43,  0.14,
     1.81, -5.51, -2.12,
     2.13, -5.94, -2.26,
    -2.29,  4.20,  2.01,
    -2.64,  4.68,  1.25,
     2.69, -2.53, -1.39,
     5.33, -7.21, -2.65,
    -2.20,  3.38,  0.59,
    -2.23,  4.35,  0.74,
     2.22, -4.41, -0.78,
     4.45, -8.75, -1.52
  )
)

plot_dt[, gene_name := factor(gene_name, levels = c("TRA2A", "SEC24C", "MAPKBP1"))]
plot_dt[, contrast := factor(contrast, levels = c("Baseline", "Sensitive VCR", "Resistant VCR", "Interaction"))]
plot_dt[, metric := factor(metric, levels = c("Scanning", "Ribosome engagement", "Collision"))]

limma_map <- data.table(
  metric = c("Scanning", "Ribosome engagement", "Collision"),
  metric_dir = c("scanning_score", "ribosome_efficiency_score", "collision_score")
)
contrast_map <- data.table(
  contrast = c("Baseline", "Sensitive VCR", "Resistant VCR", "Interaction"),
  limma_contrast = c("Resistance_baseline", "VCR_sensitive", "VCR_resistant", "Interaction")
)

pvalue_dt <- rbindlist(lapply(seq_len(nrow(limma_map)), function(i) {
  rbindlist(lapply(seq_len(nrow(contrast_map)), function(j) {
    f <- file.path(
      limma_dir,
      limma_map$metric_dir[i],
      paste0(contrast_map$limma_contrast[j], "_limma_all_genes.csv")
    )
    if (!file.exists(f)) stop("Missing limma table: ", f)
    d <- fread(f)
    d[gene_name %chin% c("TRA2A", "SEC24C", "MAPKBP1"), .(
      gene_name,
      metric = limma_map$metric[i],
      contrast = contrast_map$contrast[j],
      limma_logFC = logFC,
      raw_P = P.Value
    )]
  }), use.names = TRUE)
}), use.names = TRUE)

stars_from_p <- function(p) {
  fifelse(is.na(p), "NA",
    fifelse(p < 0.001, "***",
      fifelse(p < 0.01, "**",
        fifelse(p < 0.05, "*", "ns")
      )
    )
  )
}

plot_dt <- merge(plot_dt, pvalue_dt, by = c("gene_name", "metric", "contrast"), all.x = TRUE)
plot_dt[, gene_name := factor(gene_name, levels = c("TRA2A", "SEC24C", "MAPKBP1"))]
plot_dt[, contrast := factor(contrast, levels = c("Baseline", "Sensitive VCR", "Resistant VCR", "Interaction"))]
plot_dt[, metric := factor(metric, levels = c("Scanning", "Ribosome engagement", "Collision"))]
plot_dt[, p_stars := stars_from_p(raw_P)]
plot_dt[, direction := fifelse(median_delta > 0, "Up", fifelse(median_delta < 0, "Down", "No change"))]
plot_dt[, annotation := p_stars]
plot_dt[, label_y := median_delta + fifelse(median_delta >= 0, 0.35, -0.35)]
plot_dt[, label_vjust := fifelse(median_delta >= 0, 0, 1)]

metric_cols <- c(
  "Scanning" = "#4393C3",
  "Ribosome engagement" = "#D6604D",
  "Collision" = "#4DAC26"
)
figure_title <- "Multi-stage translational responses of representative\nTCP-seq candidates across experimental contrasts"

p <- ggplot(plot_dt, aes(x = contrast, y = median_delta, fill = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    color = "grey25",
    linewidth = 0.18
  ) +
  geom_text(
    aes(y = label_y, label = annotation, group = metric, vjust = label_vjust),
    position = position_dodge(width = 0.78),
    size = 2.6,
    lineheight = 0.86,
    color = "grey15"
  ) +
  facet_wrap(~ gene_name, nrow = 1) +
  scale_fill_manual(values = metric_cols, breaks = names(metric_cols), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.13, 0.13))) +
  labs(
    title = figure_title,
    x = NULL,
    y = "Median delta",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 12)
  )

out_png <- file.path(out_dir, "three_gene_precomputed_metric_delta_bars.png")
out_pdf <- file.path(out_dir, "three_gene_precomputed_metric_delta_bars.pdf")
out_word_png <- file.path(out_dir, "three_gene_precomputed_metric_delta_bars_word_layout.png")
out_word_pdf <- file.path(out_dir, "three_gene_precomputed_metric_delta_bars_word_layout.pdf")
out_csv <- file.path(out_dir, "three_gene_precomputed_metric_delta_values.csv")
out_p_csv <- file.path(out_dir, "three_gene_precomputed_metric_delta_values_with_limma_pvalues.csv")

ggsave(out_png, p, width = 14, height = 5, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 14, height = 5, bg = "white")

y_min <- -10
y_max <- max(6.5, max(c(plot_dt$median_delta, plot_dt$label_y), na.rm = TRUE) + 0.45)

make_gene_plot <- function(gene) {
  ggplot(plot_dt[gene_name == gene], aes(x = contrast, y = median_delta, fill = metric)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.42) +
    geom_col(
      position = position_dodge(width = 0.78),
      width = 0.68,
      color = "grey25",
      linewidth = 0.16
    ) +
    geom_text(
      aes(y = label_y, label = annotation, group = metric, vjust = label_vjust),
      position = position_dodge(width = 0.78),
      size = 2.6,
      color = "grey15"
    ) +
    scale_fill_manual(values = metric_cols, breaks = names(metric_cols), drop = FALSE) +
    scale_y_continuous(limits = c(y_min, y_max), expand = expansion(mult = c(0.02, 0.02))) +
    labs(title = gene, x = NULL, y = "Median delta", fill = NULL) +
    theme_minimal(base_size = 10.5) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
      axis.title.y = element_text(size = 10.5),
      plot.margin = margin(4, 8, 4, 8)
    )
}

p_tra2a <- make_gene_plot("TRA2A")
p_sec24c <- make_gene_plot("SEC24C")
p_mapkbp1 <- make_gene_plot("MAPKBP1")

word_layout <- (
  p_tra2a + p_sec24c +
    plot_layout(ncol = 2)
) / (
  plot_spacer() + p_mapkbp1 + plot_spacer() +
    plot_layout(ncol = 3, widths = c(0.5, 1, 0.5))
) +
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(title = figure_title) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13, lineheight = 0.95)
  )

ggsave(out_word_png, word_layout, width = 7.3, height = 8.2, dpi = 300, bg = "white")
ggsave(out_word_pdf, word_layout, width = 7.3, height = 8.2, bg = "white")
fwrite(plot_dt, out_csv)
fwrite(plot_dt, out_p_csv)

cat("\nPrecomputed metric delta bar chart complete.\n")
cat("PNG:\n", out_png, "\n", sep = "")
cat("PDF:\n", out_pdf, "\n", sep = "")
cat("Word-layout PNG:\n", out_word_png, "\n", sep = "")
cat("Word-layout PDF:\n", out_word_pdf, "\n", sep = "")
cat("Data:\n", out_csv, "\n", sep = "")
cat("Data with limma P-values:\n", out_p_csv, "\n", sep = "")
