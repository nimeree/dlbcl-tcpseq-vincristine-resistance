# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- analysis_path()
metric_file <- file.path(base_dir, "Translation_indexes_fixed", "Gene_Level_Clean", "gene_level_clean_translation_metrics_all_samples.csv")
out_dir <- file.path(base_dir, "Figure_A_B_gene_metric_and_psite_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes <- c("TAX1BP1", "TRA2A", "DEPP1", "HIVEP1", "ZNF266", "H2AC18")
metric_map <- data.table(
  metric = c("scanning_score", "ribosome_efficiency_score", "collision_score"),
  metric_label = c("Scanning", "Ribosome\nengagement", "Collision")
)

dt <- fread(metric_file)
dt[, treatment := fifelse(grepl("Vin|VCR", condition, ignore.case = TRUE), "VCR", "DMSO")]
dt[, cell_line := fifelse(grepl("^Resistant", condition), "Resistant", "Sensitive")]

plot_long <- melt(
  dt[gene_name %chin% genes & treatment == "DMSO"],
  id.vars = c("gene_name", "sample", "condition", "cell_line", "treatment"),
  measure.vars = metric_map$metric,
  variable.name = "metric",
  value.name = "value"
)
plot_long <- merge(plot_long, metric_map, by = "metric", all.x = TRUE)

med_dt <- plot_long[is.finite(value), .(
  median_value = median(value, na.rm = TRUE),
  n_replicates = .N
), by = .(gene_name, metric, metric_label, cell_line)]

wide <- dcast(
  med_dt,
  gene_name + metric + metric_label ~ cell_line,
  value.var = "median_value"
)
wide_n <- dcast(
  med_dt,
  gene_name + metric + metric_label ~ cell_line,
  value.var = "n_replicates"
)
setnames(wide_n, intersect(c("Sensitive", "Resistant"), names(wide_n)), paste0(intersect(c("Sensitive", "Resistant"), names(wide_n)), "_n"))

delta_dt <- merge(wide, wide_n, by = c("gene_name", "metric", "metric_label"), all = TRUE)
delta_dt[, median_delta := Resistant - Sensitive]

all_combos <- CJ(gene_name = genes, metric = metric_map$metric)
all_combos <- merge(all_combos, metric_map, by = "metric", all.x = TRUE)
delta_dt <- merge(all_combos, delta_dt, by = c("gene_name", "metric", "metric_label"), all.x = TRUE)
delta_dt[, gene_name := factor(gene_name, levels = genes)]
delta_dt[, metric_label := factor(metric_label, levels = metric_map$metric_label)]
delta_dt[, direction := fifelse(
  is.na(median_delta), "NA",
  fifelse(median_delta < 0, "Sensitive direction", "Resistant direction")
)]

bar_cols <- c(
  "Sensitive direction" = "#2C7FB8",
  "Resistant direction" = "#D7301F",
  "NA" = "grey82"
)

na_dt <- delta_dt[is.na(median_delta)]

p <- ggplot(delta_dt, aes(x = metric_label, y = median_delta, fill = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  geom_col(width = 0.68, color = "grey30", linewidth = 0.2, na.rm = TRUE) +
  geom_text(
    data = na_dt,
    aes(x = metric_label, y = 0, label = "NA"),
    inherit.aes = FALSE,
    color = "grey35",
    fontface = "bold",
    size = 3.4,
    vjust = -0.45
  ) +
  facet_wrap(~ gene_name, nrow = 2, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = bar_cols, breaks = c("Sensitive direction", "Resistant direction"), drop = FALSE) +
  labs(
    title = "Baseline translation metric shifts by gene",
    subtitle = "Median delta = Resistant DMSO median minus Sensitive DMSO median",
    x = NULL,
    y = "Median delta (Resistant DMSO minus Sensitive DMSO)",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 17),
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12)
  )

out_png <- file.path(out_dir, "six_gene_resistant_baseline_metric_delta_bars.png")
out_pdf <- file.path(out_dir, "six_gene_resistant_baseline_metric_delta_bars.pdf")
ggsave(out_png, p, width = 14, height = 10, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 14, height = 10, bg = "white")

fwrite(delta_dt, file.path(out_dir, "six_gene_resistant_baseline_metric_delta_values.csv"))

cat("\nSix-gene baseline metric delta bar plot complete.\n")
cat("PNG:\n", out_png, "\n", sep = "")
cat("PDF:\n", out_pdf, "\n", sep = "")
cat("\nDelta values:\n")
print(delta_dt[order(gene_name, metric_label), .(
  gene_name,
  metric_label,
  Sensitive_median = Sensitive,
  Resistant_median = Resistant,
  median_delta,
  Sensitive_n,
  Resistant_n
)])
