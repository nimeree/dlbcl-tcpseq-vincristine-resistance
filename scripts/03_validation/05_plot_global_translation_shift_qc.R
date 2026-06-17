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
OUT_DIR <- analysis_path("Translation_indexes_fixed", "Global_Shift_QC")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, width, height) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = width, height = height)
}

parse_sample_meta <- function(sample) {
  x <- data.table(sample = sample)
  x[, cell_line := fifelse(grepl("^SU8R", sample, ignore.case = TRUE), "Resistant", "Sensitive")]
  x[, treatment_raw := fifelse(grepl("-Vin-", sample, ignore.case = TRUE), "Vin", "DMSO")]
  x[, treatment := fifelse(treatment_raw == "Vin", "VCR", "DMSO")]
  x[, replicate := fifelse(grepl("Rep1", sample, ignore.case = TRUE), "Rep1",
    fifelse(grepl("Rep2", sample, ignore.case = TRUE), "Rep2", NA_character_))]
  x[, condition := paste(cell_line, treatment, sep = "_")]
  x[, sample_label := paste(condition, replicate, sep = "_")]
  x
}

metric_labels <- c(
  ribosome_efficiency_score = "Ribosome efficiency score",
  protein_output_score = "Protein output score"
)

format_p <- function(p) {
  fifelse(
    is.na(p),
    "p = NA",
    fifelse(
      p < 2.2e-16,
      "p < 2.2e-16",
      fifelse(
        p < 0.001,
        paste0("p = ", formatC(p, format = "e", digits = 2)),
        paste0("p = ", signif(p, 3))
      )
    )
  )
}

format_delta <- function(x) {
  paste0("median shift = ", sprintf("%.3f", x))
}

message("Loading data")
dt <- fread(INFILE, select = c(
  "sample", "gene_id_clean", "gene_name", "ribosome_efficiency_score", "protein_output_score"
))
for (cc in names(metric_labels)) dt[, (cc) := as.numeric(get(cc))]
dt <- dt[!is.na(gene_id_clean) & gene_id_clean != ""]

meta <- parse_sample_meta(unique(dt$sample))
dt <- merge(dt, meta, by = "sample", all.x = TRUE)
sample_levels <- meta[order(cell_line, treatment, replicate), sample_label]
condition_levels <- c("Sensitive_DMSO", "Sensitive_VCR", "Resistant_DMSO", "Resistant_VCR")
dt[, sample_label := factor(sample_label, levels = sample_levels)]
dt[, condition := factor(condition, levels = condition_levels)]
dt[, treatment := factor(treatment, levels = c("DMSO", "VCR"))]
dt[, cell_line := factor(cell_line, levels = c("Sensitive", "Resistant"))]

# Collapse transcript/isoform rows to gene x sample medians to avoid isoform pseudo-replication.
gene_sample <- dt[, .(
  ribosome_efficiency_score = median(ribosome_efficiency_score, na.rm = TRUE),
  protein_output_score = median(protein_output_score, na.rm = TRUE),
  gene_name = names(sort(table(gene_name), decreasing = TRUE))[1]
), by = .(sample, sample_label, cell_line, treatment, condition, replicate, gene_id_clean)]

long <- melt(
  gene_sample,
  id.vars = c("sample", "sample_label", "cell_line", "treatment", "condition", "replicate", "gene_id_clean", "gene_name"),
  measure.vars = names(metric_labels),
  variable.name = "metric",
  value.name = "value"
)
long <- long[is.finite(value)]
long[, metric_label := factor(metric_labels[metric], levels = metric_labels)]

message("Running Wilcoxon tests")
contrast_table <- data.table(
  contrast = c(
    "Sensitive VCR vs Sensitive DMSO",
    "Resistant VCR vs Resistant DMSO",
    "Resistant DMSO vs Sensitive DMSO",
    "Resistant VCR vs Sensitive VCR"
  ),
  group_a = c("Sensitive_VCR", "Resistant_VCR", "Resistant_DMSO", "Resistant_VCR"),
  group_b = c("Sensitive_DMSO", "Resistant_DMSO", "Sensitive_DMSO", "Sensitive_VCR"),
  biological_question = c(
    "VCR suppression in sensitive cells",
    "VCR suppression in resistant cells",
    "Baseline resistant vs sensitive difference",
    "Resistant vs sensitive difference under VCR"
  ),
  alternative = c("less", "less", "greater", "greater")
)

auc_calc <- function(a, b) {
  ranks <- rank(c(a, b), ties.method = "average")
  n1 <- as.numeric(length(a))
  n0 <- as.numeric(length(b))
  (sum(ranks[seq_len(length(a))]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_one <- function(metric_name, group_a, group_b, alternative) {
  a <- gene_sample[condition == group_a & is.finite(get(metric_name)), get(metric_name)]
  b <- gene_sample[condition == group_b & is.finite(get(metric_name)), get(metric_name)]
  p <- suppressWarnings(wilcox.test(a, b, alternative = alternative, exact = FALSE)$p.value)
  auc <- auc_calc(a, b)
  data.table(
    metric = metric_name,
    n_group_a = length(a),
    n_group_b = length(b),
    median_group_a = median(a, na.rm = TRUE),
    median_group_b = median(b, na.rm = TRUE),
    median_difference_a_minus_b = median(a, na.rm = TRUE) - median(b, na.rm = TRUE),
    wilcox_alternative = alternative,
    wilcox_p = p,
    auc_group_a_gt_group_b = auc,
    rank_biserial = 2 * auc - 1
  )
}

tests <- rbindlist(lapply(seq_len(nrow(contrast_table)), function(i) {
  row <- contrast_table[i]
  out <- rbindlist(lapply(names(metric_labels), function(metric_name) {
    test_one(metric_name, row$group_a, row$group_b, row$alternative)
  }))
  cbind(row[, .(contrast, group_a, group_b, biological_question)], out)
}), fill = TRUE)
tests[, metric_label := metric_labels[metric]]
setcolorder(tests, c("contrast", "biological_question", "metric_label", "group_a", "group_b"))

message("Plotting faceted density overlays")
density_ann <- tests[contrast %in% c(
  "Sensitive VCR vs Sensitive DMSO",
  "Resistant VCR vs Resistant DMSO"
)]
density_ann[, cell_line := fifelse(grepl("^Sensitive", group_a), "Sensitive", "Resistant")]
density_ann[, metric_label := factor(metric_label, levels = metric_labels)]
density_ann[, cell_line := factor(cell_line, levels = c("Sensitive", "Resistant"))]
density_ann[, label := paste0(
  gsub(" VCR vs .*", " VCR vs DMSO", contrast),
  "\n", format_p(wilcox_p),
  "\n", format_delta(median_difference_a_minus_b)
)]
density_pos <- long[, .(
  x = quantile(value, 0.05, na.rm = TRUE),
  y = Inf
), by = .(metric_label, cell_line)]
density_ann <- merge(density_ann, density_pos, by = c("metric_label", "cell_line"), all.x = TRUE)

p_density <- ggplot(long, aes(x = value, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.12, linewidth = 0.8) +
  geom_label(
    data = density_ann,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.08,
    size = 3.2,
    fill = "white",
    linewidth = 0.2
  ) +
  facet_grid(metric_label ~ cell_line, scales = "free") +
  scale_color_manual(values = c("DMSO" = "#2C7A7B", "VCR" = "#B8323B")) +
  scale_fill_manual(values = c("DMSO" = "#2C7A7B", "VCR" = "#B8323B")) +
  labs(
    title = "Global Translation Metric Distributions by Treatment and Cell Line",
    subtitle = "Gene-sample medians; VCR curves overlaid against DMSO within each cell line",
    x = "Metric value",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
save_plot(p_density, "density_ribosome_efficiency_and_protein_output_score_by_cell_line_treatment", 10.5, 7.2)

message("Plotting individual-sample violins/boxplots")
sample_medians <- long[, .(sample_median = median(value, na.rm = TRUE)), by = .(sample_label, condition, metric_label)]
sample_medians[, label := sprintf("%.2f", sample_median)]
p_violin <- ggplot(long, aes(x = sample_label, y = value, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.72, linewidth = 0.2) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9, linewidth = 0.22) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 7, color = "black") +
  geom_text(
    data = sample_medians,
    aes(x = sample_label, y = sample_median, label = label),
    inherit.aes = FALSE,
    size = 3.1,
    vjust = -0.75,
    color = "black"
  ) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 1) +
  labs(
    title = "Per-Sample Metric Distributions",
    subtitle = "All eight samples shown individually; black ticks mark medians",
    x = NULL,
    y = "Metric value"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none", panel.grid.minor = element_blank())
save_plot(p_violin, "violin_boxplot_ribosome_efficiency_and_protein_output_score_by_sample", 10.5, 8.2)

message("Plotting empirical CDFs")
ecdf_ann <- copy(tests)
ecdf_ann[, metric_label := factor(metric_label, levels = metric_labels)]
ecdf_ann[, panel_id := as.integer(metric_label)]
ecdf_ann[, x := -Inf]
ecdf_ann[, y := fifelse(
  contrast == "Sensitive VCR vs Sensitive DMSO", 0.22,
  fifelse(contrast == "Resistant VCR vs Resistant DMSO", 0.38,
    fifelse(contrast == "Resistant DMSO vs Sensitive DMSO", 0.54, 0.70)
  )
)]
ecdf_ann[, label := paste0(
  contrast,
  "\n", format_p(wilcox_p),
  "; ", format_delta(median_difference_a_minus_b)
)]
p_ecdf <- ggplot(long, aes(x = value, color = condition)) +
  stat_ecdf(linewidth = 0.85) +
  geom_label(
    data = ecdf_ann,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.02,
    vjust = 0.5,
    size = 2.85,
    fill = "white",
    linewidth = 0.2
  ) +
  facet_wrap(~ metric_label, scales = "free", ncol = 1) +
  scale_color_manual(values = c(
    "Sensitive_DMSO" = "#2C7A7B",
    "Sensitive_VCR" = "#B8323B",
    "Resistant_DMSO" = "#2B6CB0",
    "Resistant_VCR" = "#805AD5"
  )) +
  labs(
    title = "Empirical CDFs Reveal Global Distribution Shifts",
    subtitle = "Left-shifted curves indicate lower metric values across the gene distribution",
    x = "Metric value",
    y = "Cumulative fraction of genes",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
save_plot(p_ecdf, "ecdf_ribosome_efficiency_and_protein_output_score_by_condition", 9, 7.5)

message("Plotting DMSO-only baseline comparison")
dmso <- long[treatment == "DMSO"]
dmso_ann <- tests[contrast == "Resistant DMSO vs Sensitive DMSO"]
dmso_ann[, metric_label := factor(metric_label, levels = metric_labels)]
dmso_ann[, x := 1.5]
dmso_ann[, y := long[treatment == "DMSO", .(y = quantile(value, 0.97, na.rm = TRUE)), by = metric_label]$y]
dmso_ann[, label := paste0("Resistant vs Sensitive\n", format_p(wilcox_p), "\n", format_delta(median_difference_a_minus_b))]
p_dmso <- ggplot(dmso, aes(x = cell_line, y = value, fill = cell_line)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.72, linewidth = 0.2) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9, linewidth = 0.22) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 8, color = "black") +
  geom_label(
    data = dmso_ann,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    fill = "white",
    linewidth = 0.2
  ) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Sensitive" = "#2C7A7B", "Resistant" = "#2B6CB0")) +
  labs(
    title = "Baseline Translation Differences in DMSO",
    subtitle = "Sensitive vs resistant cells before drug exposure",
    x = NULL,
    y = "Metric value"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
save_plot(p_dmso, "dmso_baseline_sensitive_vs_resistant_boxplot", 8.5, 5.2)

fwrite(tests, file.path(OUT_DIR, "wilcoxon_global_shift_tests.csv"))

summary_by_condition <- long[, .(
  n_gene_sample_values = .N,
  median = median(value, na.rm = TRUE),
  mean = mean(value, na.rm = TRUE),
  q25 = quantile(value, 0.25, na.rm = TRUE),
  q75 = quantile(value, 0.75, na.rm = TRUE)
), by = .(metric, metric_label, condition, cell_line, treatment)]
fwrite(summary_by_condition, file.path(OUT_DIR, "metric_summary_by_condition.csv"))

sample_summary <- long[, .(
  n_gene_values = .N,
  median = median(value, na.rm = TRUE),
  mean = mean(value, na.rm = TRUE),
  q25 = quantile(value, 0.25, na.rm = TRUE),
  q75 = quantile(value, 0.75, na.rm = TRUE)
), by = .(metric, metric_label, sample, sample_label, condition, replicate)]
fwrite(sample_summary, file.path(OUT_DIR, "metric_summary_by_sample.csv"))

message("Done. Outputs written to: ", OUT_DIR)
