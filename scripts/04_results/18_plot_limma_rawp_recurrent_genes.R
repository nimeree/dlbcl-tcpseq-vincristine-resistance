# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# Recurrent gene summary across raw-p limma analyses
# - Input: Limma_translation_metrics_lfc0.7_rawP0.05 significant CSVs
# - Finds genes significant across multiple metric/contrast analyses
# - Writes tables plus recurrence plots
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
})

BASE_DIR <- analysis_path()
LIMMA_DIR <- file.path(BASE_DIR, "Limma_translation_metrics_lfc0.7_rawP0.05")
RESULT_DIR <- file.path(LIMMA_DIR, "Results")
OUT_DIR <- file.path(LIMMA_DIR, "Recurrent_gene_analysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

METRIC_ORDER <- c(
  "ribosome_efficiency_score",
  "protein_output_score",
  "collision_score",
  "scanning_score"
)

CONTRAST_ORDER <- c(
  "VCR_sensitive",
  "VCR_resistant",
  "Resistance_baseline",
  "Interaction"
)

METRIC_LABELS <- c(
  ribosome_efficiency_score = "Ribosome efficiency",
  protein_output_score = "Protein output",
  collision_score = "Collision",
  scanning_score = "Scanning"
)

CONTRAST_LABELS <- c(
  VCR_sensitive = "VCR sensitive",
  VCR_resistant = "VCR resistant",
  Resistance_baseline = "Baseline resistance",
  Interaction = "Interaction"
)

sig_files <- list.files(
  RESULT_DIR,
  pattern = "_limma_sig_rawP0\\.05_lfc0\\.7\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
stopifnot(length(sig_files) > 0)

read_sig_file <- function(path) {
  metric <- basename(dirname(path))
  contrast <- sub("_limma_sig_rawP0\\.05_lfc0\\.7\\.csv$", "", basename(path))
  dt <- fread(path)
  if (!nrow(dt)) return(data.table())
  dt[, `:=`(
    metric = metric,
    contrast = contrast,
    metric_label = METRIC_LABELS[metric],
    contrast_label = CONTRAST_LABELS[contrast],
    analysis = paste(metric, contrast, sep = "__"),
    analysis_label = paste(METRIC_LABELS[metric], CONTRAST_LABELS[contrast], sep = " | "),
    signed_score = sign(logFC) * -log10(pmax(P.Value, .Machine$double.xmin))
  )]
  dt[]
}

hits <- rbindlist(lapply(sig_files, read_sig_file), fill = TRUE)
hits <- hits[!is.na(gene_id_clean) & gene_id_clean != ""]
hits[, metric := factor(metric, levels = METRIC_ORDER)]
hits[, contrast := factor(contrast, levels = CONTRAST_ORDER)]
hits[, analysis_label := factor(
  analysis_label,
  levels = as.vector(outer(METRIC_LABELS[METRIC_ORDER], CONTRAST_LABELS[CONTRAST_ORDER], paste, sep = " | "))
)]

gene_summary <- hits[, .(
  gene_name = {
    x <- gene_name[!is.na(gene_name) & gene_name != ""]
    if (length(x)) names(sort(table(x), decreasing = TRUE))[1] else gene_id_clean[1]
  },
  n_analyses = uniqueN(analysis),
  n_metrics = uniqueN(metric),
  n_contrasts = uniqueN(contrast),
  metrics = paste(sort(unique(as.character(metric))), collapse = "; "),
  contrasts = paste(sort(unique(as.character(contrast))), collapse = "; "),
  analyses = paste(sort(unique(as.character(analysis_label))), collapse = "; "),
  max_abs_logFC = max(abs(logFC), na.rm = TRUE),
  best_raw_p = min(P.Value, na.rm = TRUE),
  best_fdr = min(adj.P.Val, na.rm = TRUE),
  mean_rank_score = mean(rank_score, na.rm = TRUE)
), by = gene_id_clean][order(-n_analyses, -n_metrics, -n_contrasts, best_raw_p)]

recurrent_genes <- gene_summary[n_analyses >= 2]
multi_metric_genes <- gene_summary[n_metrics >= 2]
multi_contrast_genes <- gene_summary[n_contrasts >= 2]

fwrite(hits, file.path(OUT_DIR, "all_significant_limma_rawP_hits_long.csv"))
fwrite(gene_summary, file.path(OUT_DIR, "gene_recurrence_summary_all_hits.csv"))
fwrite(recurrent_genes, file.path(OUT_DIR, "genes_seen_in_multiple_analyses.csv"))
fwrite(multi_metric_genes, file.path(OUT_DIR, "genes_seen_in_multiple_metrics.csv"))
fwrite(multi_contrast_genes, file.path(OUT_DIR, "genes_seen_in_multiple_contrasts.csv"))

analysis_counts <- hits[, .(gene_count = uniqueN(gene_id_clean)), by = .(metric, metric_label, contrast, contrast_label, analysis_label)]
fwrite(analysis_counts[order(metric, contrast)], file.path(OUT_DIR, "significant_gene_counts_by_analysis.csv"))

metric_overlap <- hits[, .(genes = list(unique(gene_id_clean))), by = metric]
metric_pairwise <- rbindlist(lapply(combn(as.character(METRIC_ORDER), 2, simplify = FALSE), function(pair) {
  a <- metric_overlap[metric == pair[1], genes[[1]]]
  b <- metric_overlap[metric == pair[2], genes[[1]]]
  data.table(
    metric_1 = pair[1],
    metric_2 = pair[2],
    n_metric_1 = length(a),
    n_metric_2 = length(b),
    n_overlap = length(intersect(a, b)),
    overlap_genes = paste(sort(intersect(a, b)), collapse = ";")
  )
}))
fwrite(metric_pairwise, file.path(OUT_DIR, "pairwise_metric_overlaps_all_contrasts.csv"))

contrast_overlap <- hits[, .(genes = list(unique(gene_id_clean))), by = contrast]
contrast_pairwise <- rbindlist(lapply(combn(as.character(CONTRAST_ORDER), 2, simplify = FALSE), function(pair) {
  a <- contrast_overlap[contrast == pair[1], genes[[1]]]
  b <- contrast_overlap[contrast == pair[2], genes[[1]]]
  data.table(
    contrast_1 = pair[1],
    contrast_2 = pair[2],
    n_contrast_1 = length(a),
    n_contrast_2 = length(b),
    n_overlap = length(intersect(a, b)),
    overlap_genes = paste(sort(intersect(a, b)), collapse = ";")
  )
}))
fwrite(contrast_pairwise, file.path(OUT_DIR, "pairwise_contrast_overlaps_all_metrics.csv"))

wb <- createWorkbook()
addWorksheet(wb, "Top recurrent genes")
writeDataTable(wb, "Top recurrent genes", gene_summary)
addWorksheet(wb, "Multiple analyses")
writeDataTable(wb, "Multiple analyses", recurrent_genes)
addWorksheet(wb, "Multiple metrics")
writeDataTable(wb, "Multiple metrics", multi_metric_genes)
addWorksheet(wb, "Multiple contrasts")
writeDataTable(wb, "Multiple contrasts", multi_contrast_genes)
addWorksheet(wb, "All hits long")
writeDataTable(wb, "All hits long", hits)
addWorksheet(wb, "Analysis counts")
writeDataTable(wb, "Analysis counts", analysis_counts)
addWorksheet(wb, "Metric overlaps")
writeDataTable(wb, "Metric overlaps", metric_pairwise)
addWorksheet(wb, "Contrast overlaps")
writeDataTable(wb, "Contrast overlaps", contrast_pairwise)
for (sheet in names(wb)) {
  freezePane(wb, sheet, firstRow = TRUE)
}
saveWorkbook(wb, file.path(OUT_DIR, "limma_rawP_recurrent_gene_summary.xlsx"), overwrite = TRUE)

top_n <- 35
top_genes <- gene_summary[1:min(.N, top_n), gene_id_clean]
plot_dt <- hits[gene_id_clean %in% top_genes]
plot_dt <- merge(plot_dt, gene_summary[, .(gene_id_clean, gene_display = fifelse(gene_name == "" | is.na(gene_name), gene_id_clean, gene_name), n_analyses)], by = "gene_id_clean")
plot_dt[, gene_display := factor(gene_display, levels = rev(gene_summary[gene_id_clean %in% top_genes][order(-n_analyses, best_raw_p), gene_name]))]

heat <- ggplot(plot_dt, aes(x = analysis_label, y = gene_display, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2B6CB0", mid = "white", high = "#B8323B", midpoint = 0, name = "logFC") +
  labs(
    title = "Genes Recurring Across Multiple limma Raw-P Analyses",
    subtitle = "Tiles show significant calls at raw P < 0.05 and |logFC| >= 0.7; color shows direction and effect size",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey30")
  )
ggsave(file.path(OUT_DIR, "top_recurrent_genes_analysis_heatmap.png"), heat, width = 13.5, height = 9.5, dpi = 300)
ggsave(file.path(OUT_DIR, "top_recurrent_genes_analysis_heatmap.pdf"), heat, width = 13.5, height = 9.5)

bar_dt <- gene_summary[1:min(.N, 30)]
bar_dt[, gene_display := fifelse(gene_name == "" | is.na(gene_name), gene_id_clean, gene_name)]
bar_dt[, gene_display := factor(gene_display, levels = rev(gene_display))]
bar <- ggplot(bar_dt, aes(x = n_analyses, y = gene_display, fill = n_metrics)) +
  geom_col(width = 0.72, alpha = 0.94) +
  geom_text(aes(label = paste0(n_analyses, " analyses / ", n_metrics, " metrics")), hjust = -0.05, size = 3.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)), breaks = seq(0, 16, 2)) +
  scale_fill_gradient(low = "#A7B0BA", high = "#B8323B", name = "Metrics") +
  labs(
    title = "Top Recurrent Genes Across limma Raw-P Analyses",
    subtitle = "Ranked by number of metric-contrast analyses where each gene is significant",
    x = "Number of analyses with significant call",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15))
ggsave(file.path(OUT_DIR, "top_recurrent_genes_barplot.png"), bar, width = 10.5, height = 8.5, dpi = 300)
ggsave(file.path(OUT_DIR, "top_recurrent_genes_barplot.pdf"), bar, width = 10.5, height = 8.5)

dist_dt <- gene_summary[, .N, by = n_analyses][order(n_analyses)]
dist <- ggplot(dist_dt, aes(x = factor(n_analyses), y = N)) +
  geom_col(fill = "#2C7A7B", width = 0.72, alpha = 0.95) +
  geom_text(aes(label = N), vjust = -0.25, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "How Often Significant Genes Recur",
    subtitle = "Across 16 raw-p limma analyses: 4 metrics x 4 contrasts",
    x = "Number of analyses where gene is significant",
    y = "Number of genes"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15))
ggsave(file.path(OUT_DIR, "gene_recurrence_frequency_distribution.png"), dist, width = 8, height = 5.5, dpi = 300)
ggsave(file.path(OUT_DIR, "gene_recurrence_frequency_distribution.pdf"), dist, width = 8, height = 5.5)

metric_counts <- hits[, .(gene_count = uniqueN(gene_id_clean)), by = .(metric, metric_label)]
metric_counts[, metric_label := factor(metric_label, levels = METRIC_LABELS[METRIC_ORDER])]
metric_bar <- ggplot(metric_counts, aes(x = metric_label, y = gene_count, fill = metric_label)) +
  geom_col(width = 0.68, alpha = 0.95) +
  geom_text(aes(label = gene_count), vjust = -0.35, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  scale_fill_manual(values = c(
    "Ribosome efficiency" = "#2C7A7B",
    "Protein output" = "#805AD5",
    "Collision" = "#B8323B",
    "Scanning" = "#D69E2E"
  ), guide = "none") +
  labs(
    title = "Unique Significant Genes by Metric",
    subtitle = "A gene is counted once per metric if significant in any contrast",
    x = NULL,
    y = "Unique genes"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15))
ggsave(file.path(OUT_DIR, "unique_significant_genes_by_metric.png"), metric_bar, width = 8.5, height = 5.8, dpi = 300)
ggsave(file.path(OUT_DIR, "unique_significant_genes_by_metric.pdf"), metric_bar, width = 8.5, height = 5.8)

cat("\nRecurrent gene analysis complete\n")
cat("================================\n")
cat("Total unique significant genes:", uniqueN(hits$gene_id_clean), "\n")
cat("Genes significant in >=2 analyses:", nrow(recurrent_genes), "\n")
cat("Genes significant in >=2 metrics:", nrow(multi_metric_genes), "\n")
cat("Genes significant in >=2 contrasts:", nrow(multi_contrast_genes), "\n")
cat("Outputs written to:", OUT_DIR, "\n")
