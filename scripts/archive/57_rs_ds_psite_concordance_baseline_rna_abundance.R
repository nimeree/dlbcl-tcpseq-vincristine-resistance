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
psite_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
metric_file <- file.path(base_dir, "Translation_indexes_fixed", "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")
out_dir <- file.path(psite_dir, "RS_DS_psite_concordance_baseline_RNA_abundance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_cut <- 0.05
lfc_cut <- 0.7

contrast_map <- data.table(
  contrast = c("Resistance_baseline", "VCR_sensitive", "VCR_resistant", "Interaction"),
  psite_contrast = c("Resistance_baseline", "Sensitive_Vin_vs_DMSO", "Resistant_Vin_vs_DMSO", "Interaction")
)

category_levels <- c(
  "Both sig same direction",
  "Both sig opposite direction",
  "DS only",
  "RS only",
  "Neither"
)

direction_call <- function(logfc, pval) {
  fifelse(
    !is.na(pval) & pval < p_cut & !is.na(logfc) & logfc >= lfc_cut, "Up",
    fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc <= -lfc_cut, "Down", "NS")
  )
}

read_fraction_contrast <- function(frac, psite_contrast) {
  f <- file.path(
    psite_dir,
    paste0("Fraction_", frac),
    paste0(psite_contrast, "_", frac, "_psite_limma_all_genes.csv")
  )
  if (!file.exists(f)) stop("Missing P-site limma all-genes file: ", f)

  d <- fread(f)
  if (!"gene_name" %in% names(d)) d[, gene_name := gene_id_clean]
  d[, gene_key := fifelse(!is.na(gene_id_clean) & nzchar(gene_id_clean), gene_id_clean, gene_name)]
  d[, gene_key := sub("\\.\\d+$", "", gene_key)]
  d[, direction_call := direction_call(logFC, P.Value)]
  d[, significant_call := direction_call %in% c("Up", "Down")]

  d[, .(
    gene_key,
    gene_id_clean,
    gene_name,
    logFC,
    AveExpr,
    P.Value,
    adj.P.Val,
    direction = direction_call,
    significant = significant_call
  )]
}

classify_one_contrast <- function(contrast, psite_contrast) {
  rs <- read_fraction_contrast("RS", psite_contrast)
  ds <- read_fraction_contrast("DS", psite_contrast)

  setnames(
    rs,
    c("gene_id_clean", "gene_name", "logFC", "AveExpr", "P.Value", "adj.P.Val", "direction", "significant"),
    c("RS_gene_id_clean", "RS_gene_name", "RS_logFC", "RS_AveExpr", "RS_P.Value", "RS_adj.P.Val", "RS_direction", "RS_significant")
  )
  setnames(
    ds,
    c("gene_id_clean", "gene_name", "logFC", "AveExpr", "P.Value", "adj.P.Val", "direction", "significant"),
    c("DS_gene_id_clean", "DS_gene_name", "DS_logFC", "DS_AveExpr", "DS_P.Value", "DS_adj.P.Val", "DS_direction", "DS_significant")
  )

  joined <- merge(rs, ds, by = "gene_key", all = TRUE)
  joined[is.na(RS_significant), RS_significant := FALSE]
  joined[is.na(DS_significant), DS_significant := FALSE]

  joined[, gene_name := fcoalesce(RS_gene_name, DS_gene_name)]
  joined[, gene_id_clean := fcoalesce(RS_gene_id_clean, DS_gene_id_clean, gene_key)]
  joined[, concordance_category := fifelse(
    RS_significant & DS_significant & RS_direction == DS_direction, "Both sig same direction",
    fifelse(
      RS_significant & DS_significant & RS_direction != DS_direction, "Both sig opposite direction",
      fifelse(
        DS_significant & !RS_significant, "DS only",
        fifelse(RS_significant & !DS_significant, "RS only", "Neither")
      )
    )
  )]
  joined[, concordance_category := factor(concordance_category, levels = category_levels)]
  joined[, `:=`(contrast = contrast, psite_contrast = psite_contrast)]

  setcolorder(joined, c(
    "contrast", "psite_contrast", "gene_key", "gene_id_clean", "gene_name",
    "concordance_category",
    "RS_logFC", "RS_AveExpr", "RS_P.Value", "RS_adj.P.Val", "RS_direction", "RS_significant",
    "DS_logFC", "DS_AveExpr", "DS_P.Value", "DS_adj.P.Val", "DS_direction", "DS_significant"
  ))
  joined[order(concordance_category, gene_name)]
}

message("Classifying RS/DS P-site concordance across contrasts...")
all_concordance <- rbindlist(
  Map(classify_one_contrast, contrast_map$contrast, contrast_map$psite_contrast),
  use.names = TRUE,
  fill = TRUE
)

summary_dt <- all_concordance[, .N, by = .(contrast, psite_contrast, concordance_category)]
summary_dt[, concordance_category := factor(concordance_category, levels = category_levels)]
setorder(summary_dt, contrast, concordance_category)

fwrite(all_concordance, file.path(out_dir, "RS_DS_psite_limma_concordance_all_contrasts.csv"))
fwrite(summary_dt, file.path(out_dir, "RS_DS_psite_limma_concordance_category_counts.csv"))

message("Building baseline RNA abundance table from RNA baseline CPM columns...")
metric <- fread(
  metric_file,
  select = c("gene_id_clean", "gene_name", "transcript", "baseline_sensitive_cpm", "baseline_resistant_cpm")
)
metric[, gene_key := fifelse(!is.na(gene_id_clean) & nzchar(gene_id_clean), gene_id_clean, gene_name)]
metric[, gene_key := sub("\\.\\d+$", "", gene_key)]

baseline_unique <- unique(metric[
  !is.na(gene_key) & nzchar(gene_key),
  .(gene_key, gene_id_clean, gene_name, transcript, baseline_sensitive_cpm, baseline_resistant_cpm)
])
baseline_unique[, baseline_mean_cpm := rowMeans(
  cbind(baseline_sensitive_cpm, baseline_resistant_cpm),
  na.rm = TRUE
)]
baseline_rna <- baseline_unique[
  is.finite(baseline_mean_cpm),
  .(
    baseline_RNA_CPM = median(baseline_mean_cpm, na.rm = TRUE),
    baseline_sensitive_CPM = median(baseline_sensitive_cpm, na.rm = TRUE),
    baseline_resistant_CPM = median(baseline_resistant_cpm, na.rm = TRUE)
  ),
  by = .(gene_key)
]
baseline_rna[, baseline_log2CPM := log2(baseline_RNA_CPM + 1)]

plot_dt <- merge(
  all_concordance[contrast == "Resistance_baseline"],
  baseline_rna,
  by = "gene_key",
  all.x = TRUE
)
plot_dt <- plot_dt[is.finite(baseline_log2CPM)]
plot_dt[, concordance_category := factor(as.character(concordance_category), levels = category_levels)]

baseline_median <- median(plot_dt$baseline_log2CPM, na.rm = TRUE)
category_stats <- plot_dt[, .(
  n_genes = .N,
  median_baseline_log2CPM = median(baseline_log2CPM, na.rm = TRUE),
  mean_baseline_log2CPM = mean(baseline_log2CPM, na.rm = TRUE)
), by = .(concordance_category)]
setorder(category_stats, concordance_category)

fwrite(
  plot_dt,
  file.path(out_dir, "Resistance_baseline_RS_DS_concordance_with_baseline_RNA_log2CPM.csv")
)
fwrite(
  category_stats,
  file.path(out_dir, "Resistance_baseline_baseline_RNA_log2CPM_by_concordance_summary.csv")
)

category_cols <- c(
  "Both sig same direction" = "#2F855A",
  "Both sig opposite direction" = "#C53030",
  "DS only" = "#DD6B20",
  "RS only" = "#805AD5",
  "Neither" = "grey72"
)

p <- ggplot(plot_dt, aes(concordance_category, baseline_log2CPM, fill = concordance_category)) +
  geom_violin(width = 0.86, alpha = 0.42, color = "grey35", linewidth = 0.35, trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.78, linewidth = 0.35) +
  geom_point(
    aes(color = concordance_category),
    position = position_jitter(width = 0.18, height = 0, seed = 7),
    alpha = 0.38,
    size = 0.75
  ) +
  geom_hline(yintercept = baseline_median, linetype = "dashed", color = "black", linewidth = 0.45) +
  scale_fill_manual(values = category_cols, breaks = names(category_cols), drop = FALSE) +
  scale_color_manual(values = category_cols, breaks = names(category_cols), drop = FALSE) +
  labs(
    title = "Resistance baseline RS/DS P-site concordance by baseline RNA abundance",
    subtitle = paste0(
      "Concordance: raw P < ", p_cut, " and |logFC| >= ", lfc_cut,
      "; dashed line = median of plotted genes"
    ),
    x = "RS/DS concordance category",
    y = "Baseline RNA abundance (log2 CPM + 1)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1),
    panel.grid.minor.x = element_blank()
  )

plot_png <- file.path(out_dir, "Resistance_baseline_RS_DS_concordance_baseline_RNA_log2CPM_violin.png")
plot_pdf <- file.path(out_dir, "Resistance_baseline_RS_DS_concordance_baseline_RNA_log2CPM_violin.pdf")
ggsave(plot_png, p, width = 8.5, height = 5.6, dpi = 300, bg = "white")
ggsave(plot_pdf, p, width = 8.5, height = 5.6, bg = "white")

cat("\nRS/DS P-site concordance classification complete.\n")
cat("Baseline RNA source: baseline_sensitive_cpm and baseline_resistant_cpm from:\n", metric_file, "\n", sep = "")
cat("Plotted abundance: log2(mean baseline RNA CPM + 1), collapsed to gene-level median across transcripts.\n")
cat("\nCategory counts:\n")
print(summary_dt)
cat("\nResistance_baseline baseline RNA abundance summary:\n")
print(category_stats)
cat("\nAll-gene plotted median log2CPM: ", round(baseline_median, 4), "\n", sep = "")
cat("\nSaved outputs to:\n", out_dir, "\n", sep = "")
