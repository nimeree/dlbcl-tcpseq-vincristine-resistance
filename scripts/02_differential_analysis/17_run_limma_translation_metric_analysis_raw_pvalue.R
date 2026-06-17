# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# limma differential analysis for continuous translation metrics
# - 2x2 design: Resistant/Sensitive x DMSO/VCR
# - Metrics are analyzed separately
# - Cutoffs kept consistent with DESeq2 analysis: |logFC| >= 0.7
# - Exploratory reporting threshold: raw P.Value < 0.05
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(openxlsx)
})

BASE_DIR <- analysis_path()
INFILE <- file.path(BASE_DIR, "Translation_indexes_fixed", "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")
OUT_DIR <- file.path(BASE_DIR, "Limma_translation_metrics_lfc0.7_rawP0.05")
MATRIX_DIR <- file.path(OUT_DIR, "Matrices")
RESULT_DIR <- file.path(OUT_DIR, "Results")
PLOT_DIR <- file.path(OUT_DIR, "Plots")
INTEGRATION_DIR <- file.path(OUT_DIR, "Multi_metric_integration")
dir.create(MATRIX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INTEGRATION_DIR, recursive = TRUE, showWarnings = FALSE)

P_CUT <- 0.05
LFC_CUT <- 0.7

METRICS <- c(
  "ribosome_efficiency_score",
  "protein_output_score",
  "collision_score",
  "scanning_score"
)

METRIC_LABELS <- c(
  ribosome_efficiency_score = "Ribosome efficiency score",
  protein_output_score = "Protein output score",
  collision_score = "Collision score",
  scanning_score = "Scanning score"
)

CONTRAST_LABELS <- c(
  VCR_sensitive = "VCR effect in sensitive cells",
  VCR_resistant = "VCR effect in resistant cells",
  Resistance_baseline = "Baseline resistance",
  Interaction = "Differential VCR response in resistant vs sensitive"
)

mode_first <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

parse_sample_meta <- function(samples) {
  meta <- data.table(sample = samples)
  meta[, cell_line := fifelse(grepl("^SU8R", sample, ignore.case = TRUE), "Resistant", "Sensitive")]
  meta[, treatment := fifelse(grepl("-Vin-", sample, ignore.case = TRUE), "VCR", "DMSO")]
  meta[, replicate := fifelse(grepl("Rep1", sample, ignore.case = TRUE), "Rep1",
    fifelse(grepl("Rep2", sample, ignore.case = TRUE), "Rep2", NA_character_))]
  meta[, group := factor(
    paste(cell_line, treatment, sep = "_"),
    levels = c("Resistant_DMSO", "Resistant_VCR", "Sensitive_DMSO", "Sensitive_VCR")
  )]
  meta[order(group, replicate)]
}

make_volcano <- function(res, metric, contrast) {
  plot_dt <- copy(res)
  plot_dt[, neglog10_p := -log10(pmax(P.Value, .Machine$double.xmin))]
  plot_dt[, call := fifelse(P.Value < P_CUT & logFC >= LFC_CUT, "Up",
    fifelse(P.Value < P_CUT & logFC <= -LFC_CUT, "Down", "NS"))]
  n_up <- plot_dt[call == "Up", .N]
  n_down <- plot_dt[call == "Down", .N]

  top_lab <- plot_dt[call != "NS"][order(P.Value)][1:min(.N, 12)]
  top_lab[, label := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id_clean)]

  p <- ggplot(plot_dt, aes(x = logFC, y = neglog10_p)) +
    geom_point(aes(color = call), alpha = 0.75, size = 1.25) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dotted", linewidth = 0.55) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 0.55) +
    ggrepel::geom_text_repel(
      data = top_lab,
      aes(label = label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2
    ) +
    annotate(
      "label",
      x = -Inf,
      y = Inf,
      label = paste0("Down: ", n_down),
      hjust = -0.1,
      vjust = 1.2,
      size = 3.2,
      color = "#2B6CB0",
      fill = "white",
      label.size = 0.25
    ) +
    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = paste0("Up: ", n_up),
      hjust = 1.1,
      vjust = 1.2,
      size = 3.2,
      color = "#B8323B",
      fill = "white",
      label.size = 0.25
    ) +
    scale_color_manual(values = c(Up = "#B8323B", Down = "#2B6CB0", NS = "grey72")) +
    labs(
      title = paste(METRIC_LABELS[[metric]], "-", CONTRAST_LABELS[[contrast]]),
      subtitle = paste0("limma eBayes(trend=TRUE); significant: raw P < ", P_CUT, " and |logFC| >= ", LFC_CUT),
      x = "limma logFC",
      y = "-log10(P value)",
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  out_base <- file.path(PLOT_DIR, metric, paste0(contrast, "_volcano_limma_lfc0.7_rawP0.05"))
  dir.create(dirname(out_base), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(out_base, ".png"), p, width = 7.4, height = 5.8, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 7.4, height = 5.8)
}

run_metric_limma <- function(dt, metric) {
  message("[Metric] ", metric)

  metric_dt <- dt[, .(
    gene_name = mode_first(gene_name),
    value = median(get(metric), na.rm = TRUE),
    baseline_cpm_line = median(baseline_cpm_line, na.rm = TRUE)
  ), by = .(gene_id_clean, sample)]

  metric_dt <- metric_dt[baseline_cpm_line > 0]
  if (metric == "scanning_score") {
    metric_dt[value > 10, value := NA_real_]
  }

  mat_dt <- dcast(metric_dt, gene_id_clean ~ sample, value.var = "value")
  gene_annot <- metric_dt[, .(gene_name = mode_first(gene_name)), by = gene_id_clean]

  sample_cols <- setdiff(names(mat_dt), "gene_id_clean")
  keep <- complete.cases(mat_dt[, ..sample_cols])
  mat_dt <- mat_dt[keep]

  mat <- as.matrix(mat_dt[, ..sample_cols])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- mat_dt$gene_id_clean

  meta <- parse_sample_meta(colnames(mat))
  mat <- mat[, meta$sample, drop = FALSE]

  design <- model.matrix(~ 0 + group, data = meta)
  colnames(design) <- levels(meta$group)

  contrasts_mat <- makeContrasts(
    VCR_sensitive = Sensitive_VCR - Sensitive_DMSO,
    VCR_resistant = Resistant_VCR - Resistant_DMSO,
    Resistance_baseline = Resistant_DMSO - Sensitive_DMSO,
    Interaction = (Resistant_VCR - Resistant_DMSO) - (Sensitive_VCR - Sensitive_DMSO),
    levels = design
  )

  fit <- lmFit(mat, design)
  fit2 <- contrasts.fit(fit, contrasts_mat)
  fit2 <- eBayes(fit2, trend = TRUE)

  matrix_out <- merge(
    data.table(gene_id_clean = rownames(mat), mat, check.names = FALSE),
    gene_annot,
    by = "gene_id_clean",
    all.x = TRUE
  )
  setcolorder(matrix_out, c("gene_id_clean", "gene_name", colnames(mat)))
  fwrite(matrix_out, file.path(MATRIX_DIR, paste0("limma_input_matrix_", metric, "_complete_8_samples.csv")))

  metric_result_dir <- file.path(RESULT_DIR, metric)
  dir.create(metric_result_dir, recursive = TRUE, showWarnings = FALSE)

  results <- list()
  direction_hits <- list()
  wb <- createWorkbook()
  for (contrast in colnames(contrasts_mat)) {
    res <- as.data.table(topTable(fit2, coef = contrast, number = Inf, adjust.method = "BH", sort.by = "P"), keep.rownames = "gene_id_clean")
    res <- merge(res, gene_annot, by = "gene_id_clean", all.x = TRUE)
    res[, rank_score := -log10(pmax(P.Value, .Machine$double.xmin)) * abs(logFC)]
    res[, significant_rawP0.05_lfc0.7 := P.Value < P_CUT & abs(logFC) >= LFC_CUT]
    res[, direction := fifelse(significant_rawP0.05_lfc0.7 & logFC >= LFC_CUT, "Up",
      fifelse(significant_rawP0.05_lfc0.7 & logFC <= -LFC_CUT, "Down", "NS"))]
    setcolorder(res, c(
      "gene_id_clean", "gene_name", "logFC", "AveExpr", "t", "P.Value",
      "adj.P.Val", "B", "rank_score", "significant_rawP0.05_lfc0.7", "direction"
    ))
    fwrite(res, file.path(metric_result_dir, paste0(contrast, "_limma_all_genes.csv")))
    sig_res <- res[significant_rawP0.05_lfc0.7 == TRUE]
    sig_up <- res[direction == "Up"]
    sig_down <- res[direction == "Down"]
    fwrite(sig_res, file.path(metric_result_dir, paste0(contrast, "_limma_sig_rawP0.05_lfc0.7.csv")))
    fwrite(sig_up, file.path(metric_result_dir, paste0(contrast, "_limma_sig_up_rawP0.05_lfc0.7.csv")))
    fwrite(sig_down, file.path(metric_result_dir, paste0(contrast, "_limma_sig_down_rawP0.05_lfc0.7.csv")))
    direction_hits[[contrast]] <- rbindlist(list(
      data.table(metric = metric, contrast = contrast, sig_up),
      data.table(metric = metric, contrast = contrast, sig_down)
    ), fill = TRUE)

    addWorksheet(wb, substr(contrast, 1, 31))
    writeDataTable(wb, substr(contrast, 1, 31), res)
    freezePane(wb, substr(contrast, 1, 31), firstRow = TRUE)
    setColWidths(wb, substr(contrast, 1, 31), cols = 1:ncol(res), widths = "auto")

    make_volcano(res, metric, contrast)
    results[[contrast]] <- res
  }

  metric_direction_hits <- rbindlist(direction_hits, fill = TRUE)
  fwrite(
    metric_direction_hits[direction == "Up"],
    file.path(metric_result_dir, paste0(metric, "_limma_sig_up_rawP0.05_lfc0.7_all_contrasts.csv"))
  )
  fwrite(
    metric_direction_hits[direction == "Down"],
    file.path(metric_result_dir, paste0(metric, "_limma_sig_down_rawP0.05_lfc0.7_all_contrasts.csv"))
  )
  fwrite(
    metric_direction_hits[direction == "Up"],
    file.path(MATRIX_DIR, paste0("limma_sig_up_", metric, "_rawP0.05_lfc0.7_all_contrasts.csv"))
  )
  fwrite(
    metric_direction_hits[direction == "Down"],
    file.path(MATRIX_DIR, paste0("limma_sig_down_", metric, "_rawP0.05_lfc0.7_all_contrasts.csv"))
  )

  addWorksheet(wb, "Methods")
  writeData(wb, "Methods", data.table(
    item = c("model", "variance_moderation", "raw_pvalue_cutoff", "logfc_cutoff", "genes_in_matrix"),
    value = c("limma model.matrix(~ 0 + group)", "eBayes(trend=TRUE)", P_CUT, LFC_CUT, nrow(mat))
  ))
  saveWorkbook(wb, file.path(metric_result_dir, paste0(metric, "_limma_results_lfc0.7_rawP0.05.xlsx")), overwrite = TRUE)

  list(
    metric = metric,
    n_complete_genes = nrow(mat),
    sample_metadata = meta,
    design = design,
    contrasts = contrasts_mat,
    results = results
  )
}

message("[Read] ", INFILE)
cols <- unique(c("gene_id_clean", "gene_name", "sample", "baseline_cpm_line", METRICS))
dt <- fread(INFILE, select = cols)
dt <- dt[!is.na(gene_id_clean) & gene_id_clean != ""]
for (cc in c("baseline_cpm_line", METRICS)) {
  dt[, (cc) := as.numeric(get(cc))]
}

analysis <- lapply(METRICS, function(metric) run_metric_limma(dt, metric))
names(analysis) <- METRICS

summary_dt <- rbindlist(lapply(analysis, function(x) {
  rbindlist(lapply(names(x$results), function(contrast) {
    res <- x$results[[contrast]]
    data.table(
      metric = x$metric,
      metric_label = METRIC_LABELS[[x$metric]],
      contrast = contrast,
      contrast_label = CONTRAST_LABELS[[contrast]],
      complete_genes_tested = x$n_complete_genes,
      significant_rawP0.05_lfc0.7 = res[significant_rawP0.05_lfc0.7 == TRUE, .N],
      up = res[direction == "Up", .N],
      down = res[direction == "Down", .N],
      min_p = min(res$P.Value, na.rm = TRUE),
      min_adj_p = min(res$adj.P.Val, na.rm = TRUE)
    )
  }))
}))
fwrite(summary_dt, file.path(OUT_DIR, "limma_metric_contrast_summary_rawP0.05.csv"))

interaction_hits <- rbindlist(lapply(METRICS, function(metric) {
  res <- analysis[[metric]]$results[["Interaction"]]
  res[significant_rawP0.05_lfc0.7 == TRUE, .(
    gene_id_clean,
    gene_name,
    metric,
    logFC,
    P.Value,
    adj.P.Val,
    rank_score,
    direction
  )]
}), fill = TRUE)
fwrite(interaction_hits, file.path(INTEGRATION_DIR, "interaction_significant_genes_by_metric_long.csv"))

if (nrow(interaction_hits)) {
  membership <- interaction_hits[, .(
    gene_name = mode_first(gene_name),
    n_metrics_hit = uniqueN(metric),
    metrics_hit = paste(sort(unique(metric)), collapse = "; ")
  ), by = gene_id_clean][order(-n_metrics_hit, gene_name)]

  wide_logfc <- dcast(interaction_hits, gene_id_clean ~ metric, value.var = "logFC")
  setnames(wide_logfc, setdiff(names(wide_logfc), "gene_id_clean"), paste0("Interaction_logFC_", setdiff(names(wide_logfc), "gene_id_clean")))
  wide_p <- dcast(interaction_hits, gene_id_clean ~ metric, value.var = "P.Value")
  setnames(wide_p, setdiff(names(wide_p), "gene_id_clean"), paste0("Interaction_rawP_", setdiff(names(wide_p), "gene_id_clean")))
  wide_fdr <- dcast(interaction_hits, gene_id_clean ~ metric, value.var = "adj.P.Val")
  setnames(wide_fdr, setdiff(names(wide_fdr), "gene_id_clean"), paste0("Interaction_FDR_", setdiff(names(wide_fdr), "gene_id_clean")))
  membership <- Reduce(function(x, y) merge(x, y, by = "gene_id_clean", all.x = TRUE), list(membership, wide_logfc, wide_p, wide_fdr))
  fwrite(membership, file.path(INTEGRATION_DIR, "interaction_multi_metric_gene_membership.csv"))

  pairwise_overlap <- rbindlist(lapply(combn(METRICS, 2, simplify = FALSE), function(pair) {
    a <- interaction_hits[metric == pair[1], unique(gene_id_clean)]
    b <- interaction_hits[metric == pair[2], unique(gene_id_clean)]
    data.table(
      metric_1 = pair[1],
      metric_2 = pair[2],
      n_metric_1 = length(a),
      n_metric_2 = length(b),
      n_overlap = length(intersect(a, b)),
      overlap_genes = paste(sort(intersect(a, b)), collapse = ";")
    )
  }))
  fwrite(pairwise_overlap, file.path(INTEGRATION_DIR, "interaction_pairwise_metric_overlaps.csv"))
} else {
  fwrite(data.table(), file.path(INTEGRATION_DIR, "interaction_multi_metric_gene_membership.csv"))
  fwrite(data.table(), file.path(INTEGRATION_DIR, "interaction_pairwise_metric_overlaps.csv"))
}

method_note <- c(
  "Continuous translation metrics were analyzed separately with limma.",
  "For each metric, transcript rows were collapsed to gene-sample medians.",
  "Genes with baseline_cpm_line <= 0 were removed before testing.",
  "Only genes with complete values in all eight samples were retained for each metric.",
  "The model used a no-intercept 2x2 group design: Resistant_DMSO, Resistant_VCR, Sensitive_DMSO, Sensitive_VCR.",
  "Contrasts were VCR_sensitive, VCR_resistant, Resistance_baseline, and Interaction.",
  "Empirical Bayes variance moderation used eBayes(trend=TRUE).",
  "Significant genes are reported at raw P.Value < 0.05 and |logFC| >= 0.7.",
  "Because there are n=2 biological replicates per condition, results should be treated as hypothesis-generating."
)
writeLines(method_note, file.path(OUT_DIR, "limma_methods_note.txt"))

message("Done. limma raw-p outputs written to: ", OUT_DIR)
