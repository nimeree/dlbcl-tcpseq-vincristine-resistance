# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# GO Biological Process-only g:Profiler ORA for baseline resistance
# Combined Up + Down significant genes.
#
# P-site fraction limma baseline is generated here as:
#   Resistant_DMSO - Sensitive_DMSO
#
# Translation metric limma baseline uses existing:
#   Resistance_baseline_limma_all_genes.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(gprofiler2)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(openxlsx)
})

BASE_DIR <- analysis_path()
PSITE_DIR <- file.path(BASE_DIR, "Psite_fraction_limma_lfc0.7_rawP0.05")
METRIC_DIR <- file.path(BASE_DIR, "Limma_translation_metrics_lfc0.7_rawP0.05")

PSITE_OUT <- file.path(PSITE_DIR, "GO_BP_ORA_Baseline_Combined")
METRIC_OUT <- file.path(METRIC_DIR, "GO_BP_ORA_Baseline_Combined")

SOURCE <- "GO:BP"
USER_THRESHOLD <- 0.05
MAX_TERM_SIZE <- 500
MIN_INTERSECTION_SIZE <- 5
TOP_N <- 15
P_CUT <- 0.05
LFC_CUT <- 0.7

nice_names <- c(
  SSU = "SSU",
  RS = "RS",
  DS = "DS",
  scanning_score = "Scanning",
  ribosome_efficiency_score = "Ribosome_engagement",
  collision_score = "Collision",
  protein_output_score = "Protein_output"
)

display_names <- c(
  SSU = "SSU P-site count baseline",
  RS = "RS P-site count baseline",
  DS = "DS P-site count baseline",
  scanning_score = "Scanning score",
  ribosome_efficiency_score = "Ribosome engagement score",
  collision_score = "Collision score",
  protein_output_score = "Protein output score"
)

split_intersection <- function(x) {
  if (is.na(x) || x == "") character(0) else unlist(strsplit(x, ",", fixed = TRUE))
}

direction_call <- function(logfc, pval) {
  fifelse(!is.na(pval) & pval < P_CUT & !is.na(logfc) & logfc >= LFC_CUT, "Up",
          fifelse(!is.na(pval) & pval < P_CUT & !is.na(logfc) & logfc <= -LFC_CUT, "Down", "NS"))
}

ensure_psite_baseline_limma <- function() {
  gene_counts_file <- file.path(PSITE_DIR, "psite_gene_counts_long_by_fraction_sample.csv")
  meta_file <- file.path(PSITE_DIR, "psite_fraction_limma_sample_metadata.csv")
  if (!file.exists(gene_counts_file) || !file.exists(meta_file)) {
    stop("Missing P-site limma count inputs in: ", PSITE_DIR)
  }

  gene_counts_long <- fread(gene_counts_file)
  sample_meta <- fread(meta_file)
  sample_meta <- sample_meta[condition %in% c("Sensitive_DMSO", "Resistant_DMSO")]
  sample_meta[, condition := factor(condition, levels = c("Sensitive_DMSO", "Resistant_DMSO"))]

  summary_rows <- list()

  for (frac in c("SSU", "RS", "DS")) {
    frac_dir <- file.path(PSITE_DIR, paste0("Fraction_", frac))
    dir.create(frac_dir, recursive = TRUE, showWarnings = FALSE)
    out_all <- file.path(frac_dir, paste0("Resistance_baseline_", frac, "_psite_limma_all_genes.csv"))
    out_sig <- file.path(frac_dir, paste0("Resistance_baseline_", frac, "_psite_limma_sig_rawP0.05_lfc0.7.csv"))
    out_up <- file.path(frac_dir, paste0("Resistance_baseline_", frac, "_psite_limma_sig_up_rawP0.05_lfc0.7.csv"))
    out_down <- file.path(frac_dir, paste0("Resistance_baseline_", frac, "_psite_limma_sig_down_rawP0.05_lfc0.7.csv"))

    d <- gene_counts_long[fraction == frac & sample %in% sample_meta[fraction == frac, sample]]
    meta <- unique(sample_meta[fraction == frac][order(condition, replicate)])
    sample_levels <- meta$sample

    mat_dt <- dcast(d, gene_id_clean + gene_name ~ sample, value.var = "psite_count", fill = 0)
    gene_info <- mat_dt[, .(gene_id_clean, gene_name)]
    mat <- as.matrix(mat_dt[, ..sample_levels])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- make.unique(gene_info$gene_id_clean)

    keep <- rowSums(mat >= 10) >= 2
    mat <- mat[keep, , drop = FALSE]
    gene_info <- gene_info[keep]
    lib_sizes <- colSums(mat)
    log_cpm <- log2(t(t(mat + 0.5) / pmax(lib_sizes, 1)) * 1e6)

    group <- factor(meta$condition, levels = c("Sensitive_DMSO", "Resistant_DMSO"))
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    fit <- lmFit(log_cpm, design)
    fit <- eBayes(fit, trend = TRUE)
    contrast_matrix <- makeContrasts(Resistance_baseline = Resistant_DMSO - Sensitive_DMSO, levels = design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2, trend = TRUE)

    tt <- as.data.table(topTable(fit2, coef = "Resistance_baseline", number = Inf, sort.by = "none"), keep.rownames = "row_key")
    tt[, gene_id_clean := gene_info$gene_id_clean]
    tt[, gene_name := gene_info$gene_name]
    tt[, significant_rawP0.05_lfc0.7 := !is.na(P.Value) & P.Value < P_CUT & abs(logFC) >= LFC_CUT]
    tt[, direction := direction_call(logFC, P.Value)]
    setcolorder(tt, c("gene_id_clean", "gene_name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "significant_rawP0.05_lfc0.7", "direction"))
    tt <- tt[order(P.Value)]

    fwrite(tt, out_all)
    fwrite(tt[significant_rawP0.05_lfc0.7 == TRUE], out_sig)
    fwrite(tt[direction == "Up"], out_up)
    fwrite(tt[direction == "Down"], out_down)

    summary_rows[[frac]] <- data.table(
      fraction = frac,
      genes_tested = nrow(tt),
      n_sig = tt[significant_rawP0.05_lfc0.7 == TRUE, .N],
      n_up = tt[direction == "Up", .N],
      n_down = tt[direction == "Down", .N],
      median_logFC = median(tt$logFC, na.rm = TRUE)
    )
  }

  summary_dt <- rbindlist(summary_rows)
  fwrite(summary_dt, file.path(PSITE_DIR, "psite_fraction_limma_baseline_summary.csv"))
  summary_dt
}

read_psite_query <- function(frac) {
  f <- file.path(PSITE_DIR, paste0("Fraction_", frac), paste0("Resistance_baseline_", frac, "_psite_limma_sig_rawP0.05_lfc0.7.csv"))
  if (!file.exists(f)) stop("Missing P-site baseline significant table: ", f)
  dt <- fread(f)
  if ("gene_id_clean" %in% names(dt)) setnames(dt, "gene_id_clean", "gene_id")
  if ("logFC" %in% names(dt)) setnames(dt, "logFC", "log2FoldChange")
  if (!"gene_name" %in% names(dt)) dt[, gene_name := NA_character_]
  if (!"direction" %in% names(dt)) dt[, direction := fifelse(log2FoldChange >= 0, "Up", "Down")]
  dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
  unique(dt[!is.na(gene_id) & gene_id != "", .(gene_id, gene_name, log2FoldChange, direction)], by = "gene_id")
}

read_metric_query <- function(metric) {
  f <- file.path(METRIC_DIR, "Results", metric, "Resistance_baseline_limma_all_genes.csv")
  if (!file.exists(f)) stop("Missing metric baseline limma table: ", f)
  dt <- fread(f)
  if ("gene_id_clean" %in% names(dt)) setnames(dt, "gene_id_clean", "gene_id")
  if ("logFC" %in% names(dt)) setnames(dt, "logFC", "log2FoldChange")
  if (!"gene_name" %in% names(dt)) dt[, gene_name := NA_character_]
  dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
  dt[, direction := direction_call(log2FoldChange, P.Value)]
  unique(dt[direction %in% c("Up", "Down") & !is.na(gene_id) & gene_id != "",
            .(gene_id, gene_name, log2FoldChange, direction)], by = "gene_id")
}

add_hit_metadata <- function(ora_dt, query_dt) {
  symbol_map <- unique(query_dt[, .(gene_id, gene_name, log2FoldChange, direction)])
  ora_dt[, hit_gene_ids := intersection]
  ora_dt[, hit_gene_symbols := vapply(intersection, function(x) {
    ids <- split_intersection(x)
    hit <- symbol_map[gene_id %in% ids]
    labels <- ifelse(!is.na(hit$gene_name) & hit$gene_name != "", hit$gene_name, hit$gene_id)
    paste(sort(unique(labels)), collapse = ";")
  }, character(1))]
  ora_dt[, mean_hit_log2FC := vapply(intersection, function(x) {
    ids <- split_intersection(x)
    vals <- symbol_map[gene_id %in% ids]$log2FoldChange
    if (length(vals) == 0) NA_real_ else mean(vals, na.rm = TRUE)
  }, numeric(1))]
  ora_dt[, n_hit_up := vapply(intersection, function(x) {
    ids <- split_intersection(x)
    symbol_map[gene_id %in% ids & direction == "Up", .N]
  }, integer(1))]
  ora_dt[, n_hit_down := vapply(intersection, function(x) {
    ids <- split_intersection(x)
    symbol_map[gene_id %in% ids & direction == "Down", .N]
  }, integer(1))]
  ora_dt[, hit_up_genes := vapply(intersection, function(x) {
    ids <- split_intersection(x)
    hit <- symbol_map[gene_id %in% ids & direction == "Up"]
    labels <- ifelse(!is.na(hit$gene_name) & hit$gene_name != "", hit$gene_name, hit$gene_id)
    paste(sort(unique(labels)), collapse = ";")
  }, character(1))]
  ora_dt[, hit_down_genes := vapply(intersection, function(x) {
    ids <- split_intersection(x)
    hit <- symbol_map[gene_id %in% ids & direction == "Down"]
    labels <- ifelse(!is.na(hit$gene_name) & hit$gene_name != "", hit$gene_name, hit$gene_id)
    paste(sort(unique(labels)), collapse = ";")
  }, character(1))]
  ora_dt[]
}

select_terms_for_plot <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  plot_dt <- dt[order(p_value, -intersection_size)]
  plot_dt <- head(plot_dt, TOP_N)
  plot_dt[, term_label := stringr::str_wrap(term_name, width = 48)]
  plot_dt[, neglog10_gscs := -log10(p_value)]
  plot_dt[]
}

make_plot <- function(dt, title, out_base) {
  plot_dt <- select_terms_for_plot(dt)
  if (nrow(plot_dt) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No GO Biological Process terms after filters", size = 5) +
      labs(title = title) +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave(paste0(out_base, ".png"), p, width = 10, height = 5.5, dpi = 300, bg = "white")
    ggsave(paste0(out_base, ".pdf"), p, width = 10, height = 5.5, bg = "white")
    return(invisible(NULL))
  }

  p <- ggplot(plot_dt, aes(
    x = neglog10_gscs,
    y = forcats::fct_reorder(term_label, neglog10_gscs)
  )) +
    geom_segment(aes(
      x = 0, xend = neglog10_gscs,
      y = forcats::fct_reorder(term_label, neglog10_gscs),
      yend = forcats::fct_reorder(term_label, neglog10_gscs)
    ), linewidth = 0.8, color = "grey78") +
    geom_point(aes(size = intersection_size, color = mean_hit_log2FC), alpha = 0.95) +
    scale_color_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
    scale_size_continuous(range = c(3, 8)) +
    labs(
      title = title,
      subtitle = paste0("Baseline; GO Biological Process only; Up + Down combined; g:SCS < ", USER_THRESHOLD,
                        "; term size <= ", MAX_TERM_SIZE,
                        "; hit genes >= ", MIN_INTERSECTION_SIZE),
      x = "-log10(g:SCS corrected p-value)",
      y = NULL,
      size = "Hit genes",
      color = "Mean hit\nlogFC"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      axis.text.y = element_text(color = "black", size = 9.5),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )

  height <- max(6.3, min(10.5, 2.4 + 0.36 * nrow(plot_dt)))
  ggsave(paste0(out_base, ".png"), p, width = 11.5, height = height, dpi = 300, limitsize = FALSE, bg = "white")
  ggsave(paste0(out_base, ".pdf"), p, width = 11.5, height = height, limitsize = FALSE, bg = "white")
}

run_one_query <- function(dataset, key, query_dt, out_root) {
  nice <- nice_names[[key]]
  label <- display_names[[key]]
  query <- unique(query_dt$gene_id)
  message("[Baseline GO:BP ORA] ", dataset, " / ", nice, " (", length(query), " genes)")

  tables_dir <- file.path(out_root, "Tables")
  plots_dir <- file.path(out_root, "Plots", nice)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(query_dt, file.path(tables_dir, paste0(nice, "_query_genes.csv")))

  if (length(query) == 0) {
    all_dt <- data.table()
  } else {
    res <- gprofiler2::gost(
      query = query,
      organism = "hsapiens",
      sources = SOURCE,
      correction_method = "g_SCS",
      domain_scope = "annotated",
      user_threshold = USER_THRESHOLD,
      evcodes = TRUE
    )
    if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
      all_dt <- data.table()
    } else {
      all_dt <- as.data.table(res$result)
      all_dt <- add_hit_metadata(all_dt, query_dt)
      all_dt[, `:=`(
        dataset = dataset,
        contrast = "Baseline",
        analysis = key,
        analysis_label = label,
        direction = "Combined",
        query_gene_count_input = length(query),
        query_up_genes = query_dt[direction == "Up", .N],
        query_down_genes = query_dt[direction == "Down", .N],
        correction_method = "g_SCS",
        domain_scope = "annotated",
        user_threshold = USER_THRESHOLD,
        max_term_size_filter = MAX_TERM_SIZE,
        min_intersection_filter = MIN_INTERSECTION_SIZE
      )]
    }
    Sys.sleep(0.5)
  }

  filtered_dt <- if (nrow(all_dt) == 0 || !"significant" %in% names(all_dt)) {
    data.table()
  } else {
    all_dt[
      significant == TRUE &
        p_value < USER_THRESHOLD &
        term_size <= MAX_TERM_SIZE &
        intersection_size >= MIN_INTERSECTION_SIZE
    ]
  }
  plotted_dt <- select_terms_for_plot(filtered_dt)

  fwrite(all_dt, file.path(tables_dir, paste0(nice, "_all_gSCS_significant_terms.csv")))
  fwrite(filtered_dt, file.path(tables_dir, paste0(nice, "_filtered_terms.csv")))
  fwrite(plotted_dt, file.path(tables_dir, paste0(nice, "_terms_shown_in_plot.csv")))
  make_plot(filtered_dt, paste(label, "baseline GO:BP ORA"), file.path(plots_dir, paste0(nice, "_GO_BP_ORA_baseline_combined")))

  list(
    all = all_dt,
    filtered = filtered_dt,
    plotted = plotted_dt,
    query = query_dt[, .(dataset, contrast = "Baseline", analysis = key, analysis_label = label, gene_id, gene_name, log2FoldChange, direction)],
    summary = data.table(
      dataset = dataset,
      contrast = "Baseline",
      analysis = key,
      analysis_label = label,
      query_gene_count = length(query),
      query_up_genes = query_dt[direction == "Up", .N],
      query_down_genes = query_dt[direction == "Down", .N],
      gprofiler_significant_terms = nrow(all_dt),
      filtered_terms = nrow(filtered_dt),
      plotted_terms = nrow(plotted_dt)
    )
  )
}

write_dataset_workbook <- function(out_root, results, workbook_name) {
  tables_dir <- file.path(out_root, "Tables")
  all_dt <- rbindlist(lapply(results, `[[`, "all"), fill = TRUE)
  filtered_dt <- rbindlist(lapply(results, `[[`, "filtered"), fill = TRUE)
  plotted_dt <- rbindlist(lapply(results, `[[`, "plotted"), fill = TRUE)
  query_dt <- rbindlist(lapply(results, `[[`, "query"), fill = TRUE)
  summary_dt <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)

  fwrite(summary_dt, file.path(tables_dir, "summary.csv"))
  fwrite(query_dt, file.path(tables_dir, "combined_query_genes.csv"))
  fwrite(all_dt, file.path(tables_dir, "all_gSCS_significant_terms.csv"))
  fwrite(filtered_dt, file.path(tables_dir, "filtered_terms.csv"))
  fwrite(plotted_dt, file.path(tables_dir, "terms_shown_in_plots.csv"))

  wb <- createWorkbook()
  addWorksheet(wb, "Summary")
  writeDataTable(wb, "Summary", summary_dt, tableStyle = "TableStyleMedium2")
  freezePane(wb, "Summary", firstRow = TRUE)
  setColWidths(wb, "Summary", cols = 1:ncol(summary_dt), widths = "auto")

  addWorksheet(wb, "Terms_shown_in_plots")
  writeDataTable(wb, "Terms_shown_in_plots", plotted_dt, tableStyle = "TableStyleMedium5")
  freezePane(wb, "Terms_shown_in_plots", firstRow = TRUE)
  if (ncol(plotted_dt) > 0) setColWidths(wb, "Terms_shown_in_plots", cols = 1:ncol(plotted_dt), widths = "auto")

  addWorksheet(wb, "Filtered_terms")
  writeDataTable(wb, "Filtered_terms", filtered_dt, tableStyle = "TableStyleMedium4")
  freezePane(wb, "Filtered_terms", firstRow = TRUE)
  if (ncol(filtered_dt) > 0) setColWidths(wb, "Filtered_terms", cols = 1:ncol(filtered_dt), widths = "auto")

  addWorksheet(wb, "Query_genes")
  writeDataTable(wb, "Query_genes", query_dt, tableStyle = "TableStyleMedium9")
  freezePane(wb, "Query_genes", firstRow = TRUE)
  if (ncol(query_dt) > 0) setColWidths(wb, "Query_genes", cols = 1:ncol(query_dt), widths = "auto")

  saveWorkbook(wb, file.path(out_root, workbook_name), overwrite = TRUE)
  summary_dt
}

baseline_count_summary <- ensure_psite_baseline_limma()

psite_results <- list(
  SSU = run_one_query("P-site fraction limma", "SSU", read_psite_query("SSU"), PSITE_OUT),
  RS = run_one_query("P-site fraction limma", "RS", read_psite_query("RS"), PSITE_OUT),
  DS = run_one_query("P-site fraction limma", "DS", read_psite_query("DS"), PSITE_OUT)
)
psite_summary <- write_dataset_workbook(PSITE_OUT, psite_results, "GO_BP_ORA_Psite_fraction_limma_baseline_combined.xlsx")

metric_keys <- c("scanning_score", "ribosome_efficiency_score", "collision_score", "protein_output_score")
metric_results <- setNames(lapply(metric_keys, function(metric) {
  run_one_query("Translation metric limma", metric, read_metric_query(metric), METRIC_OUT)
}), metric_keys)
metric_summary <- write_dataset_workbook(METRIC_OUT, metric_results, "GO_BP_ORA_translation_metric_limma_baseline_combined.xlsx")

cat("\nDone. Baseline GO:BP-only combined ORA outputs saved to:\n")
cat("P-site fraction limma: ", PSITE_OUT, "\n", sep = "")
cat("Translation metric limma: ", METRIC_OUT, "\n\n", sep = "")
cat("P-site baseline limma summary:\n")
print(baseline_count_summary)
cat("\nP-site ORA summary:\n")
print(psite_summary)
cat("\nMetric ORA summary:\n")
print(metric_summary)
