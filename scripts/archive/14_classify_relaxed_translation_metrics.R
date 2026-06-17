# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# Relaxed translation-state metrics from t2g_v3 DESeq2 outputs
# - Efficient translation: SSU Up or NS, RS Up, DS Down or NS
# - Scanning bottleneck: SSU Up, RS Down or NS
# - Collision/stalling: SSU Up, RS Down or NS, DS Up
# - Metrics are non-exclusive: a gene may satisfy more than one metric.
# - Threshold: pvalue < 0.05 and abs(log2FoldChange) >= 0.7
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
})

BASE_DIR <- analysis_path()
OUT_DIR <- file.path(BASE_DIR, "Translation_Metrics_Relaxed_lfc0.7_p0.05")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

P_CUT <- 0.05
LFC_CUT <- 0.7
FRACTIONS <- c("SSU", "RS", "DS")

COMPARISONS <- c(
  "Sensitive_Vin_vs_DMSO",
  "Resistant_Vin_vs_DMSO",
  "Vin_Resistant_vs_Sensitive"
)

COMPARISON_TITLES <- c(
  Sensitive_Vin_vs_DMSO = "Sensitive Vin vs DMSO",
  Resistant_Vin_vs_DMSO = "Resistant Vin vs DMSO",
  Vin_Resistant_vs_Sensitive = "Vin Resistant vs Sensitive"
)

METRIC_LEVELS <- c(
  "Efficient translation",
  "Scanning bottleneck",
  "Collision / stalling"
)

METRIC_DEFINITIONS <- data.table(
  translation_metric = METRIC_LEVELS,
  logic = c(
    "SSU Up or NS, RS Up, DS Down or NS",
    "SSU Up, RS Down or NS",
    "SSU Up, RS Down or NS, DS Up"
  ),
  interpretation = c(
    "Increased 80S/translating ribosome signal without increased disome burden.",
    "Increased small-subunit signal without increased 80S/translating ribosome signal.",
    "Increased small-subunit and disome/collision signal without increased 80S/translating ribosome signal."
  )
)

status_call <- function(pvalue, log2fc) {
  fifelse(
    !is.na(pvalue) & !is.na(log2fc) & pvalue < P_CUT & log2fc >= LFC_CUT,
    "Up",
    fifelse(
      !is.na(pvalue) & !is.na(log2fc) & pvalue < P_CUT & log2fc <= -LFC_CUT,
      "Down",
      "NS"
    )
  )
}

first_or_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[1]
}

read_fraction_result <- function(comparison, fraction) {
  f <- file.path(
    BASE_DIR,
    paste0("Fraction_", fraction),
    paste0(comparison, "_", fraction, "_results_all.csv")
  )
  if (!file.exists(f)) stop("Missing DESeq2 file: ", f)

  dt <- fread(f)
  dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
  dt[, status := status_call(pvalue, log2FoldChange)]
  dt[, .(gene_id, log2FoldChange, pvalue, padj, status)]
}

read_annotation <- function() {
  xlsx_files <- list.files(
    BASE_DIR,
    pattern = "_top500_differential_genes\\.xlsx$",
    recursive = TRUE,
    full.names = TRUE
  )
  xlsx_files <- xlsx_files[!grepl("/~\\$", xlsx_files)]
  if (length(xlsx_files) == 0) {
    return(data.table(
      gene_id = character(),
      gene_name = character(),
      gene_function = character(),
      gene_summary = character()
    ))
  }

  ann <- rbindlist(lapply(xlsx_files, function(f) {
    d <- as.data.table(read.xlsx(f))
    needed <- c("gene_id", "gene_name", "gene_function", "gene_summary")
    if (!all(needed %in% names(d))) return(NULL)
    d[, gene_id := sub("\\.\\d+$", "", gene_id)]
    d[, ..needed]
  }), fill = TRUE)

  ann <- ann[!is.na(gene_id) & gene_id != ""]
  ann[, .(
    gene_name = first_or_na(gene_name),
    gene_function = first_or_na(gene_function),
    gene_summary = first_or_na(gene_summary)
  ), by = gene_id]
}

make_metric_rows <- function(dt) {
  rows <- list()

  efficient <- dt[SSU_status %in% c("Up", "NS") & RS_status == "Up" & DS_status %in% c("Down", "NS")]
  if (nrow(efficient) > 0) {
    efficient[, translation_metric := "Efficient translation"]
    rows[["Efficient translation"]] <- efficient
  }

  scanning <- dt[SSU_status == "Up" & RS_status %in% c("Down", "NS")]
  if (nrow(scanning) > 0) {
    scanning[, translation_metric := "Scanning bottleneck"]
    rows[["Scanning bottleneck"]] <- scanning
  }

  collision <- dt[SSU_status == "Up" & RS_status %in% c("Down", "NS") & DS_status == "Up"]
  if (nrow(collision) > 0) {
    collision[, translation_metric := "Collision / stalling"]
    rows[["Collision / stalling"]] <- collision
  }

  rbindlist(rows, fill = TRUE)
}

annotation <- read_annotation()
all_metric_rows <- list()

for (comparison in COMPARISONS) {
  frac_tables <- lapply(FRACTIONS, function(frac) {
    d <- read_fraction_result(comparison, frac)
    setnames(
      d,
      c("log2FoldChange", "pvalue", "padj", "status"),
      paste0(frac, c("_log2FC", "_pvalue", "_padj", "_status"))
    )
    d
  })

  merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), frac_tables)
  for (frac in FRACTIONS) {
    status_col <- paste0(frac, "_status")
    merged[is.na(get(status_col)), (status_col) := "NS"]
  }

  metric_rows <- make_metric_rows(merged)
  if (nrow(metric_rows) == 0) next
  metric_rows[, comparison := comparison]
  metric_rows <- merge(metric_rows, annotation, by = "gene_id", all.x = TRUE)
  metric_rows[, gene_name := fifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)]
  all_metric_rows[[comparison]] <- metric_rows
}

metric_dt <- rbindlist(all_metric_rows, fill = TRUE)
metric_dt[, translation_metric := factor(translation_metric, levels = METRIC_LEVELS)]
metric_dt[, comparison_title := COMPARISON_TITLES[comparison]]
metric_dt[, comparison_title := factor(comparison_title, levels = COMPARISON_TITLES[COMPARISONS])]

count_dt <- metric_dt[, .(gene_count = uniqueN(gene_id)), by = .(
  comparison,
  comparison_title,
  translation_metric
)]
all_grid <- CJ(comparison = COMPARISONS, translation_metric = METRIC_LEVELS)
all_grid[, comparison_title := COMPARISON_TITLES[comparison]]
all_grid[, comparison_title := factor(comparison_title, levels = COMPARISON_TITLES[COMPARISONS])]
all_grid[, translation_metric := factor(translation_metric, levels = METRIC_LEVELS)]
count_dt <- merge(all_grid, count_dt, by = c("comparison", "comparison_title", "translation_metric"), all.x = TRUE)
count_dt[is.na(gene_count), gene_count := 0L]

metric_colors <- c(
  "Efficient translation" = "#2C7A7B",
  "Scanning bottleneck" = "#D69E2E",
  "Collision / stalling" = "#B8323B"
)

p <- ggplot(count_dt, aes(x = translation_metric, y = gene_count, fill = translation_metric)) +
  geom_col(width = 0.68, alpha = 0.94) +
  geom_text(aes(label = gene_count), vjust = -0.35, size = 4.8, fontface = "bold") +
  facet_wrap(~ comparison_title, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = metric_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = "Relaxed Translation Metrics from DESeq2 Fraction Patterns",
    subtitle = paste0("Thresholds: p < ", P_CUT, " and abs(log2FC) >= ", LFC_CUT, "; metrics are non-exclusive"),
    x = NULL,
    y = "Number of genes"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 12, color = "black", angle = 25, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 18, 12, 18)
  )

plot_base <- file.path(OUT_DIR, "relaxed_translation_metric_count_barplot_lfc0.7_p0.05")
ggsave(paste0(plot_base, ".png"), p, width = 12, height = 6.8, dpi = 300)
ggsave(paste0(plot_base, ".pdf"), p, width = 12, height = 6.8)

fwrite(count_dt, file.path(OUT_DIR, "relaxed_translation_metric_counts.csv"))
fwrite(metric_dt, file.path(OUT_DIR, "relaxed_translation_metric_gene_list.csv"))

wb <- createWorkbook()
addWorksheet(wb, "Metric definitions")
writeData(wb, "Metric definitions", METRIC_DEFINITIONS)

addWorksheet(wb, "Counts")
writeData(wb, "Counts", count_dt[order(comparison, translation_metric)])

excel_cols <- c(
  "comparison", "translation_metric", "gene_id", "gene_name",
  "SSU_status", "SSU_log2FC", "SSU_pvalue", "SSU_padj",
  "RS_status", "RS_log2FC", "RS_pvalue", "RS_padj",
  "DS_status", "DS_log2FC", "DS_pvalue", "DS_padj",
  "gene_function", "gene_summary"
)

for (metric in METRIC_LEVELS) {
  sheet_name <- substr(gsub("[/]", "-", metric), 1, 31)
  addWorksheet(wb, sheet_name)
  out <- metric_dt[as.character(translation_metric) == metric]
  out <- out[order(comparison, gene_name)]
  writeData(wb, sheet_name, out[, ..excel_cols])
}

for (comp in COMPARISONS) {
  sheet_name <- substr(gsub("_", " ", comp), 1, 31)
  addWorksheet(wb, sheet_name)
  out <- metric_dt[comparison == comp]
  out <- out[order(translation_metric, gene_name)]
  writeData(wb, sheet_name, out[, ..excel_cols])
}

saveWorkbook(
  wb,
  file.path(OUT_DIR, "relaxed_translation_metric_gene_lists_lfc0.7_p0.05.xlsx"),
  overwrite = TRUE
)

message("Done. Relaxed translation metric outputs saved to: ", OUT_DIR)
