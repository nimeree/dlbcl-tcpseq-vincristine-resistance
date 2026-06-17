# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# Common/uncommon gene workbook for t2g_v3 DESeq2 analyses
# - Uses the same Venn thresholds:
#   pvalue < 0.05 and abs(log2FoldChange) >= 0.7
# - One workbook with:
#   * Summary region counts
#   * One detail sheet per comparison/direction
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

IN_DIR <- analysis_path()
OUT_DIR <- file.path(IN_DIR, "Venn_Diagrams_lfc0.7")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_XLSX <- file.path(OUT_DIR, "DESeq2_common_uncommon_genes_SSU_RS_DS_lfc0.7.xlsx")

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

SHEET_PREFIX <- c(
  Sensitive_Vin_vs_DMSO = "Sens",
  Resistant_Vin_vs_DMSO = "Res",
  Vin_Resistant_vs_Sensitive = "VRS"
)

read_result <- function(comparison, fraction) {
  f <- file.path(IN_DIR, paste0("Fraction_", fraction), paste0(comparison, "_", fraction, "_results_all.csv"))
  if (!file.exists(f)) stop("Missing input file: ", f)
  dt <- fread(f)
  dt[, fraction := fraction]
  dt[, comparison := comparison]
  dt
}

first_non_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[1]
}

make_detail <- function(comparison, direction) {
  per_fraction <- setNames(lapply(FRACTIONS, function(frac) read_result(comparison, frac)), FRACTIONS)

  sig_sets <- lapply(per_fraction, function(dt) {
    dt <- dt[!is.na(pvalue) & !is.na(log2FoldChange)]
    if (direction == "Up") {
      dt <- dt[pvalue < P_CUT & log2FoldChange >= LFC_CUT]
    } else {
      dt <- dt[pvalue < P_CUT & log2FoldChange <= -LFC_CUT]
    }
    unique(dt$gene_id)
  })

  union_genes <- sort(unique(unlist(sig_sets)))
  detail <- data.table(gene_id = union_genes)
  detail[, in_SSU := gene_id %in% sig_sets$SSU]
  detail[, in_RS  := gene_id %in% sig_sets$RS]
  detail[, in_DS  := gene_id %in% sig_sets$DS]
  detail[, n_fractions := in_SSU + in_RS + in_DS]
  detail[, region := fcase(
    in_SSU & in_RS & in_DS, "Common_all_three",
    in_SSU & in_RS, "Shared_SSU_RS",
    in_SSU & in_DS, "Shared_SSU_DS",
    in_RS & in_DS, "Shared_RS_DS",
    in_SSU, "Unique_SSU",
    in_RS, "Unique_RS",
    in_DS, "Unique_DS"
  )]
  detail[, common_status := fifelse(region == "Common_all_three", "Common", "Uncommon")]

  all_rows <- rbindlist(per_fraction, fill = TRUE)
  annot <- all_rows[, .(
    gene_name = first_non_na(gene_name),
    gene_function = first_non_na(gene_function),
    gene_type = first_non_na(gene_type)
  ), by = gene_id]
  detail <- merge(detail, annot, by = "gene_id", all.x = TRUE)

  for (frac in FRACTIONS) {
    cols <- per_fraction[[frac]][, .(
      gene_id,
      log2FoldChange,
      pvalue,
      padj
    )]
    if (direction == "Up") {
      cols[, significant := !is.na(pvalue) & !is.na(log2FoldChange) & pvalue < P_CUT & log2FoldChange >= LFC_CUT]
    } else {
      cols[, significant := !is.na(pvalue) & !is.na(log2FoldChange) & pvalue < P_CUT & log2FoldChange <= -LFC_CUT]
    }
    setnames(cols, c("log2FoldChange", "pvalue", "padj", "significant"),
             paste0(c("log2FC_", "pvalue_", "padj_", "significant_"), frac))
    detail <- merge(detail, cols, by = "gene_id", all.x = TRUE)
  }

  region_order <- c(
    "Common_all_three",
    "Shared_SSU_RS",
    "Shared_SSU_DS",
    "Shared_RS_DS",
    "Unique_SSU",
    "Unique_RS",
    "Unique_DS"
  )
  detail[, region := factor(region, levels = region_order)]
  setorder(detail, region, gene_name, gene_id)
  detail[, region := as.character(region)]

  setcolorder(detail, c(
    "common_status", "region", "n_fractions",
    "gene_id", "gene_name", "gene_function", "gene_type",
    "in_SSU", "in_RS", "in_DS",
    "log2FC_SSU", "pvalue_SSU", "padj_SSU", "significant_SSU",
    "log2FC_RS", "pvalue_RS", "padj_RS", "significant_RS",
    "log2FC_DS", "pvalue_DS", "padj_DS", "significant_DS"
  ))

  detail
}

wb <- createWorkbook()

summary_rows <- list()
detail_sheets <- list()

for (comparison in COMPARISONS) {
  for (direction in c("Up", "Down")) {
    detail <- make_detail(comparison, direction)
    sheet_name <- paste0(SHEET_PREFIX[[comparison]], "_", direction)
    detail_sheets[[sheet_name]] <- detail

    counts <- detail[, .N, by = .(common_status, region)]
    counts[, comparison := COMPARISON_TITLES[[comparison]]]
    counts[, direction := direction]
    setcolorder(counts, c("comparison", "direction", "common_status", "region", "N"))
    summary_rows[[paste(comparison, direction, sep = "_")]] <- counts
  }
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
setorder(summary_dt, comparison, direction, common_status, region)

header_style <- createStyle(
  fgFill = "#1F4E78",
  fontColour = "white",
  textDecoration = "bold",
  halign = "center",
  border = "Bottom"
)
common_style <- createStyle(fgFill = "#D9EAD3")
uncommon_style <- createStyle(fgFill = "#FFF2CC")

addWorksheet(wb, "Summary")
writeDataTable(wb, "Summary", summary_dt, tableStyle = "TableStyleMedium2")
freezePane(wb, "Summary", firstRow = TRUE)
setColWidths(wb, "Summary", cols = 1:ncol(summary_dt), widths = "auto")
addStyle(wb, "Summary", header_style, rows = 1, cols = 1:ncol(summary_dt), gridExpand = TRUE)

for (sheet_name in names(detail_sheets)) {
  dt <- detail_sheets[[sheet_name]]
  addWorksheet(wb, sheet_name)
  writeDataTable(wb, sheet_name, dt, tableStyle = "TableStyleMedium9")
  freezePane(wb, sheet_name, firstRow = TRUE, firstCol = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(dt), widths = "auto")
  addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(dt), gridExpand = TRUE)

  common_rows <- which(dt$common_status == "Common") + 1
  uncommon_rows <- which(dt$common_status == "Uncommon") + 1
  if (length(common_rows) > 0) {
    addStyle(wb, sheet_name, common_style, rows = common_rows, cols = 1, gridExpand = TRUE, stack = TRUE)
  }
  if (length(uncommon_rows) > 0) {
    addStyle(wb, sheet_name, uncommon_style, rows = uncommon_rows, cols = 1, gridExpand = TRUE, stack = TRUE)
  }
}

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
message("Saved: ", OUT_XLSX)
