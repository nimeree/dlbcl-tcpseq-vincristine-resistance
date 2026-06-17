# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# Translation-state classification from t2g_v3 DESeq2 outputs
# - Uses SSU, RS, DS DESeq2 log2FC/pvalue patterns
# - Threshold: pvalue < 0.05 and abs(log2FoldChange) >= 0.7
# - Outputs count barplot and Excel gene lists
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
})

BASE_DIR <- analysis_path()
OUT_DIR <- file.path(BASE_DIR, "Translation_State_Classification_lfc0.7_p0.05")
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

CLASS_LEVELS <- c(
  "Efficient translation",
  "Scanning bottleneck",
  "Collision / slow translation",
  "Translation repression",
  "Disome-specific stalling",
  "Reduced ribosome occupancy",
  "Other significant pattern"
)

CLASS_DESCRIPTIONS <- data.table(
  translation_class = CLASS_LEVELS,
  logic = c(
    "SSU Up, RS Up, DS Down",
    "SSU Up, RS Down",
    "RS Up, DS Up",
    "SSU Down, RS Down",
    "DS Up, without SSU Up or RS Up",
    "RS Down, DS Down",
    "At least one significant SSU/RS/DS change, but not one of the defined classes"
  ),
  interpretation = c(
    "More small-subunit loading and 80S ribosome signal, with fewer disomes/collisions.",
    "Small-subunit signal accumulates but does not progress into translating ribosomes.",
    "More translating ribosomes together with more disomes, consistent with queuing or slow elongation.",
    "Reduced initiation/loading and reduced translating ribosome occupancy.",
    "Disome/collision signal without broad SSU/RS upregulation.",
    "Reduced elongating ribosome signal and reduced disome signal.",
    "Significant translational change with a mixed or less interpretable SSU/RS/DS pattern."
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
  dt[, .(
    gene_id,
    log2FoldChange,
    pvalue,
    padj,
    status
  )]
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
    missing <- setdiff(needed, names(d))
    if (length(missing) > 0) return(NULL)
    d[, gene_id := sub("\\.\\d+$", "", gene_id)]
    d[, ..needed]
  }), fill = TRUE)

  ann <- ann[!is.na(gene_id) & gene_id != ""]
  ann[, .(
    gene_name = first(na.omit(gene_name)),
    gene_function = first(na.omit(gene_function)),
    gene_summary = first(na.omit(gene_summary))
  ), by = gene_id]
}

first_or_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[1]
}

classify_patterns <- function(dt) {
  dt[, translation_class := fcase(
    SSU_status == "Up" & RS_status == "Up" & DS_status == "Down",
    "Efficient translation",

    SSU_status == "Up" & RS_status == "Down",
    "Scanning bottleneck",

    RS_status == "Up" & DS_status == "Up",
    "Collision / slow translation",

    SSU_status == "Down" & RS_status == "Down",
    "Translation repression",

    DS_status == "Up" & SSU_status != "Up" & RS_status != "Up",
    "Disome-specific stalling",

    RS_status == "Down" & DS_status == "Down",
    "Reduced ribosome occupancy",

    SSU_status != "NS" | RS_status != "NS" | DS_status != "NS",
    "Other significant pattern",

    default = NA_character_
  )]
  dt[!is.na(translation_class)]
}

annotation <- read_annotation()
all_classified <- list()

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

  classified <- classify_patterns(merged)
  classified[, comparison := comparison]
  classified <- merge(classified, annotation, by = "gene_id", all.x = TRUE)
  classified[, gene_name := fifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)]
  all_classified[[comparison]] <- classified
}

classified_dt <- rbindlist(all_classified, fill = TRUE)
classified_dt[, translation_class := factor(translation_class, levels = CLASS_LEVELS)]
classified_dt[, comparison_title := COMPARISON_TITLES[comparison]]
classified_dt[, comparison_title := factor(comparison_title, levels = COMPARISON_TITLES[COMPARISONS])]

count_dt <- classified_dt[, .(gene_count = .N), by = .(comparison, comparison_title, translation_class)]
all_grid <- CJ(comparison = COMPARISONS, translation_class = CLASS_LEVELS)
all_grid[, comparison_title := COMPARISON_TITLES[comparison]]
all_grid[, comparison_title := factor(comparison_title, levels = COMPARISON_TITLES[COMPARISONS])]
all_grid[, translation_class := factor(translation_class, levels = CLASS_LEVELS)]
count_dt <- merge(all_grid, count_dt, by = c("comparison", "comparison_title", "translation_class"), all.x = TRUE)
count_dt[is.na(gene_count), gene_count := 0L]

class_colors <- c(
  "Efficient translation" = "#2C7A7B",
  "Scanning bottleneck" = "#D69E2E",
  "Collision / slow translation" = "#B8323B",
  "Translation repression" = "#4A5568",
  "Disome-specific stalling" = "#805AD5",
  "Reduced ribosome occupancy" = "#2B6CB0",
  "Other significant pattern" = "#718096"
)

p <- ggplot(count_dt, aes(x = translation_class, y = gene_count, fill = translation_class)) +
  geom_col(width = 0.72, alpha = 0.94) +
  geom_text(aes(label = gene_count), vjust = -0.35, size = 4.2, fontface = "bold") +
  facet_wrap(~ comparison_title, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = class_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = "Translation-State Classification of Differentially Translated Genes",
    subtitle = paste0("DESeq2 thresholds: p < ", P_CUT, " and abs(log2FC) >= ", LFC_CUT),
    x = NULL,
    y = "Number of genes"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 11, color = "black", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 18, 12, 18)
  )

plot_base <- file.path(OUT_DIR, "translation_state_classification_count_barplot_lfc0.7_p0.05")
ggsave(paste0(plot_base, ".png"), p, width = 14, height = 7.5, dpi = 300)
ggsave(paste0(plot_base, ".pdf"), p, width = 14, height = 7.5)

fwrite(count_dt, file.path(OUT_DIR, "translation_state_classification_counts.csv"))
fwrite(classified_dt, file.path(OUT_DIR, "translation_state_classified_gene_list.csv"))

wb <- createWorkbook()
addWorksheet(wb, "Class definitions")
writeData(wb, "Class definitions", CLASS_DESCRIPTIONS)

addWorksheet(wb, "Counts")
writeData(wb, "Counts", count_dt[order(comparison, translation_class)])

excel_cols <- c(
  "comparison", "translation_class", "gene_id", "gene_name",
  "SSU_status", "SSU_log2FC", "SSU_pvalue", "SSU_padj",
  "RS_status", "RS_log2FC", "RS_pvalue", "RS_padj",
  "DS_status", "DS_log2FC", "DS_pvalue", "DS_padj",
  "gene_function", "gene_summary"
)

for (comp in COMPARISONS) {
  sheet_name <- gsub("_", " ", comp)
  addWorksheet(wb, sheet_name)
  out <- classified_dt[comparison == comp]
  out <- out[order(translation_class, gene_name)]
  writeData(wb, sheet_name, out[, ..excel_cols])
}

for (class_name in CLASS_LEVELS) {
  sheet_name <- substr(gsub("[/]", "-", class_name), 1, 31)
  addWorksheet(wb, sheet_name)
  out <- classified_dt[as.character(translation_class) == class_name]
  out <- out[order(comparison, gene_name)]
  writeData(wb, sheet_name, out[, ..excel_cols])
}

saveWorkbook(
  wb,
  file.path(OUT_DIR, "translation_state_classified_gene_lists_lfc0.7_p0.05.xlsx"),
  overwrite = TRUE
)

message("Done. Translation-state classification outputs saved to: ", OUT_DIR)
