# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# Venn diagrams for t2g_v3 DESeq2 genome length-filtered outputs
# - 3-set Venns: SSU, RS, DS
# - Separate Up and Down diagrams for each comparison
# - Threshold: pvalue < 0.05 and abs(log2FoldChange) >= 0.7
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggVennDiagram)
  library(openxlsx)
})

IN_DIR <- analysis_path()
OUT_DIR <- file.path(IN_DIR, "Venn_Diagrams_lfc0.7")
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

read_gene_set <- function(comparison, fraction, direction) {
  f <- file.path(IN_DIR, paste0("Fraction_", fraction), paste0(comparison, "_", fraction, "_results_all.csv"))
  if (!file.exists(f)) stop("Missing input file: ", f)

  dt <- fread(f)
  dt <- dt[!is.na(pvalue) & !is.na(log2FoldChange)]

  if (direction == "Up") {
    dt <- dt[pvalue < P_CUT & log2FoldChange >= LFC_CUT]
  } else if (direction == "Down") {
    dt <- dt[pvalue < P_CUT & log2FoldChange <= -LFC_CUT]
  } else {
    stop("direction must be Up or Down")
  }

  unique(dt$gene_id)
}

classify_regions <- function(sets) {
  ssu <- sets$SSU
  rs <- sets$RS
  ds <- sets$DS

  list(
    SSU_only = setdiff(ssu, union(rs, ds)),
    RS_only = setdiff(rs, union(ssu, ds)),
    DS_only = setdiff(ds, union(ssu, rs)),
    SSU_RS_only = setdiff(intersect(ssu, rs), ds),
    SSU_DS_only = setdiff(intersect(ssu, ds), rs),
    RS_DS_only = setdiff(intersect(rs, ds), ssu),
    SSU_RS_DS = Reduce(intersect, list(ssu, rs, ds))
  )
}

make_venn_plot <- function(sets, title, subtitle) {
  ggVennDiagram::ggVennDiagram(
    sets,
    label = "count",
    label_alpha = 0,
    edge_size = 0.8,
    set_size = 5
  ) +
    scale_fill_gradient(low = "#F7FBFF", high = "#3182BD") +
    labs(
      title = title,
      subtitle = subtitle
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
      plot.margin = margin(14, 14, 14, 14)
    )
}

summary_rows <- list()

for (comparison in COMPARISONS) {
  for (direction in c("Up", "Down")) {
    sets <- setNames(lapply(FRACTIONS, function(frac) read_gene_set(comparison, frac, direction)), FRACTIONS)
    regions <- classify_regions(sets)
    region_counts <- vapply(regions, length, integer(1))

    title <- paste(COMPARISON_TITLES[[comparison]], direction, "Genes")
    subtitle <- paste0("p < ", P_CUT, " and ",
                       ifelse(direction == "Up", "log2FC >= ", "log2FC <= -"),
                       LFC_CUT)

    p <- make_venn_plot(sets, title, subtitle)

    out_base <- paste0(comparison, "_", direction, "_SSU_RS_DS_venn_lfc0.7")
    ggsave(file.path(OUT_DIR, paste0(out_base, ".png")), p, width = 7.2, height = 6.4, dpi = 300)
    ggsave(file.path(OUT_DIR, paste0(out_base, ".pdf")), p, width = 7.2, height = 6.4)

    wb <- createWorkbook()
    for (nm in names(regions)) {
      addWorksheet(wb, nm)
      writeData(wb, nm, data.table(gene_id = sort(regions[[nm]])))
    }
    saveWorkbook(wb, file.path(OUT_DIR, paste0(out_base, "_gene_lists.xlsx")), overwrite = TRUE)

    summary_rows[[paste(comparison, direction, sep = "_")]] <- data.table(
      comparison = comparison,
      direction = direction,
      total_SSU = length(sets$SSU),
      total_RS = length(sets$RS),
      total_DS = length(sets$DS),
      region = names(region_counts),
      count = as.integer(region_counts)
    )
  }
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, file.path(OUT_DIR, "venn_region_counts_summary.csv"))

message("Done. Venn diagrams and region gene lists saved to: ", OUT_DIR)
