# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# GO Biological Process-only g:Profiler ORA
# P-site fraction-specific limma count analysis
# Up + Down significant genes combined per fraction and contrast.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(gprofiler2)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(openxlsx)
})

BASE_DIR <- analysis_path()
PSITE_DIR <- file.path(BASE_DIR, "Psite_fraction_limma_lfc0.7_rawP0.05")
OUT_DIR <- file.path(PSITE_DIR, "GO_BP_ORA_All_Contrasts_Combined")

SOURCE <- "GO:BP"
USER_THRESHOLD <- 0.05
MAX_TERM_SIZE <- 500
MIN_INTERSECTION_SIZE <- 5
TOP_N <- 15

CONTRASTS <- c(
  "Resistance_baseline",
  "Interaction"
)

CONTRAST_LABELS <- c(
  Sensitive_Vin_vs_DMSO = "Sensitive Vin vs DMSO",
  Resistant_Vin_vs_DMSO = "Resistant Vin vs DMSO",
  Vin_Resistant_vs_Sensitive = "Resistant vs sensitive under Vin",
  Resistance_baseline = "Baseline resistance",
  Interaction = "True interaction"
)

CONTRAST_NICE <- c(
  Sensitive_Vin_vs_DMSO = "Sensitive_Vin_vs_DMSO",
  Resistant_Vin_vs_DMSO = "Resistant_Vin_vs_DMSO",
  Vin_Resistant_vs_Sensitive = "Resistant_vs_Sensitive_under_Vin",
  Resistance_baseline = "Baseline_resistance",
  Interaction = "True_interaction"
)

FRACTIONS <- c("SSU", "RS", "DS")

split_intersection <- function(x) {
  if (is.na(x) || x == "") character(0) else unlist(strsplit(x, ",", fixed = TRUE))
}

sig_file <- function(contrast, frac) {
  file.path(
    PSITE_DIR,
    paste0("Fraction_", frac),
    paste0(contrast, "_", frac, "_psite_limma_sig_rawP0.05_lfc0.7.csv")
  )
}

read_query <- function(contrast, frac) {
  f <- sig_file(contrast, frac)
  if (!file.exists(f)) stop("Missing significant limma file: ", f)
  dt <- fread(f)
  if ("gene_id_clean" %in% names(dt)) setnames(dt, "gene_id_clean", "gene_id")
  if ("logFC" %in% names(dt)) setnames(dt, "logFC", "log2FoldChange")
  if (!"gene_name" %in% names(dt)) dt[, gene_name := NA_character_]
  if (!"direction" %in% names(dt)) dt[, direction := fifelse(log2FoldChange >= 0, "Up", "Down")]
  dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
  unique(dt[!is.na(gene_id) & gene_id != "", .(gene_id, gene_name, log2FoldChange, direction)], by = "gene_id")
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
      subtitle = paste0("GO Biological Process only; Up + Down combined; g:SCS < ", USER_THRESHOLD,
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

run_one <- function(contrast, frac) {
  query_dt <- read_query(contrast, frac)
  query <- unique(query_dt$gene_id)
  contrast_nice <- CONTRAST_NICE[[contrast]]
  contrast_label <- CONTRAST_LABELS[[contrast]]

  message("[GO:BP ORA fraction limma] ", contrast_nice, " / ", frac, " (", length(query), " genes)")

  tables_dir <- file.path(OUT_DIR, contrast_nice, "Tables")
  plots_dir <- file.path(OUT_DIR, contrast_nice, "Plots", frac)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(query_dt, file.path(tables_dir, paste0(frac, "_query_genes.csv")))

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
        dataset = "P-site fraction limma",
        contrast = contrast,
        contrast_label = contrast_label,
        fraction = frac,
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

  fwrite(all_dt, file.path(tables_dir, paste0(frac, "_all_gSCS_significant_terms.csv")))
  fwrite(filtered_dt, file.path(tables_dir, paste0(frac, "_filtered_terms.csv")))
  fwrite(plotted_dt, file.path(tables_dir, paste0(frac, "_terms_shown_in_plot.csv")))

  make_plot(
    filtered_dt,
    paste(contrast_label, frac, "GO:BP ORA"),
    file.path(plots_dir, paste0(frac, "_GO_BP_ORA_combined_direction"))
  )

  list(
    all = all_dt,
    filtered = filtered_dt,
    plotted = plotted_dt,
    query = query_dt[, .(dataset = "P-site fraction limma", contrast, contrast_label, fraction = frac, gene_id, gene_name, log2FoldChange, direction)],
    summary = data.table(
      dataset = "P-site fraction limma",
      contrast = contrast,
      contrast_label = contrast_label,
      fraction = frac,
      query_gene_count = length(query),
      query_up_genes = query_dt[direction == "Up", .N],
      query_down_genes = query_dt[direction == "Down", .N],
      gprofiler_significant_terms = nrow(all_dt),
      filtered_terms = nrow(filtered_dt),
      plotted_terms = nrow(plotted_dt)
    )
  )
}

all_results <- list()
for (contrast in CONTRASTS) {
  for (frac in FRACTIONS) {
    all_results[[paste(contrast, frac, sep = "_")]] <- run_one(contrast, frac)
  }
}

all_dt <- rbindlist(lapply(all_results, `[[`, "all"), fill = TRUE)
filtered_dt <- rbindlist(lapply(all_results, `[[`, "filtered"), fill = TRUE)
plotted_dt <- rbindlist(lapply(all_results, `[[`, "plotted"), fill = TRUE)
query_dt <- rbindlist(lapply(all_results, `[[`, "query"), fill = TRUE)
summary_dt <- rbindlist(lapply(all_results, `[[`, "summary"), fill = TRUE)

dir.create(file.path(OUT_DIR, "All_Tables"), recursive = TRUE, showWarnings = FALSE)
fwrite(summary_dt, file.path(OUT_DIR, "All_Tables", "summary.csv"))
fwrite(query_dt, file.path(OUT_DIR, "All_Tables", "combined_query_genes.csv"))
fwrite(all_dt, file.path(OUT_DIR, "All_Tables", "all_gSCS_significant_terms.csv"))
fwrite(filtered_dt, file.path(OUT_DIR, "All_Tables", "filtered_terms.csv"))
fwrite(plotted_dt, file.path(OUT_DIR, "All_Tables", "terms_shown_in_plots.csv"))

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

saveWorkbook(wb, file.path(OUT_DIR, "GO_BP_ORA_Psite_fraction_limma_all_contrasts_combined.xlsx"), overwrite = TRUE)

cat("\nDone. GO:BP-only fraction-specific limma ORA outputs saved to:\n", OUT_DIR, "\n\n", sep = "")
cat("Summary:\n")
print(summary_dt)
