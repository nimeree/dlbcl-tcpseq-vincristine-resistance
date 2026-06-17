# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# g:Profiler ORA for P-site fraction limma true-interaction lists
# Up and Down genes combined into one query per fraction.
# Sources: Reactome, GO Biological Process, GO Cellular Component
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
PSITE_LIMMA_DIR <- file.path(BASE_DIR, "Psite_fraction_limma_lfc0.7_rawP0.05")
OUT_DIR <- file.path(PSITE_LIMMA_DIR, "Pathway_gProfiler_Clean_REAC_GOBP_GOCC_true_interaction_combined_direction_lfc0.7_rawP0.05")
TABLE_DIR <- file.path(OUT_DIR, "Tables")
PLOT_DIR <- file.path(OUT_DIR, "Plots")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

SOURCES <- c("REAC", "GO:BP", "GO:CC")
SOURCE_LABELS <- c(REAC = "Reactome", `GO:BP` = "GO Biological Process", `GO:CC` = "GO Cellular Component")
USER_THRESHOLD <- 0.05
MAX_TERM_SIZE <- 500
MIN_INTERSECTION_SIZE <- 5
TOP_N_PER_SOURCE <- 5

FRACTION_LABELS <- c(
  SSU = "SSU P-site count interaction",
  RS = "RS P-site count interaction",
  DS = "DS P-site count interaction"
)

split_intersection <- function(x) {
  if (is.na(x) || x == "") character(0) else unlist(strsplit(x, ",", fixed = TRUE))
}

read_query <- function(frac) {
  f <- file.path(
    PSITE_LIMMA_DIR,
    paste0("Fraction_", frac),
    paste0("Interaction_", frac, "_psite_limma_sig_rawP0.05_lfc0.7.csv")
  )
  if (!file.exists(f)) stop("Missing significant interaction file: ", f)
  dt <- fread(f)
  if ("gene_id_clean" %in% names(dt)) setnames(dt, "gene_id_clean", "gene_id")
  if ("logFC" %in% names(dt)) setnames(dt, "logFC", "log2FoldChange")
  if (!"gene_name" %in% names(dt)) dt[, gene_name := NA_character_]
  if (!"log2FoldChange" %in% names(dt)) dt[, log2FoldChange := NA_real_]
  if (!"direction" %in% names(dt)) {
    dt[, direction := fifelse(log2FoldChange >= 0, "Up", "Down")]
  }
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
  ora_dt[, n_hit_up := vapply(intersection, function(x) {
    ids <- split_intersection(x)
    symbol_map[gene_id %in% ids & direction == "Up", .N]
  }, integer(1))]
  ora_dt[, n_hit_down := vapply(intersection, function(x) {
    ids <- split_intersection(x)
    symbol_map[gene_id %in% ids & direction == "Down", .N]
  }, integer(1))]
  ora_dt[, source_label := SOURCE_LABELS[source]]
  ora_dt[]
}

empty_plot <- function(title, subtitle, out_base) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = subtitle, size = 5) +
    labs(title = title) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(paste0(out_base, ".png"), p, width = 10, height = 5.5, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 10, height = 5.5)
}

select_terms_for_plot <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  plot_dt <- dt[!is.na(source) & !is.na(term_name) & !is.na(p_value)]
  plot_dt <- plot_dt[order(source, p_value, -intersection_size)]
  plot_dt <- plot_dt[, head(.SD, TOP_N_PER_SOURCE), by = source]
  plot_dt[, source := factor(source, levels = SOURCES)]
  plot_dt[, source_label := factor(source_label, levels = SOURCE_LABELS[SOURCES])]
  plot_dt[, term_label := stringr::str_wrap(term_name, width = 42)]
  plot_dt[, neglog10_gscs := -log10(p_value)]
  plot_dt[, gene_ratio := intersection_size / query_size]
  plot_dt[]
}

make_combined_plot <- function(dt, title, out_base) {
  plot_dt <- select_terms_for_plot(dt)
  if (nrow(plot_dt) == 0) {
    empty_plot(title, "No Reactome, GO:BP, or GO:CC terms after filters", out_base)
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
    facet_grid(source_label ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
    scale_size_continuous(range = c(3, 8)) +
    labs(
      title = title,
      subtitle = paste0("Up and down significant genes combined; g:Profiler ORA, g:SCS < ", USER_THRESHOLD,
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
      strip.text.y = element_text(face = "bold", angle = 0),
      strip.background = element_rect(fill = "grey92", color = "grey60"),
      axis.text.y = element_text(color = "black", size = 9.5),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )

  n_terms <- nrow(plot_dt)
  height <- max(6.5, min(11, 2.2 + 0.34 * n_terms))
  ggsave(paste0(out_base, ".png"), p, width = 11.5, height = height, dpi = 300, limitsize = FALSE)
  ggsave(paste0(out_base, ".pdf"), p, width = 11.5, height = height, limitsize = FALSE)
}

all_results <- list()
filtered_results <- list()
plotted_results <- list()
summary_rows <- list()
query_rows <- list()

for (frac in c("SSU", "RS", "DS")) {
  query_dt <- read_query(frac)
  query <- unique(query_dt$gene_id)
  key <- paste("Interaction", frac, "CombinedDirection", sep = "_")
  message("[g:Profiler combined direction] ", key, " (", length(query), " genes)")

  plot_subdir <- file.path(PLOT_DIR, frac)
  dir.create(plot_subdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(query_dt, file.path(TABLE_DIR, paste0(key, "_query_genes.csv")))

  gost_result <- gprofiler2::gost(
    query = query,
    organism = "hsapiens",
    sources = SOURCES,
    correction_method = "g_SCS",
    domain_scope = "annotated",
    user_threshold = USER_THRESHOLD,
    evcodes = TRUE
  )

  if (is.null(gost_result) || is.null(gost_result$result) || nrow(gost_result$result) == 0) {
    gost_dt <- data.table()
  } else {
    gost_dt <- as.data.table(gost_result$result)
    gost_dt <- add_hit_metadata(gost_dt, query_dt)
    gost_dt[, `:=`(
      comparison = "Interaction",
      comparison_label = "True interaction: resistant VCR response vs sensitive VCR response",
      fraction = frac,
      fraction_label = FRACTION_LABELS[[frac]],
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

  if (nrow(gost_dt) == 0 || !"significant" %in% names(gost_dt)) {
    filtered_dt <- data.table()
  } else {
    filtered_dt <- gost_dt[
      significant == TRUE &
        p_value < USER_THRESHOLD &
        term_size <= MAX_TERM_SIZE &
        intersection_size >= MIN_INTERSECTION_SIZE
    ]
  }

  plotted_dt <- select_terms_for_plot(filtered_dt)
  title <- paste(FRACTION_LABELS[[frac]], "combined-direction pathway ORA")
  make_combined_plot(filtered_dt, title, file.path(plot_subdir, paste0(key, "_Reactome_GOBP_GOCC_combined_direction")))

  all_results[[key]] <- gost_dt
  filtered_results[[key]] <- filtered_dt
  plotted_results[[key]] <- plotted_dt
  query_rows[[key]] <- query_dt[, .(
    comparison = "Interaction",
    fraction = frac,
    gene_id,
    gene_name,
    log2FoldChange,
    direction
  )]

  plotted_source_count <- function(src) {
    if (nrow(plotted_dt) == 0 || !"source" %in% names(plotted_dt)) 0L else nrow(plotted_dt[source == src])
  }
  summary_rows[[key]] <- data.table(
    comparison = "Interaction",
    fraction = frac,
    fraction_label = FRACTION_LABELS[[frac]],
    direction = "Combined",
    query_gene_count = length(query),
    query_up_genes = query_dt[direction == "Up", .N],
    query_down_genes = query_dt[direction == "Down", .N],
    gprofiler_significant_terms = nrow(gost_dt),
    filtered_terms = nrow(filtered_dt),
    plotted_terms = nrow(plotted_dt),
    plotted_reactome_terms = plotted_source_count("REAC"),
    plotted_gobp_terms = plotted_source_count("GO:BP"),
    plotted_gocc_terms = plotted_source_count("GO:CC")
  )
}

all_dt <- rbindlist(all_results, fill = TRUE)
filtered_dt <- rbindlist(filtered_results, fill = TRUE)
plotted_dt <- rbindlist(plotted_results, fill = TRUE)
summary_dt <- rbindlist(summary_rows, fill = TRUE)
query_dt_all <- rbindlist(query_rows, fill = TRUE)

preferred_cols <- c(
  "comparison", "comparison_label", "fraction", "fraction_label", "direction",
  "query_gene_count_input", "query_up_genes", "query_down_genes",
  "source", "source_label", "term_id", "term_name", "p_value", "significant",
  "term_size", "query_size", "intersection_size", "precision", "recall",
  "effective_domain_size", "hit_gene_symbols", "hit_gene_ids",
  "mean_hit_log2FC", "n_hit_up", "n_hit_down", "hit_up_genes", "hit_down_genes",
  "correction_method", "domain_scope", "user_threshold",
  "max_term_size_filter", "min_intersection_filter"
)
if (nrow(all_dt) > 0) setcolorder(all_dt, intersect(preferred_cols, names(all_dt)))
if (nrow(filtered_dt) > 0) setcolorder(filtered_dt, intersect(preferred_cols, names(filtered_dt)))
if (nrow(plotted_dt) > 0) setcolorder(plotted_dt, intersect(preferred_cols, names(plotted_dt)))

fwrite(query_dt_all, file.path(TABLE_DIR, "combined_direction_query_genes.csv"))
fwrite(all_dt, file.path(TABLE_DIR, "gprofiler_combined_direction_all_gSCS_significant_terms.csv"))
fwrite(filtered_dt, file.path(TABLE_DIR, "gprofiler_combined_direction_filtered_terms_termSizeLE500_hitGE5.csv"))
fwrite(plotted_dt, file.path(TABLE_DIR, "gprofiler_combined_direction_terms_shown_in_plots.csv"))
fwrite(summary_dt, file.path(TABLE_DIR, "gprofiler_combined_direction_summary_counts.csv"))

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

addWorksheet(wb, "Combined_query_genes")
writeDataTable(wb, "Combined_query_genes", query_dt_all, tableStyle = "TableStyleMedium9")
freezePane(wb, "Combined_query_genes", firstRow = TRUE)
setColWidths(wb, "Combined_query_genes", cols = 1:ncol(query_dt_all), widths = "auto")

saveWorkbook(wb, file.path(OUT_DIR, "gProfiler_REAC_GOBP_GOCC_true_interaction_combined_direction_lfc0.7_rawP0.05.xlsx"), overwrite = TRUE)

cat("\nDone. Combined-direction pathway outputs saved to:\n", OUT_DIR, "\n\n", sep = "")
cat("Summary:\n")
print(summary_dt)
