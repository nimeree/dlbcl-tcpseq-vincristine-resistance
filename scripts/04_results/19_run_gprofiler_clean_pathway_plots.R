# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# Clean g:Profiler ORA pathway plots for limma metric up/down lists
# Sources: Reactome, GO Biological Process, GO Cellular Component
# Statistics: gprofiler2::gost() with g:SCS correction
# Background: domain_scope = "annotated"
# Filters: p < 0.05, term_size <= 500, intersection_size >= 5
# Plot: one combined figure per list, top terms per source
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
LIMMA_DIR <- file.path(BASE_DIR, "Limma_translation_metrics_lfc0.7_rawP0.05")
LIMMA_RESULTS_DIR <- file.path(LIMMA_DIR, "Results")
OUT_DIR <- file.path(LIMMA_DIR, "Pathway_gProfiler_Clean_REAC_GOBP_GOCC_lfc0.7_rawP0.05")
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

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

parse_input_file <- function(f) {
  b <- basename(f)
  metric <- basename(dirname(f))
  direction <- sub("^.*_limma_sig_(up|down)_rawP0\\.05_lfc0\\.7\\.csv$", "\\1", b)
  contrast <- sub("_limma_sig_(up|down)_rawP0\\.05_lfc0\\.7\\.csv$", "", b)
  data.table(
    file = f,
    comparison = contrast,
    comparison_label = if (contrast %in% names(CONTRAST_LABELS)) CONTRAST_LABELS[[contrast]] else contrast,
    fraction = metric,
    fraction_label = if (metric %in% names(METRIC_LABELS)) METRIC_LABELS[[metric]] else metric,
    direction = stringr::str_to_title(direction)
  )
}

read_query <- function(f) {
  dt <- fread(f)
  if ("gene_id_clean" %in% names(dt)) {
    setnames(dt, "gene_id_clean", "gene_id")
  }
  if ("logFC" %in% names(dt)) {
    setnames(dt, "logFC", "log2FoldChange")
  }
  if (!"gene_name" %in% names(dt)) {
    dt[, gene_name := NA_character_]
  }
  if (!"log2FoldChange" %in% names(dt)) {
    dt[, log2FoldChange := NA_real_]
  }
  dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
  unique(dt[!is.na(gene_id) & gene_id != ""], by = "gene_id")
}

split_intersection <- function(x) {
  if (is.na(x) || x == "") character(0) else unlist(strsplit(x, ",", fixed = TRUE))
}

add_hit_metadata <- function(ora_dt, query_dt) {
  symbol_map <- unique(query_dt[, .(gene_id, gene_name, log2FoldChange)])
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
    geom_segment(aes(x = 0, xend = neglog10_gscs,
                     y = forcats::fct_reorder(term_label, neglog10_gscs),
                     yend = forcats::fct_reorder(term_label, neglog10_gscs)),
                 linewidth = 0.8, color = "grey78") +
    geom_point(aes(size = intersection_size, color = gene_ratio), alpha = 0.95) +
    facet_grid(source_label ~ ., scales = "free_y", space = "free_y") +
    scale_color_viridis_c(option = "plasma", direction = -1) +
    scale_size_continuous(range = c(3, 8)) +
    labs(
      title = title,
      subtitle = paste0("g:Profiler ORA, g:SCS < ", USER_THRESHOLD,
                        "; term size <= ", MAX_TERM_SIZE,
                        "; hit genes >= ", MIN_INTERSECTION_SIZE,
                        "; up to ", TOP_N_PER_SOURCE, " terms per source"),
      x = "-log10(g:SCS corrected p-value)",
      y = NULL,
      size = "Hit genes",
      color = "Gene ratio"
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

input_files <- list.files(
  LIMMA_RESULTS_DIR,
  pattern = "_limma_sig_(up|down)_rawP0\\.05_lfc0\\.7\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
input_files <- input_files[!grepl("_all_contrasts\\.csv$", input_files)]
if (length(input_files) == 0) stop("No limma significant up/down CSV files found in: ", LIMMA_RESULTS_DIR)
inputs <- rbindlist(lapply(input_files, parse_input_file), fill = TRUE)
setorder(inputs, comparison, fraction, direction)

all_results <- list()
filtered_results <- list()
plotted_results <- list()
summary_rows <- list()

for (i in seq_len(nrow(inputs))) {
  info <- inputs[i]
  query_dt <- read_query(info$file)
  query <- unique(query_dt$gene_id)
  key <- paste(info$comparison, info$fraction, info$direction, sep = "_")
  message("[g:Profiler clean] ", key, " (", length(query), " genes)")

  plot_subdir <- file.path(PLOT_DIR, safe_name(info$comparison), info$fraction, info$direction)
  dir.create(plot_subdir, recursive = TRUE, showWarnings = FALSE)

  if (length(query) == 0) {
    gost_dt <- data.table()
  } else {
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
        comparison = info$comparison,
        comparison_label = info$comparison_label,
        fraction = info$fraction,
        fraction_label = info$fraction_label,
        direction = info$direction,
        query_gene_count_input = length(query),
        correction_method = "g_SCS",
        domain_scope = "annotated",
        user_threshold = USER_THRESHOLD,
        max_term_size_filter = MAX_TERM_SIZE,
        min_intersection_filter = MIN_INTERSECTION_SIZE
      )]
    }
    Sys.sleep(0.5)
  }

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

  title <- paste(info$comparison_label, info$fraction_label, info$direction)
  make_combined_plot(filtered_dt, paste(title, "pathway ORA"), file.path(plot_subdir, paste0(key, "_Reactome_GOBP_GOCC_combined")))

  all_results[[key]] <- gost_dt
  filtered_results[[key]] <- filtered_dt
  plotted_results[[key]] <- plotted_dt
  plotted_source_count <- function(src) {
    if (nrow(plotted_dt) == 0 || !"source" %in% names(plotted_dt)) 0L else nrow(plotted_dt[source == src])
  }
  summary_rows[[key]] <- data.table(
    comparison = info$comparison,
    comparison_label = info$comparison_label,
    fraction = info$fraction,
    fraction_label = info$fraction_label,
    direction = info$direction,
    query_gene_count = length(query),
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

preferred_cols <- c(
  "comparison", "comparison_label", "fraction", "fraction_label", "direction", "query_gene_count_input",
  "source", "source_label", "term_id", "term_name", "p_value", "significant",
  "term_size", "query_size", "intersection_size", "precision", "recall",
  "effective_domain_size", "gene_ratio", "hit_gene_symbols", "hit_gene_ids",
  "mean_hit_log2FC", "correction_method", "domain_scope", "user_threshold",
  "max_term_size_filter", "min_intersection_filter"
)
if (nrow(all_dt) > 0) setcolorder(all_dt, intersect(preferred_cols, names(all_dt)))
if (nrow(filtered_dt) > 0) setcolorder(filtered_dt, intersect(preferred_cols, names(filtered_dt)))
if (nrow(plotted_dt) > 0) setcolorder(plotted_dt, intersect(preferred_cols, names(plotted_dt)))

fwrite(all_dt, file.path(TABLE_DIR, "gprofiler_clean_all_gSCS_significant_terms.csv"))
fwrite(filtered_dt, file.path(TABLE_DIR, "gprofiler_clean_filtered_terms_termSizeLE500_hitGE5.csv"))
fwrite(plotted_dt, file.path(TABLE_DIR, "gprofiler_clean_terms_shown_in_plots.csv"))
fwrite(summary_dt, file.path(TABLE_DIR, "gprofiler_clean_summary_counts.csv"))

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

addWorksheet(wb, "All_gSCS_terms")
writeDataTable(wb, "All_gSCS_terms", all_dt, tableStyle = "TableStyleMedium9")
freezePane(wb, "All_gSCS_terms", firstRow = TRUE)
if (ncol(all_dt) > 0) setColWidths(wb, "All_gSCS_terms", cols = 1:ncol(all_dt), widths = "auto")

saveWorkbook(wb, file.path(OUT_DIR, "gProfiler_clean_REAC_GOBP_GOCC_lfc0.7_rawP0.05.xlsx"), overwrite = TRUE)
message("Done. Clean pathway outputs saved to: ", OUT_DIR)
