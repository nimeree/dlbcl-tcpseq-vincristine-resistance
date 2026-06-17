# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# ORA pathway analysis for t2g_v3 DESeq2 outputs
# Databases: Reactome, GO Biological Process, Hallmark
# Test: one-sided hypergeometric / Fisher exact ORA
# DE gene cutoff: raw pvalue < 0.05 and abs(log2FoldChange) >= 0.7
# Pathway significance shown by raw pathway pvalue.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(openxlsx)
})

BASE_DIR <- analysis_path()
OUT_DIR <- file.path(BASE_DIR, "Pathway_ORA_lfc0.7_p0.05_rawP")
PLOT_DIR <- file.path(OUT_DIR, "Plots")
TABLE_DIR <- file.path(OUT_DIR, "Tables")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

P_CUT <- 0.05
LFC_CUT <- 0.7
MIN_PATHWAY_SIZE <- 10
MAX_PATHWAY_SIZE <- 500
TOP_N_PLOT <- 15

COMPARISONS <- c(
  "Sensitive_Vin_vs_DMSO",
  "Resistant_Vin_vs_DMSO",
  "Vin_Resistant_vs_Sensitive"
)
FRACTIONS <- c("SSU", "RS", "DS")
DATABASES <- c("Reactome", "GOBP", "Hallmark")

COMPARISON_TITLES <- c(
  Sensitive_Vin_vs_DMSO = "Sensitive Vin vs DMSO",
  Resistant_Vin_vs_DMSO = "Resistant Vin vs DMSO",
  Vin_Resistant_vs_Sensitive = "Vin Resistant vs Sensitive"
)

clean_pathway_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  stringr::str_to_sentence(tolower(x))
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

first_non_na <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) NA_character_ else x[1]
}

message("[MSigDB] Loading gene sets from msigdbr...")
msig <- as.data.table(msigdbr(species = "Homo sapiens"))
msig[, entrez_id := as.character(ncbi_gene)]

gene_sets <- rbindlist(list(
  msig[gs_collection == "C2" & gs_subcollection == "CP:REACTOME",
       .(database = "Reactome", pathway_id = gs_id, pathway_name = gs_name,
         pathway_description = gs_description, entrez_id, gene_symbol)],
  msig[gs_collection == "C5" & gs_subcollection == "GO:BP",
       .(database = "GOBP", pathway_id = gs_id, pathway_name = gs_name,
         pathway_description = gs_description, entrez_id, gene_symbol)],
  msig[gs_collection == "H",
       .(database = "Hallmark", pathway_id = gs_id, pathway_name = gs_name,
         pathway_description = gs_description, entrez_id, gene_symbol)]
), fill = TRUE)
gene_sets <- unique(gene_sets[!is.na(entrez_id) & entrez_id != ""])
gene_sets[, pathway_label := clean_pathway_name(pathway_name)]

message("[Annot] Building Ensembl -> Entrez/Symbol map...")
all_result_files <- list.files(BASE_DIR, pattern = "_results_all\\.csv$", recursive = TRUE, full.names = TRUE)
all_gene_ids <- unique(unlist(lapply(all_result_files, function(f) {
  dt <- fread(f, select = "gene_id")
  sub("\\.\\d+$", "", dt$gene_id)
})))

id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_gene_ids,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)
id_map <- as.data.table(id_map)
setnames(id_map, "ENSEMBL", "gene_id")
id_map[, gene_id := sub("\\.\\d+$", "", gene_id)]
id_map[, entrez_id := as.character(ENTREZID)]
id_map <- id_map[!is.na(entrez_id) & entrez_id != ""]
id_map <- id_map[, .(
  entrez_id = first_non_na(entrez_id),
  symbol = first_non_na(SYMBOL)
), by = gene_id]

read_result <- function(comparison, fraction) {
  f <- file.path(BASE_DIR, paste0("Fraction_", fraction), paste0(comparison, "_", fraction, "_results_all.csv"))
  if (!file.exists(f)) stop("Missing result file: ", f)
  dt <- fread(f)
  dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
  dt <- merge(dt, id_map, by = "gene_id", all.x = TRUE)
  dt <- dt[!is.na(entrez_id) & entrez_id != ""]
  dt
}

ora_one <- function(query_dt, universe_entrez, db_sets) {
  query_dt <- unique(query_dt[!is.na(entrez_id) & entrez_id != ""], by = "entrez_id")
  query_entrez <- unique(query_dt$entrez_id)
  universe_entrez <- unique(universe_entrez)
  db_sets <- db_sets[entrez_id %in% universe_entrez]

  set_dt <- db_sets[, .(
    pathway_genes = list(unique(entrez_id)),
    pathway_symbols = list(unique(gene_symbol[!is.na(gene_symbol) & gene_symbol != ""])),
    pathway_size = uniqueN(entrez_id),
    pathway_label = first_non_na(pathway_label),
    pathway_description = first_non_na(pathway_description)
  ), by = .(database, pathway_id, pathway_name)]
  set_dt <- set_dt[pathway_size >= MIN_PATHWAY_SIZE & pathway_size <= MAX_PATHWAY_SIZE]

  M <- length(universe_entrez)
  N <- length(query_entrez)
  if (N == 0 || M == 0 || nrow(set_dt) == 0) return(data.table())

  out <- set_dt[, {
    hits <- intersect(query_entrez, pathway_genes[[1]])
    hit_dt <- query_dt[entrez_id %in% hits]
    x <- length(hits)
    K <- pathway_size
    p <- if (x > 0) phyper(x - 1, K, M - K, N, lower.tail = FALSE) else 1
    list(
      raw_pvalue = p,
      hit_count = x,
      query_gene_count = N,
      pathway_gene_count = K,
      universe_gene_count = M,
      gene_ratio = x / N,
      background_ratio = K / M,
      mean_hit_log2FC = if (x > 0) mean(hit_dt$log2FoldChange, na.rm = TRUE) else NA_real_,
      median_hit_log2FC = if (x > 0) median(hit_dt$log2FoldChange, na.rm = TRUE) else NA_real_,
      min_hit_log2FC = if (x > 0) min(hit_dt$log2FoldChange, na.rm = TRUE) else NA_real_,
      max_hit_log2FC = if (x > 0) max(hit_dt$log2FoldChange, na.rm = TRUE) else NA_real_,
      hit_entrez_ids = paste(sort(hits), collapse = ";"),
      hit_gene_symbols = paste(sort(unique(hit_dt$symbol[!is.na(hit_dt$symbol) & hit_dt$symbol != ""])), collapse = ";")
    )
  }, by = .(database, pathway_id, pathway_name, pathway_label, pathway_description)]

  out[, p_adjust_bh := p.adjust(raw_pvalue, method = "BH"), by = database]
  setorder(out, raw_pvalue, -hit_count)
  out
}

make_plot <- function(dt, title, out_base) {
  sig <- dt[raw_pvalue < 0.05 & hit_count > 0]

  if (nrow(sig) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No pathways with raw p < 0.05", size = 5) +
      labs(title = title, subtitle = "ORA: hypergeometric/Fisher exact test") +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  } else {
    plot_dt <- sig[order(raw_pvalue)][1:min(TOP_N_PLOT, .N)]
    plot_dt[, pathway_plot := stringr::str_wrap(pathway_label, width = 42)]
    plot_dt[, neglog10_raw_p := -log10(raw_pvalue)]

    p <- ggplot(plot_dt, aes(
      x = neglog10_raw_p,
      y = forcats::fct_reorder(pathway_plot, neglog10_raw_p)
    )) +
      geom_segment(
        aes(x = 0, xend = neglog10_raw_p,
            y = forcats::fct_reorder(pathway_plot, neglog10_raw_p),
            yend = forcats::fct_reorder(pathway_plot, neglog10_raw_p),
            color = mean_hit_log2FC),
        linewidth = 1.25,
        alpha = 0.9,
        lineend = "round"
      ) +
      geom_point(aes(size = hit_count, color = mean_hit_log2FC), alpha = 0.98) +
      scale_color_gradient2(
        low = "#2B6CB0",
        mid = "grey92",
        high = "#B8323B",
        midpoint = 0,
        na.value = "grey60"
      ) +
      scale_size_continuous(range = c(2.8, 8)) +
      labs(
        title = title,
        subtitle = "Lollipop length = raw-p significance; color = mean log2FC of pathway hit genes",
        x = "-log10(raw pathway p-value)",
        y = NULL,
        size = "Hit genes",
        color = "Mean hit\nlog2FC"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 12),
        axis.text.y = element_text(size = 14, color = "black", lineheight = 0.95),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 13, margin = margin(t = 8)),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = "right",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  }

  ggsave(paste0(out_base, ".png"), p, width = 10.5, height = 7.2, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 10.5, height = 7.2)
}

all_results <- list()
sig_results <- list()
summary_rows <- list()

for (comparison in COMPARISONS) {
  for (fraction in FRACTIONS) {
    res <- read_result(comparison, fraction)
    universe <- unique(res$entrez_id)

    for (direction in c("Up", "Down")) {
      if (direction == "Up") {
        query <- res[pvalue < P_CUT & log2FoldChange >= LFC_CUT,
                     .(entrez_id, symbol, log2FoldChange)]
      } else {
        query <- res[pvalue < P_CUT & log2FoldChange <= -LFC_CUT,
                     .(entrez_id, symbol, log2FoldChange)]
      }

      for (db_name in DATABASES) {
        db_sets <- gene_sets[database == db_name]
        ora <- ora_one(query, universe, db_sets)
        if (nrow(ora) == 0) {
          ora <- data.table(
            database = db_name,
            pathway_id = NA_character_,
            pathway_name = NA_character_,
            pathway_label = NA_character_,
            pathway_description = NA_character_,
            raw_pvalue = NA_real_,
            p_adjust_bh = NA_real_,
            hit_count = 0L,
            query_gene_count = uniqueN(query$entrez_id),
            pathway_gene_count = NA_integer_,
            universe_gene_count = length(universe),
            gene_ratio = NA_real_,
            background_ratio = NA_real_,
            mean_hit_log2FC = NA_real_,
            median_hit_log2FC = NA_real_,
            min_hit_log2FC = NA_real_,
            max_hit_log2FC = NA_real_,
            hit_entrez_ids = NA_character_,
            hit_gene_symbols = NA_character_
          )
        }

        ora[, `:=`(
          comparison = comparison,
          comparison_label = COMPARISON_TITLES[[comparison]],
          fraction = fraction,
          direction = direction,
          de_pvalue_cutoff = P_CUT,
          de_abs_log2fc_cutoff = LFC_CUT,
          pathway_raw_p_significant = !is.na(raw_pvalue) & raw_pvalue < 0.05
        )]

        setcolorder(ora, c(
          "comparison", "comparison_label", "fraction", "direction", "database",
          "pathway_id", "pathway_name", "pathway_label", "pathway_description",
          "raw_pvalue", "p_adjust_bh", "pathway_raw_p_significant",
          "hit_count", "query_gene_count", "pathway_gene_count", "universe_gene_count",
          "gene_ratio", "background_ratio",
          "mean_hit_log2FC", "median_hit_log2FC", "min_hit_log2FC", "max_hit_log2FC",
          "hit_gene_symbols", "hit_entrez_ids",
          "de_pvalue_cutoff", "de_abs_log2fc_cutoff"
        ))

        key <- paste(comparison, fraction, direction, db_name, sep = "_")
        all_results[[key]] <- ora
        sig_results[[key]] <- ora[pathway_raw_p_significant == TRUE & hit_count > 0]

        summary_rows[[key]] <- data.table(
          comparison = comparison,
          fraction = fraction,
          direction = direction,
          database = db_name,
          query_gene_count = uniqueN(query$entrez_id),
          universe_gene_count = length(universe),
          significant_pathway_count_raw_p05 = nrow(sig_results[[key]])
        )

        title <- paste(COMPARISON_TITLES[[comparison]], fraction, direction, db_name)
        plot_subdir <- file.path(PLOT_DIR, db_name, comparison, fraction)
        dir.create(plot_subdir, recursive = TRUE, showWarnings = FALSE)
        out_base <- file.path(plot_subdir, paste0(safe_name(comparison), "_", fraction, "_", direction, "_", db_name, "_ORA_rawP"))
        make_plot(ora, title, out_base)
      }
    }
  }
}

all_dt <- rbindlist(all_results, fill = TRUE)
sig_dt <- rbindlist(sig_results, fill = TRUE)
summary_dt <- rbindlist(summary_rows, fill = TRUE)

fwrite(all_dt, file.path(TABLE_DIR, "all_pathway_ORA_results_raw_pvalues.csv"))
fwrite(sig_dt, file.path(TABLE_DIR, "significant_pathway_ORA_results_raw_pvalue_lt_0.05.csv"))
fwrite(summary_dt, file.path(TABLE_DIR, "pathway_ORA_summary_counts.csv"))

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeDataTable(wb, "Summary", summary_dt, tableStyle = "TableStyleMedium2")
freezePane(wb, "Summary", firstRow = TRUE)
setColWidths(wb, "Summary", cols = 1:ncol(summary_dt), widths = "auto")

addWorksheet(wb, "Significant_raw_p_lt_0.05")
writeDataTable(wb, "Significant_raw_p_lt_0.05", sig_dt, tableStyle = "TableStyleMedium4")
freezePane(wb, "Significant_raw_p_lt_0.05", firstRow = TRUE)
setColWidths(wb, "Significant_raw_p_lt_0.05", cols = 1:ncol(sig_dt), widths = "auto")

addWorksheet(wb, "All_ORA_results")
writeDataTable(wb, "All_ORA_results", all_dt, tableStyle = "TableStyleMedium9")
freezePane(wb, "All_ORA_results", firstRow = TRUE)
setColWidths(wb, "All_ORA_results", cols = 1:ncol(all_dt), widths = "auto")

saveWorkbook(wb, file.path(OUT_DIR, "Pathway_ORA_Reactome_GOBP_Hallmark_lfc0.7_p0.05_rawP.xlsx"), overwrite = TRUE)

message("Done. Pathway ORA outputs saved to: ", OUT_DIR)
