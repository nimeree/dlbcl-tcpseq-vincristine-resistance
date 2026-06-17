# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- analysis_path()
psite_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
metric_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05")
out_dir <- file.path(base_dir, "Cell_cycle_deep_dive")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_cut <- 0.05
lfc_cut <- 0.7

cell_cycle_pattern <- paste(
  c(
    "cell cycle", "mitotic", "mitosis", "chromosome segregation", "nuclear division",
    "DNA replication", "M phase", "G1", "G2", "S phase", "spindle", "centromere",
    "chromatid", "cytokinesis"
  ),
  collapse = "|"
)

split_ids <- function(x) {
  if (is.na(x) || x == "") character(0) else unlist(strsplit(x, ",", fixed = TRUE))
}

direction_call <- function(logfc, pval) {
  fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc >= lfc_cut, "Up",
          fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc <= -lfc_cut, "Down", "NS"))
}

majority_direction <- function(up, down) {
  fifelse(up > down, "Mostly Up", fifelse(down > up, "Mostly Down", "Mixed"))
}

read_go_tables <- function() {
  psite_file <- file.path(psite_dir, "GO_BP_ORA_All_Contrasts_Combined", "All_Tables", "filtered_terms.csv")
  metric_baseline_file <- file.path(metric_dir, "GO_BP_ORA_Baseline_Combined", "Tables", "filtered_terms.csv")
  metric_interaction_file <- file.path(metric_dir, "GO_BP_ORA_Interaction_Combined", "Tables", "filtered_terms.csv")

  psite <- fread(psite_file)
  psite[, dataset := "Fraction limma counts"]
  psite[, contrast_group := fifelse(contrast == "Resistance_baseline", "Baseline", "Interaction")]
  psite[, analysis := fraction]
  psite[, analysis_label := fraction]

  metric_base <- fread(metric_baseline_file)
  metric_base[, dataset := "Translation metric limma"]
  metric_base[, contrast_group := "Baseline"]
  if (!"analysis" %in% names(metric_base)) metric_base[, analysis := NA_character_]

  metric_int <- fread(metric_interaction_file)
  metric_int[, dataset := "Translation metric limma"]
  metric_int[, contrast_group := "Interaction"]
  metric_int[, contrast := "Interaction"]
  if (!"analysis" %in% names(metric_int)) metric_int[, analysis := NA_character_]

  go <- rbindlist(list(psite, metric_base, metric_int), fill = TRUE)
  go[, term_majority := majority_direction(n_hit_up, n_hit_down)]
  go[, cell_cycle_related := grepl(cell_cycle_pattern, term_name, ignore.case = TRUE)]
  go[]
}

extract_term_hit_genes <- function(go_dt) {
  cc <- go_dt[cell_cycle_related == TRUE]
  if (nrow(cc) == 0) return(data.table())
  rbindlist(lapply(seq_len(nrow(cc)), function(i) {
    ids <- unique(split_ids(cc$intersection[i]))
    if (length(ids) == 0) return(data.table())
    data.table(
      dataset = cc$dataset[i],
      contrast_group = cc$contrast_group[i],
      analysis = cc$analysis[i],
      analysis_label = cc$analysis_label[i],
      term_id = cc$term_id[i],
      term_name = cc$term_name[i],
      term_majority = cc$term_majority[i],
      gene_id = ids
    )
  }), fill = TRUE)
}

read_psite_all <- function(frac, contrast_group) {
  contrast <- if (contrast_group == "Baseline") "Resistance_baseline" else "Interaction"
  f <- file.path(psite_dir, paste0("Fraction_", frac), paste0(contrast, "_", frac, "_psite_limma_all_genes.csv"))
  d <- fread(f)
  d[, gene_id := sub("\\.\\d+$", "", gene_id_clean)]
  d[, gene_symbol := gene_name]
  d[, direction := direction_call(logFC, P.Value)]
  d[, .(gene_id, gene_symbol, logFC, P.Value, direction)]
}

metric_files <- data.table(
  metric = c("scanning", "ribosome_engagement", "collision", "protein_output"),
  folder = c("scanning_score", "ribosome_efficiency_score", "collision_score", "protein_output_score"),
  label = c("Scanning", "Ribosome engagement", "Collision", "Protein output")
)

read_metric_all <- function(metric_folder, contrast_group) {
  contrast <- if (contrast_group == "Baseline") "Resistance_baseline" else "Interaction"
  f <- file.path(metric_dir, "Results", metric_folder, paste0(contrast, "_limma_all_genes.csv"))
  d <- fread(f)
  d[, gene_id := sub("\\.\\d+$", "", gene_id_clean)]
  d[, gene_symbol := gene_name]
  d[, direction := direction_call(logFC, P.Value)]
  d[, .(gene_id, gene_symbol, logFC, P.Value, direction)]
}

summarise_lfc_for_gene_set <- function(dt, cell_cycle_ids) {
  dt[, is_cell_cycle := gene_id %in% cell_cycle_ids]
  cell <- dt[is_cell_cycle == TRUE]
  bg <- dt[is_cell_cycle == FALSE]
  wt <- if (nrow(cell) > 1 && nrow(bg) > 1) {
    suppressWarnings(wilcox.test(cell$logFC, bg$logFC))
  } else {
    list(p.value = NA_real_)
  }
  data.table(
    tested_cell_cycle_genes = nrow(cell),
    sig_up = cell[direction == "Up", .N],
    sig_down = cell[direction == "Down", .N],
    sig_any = cell[direction != "NS", .N],
    pct_sig_up = ifelse(nrow(cell) == 0, NA_real_, 100 * cell[direction == "Up", .N] / nrow(cell)),
    pct_sig_down = ifelse(nrow(cell) == 0, NA_real_, 100 * cell[direction == "Down", .N] / nrow(cell)),
    median_logFC_cell_cycle = median(cell$logFC, na.rm = TRUE),
    median_logFC_background = median(bg$logFC, na.rm = TRUE),
    wilcox_cell_cycle_vs_background_p = wt$p.value
  )
}

go <- read_go_tables()
cc_hits <- extract_term_hit_genes(go)
cell_cycle_ids <- sort(unique(cc_hits$gene_id))

pathway_summary <- go[, .(
  total_gobp_terms = .N,
  cell_cycle_terms = sum(cell_cycle_related),
  pct_cell_cycle_terms = 100 * sum(cell_cycle_related) / .N,
  cell_cycle_mostly_up = sum(cell_cycle_related & term_majority == "Mostly Up"),
  cell_cycle_mostly_down = sum(cell_cycle_related & term_majority == "Mostly Down"),
  cell_cycle_mixed = sum(cell_cycle_related & term_majority == "Mixed")
), by = .(dataset, contrast_group, analysis_label)]
setorder(pathway_summary, dataset, contrast_group, analysis_label)

cc_term_details <- go[cell_cycle_related == TRUE, .(
  dataset, contrast_group, analysis_label, term_id, term_name,
  p_value, intersection_size, n_hit_up, n_hit_down, term_majority,
  mean_hit_log2FC
)]
setorder(cc_term_details, dataset, contrast_group, analysis_label, p_value)

hit_gene_summary <- cc_hits[, .(
  unique_cell_cycle_hit_genes = uniqueN(gene_id),
  cell_cycle_terms_contributing = uniqueN(term_id)
), by = .(dataset, contrast_group, analysis_label)]
setorder(hit_gene_summary, dataset, contrast_group, analysis_label)

fraction_gene_summary <- rbindlist(lapply(c("Baseline", "Interaction"), function(ct) {
  rbindlist(lapply(c("SSU", "RS", "DS"), function(frac) {
    d <- read_psite_all(frac, ct)
    cbind(
      data.table(dataset = "Fraction limma counts", contrast_group = ct, analysis = frac),
      summarise_lfc_for_gene_set(d, cell_cycle_ids)
    )
  }))
}), fill = TRUE)

metric_gene_summary <- rbindlist(lapply(c("Baseline", "Interaction"), function(ct) {
  rbindlist(lapply(seq_len(nrow(metric_files)), function(i) {
    d <- read_metric_all(metric_files$folder[i], ct)
    cbind(
      data.table(dataset = "Translation metric limma", contrast_group = ct, analysis = metric_files$label[i]),
      summarise_lfc_for_gene_set(d, cell_cycle_ids)
    )
  }))
}), fill = TRUE)

metric_pattern_summary <- rbindlist(lapply(c("Baseline", "Interaction"), function(ct) {
  metric_dt <- Reduce(function(x, y) merge(x, y, by = c("gene_id", "gene_symbol"), all = TRUE), lapply(seq_len(nrow(metric_files)), function(i) {
    d <- read_metric_all(metric_files$folder[i], ct)
    setnames(d, c("logFC", "P.Value", "direction"), paste0(metric_files$metric[i], c("_logFC", "_P.Value", "_direction")))
    d
  }))
  metric_dt[, is_cell_cycle := gene_id %in% cell_cycle_ids]

  summarise_patterns <- function(x, label) {
    data.table(
      contrast_group = ct,
      gene_set = label,
      n_genes = nrow(x),
      scanning_up_bottleneck = x[scanning_direction == "Up", .N],
      scanning_down_less_bottleneck = x[scanning_direction == "Down", .N],
      ribosome_engagement_up = x[ribosome_engagement_direction == "Up", .N],
      ribosome_engagement_down = x[ribosome_engagement_direction == "Down", .N],
      collision_up_slow = x[collision_direction == "Up", .N],
      collision_down_less_slow = x[collision_direction == "Down", .N],
      productive_pattern_scanDown_riboUp_collDown = x[
        scanning_direction == "Down" &
          ribosome_engagement_direction == "Up" &
          collision_direction == "Down", .N],
      slowdown_pattern_scanUp_or_riboDown_or_collUp = x[
        scanning_direction == "Up" |
          ribosome_engagement_direction == "Down" |
          collision_direction == "Up", .N],
      riboUp_collisionDown = x[
        ribosome_engagement_direction == "Up" &
          collision_direction == "Down", .N],
      riboDown_collisionUp = x[
        ribosome_engagement_direction == "Down" &
          collision_direction == "Up", .N]
    )
  }

  rbindlist(list(
    summarise_patterns(metric_dt[is_cell_cycle == TRUE], "cell_cycle_GO_hit_genes"),
    summarise_patterns(metric_dt[is_cell_cycle == FALSE], "background_genes")
  ))
}), fill = TRUE)

fwrite(pathway_summary, file.path(out_dir, "cell_cycle_pathway_summary.csv"))
fwrite(cc_term_details, file.path(out_dir, "cell_cycle_term_details.csv"))
fwrite(hit_gene_summary, file.path(out_dir, "cell_cycle_hit_gene_summary.csv"))
fwrite(fraction_gene_summary, file.path(out_dir, "cell_cycle_gene_level_fraction_limma_summary.csv"))
fwrite(metric_gene_summary, file.path(out_dir, "cell_cycle_gene_level_metric_limma_summary.csv"))
fwrite(metric_pattern_summary, file.path(out_dir, "cell_cycle_metric_pattern_summary.csv"))

cat("\nCell-cycle GO:BP term recurrence:\n")
print(pathway_summary[, .(
  dataset, contrast_group, analysis_label,
  total_gobp_terms,
  cell_cycle_terms,
  pct_cell_cycle_terms = round(pct_cell_cycle_terms, 1),
  mostly_up = cell_cycle_mostly_up,
  mostly_down = cell_cycle_mostly_down,
  mixed = cell_cycle_mixed
)])

cat("\nUnique cell-cycle GO-hit genes represented in ORA terms:\n")
print(hit_gene_summary)

cat("\nGene-level direction for cell-cycle GO-hit genes in fraction-specific limma:\n")
print(fraction_gene_summary[, .(
  contrast_group, analysis,
  tested_cell_cycle_genes,
  sig_up,
  sig_down,
  sig_any,
  pct_sig_up = round(pct_sig_up, 1),
  pct_sig_down = round(pct_sig_down, 1),
  median_logFC_cell_cycle = round(median_logFC_cell_cycle, 3),
  median_logFC_background = round(median_logFC_background, 3),
  wilcox_p = signif(wilcox_cell_cycle_vs_background_p, 3)
)])

cat("\nGene-level direction for cell-cycle GO-hit genes in translation metrics:\n")
print(metric_gene_summary[, .(
  contrast_group, analysis,
  tested_cell_cycle_genes,
  sig_up,
  sig_down,
  sig_any,
  pct_sig_up = round(pct_sig_up, 1),
  pct_sig_down = round(pct_sig_down, 1),
  median_logFC_cell_cycle = round(median_logFC_cell_cycle, 3),
  median_logFC_background = round(median_logFC_background, 3),
  wilcox_p = signif(wilcox_cell_cycle_vs_background_p, 3)
)])

cat("\nTranslation-metric pattern counts for cell-cycle GO-hit genes vs background:\n")
print(metric_pattern_summary)

cat("\nCell-cycle term details:\n")
print(cc_term_details[, .(
  dataset, contrast_group, analysis_label, term_name,
  p_value = signif(p_value, 3),
  intersection_size,
  n_hit_up,
  n_hit_down,
  term_majority,
  mean_hit_log2FC = round(mean_hit_log2FC, 3)
)])

cat("\nSaved tables to: ", out_dir, "\n", sep = "")
