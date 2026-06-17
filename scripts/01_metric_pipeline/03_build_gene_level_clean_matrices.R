# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

INFILE <- analysis_path("Translation_indexes_fixed", "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")
OUT_DIR <- analysis_path("Translation_indexes_fixed", "Gene_Level_Clean")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

mode_first <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

parse_sample_meta <- function(sample) {
  data.table(sample = sample)[, `:=`(
    cell_line_from_sample = fifelse(grepl("^SU8R", sample, ignore.case = TRUE), "Resistant", "Sensitive"),
    treatment = fifelse(grepl("-Vin-", sample, ignore.case = TRUE), "VCR", "DMSO"),
    replicate = fifelse(grepl("Rep1", sample, ignore.case = TRUE), "Rep1",
      fifelse(grepl("Rep2", sample, ignore.case = TRUE), "Rep2", NA_character_))
  )][, condition := paste(cell_line_from_sample, treatment, sep = "_")]
}

dt <- fread(INFILE)
dt <- dt[!is.na(gene_id_clean) & gene_id_clean != ""]

id_cols <- c("gene_id_clean", "gene_name", "sample", "cell_line")
numeric_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
metric_cols <- setdiff(numeric_cols, character())

gene_df <- dt[, c(
  .(gene_name = mode_first(gene_name)),
  lapply(.SD, function(x) as.numeric(median(x, na.rm = TRUE)))
), by = .(gene_id_clean, sample, cell_line), .SDcols = metric_cols]

meta <- parse_sample_meta(unique(gene_df$sample))
gene_df <- merge(gene_df, meta, by = "sample", all.x = TRUE)
gene_df[, sample_label := paste(condition, replicate, sep = "_")]

# Supervisor-suggested cleaning.
gene_df <- gene_df[baseline_cpm_line > 0]
gene_df[scanning_score > 10, scanning_score := NA_real_]
gene_df[total_translation_rate_proxy_stable > 1e6, total_translation_rate_proxy_stable := NA_real_]

all_samples <- sort(unique(gene_df$sample))
n_samples <- length(all_samples)

complete_gene_set <- function(required_cols) {
  gene_df[, .(
    n_samples_present = uniqueN(sample),
    complete_all_required = all(vapply(required_cols, function(cc) {
      cc %in% names(.SD) && all(is.finite(.SD[[cc]]))
    }, logical(1)))
  ), by = gene_id_clean, .SDcols = required_cols][
    n_samples_present == n_samples & complete_all_required == TRUE,
    gene_id_clean
  ]
}

te_protein_genes <- complete_gene_set(c("ribosome_efficiency_score", "protein_output_score"))
collision_genes <- complete_gene_set("collision_score")

gene_df[, complete_te_protein_set := gene_id_clean %chin% te_protein_genes]
gene_df[, complete_collision_set := gene_id_clean %chin% collision_genes]

fwrite(gene_df, file.path(OUT_DIR, "gene_level_clean_translation_metrics_all_samples.csv"))
fwrite(gene_df[complete_te_protein_set == TRUE], file.path(OUT_DIR, "gene_level_clean_TE_protein_output_complete_8_samples.csv"))
fwrite(gene_df[complete_collision_set == TRUE], file.path(OUT_DIR, "gene_level_clean_collision_complete_8_samples.csv"))

make_wide <- function(metric, gene_ids, file_name) {
  x <- gene_df[gene_id_clean %chin% gene_ids & is.finite(get(metric)),
    .(value = median(get(metric), na.rm = TRUE)),
    by = .(gene_id_clean, gene_name, sample_label)
  ]
  wide <- dcast(x, gene_id_clean + gene_name ~ sample_label, value.var = "value")
  fwrite(wide, file.path(OUT_DIR, file_name))
  invisible(wide)
}

make_wide("ribosome_efficiency_score", te_protein_genes, "matrix_ribosome_efficiency_score_complete_8_samples.csv")
make_wide("protein_output_score", te_protein_genes, "matrix_protein_output_score_complete_8_samples.csv")
make_wide("collision_score", collision_genes, "matrix_collision_score_complete_8_samples.csv")

summary_dt <- rbindlist(list(
  data.table(metric_set = "all_gene_sample_rows_after_baseline_filter", genes = uniqueN(gene_df$gene_id_clean), gene_sample_rows = nrow(gene_df)),
  data.table(metric_set = "complete_TE_protein_output_8_samples", genes = length(te_protein_genes), gene_sample_rows = nrow(gene_df[complete_te_protein_set == TRUE])),
  data.table(metric_set = "complete_collision_score_8_samples", genes = length(collision_genes), gene_sample_rows = nrow(gene_df[complete_collision_set == TRUE])),
  data.table(metric_set = "scanning_score_gt_10_removed", genes = NA_integer_, gene_sample_rows = dt[, .(
    gene_id_clean, sample, scan = median(scanning_score, na.rm = TRUE)
  ), by = .(gene_id_clean, sample)][scan > 10, .N]),
  data.table(metric_set = "total_translation_rate_proxy_stable_gt_1e6_removed", genes = NA_integer_, gene_sample_rows = dt[, .(
    gene_id_clean, sample, proxy = median(total_translation_rate_proxy_stable, na.rm = TRUE)
  ), by = .(gene_id_clean, sample)][proxy > 1e6, .N])
), fill = TRUE)
fwrite(summary_dt, file.path(OUT_DIR, "gene_level_clean_matrix_summary.csv"))

cat("\nGene-level clean matrix summary\n")
cat("===============================\n")
print(summary_dt)
cat("\nOutputs written to:", OUT_DIR, "\n")
