# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

f <- analysis_path("Translation_indexes_fixed", "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")
backup <- sub("[.]csv$", "_pre_rsrate_stablemask_backup.csv", f)

if (!file.exists(backup)) {
  file.copy(f, backup)
}

dt <- fread(f)
req <- c("n_cds", "n_core", "rs_rate", "baseline_cpm_line")
missing <- setdiff(req, names(dt))
if (length(missing)) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}
numeric_needed <- intersect(c(
  req,
  "initiation_rate_index", "elongation_rate_index", "total_translation_rate_proxy",
  "initiation_rate_index_stable", "elongation_rate_index_stable", "total_translation_rate_proxy_stable"
), names(dt))
for (cc in numeric_needed) {
  dt[, (cc) := as.numeric(get(cc))]
}

dt[, stable_mask := n_cds >= 20 & n_core >= 20]

dt[, ribosome_efficiency_score := NA_real_]
dt[
  !is.na(baseline_cpm_line) & stable_mask == TRUE,
  ribosome_efficiency_score := log2((rs_rate + 1) / (baseline_cpm_line + 1))
]

dt[, protein_output_score := NA_real_]
dt[
  !is.na(baseline_cpm_line) & stable_mask == TRUE,
  protein_output_score := log2((baseline_cpm_line + 1) * (rs_rate + 1))
]

eps <- 1e-6
if ("initiation_rate_index" %in% names(dt)) dt[, initiation_rate_TEproxy := initiation_rate_index * pmax(ribosome_efficiency_score, eps)]
if ("elongation_rate_index" %in% names(dt)) dt[, elongation_rate_TEproxy := elongation_rate_index * pmax(ribosome_efficiency_score, eps)]
if ("total_translation_rate_proxy" %in% names(dt)) dt[, total_translation_rate_TEproxy := total_translation_rate_proxy * pmax(ribosome_efficiency_score, eps)]
if ("initiation_rate_index_stable" %in% names(dt)) dt[, initiation_rate_TEproxy_stable := initiation_rate_index_stable * pmax(ribosome_efficiency_score, eps)]
if ("elongation_rate_index_stable" %in% names(dt)) dt[, elongation_rate_TEproxy_stable := elongation_rate_index_stable * pmax(ribosome_efficiency_score, eps)]
if ("total_translation_rate_proxy_stable" %in% names(dt)) dt[, total_translation_rate_TEproxy_stable := total_translation_rate_proxy_stable * pmax(ribosome_efficiency_score, eps)]

fwrite(dt, f)

cat("Updated:", f, "\n")
cat("Backup:", backup, "\n")
cat("Rows:", nrow(dt), "\n")
cat("Stable rows:", sum(dt$stable_mask, na.rm = TRUE), "\n")
cat("Finite ribosome_efficiency_score:", sum(is.finite(dt$ribosome_efficiency_score)), "\n")
cat("Finite protein_output_score:", sum(is.finite(dt$protein_output_score)), "\n")
cat("Genes with at least one finite ribosome_efficiency_score:", uniqueN(dt[is.finite(ribosome_efficiency_score), gene_id_clean]), "\n")
