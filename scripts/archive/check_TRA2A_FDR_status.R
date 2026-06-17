# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)

base <- analysis_path()
files <- list.files(base, recursive = TRUE, full.names = TRUE, pattern = "[.]csv$")

wanted_cols <- c(
  "metric", "contrast", "gene_id", "gene_name", "gene", "gene_symbol",
  "logFC", "log2FoldChange", "AveExpr", "t", "P.Value", "adj.P.Val",
  "pvalue", "padj", "significant", "sig_p05_lfc0.7", "direction",
  "regulation", "row_id", "transcript_name", "chi_square_pvalue",
  "max_abs_proportion_shift"
)

hits <- rbindlist(lapply(files, function(f) {
  header <- tryCatch(fread(f, nrows = 0), error = function(e) NULL)
  if (is.null(header)) return(NULL)
  gene_cols <- intersect(c("gene_name", "gene", "gene_symbol"), names(header))
  if (!length(gene_cols)) return(NULL)

  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  idx <- rep(FALSE, nrow(d))
  for (gc in gene_cols) idx <- idx | d[[gc]] == "TRA2A"
  hit <- d[idx]
  if (!nrow(hit)) return(NULL)

  cols <- intersect(wanted_cols, names(hit))
  out <- hit[, ..cols]
  out[, file := sub(base, "", f, fixed = TRUE)]
  out
}), fill = TRUE)

if (!nrow(hits)) {
  cat("No TRA2A rows found.\n")
} else {
  print(hits, nrows = 250, width = 220)
}
