# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(gprofiler2)
})

limma_dir <- analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05")
res_dir <- file.path(limma_dir, "Results")
follow_dir <- file.path(limma_dir, "Multi_metric_integration", "LPXN_HNRNPD_followup_tables")
out_file <- file.path(follow_dir, "productive_translation_38gene_gprofiler_REAC_GOBP_GOCC.csv")
dir.create(follow_dir, recursive = TRUE, showWarnings = FALSE)

read_metric <- function(metric) {
  d <- fread(file.path(res_dir, metric, "Interaction_limma_all_genes.csv"))
  d <- d[!is.na(gene_name) & gene_name != ""]
  d[, direction := fifelse(
    P.Value < 0.05 & logFC >= 0.7, "Up",
    fifelse(P.Value < 0.05 & logFC <= -0.7, "Down", "NS")
  )]
  d[, .(gene_id_clean, gene_name, logFC, P.Value, adj.P.Val, direction)]
}

rib <- read_metric("ribosome_efficiency_score")
setnames(rib, c("logFC", "P.Value", "adj.P.Val", "direction"), c("ribosome_logFC", "ribosome_P", "ribosome_FDR", "ribosome_dir"))
collision <- read_metric("collision_score")
setnames(collision, c("logFC", "P.Value", "adj.P.Val", "direction"), c("collision_logFC", "collision_P", "collision_FDR", "collision_dir"))

productive <- merge(rib, collision, by = c("gene_id_clean", "gene_name"))
productive <- productive[ribosome_dir == "Up" & collision_dir == "Down"]

# The old "38-gene productive" bar was the productive ribosome/collision
# set excluding the three-gene full-cycle subset highlighted separately.
full_cycle_old <- c("MAPKBP1", "SEC24C", "TRA2A")
query_dt <- productive[!gene_name %in% full_cycle_old]
query <- unique(query_dt$gene_id_clean)

cat("Input genes:", nrow(query_dt), "\n")
cat("Symbols:", paste(sort(query_dt$gene_name), collapse = ", "), "\n")
write.csv(query_dt, file.path(follow_dir, "productive_translation_38gene_query_used_for_gprofiler.csv"), row.names = FALSE)

res <- gost(
  query = query,
  organism = "hsapiens",
  sources = c("REAC", "GO:BP", "GO:CC"),
  correction_method = "g_SCS",
  domain_scope = "annotated",
  user_threshold = 0.05,
  evcodes = TRUE
)

if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
  out <- data.table()
  cat("\nNo significant g:Profiler terms returned.\n")
} else {
  out <- as.data.table(res$result)
  out <- out[significant == TRUE]
  out <- out[order(p_value, -intersection_size)]
  list_cols <- names(out)[vapply(out, is.list, logical(1))]
  for (lc in list_cols) {
    out[, (lc) := vapply(get(lc), function(x) paste(unlist(x), collapse = ";"), character(1))]
  }
  write.csv(out, out_file, row.names = FALSE)

  filtered <- out[term_size <= 500 & intersection_size >= 3]
  cat("\nSignificant terms after term_size <= 500 and hit genes >= 3:\n")
  print(filtered[, .(source, term_name, p_value, term_size, intersection_size, intersection)])
  cat("\nAll significant terms saved to:\n", out_file, "\n", sep = "")
}
