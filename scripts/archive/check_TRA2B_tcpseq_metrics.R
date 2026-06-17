# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

base <- analysis_path()

metrics <- c(
  ribosome_engagement = "ribosome_efficiency_score",
  protein_output = "protein_output_score",
  scanning = "scanning_score",
  collision = "collision_score"
)
contrasts <- c("Resistance_baseline", "Interaction", "VCR_sensitive", "VCR_resistant")

for (metric_label in names(metrics)) {
  metric <- metrics[[metric_label]]
  cat("\n====", metric_label, "(", metric, ") ====\n")
  for (con in contrasts) {
    f <- file.path(
      base,
      "Limma_translation_metrics_lfc0.7_rawP0.05",
      "Results",
      metric,
      paste0(con, "_limma_all_genes.csv")
    )
    d <- fread(f)
    x <- d[gene_name == "TRA2B" | gene_id_clean == "ENSG00000136527"]
    if (nrow(x)) {
      x[, contrast := con]
      print(x[, .(
        contrast, gene_id_clean, gene_name, logFC, AveExpr, t,
        P.Value, adj.P.Val, significant_rawP0.05_lfc0.7, direction
      )])
    } else {
      cat(con, ": not found\n")
    }
  }
}

cat("\n==== gene-level metric values ====\n")
gf <- file.path(
  base,
  "Translation_indexes_fixed",
  "Gene_Level_Clean",
  "gene_level_clean_translation_metrics_all_samples.csv"
)
g <- fread(gf)
x <- g[gene_name == "TRA2B" | gene_id_clean == "ENSG00000136527"]
print(x[, .(
  sample_label, cell_line_from_sample, treatment, replicate,
  n_core, n_utr5, rs_core_cpm, baseline_cpm_line,
  ribosome_efficiency_score, protein_output_score,
  collision_score, scanning_score
)][order(cell_line_from_sample, treatment, replicate)])

cat("\nSummary by condition\n")
print(x[, .(
  n = .N,
  mean_n_core = mean(n_core),
  mean_rs_core_cpm = mean(rs_core_cpm),
  mean_baseline_cpm = mean(baseline_cpm_line),
  mean_ribo_engagement = mean(ribosome_efficiency_score),
  mean_protein_output = mean(protein_output_score),
  mean_collision = mean(collision_score, na.rm = TRUE),
  mean_scanning = mean(scanning_score, na.rm = TRUE)
), by = .(condition)][order(condition)])
