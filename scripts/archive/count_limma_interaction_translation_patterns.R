# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

base <- analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05", "Results")
contrast <- "Interaction"
p_cut <- 0.05
lfc_cut <- 0.7

metrics <- c(
  scanning = "scanning_score",
  ribosome = "ribosome_efficiency_score",
  protein = "protein_output_score",
  collision = "collision_score"
)

read_res <- function(metric_key) {
  metric <- metrics[[metric_key]]
  d <- fread(file.path(base, metric, paste0(contrast, "_limma_all_genes.csv")))
  d[, direction := fifelse(
    P.Value < p_cut & logFC >= lfc_cut,
    "Up",
    fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS")
  )]
  keep <- d[, .(gene_id_clean, gene_name, logFC, P.Value, adj.P.Val, direction)]
  setnames(
    keep,
    c("logFC", "P.Value", "adj.P.Val", "direction"),
    paste(metric_key, c("logFC", "P", "FDR", "dir"), sep = "_")
  )
  keep
}

wide <- Reduce(
  function(x, y) merge(x, y, by = c("gene_id_clean", "gene_name"), all = FALSE),
  lapply(names(metrics), read_res)
)

patterns <- list(
  "Full-cycle translation" = wide[
    scanning_dir == "Up" &
      ribosome_dir == "Up" &
      protein_dir == "Up" &
      collision_dir == "Down"
  ],
  "Collision stress" = wide[
    ribosome_dir == "Down" &
      collision_dir == "Up" &
      protein_dir == "Down"
  ],
  "Scanning bottleneck" = wide[
    scanning_dir == "Up" &
      ribosome_dir == "Down" &
      protein_dir == "Down"
  ]
)

for (nm in names(patterns)) {
  cat("\n", nm, "\n", sep = "")
  cat("Count:", nrow(patterns[[nm]]), "\n")
  if (nrow(patterns[[nm]]) > 0) {
    cat("Genes:", paste(sort(patterns[[nm]]$gene_name), collapse = ", "), "\n")
  }
}

details <- rbindlist(lapply(names(patterns), function(nm) {
  x <- copy(patterns[[nm]])
  if (!nrow(x)) return(NULL)
  x[, pattern := nm]
  x
}), fill = TRUE)

cat("\nDetailed table\n")
if (nrow(details)) {
  setcolorder(details, c("pattern", "gene_id_clean", "gene_name"))
  print(details[order(pattern, gene_name), .(
    pattern, gene_name,
    scanning_logFC, ribosome_logFC, protein_logFC, collision_logFC,
    scanning_P, ribosome_P, protein_P, collision_P,
    scanning_FDR, ribosome_FDR, protein_FDR, collision_FDR
  )])
} else {
  cat("No hits in any pattern\n")
}
