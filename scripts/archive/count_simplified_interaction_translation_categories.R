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
  collision = "collision_score"
)

read_res <- function(k) {
  d <- fread(file.path(base, metrics[[k]], paste0(contrast, "_limma_all_genes.csv")))
  d[, direction := fifelse(
    P.Value < p_cut & logFC >= lfc_cut,
    "Up",
    fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS")
  )]
  x <- d[, .(gene_id_clean, gene_name, logFC, P.Value, adj.P.Val, direction)]
  setnames(
    x,
    c("logFC", "P.Value", "adj.P.Val", "direction"),
    paste(k, c("logFC", "P", "FDR", "dir"), sep = "_")
  )
  x
}

wide <- Reduce(
  function(x, y) merge(x, y, by = c("gene_id_clean", "gene_name"), all = FALSE),
  lapply(names(metrics), read_res)
)

patterns <- list(
  "Scanning bottleneck (Scanning Up)" = wide[scanning_dir == "Up"],
  "Collision stress (Collision Up)" = wide[collision_dir == "Up"],
  "Productive/prioritised translation (Ribosome Up + Collision Down)" =
    wide[ribosome_dir == "Up" & collision_dir == "Down"],
  "Full-cycle prioritisation subset (Scanning Up + Ribosome Up + Collision Down)" =
    wide[scanning_dir == "Up" & ribosome_dir == "Up" & collision_dir == "Down"]
)

for (nm in names(patterns)) {
  cat("\n", nm, "\n", sep = "")
  cat("Count:", nrow(patterns[[nm]]), "\n")
  if (nrow(patterns[[nm]]) > 0) {
    genes <- sort(patterns[[nm]]$gene_name)
    cat("Genes:", paste(genes, collapse = ", "), "\n")
  }
}

cat("\nOverlap notes\n")
cat("Scanning Up AND Collision Up:", nrow(wide[scanning_dir == "Up" & collision_dir == "Up"]), "\n")
cat("Scanning Up AND Ribosome Up:", nrow(wide[scanning_dir == "Up" & ribosome_dir == "Up"]), "\n")
cat("Collision Up AND Ribosome Down:", nrow(wide[collision_dir == "Up" & ribosome_dir == "Down"]), "\n")
