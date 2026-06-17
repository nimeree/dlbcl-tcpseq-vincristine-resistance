# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(gprofiler2)
})

base_dir <- analysis_path()
count_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
metric_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05", "Results")

p_cut <- 0.05
lfc_cut <- 0.7

direction_call <- function(logfc, pval) {
  fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc >= lfc_cut, "Up",
          fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc <= -lfc_cut, "Down", "NS"))
}

metric <- fread(file.path(metric_dir, "collision_score", "Interaction_limma_all_genes.csv"))
count <- fread(file.path(count_dir, "Fraction_DS", "Interaction_DS_psite_limma_all_genes.csv"))

metric[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
count[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
if (!"gene_name" %in% names(metric)) metric[, gene_name := gene_id_clean]
if (!"gene_name" %in% names(count)) count[, gene_name := gene_id_clean]

metric[, metric_direction := direction_call(logFC, P.Value)]
count[, count_direction := direction_call(logFC, P.Value)]

joined <- merge(
  metric[, .(
    gene_id_clean,
    gene_name_metric = gene_name,
    metric_logFC = logFC,
    metric_P = P.Value,
    metric_direction
  )],
  count[, .(
    gene_id_clean,
    gene_name_count = gene_name,
    DS_logFC = logFC,
    DS_P = P.Value,
    count_direction
  )],
  by = "gene_id_clean"
)
joined[, gene_name := fifelse(!is.na(gene_name_metric) & nzchar(gene_name_metric), gene_name_metric, gene_name_count)]

conv <- joined[metric_direction == "Down" & count_direction == "Down"]
cat("Convergent Interaction collision/DS Down genes:", nrow(conv), "\n\n")

res <- gost(
  query = unique(conv$gene_id_clean),
  organism = "hsapiens",
  sources = c("GO:BP"),
  correction_method = "g_SCS",
  domain_scope = "annotated",
  user_threshold = 0.05,
  evcodes = TRUE
)

dt <- as.data.table(res$result)
target_patterns <- c(
  "mRNA transport",
  "mRNA export from nucleus",
  "nucleic acid transport",
  "RNA transport",
  "establishment of RNA localization",
  "nuclear export"
)
pattern <- paste(target_patterns, collapse = "|")
terms <- dt[
  significant == TRUE &
    term_size <= 500 &
    intersection_size >= 5 &
    grepl(pattern, term_name, ignore.case = TRUE)
]
setorder(terms, p_value)

cat("Target mRNA/RNA transport-export GO terms:\n")
print(terms[, .(term_id, term_name, p_value, term_size, intersection_size, intersection)])

cat("\nHit genes per term, mapped back to symbols and logFCs:\n")
for (i in seq_len(nrow(terms))) {
  ids <- unlist(strsplit(terms$intersection[i], ","))
  ids <- sub("\\.\\d+$", "", ids)
  hits <- conv[gene_id_clean %in% ids, .(
    gene_name,
    gene_id_clean,
    collision_metric_logFC = round(metric_logFC, 3),
    collision_metric_P = signif(metric_P, 3),
    DS_psite_logFC = round(DS_logFC, 3),
    DS_psite_P = signif(DS_P, 3)
  )]
  setorder(hits, gene_name)
  cat("\n", terms$term_name[i], " (", terms$term_id[i], "), hits=", nrow(hits), "\n", sep = "")
  print(hits)
}

all_ids <- unique(unlist(strsplit(terms$intersection, ",")))
all_ids <- sub("\\.\\d+$", "", all_ids)
all_hits <- conv[gene_id_clean %in% all_ids, .(
  gene_name,
  gene_id_clean,
  collision_metric_logFC = round(metric_logFC, 3),
  collision_metric_P = signif(metric_P, 3),
  DS_psite_logFC = round(DS_logFC, 3),
  DS_psite_P = signif(DS_P, 3)
)]
setorder(all_hits, gene_name)

cat("\nUnique genes across these mRNA/RNA transport-export terms:", nrow(all_hits), "\n")
print(all_hits)
