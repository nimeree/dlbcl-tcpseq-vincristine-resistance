# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- analysis_path()
limma_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05")
res_dir <- file.path(limma_dir, "Results")
out_dir <- file.path(limma_dir, "Multi_metric_integration", "LPXN_HNRNPD_followup_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_cut <- 0.05
lfc_cut <- 0.7
contrasts <- c("Resistance_baseline", "VCR_sensitive", "VCR_resistant", "Interaction")
metrics <- c(
  scanning = "scanning_score",
  ribosome_engagement = "ribosome_efficiency_score",
  protein_output = "protein_output_score",
  collision = "collision_score"
)

read_metric <- function(metric_key, contrast) {
  f <- file.path(res_dir, metrics[[metric_key]], paste0(contrast, "_limma_all_genes.csv"))
  d <- fread(f)
  if ("gene_id" %in% names(d) && !"gene_id_clean" %in% names(d)) {
    d[, gene_id_clean := sub("\\.\\d+$", "", gene_id)]
  }
  d[, direction := fifelse(
    P.Value < p_cut & logFC >= lfc_cut, "Up",
    fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS")
  )]
  x <- d[, .(gene_id_clean, gene_name, logFC, P.Value, adj.P.Val, direction)]
  setnames(
    x,
    c("logFC", "P.Value", "adj.P.Val", "direction"),
    paste(metric_key, c("logFC", "P", "FDR", "dir"), sep = "_")
  )
  x
}

wide_for_contrast <- function(contrast) {
  Reduce(
    function(x, y) merge(x, y, by = c("gene_id_clean", "gene_name"), all = FALSE),
    lapply(names(metrics), read_metric, contrast = contrast)
  )
}

interaction <- wide_for_contrast("Interaction")
productive <- interaction[ribosome_engagement_dir == "Up" & collision_dir == "Down"]
productive <- productive[order(gene_name)]

write.csv(productive, file.path(out_dir, "interaction_productive_translation_38gene_set.csv"), row.names = FALSE)

target_genes <- c(
  "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2",
  "ATF2", "ATF3", "ATF4", "ELK1", "EGR1", "EGR2", "EGR3",
  "DUSP1", "DUSP2", "DUSP4", "DUSP5", "DUSP6", "DUSP8", "DUSP10",
  "GADD45A", "GADD45B", "GADD45G",
  "MAPKAPK2", "MAPKAPK3", "HSPB1", "ZFP36", "CREB1",
  "DDIT3", "PPP1R15A", "NFKBIA", "IER3",
  "HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "DNAJB1",
  "IL6", "CXCL8", "PTGS2", "TNF"
)

run_target_test <- function(metric_key) {
  logfc_col <- paste0(metric_key, "_logFC")
  p_col <- paste0(metric_key, "_P")
  dir_col <- paste0(metric_key, "_dir")
  dt <- copy(interaction)
  dt[, in_target_set := gene_name %in% target_genes]
  present <- sort(unique(dt[in_target_set == TRUE, gene_name]))
  test_dt <- dt[!is.na(get(logfc_col))]
  wt <- wilcox.test(
    test_dt[[logfc_col]][test_dt$in_target_set],
    test_dt[[logfc_col]][!test_dt$in_target_set],
    exact = FALSE
  )
  data.table(
    metric = metric_key,
    n_target_input = length(target_genes),
    n_target_present = length(present),
    target_present = paste(present, collapse = ";"),
    median_target_logFC = median(test_dt[[logfc_col]][test_dt$in_target_set], na.rm = TRUE),
    median_background_logFC = median(test_dt[[logfc_col]][!test_dt$in_target_set], na.rm = TRUE),
    mean_target_logFC = mean(test_dt[[logfc_col]][test_dt$in_target_set], na.rm = TRUE),
    mean_background_logFC = mean(test_dt[[logfc_col]][!test_dt$in_target_set], na.rm = TRUE),
    wilcoxon_p = wt$p.value,
    n_target_up = sum(test_dt$in_target_set & test_dt[[dir_col]] == "Up", na.rm = TRUE),
    n_target_down = sum(test_dt$in_target_set & test_dt[[dir_col]] == "Down", na.rm = TRUE),
    n_target_ns = sum(test_dt$in_target_set & test_dt[[dir_col]] == "NS", na.rm = TRUE)
  )
}

jnk_tests <- rbindlist(lapply(c("ribosome_engagement", "collision", "scanning", "protein_output"), run_target_test))

jnk_gene_profile <- interaction[gene_name %in% target_genes]
jnk_gene_profile <- jnk_gene_profile[order(gene_name)]

write.csv(jnk_tests, file.path(out_dir, "JNK_p38_target_interaction_wilcoxon_tests.csv"), row.names = FALSE)
write.csv(jnk_gene_profile, file.path(out_dir, "JNK_p38_target_interaction_gene_profiles.csv"), row.names = FALSE)

gene_metric_profile <- rbindlist(lapply(contrasts, function(contrast) {
  w <- wide_for_contrast(contrast)
  x <- w[gene_name %in% c("LPXN", "HNRNPD")]
  x[, contrast := contrast]
  x
}), fill = TRUE)

setcolorder(gene_metric_profile, c("contrast", setdiff(names(gene_metric_profile), "contrast")))
write.csv(gene_metric_profile, file.path(out_dir, "LPXN_HNRNPD_all_contrast_metric_profile.csv"), row.names = FALSE)

cat("\nProductive translation set: ribosome engagement Up + collision Down in Interaction\n")
cat("Count:", nrow(productive), "\n")
cat("Genes:", paste(productive$gene_name, collapse = ", "), "\n")

cat("\nJNK/p38 target-set Wilcoxon tests vs all other genes, Interaction logFC\n")
print(jnk_tests)

cat("\nJNK/p38 target genes present in limma matrix\n")
cat(jnk_tests[metric == "ribosome_engagement", target_present], "\n")

cat("\nLPXN and HNRNPD metric profile across contrasts\n")
profile_long <- melt(
  gene_metric_profile,
  id.vars = c("contrast", "gene_id_clean", "gene_name"),
  measure.vars = patterns("_logFC$", "_P$", "_FDR$", "_dir$"),
  variable.name = "metric_field"
)
print(gene_metric_profile[
  order(gene_name, match(contrast, contrasts)),
  .(
    gene_name,
    contrast,
    scanning_logFC, scanning_P, scanning_dir,
    ribosome_engagement_logFC, ribosome_engagement_P, ribosome_engagement_dir,
    protein_output_logFC, protein_output_P, protein_output_dir,
    collision_logFC, collision_P, collision_dir
  )
])

cat("\nOutput directory:\n")
cat(out_dir, "\n")
