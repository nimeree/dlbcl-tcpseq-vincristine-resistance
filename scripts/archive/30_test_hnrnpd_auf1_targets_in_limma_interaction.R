# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05")
results_dir <- file.path(base_dir, "Results")
out_dir <- file.path(base_dir, "Multi_metric_integration", "HNRNPD_AUF1_targets")
target_file <- file.path(out_dir, "Yoon2014_AUF1_PARCLIP_gene_level_targets.csv")

p_cut <- 0.05
lfc_cut <- 0.7

metrics <- c(
  scanning = "scanning_score",
  ribosome_engagement = "ribosome_efficiency_score",
  protein_output = "protein_output_score",
  collision = "collision_score"
)

read_metric <- function(metric_key) {
  f <- file.path(results_dir, metrics[[metric_key]], "Interaction_limma_all_genes.csv")
  d <- fread(f)
  if ("gene_id" %in% names(d) && !"gene_id_clean" %in% names(d)) {
    d[, gene_id_clean := sub("\\.\\d+$", "", gene_id)]
  }
  d <- d[!is.na(gene_name) & gene_name != ""]
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

wide <- Reduce(
  function(x, y) merge(x, y, by = c("gene_id_clean", "gene_name"), all = FALSE),
  lapply(names(metrics), read_metric)
)

targets <- fread(target_file)
targets[, gene_name := as.character(gene_name)]
targets[, AUF1_any := total_auf1_target_sites > 0]

dat <- merge(
  wide,
  targets[, .(
    gene_name,
    AUF1_any,
    n_target_transcripts,
    total_auf1_target_sites,
    total_auf1_t_to_c,
    total_hur_target_sites,
    total_hur_t_to_c
  )],
  by = "gene_name",
  all.x = TRUE
)
dat[is.na(AUF1_any), `:=`(
  AUF1_any = FALSE,
  n_target_transcripts = 0,
  total_auf1_target_sites = 0,
  total_auf1_t_to_c = 0,
  total_hur_target_sites = 0,
  total_hur_t_to_c = 0
)]

# HNRNPD itself is the regulator, not a downstream target for this test.
dat[, test_universe := gene_name != "HNRNPD"]
target_site_cut <- quantile(dat[test_universe == TRUE & AUF1_any == TRUE, total_auf1_target_sites], 0.75, na.rm = TRUE)
tc_cut <- quantile(dat[test_universe == TRUE & AUF1_any == TRUE, total_auf1_t_to_c], 0.75, na.rm = TRUE)
dat[, AUF1_high_sites := AUF1_any & total_auf1_target_sites >= target_site_cut]
dat[, AUF1_high_tc := AUF1_any & total_auf1_t_to_c >= tc_cut]

run_test <- function(metric_key, set_col, alternative) {
  logfc_col <- paste0(metric_key, "_logFC")
  dir_col <- paste0(metric_key, "_dir")
  d <- dat[test_universe == TRUE]
  in_set <- d[[set_col]] == TRUE
  x <- d[[logfc_col]][in_set]
  y <- d[[logfc_col]][!in_set]
  w_two <- wilcox.test(x, y, exact = FALSE, alternative = "two.sided")
  w_one <- wilcox.test(x, y, exact = FALSE, alternative = alternative)
  data.table(
    target_set = set_col,
    metric = metric_key,
    expected_direction = fifelse(alternative == "less", "lower in AUF1 targets", "higher in AUF1 targets"),
    n_targets_present = sum(in_set),
    n_background = sum(!in_set),
    median_target_logFC = median(x, na.rm = TRUE),
    median_background_logFC = median(y, na.rm = TRUE),
    mean_target_logFC = mean(x, na.rm = TRUE),
    mean_background_logFC = mean(y, na.rm = TRUE),
    wilcoxon_two_sided_p = w_two$p.value,
    wilcoxon_expected_direction_p = w_one$p.value,
    n_target_up = sum(in_set & d[[dir_col]] == "Up", na.rm = TRUE),
    n_target_down = sum(in_set & d[[dir_col]] == "Down", na.rm = TRUE),
    n_target_ns = sum(in_set & d[[dir_col]] == "NS", na.rm = TRUE)
  )
}

tests <- rbindlist(list(
  run_test("ribosome_engagement", "AUF1_any", "less"),
  run_test("protein_output", "AUF1_any", "less"),
  run_test("collision", "AUF1_any", "greater"),
  run_test("scanning", "AUF1_any", "greater"),
  run_test("ribosome_engagement", "AUF1_high_sites", "less"),
  run_test("protein_output", "AUF1_high_sites", "less"),
  run_test("collision", "AUF1_high_sites", "greater"),
  run_test("scanning", "AUF1_high_sites", "greater"),
  run_test("ribosome_engagement", "AUF1_high_tc", "less"),
  run_test("protein_output", "AUF1_high_tc", "less"),
  run_test("collision", "AUF1_high_tc", "greater"),
  run_test("scanning", "AUF1_high_tc", "greater")
), fill = TRUE)
tests[, BH_expected_direction_p := p.adjust(wilcoxon_expected_direction_p, method = "BH")]
tests[, BH_two_sided_p := p.adjust(wilcoxon_two_sided_p, method = "BH")]

gene_profiles <- dat[AUF1_any == TRUE & test_universe == TRUE][order(-total_auf1_target_sites)]

write.csv(dat, file.path(out_dir, "AUF1_target_annotation_interaction_limma_all_genes.csv"), row.names = FALSE)
write.csv(gene_profiles, file.path(out_dir, "AUF1_targets_present_in_interaction_limma_gene_profiles.csv"), row.names = FALSE)
write.csv(tests, file.path(out_dir, "AUF1_target_interaction_metric_wilcoxon_tests.csv"), row.names = FALSE)

cat("\nAUF1/HNRNPD target source: Yoon et al. 2014 AUF1 PAR-CLIP Supplementary Table 1\n")
cat("ENCODE released HNRNPD eCLIP experiments found: 0\n")
cat("Total AUF1 PAR-CLIP target genes in source:", targets[AUF1_any == TRUE, .N], "\n")
cat("Targets present in complete limma interaction metric universe, excluding HNRNPD:", dat[test_universe == TRUE & AUF1_any == TRUE, .N], "\n")
cat("Background genes:", dat[test_universe == TRUE & AUF1_any == FALSE, .N], "\n")
cat("High-site threshold, top quartile target-site count:", target_site_cut, "\n")
cat("High-T-to-C threshold, top quartile T-to-C count:", tc_cut, "\n")

cat("\nWilcoxon tests vs background:\n")
print(tests[order(target_set, metric)])

cat("\nDirection counts for all AUF1 targets:\n")
print(tests[target_set == "AUF1_any", .(
  metric, n_targets_present, n_target_up, n_target_down, n_target_ns,
  median_target_logFC, median_background_logFC,
  wilcoxon_expected_direction_p, BH_expected_direction_p
)])

cat("\nTop AUF1 target genes present in limma by target-site count:\n")
print(head(gene_profiles[, .(
  gene_name,
  total_auf1_target_sites,
  total_auf1_t_to_c,
  ribosome_engagement_logFC,
  ribosome_engagement_P,
  ribosome_engagement_dir,
  collision_logFC,
  collision_P,
  collision_dir
)], 30))

cat("\nOutput directory:\n")
cat(out_dir, "\n")
