# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- analysis_path()

count_path <- file.path(
  base_dir,
  "Psite_fraction_limma_lfc0.7_rawP0.05",
  "Pathway_gProfiler_Clean_REAC_GOBP_GOCC_true_interaction_lfc0.7_rawP0.05",
  "Tables",
  "gprofiler_clean_filtered_terms_termSizeLE500_hitGE5.csv"
)
metric_path <- file.path(
  base_dir,
  "Limma_translation_metrics_lfc0.7_rawP0.05",
  "Pathway_gProfiler_Clean_REAC_GOBP_GOCC_lfc0.7_rawP0.05",
  "Tables",
  "gprofiler_clean_filtered_terms_termSizeLE500_hitGE5.csv"
)
dseq_path <- file.path(
  base_dir,
  "Pathway_gProfiler_Clean_REAC_GOBP_GOCC_lfc0.7_p0.05",
  "Tables",
  "gprofiler_clean_filtered_terms_termSizeLE500_hitGE5.csv"
)

out_dir <- file.path(
  base_dir,
  "Psite_fraction_limma_lfc0.7_rawP0.05",
  "Pathway_concordance_with_metric_and_DESeq"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

count <- fread(count_path)
metric <- fread(metric_path)
dseq <- fread(dseq_path)

term_key_cols <- c("source", "term_id")

metric_map <- data.table(
  fraction = c("SSU", "RS", "DS"),
  metric_fraction = c("scanning_score", "ribosome_efficiency_score", "collision_score"),
  label = c("SSU vs scanning", "RS vs ribosome engagement", "DS vs collision")
)

compare_sets <- function(a, b, by_cols = term_key_cols) {
  a_terms <- unique(a[, ..by_cols])
  b_terms <- unique(b[, ..by_cols])
  a_terms[, in_a := TRUE]
  b_terms[, in_b := TRUE]
  joined <- merge(a_terms, b_terms, by = by_cols, all = TRUE)
  n_a <- nrow(a_terms)
  n_b <- nrow(b_terms)
  n_overlap <- joined[in_a == TRUE & in_b == TRUE, .N]
  n_union <- joined[, .N]
  data.table(
    n_count_terms = n_a,
    n_reference_terms = n_b,
    n_overlap = n_overlap,
    jaccard = ifelse(n_union == 0, NA_real_, n_overlap / n_union),
    count_terms_recovered_pct = ifelse(n_a == 0, NA_real_, 100 * n_overlap / n_a),
    reference_terms_recovered_pct = ifelse(n_b == 0, NA_real_, 100 * n_overlap / n_b)
  )
}

extract_common_terms <- function(a, b, ref_name) {
  common <- merge(
    a,
    b[, .(
      source, term_id,
      reference_term_name = term_name,
      reference_p_value = p_value,
      reference_intersection_size = intersection_size,
      reference_query_size = query_size,
      reference_hit_gene_symbols = hit_gene_symbols
    )],
    by = term_key_cols,
    all = FALSE,
    suffixes = c("_count", "_reference")
  )
  if (nrow(common) == 0) return(data.table())
  common[, reference := ref_name]
  common[, .(
    reference, fraction, direction, source, term_id,
    term_name,
    count_p_value = p_value,
    reference_p_value,
    count_intersection_size = intersection_size,
    reference_intersection_size,
    count_hit_gene_symbols = hit_gene_symbols,
    reference_hit_gene_symbols
  )]
}

summary_rows <- list()
common_rows <- list()

for (i in seq_len(nrow(metric_map))) {
  frac <- metric_map$fraction[i]
  metric_frac <- metric_map$metric_fraction[i]
  for (dirn in c("Up", "Down")) {
    count_set <- count[comparison == "Interaction" & fraction == frac & direction == dirn]
    metric_set <- metric[comparison == "Interaction" & fraction == metric_frac & direction == dirn]
    dseq_set <- dseq[comparison == "Vin_Resistant_vs_Sensitive" & fraction == frac & direction == dirn]

    s_metric <- compare_sets(count_set, metric_set)
    s_metric[, `:=`(
      reference = "Metric limma interaction",
      comparison_basis = metric_map$label[i],
      fraction = frac,
      direction = dirn
    )]
    summary_rows[[length(summary_rows) + 1L]] <- s_metric
    common_rows[[length(common_rows) + 1L]] <- extract_common_terms(count_set, metric_set, "Metric limma interaction")

    s_dseq <- compare_sets(count_set, dseq_set)
    s_dseq[, `:=`(
      reference = "DESeq fraction resistant-vs-sensitive under VCR",
      comparison_basis = paste0(frac, " vs ", frac),
      fraction = frac,
      direction = dirn
    )]
    summary_rows[[length(summary_rows) + 1L]] <- s_dseq
    common_rows[[length(common_rows) + 1L]] <- extract_common_terms(count_set, dseq_set, "DESeq fraction resistant-vs-sensitive under VCR")
  }
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
setcolorder(summary_dt, c(
  "reference", "comparison_basis", "fraction", "direction",
  "n_count_terms", "n_reference_terms", "n_overlap", "jaccard",
  "count_terms_recovered_pct", "reference_terms_recovered_pct"
))
setorder(summary_dt, reference, fraction, direction)

common_dt <- rbindlist(common_rows, fill = TRUE)
if (nrow(common_dt) > 0) setorder(common_dt, reference, fraction, direction, count_p_value)

fwrite(summary_dt, file.path(out_dir, "pathway_concordance_summary_exact_term_overlap.csv"))
fwrite(common_dt, file.path(out_dir, "pathway_concordance_common_terms_exact_term_overlap.csv"))

plot_dt <- copy(summary_dt)
plot_dt[, direction := factor(direction, levels = c("Up", "Down"))]
plot_dt[, fraction := factor(fraction, levels = c("SSU", "RS", "DS"))]
plot_dt[, label := paste0(n_overlap, "/", n_count_terms)]

p <- ggplot(plot_dt, aes(direction, fraction, fill = count_terms_recovered_pct)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 4, fontface = "bold") +
  facet_wrap(~reference, nrow = 1) +
  scale_fill_gradient(low = "grey95", high = "#2C7BB6", na.value = "grey95", limits = c(0, 100)) +
  labs(
    title = "Pathway concordance with P-site fraction limma true interaction",
    subtitle = "Exact overlap of filtered g:Profiler Reactome/GO:BP/GO:CC terms; tile label = shared/count-limma terms",
    x = NULL,
    y = NULL,
    fill = "% count terms\nrecovered"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(file.path(out_dir, "pathway_concordance_exact_term_overlap_heatmap.png"), p, width = 11, height = 4.6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "pathway_concordance_exact_term_overlap_heatmap.pdf"), p, width = 11, height = 4.6, bg = "white")

cat("\nPathway concordance summary:\n")
print(summary_dt[, .(
  reference, fraction, direction,
  count_terms = n_count_terms,
  reference_terms = n_reference_terms,
  overlap = n_overlap,
  jaccard = round(jaccard, 3),
  count_recovered_pct = round(count_terms_recovered_pct, 1)
)])

cat("\nCommon terms:\n")
if (nrow(common_dt) == 0) {
  cat("No exact term overlaps found.\n")
} else {
  print(common_dt[, .(
    reference, fraction, direction, source, term_name,
    count_p_value = signif(count_p_value, 3),
    reference_p_value = signif(reference_p_value, 3)
  )])
}

cat("\nSaved pathway concordance outputs to:\n", out_dir, "\n", sep = "")
