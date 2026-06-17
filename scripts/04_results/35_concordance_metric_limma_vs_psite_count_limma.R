# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- analysis_path()
metric_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05", "Results")
count_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
out_dir <- file.path(count_dir, "Concordance_with_metric_limma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_cut <- 0.05
lfc_cut <- 0.7

metric_map <- data.table(
  metric = c("Scanning", "Ribosome engagement", "Collision"),
  metric_folder = c("scanning_score", "ribosome_efficiency_score", "collision_score"),
  fraction = c("SSU", "RS", "DS")
)

contrast_map <- data.table(
  metric_contrast = c("VCR_sensitive", "VCR_resistant", "Interaction"),
  count_contrast = c("Sensitive_Vin_vs_DMSO", "Resistant_Vin_vs_DMSO", "Vin_Resistant_vs_Sensitive"),
  contrast_label = c(
    "Sensitive VCR vs DMSO",
    "Resistant VCR vs DMSO",
    "Resistant vs sensitive under VCR"
  )
)

metric_file <- function(folder, contrast) {
  file.path(metric_dir, folder, paste0(contrast, "_limma_all_genes.csv"))
}

count_file <- function(fraction, contrast) {
  file.path(count_dir, paste0("Fraction_", fraction),
            paste0(contrast, "_", fraction, "_psite_limma_all_genes.csv"))
}

direction_call <- function(logfc, pval) {
  fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc >= lfc_cut, "Up",
          fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc <= -lfc_cut, "Down", "NS"))
}

read_metric <- function(folder, contrast) {
  f <- metric_file(folder, contrast)
  if (!file.exists(f)) stop("Missing metric limma file: ", f)
  d <- fread(f)
  if (!"gene_name" %in% names(d)) d[, gene_name := gene_id_clean]
  d[, gene_key := fifelse(!is.na(gene_id_clean) & nzchar(gene_id_clean), gene_id_clean, gene_name)]
  d[, gene_key := sub("\\.\\d+$", "", gene_key)]
  d[, .(
    gene_key,
    gene_name_metric = gene_name,
    metric_logFC = logFC,
    metric_P.Value = P.Value,
    metric_adj.P.Val = adj.P.Val,
    metric_direction = direction_call(logFC, P.Value)
  )]
}

read_count <- function(fraction, contrast) {
  f <- count_file(fraction, contrast)
  if (!file.exists(f)) stop("Missing count limma file: ", f)
  d <- fread(f)
  if (!"gene_name" %in% names(d)) d[, gene_name := gene_id_clean]
  d[, gene_key := fifelse(!is.na(gene_id_clean) & nzchar(gene_id_clean), gene_id_clean, gene_name)]
  d[, gene_key := sub("\\.\\d+$", "", gene_key)]
  d[, .(
    gene_key,
    gene_name_count = gene_name,
    count_logFC = logFC,
    count_P.Value = P.Value,
    count_adj.P.Val = adj.P.Val,
    count_direction = direction_call(logFC, P.Value)
  )]
}

all_pairs <- list()
for (i in seq_len(nrow(metric_map))) {
  for (j in seq_len(nrow(contrast_map))) {
    m <- metric_map[i]
    c <- contrast_map[j]
    metric_dt <- read_metric(m$metric_folder, c$metric_contrast)
    count_dt <- read_count(m$fraction, c$count_contrast)
    joined <- merge(metric_dt, count_dt, by = "gene_key", all = FALSE)
    joined[, `:=`(
      gene_name = fifelse(!is.na(gene_name_metric) & nzchar(gene_name_metric), gene_name_metric, gene_name_count),
      metric = m$metric,
      fraction = m$fraction,
      metric_contrast = c$metric_contrast,
      count_contrast = c$count_contrast,
      contrast_label = c$contrast_label
    )]
    all_pairs[[length(all_pairs) + 1L]] <- joined
  }
}

all <- rbindlist(all_pairs, use.names = TRUE)
all[, metric_sig := metric_direction %in% c("Up", "Down")]
all[, count_sig := count_direction %in% c("Up", "Down")]
all[, direction_relation := fifelse(metric_sig & count_sig & metric_direction == count_direction, "same_direction",
                             fifelse(metric_sig & count_sig & metric_direction != count_direction, "opposite_direction",
                             fifelse(metric_sig & !count_sig, "metric_only",
                             fifelse(!metric_sig & count_sig, "count_only", "neither_sig"))))]

cor_stats <- all[, {
  ok <- is.finite(metric_logFC) & is.finite(count_logFC)
  pear <- suppressWarnings(cor.test(metric_logFC[ok], count_logFC[ok], method = "pearson"))
  spear <- suppressWarnings(cor.test(metric_logFC[ok], count_logFC[ok], method = "spearman", exact = FALSE))
  .(
    n_genes = sum(ok),
    pearson_r = unname(pear$estimate),
    pearson_p = pear$p.value,
    spearman_rho = unname(spear$estimate),
    spearman_p = spear$p.value,
    metric_sig_n = sum(metric_sig, na.rm = TRUE),
    count_sig_n = sum(count_sig, na.rm = TRUE),
    both_sig_n = sum(metric_sig & count_sig, na.rm = TRUE),
    same_direction_n = sum(metric_sig & count_sig & metric_direction == count_direction, na.rm = TRUE),
    opposite_direction_n = sum(metric_sig & count_sig & metric_direction != count_direction, na.rm = TRUE)
  )
}, by = .(metric, fraction, metric_contrast, count_contrast, contrast_label)]

concordance_by_metric_sig <- all[metric_sig == TRUE, .(
  metric_sig_n = .N,
  count_same_direction = sum(count_sig & metric_direction == count_direction, na.rm = TRUE),
  count_opposite_direction = sum(count_sig & metric_direction != count_direction, na.rm = TRUE),
  count_not_sig = sum(!count_sig, na.rm = TRUE),
  count_same_direction_pct = 100 * sum(count_sig & metric_direction == count_direction, na.rm = TRUE) / .N
), by = .(metric, fraction, metric_contrast, count_contrast, contrast_label)]

concordance_by_count_sig <- all[count_sig == TRUE, .(
  count_sig_n = .N,
  metric_same_direction = sum(metric_sig & metric_direction == count_direction, na.rm = TRUE),
  metric_opposite_direction = sum(metric_sig & metric_direction != count_direction, na.rm = TRUE),
  metric_not_sig = sum(!metric_sig, na.rm = TRUE),
  metric_same_direction_pct = 100 * sum(metric_sig & metric_direction == count_direction, na.rm = TRUE) / .N
), by = .(metric, fraction, metric_contrast, count_contrast, contrast_label)]

overall <- all[, {
  ok <- is.finite(metric_logFC) & is.finite(count_logFC)
  pear <- suppressWarnings(cor.test(metric_logFC[ok], count_logFC[ok], method = "pearson"))
  spear <- suppressWarnings(cor.test(metric_logFC[ok], count_logFC[ok], method = "spearman", exact = FALSE))
  .(
    n_rows = sum(ok),
    pearson_r = unname(pear$estimate),
    pearson_p = pear$p.value,
    spearman_rho = unname(spear$estimate),
    spearman_p = spear$p.value
  )
}]

key_genes <- c("TRA2A", "LPXN", "HNRNPD", "MAPKBP1", "SEC24C", "HSP90AB1")
key_table <- all[gene_name %in% key_genes, .(
  gene_name, metric, fraction, contrast_label,
  metric_logFC, metric_P.Value, metric_direction,
  count_logFC, count_P.Value, count_direction
)]
setorder(key_table, gene_name, metric, contrast_label)

fwrite(all, file.path(out_dir, "metric_limma_vs_psite_count_limma_all_pairs.csv"))
fwrite(cor_stats, file.path(out_dir, "metric_limma_vs_psite_count_limma_correlations.csv"))
fwrite(concordance_by_metric_sig, file.path(out_dir, "metric_sig_gene_concordance_with_count_limma.csv"))
fwrite(concordance_by_count_sig, file.path(out_dir, "count_sig_gene_concordance_with_metric_limma.csv"))
fwrite(key_table, file.path(out_dir, "key_gene_metric_vs_psite_count_limma.csv"))

cat("\nOverall logFC correlation across all metric/fraction/contrast pairs:\n")
print(overall)

cat("\nPer metric/fraction/contrast logFC correlations:\n")
print(cor_stats[, .(
  metric, fraction, contrast_label, n_genes,
  pearson_r = round(pearson_r, 3),
  spearman_rho = round(spearman_rho, 3),
  metric_sig_n, count_sig_n, both_sig_n,
  same_direction_n, opposite_direction_n
)])

cat("\nAmong metric-significant genes, how often are the P-site count limma calls also significant in the same direction?\n")
print(concordance_by_metric_sig[, .(
  metric, fraction, contrast_label, metric_sig_n,
  count_same_direction, count_opposite_direction, count_not_sig,
  count_same_direction_pct = round(count_same_direction_pct, 1)
)])

cat("\nAmong P-site count-significant genes, how often are the metric limma calls also significant in the same direction?\n")
print(concordance_by_count_sig[, .(
  metric, fraction, contrast_label, count_sig_n,
  metric_same_direction, metric_opposite_direction, metric_not_sig,
  metric_same_direction_pct = round(metric_same_direction_pct, 1)
)])

cat("\nKey genes:\n")
print(key_table[, .(
  gene_name, metric, fraction, contrast_label,
  metric_logFC = round(metric_logFC, 3), metric_P.Value = signif(metric_P.Value, 3), metric_direction,
  count_logFC = round(count_logFC, 3), count_P.Value = signif(count_P.Value, 3), count_direction
)])

cat("\nSaved concordance tables to:\n", out_dir, "\n", sep = "")
