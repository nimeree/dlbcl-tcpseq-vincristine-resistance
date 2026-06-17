# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

base_dir <- analysis_path()
count_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
metric_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05", "Results")
out_dir <- file.path(count_dir, "True_interaction_concordance_with_metric_limma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_cut <- 0.05
lfc_cut <- 0.7

gene_counts_long <- fread(file.path(count_dir, "psite_gene_counts_long_by_fraction_sample.csv"))
sample_meta <- fread(file.path(count_dir, "psite_fraction_limma_sample_metadata.csv"))
sample_meta[, condition := factor(condition, levels = c("Sensitive_DMSO", "Sensitive_Vin", "Resistant_DMSO", "Resistant_Vin"))]

direction_call <- function(logfc, pval) {
  fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc >= lfc_cut, "Up",
          fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc <= -lfc_cut, "Down", "NS"))
}

run_true_interaction_fraction <- function(frac) {
  frac_dir <- file.path(count_dir, paste0("Fraction_", frac))
  d <- gene_counts_long[fraction == frac]
  meta <- unique(sample_meta[fraction == frac][order(condition, replicate)])
  sample_levels <- meta$sample

  mat_dt <- dcast(d, gene_id_clean + gene_name ~ sample, value.var = "psite_count", fill = 0)
  gene_info <- mat_dt[, .(gene_id_clean, gene_name)]
  mat <- as.matrix(mat_dt[, ..sample_levels])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(gene_info$gene_id_clean)

  keep <- rowSums(mat >= 10) >= 2
  mat <- mat[keep, , drop = FALSE]
  gene_info <- gene_info[keep]
  lib_sizes <- colSums(mat)
  log_cpm <- log2(t(t(mat + 0.5) / pmax(lib_sizes, 1)) * 1e6)

  group <- factor(meta$condition, levels = c("Sensitive_DMSO", "Sensitive_Vin", "Resistant_DMSO", "Resistant_Vin"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- lmFit(log_cpm, design)
  fit <- eBayes(fit, trend = TRUE)
  contrast_matrix <- makeContrasts(
    Interaction = (Resistant_Vin - Resistant_DMSO) - (Sensitive_Vin - Sensitive_DMSO),
    levels = design
  )
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE)

  tt <- as.data.table(topTable(fit2, coef = "Interaction", number = Inf, sort.by = "none"), keep.rownames = "row_key")
  tt[, gene_id_clean := gene_info$gene_id_clean]
  tt[, gene_name := gene_info$gene_name]
  tt[, significant_rawP0.05_lfc0.7 := !is.na(P.Value) & P.Value < p_cut & abs(logFC) >= lfc_cut]
  tt[, direction := direction_call(logFC, P.Value)]
  setcolorder(tt, c("gene_id_clean", "gene_name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "significant_rawP0.05_lfc0.7", "direction"))
  tt <- tt[order(P.Value)]

  out_all <- file.path(frac_dir, paste0("Interaction_", frac, "_psite_limma_all_genes.csv"))
  out_sig <- file.path(frac_dir, paste0("Interaction_", frac, "_psite_limma_sig_rawP0.05_lfc0.7.csv"))
  fwrite(tt, out_all)
  fwrite(tt[significant_rawP0.05_lfc0.7 == TRUE], out_sig)
  tt
}

metric_map <- data.table(
  metric = c("Scanning", "Ribosome engagement", "Collision"),
  metric_folder = c("scanning_score", "ribosome_efficiency_score", "collision_score"),
  fraction = c("SSU", "RS", "DS")
)

read_metric_interaction <- function(folder) {
  f <- file.path(metric_dir, folder, "Interaction_limma_all_genes.csv")
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

read_count_interaction <- function(frac) {
  f <- file.path(count_dir, paste0("Fraction_", frac), paste0("Interaction_", frac, "_psite_limma_all_genes.csv"))
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

summary_counts <- rbindlist(lapply(metric_map$fraction, function(frac) {
  tt <- run_true_interaction_fraction(frac)
  data.table(
    fraction = frac,
    genes_tested = nrow(tt),
    n_sig = tt[significant_rawP0.05_lfc0.7 == TRUE, .N],
    n_up = tt[direction == "Up", .N],
    n_down = tt[direction == "Down", .N],
    median_logFC = median(tt$logFC, na.rm = TRUE)
  )
}))
fwrite(summary_counts, file.path(out_dir, "psite_fraction_limma_true_interaction_summary.csv"))

all_pairs <- rbindlist(lapply(seq_len(nrow(metric_map)), function(i) {
  m <- metric_map[i]
  joined <- merge(read_metric_interaction(m$metric_folder), read_count_interaction(m$fraction), by = "gene_key", all = FALSE)
  joined[, `:=`(
    gene_name = fifelse(!is.na(gene_name_metric) & nzchar(gene_name_metric), gene_name_metric, gene_name_count),
    metric = m$metric,
    fraction = m$fraction
  )]
  joined
}), use.names = TRUE)

all_pairs[, metric_sig := metric_direction %in% c("Up", "Down")]
all_pairs[, count_sig := count_direction %in% c("Up", "Down")]

cor_stats <- all_pairs[, {
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
}, by = .(metric, fraction)]

concordance_metric_sig <- all_pairs[metric_sig == TRUE, .(
  metric_sig_n = .N,
  count_same_direction = sum(count_sig & metric_direction == count_direction, na.rm = TRUE),
  count_opposite_direction = sum(count_sig & metric_direction != count_direction, na.rm = TRUE),
  count_not_sig = sum(!count_sig, na.rm = TRUE),
  count_same_direction_pct = 100 * sum(count_sig & metric_direction == count_direction, na.rm = TRUE) / .N
), by = .(metric, fraction)]

concordance_count_sig <- all_pairs[count_sig == TRUE, .(
  count_sig_n = .N,
  metric_same_direction = sum(metric_sig & metric_direction == count_direction, na.rm = TRUE),
  metric_opposite_direction = sum(metric_sig & metric_direction != count_direction, na.rm = TRUE),
  metric_not_sig = sum(!metric_sig, na.rm = TRUE),
  metric_same_direction_pct = 100 * sum(metric_sig & metric_direction == count_direction, na.rm = TRUE) / .N
), by = .(metric, fraction)]

key_genes <- c("TRA2A", "LPXN", "HNRNPD", "MAPKBP1", "SEC24C", "HSP90AB1")
key_table <- all_pairs[gene_name %in% key_genes, .(
  gene_name, metric, fraction,
  metric_logFC, metric_P.Value, metric_direction,
  count_logFC, count_P.Value, count_direction
)]
setorder(key_table, gene_name, metric)

fwrite(all_pairs, file.path(out_dir, "true_interaction_metric_limma_vs_psite_count_limma_all_pairs.csv"))
fwrite(cor_stats, file.path(out_dir, "true_interaction_correlations.csv"))
fwrite(concordance_metric_sig, file.path(out_dir, "true_interaction_metric_sig_concordance.csv"))
fwrite(concordance_count_sig, file.path(out_dir, "true_interaction_count_sig_concordance.csv"))
fwrite(key_table, file.path(out_dir, "true_interaction_key_genes.csv"))

plot_dt <- merge(
  all_pairs,
  cor_stats[, .(metric, fraction, pearson_r, spearman_rho)],
  by = c("metric", "fraction"),
  all.x = TRUE
)
plot_dt[, metric := factor(metric, levels = c("Scanning", "Ribosome engagement", "Collision"))]
plot_dt[, sig_group := fifelse(metric_sig & count_sig & metric_direction == count_direction, "Both sig, same direction",
                        fifelse(metric_sig & count_sig & metric_direction != count_direction, "Both sig, opposite direction",
                        fifelse(metric_sig, "Metric only",
                        fifelse(count_sig, "Count only", "Neither"))))]
plot_dt[, panel_label := paste0(metric, " vs ", fraction,
                                "\nPearson r=", round(pearson_r, 2),
                                "; Spearman rho=", round(spearman_rho, 2))]

plot_cols <- c(
  "Both sig, same direction" = "#2F855A",
  "Both sig, opposite direction" = "#C53030",
  "Metric only" = "#805AD5",
  "Count only" = "#DD6B20",
  "Neither" = "grey78"
)

scatter_base <- ggplot(plot_dt, aes(metric_logFC, count_logFC)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_point(aes(color = sig_group), alpha = 0.65, size = 0.95) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.45) +
  scale_color_manual(values = plot_cols, breaks = names(plot_cols)) +
  labs(
    title = "True interaction concordance: metric limma vs P-site count limma",
    subtitle = "Count interaction = (Resistant Vin - Resistant DMSO) - (Sensitive Vin - Sensitive DMSO)",
    x = "Metric limma interaction logFC",
    y = "P-site count limma interaction logFC",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

combined_plot <- scatter_base +
  facet_wrap(~panel_label, scales = "free", nrow = 1)

ggsave(
  file.path(out_dir, "true_interaction_metric_vs_psite_count_limma_scatter_combined.png"),
  combined_plot,
  width = 12.5,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

for (m in levels(plot_dt$metric)) {
  one <- plot_dt[metric == m]
  p <- scatter_base %+% one +
    facet_wrap(~panel_label, scales = "free") +
    theme(legend.position = "bottom")
  safe_name <- gsub("[^A-Za-z0-9]+", "_", tolower(m))
  ggsave(
    file.path(out_dir, paste0("true_interaction_", safe_name, "_vs_", unique(one$fraction), "_scatter.png")),
    p,
    width = 6.2,
    height = 5.4,
    dpi = 300,
    bg = "white"
  )
}

cat("\nP-site count limma true interaction summary:\n")
print(summary_counts)

cat("\nTrue interaction concordance: metric limma interaction vs P-site count limma interaction\n")
print(cor_stats[, .(
  metric, fraction, n_genes,
  pearson_r = round(pearson_r, 3),
  spearman_rho = round(spearman_rho, 3),
  metric_sig_n, count_sig_n, both_sig_n, same_direction_n, opposite_direction_n
)])

cat("\nAmong metric-interaction significant genes:\n")
print(concordance_metric_sig[, .(
  metric, fraction, metric_sig_n, count_same_direction, count_opposite_direction, count_not_sig,
  count_same_direction_pct = round(count_same_direction_pct, 1)
)])

cat("\nAmong count-interaction significant genes:\n")
print(concordance_count_sig[, .(
  metric, fraction, count_sig_n, metric_same_direction, metric_opposite_direction, metric_not_sig,
  metric_same_direction_pct = round(metric_same_direction_pct, 1)
)])

cat("\nKey genes:\n")
print(key_table[, .(
  gene_name, metric, fraction,
  metric_logFC = round(metric_logFC, 3), metric_P.Value = signif(metric_P.Value, 3), metric_direction,
  count_logFC = round(count_logFC, 3), count_P.Value = signif(count_P.Value, 3), count_direction
)])

cat("\nSaved true interaction concordance tables to:\n", out_dir, "\n", sep = "")
