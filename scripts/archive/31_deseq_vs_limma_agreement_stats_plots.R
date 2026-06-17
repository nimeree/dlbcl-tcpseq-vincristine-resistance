# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

base_dir <- analysis_path()
limma_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05")
limma_results_dir <- file.path(limma_dir, "Results")
out_dir <- file.path(limma_dir, "DESeq2_vs_limma_agreement")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_cut <- 0.05
lfc_cut <- 0.7

metrics <- c(
  scanning = "scanning_score",
  ribosome_engagement = "ribosome_efficiency_score",
  protein_output = "protein_output_score",
  collision = "collision_score"
)

metric_labels <- c(
  scanning = "Scanning score",
  ribosome_engagement = "Ribosome engagement",
  protein_output = "Protein output",
  collision = "Collision score"
)

fractions <- c("SSU", "RS", "DS")
metric_fraction_map <- data.table(
  metric = c("scanning", "ribosome_engagement", "protein_output", "collision"),
  fraction = c("SSU", "RS", "RS", "DS"),
  expectation = c(
    "Scanning metric should most directly reflect SSU relative to RS",
    "Ribosome engagement should most directly reflect RS-associated signal",
    "Protein output used the same fraction-level signal as ribosome engagement in this implementation",
    "Collision metric should most directly reflect DS relative to RS"
  )
)

highlight_genes <- c("TRA2A", "LPXN", "HNRNPD", "MAPKBP1", "SEC24C")

read_limma <- function(metric_key, contrast) {
  f <- file.path(limma_results_dir, metrics[[metric_key]], paste0(contrast, "_limma_all_genes.csv"))
  d <- fread(f)
  d <- d[!is.na(gene_name) & gene_name != ""]
  d[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  d[, limma_dir := fifelse(
    P.Value < p_cut & logFC >= lfc_cut, "Up",
    fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS")
  )]
  d[, .(
    gene_id = gene_id_clean,
    gene_name,
    metric = metric_key,
    contrast,
    limma_logFC = logFC,
    limma_p = P.Value,
    limma_fdr = adj.P.Val,
    limma_dir
  )]
}

read_deseq <- function(fraction, contrast) {
  f <- file.path(base_dir, paste0("Fraction_", fraction), paste0(contrast, "_", fraction, "_results_all.csv"))
  d <- fread(f)
  d <- d[!is.na(gene_name) & gene_name != ""]
  d[, gene_id := sub("\\.\\d+$", "", gene_id)]
  d[, deseq_dir := fifelse(
    pvalue < p_cut & abs(log2FoldChange) >= lfc_cut & log2FoldChange >= 0, "Up",
    fifelse(pvalue < p_cut & abs(log2FoldChange) >= lfc_cut & log2FoldChange < 0, "Down", "NS")
  )]
  d[, .(
    gene_id,
    gene_name,
    fraction,
    deseq_contrast = contrast,
    deseq_log2FC = log2FoldChange,
    deseq_p = pvalue,
    deseq_fdr = padj,
    deseq_dir,
    baseMean
  )]
}

limma_dt <- rbindlist(lapply(names(metrics), function(m) {
  rbindlist(lapply(c("VCR_sensitive", "VCR_resistant", "Interaction"), read_limma, metric_key = m))
}), fill = TRUE)

deseq_pairwise <- rbindlist(lapply(fractions, function(fr) {
  rbindlist(lapply(c("Sensitive_Vin_vs_DMSO", "Resistant_Vin_vs_DMSO"), read_deseq, fraction = fr))
}), fill = TRUE)

make_delta <- function(fr) {
  sens <- read_deseq(fr, "Sensitive_Vin_vs_DMSO")
  res <- read_deseq(fr, "Resistant_Vin_vs_DMSO")
  setnames(sens, names(sens), paste0("sens_", names(sens)))
  setnames(res, names(res), paste0("res_", names(res)))
  d <- merge(
    res,
    sens,
    by.x = c("res_gene_id", "res_gene_name"),
    by.y = c("sens_gene_id", "sens_gene_name")
  )
  d[, .(
    gene_id = res_gene_id,
    gene_name = res_gene_name,
    fraction = fr,
    deseq_contrast = "Interaction_delta_ResistantMinusSensitive",
    deseq_log2FC = res_deseq_log2FC - sens_deseq_log2FC,
    resistant_deseq_log2FC = res_deseq_log2FC,
    sensitive_deseq_log2FC = sens_deseq_log2FC,
    resistant_p = res_deseq_p,
    sensitive_p = sens_deseq_p,
    resistant_dir = res_deseq_dir,
    sensitive_dir = sens_deseq_dir,
    deseq_dir = fifelse(
      abs(res_deseq_log2FC - sens_deseq_log2FC) >= lfc_cut & (res_deseq_log2FC - sens_deseq_log2FC) >= 0,
      "Up",
      fifelse(abs(res_deseq_log2FC - sens_deseq_log2FC) >= lfc_cut, "Down", "NS")
    ),
    baseMean = (res_baseMean + sens_baseMean) / 2
  )]
}

deseq_delta <- rbindlist(lapply(fractions, make_delta), fill = TRUE)

compare_sets <- list(
  VCR_sensitive = "Sensitive_Vin_vs_DMSO",
  VCR_resistant = "Resistant_Vin_vs_DMSO",
  Interaction = "Interaction_delta_ResistantMinusSensitive"
)

comparison_dt <- rbindlist(lapply(names(compare_sets), function(limma_contrast) {
  deseq_contrast <- compare_sets[[limma_contrast]]
  dseq <- if (deseq_contrast == "Interaction_delta_ResistantMinusSensitive") {
    deseq_delta
  } else {
    contrast_name <- deseq_contrast
    deseq_pairwise[deseq_contrast == contrast_name]
  }
  merge(
    limma_dt[contrast == limma_contrast],
    dseq,
    by = c("gene_id", "gene_name"),
    allow.cartesian = TRUE
  )
}), fill = TRUE)

comparison_dt <- merge(comparison_dt, metric_fraction_map, by = c("metric", "fraction"), all.x = TRUE)
comparison_dt[, direct_metric_fraction_pair := !is.na(expectation)]
comparison_dt[, metric_label := metric_labels[metric]]
comparison_dt[, metric_label := factor(metric_label, levels = metric_labels)]
comparison_dt[, fraction := factor(fraction, levels = fractions)]
comparison_dt[, contrast_label := fifelse(
  contrast == "VCR_sensitive", "Sensitive VCR vs DMSO",
  fifelse(contrast == "VCR_resistant", "Resistant VCR vs DMSO", "Interaction: resistant response - sensitive response")
)]

cor_stats <- comparison_dt[
  ,
  {
    ok <- is.finite(limma_logFC) & is.finite(deseq_log2FC)
    ct <- suppressWarnings(cor.test(limma_logFC[ok], deseq_log2FC[ok], method = "spearman", exact = FALSE))
    list(
      n_genes = sum(ok),
      spearman_rho = unname(ct$estimate),
      p_value = ct$p.value,
      pearson_r = suppressWarnings(cor(limma_logFC[ok], deseq_log2FC[ok], method = "pearson"))
    )
  },
  by = .(contrast, contrast_label, metric, metric_label, fraction, direct_metric_fraction_pair)
]
cor_stats[, BH_p_value := p.adjust(p_value, method = "BH")]
setorder(cor_stats, contrast, metric, fraction)

direct_pairs <- comparison_dt[direct_metric_fraction_pair == TRUE]
direct_pairs[, sign_class := fifelse(
  limma_dir == "NS", "limma NS",
  fifelse(deseq_dir == "NS", "limma-only",
          fifelse(limma_dir == deseq_dir, "same direction", "opposite direction"))
)]

concordance_stats <- direct_pairs[
  limma_dir != "NS",
  .(
    n_limma_sig = .N,
    n_same_direction = sum(sign_class == "same direction"),
    n_opposite_direction = sum(sign_class == "opposite direction"),
    n_limma_only = sum(sign_class == "limma-only"),
    pct_same_among_deseq_called = ifelse(
      sum(sign_class %in% c("same direction", "opposite direction")) > 0,
      100 * sum(sign_class == "same direction") / sum(sign_class %in% c("same direction", "opposite direction")),
      NA_real_
    )
  ),
  by = .(contrast, contrast_label, metric, metric_label, fraction)
]

discordant_genes <- direct_pairs[
  limma_dir != "NS" & sign_class %in% c("opposite direction", "limma-only"),
  .(
    gene_id,
    gene_name,
    contrast,
    contrast_label,
    metric,
    metric_label,
    fraction,
    limma_logFC,
    limma_p,
    limma_fdr,
    limma_dir,
    deseq_log2FC,
    deseq_p,
    deseq_fdr,
    deseq_dir,
    sign_class,
    resistant_deseq_log2FC,
    sensitive_deseq_log2FC,
    resistant_dir,
    sensitive_dir
  )
][order(contrast, metric, sign_class, gene_name)]

key_gene_profile <- direct_pairs[
  gene_name %in% highlight_genes,
  .(
    gene_name,
    contrast,
    contrast_label,
    metric,
    metric_label,
    fraction,
    limma_logFC,
    limma_p,
    limma_dir,
    deseq_log2FC,
    deseq_p,
    deseq_dir,
    resistant_deseq_log2FC,
    sensitive_deseq_log2FC,
    sign_class
  )
][order(gene_name, contrast, metric)]

write.csv(comparison_dt, file.path(out_dir, "deseq2_limma_all_metric_fraction_comparisons.csv"), row.names = FALSE)
write.csv(cor_stats, file.path(out_dir, "deseq2_limma_spearman_correlation_stats.csv"), row.names = FALSE)
write.csv(concordance_stats, file.path(out_dir, "deseq2_limma_direct_pair_concordance_counts.csv"), row.names = FALSE)
write.csv(discordant_genes, file.path(out_dir, "deseq2_limma_direct_pair_discordant_or_limma_only_genes.csv"), row.names = FALSE)
write.csv(key_gene_profile, file.path(out_dir, "key_gene_deseq2_limma_profiles.csv"), row.names = FALSE)
write.csv(metric_fraction_map, file.path(out_dir, "metric_to_fraction_interpretation_map.csv"), row.names = FALSE)

heat_dt <- cor_stats
plot_metric_order <- c("Scanning score", "Ribosome engagement", "Collision score")
heat_dt <- heat_dt[as.character(metric_label) %in% plot_metric_order]
heat_dt[, metric_label := factor(as.character(metric_label), levels = plot_metric_order)]
heat_dt[, rho_label := sprintf("%.2f", spearman_rho)]
heat_dt[, comparison_label := paste(metric_label, "vs", fraction)]

p_heat <- ggplot(heat_dt, aes(fraction, metric_label, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = rho_label), size = 3.1) +
  facet_wrap(~ contrast_label, ncol = 1) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C",
    midpoint = 0, limits = c(-0.45, 0.45), oob = squish,
    name = "Spearman\nrho"
  ) +
  labs(
    title = "limma metric agreement with fraction-specific DESeq2",
    subtitle = "Values are Spearman correlations between limma metric logFC and DESeq2 fraction log2FC across genes",
    x = "DESeq2 TCP-seq fraction",
    y = "limma translation metric"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "A_deseq2_limma_correlation_heatmap.png"), p_heat, width = 8.5, height = 8.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "A_deseq2_limma_correlation_heatmap.pdf"), p_heat, width = 8.5, height = 8.2, bg = "white")

scatter_dt <- direct_pairs[contrast == "Interaction" & as.character(metric_label) %in% plot_metric_order]
scatter_dt[, metric_label := factor(as.character(metric_label), levels = plot_metric_order)]
scatter_dt[, label_gene := fifelse(gene_name %in% highlight_genes, gene_name, "")]
scatter_dt[, point_group := fifelse(gene_name %in% highlight_genes, gene_name, "Other genes")]

p_scatter <- ggplot(scatter_dt, aes(deseq_log2FC, limma_logFC)) +
  geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", color = "grey45") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", color = "grey45") +
  geom_point(aes(color = sign_class), alpha = 0.62, size = 1.35) +
  geom_text(
    data = scatter_dt[label_gene != ""],
    aes(label = label_gene),
    size = 2.8,
    vjust = -0.7,
    check_overlap = TRUE,
    color = "black"
  ) +
  facet_wrap(~ metric_label, scales = "free", ncol = 3) +
  scale_color_manual(
    values = c(
      "same direction" = "#2563EB",
      "opposite direction" = "#DC2626",
      "limma-only" = "#6B7280",
      "limma NS" = "#D1D5DB"
    ),
    name = NULL
  ) +
  labs(
    title = "Interaction contrast: limma metric changes versus DESeq2 fraction-response changes",
    subtitle = "DESeq2 interaction proxy = Resistant VCR response minus Sensitive VCR response for the mapped fraction",
    x = "DESeq2 fraction delta log2FC",
    y = "limma metric interaction logFC"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "B_interaction_direct_pair_scatterplots.png"), p_scatter, width = 12, height = 4.8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "B_interaction_direct_pair_scatterplots.pdf"), p_scatter, width = 12, height = 4.8, bg = "white")

bar_dt <- direct_pairs[limma_dir != "NS" & as.character(metric_label) %in% plot_metric_order, .N, by = .(contrast_label, metric_label, sign_class)]
bar_dt[, metric_label := factor(as.character(metric_label), levels = plot_metric_order)]
bar_dt[, sign_class := factor(sign_class, levels = c("same direction", "opposite direction", "limma-only"))]

p_bar <- ggplot(bar_dt, aes(metric_label, N, fill = sign_class)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5), size = 2.8, color = "white") +
  facet_wrap(~ contrast_label, ncol = 1) +
  scale_fill_manual(
    values = c(
      "same direction" = "#2563EB",
      "opposite direction" = "#DC2626",
      "limma-only" = "#6B7280"
    ),
    name = NULL
  ) +
  labs(
    title = "How often limma-significant genes are supported by the matching DESeq2 fraction",
    subtitle = "Same direction/opposite direction use the matching fraction; limma-only means the mapped DESeq2 fraction did not pass the same raw p and |logFC| rule",
    x = NULL,
    y = "limma-significant genes"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "C_limma_significant_concordance_barplot.png"), p_bar, width = 8.4, height = 8.0, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "C_limma_significant_concordance_barplot.pdf"), p_bar, width = 8.4, height = 8.0, bg = "white")

key_dt <- key_gene_profile[contrast == "Interaction" & as.character(metric_label) %in% plot_metric_order]
key_dt[, metric_label := factor(as.character(metric_label), levels = plot_metric_order)]
key_long <- melt(
  key_dt,
  id.vars = c("gene_name", "metric_label", "fraction", "sign_class"),
  measure.vars = c("limma_logFC", "deseq_log2FC"),
  variable.name = "analysis",
  value.name = "logFC"
)
key_long[, analysis := fifelse(analysis == "limma_logFC", "limma metric", "DESeq2 mapped fraction")]

p_key <- ggplot(key_long, aes(metric_label, logFC, fill = analysis)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  facet_wrap(~ gene_name, ncol = 3) +
  scale_fill_manual(values = c("limma metric" = "#111827", "DESeq2 mapped fraction" = "#9CA3AF"), name = NULL) +
  labs(
    title = "Key genes: limma interaction metric versus mapped DESeq2 fraction delta",
    subtitle = "This highlights why candidate genes can be metric-driven even when fraction-level DESeq2 is weak",
    x = NULL,
    y = "logFC"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "D_key_gene_limma_vs_deseq2_barplot.png"), p_key, width = 10.5, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "D_key_gene_limma_vs_deseq2_barplot.pdf"), p_key, width = 10.5, height = 6.2, bg = "white")

cat("\nOutput directory:\n", out_dir, "\n", sep = "")
cat("\nDirect metric-to-fraction correlation stats:\n")
print(cor_stats[direct_metric_fraction_pair == TRUE][
  ,
  .(contrast, metric, fraction, n_genes, spearman_rho, p_value, BH_p_value)
][order(contrast, metric)])

cat("\nConcordance counts for limma-significant genes against mapped DESeq2 fractions:\n")
print(concordance_stats[order(contrast, metric)])

cat("\nKey gene interaction profiles:\n")
print(key_gene_profile[contrast == "Interaction"][order(gene_name, metric)])
