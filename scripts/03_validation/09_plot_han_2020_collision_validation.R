# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

HAN_FILE <- external_path("Downloads", "NIHMS1591487-supplement-Table_S1.csv")
GENE_METRICS <- analysis_path("Translation_indexes_fixed", "Gene_Level_Clean", "gene_level_clean_collision_complete_8_samples.csv")
OUT_DIR <- analysis_path("Translation_indexes_fixed", "Validation_Plots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, width, height) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = width, height = height)
}

format_p <- function(p) {
  ifelse(p < 1e-4, "p < 0.0001", paste0("p = ", signif(p, 3)))
}

han <- fread(HAN_FILE)
han[, gene := toupper(gene)]
han_gene <- han[, .(
  han_n_pause_sites = .N,
  han_max_pause_score = max(pause_score, na.rm = TRUE),
  han_median_pause_score = median(pause_score, na.rm = TRUE)
), by = gene][order(-han_max_pause_score)]
han_gene[, han_rank := seq_len(.N)]
han_top200 <- han_gene[han_rank <= 200, gene]

metrics <- fread(GENE_METRICS)
gene_med <- metrics[, .(
  gene_name = toupper(names(sort(table(gene_name), decreasing = TRUE))[1]),
  collision_score = median(collision_score, na.rm = TRUE),
  collision_score_DMSO = median(collision_score[treatment == "DMSO"], na.rm = TRUE),
  collision_score_VCR = median(collision_score[treatment == "VCR"], na.rm = TRUE),
  rs_core_cpm = median(rs_core_cpm, na.rm = TRUE)
), by = gene_id_clean]

gene_med[, han_any_collision_gene := gene_name %chin% han_gene$gene]
gene_med[, han_top200_gene := gene_name %chin% han_top200]
gene_med <- merge(
  gene_med,
  han_gene,
  by.x = "gene_name",
  by.y = "gene",
  all.x = TRUE
)

run_wilcox <- function(flag_col, y_col, label) {
  x <- gene_med[get(flag_col) == TRUE & is.finite(get(y_col)), get(y_col)]
  bg <- gene_med[get(flag_col) == FALSE & is.finite(get(y_col)), get(y_col)]
  data.table(
    comparison = label,
    metric = y_col,
    n_han = length(x),
    n_background = length(bg),
    median_han = median(x, na.rm = TRUE),
    median_background = median(bg, na.rm = TRUE),
    median_difference = median(x, na.rm = TRUE) - median(bg, na.rm = TRUE),
    wilcox_p_han_greater = suppressWarnings(wilcox.test(x, bg, alternative = "greater", exact = FALSE)$p.value)
  )
}

tests <- rbindlist(list(
  run_wilcox("han_top200_gene", "collision_score", "Han top 200 genes by max pause score"),
  run_wilcox("han_top200_gene", "collision_score_DMSO", "Han top 200 genes by max pause score"),
  run_wilcox("han_any_collision_gene", "collision_score", "Han any collision gene"),
  run_wilcox("han_any_collision_gene", "collision_score_DMSO", "Han any collision gene")
))
fwrite(tests, file.path(OUT_DIR, "han_2020_collision_rank_validation_stats.csv"))
fwrite(gene_med, file.path(OUT_DIR, "han_2020_collision_rank_validation_joined_data.csv"))

plot_dt <- rbindlist(list(
  gene_med[is.finite(collision_score), .(
    gene_name,
    collision_score,
    group = fifelse(han_top200_gene, "Han top 200 collision genes", "Background"),
    metric_label = "All samples median collision_score"
  )],
  gene_med[is.finite(collision_score_DMSO), .(
    gene_name,
    collision_score = collision_score_DMSO,
    group = fifelse(han_top200_gene, "Han top 200 collision genes", "Background"),
    metric_label = "DMSO median collision_score"
  )]
))
plot_dt[, group := factor(group, levels = c("Background", "Han top 200 collision genes"))]
plot_dt[, metric_label := factor(metric_label, levels = c("All samples median collision_score", "DMSO median collision_score"))]

ann <- tests[comparison == "Han top 200 genes by max pause score"]
ann[, metric_label := fifelse(metric == "collision_score", "All samples median collision_score", "DMSO median collision_score")]
ann[, metric_label := factor(metric_label, levels = levels(plot_dt$metric_label))]
ann_pos <- plot_dt[, .(
  x = 1.5,
  y = quantile(collision_score, 0.98, na.rm = TRUE)
), by = metric_label]
ann <- merge(ann, ann_pos, by = "metric_label", all.x = TRUE)
ann[, label := paste0(
  "Han median = ", sprintf("%.3f", median_han),
  "\nBackground median = ", sprintf("%.3f", median_background),
  "\n", format_p(wilcox_p_han_greater),
  "\nn = ", n_han, " vs ", n_background
)]

p_violin <- ggplot(plot_dt, aes(x = group, y = collision_score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.72, linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.92, linewidth = 0.25) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 8, color = "black") +
  geom_label(
    data = ann,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.4,
    fill = "white",
    linewidth = 0.2
  ) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Background" = "#A7B0BA", "Han top 200 collision genes" = "#B8323B")) +
  labs(
    title = "Han 2020 Collision Genes vs This Dataset's Collision Score",
    subtitle = "Han top 200 genes ranked by maximum disome pause score; clean collision-complete gene set",
    x = NULL,
    y = "Collision score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 18, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 15)
  )
save_plot(p_violin, "han_2020_top200_collision_score_validation", 10.8, 5.8)

rank_dt <- gene_med[is.finite(collision_score)]
rank_dt[, our_rank := frank(-collision_score, ties.method = "average")]
rank_dt[, han_group := fifelse(han_top200_gene, "Han top 200", "Background")]
rank_dt[, han_group := factor(han_group, levels = c("Background", "Han top 200"))]

p_rank <- ggplot(rank_dt, aes(x = han_group, y = our_rank, fill = han_group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.72, linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.92, linewidth = 0.25) +
  scale_y_reverse() +
  scale_fill_manual(values = c("Background" = "#A7B0BA", "Han top 200" = "#B8323B")) +
  labs(
    title = "Rank Distribution of Han 2020 Top Collision Genes",
    subtitle = "Higher placement means smaller rank number in this dataset",
    x = NULL,
    y = "Rank by this dataset's collision_score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 18, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 15)
  )
save_plot(p_rank, "han_2020_top200_collision_rank_distribution", 7.5, 5.8)

overlap_summary <- data.table(
  han_pause_sites = nrow(han),
  han_unique_genes = nrow(han_gene),
  han_top200_unique_genes = length(han_top200),
  clean_collision_genes = nrow(gene_med),
  han_any_overlap_clean_genes = gene_med[han_any_collision_gene == TRUE, .N],
  han_top200_overlap_clean_genes = gene_med[han_top200_gene == TRUE, .N]
)
fwrite(overlap_summary, file.path(OUT_DIR, "han_2020_collision_overlap_summary.csv"))

cat("\nHan 2020 collision gene validation\n")
cat("==================================\n")
print(overlap_summary)
cat("\nRank-sum tests\n")
print(tests[, .(
  comparison,
  metric,
  n_han,
  n_background,
  median_han = round(median_han, 4),
  median_background = round(median_background, 4),
  median_difference = round(median_difference, 4),
  wilcox_p_han_greater = signif(wilcox_p_han_greater, 3)
)])
cat("\nOutputs written to:", OUT_DIR, "\n")
