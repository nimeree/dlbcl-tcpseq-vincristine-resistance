# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# Plot robust cross-framework genes as two heatmap panels:
# A) DS p-site + collision score robust genes
# B) RS p-site + ribosome engagement score robust genes

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE_DIR <- analysis_path()
PSITE_DIR <- file.path(BASE_DIR, "Psite_fraction_limma_lfc0.7_rawP0.05")
METRIC_DIR <- file.path(BASE_DIR, "Limma_translation_metrics_lfc0.7_rawP0.05")
OUT_DIR <- file.path(METRIC_DIR, "Cross_framework_synthesis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

P_CUT <- 0.05
LFC_CUT <- 0.7

read_psite <- function(frac, contrast) {
  f <- file.path(
    PSITE_DIR,
    paste0("Fraction_", frac),
    paste0(contrast, "_", frac, "_psite_limma_all_genes.csv")
  )
  d <- fread(f)
  if ("logFC" %in% names(d)) setnames(d, "logFC", "psite_logFC")
  if ("P.Value" %in% names(d)) setnames(d, "P.Value", "psite_P.Value")
  d[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  d[, psite_sig := psite_P.Value < P_CUT & abs(psite_logFC) >= LFC_CUT]
  d[, psite_dir := fifelse(psite_sig & psite_logFC > 0, "Up",
                    fifelse(psite_sig & psite_logFC < 0, "Down", "NS"))]
  unique(d[, .(gene_id_clean, gene_name, psite_logFC, psite_P.Value, psite_sig, psite_dir)], by = "gene_id_clean")
}

read_metric <- function(metric, contrast) {
  f <- file.path(METRIC_DIR, "Results", metric, paste0(contrast, "_limma_all_genes.csv"))
  d <- fread(f)
  if ("logFC" %in% names(d)) setnames(d, "logFC", "metric_logFC")
  if ("P.Value" %in% names(d)) setnames(d, "P.Value", "metric_P.Value")
  d[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  d[, metric_sig := metric_P.Value < P_CUT & abs(metric_logFC) >= LFC_CUT]
  d[, metric_dir := fifelse(metric_sig & metric_logFC > 0, "Up",
                     fifelse(metric_sig & metric_logFC < 0, "Down", "NS"))]
  unique(d[, .(gene_id_clean, gene_name, metric_logFC, metric_P.Value, metric_sig, metric_dir)], by = "gene_id_clean")
}

merge_pair <- function(frac, metric, contrast) {
  x <- merge(read_psite(frac, contrast), read_metric(metric, contrast),
             by = c("gene_id_clean", "gene_name"), all = TRUE)
  x[is.na(psite_sig), psite_sig := FALSE]
  x[is.na(metric_sig), metric_sig := FALSE]
  x[is.na(psite_dir), psite_dir := "NS"]
  x[is.na(metric_dir), metric_dir := "NS"]
  x[, same_dir := psite_sig & metric_sig & psite_dir == metric_dir]
  x[]
}

same_genes <- function(frac, metric, contrast) {
  merge_pair(frac, metric, contrast)[same_dir == TRUE, .(gene_id_clean, gene_name)]
}

make_matrix <- function(ids, frac, metric, psite_label, metric_label) {
  b <- merge_pair(frac, metric, "Resistance_baseline")[gene_id_clean %in% ids]
  i <- merge_pair(frac, metric, "Interaction")[gene_id_clean %in% ids]
  b <- b[, .(
    gene_id_clean, gene_name,
    baseline_psite = psite_logFC,
    baseline_metric = metric_logFC
  )]
  i <- i[, .(
    gene_id_clean, gene_name,
    interaction_psite = psite_logFC,
    interaction_metric = metric_logFC
  )]
  wide <- merge(b, i, by = c("gene_id_clean", "gene_name"), all = TRUE)
  setnames(
    wide,
    c("baseline_psite", "baseline_metric", "interaction_psite", "interaction_metric"),
    c(
      paste0("Baseline\n", psite_label),
      paste0("Baseline\n", metric_label),
      paste0("Interaction\n", psite_label),
      paste0("Interaction\n", metric_label)
    )
  )
  wide
}

robust_collision_ids <- intersect(
  same_genes("DS", "collision_score", "Resistance_baseline")$gene_id_clean,
  same_genes("DS", "collision_score", "Interaction")$gene_id_clean
)
robust_re_ids <- intersect(
  same_genes("RS", "ribosome_efficiency_score", "Resistance_baseline")$gene_id_clean,
  same_genes("RS", "ribosome_efficiency_score", "Interaction")$gene_id_clean
)

collision_wide <- make_matrix(
  robust_collision_ids,
  "DS",
  "collision_score",
  "DS p-site",
  "Collision score"
)
collision_wide <- collision_wide[order(-`Baseline\nDS p-site`)]

re_wide <- make_matrix(
  robust_re_ids,
  "RS",
  "ribosome_efficiency_score",
  "RS p-site",
  "Ribosome engagement"
)
re_wide[, ilf2_flag := gene_name == "ILF2"]
re_wide <- re_wide[order(ilf2_flag, -`Baseline\nRS p-site`)]
re_wide[, ilf2_flag := NULL]

to_long <- function(wide, cap) {
  value_cols <- setdiff(names(wide), c("gene_id_clean", "gene_name"))
  long <- melt(
    wide,
    id.vars = c("gene_id_clean", "gene_name"),
    measure.vars = value_cols,
    variable.name = "measure",
    value.name = "logFC"
  )
  long[, measure := factor(measure, levels = value_cols)]
  long[, gene_name := factor(gene_name, levels = rev(wide$gene_name))]
  long[, capped_logFC := pmax(pmin(logFC, cap), -cap)]
  long
}

plot_panel <- function(long, title, cap, show_legend = TRUE) {
  p <- ggplot(long, aes(x = measure, y = gene_name, fill = capped_logFC)) +
    geom_tile(color = "white", linewidth = 0.55) +
    geom_vline(xintercept = 2.5, linewidth = 0.9, color = "grey20") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-cap, cap),
      breaks = c(-cap, 0, cap),
      name = "limma\nlogFC"
    ) +
    scale_x_discrete(position = "top") +
    labs(title = title, x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 16),
      axis.text.x = element_text(size = 11, color = "black", lineheight = 0.95),
      axis.text.y = element_text(size = 10.5, color = "black"),
      panel.grid = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      plot.margin = margin(8, 12, 8, 8)
    )
  p
}

collision_long <- to_long(collision_wide, cap = 3)
re_long <- to_long(re_wide, cap = 5)

p_collision <- plot_panel(
  collision_long,
  "A. Robust collision genes",
  cap = 3,
  show_legend = TRUE
)
p_re <- plot_panel(
  re_long,
  "B. Robust ribosome engagement genes",
  cap = 5,
  show_legend = TRUE
)

ggsave(file.path(OUT_DIR, "robust_collision_genes_heatmap.png"), p_collision,
       width = 7.8, height = 8.2, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "robust_collision_genes_heatmap.pdf"), p_collision,
       width = 7.8, height = 8.2, bg = "white")

ggsave(file.path(OUT_DIR, "robust_ribosome_engagement_genes_heatmap.png"), p_re,
       width = 7.8, height = 4.6, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "robust_ribosome_engagement_genes_heatmap.pdf"), p_re,
       width = 7.8, height = 4.6, bg = "white")

combined <- p_collision / p_re + plot_layout(heights = c(2.1, 1))
ggsave(file.path(OUT_DIR, "robust_cross_framework_genes_heatmap_panels.png"), combined,
       width = 8.3, height = 12.5, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "robust_cross_framework_genes_heatmap_panels.pdf"), combined,
       width = 8.3, height = 12.5, bg = "white")

cat("Saved heatmaps to:\n", OUT_DIR, "\n", sep = "")
cat("Collision genes:", nrow(collision_wide), "\n")
cat("Ribosome engagement genes:", nrow(re_wide), "\n")
