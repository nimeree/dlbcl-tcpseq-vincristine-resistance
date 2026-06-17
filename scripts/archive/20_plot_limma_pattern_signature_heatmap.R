# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

base_dir <- analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05")
results_dir <- file.path(base_dir, "Results")
out_dir <- file.path(base_dir, "Multi_metric_integration", "Pattern_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

contrast <- "Interaction"
lfc_cutoff <- 0.7
p_cutoff <- 0.05

metrics <- c(
  scanning_score = "Scanning",
  ribosome_efficiency_score = "Ribosome engagement",
  collision_score = "Collision"
)

read_metric <- function(metric) {
  path <- file.path(results_dir, metric, paste0(contrast, "_limma_all_genes.csv"))
  dat <- read.csv(path, check.names = FALSE)
  dat$metric <- metric
  dat
}

all_results <- do.call(rbind, lapply(names(metrics), read_metric))
all_results <- all_results[!is.na(all_results$gene_name) & all_results$gene_name != "", ]
all_results$direction <- ifelse(
  all_results$P.Value < p_cutoff & all_results$logFC >= lfc_cutoff, "Up",
  ifelse(all_results$P.Value < p_cutoff & all_results$logFC <= -lfc_cutoff, "Down", "NS")
)

get_genes <- function(metric, direction) {
  unique(all_results$gene_name[all_results$metric == metric & all_results$direction == direction])
}

pattern1_genes <- Reduce(intersect, list(
  get_genes("ribosome_efficiency_score", "Up"),
  get_genes("scanning_score", "Down"),
  get_genes("collision_score", "Down")
))

scanning_bottleneck_genes <- get_genes("scanning_score", "Up")
collision_stress_genes <- get_genes("collision_score", "Up")

pattern_sets <- list(
  "Scanning bottleneck" = sort(scanning_bottleneck_genes),
  "Full-cycle translation" = sort(pattern1_genes),
  "Collision stress" = sort(collision_stress_genes)
)

heatmap_dat <- do.call(rbind, lapply(names(pattern_sets), function(pattern_label) {
  genes <- pattern_sets[[pattern_label]]
  subset <- all_results[all_results$gene_name %in% genes, c("gene_name", "metric", "logFC")]
  subset$pattern_label <- pattern_label
  subset
}))

heatmap_dat$metric_label <- metrics[heatmap_dat$metric]
heatmap_dat$metric_label <- factor(heatmap_dat$metric_label, levels = unname(metrics))
heatmap_dat$logFC_cap <- pmax(pmin(heatmap_dat$logFC, 3), -3)
heatmap_dat$pattern_strip <- factor(
  heatmap_dat$pattern_label,
  levels = names(pattern_sets)
)

gene_order <- unlist(lapply(names(pattern_sets), function(pattern_label) {
  genes <- pattern_sets[[pattern_label]]
  order_metric <- if (pattern_label == "Collision stress") "collision_score" else if (pattern_label == "Scanning bottleneck") "scanning_score" else "ribosome_efficiency_score"
  sub <- heatmap_dat[heatmap_dat$pattern_label == pattern_label & heatmap_dat$metric == order_metric, ]
  decreasing <- pattern_label != "Collision stress"
  genes[order(sub$logFC[match(genes, sub$gene_name)], decreasing = decreasing, na.last = TRUE)]
}))
heatmap_dat$gene_name_plot <- paste(heatmap_dat$pattern_label, heatmap_dat$gene_name, sep = "__")
gene_levels <- unlist(lapply(names(pattern_sets), function(pattern_label) {
  paste(pattern_label, rev(gene_order[gene_order %in% pattern_sets[[pattern_label]]]), sep = "__")
}))
heatmap_dat$gene_name_plot <- factor(heatmap_dat$gene_name_plot, levels = unique(gene_levels))

heatmap_scale <- scale_fill_gradient2(
  low = "#2C7BB6", mid = "white", high = "#D7191C",
  midpoint = 0, limits = c(-3, 3), name = "limma\nlogFC"
)

combined <- ggplot(heatmap_dat, aes(metric_label, gene_name_plot, fill = logFC_cap)) +
  geom_tile(color = "white", linewidth = 0.35) +
  heatmap_scale +
  facet_grid(pattern_strip ~ ., scales = "free_y", switch = "y") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  labs(
    title = paste0("Interaction multi-metric signatures (raw P < ", p_cutoff, ", |logFC| >= ", lfc_cutoff, ")"),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.25, "in"),
    strip.placement = "outside",
    strip.background.y = element_rect(fill = "#F1F3F4", color = NA),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 8.5, color = "#222222"),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_text(size = 9, color = "black", angle = 35, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

png_path <- file.path(out_dir, "pattern1_pattern3_metric_logfc_heatmap.png")
pdf_path <- file.path(out_dir, "pattern1_pattern3_metric_logfc_heatmap.pdf")

ggsave(png_path, combined, width = 7.2, height = 14.2, dpi = 300, bg = "white")
ggsave(pdf_path, combined, width = 7.2, height = 14.2, bg = "white")

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("Scanning bottleneck genes: ", length(scanning_bottleneck_genes))
message("Full-cycle translation genes: ", length(pattern1_genes))
message("Collision stress genes: ", length(collision_stress_genes))
