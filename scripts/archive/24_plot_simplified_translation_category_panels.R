# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
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
  all_results$P.Value < p_cutoff & all_results$logFC >= lfc_cutoff,
  "Up",
  ifelse(all_results$P.Value < p_cutoff & all_results$logFC <= -lfc_cutoff, "Down", "NS")
)

get_genes <- function(metric, direction) {
  unique(all_results$gene_name[all_results$metric == metric & all_results$direction == direction])
}

category_sets <- list(
  "A" = list(
    title = "A. Scanning bottleneck",
    subtitle = "Scanning Up",
    genes = sort(get_genes("scanning_score", "Up")),
    order_metric = "scanning_score",
    decreasing = TRUE,
    file_stub = "A_scanning_bottleneck_heatmap"
  ),
  "B" = list(
    title = "B. Full-cycle translation",
    subtitle = "Scanning Down + Ribosome engagement Up + Collision Down",
    genes = sort(Reduce(intersect, list(
      get_genes("scanning_score", "Down"),
      get_genes("ribosome_efficiency_score", "Up"),
      get_genes("collision_score", "Down")
    ))),
    order_metric = "ribosome_efficiency_score",
    decreasing = TRUE,
    file_stub = "B_full_cycle_translation_heatmap"
  ),
  "C" = list(
    title = "C. Collision stress",
    subtitle = "Collision Up",
    genes = sort(get_genes("collision_score", "Up")),
    order_metric = "collision_score",
    decreasing = TRUE,
    file_stub = "C_collision_stress_heatmap"
  )
)

heatmap_scale <- scale_fill_gradient2(
  low = "#2C7BB6", mid = "white", high = "#D7191C",
  midpoint = 0, limits = c(-3, 3), name = "limma\nlogFC"
)

make_panel <- function(panel) {
  genes <- panel$genes
  dat <- all_results[all_results$gene_name %in% genes, c("gene_name", "metric", "logFC")]
  dat$metric_label <- metrics[dat$metric]
  dat$metric_label <- factor(dat$metric_label, levels = unname(metrics))
  dat$logFC_cap <- pmax(pmin(dat$logFC, 3), -3)

  ord <- dat[dat$metric == panel$order_metric, ]
  gene_order <- genes[order(ord$logFC[match(genes, ord$gene_name)], decreasing = panel$decreasing, na.last = TRUE)]
  dat$gene_name <- factor(dat$gene_name, levels = rev(unique(gene_order)))

  height <- max(3.2, min(10.5, 1.4 + 0.17 * length(genes)))
  width <- 5.9

  p <- ggplot(dat, aes(metric_label, gene_name, fill = logFC_cap)) +
    geom_tile(color = "white", linewidth = 0.35) +
    heatmap_scale +
    labs(
      title = panel$title,
      subtitle = paste0(panel$subtitle, " (n = ", length(genes), ")"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 7.5, color = "black"),
      axis.text.x = element_text(size = 9, color = "black", angle = 35, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 9.5, color = "#374151", hjust = 0),
      plot.margin = margin(6, 8, 6, 8),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  png_path <- file.path(out_dir, paste0(panel$file_stub, ".png"))
  pdf_path <- file.path(out_dir, paste0(panel$file_stub, ".pdf"))
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_path, p, width = width, height = height, bg = "white")
  message("Saved: ", png_path)
  message(panel$title, " genes: ", length(genes))
  invisible(p)
}

invisible(lapply(category_sets, make_panel))
