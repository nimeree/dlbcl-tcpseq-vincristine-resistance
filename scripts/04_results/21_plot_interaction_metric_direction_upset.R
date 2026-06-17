# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(cowplot)
})

base_dir <- analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05")
results_dir <- file.path(base_dir, "Results")
out_dir <- file.path(base_dir, "Multi_metric_integration", "Pattern_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

contrast <- "Interaction"
lfc_cutoff <- 0.7
p_cutoff <- 0.05

sets <- list(
  "Scanning Up" = list(metric = "scanning_score", direction = "Up"),
  "Ribosome engagement Up" = list(metric = "ribosome_efficiency_score", direction = "Up"),
  "Collision Down" = list(metric = "collision_score", direction = "Down"),
  "Collision Up" = list(metric = "collision_score", direction = "Up")
)

read_metric <- function(metric) {
  path <- file.path(results_dir, metric, paste0(contrast, "_limma_all_genes.csv"))
  dat <- read.csv(path, check.names = FALSE)
  dat <- dat[!is.na(dat$gene_name) & dat$gene_name != "", ]
  dat$metric <- metric
  dat$direction <- ifelse(
    dat$P.Value < p_cutoff & dat$logFC >= lfc_cutoff, "Up",
    ifelse(dat$P.Value < p_cutoff & dat$logFC <= -lfc_cutoff, "Down", "NS")
  )
  dat
}

all_results <- do.call(rbind, lapply(unique(vapply(sets, `[[`, character(1), "metric")), read_metric))

gene_sets <- lapply(sets, function(x) {
  unique(all_results$gene_name[all_results$metric == x$metric & all_results$direction == x$direction])
})

all_genes <- sort(unique(unlist(gene_sets)))
membership <- data.frame(gene_name = all_genes, stringsAsFactors = FALSE)
for (set_name in names(gene_sets)) {
  membership[[set_name]] <- membership$gene_name %in% gene_sets[[set_name]]
}
membership$n_sets <- rowSums(membership[names(gene_sets)])
membership <- membership[membership$n_sets > 0, ]

combo_key <- apply(membership[names(gene_sets)], 1, function(x) paste(as.integer(x), collapse = ""))
combo_counts <- as.data.frame(table(combo_key), stringsAsFactors = FALSE)
names(combo_counts) <- c("combo_key", "count")
combo_counts <- combo_counts[combo_counts$count > 0, ]
combo_counts$n_sets <- rowSums(do.call(rbind, strsplit(combo_counts$combo_key, "")) == "1")
priority_keys <- c("1110", "0110", "1000", "0001")
combo_counts$priority <- match(combo_counts$combo_key, priority_keys)
combo_counts$priority[is.na(combo_counts$priority)] <- length(priority_keys) + 1
combo_counts <- combo_counts[order(combo_counts$priority, -combo_counts$n_sets, -combo_counts$count), ]
combo_counts$combo_id <- seq_len(nrow(combo_counts))
combo_counts$display_count <- combo_counts$count

matrix_dat <- do.call(rbind, lapply(seq_len(nrow(combo_counts)), function(i) {
  present <- strsplit(combo_counts$combo_key[i], "")[[1]] == "1"
  data.frame(
    combo_id = combo_counts$combo_id[i],
    set_name = names(gene_sets),
    present = present,
    stringsAsFactors = FALSE
  )
}))
matrix_dat$set_name <- factor(matrix_dat$set_name, levels = rev(names(gene_sets)))

set_size_dat <- data.frame(
  set_name = factor(names(gene_sets), levels = rev(names(gene_sets))),
  size = vapply(gene_sets, length, integer(1))
)

combo_counts$label <- factor(combo_counts$combo_id, levels = combo_counts$combo_id)
matrix_dat$label <- factor(matrix_dat$combo_id, levels = combo_counts$combo_id)

highlight_keys <- c(
  "Productive/prioritised translation" = paste(as.integer(names(gene_sets) %in% c(
    "Ribosome engagement Up", "Collision Down"
  )), collapse = ""),
  "Full-cycle prioritisation subset" = paste(as.integer(names(gene_sets) %in% c(
    "Scanning Up", "Ribosome engagement Up", "Collision Down"
  )), collapse = "")
)
combo_counts$collision_up_present <- substr(combo_counts$combo_key, 4, 4) == "1"
combo_counts$signature <- "Other"
combo_counts$signature[combo_counts$collision_up_present] <- "Collision stress"
combo_counts$signature[combo_counts$combo_key == "1000"] <- "Scanning bottleneck"
combo_counts$signature[combo_counts$combo_key == highlight_keys["Productive/prioritised translation"]] <- "Productive/prioritised translation"
combo_counts$signature[combo_counts$combo_key == highlight_keys["Full-cycle prioritisation subset"]] <- "Full-cycle prioritisation subset"

signature_cols <- c(
  "Scanning bottleneck" = "#D89000",
  "Collision stress" = "#7C3AED",
  "Productive/prioritised translation" = "#C0392B",
  "Full-cycle prioritisation subset" = "#8E44AD",
  "Other" = "#6B7280"
)

bar_plot <- ggplot(combo_counts, aes(label, display_count, fill = signature)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = display_count), vjust = -0.35, size = 3.2) +
  scale_fill_manual(values = signature_cols, name = NULL, guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(
    title = paste0("Interaction multi-metric overlaps (raw P < ", p_cutoff, ", |logFC| >= ", lfc_cutoff, ")"),
    subtitle = "Set sizes show simplified category totals; bars show exclusive intersections.",
    x = NULL,
    y = "Genes"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.2, color = "#374151", margin = margin(t = 3, b = 3)),
    legend.position = "top",
    legend.box = "vertical",
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(0.18, "in"),
    plot.margin = margin(4, 8, 4, 8)
  )

matrix_plot <- ggplot(matrix_dat, aes(label, set_name)) +
  geom_point(aes(color = present), size = 3.1) +
  geom_line(
    data = matrix_dat[matrix_dat$present, ],
    aes(group = label),
    color = "#374151",
    linewidth = 0.45
  ) +
  scale_color_manual(values = c(`TRUE` = "#111827", `FALSE` = "#D1D5DB"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 9),
    plot.margin = margin(0, 5, 2, 5)
  )

set_size_plot <- ggplot(set_size_dat, aes(size, set_name)) +
  geom_col(width = 0.62, fill = "#4B5563") +
  geom_text(aes(label = size), hjust = -0.08, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.32))) +
  coord_cartesian(clip = "off") +
  labs(x = "Set size", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 18, 2, 22)
  )

aligned <- align_plots(bar_plot, matrix_plot, align = "v", axis = "lr")
bar_aligned <- aligned[[1]]
matrix_aligned <- aligned[[2]]

bottom_panel <- plot_grid(
  set_size_plot,
  matrix_aligned,
  ncol = 2,
  rel_widths = c(1.55, 5.25),
  align = "h",
  axis = "tb"
)

top_panel <- plot_grid(
  ggdraw(),
  bar_aligned,
  ncol = 2,
  rel_widths = c(1.55, 5.25),
  align = "h",
  axis = "tb"
)

combined <- plot_grid(
  top_panel,
  bottom_panel,
  ncol = 1,
  rel_heights = c(3.1, 1.65),
  align = "v",
  axis = "lr"
)

png_path <- file.path(out_dir, "interaction_metric_direction_upset.png")
pdf_path <- file.path(out_dir, "interaction_metric_direction_upset.pdf")

ggsave(png_path, combined, width = 11.6, height = 5.4, dpi = 300, bg = "white")
ggsave(pdf_path, combined, width = 11.6, height = 5.4, bg = "white")

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("Interaction genes in union: ", nrow(membership))
message("Scanning bottleneck genes: ", sum(membership[["Scanning Up"]]))
message("Collision stress genes: ", sum(membership[["Collision Up"]]))
message("Productive/prioritised translation genes: ", sum(membership[["Ribosome engagement Up"]] & membership[["Collision Down"]]))
message("Full-cycle prioritisation subset genes: ", sum(membership[["Scanning Up"]] & membership[["Ribosome engagement Up"]] & membership[["Collision Down"]]))
