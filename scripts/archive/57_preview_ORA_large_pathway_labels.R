# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# Preview one Word-friendly ORA plot with larger pathway labels.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(stringr)
})

in_file <- analysis_path("Psite_fraction_limma_lfc0.7_rawP0.05", "GO_BP_ORA_Baseline_Combined", "Tables", "SSU_terms_shown_in_plot.csv")
out_base <- analysis_path("Psite_fraction_limma_lfc0.7_rawP0.05", "GO_BP_ORA_Baseline_Combined", "Plots", "SSU", "SSU_GO_BP_ORA_baseline_combined_large_labels_preview")

dt <- fread(in_file)
dt <- dt[order(p_value, -intersection_size)]
dt[, term_label := stringr::str_wrap(term_name, width = 34)]
dt[, neglog10_gscs := -log10(p_value)]

p <- ggplot(dt, aes(
  x = neglog10_gscs,
  y = forcats::fct_reorder(term_label, neglog10_gscs)
)) +
  geom_segment(aes(
    x = 0, xend = neglog10_gscs,
    y = forcats::fct_reorder(term_label, neglog10_gscs),
    yend = forcats::fct_reorder(term_label, neglog10_gscs)
  ), linewidth = 0.9, color = "grey76") +
  geom_point(aes(size = intersection_size, color = mean_hit_log2FC), alpha = 0.95) +
  scale_color_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
  scale_size_continuous(range = c(5, 10)) +
  labs(
    title = "SSU P-site count baseline GO:BP ORA",
    subtitle = "Baseline; GO Biological Process only; Up + Down combined; g:SCS < 0.05; term size <= 500; hit genes >= 5",
    x = "-log10(g:SCS corrected p-value)",
    y = NULL,
    size = "Hit genes",
    color = "Mean hit\nlogFC"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 13),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(color = "black", size = 15, lineheight = 0.95),
    panel.grid.major.y = element_blank(),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(t = 10, r = 18, b = 10, l = 16)
  )

ggsave(paste0(out_base, ".png"), p, width = 13.5, height = 7.2, dpi = 300, limitsize = FALSE, bg = "white")
ggsave(paste0(out_base, ".pdf"), p, width = 13.5, height = 7.2, limitsize = FALSE, bg = "white")

cat("Saved preview:\n", paste0(out_base, ".png"), "\n", sep = "")
