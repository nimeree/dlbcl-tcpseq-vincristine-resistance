# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

base_dir <- analysis_path()
out_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")

p_cut <- 0.05
lfc_cut <- 0.7
top_n_each <- 10

make_volcano <- function(dt, title, out_png) {
  df <- copy(dt)
  df <- df[!is.na(P.Value) & !is.na(logFC)]
  df[, neglog10_p := -log10(P.Value)]
  df[, direction := fifelse(
    P.Value < p_cut & logFC >= lfc_cut, "Up",
    fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS")
  )]
  n_up <- df[direction == "Up", .N]
  n_down <- df[direction == "Down", .N]

  lab <- rbindlist(list(
    df[direction == "Up"][order(P.Value)][1:min(top_n_each, .N)],
    df[direction == "Down"][order(P.Value)][1:min(top_n_each, .N)]
  ), fill = TRUE)

  p <- ggplot(df, aes(logFC, neglog10_p)) +
    geom_point(aes(color = direction), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c(Up = "#D7191C", Down = "#2C7BB6", NS = "grey72")) +
    geom_hline(yintercept = -log10(p_cut), linetype = "dotted", linewidth = 0.7) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", linewidth = 0.7) +
    ggrepel::geom_text_repel(
      data = lab,
      aes(label = gene_name),
      size = 2.8,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25
    ) +
    annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.25,
             label = paste0("Up: ", n_up), color = "#D7191C", size = 3.2, fill = "white") +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.25,
             label = paste0("Down: ", n_down), color = "#2C7BB6", size = 3.2, fill = "white") +
    labs(
      title = title,
      subtitle = paste0(
        "P-site offset gene counts; limma baseline = Resistant DMSO - Sensitive DMSO; ",
        "raw P < ", p_cut, ", |logFC| >= ", lfc_cut
      ),
      x = "limma baseline logFC",
      y = "-log10(P value)",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )

  ggsave(out_png, p, width = 7.2, height = 5.5, dpi = 300, bg = "white")
}

summary_rows <- list()
for (frac in c("SSU", "RS", "DS")) {
  frac_dir <- file.path(out_dir, paste0("Fraction_", frac))
  in_file <- file.path(frac_dir, paste0("Resistance_baseline_", frac, "_psite_limma_all_genes.csv"))
  if (!file.exists(in_file)) stop("Missing baseline table: ", in_file)

  dt <- fread(in_file)
  dt[, significant_rawP0.05_lfc0.7 := !is.na(P.Value) & P.Value < p_cut & abs(logFC) >= lfc_cut]
  dt[, direction := fifelse(
    P.Value < p_cut & logFC >= lfc_cut, "Up",
    fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS")
  )]

  fwrite(dt[significant_rawP0.05_lfc0.7 == TRUE],
         file.path(frac_dir, paste0("Resistance_baseline_", frac, "_psite_limma_sig_rawP0.05_lfc0.7.csv")))
  fwrite(dt[direction == "Up"],
         file.path(frac_dir, paste0("Resistance_baseline_", frac, "_psite_limma_sig_up_rawP0.05_lfc0.7.csv")))
  fwrite(dt[direction == "Down"],
         file.path(frac_dir, paste0("Resistance_baseline_", frac, "_psite_limma_sig_down_rawP0.05_lfc0.7.csv")))

  plot_file <- file.path(frac_dir, paste0("Resistance_baseline_", frac, "_psite_limma_volcano.png"))
  make_volcano(
    dt,
    title = paste0(frac, " P-site limma: baseline resistance"),
    out_png = plot_file
  )

  summary_rows[[frac]] <- data.table(
    fraction = frac,
    genes_tested = nrow(dt),
    n_sig = dt[significant_rawP0.05_lfc0.7 == TRUE, .N],
    n_up = dt[direction == "Up", .N],
    n_down = dt[direction == "Down", .N],
    plot = plot_file
  )
}

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(out_dir, "psite_fraction_limma_baseline_plot_summary.csv"))

cat("\nGenerated baseline P-site fraction limma volcano plots:\n")
print(summary_dt)
