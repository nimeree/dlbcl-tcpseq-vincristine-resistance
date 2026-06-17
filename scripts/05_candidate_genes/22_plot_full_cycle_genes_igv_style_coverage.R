# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- analysis_path()
index_dir <- file.path(base_dir, "Translation_indexes_fixed")
out_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05", "Multi_metric_integration", "Full_cycle_IGV_style")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes <- c("MAPKBP1", "SEC24C", "TRA2A")

psite_file <- file.path(index_dir, "transcript_psite_matrix_long_ALL_samples.csv")
metric_file <- file.path(index_dir, "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")

psite <- fread(psite_file)
tx_map <- unique(fread(metric_file, select = c("transcript", "gene_id_clean", "gene_name"))[gene_name %in% genes])
psite <- merge(psite, tx_map, by = "transcript")

top_tx <- psite[, .(
  total_psites = sum(psite_count, na.rm = TRUE),
  n_positions = uniqueN(codon_pos)
), by = .(gene_name, transcript)][
  order(gene_name, -total_psites)
][, .SD[1], by = gene_name]

fwrite(top_tx, file.path(out_dir, "full_cycle_genes_top_transcripts_for_coverage.csv"))

plot_dat <- merge(psite, top_tx[, .(gene_name, transcript)], by = c("gene_name", "transcript"))

plot_dat[, condition := fifelse(grepl("SU8R-DMSO|SU8-R-DMSO", sample), "Resistant DMSO",
                         fifelse(grepl("SU8R-Vin|SU8-R-Vin", sample), "Resistant VCR",
                         fifelse(grepl("SU8-DMSO", sample), "Sensitive DMSO",
                         fifelse(grepl("SU8-Vin", sample), "Sensitive VCR", NA_character_))))]

plot_dat <- plot_dat[!is.na(condition)]
plot_dat[, condition := factor(condition, levels = c("Sensitive DMSO", "Sensitive VCR", "Resistant DMSO", "Resistant VCR"))]
plot_dat[, fraction := factor(fraction, levels = c("SSU", "RS", "DS"))]

mean_cov <- plot_dat[, .(
  mean_cpm = mean(psite_cpm, na.rm = TRUE),
  mean_count = mean(psite_count, na.rm = TRUE)
), by = .(gene_name, transcript, fraction, condition, codon_pos)]

all_positions <- mean_cov[, .(codon_pos = seq(min(codon_pos), max(codon_pos))), by = .(gene_name, transcript, fraction, condition)]
mean_cov <- merge(all_positions, mean_cov, by = c("gene_name", "transcript", "fraction", "condition", "codon_pos"), all.x = TRUE)
mean_cov[is.na(mean_cpm), mean_cpm := 0]
mean_cov[is.na(mean_count), mean_count := 0]

track_summary <- mean_cov[, .(
  total_mean_cpm = sum(mean_cpm),
  max_mean_cpm = max(mean_cpm),
  covered_positions = sum(mean_cpm > 0)
), by = .(gene_name, transcript, fraction, condition)]
fwrite(track_summary, file.path(out_dir, "full_cycle_genes_igv_style_track_summary.csv"))

make_gene_plot <- function(gene) {
  d <- mean_cov[gene_name == gene]
  tx <- unique(d$transcript)

  d[, track := paste(fraction, condition, sep = " | ")]
  track_levels <- as.vector(outer(
    rev(levels(d$fraction)),
    rev(levels(d$condition)),
    paste,
    sep = " | "
  ))
  d[, track := factor(track, levels = track_levels)]

  p <- ggplot(d, aes(codon_pos, mean_cpm)) +
    geom_col(width = 1, fill = "#3F4A5A", color = NA) +
    facet_grid(track ~ ., scales = "free_y", switch = "y") +
    labs(
      title = paste0(gene, " IGV-style P-site coverage"),
      subtitle = paste0("Transcript: ", tx, " | replicate mean P-site CPM | transcript codon coordinates"),
      x = "Codon position on transcript",
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing.y = unit(0.06, "in"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, size = 7.2, color = "black"),
      strip.background = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "#4B5563"),
      plot.margin = margin(5, 5, 5, 5)
    )

  png_path <- file.path(out_dir, paste0(gene, "_igv_style_psite_coverage.png"))
  pdf_path <- file.path(out_dir, paste0(gene, "_igv_style_psite_coverage.pdf"))
  ggsave(png_path, p, width = 8.6, height = 7.2, dpi = 300, bg = "white")
  ggsave(pdf_path, p, width = 8.6, height = 7.2, bg = "white")
  list(png = png_path, pdf = pdf_path)
}

paths <- lapply(genes, make_gene_plot)

message("Saved IGV-style coverage plots to: ", out_dir)
print(top_tx)
print(track_summary[order(gene_name, fraction, condition)])
