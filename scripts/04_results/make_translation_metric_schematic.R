# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

out_dir <- analysis_path("Translation_metric_schematic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

box <- function(xmin, xmax, ymin, ymax, label, fill, color = "#253041", size = 3.4, fontface = "plain") {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      fill = fill, color = color, linewidth = 0.45),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
      label = label, size = size, fontface = fontface, lineheight = 0.92)
  )
}

arrow_seg <- function(x, xend, y, yend = y, color = "#253041") {
  annotate("segment", x = x, xend = xend, y = y, yend = yend,
    arrow = arrow(length = unit(0.16, "in"), type = "closed"),
    linewidth = 0.55, color = color)
}

metric_rows <- data.frame(
  y = c(6.25, 4.75, 3.25, 1.75),
  title = c("Scanning score", "Collision score", "Ribosome efficiency score", "Protein output score"),
  formula = c(
    "(SSU 5'UTR density + 1e-3) / (RS CDS density + 1e-3)",
    "(DS CDS density + 1e-3) / (RS CDS density + 1e-3)",
    "log2((RS core CPM + 1) / (RNA baseline CPM + 1))",
    "log2((RNA baseline CPM + 1) x (RS core CPM + 1))"
  ),
  meaning = c(
    "Relative 5'UTR small-subunit loading; high values suggest scanning/initiation bottleneck relative to 80S signal.",
    "Disome enrichment relative to 80S signal; high values suggest more collided or stalled ribosomes.",
    "80S ribosome occupancy relative to RNA abundance; high values mean more ribosome loading per transcript.",
    "Composite output proxy; high values require both RNA abundance and 80S ribosome occupancy."
  ),
  fill = c("#E8F1FB", "#FBEAEA", "#EAF5EC", "#FFF4D9")
)

p <- ggplot() +
  coord_cartesian(xlim = c(0, 14.4), ylim = c(0, 8), clip = "off") +
  theme_void(base_size = 12) +
  annotate("text", x = 0.05, y = 7.65, hjust = 0, label = "TCP-seq fractions", fontface = "bold", size = 4.3) +
  annotate("text", x = 4.2, y = 7.65, hjust = 0, label = "Derived limma metric matrices", fontface = "bold", size = 4.3) +
  annotate("text", x = 11.15, y = 7.65, hjust = 0, label = "Fraction-specific DESeq2", fontface = "bold", size = 4.3) +
  box(0.25, 2.55, 6.55, 7.15, "SSU fraction\nsmall subunit reads", "#D8E9FF", size = 3.15)[[1]] +
  box(0.25, 2.55, 6.55, 7.15, "SSU fraction\nsmall subunit reads", "#D8E9FF", size = 3.15)[[2]] +
  box(0.25, 2.55, 5.45, 6.05, "RS fraction\n80S / monosome reads", "#DDF0DF", size = 3.15)[[1]] +
  box(0.25, 2.55, 5.45, 6.05, "RS fraction\n80S / monosome reads", "#DDF0DF", size = 3.15)[[2]] +
  box(0.25, 2.55, 4.35, 4.95, "DS fraction\ndisome reads", "#F6DCDC", size = 3.15)[[1]] +
  box(0.25, 2.55, 4.35, 4.95, "DS fraction\ndisome reads", "#F6DCDC", size = 3.15)[[2]] +
  box(0.25, 2.55, 2.95, 3.55, "RNA baseline\ncell-line mRNA CPM", "#EFE7D5", size = 3.15)[[1]] +
  box(0.25, 2.55, 2.95, 3.55, "RNA baseline\ncell-line mRNA CPM", "#EFE7D5", size = 3.15)[[2]] +
  arrow_seg(2.75, 3.75, 6.85) +
  arrow_seg(2.75, 3.75, 5.75) +
  arrow_seg(2.75, 3.75, 4.65) +
  arrow_seg(2.75, 3.75, 3.25) +
  annotate("text", x = 1.4, y = 2.15, label = "Transcript-level metrics were\ncollapsed to gene level by median\nbefore limma.", size = 2.95, lineheight = 0.95, color = "#374151") +
  annotate("rect", xmin = 3.85, xmax = 10.75, ymin = 0.75, ymax = 7.15, fill = "#FFFFFF", color = "#CBD5E1", linewidth = 0.45) +
  lapply(seq_len(nrow(metric_rows)), function(i) {
    y <- metric_rows$y[i]
    list(
      annotate("rect", xmin = 4.15, xmax = 10.45, ymin = y - 0.53, ymax = y + 0.53,
        fill = metric_rows$fill[i], color = "#334155", linewidth = 0.35),
      annotate("text", x = 4.38, y = y + 0.27, hjust = 0, label = metric_rows$title[i],
        fontface = "bold", size = 3.25),
      annotate("text", x = 4.38, y = y, hjust = 0, label = metric_rows$formula[i],
        family = "mono", size = 2.55),
      annotate("text", x = 4.38, y = y - 0.28, hjust = 0, label = metric_rows$meaning[i],
        size = 2.55, lineheight = 0.9, color = "#374151")
    )
  }) +
  annotate("text", x = 7.3, y = 0.35,
    label = "limma asks: does this continuous score change by condition?\nContrasts: VCR sensitive, VCR resistant, baseline resistance, interaction.",
    size = 3.05, lineheight = 0.95, fontface = "bold", color = "#111827") +
  arrow_seg(10.95, 11.55, 5.75) +
  box(11.75, 14.05, 6.2, 6.85, "SSU count matrix", "#D8E9FF", size = 3.05)[[1]] +
  box(11.75, 14.05, 6.2, 6.85, "SSU count matrix", "#D8E9FF", size = 3.05)[[2]] +
  box(11.75, 14.05, 5.25, 5.9, "RS count matrix", "#DDF0DF", size = 3.05)[[1]] +
  box(11.75, 14.05, 5.25, 5.9, "RS count matrix", "#DDF0DF", size = 3.05)[[2]] +
  box(11.75, 14.05, 4.3, 4.95, "DS count matrix", "#F6DCDC", size = 3.05)[[1]] +
  box(11.75, 14.05, 4.3, 4.95, "DS count matrix", "#F6DCDC", size = 3.05)[[2]] +
  annotate("text", x = 12.9, y = 3.45,
    label = "DESeq2 asks: does a gene's\nraw count change within one fraction?\n\nSSU up = more small-subunit signal\nRS up = more 80S occupancy\nDS up = more disome signal\n\nIt does not directly test ratios\nsuch as DS/RS or SSU/RS.",
    size = 3.0, lineheight = 0.95, color = "#111827") +
  annotate("rect", xmin = 0.15, xmax = 14.1, ymin = 0.05, ymax = 0.18, fill = "#253041", color = NA) +
  annotate("text", x = 7.1, y = -0.12,
    label = "Practical interpretation: fraction-specific DESeq2 tells which ribosome population changes; limma metric scores tell whether relationships between populations change.",
    size = 3.1, fontface = "bold", color = "#253041")

png_path <- file.path(out_dir, "translation_metric_schematic_limma_vs_fraction_deseq.png")
pdf_path <- file.path(out_dir, "translation_metric_schematic_limma_vs_fraction_deseq.pdf")
ggsave(png_path, p, width = 14, height = 8, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(pdf_path, p, width = 14, height = 8, bg = "white", limitsize = FALSE)

writeLines(c(
  "Translation metric schematic",
  "",
  "Ribosome efficiency score = log2((RS core CPM + 1) / (RNA baseline CPM + 1)).",
  "Protein output score = log2((RNA baseline CPM + 1) x (RS core CPM + 1)).",
  "Collision score = (DS core density per kb + 1e-3) / (RS core density per kb + 1e-3).",
  "Scanning score = (SSU 5'UTR density per kb + 1e-3) / (RS core density per kb + 1e-3).",
  "",
  "limma was run on gene-level continuous metric matrices after transcript-level values were collapsed by median.",
  "Fraction-specific DESeq2 was run on count matrices for SSU, RS, and DS separately; it tests abundance within a fraction rather than cross-fraction ratios."
), file.path(out_dir, "translation_metric_schematic_notes.txt"))

message("Saved schematic: ", png_path)
