# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)
library(ggplot2)
library(patchwork)

out_dir <- analysis_path("LongRead_TRA2A", "TRA2A_protein_isoform_comparison")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

seqs <- c(
  "TRA2A-201" = "MSDVEENNFEGRESRSQSKSPTGTPARVKSESRSGSRSPSRVSKHSESHSRSRSKSRSRSRRHSHRRYTRSRSHSHSHRRRSRSRSYTPEYRRRRSRSHSPMSNRRRHTGSRANPDPNTCLGVFGLSLYTTERDLREVFSRYGPLSGVNVVYDQRTGRSRGFAFVYFERIDDSKEAMERANGMELDGRRIRVDYSITKRAHTPTPGIYMGRPTHSGGGGGGGGGGGGGGGGRRRDSYYDRGYDRGYDRYEDYDYRYRRRSPSPYYSRYRSRSRSRSYSPRRY",
  "TRA2A-202" = "MSNRRRHTGSRANPDPNTCLGVFGLSLYTTERDLREVFSRYGPLSGVNVVYDQRTGRSRGFAFVYFERIDDSKEAMERANGMELDGRRIRVDYSITKRAHTPTPGIYMGRPTHSGGGGGGGGGGGGGGGGRRRDSYYDRGYDRGYDRYEDYDYRYRRSPSPYYSRYRSRSRSRSYSPRRY",
  "TRA2A-212" = "MSNRRRHTGSRANPDPNTCLGVFGLSLYTTERDLREVFSRYGPLSGVNVVYDQRTGRSRGFAFVYFERIDDSKEAMERANGMELDGRRIRVDYSITKRAHTPTPGIYMGRPTHSGGGGGGGGGGGGGGGGRRRDSYYDRGYDRGYDRYEDYDYRYRRRSPSPYYSRYRSRSRSRSYSPRRY"
)

tx_ids <- c(
  "TRA2A-201" = "ENST00000297071.9",
  "TRA2A-202" = "ENST00000392502.8",
  "TRA2A-212" = "ENST00000621813.4"
)

protein_ids <- c(
  "TRA2A-201" = "ENSP00000297071.4",
  "TRA2A-202" = "ENSP00000376290.4",
  "TRA2A-212" = "ENSP00000480822.1"
)

iso_levels <- c("TRA2A-201", "TRA2A-212", "TRA2A-202")

diff_202_212 <- data.table(
  pos_202 = seq_len(nchar(seqs[["TRA2A-202"]])),
  aa_202 = strsplit(seqs[["TRA2A-202"]], "")[[1]]
)
diff_202_212[, aa_212_same_index := strsplit(seqs[["TRA2A-212"]], "")[[1]][seq_len(.N)]]
first_diff <- diff_202_212[aa_202 != aa_212_same_index][1]

bars <- data.table(
  isoform = factor(
    c(
      rep("TRA2A-201", 4),
      rep("TRA2A-212", 3),
      rep("TRA2A-202", 3)
    ),
    levels = rev(iso_levels)
  ),
  domain = c(
    "RS1", "Linker", "RRM", "RS2",
    "Linker", "RRM", "RS2",
    "Linker", "RRM", "RS2"
  ),
  start = c(
    1, 102, 115, 195,
    102, 115, 195,
    102, 115, 195
  ),
  end = c(
    101, 114, 194, 282,
    114, 194, 282,
    114, 194, 281
  )
)

label_dt <- data.table(
  isoform = factor(iso_levels, levels = rev(iso_levels)),
  x = c(282, 282, 281),
  label = c(
    "282 aa | TRA2A-201\nENSP00000297071.4",
    "181 aa | TRA2A-212\nENSP00000480822.1",
    "180 aa | TRA2A-202\nENSP00000376290.4"
  )
)

plot_main <- ggplot(bars, aes(y = isoform)) +
  geom_rect(aes(xmin = start, xmax = end, ymin = as.numeric(isoform) - 0.23, ymax = as.numeric(isoform) + 0.23, fill = domain), color = "white", linewidth = 0.4) +
  geom_vline(xintercept = c(101.5, 114.5, 194.5), linetype = "dashed", color = "#4B5563", linewidth = 0.35) +
  geom_text(data = label_dt, aes(x = x + 4, y = isoform, label = label), inherit.aes = FALSE, hjust = 0, size = 3.1, lineheight = 0.95) +
  annotate("segment", x = 1, xend = 101, y = 3.28, yend = 3.28, linewidth = 0.35, arrow = arrow(ends = "both", length = unit(0.06, "in"))) +
  annotate("text", x = 51, y = 3.42, label = "RS1 / N-terminal 101 aa absent from 202/212", size = 3.25, fontface = "bold") +
  annotate("text", x = 51, y = 3, label = "RS1", size = 3.8, fontface = "bold", color = "white") +
  annotate("text", x = 108, y = 2.72, label = "linker", angle = 90, size = 2.6, color = "#374151") +
  annotate("text", x = 154.5, y = 3, label = "RRM", size = 3.8, fontface = "bold", color = "white") +
  annotate("text", x = 238.5, y = 3, label = "RS2", size = 3.8, fontface = "bold", color = "white") +
  scale_x_continuous(limits = c(1, 340), breaks = c(1, 50, 101, 115, 150, 194, 250, 282), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(
    "RS1" = "#6B7280",
    "Linker" = "#F2C94C",
    "RRM" = "#2F80ED",
    "RS2" = "#58A55C"
  )) +
  labs(
    title = "TRA2A isoform switch removes the N-terminal RS1 domain",
    subtitle = "TRA2A-202/212 retain the RRM and RS2 regions but lack the TRA2A-201 N-terminal RS1 segment",
    x = "Amino acid coordinate aligned to TRA2A-201",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(24, 12, 10, 8)
  )

zoom_seq <- data.table(
  isoform = factor(c("TRA2A-201 / 212", "TRA2A-202"), levels = rev(c("TRA2A-201 / 212", "TRA2A-202"))),
  sequence = c("YEDYDYRYRRRSPSPYYSRYRSRSRS", "YEDYDYRYRRSPSPYYSRYRSRSRS"),
  x = 1
)

zoom_letters <- rbindlist(lapply(seq_len(nrow(zoom_seq)), function(i) {
  data.table(
    isoform = zoom_seq$isoform[i],
    x = seq_len(nchar(zoom_seq$sequence[i])),
    aa = strsplit(zoom_seq$sequence[i], "")[[1]]
  )
}))
zoom_letters[, highlight := isoform == "TRA2A-201 / 212" & x == 10]

plot_zoom <- ggplot(zoom_letters, aes(x, isoform)) +
  geom_tile(aes(fill = highlight), width = 0.95, height = 0.72, color = "white") +
  geom_text(aes(label = aa), family = "mono", size = 4) +
  scale_fill_manual(values = c("TRUE" = "#F59E0B", "FALSE" = "#E5E7EB"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = "TRA2A-202 differs from TRA2A-212 by one arginine in the C-terminal R-rich region",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

combined <- plot_main / plot_zoom + plot_layout(heights = c(2.4, 1))

png_path <- file.path(out_dir, "TRA2A_201_202_212_protein_isoform_comparison.png")
pdf_path <- file.path(out_dir, "TRA2A_201_202_212_protein_isoform_comparison.pdf")
ggsave(png_path, combined, width = 11.5, height = 7.8, dpi = 300, bg = "white")
ggsave(pdf_path, combined, width = 11.5, height = 7.8, bg = "white")

summary_dt <- data.table(
  isoform = names(seqs),
  transcript_id = tx_ids[names(seqs)],
  protein_id = protein_ids[names(seqs)],
  aa_length = nchar(seqs),
  interpretation = c(
    "Full-length Ensembl protein isoform with RS1, RRM, and RS2 regions",
    "Shorter protein lacking RS1/N-terminal 101 aa; retains RRM/RS2; differs from 212 by one fewer arginine in C-terminal R-rich region",
    "Shorter protein corresponding to TRA2A-201 amino acids 102-282; lacks RS1 but retains RRM/RS2"
  )
)
fwrite(summary_dt, file.path(out_dir, "TRA2A_201_202_212_protein_isoform_summary.csv"))

writeLines(c(
  "TRA2A protein isoform comparison generated from Ensembl GRCh38 release 114 CDS-derived amino acid sequences.",
  "TRA2A-201: ENST00000297071.9 / ENSP00000297071.4, 282 aa.",
  "TRA2A-202: ENST00000392502.8 / ENSP00000376290.4, 180 aa.",
  "TRA2A-212: ENST00000621813.4 / ENSP00000480822.1, 181 aa.",
  "Domain schematic used in the figure: RS1 1-101; RRM 115-194; RS2 195-282, with amino acid coordinates aligned to TRA2A-201.",
  "TRA2A-202/212 lack the N-terminal 101 aa present in TRA2A-201, corresponding to the RS1 segment shown in the figure.",
  "TRA2A-212 maps exactly to TRA2A-201 amino acids 102-282.",
  "TRA2A-202 is nearly identical to TRA2A-212 but has one fewer arginine in the C-terminal R-rich region."
), file.path(out_dir, "README_TRA2A_protein_isoform_comparison.txt"))

message("Saved: ", png_path)
print(summary_dt)
