# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(zoo)
})

base_dir <- analysis_path()
index_dir <- file.path(base_dir, "Translation_indexes_fixed")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
out_dir <- file.path(base_dir, "LongRead_TRA2A", "TRA2A_ribosome_locus_coverage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene <- "TRA2A"
canonical_tx <- "ENST00000297071.9"
psite_file <- file.path(index_dir, "transcript_psite_matrix_long_ALL_samples.csv")
metric_file <- file.path(index_dir, "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

standard_chr <- function(x) {
  x <- sub("^chr", "", x)
  fifelse(x == "M", "MT", x)
}

message("Finding TRA2A transcripts in metric table...")
tx_map <- unique(fread(metric_file, select = c("transcript", "gene_id_clean", "gene_name"))[gene_name == gene])

message("Reading TRA2A P-site rows...")
psite <- fread(
  psite_file,
  select = c("sample", "fraction", "transcript", "codon_pos", "psite_count", "psite_cpm")
)[transcript %in% tx_map$transcript]
psite <- merge(psite, tx_map, by = "transcript")

tx_totals <- psite[, .(
  total_psites = sum(psite_count, na.rm = TRUE),
  total_cpm = sum(psite_cpm, na.rm = TRUE),
  min_pos = min(codon_pos),
  max_pos = max(codon_pos),
  n_positions = uniqueN(codon_pos)
), by = .(gene_name, transcript)][order(-total_psites)]

fwrite(tx_totals, file.path(out_dir, "TRA2A_transcripts_psite_coverage_summary.csv"))

plot_transcripts <- unique(c(canonical_tx, tx_totals$transcript[1]))
plot_transcripts <- plot_transcripts[plot_transcripts %in% tx_totals$transcript]

message("Parsing GTF transcript coordinates...")
gtf_dt <- fread(
  gtf,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
)

get_tx_features <- function(tx_versioned) {
  tx_id <- sub("[.][0-9]+$", "", tx_versioned)
  dt <- gtf_dt[grepl(paste0('transcript_id "', tx_id, '"'), attributes)]
  attrs <- dt$attributes
  dt[, transcript_id := grab_attr(attrs, "transcript_id")]
  dt[, transcript_version := grab_attr(attrs, "transcript_version")]
  dt[, row_id := paste0(transcript_id, ".", transcript_version)]
  dt <- dt[row_id == tx_versioned]
  if (!nrow(dt)) return(NULL)

  strand <- unique(dt[feature == "transcript", strand])
  exons <- dt[feature == "exon", .(
    chr = standard_chr(chr),
    start = as.integer(start),
    end = as.integer(end),
    strand = strand,
    exon_number = suppressWarnings(as.integer(grab_attr(attributes, "exon_number")))
  )]
  if (strand == "+") {
    setorder(exons, start, end)
  } else {
    setorder(exons, -start, -end)
  }
  exons[, exon_len := end - start + 1L]
  exons[, tx_start := cumsum(shift(exon_len, fill = 0L)) + 1L]
  exons[, tx_end := tx_start + exon_len - 1L]

  map_interval <- function(g_start, g_end) {
    hits <- exons[start <= g_end & end >= g_start]
    if (!nrow(hits)) return(NULL)
    rbindlist(lapply(seq_len(nrow(hits)), function(i) {
      h <- hits[i]
      ov_start <- max(g_start, h$start)
      ov_end <- min(g_end, h$end)
      if (strand == "+") {
        tx_s <- h$tx_start + (ov_start - h$start)
        tx_e <- h$tx_start + (ov_end - h$start)
      } else {
        tx_s <- h$tx_start + (h$end - ov_end)
        tx_e <- h$tx_start + (h$end - ov_start)
      }
      data.table(start = min(tx_s, tx_e), end = max(tx_s, tx_e))
    }))
  }

  cds <- dt[feature == "CDS"]
  cds_tx <- if (nrow(cds)) {
    rbindlist(lapply(seq_len(nrow(cds)), function(i) map_interval(as.integer(cds$start[i]), as.integer(cds$end[i]))))
  } else {
    NULL
  }

  transcript_len <- max(exons$tx_end)
  if (!is.null(cds_tx) && nrow(cds_tx)) {
    cds_start <- min(cds_tx$start)
    cds_end <- max(cds_tx$end)
    regions <- rbindlist(list(
      data.table(region = "5'UTR", start = 1L, end = cds_start - 1L),
      data.table(region = "CDS", start = cds_start, end = cds_end),
      data.table(region = "3'UTR", start = cds_end + 1L, end = transcript_len)
    ))
    regions <- regions[start <= end]
  } else {
    regions <- data.table(region = "Transcript", start = 1L, end = transcript_len)
  }
  regions[, transcript := tx_versioned]
  regions[]
}

regions <- rbindlist(lapply(plot_transcripts, get_tx_features), fill = TRUE)
fwrite(regions, file.path(out_dir, "TRA2A_plot_transcript_regions_transcript_coordinates.csv"))

psite[, condition := fifelse(grepl("SU8R-DMSO|SU8-R-DMSO", sample), "Resistant DMSO",
                      fifelse(grepl("SU8R-Vin|SU8-R-Vin", sample), "Resistant VCR",
                      fifelse(grepl("SU8-DMSO", sample), "Sensitive DMSO",
                      fifelse(grepl("SU8-Vin", sample), "Sensitive VCR", NA_character_))))]
psite <- psite[!is.na(condition)]
psite[, condition := factor(condition, levels = c("Sensitive DMSO", "Sensitive VCR", "Resistant DMSO", "Resistant VCR"))]
psite[, fraction := factor(fraction, levels = c("SSU", "RS", "DS"))]

mean_cov <- psite[transcript %in% plot_transcripts, .(
  mean_cpm = mean(psite_cpm, na.rm = TRUE),
  sem_cpm = sd(psite_cpm, na.rm = TRUE) / sqrt(.N),
  mean_count = mean(psite_count, na.rm = TRUE)
), by = .(gene_name, transcript, fraction, condition, codon_pos)]

all_positions <- mean_cov[, .(codon_pos = seq(min(codon_pos), max(codon_pos))), by = .(gene_name, transcript, fraction, condition)]
mean_cov <- merge(all_positions, mean_cov, by = c("gene_name", "transcript", "fraction", "condition", "codon_pos"), all.x = TRUE)
mean_cov[is.na(mean_cpm), mean_cpm := 0]
mean_cov[is.na(mean_count), mean_count := 0]
mean_cov[, smooth_cpm := zoo::rollmean(mean_cpm, k = 9, fill = NA, align = "center"), by = .(transcript, fraction, condition)]

track_summary <- mean_cov[, .(
  total_mean_cpm = sum(mean_cpm),
  max_mean_cpm = max(mean_cpm),
  covered_positions = sum(mean_cpm > 0),
  utr5_mean_cpm = {
    r <- regions[transcript == .BY$transcript & region == "5'UTR"]
    if (nrow(r)) mean(mean_cpm[codon_pos >= r$start & codon_pos <= r$end], na.rm = TRUE) else NA_real_
  },
  cds_mean_cpm = {
    r <- regions[transcript == .BY$transcript & region == "CDS"]
    if (nrow(r)) mean(mean_cpm[codon_pos >= r$start & codon_pos <= r$end], na.rm = TRUE) else NA_real_
  }
), by = .(gene_name, transcript, fraction, condition)]
fwrite(track_summary, file.path(out_dir, "TRA2A_locus_coverage_track_summary.csv"))
fwrite(mean_cov, file.path(out_dir, "TRA2A_locus_coverage_mean_by_position.csv"))

region_cols <- c("5'UTR" = "#EAF2F8", "CDS" = "#FDEDEC", "3'UTR" = "#EAF7EA", "Transcript" = "#F2F2F2")
fraction_cols <- c("SSU" = "#4C78A8", "RS" = "#54A24B", "DS" = "#E45756")

make_coverage_grid <- function(tx) {
  d <- mean_cov[transcript == tx]
  r <- regions[transcript == tx]
  ymax <- max(d$mean_cpm, na.rm = TRUE)
  bg <- r[, .(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, region)]

  p <- ggplot(d, aes(codon_pos, mean_cpm)) +
    geom_rect(
      data = bg,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
      inherit.aes = FALSE,
      alpha = 0.45
    ) +
    geom_col(aes(color = fraction), width = 1, fill = "#2F3A4A", alpha = 0.65, linewidth = 0) +
    geom_line(aes(y = smooth_cpm, color = fraction), linewidth = 0.45, na.rm = TRUE) +
    facet_grid(fraction ~ condition, scales = "free_y") +
    scale_fill_manual(values = region_cols, drop = FALSE) +
    scale_color_manual(values = fraction_cols, guide = "none") +
    labs(
      title = paste0("TRA2A ribosome/P-site coverage across conditions and fractions"),
      subtitle = paste0("Transcript: ", tx, " | replicate mean P-site CPM | shaded regions from GRCh38.114 GTF"),
      x = "Transcript coordinate",
      y = "Mean P-site CPM",
      fill = "Region"
    ) +
    coord_cartesian(ylim = c(0, ymax * 1.05)) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      legend.position = "bottom"
    )

  png_path <- file.path(out_dir, paste0("TRA2A_", tx, "_coverage_fraction_condition_grid.png"))
  pdf_path <- file.path(out_dir, paste0("TRA2A_", tx, "_coverage_fraction_condition_grid.pdf"))
  ggsave(png_path, p, width = 13, height = 8, dpi = 300, bg = "white")
  ggsave(pdf_path, p, width = 13, height = 8, bg = "white")

  for (fr in levels(d$fraction)) {
    df <- d[fraction == fr]
    pf <- ggplot(df, aes(codon_pos, mean_cpm)) +
      geom_rect(
        data = bg,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
        inherit.aes = FALSE,
        alpha = 0.45
      ) +
      geom_col(width = 1, fill = fraction_cols[[fr]], alpha = 0.65, linewidth = 0) +
      geom_line(aes(y = smooth_cpm), color = "#111827", linewidth = 0.45, na.rm = TRUE) +
      facet_grid(condition ~ ., scales = "free_y") +
      scale_fill_manual(values = region_cols, drop = FALSE) +
      labs(
        title = paste0("TRA2A ", fr, " P-site coverage"),
        subtitle = paste0("Transcript: ", tx, " | replicate mean P-site CPM"),
        x = "Transcript coordinate",
        y = "Mean P-site CPM",
        fill = "Region"
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey92", color = "grey70"),
        legend.position = "bottom"
      )
    ggsave(file.path(out_dir, paste0("TRA2A_", tx, "_", fr, "_coverage_by_condition.png")), pf, width = 9, height = 6.5, dpi = 300, bg = "white")
    ggsave(file.path(out_dir, paste0("TRA2A_", tx, "_", fr, "_coverage_by_condition.pdf")), pf, width = 9, height = 6.5, bg = "white")
  }

  r5 <- r[region == "5'UTR"]
  if (nrow(r5)) {
    d5 <- d[fraction == "SSU" & codon_pos >= r5$start & codon_pos <= r5$end]
    if (nrow(d5) && sum(d5$mean_cpm > 0, na.rm = TRUE) > 0) {
      p5 <- ggplot(d5, aes(codon_pos, mean_cpm, color = condition)) +
        geom_col(aes(fill = condition), width = 1, alpha = 0.45, position = "identity", color = NA) +
        geom_line(aes(y = smooth_cpm), linewidth = 0.65, na.rm = TRUE) +
        scale_color_manual(values = c("Sensitive DMSO" = "#4C78A8", "Sensitive VCR" = "#72B7B2", "Resistant DMSO" = "#F58518", "Resistant VCR" = "#E45756")) +
        scale_fill_manual(values = c("Sensitive DMSO" = "#4C78A8", "Sensitive VCR" = "#72B7B2", "Resistant DMSO" = "#F58518", "Resistant VCR" = "#E45756")) +
        labs(
          title = "TRA2A 5'UTR SSU signal",
          subtitle = paste0("Transcript: ", tx, " | 5'UTR coordinates ", r5$start, "-", r5$end),
          x = "Transcript coordinate within 5'UTR",
          y = "Mean SSU P-site CPM",
          color = "Condition",
          fill = "Condition"
        ) +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
      ggsave(file.path(out_dir, paste0("TRA2A_", tx, "_SSU_5UTR_focused.png")), p5, width = 9, height = 4.8, dpi = 300, bg = "white")
      ggsave(file.path(out_dir, paste0("TRA2A_", tx, "_SSU_5UTR_focused.pdf")), p5, width = 9, height = 4.8, bg = "white")
    }
  }
}

lapply(plot_transcripts, make_coverage_grid)

readme <- c(
  "TRA2A ribosome/P-site locus coverage plots.",
  paste0("Input P-site matrix: ", psite_file),
  paste0("Annotation GTF: ", gtf),
  paste0("Transcripts plotted: ", paste(plot_transcripts, collapse = ", ")),
  "Coverage is replicate mean P-site CPM by transcript coordinate.",
  "Conditions: Sensitive DMSO, Sensitive VCR, Resistant DMSO, Resistant VCR.",
  "Fractions: SSU, RS, DS.",
  "Shaded regions mark 5'UTR, CDS, and 3'UTR inferred from GRCh38.114 transcript annotation.",
  "A 9-position rolling mean line is overlaid on raw per-position CPM bars.",
  "Note: this is transcript-coordinate coverage from the transcriptome-derived P-site matrix, not genome-browser BAM coverage."
)
writeLines(readme, file.path(out_dir, "README_TRA2A_locus_coverage_plots.txt"))

message("Saved TRA2A coverage plots to: ", out_dir)
print(tx_totals)
print(track_summary[order(transcript, fraction, condition)])
