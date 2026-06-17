# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)
library(GenomicRanges)
library(ggplot2)
library(scales)

rdata <- input_path("SUDHL.RData")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
out_dir <- analysis_path("LongRead_TRA2A", "Step2_TRA2A_targets_isoform_usage")
peak_dir <- file.path(out_dir, "ENCODE_TRA2A_eCLIP")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

standard_chr <- function(x) {
  x <- sub("^chr", "", x)
  fifelse(x == "M", "MT", x)
}

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

message("Loading long-read matrices...")
e <- new.env()
load(rdata, envir = e)
annotation <- as.data.table(e$Annotation)[Condition == "SUDHL8"]
samples <- annotation$Sample

message("Parsing GTF genes and transcripts...")
gtf_dt <- fread(
  gtf,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
)

gene_dt <- gtf_dt[feature == "gene"]
gene_attrs <- gene_dt$attributes
genes <- data.table(
  chr = standard_chr(gene_dt$chr),
  start = as.integer(gene_dt$start),
  end = as.integer(gene_dt$end),
  strand = gene_dt$strand,
  gene_id = grab_attr(gene_attrs, "gene_id"),
  gene_name = grab_attr(gene_attrs, "gene_name"),
  gene_biotype = grab_attr(gene_attrs, "gene_biotype")
)
genes[, gene_id_versionless := sub("[.][0-9]+$", "", gene_id)]
gene_gr <- GRanges(
  seqnames = genes$chr,
  ranges = IRanges(genes$start, genes$end),
  strand = genes$strand,
  gene_id = genes$gene_id_versionless,
  gene_name = genes$gene_name,
  gene_biotype = genes$gene_biotype
)

tx_dt <- gtf_dt[feature == "transcript"]
tx_attrs <- tx_dt$attributes
tx_map <- data.table(
  chr = standard_chr(tx_dt$chr),
  start = as.integer(tx_dt$start),
  end = as.integer(tx_dt$end),
  strand = tx_dt$strand,
  gene_id = grab_attr(tx_attrs, "gene_id"),
  gene_name = grab_attr(tx_attrs, "gene_name"),
  gene_biotype = grab_attr(tx_attrs, "gene_biotype"),
  transcript_id = grab_attr(tx_attrs, "transcript_id"),
  transcript_version = grab_attr(tx_attrs, "transcript_version"),
  transcript_name = grab_attr(tx_attrs, "transcript_name"),
  transcript_biotype = grab_attr(tx_attrs, "transcript_biotype")
)
tx_map[, gene_id := sub("[.][0-9]+$", "", gene_id)]
tx_map[, row_id := paste0(transcript_id, ".", transcript_version)]
tx_map <- tx_map[row_id %in% rownames(e$count_matrix)]

message("Reading ENCODE TRA2A eCLIP peaks...")
peak_files <- list.files(peak_dir, pattern = "[.]bed[.]gz$", full.names = TRUE)
stopifnot(length(peak_files) > 0)

read_peak_file <- function(f) {
  dt <- fread(f, header = FALSE)
  if (ncol(dt) < 3) stop("Peak file has fewer than 3 BED columns: ", f)
  setnames(dt, paste0("V", seq_len(ncol(dt))))
  dt <- dt[, .(
    chr = standard_chr(V1),
    start = as.integer(V2) + 1L,
    end = as.integer(V3),
    peak_id = if ("V4" %in% names(dt)) V4 else paste0(basename(f), "_", seq_len(.N)),
    score = if ("V5" %in% names(dt)) suppressWarnings(as.numeric(V5)) else NA_real_,
    strand = if ("V6" %in% names(dt)) V6 else "*"
  )]
  dt[, source_file := basename(f)]
  dt[, source_cell := fifelse(grepl("K562", source_file), "K562",
                              fifelse(grepl("HepG2", source_file), "HepG2", "unknown"))]
  dt
}

peaks <- rbindlist(lapply(peak_files, read_peak_file), fill = TRUE)
peak_gr <- GRanges(seqnames = peaks$chr, ranges = IRanges(peaks$start, peaks$end), source_cell = peaks$source_cell)
hits <- findOverlaps(peak_gr, gene_gr, ignore.strand = TRUE)
target_map <- unique(data.table(
  source_cell = mcols(peak_gr)$source_cell[queryHits(hits)],
  gene_id = mcols(gene_gr)$gene_id[subjectHits(hits)],
  gene_name = mcols(gene_gr)$gene_name[subjectHits(hits)],
  gene_biotype = mcols(gene_gr)$gene_biotype[subjectHits(hits)]
))
target_summary <- target_map[
  ,
  .(
    n_source_cell_lines = uniqueN(source_cell),
    source_cell_lines = paste(sort(unique(source_cell)), collapse = ";")
  ),
  by = .(gene_id, gene_name, gene_biotype)
]
target_summary[, TRA2A_eCLIP_target := TRUE]

message("Building SUDHL8 isoform proportions...")
counts <- e$count_matrix[tx_map$row_id, samples, drop = FALSE]
counts[is.na(counts)] <- 0
count_dt <- as.data.table(as.table(counts))
setnames(count_dt, c("row_id", "Sample", "count"))
count_dt[, count := as.numeric(count)]
count_dt <- merge(count_dt, tx_map, by = "row_id", all.x = TRUE)
count_dt <- merge(count_dt, annotation, by = "Sample", all.x = TRUE)
count_dt <- count_dt[!is.na(gene_id)]

gene_sample_totals <- count_dt[, .(gene_total_count = sum(count)), by = .(gene_id, gene_name, Sample, Type, Replicate)]
count_dt <- merge(count_dt, gene_sample_totals, by = c("gene_id", "gene_name", "Sample", "Type", "Replicate"))
count_dt[, proportion := fifelse(gene_total_count > 0, count / gene_total_count, NA_real_)]

expressed_gene_summary <- gene_sample_totals[
  ,
  .(
    mean_total_count = mean(gene_total_count),
    total_count = sum(gene_total_count),
    n_samples_detected = sum(gene_total_count > 0)
  ),
  by = .(gene_id, gene_name)
]

gene_tx_n <- unique(count_dt[count > 0, .(gene_id, row_id)])[, .(n_detected_transcripts = .N), by = gene_id]
expressed_gene_summary <- merge(expressed_gene_summary, gene_tx_n, by = "gene_id", all.x = TRUE)
expressed_gene_summary[is.na(n_detected_transcripts), n_detected_transcripts := 0L]
expressed_gene_summary <- merge(expressed_gene_summary, target_summary, by = c("gene_id", "gene_name"), all.x = TRUE)
expressed_gene_summary[is.na(TRA2A_eCLIP_target), TRA2A_eCLIP_target := FALSE]

test_gene <- function(gid) {
  dt <- count_dt[gene_id == gid]
  if (uniqueN(dt$row_id[dt$count > 0]) < 2) return(NULL)
  if (sum(dt$count) < 20) return(NULL)

  agg <- dcast(dt, Type ~ row_id, value.var = "count", fun.aggregate = sum)
  mat <- as.matrix(agg[, -"Type"])
  rownames(mat) <- agg$Type
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (ncol(mat) < 2 || nrow(mat) < 2) return(NULL)

  p_value <- tryCatch({
    suppressWarnings(chisq.test(mat, simulate.p.value = any(chisq.test(mat)$expected < 5), B = 10000)$p.value)
  }, error = function(e) NA_real_)

  prop_summary <- dt[
    ,
    .(mean_prop = mean(proportion, na.rm = TRUE), mean_count = mean(count)),
    by = .(Type, row_id, transcript_name, transcript_biotype)
  ]
  wide <- dcast(prop_summary, row_id + transcript_name + transcript_biotype ~ Type, value.var = "mean_prop")
  if (!"Original" %in% names(wide)) wide[, Original := 0]
  if (!"Resistant" %in% names(wide)) wide[, Resistant := 0]
  wide[is.na(Original), Original := 0]
  wide[is.na(Resistant), Resistant := 0]
  wide[, resistant_minus_original := Resistant - Original]

  top_original <- wide[which.max(Original), transcript_name]
  top_resistant <- wide[which.max(Resistant), transcript_name]
  top_shift <- wide[which.max(abs(resistant_minus_original))]
  totals <- gene_sample_totals[gene_id == gid, .(mean_total = mean(gene_total_count)), by = Type]
  original_total <- totals[Type == "Original", mean_total]
  resistant_total <- totals[Type == "Resistant", mean_total]

  data.table(
    gene_id = gid,
    gene_name = dt$gene_name[1],
    n_detected_transcripts = uniqueN(dt$row_id[dt$count > 0]),
    total_count = sum(dt$count),
    original_mean_total_count = ifelse(length(original_total), original_total, NA_real_),
    resistant_mean_total_count = ifelse(length(resistant_total), resistant_total, NA_real_),
    top_original_isoform = top_original,
    top_resistant_isoform = top_resistant,
    dominant_isoform_switch = top_original != top_resistant,
    max_abs_proportion_shift = max(abs(wide$resistant_minus_original), na.rm = TRUE),
    max_shift_isoform = top_shift$transcript_name,
    max_shift_transcript_id = top_shift$row_id,
    max_shift_original_prop = top_shift$Original,
    max_shift_resistant_prop = top_shift$Resistant,
    max_shift_resistant_minus_original = top_shift$resistant_minus_original,
    chi_square_pvalue = p_value
  )
}

eligible_target_genes <- expressed_gene_summary[
  TRA2A_eCLIP_target == TRUE &
    n_detected_transcripts >= 2 &
    total_count >= 20,
  gene_id
]

target_tests <- rbindlist(lapply(eligible_target_genes, test_gene), fill = TRUE)
if (nrow(target_tests) > 0) {
  target_tests[, padj := p.adjust(chi_square_pvalue, method = "BH")]
  target_tests <- merge(
    target_tests,
    target_summary[, .(gene_id, source_cell_lines, n_source_cell_lines)],
    by = "gene_id",
    all.x = TRUE
  )
  setorder(target_tests, chi_square_pvalue, -max_abs_proportion_shift)
}

target_transcript_shifts <- count_dt[gene_id %in% eligible_target_genes, .(
  mean_prop = mean(proportion, na.rm = TRUE),
  mean_count = mean(count)
), by = .(gene_id, gene_name, row_id, transcript_name, transcript_biotype, Type)]
target_transcript_shifts <- dcast(
  target_transcript_shifts,
  gene_id + gene_name + row_id + transcript_name + transcript_biotype ~ Type,
  value.var = c("mean_prop", "mean_count")
)
for (nm in c("mean_prop_Original", "mean_prop_Resistant", "mean_count_Original", "mean_count_Resistant")) {
  if (!nm %in% names(target_transcript_shifts)) target_transcript_shifts[, (nm) := 0]
  target_transcript_shifts[is.na(get(nm)), (nm) := 0]
}
target_transcript_shifts[, resistant_minus_original_prop := mean_prop_Resistant - mean_prop_Original]
target_transcript_shifts[, abs_resistant_minus_original_prop := abs(resistant_minus_original_prop)]
setorder(target_transcript_shifts, -abs_resistant_minus_original_prop)

write.csv(peaks, file.path(out_dir, "ENCODE_TRA2A_eCLIP_GRCh38_peaks_used.csv"), row.names = FALSE)
write.csv(target_summary[order(gene_name)], file.path(out_dir, "TRA2A_eCLIP_target_genes_from_peak_gene_overlap.csv"), row.names = FALSE)
write.csv(expressed_gene_summary[order(-TRA2A_eCLIP_target, -total_count)], file.path(out_dir, "SUDHL8_longread_gene_expression_target_annotation.csv"), row.names = FALSE)
write.csv(target_tests, file.path(out_dir, "TRA2A_target_genes_SUDHL8_isoform_usage_tests.csv"), row.names = FALSE)
write.csv(target_transcript_shifts, file.path(out_dir, "TRA2A_target_transcript_level_isoform_shifts.csv"), row.names = FALSE)

top_plot_genes <- head(target_tests[!is.na(chi_square_pvalue)][order(chi_square_pvalue, -max_abs_proportion_shift), gene_id], 12)
plot_dt <- count_dt[gene_id %in% top_plot_genes]
plot_dt[, Type := factor(Type, levels = c("Original", "Resistant"))]
plot_sum <- plot_dt[, .(prop = sum(count) / unique(gene_total_count)), by = .(Sample, Type, gene_name, transcript_name)]
plot_sum[, gene_name := factor(gene_name, levels = target_tests[gene_id %in% top_plot_genes][order(chi_square_pvalue, -max_abs_proportion_shift), gene_name])]

if (nrow(plot_sum) > 0) {
  p <- ggplot(plot_sum, aes(x = Type, y = prop, fill = transcript_name)) +
    geom_col(width = 0.78, color = "white", linewidth = 0.1, position = "fill") +
    facet_wrap(~ gene_name, scales = "free_y") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Top TRA2A eCLIP target genes with SUDHL8 isoform usage shifts",
      subtitle = "Targets defined by ENCODE TRA2A eCLIP peak overlap with GRCh38 genes",
      x = NULL,
      y = "Isoform proportion",
      fill = "Transcript"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )
  ggsave(file.path(out_dir, "top_TRA2A_target_isoform_usage_shifts_SUDHL8.png"), p, width = 11, height = 8, dpi = 300)
}

source_note <- c(
  "Step 2 target-gene definition:",
  "TRA2A target genes were defined by overlapping released ENCODE TRA2A eCLIP GRCh38 merged peak BED files with GRCh38.114 gene coordinates.",
  "ENCODE experiments used:",
  "ENCSR314UMJ HepG2 TRA2A eCLIP, merged GRCh38 peak file ENCFF766OCH.",
  "ENCSR365NVO K562 TRA2A eCLIP, merged GRCh38 peak file ENCFF726PFJ.",
  "This is a focused exploratory isoform-usage screen because DRIMSeq/DEXSeq/IsoformSwitchAnalyzeR were not installed in the current R environment.",
  "Per-gene p-values are chi-square tests on aggregated transcript counts in SUDHL8 Original vs Resistant; effect sizes are mean within-gene isoform proportion shifts across replicates."
)
writeLines(source_note, file.path(out_dir, "README_step2_source_and_method.txt"))

cat("\nENCODE TRA2A eCLIP target genes from peak-gene overlap:", nrow(target_summary), "\n")
cat("Targets expressed in SUDHL8 long-read data:", expressed_gene_summary[TRA2A_eCLIP_target == TRUE & total_count > 0, .N], "\n")
cat("Targets eligible for isoform-usage testing:", length(eligible_target_genes), "\n")
cat("\nTop TRA2A target isoform-usage test results:\n")
print(head(target_tests, 20))
cat("\nTop transcript-level shifts among TRA2A targets:\n")
print(head(target_transcript_shifts, 20))
cat("\nOutput directory:\n")
cat(out_dir, "\n")
