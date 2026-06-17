# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(GenomicRanges)
})

bam_dir <- input_path("cDNA", "t2g_v3_lenFiltered")
gtf_file <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
out_dir <- analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05", "Multi_metric_integration", "Full_cycle_IGV_style")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes <- c("MAPKBP1", "SEC24C", "TRA2A")

extract_attr <- function(x, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  hit <- regexec(pattern, x)
  out <- regmatches(x, hit)
  vapply(out, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
}

gtf <- fread(
  gtf_file,
  sep = "\t",
  header = FALSE,
  skip = "#",
  quote = "",
  col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr")
)
gtf[, gene_name := extract_attr(attr, "gene_name")]
gtf[, gene_id := extract_attr(attr, "gene_id")]

gene_coords <- gtf[feature == "gene" & gene_name %in% genes,
  .(gene_name, gene_id, seqname, start, end, strand, width = end - start + 1)
]
fwrite(gene_coords, file.path(out_dir, "full_cycle_genes_genome_coordinates.csv"))

gene_gr <- GRanges(
  seqnames = gene_coords$seqname,
  ranges = IRanges(gene_coords$start, gene_coords$end),
  strand = "*",
  gene_name = gene_coords$gene_name
)
gene_gr_by_gene <- split(gene_gr, gene_gr$gene_name)

bams <- sort(list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE))

sample_meta <- data.table(bam = bams, bam_name = basename(bams))
sample_meta[, fraction := fifelse(grepl("-SSU_", bam_name), "SSU",
                           fifelse(grepl("-RS_", bam_name), "RS",
                           fifelse(grepl("-DS_", bam_name), "DS", NA_character_)))]
sample_meta[, cell_line := fifelse(grepl("^SU8R-", bam_name), "Resistant", "Sensitive")]
sample_meta[, treatment := fifelse(grepl("-Vin-", bam_name), "VCR", "DMSO")]
sample_meta[, replicate := fifelse(grepl("Rep1", bam_name), "Rep1", "Rep2")]
sample_meta[, condition := paste(cell_line, treatment, sep = "_")]

count_gene <- function(bam, gr) {
  param <- ScanBamParam(which = gr, what = "qname", flag = scanBamFlag(isUnmappedQuery = FALSE))
  aln <- scanBam(bam, param = param)
  sum(vapply(aln, function(x) length(x$qname), integer(1)))
}

total_mapped <- function(bam) {
  stats <- idxstatsBam(bam)
  sum(stats$mapped)
}

totals <- sample_meta[, .(mapped_reads = total_mapped(bam)), by = .(bam, bam_name)]
sample_meta <- merge(sample_meta, totals, by = c("bam", "bam_name"))

counts <- rbindlist(lapply(seq_len(nrow(sample_meta)), function(i) {
  sm <- sample_meta[i]
  rbindlist(lapply(names(gene_gr_by_gene), function(g) {
    n <- count_gene(sm$bam, gene_gr_by_gene[[g]])
    data.table(
      bam_name = sm$bam_name,
      gene_name = g,
      raw_reads = n,
      mapped_reads = sm$mapped_reads,
      rpm = n / sm$mapped_reads * 1e6
    )
  }))
}))

counts <- merge(counts, sample_meta[, .(bam_name, fraction, cell_line, treatment, replicate, condition)], by = "bam_name")
fwrite(counts, file.path(out_dir, "full_cycle_genes_genome_bam_gene_window_counts_by_sample.csv"))

summary_by_condition <- counts[, .(
  mean_raw_reads = mean(raw_reads),
  mean_rpm = mean(rpm),
  rep1_raw_reads = raw_reads[replicate == "Rep1"][1],
  rep2_raw_reads = raw_reads[replicate == "Rep2"][1],
  rep1_rpm = rpm[replicate == "Rep1"][1],
  rep2_rpm = rpm[replicate == "Rep2"][1]
), by = .(gene_name, fraction, condition)]
fwrite(summary_by_condition, file.path(out_dir, "full_cycle_genes_genome_bam_gene_window_counts_by_condition.csv"))

delta <- dcast(
  summary_by_condition,
  gene_name + fraction ~ condition,
  value.var = "mean_rpm"
)
delta[, resistant_vcr_delta := Resistant_VCR - Resistant_DMSO]
delta[, sensitive_vcr_delta := Sensitive_VCR - Sensitive_DMSO]
delta[, interaction_delta := resistant_vcr_delta - sensitive_vcr_delta]
fwrite(delta, file.path(out_dir, "full_cycle_genes_genome_bam_gene_window_rpm_interaction_delta.csv"))

cat("\nGene coordinates:\n")
print(gene_coords)
cat("\nMean exonic BAM reads/RPM by condition:\n")
print(summary_by_condition[order(gene_name, fraction, condition)])
cat("\nRPM interaction deltas: (Resistant VCR-DMSO) - (Sensitive VCR-DMSO)\n")
print(delta[order(gene_name, fraction)])
