# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(gprofiler2)
})

rdata <- input_path("SUDHL.RData")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
base_dir <- analysis_path()
count_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
metric_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05", "Results")

p_cut <- 0.05
lfc_cut <- 0.7

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

direction_call <- function(logfc, pval) {
  fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc >= lfc_cut, "Up",
          fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc <= -lfc_cut, "Down", "NS"))
}

load_longread <- function() {
  e <- new.env()
  load(rdata, envir = e)
  e
}

build_tx_map <- function(count_matrix) {
  gtf_dt <- fread(
    gtf,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "#",
    col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  )
  tx_dt <- gtf_dt[feature == "transcript"]
  attrs <- tx_dt$attributes
  tx_map <- data.table(
    gene_id = grab_attr(attrs, "gene_id"),
    gene_name = grab_attr(attrs, "gene_name"),
    gene_biotype = grab_attr(attrs, "gene_biotype"),
    transcript_id = grab_attr(attrs, "transcript_id"),
    transcript_version = grab_attr(attrs, "transcript_version"),
    transcript_name = grab_attr(attrs, "transcript_name"),
    transcript_biotype = grab_attr(attrs, "transcript_biotype")
  )
  tx_map[, gene_id := sub("[.][0-9]+$", "", gene_id)]
  tx_map[, row_id := paste0(transcript_id, ".", transcript_version)]
  tx_map[row_id %in% rownames(count_matrix)]
}

get_baseline_splicing_hits <- function() {
  metric <- fread(file.path(metric_dir, "ribosome_efficiency_score", "Resistance_baseline_limma_all_genes.csv"))
  count <- fread(file.path(count_dir, "Fraction_RS", "Resistance_baseline_RS_psite_limma_all_genes.csv"))
  metric[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  count[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  if (!"gene_name" %in% names(metric)) metric[, gene_name := gene_id_clean]
  if (!"gene_name" %in% names(count)) count[, gene_name := gene_id_clean]
  metric[, metric_direction := direction_call(logFC, P.Value)]
  count[, count_direction := direction_call(logFC, P.Value)]
  joined <- merge(
    metric[, .(gene_id_clean, gene_name_metric = gene_name, metric_logFC = logFC, metric_P = P.Value, metric_direction)],
    count[, .(gene_id_clean, gene_name_count = gene_name, RS_logFC = logFC, RS_P = P.Value, count_direction)],
    by = "gene_id_clean"
  )
  joined[, gene_name := fifelse(!is.na(gene_name_metric) & nzchar(gene_name_metric), gene_name_metric, gene_name_count)]
  conv_down <- joined[metric_direction == "Down" & count_direction == "Down"]

  res <- gost(
    query = unique(conv_down$gene_id_clean),
    organism = "hsapiens",
    sources = c("GO:BP", "GO:CC"),
    correction_method = "g_SCS",
    domain_scope = "annotated",
    user_threshold = 0.05,
    evcodes = TRUE
  )
  dt <- as.data.table(res$result)
  splicing_term <- dt[
    significant == TRUE &
      term_size <= 500 &
      intersection_size >= 5 &
      source == "GO:BP" &
      term_name == "regulation of RNA splicing"
  ]
  ids <- unique(unlist(strsplit(splicing_term$intersection[1], ",")))
  ids <- sub("\\.\\d+$", "", ids)
  conv_down[gene_id_clean %in% ids, .(
    gene_id_clean,
    gene_name,
    metric_logFC = round(metric_logFC, 3),
    metric_P = signif(metric_P, 3),
    RS_logFC = round(RS_logFC, 3),
    RS_P = signif(RS_P, 3)
  )][order(gene_name)]
}

e <- load_longread()
annotation <- as.data.table(e$Annotation)[Condition == "SUDHL8"]
samples <- annotation$Sample
tx_map <- build_tx_map(e$count_matrix)

splicing_hits <- get_baseline_splicing_hits()
cat("\nConvergent baseline ribosome engagement/RS Down genes driving GO: regulation of RNA splicing:\n")
print(splicing_hits)

requested_splicing <- unique(c("LSM4", "SON", splicing_hits$gene_name))
requested_export <- c("CHTOP", "ZC3H11A", "ALKBH5", "AGFG1", "SENP2")
all_genes <- unique(c(requested_splicing, requested_export))

cat("\nGene groups requested:\n")
cat("Splicing/conclusion 1 genes:", paste(requested_splicing, collapse = ", "), "\n")
cat("mRNA export/conclusion 3 genes:", paste(requested_export, collapse = ", "), "\n")

gene_tx <- tx_map[toupper(gene_name) %in% toupper(all_genes)]
missing_genes <- setdiff(toupper(all_genes), unique(toupper(gene_tx$gene_name)))
if (length(missing_genes)) {
  cat("\nGenes not found in long-read transcript count matrix:", paste(missing_genes, collapse = ", "), "\n")
}

counts <- e$count_matrix[gene_tx$row_id, samples, drop = FALSE]
counts[is.na(counts)] <- 0
count_dt <- as.data.table(as.table(counts))
setnames(count_dt, c("row_id", "Sample", "count"))
count_dt[, count := as.numeric(count)]
count_dt <- merge(count_dt, gene_tx[, .(row_id, gene_id, gene_name, transcript_name, transcript_biotype)], by = "row_id", all.x = TRUE)
count_dt <- merge(count_dt, annotation[, .(Sample, Condition, Type, Replicate)], by = "Sample", all.x = TRUE)

gene_sample <- count_dt[, .(
  nanopore_gene_total_count = sum(count),
  n_detected_transcripts = uniqueN(row_id[count > 0])
), by = .(gene_name, gene_id, Sample, Type, Replicate)]
setorder(gene_sample, gene_name, Type, Replicate)

gene_summary <- gene_sample[, .(
  counts_by_replicate = paste(nanopore_gene_total_count, collapse = ","),
  mean_count = round(mean(nanopore_gene_total_count), 2),
  median_count = median(nanopore_gene_total_count),
  min_count = min(nanopore_gene_total_count),
  max_count = max(nanopore_gene_total_count),
  samples_detected = sum(nanopore_gene_total_count > 0)
), by = .(gene_name, Type)]
setorder(gene_summary, gene_name, Type)

gene_wide <- dcast(
  gene_sample,
  gene_name + gene_id ~ Sample,
  value.var = "nanopore_gene_total_count",
  fill = 0
)
setorder(gene_wide, gene_name)

tx_summary <- count_dt[count > 0, .(
  total_count = sum(count),
  samples_detected = uniqueN(Sample)
), by = .(gene_name, transcript_name, transcript_biotype)]
setorder(tx_summary, gene_name, -total_count)

cat("\nPer-sample total Nanopore gene counts:\n")
print(gene_sample[, .(gene_name, Sample, Type, Replicate, nanopore_gene_total_count, n_detected_transcripts)])

cat("\nGroup summary of total Nanopore gene counts:\n")
print(gene_summary)

cat("\nWide per-sample count table:\n")
print(gene_wide)

cat("\nDetected transcript-level support summary, for context only:\n")
print(tx_summary)
