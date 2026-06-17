# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
gene_query <- if (length(args)) args[1] else "MAPKBP1"

rdata <- input_path("SUDHL.RData")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

e <- new.env()
load(rdata, envir = e)

annotation <- as.data.table(e$Annotation)[Condition == "SUDHL8"]
samples <- annotation$Sample

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
tx_map <- tx_map[row_id %in% rownames(e$count_matrix)]

tx_gene <- tx_map[toupper(gene_name) == toupper(gene_query)]
if (!nrow(tx_gene)) {
  stop("No long-read count-matrix transcripts found for gene query: ", gene_query)
}

counts <- e$count_matrix[tx_gene$row_id, samples, drop = FALSE]
counts[is.na(counts)] <- 0

count_dt <- as.data.table(as.table(counts))
setnames(count_dt, c("row_id", "Sample", "count"))
count_dt[, count := as.numeric(count)]
count_dt <- merge(count_dt, tx_gene, by = "row_id", all.x = TRUE)
count_dt <- merge(count_dt, annotation, by = "Sample", all.x = TRUE)

sample_tot <- count_dt[, .(gene_total = sum(count)), by = .(Sample, Type, Replicate)]
count_dt <- merge(count_dt, sample_tot, by = c("Sample", "Type", "Replicate"))
count_dt[, proportion := fifelse(gene_total > 0, count / gene_total, NA_real_)]

cat("\nGene query:", gene_query, "\n")
cat("\nTranscripts found in long-read matrix\n")
print(tx_gene[, .(row_id, gene_id, gene_name, gene_biotype, transcript_name, transcript_biotype)])

cat("\nPer-sample total long-read counts\n")
print(sample_tot[order(Type, Replicate)])

cat("\nTotal count summary by group\n")
print(sample_tot[, .(
  n = .N,
  counts = paste(gene_total, collapse = ","),
  mean = mean(gene_total),
  median = median(gene_total),
  sd = sd(gene_total),
  samples_detected = sum(gene_total > 0)
), by = Type])

cat("\nDetected isoform counts per sample\n")
print(count_dt[count > 0, .(
  Sample, Type, Replicate, transcript_name, row_id,
  transcript_biotype, count, gene_total, proportion
)][order(Type, Replicate, -count)])

cat("\nAggregated isoform counts by group\n")
agg <- dcast(count_dt, Type ~ transcript_name, value.var = "count", fun.aggregate = sum)
print(agg)

cat("\nTotal count tests\n")
if (uniqueN(sample_tot$Type) == 2) {
  cat("Wilcoxon p =", tryCatch(wilcox.test(gene_total ~ Type, data = sample_tot, exact = FALSE)$p.value, error = function(e) NA_real_), "\n")
  cat("t-test p =", tryCatch(t.test(gene_total ~ Type, data = sample_tot)$p.value, error = function(e) NA_real_), "\n")
}

mat <- xtabs(count ~ Type + transcript_name, data = count_dt)
mat <- mat[, colSums(mat) > 0, drop = FALSE]
cat("\nIsoform usage contingency matrix\n")
print(mat)

if (nrow(mat) >= 2 && ncol(mat) >= 2) {
  cat("\nAggregated isoform usage tests\n")
  cat("Chi-square p =", suppressWarnings(chisq.test(mat)$p.value), "\n")
  set.seed(1)
  cat("Monte Carlo chi-square p =", chisq.test(mat, simulate.p.value = TRUE, B = 100000)$p.value, "\n")
}
