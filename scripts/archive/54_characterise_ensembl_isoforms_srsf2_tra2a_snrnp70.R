# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
})

gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
cds_fasta <- input_path("Homo_sapiens.GRCh38.cds.all.fa.gz")

target_tx <- data.table(
  gene_name = c("SRSF2", "SRSF2", "TRA2A", "TRA2A", "SNRNP70"),
  transcript_name = c("SRSF2-202", "SRSF2-203", "TRA2A-201", "TRA2A-202", "SNRNP70-209")
)

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

gtf_dt <- fread(
  gtf,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attr")
)
gtf_dt[, `:=`(start = as.integer(start), end = as.integer(end))]

gtf_dt[, `:=`(
  gene_id = grab_attr(attr, "gene_id"),
  gene_name = grab_attr(attr, "gene_name"),
  gene_biotype = grab_attr(attr, "gene_biotype"),
  transcript_id = grab_attr(attr, "transcript_id"),
  transcript_version = grab_attr(attr, "transcript_version"),
  transcript_name = grab_attr(attr, "transcript_name"),
  transcript_biotype = grab_attr(attr, "transcript_biotype"),
  protein_id = grab_attr(attr, "protein_id"),
  protein_version = grab_attr(attr, "protein_version"),
  ccds_id = grab_attr(attr, "ccds_id"),
  exon_number = grab_attr(attr, "exon_number")
)]
gtf_dt[, row_id := paste0(transcript_id, ".", transcript_version)]
gtf_dt[, protein_row_id := fifelse(!is.na(protein_id), paste0(protein_id, ".", protein_version), NA_character_)]

target_rows <- merge(target_tx, unique(gtf_dt[feature == "transcript", .(
  gene_name, transcript_name, row_id, transcript_id, transcript_version,
  gene_biotype, transcript_biotype, chr, start, end, strand
)]), by = c("gene_name", "transcript_name"), all.x = TRUE)

feature_summary <- gtf_dt[row_id %in% target_rows$row_id, .(
  n_exons = sum(feature == "exon"),
  exon_bp = sum(ifelse(feature == "exon", end - start + 1L, 0L)),
  n_cds_exons = sum(feature == "CDS"),
  cds_bp_from_gtf = sum(ifelse(feature == "CDS", end - start + 1L, 0L)),
  n_start_codons = sum(feature == "start_codon"),
  n_stop_codons = sum(feature == "stop_codon"),
  protein_id = na.omit(protein_row_id)[1],
  ccds_id = na.omit(ccds_id)[1]
), by = .(row_id, transcript_name)]

cds <- readDNAStringSet(cds_fasta)
cds_names <- names(cds)
names(cds) <- sub(" .*", "", cds_names)
cds_sub <- cds[intersect(target_rows$row_id, names(cds))]
aa <- translate(cds_sub, if.fuzzy.codon = "X")
aa_seq <- as.character(aa)
aa_seq <- sub("\\*$", "", aa_seq)

seq_summary <- data.table(
  row_id = names(cds_sub),
  cds_nt_from_fasta = width(cds_sub),
  aa_length = nchar(aa_seq),
  protein_sequence = aa_seq
)

summary_dt <- merge(target_rows, feature_summary, by = c("row_id", "transcript_name"), all.x = TRUE)
summary_dt <- merge(summary_dt, seq_summary[, .(row_id, cds_nt_from_fasta, aa_length)], by = "row_id", all.x = TRUE)
setcolorder(summary_dt, c(
  "gene_name", "transcript_name", "row_id", "transcript_biotype", "protein_id",
  "aa_length", "cds_nt_from_fasta", "cds_bp_from_gtf", "n_exons", "exon_bp",
  "n_cds_exons", "n_start_codons", "n_stop_codons", "ccds_id"
))

cat("\n=== Ensembl v114 transcript/protein summary ===\n")
print(summary_dt[, .(
  gene_name, transcript_name, row_id, transcript_biotype, protein_id,
  aa_length, cds_nt_from_fasta, cds_bp_from_gtf, n_exons, exon_bp,
  n_cds_exons, ccds_id
)])

compare_pair <- function(a_name, b_name) {
  a_id <- summary_dt[transcript_name == a_name, row_id][1]
  b_id <- summary_dt[transcript_name == b_name, row_id][1]
  a_seq <- aa_seq[[a_id]]
  b_seq <- aa_seq[[b_id]]
  cat("\n=== Protein comparison:", a_name, "vs", b_name, "===\n")
  if (is.null(a_seq) || is.null(b_seq) || is.na(a_seq) || is.na(b_seq)) {
    cat("At least one transcript has no CDS/protein sequence in Ensembl CDS FASTA.\n")
    return(invisible(NULL))
  }
  a_chars <- strsplit(a_seq, "")[[1]]
  b_chars <- strsplit(b_seq, "")[[1]]
  min_len <- min(length(a_chars), length(b_chars))
  same_index <- sum(a_chars[seq_len(min_len)] == b_chars[seq_len(min_len)])
  diffs <- data.table(position = seq_len(min_len), a = a_chars[seq_len(min_len)], b = b_chars[seq_len(min_len)])
  diffs <- diffs[a != b]
  cat(a_name, "length:", length(a_chars), "aa\n")
  cat(b_name, "length:", length(b_chars), "aa\n")
  cat("Same-index identity over shared length:", round(100 * same_index / min_len, 2), "%\n")
  cat("First same-index difference:", if (nrow(diffs)) diffs$position[1] else "none", "\n")
  cat("Number same-index differences:", nrow(diffs), "\n")
  cat("Length difference:", length(b_chars) - length(a_chars), "aa\n")
  if (length(a_chars) != length(b_chars)) {
    cat("Longer unique tail/head cannot be interpreted as alignment-aware domain loss without protein alignment.\n")
  }
  cat("\nFirst 120 aa:\n")
  cat(a_name, ":", substr(a_seq, 1, 120), "\n")
  cat(b_name, ":", substr(b_seq, 1, 120), "\n")
  cat("\nLast 80 aa:\n")
  cat(a_name, ":", substr(a_seq, max(1, nchar(a_seq) - 79), nchar(a_seq)), "\n")
  cat(b_name, ":", substr(b_seq, max(1, nchar(b_seq) - 79), nchar(b_seq)), "\n")
}

compare_pair("SRSF2-202", "SRSF2-203")
compare_pair("TRA2A-201", "TRA2A-202")

cat("\n=== SNRNP70-209 annotation check ===\n")
print(summary_dt[transcript_name == "SNRNP70-209", .(
  gene_name, transcript_name, row_id, transcript_biotype, protein_id,
  aa_length, cds_nt_from_fasta, n_cds_exons, n_exons
)])

cat("\n=== Protein sequences for selected coding transcripts ===\n")
for (tx in c("SRSF2-202", "SRSF2-203", "TRA2A-201", "TRA2A-202")) {
  id <- summary_dt[transcript_name == tx, row_id][1]
  cat(">", tx, "|", id, "|", summary_dt[transcript_name == tx, protein_id][1], "\n", sep = "")
  cat(aa_seq[[id]], "\n")
}
