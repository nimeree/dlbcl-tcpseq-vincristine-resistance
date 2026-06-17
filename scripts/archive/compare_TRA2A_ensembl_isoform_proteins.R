# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)
library(Biostrings)

gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
cds_fasta <- input_path("Homo_sapiens.GRCh38.cds.all.fa.gz")

txs <- c(
  "TRA2A-201" = "ENST00000297071.9",
  "TRA2A-202" = "ENST00000392502.8",
  "TRA2A-212" = "ENST00000621813.4"
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

tra2a <- gtf_dt[grepl('gene_name "TRA2A"', attr)]
tra2a[, transcript_id := grab_attr(attr, "transcript_id")]
tra2a[, transcript_version := grab_attr(attr, "transcript_version")]
tra2a[, row_id := paste0(transcript_id, ".", transcript_version)]
tra2a[, transcript_name := grab_attr(attr, "transcript_name")]
tra2a[, transcript_biotype := grab_attr(attr, "transcript_biotype")]
tra2a[, protein_id := grab_attr(attr, "protein_id")]
tra2a[, protein_version := grab_attr(attr, "protein_version")]
tra2a[, protein_row_id := paste0(protein_id, ".", protein_version)]
tra2a[, ccds_id := grab_attr(attr, "ccds_id")]
tra2a[, tag := attr]

tx_info <- unique(tra2a[row_id %in% txs & feature == "transcript", .(
  row_id, transcript_name, transcript_biotype, chr, start, end, strand, ccds_id, tag
)])

cds_info <- tra2a[row_id %in% txs & feature == "CDS", .(
  cds_bp = sum(end - start + 1L),
  n_cds_exons = .N,
  protein_id = na.omit(protein_row_id)[1],
  ccds_id = na.omit(ccds_id)[1]
), by = .(row_id, transcript_name)]

cds <- readDNAStringSet(cds_fasta)
names(cds) <- sub(" .*", "", names(cds))
cds <- cds[txs]
aa <- translate(cds, if.fuzzy.codon = "X")
aa_seq <- as.character(aa)
aa_seq <- sub("\\*$", "", aa_seq)
names(aa_seq) <- names(txs)

summary <- data.table(
  transcript_name = names(txs),
  row_id = unname(txs),
  cds_nt = width(cds),
  aa_length = nchar(aa_seq),
  protein_sequence = aa_seq
)

compare_to_ref <- function(query_name, ref_name = "TRA2A-201") {
  ref <- strsplit(aa_seq[[ref_name]], "")[[1]]
  qry <- strsplit(aa_seq[[query_name]], "")[[1]]
  max_len <- max(length(ref), length(qry))
  min_len <- min(length(ref), length(qry))
  aligned_same <- sum(ref[seq_len(min_len)] == qry[seq_len(min_len)])
  diffs <- data.table(position = seq_len(min_len), ref = ref[seq_len(min_len)], query = qry[seq_len(min_len)])
  diffs <- diffs[ref != query]
  data.table(
    query = query_name,
    ref = ref_name,
    ref_length = length(ref),
    query_length = length(qry),
    identical_prefix_positions = aligned_same,
    shared_positions = min_len,
    sequence_identity_over_shared_positions = aligned_same / min_len,
    exact_same_sequence = length(ref) == length(qry) && aligned_same == min_len,
    first_difference = if (nrow(diffs)) diffs$position[1] else NA_integer_,
    n_different_positions_shared = nrow(diffs),
    length_difference = length(qry) - length(ref)
  )
}

comparisons <- rbindlist(lapply(setdiff(names(txs), "TRA2A-201"), compare_to_ref))

cat("\nEnsembl transcript/CDS metadata:\n")
print(merge(tx_info[, !"tag"], cds_info, by = c("row_id", "transcript_name"), all = TRUE))

cat("\nProtein sequence summary translated from Ensembl CDS FASTA:\n")
print(summary[, .(transcript_name, row_id, cds_nt, aa_length)])

cat("\nComparison to TRA2A-201:\n")
print(comparisons)

cat("\nProtein sequences:\n")
for (nm in names(aa_seq)) {
  cat(">", nm, "|", txs[[nm]], "\n", aa_seq[[nm]], "\n", sep = "")
}
