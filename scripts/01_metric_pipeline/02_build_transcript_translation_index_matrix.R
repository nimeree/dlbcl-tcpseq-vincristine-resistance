# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(riboWaltz)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(rtracklayer)
  library(Rsamtools)
  library(tools)
})

# -----------------------------
# User-provided paths
# -----------------------------
BAM_DIR <- normalizePath(input_path("cDNA", "t2g_v3", "Trasncript"), winslash = "/")
GTF <- normalizePath(input_path("Homo_sapiens.GRCh38.114.chr.gtf"), winslash = "/")
OFFCSV <- normalizePath(external_path("Thesis", "Analysis", "t2g_v2", "01_psite_offsets", "psite_offsets_ALL_samples.csv"), winslash = "/")
OUTDIR <- analysis_path("Translation_indexes_fixed")

BASELINE_FILE <- external_path("Thesis", "Analysis", "t2g_v1", "Version_5", "baseline_RNA_CPM_by_gene_SUDHL8.csv")

OUT_MATRIX <- file.path(OUTDIR, "transcript_psite_matrix_long_ALL_samples.csv")
OUT_METRICS <- file.path(OUTDIR, "transcript_translation_metrics_ALL_samples.csv")
OUT_METRICS_RNA <- file.path(OUTDIR, "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")
OUT_QA <- file.path(OUTDIR, "translation_index_QA_summary.csv")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(BAM_DIR), file.exists(GTF), file.exists(OFFCSV))

# -----------------------------
# Parameters
# -----------------------------
LEN_SSU <- 15:33
LEN_RS <- 15:33
LEN_DS <- 40:65
START_WIN_NT <- 90
STOP_WIN_NT <- 60
eps <- 1e-6
MIN_UTR5_FOR_METRICS <- 50
MIN_CDS_FOR_CORE_METRICS <- START_WIN_NT + STOP_WIN_NT + 1
MIN_CDS_COUNTS_STABLE <- 20
MIN_CORE_COUNTS_STABLE <- 20

# -----------------------------
# Helpers
# -----------------------------
fraction_from_name <- function(x) {
  u <- toupper(basename(x))
  if (grepl("SSU", u)) return("SSU")
  if (grepl("RS", u)) return("RS")
  if (grepl("DS", u)) return("DS")
  "NA"
}

length_range_for_fraction <- function(frac) {
  if (frac == "SSU") return(LEN_SSU)
  if (frac == "RS") return(LEN_RS)
  if (frac == "DS") return(LEN_DS)
  integer(0)
}

sample_pair_key <- function(x) {
  y <- tools::file_path_sans_ext(basename(x))
  # Remove fraction+lane token so RS/DS/SSU replicates share the same pair key.
  y <- gsub("-(SSU|RS|DS)_S[0-9]+", "", y, ignore.case = TRUE, perl = TRUE)
  y <- gsub("_(SSU|RS|DS)_S[0-9]+", "", y, ignore.case = TRUE, perl = TRUE)
  y <- gsub("(^|[-_])(SSU|RS|DS)($|[-_])", "\\1\\3", y, ignore.case = TRUE, perl = TRUE)
  y <- gsub("__+", "_", y)
  y <- gsub("--+", "-", y)
  y <- gsub("^[-_]+|[-_]+$", "", y)
  trimws(y, which = "both")
}

clean_id <- function(x) {
  y <- as.character(x)
  y <- trimws(y)
  y <- sub("^gene:", "", y)
  y <- sub("\\.\\d+$", "", y)
  y
}

link_or_copy <- function(src, dst) {
  ok <- tryCatch(file.link(src, dst), warning = function(w) FALSE, error = function(e) FALSE)
  if (!ok) ok <- file.copy(src, dst, overwrite = TRUE)
  ok
}

pick_offset_column <- function(dt) {
  nms <- names(dt)
  candidates <- c(
    "psite_offset_median", "psite_offset", "offset_from_5",
    "corrected_offset_from_5", "offset", "offset_5p"
  )
  exact <- candidates[candidates %in% nms]
  if (length(exact)) return(exact[1])
  fuzzy <- grep("psite.*off|offset", nms, ignore.case = TRUE, value = TRUE)
  if (!length(fuzzy)) stop("Could not find an offset column in OFFCSV.")
  fuzzy[1]
}

pick_length_column <- function(dt) {
  nms <- names(dt)
  candidates <- c("length", "read_length", "rlength", "readlen")
  exact <- candidates[candidates %in% nms]
  if (length(exact)) return(exact[1])
  fuzzy <- grep("length|read.*len", nms, ignore.case = TRUE, value = TRUE)
  if (!length(fuzzy)) stop("Could not find a read-length column in OFFCSV.")
  fuzzy[1]
}

pick_sample_column <- function(dt) {
  nms <- names(dt)
  candidates <- c("sample", "sample_name", "bam", "library")
  exact <- candidates[candidates %in% nms]
  if (length(exact)) return(exact[1])
  fuzzy <- grep("sample|library|bam", nms, ignore.case = TRUE, value = TRUE)
  if (!length(fuzzy)) return(NA_character_)
  fuzzy[1]
}

offsets_for_sample <- function(offset_dt, sname) {
  len_col <- pick_length_column(offset_dt)
  off_col <- pick_offset_column(offset_dt)
  smp_col <- pick_sample_column(offset_dt)
  dt <- copy(offset_dt)
  setnames(dt, len_col, "length")
  setnames(dt, off_col, "offset_raw")
  dt[, length := suppressWarnings(as.integer(length))]
  dt[, offset_raw := suppressWarnings(as.numeric(offset_raw))]
  dt <- dt[!is.na(length) & !is.na(offset_raw)]
  if (!is.na(smp_col) && smp_col %in% names(dt)) {
    dt[, sample_col := as.character(get(smp_col))]
    hit <- dt[toupper(sample_col) == toupper(sname)]
    if (nrow(hit)) dt <- hit
  }
  out <- dt[, .(psite_offset_median = median(offset_raw, na.rm = TRUE)), by = length]
  out
}

import_one_bam <- function(bam, annot_for_rw, length_range) {
  sname <- tools::file_path_sans_ext(basename(bam))
  td <- file.path(tempdir(), paste0("rw_", sname))
  if (dir.exists(td)) unlink(td, recursive = TRUE, force = TRUE)
  dir.create(td, showWarnings = FALSE)
  bam_dst <- file.path(td, basename(bam))
  bai_src <- paste0(bam, ".bai")
  bai_dst <- paste0(bam_dst, ".bai")
  if (!link_or_copy(bam, bam_dst)) stop("Could not stage BAM: ", bam)
  if (file.exists(bai_src)) link_or_copy(bai_src, bai_dst)

  rl <- riboWaltz::bamtolist(bamfolder = td, annotation = annot_for_rw)
  names(rl) <- sname
  unlink(td, recursive = TRUE, force = TRUE)
  rl[[1]] <- rl[[1]] %>%
    dplyr::filter(dplyr::between(length, min(length_range), max(length_range)))
  rl
}

apply_offsets <- function(rl, off_by_len) {
  sname <- names(rl)[1]
  goff <- as.data.table(off_by_len)
  goff[, length := suppressWarnings(as.integer(length))]
  goff[, psite_offset_median := suppressWarnings(as.numeric(psite_offset_median))]
  goff <- goff[!is.na(length) & !is.na(psite_offset_median), .(length, psite_offset_median)]
  goff[, corrected_offset_from_5 := psite_offset_median]
  goff[, corrected_offset_from_3 := pmax(length - corrected_offset_from_5 - 1, 0)]
  goff[, `:=`(
    offset_from_5 = corrected_offset_from_5,
    offset_from_3 = corrected_offset_from_3,
    psite_offset = corrected_offset_from_5,
    offset = corrected_offset_from_5,
    sample = sname
  )]
  riboWaltz::psite_info(rl, goff)
}

count_psites_by_regions <- function(rl, annot_df) {
  dt <- as.data.table(rl[[1]])
  need <- c("transcript", "psite_from_start", "psite_from_stop")
  if (!all(need %in% names(dt))) {
    stop("psite_info columns missing; got: ", paste(names(dt), collapse = ", "))
  }
  dt <- merge(
    dt,
    as.data.table(annot_df[, c("transcript", "cds_length", "utr5_length", "utr3_length")]),
    by = "transcript", all.x = TRUE
  )
  dt[, in_utr5 := (psite_from_start < 0)]
  dt[, in_start := (psite_from_start >= 0 & psite_from_start < START_WIN_NT)]
  dt[, in_stop := (psite_from_stop < STOP_WIN_NT & psite_from_stop >= 0)]
  dt[, in_cds := (psite_from_start >= 0 & psite_from_stop >= 0)]
  dt[, in_core := (in_cds & !in_start & !in_stop)]
  region_qa <- dt[, .(
    reads_after_length_filter = .N,
    psite_from_start_min = min(psite_from_start, na.rm = TRUE),
    psite_from_start_median = median(psite_from_start, na.rm = TRUE),
    psite_from_start_max = max(psite_from_start, na.rm = TRUE),
    psite_from_stop_min = min(psite_from_stop, na.rm = TRUE),
    psite_from_stop_median = median(psite_from_stop, na.rm = TRUE),
    psite_from_stop_max = max(psite_from_stop, na.rm = TRUE),
    frac_utr5 = mean(in_utr5, na.rm = TRUE),
    frac_start = mean(in_start, na.rm = TRUE),
    frac_core = mean(in_core, na.rm = TRUE),
    frac_stop = mean(in_stop, na.rm = TRUE),
    frac_cds = mean(in_cds, na.rm = TRUE)
  )]
  sums <- dt[, .(
    n_utr5 = sum(in_utr5, na.rm = TRUE),
    n_start = sum(in_start, na.rm = TRUE),
    n_core = sum(in_core, na.rm = TRUE),
    n_stop = sum(in_stop, na.rm = TRUE),
    n_cds = sum(in_cds, na.rm = TRUE)
  ), by = .(transcript)]
  txlens <- as.data.table(annot_df[, c("transcript", "cds_length", "utr5_length")])
  txlens[is.na(cds_length), cds_length := 0L]
  txlens[is.na(utr5_length), utr5_length := 0L]
  sums <- merge(sums, txlens, by = "transcript", all.x = TRUE)
  sums[, core_len_nt := pmax(cds_length - START_WIN_NT - STOP_WIN_NT, 1L)]
  sums[, utr5_len_nt := pmax(utr5_length, 1L)]
  attr(sums, "region_qa") <- region_qa
  sums
}

psite_codon_matrix_long <- function(rl) {
  dt <- as.data.table(rl[[1]])
  need <- c("transcript", "psite_from_start", "psite_from_stop")
  if (!all(need %in% names(dt))) return(data.table())
  dt <- dt[psite_from_start >= 0 & psite_from_stop >= 0]
  if (!nrow(dt)) return(data.table())
  dt[, codon_pos := as.integer(floor(psite_from_start / 3) + 1L)]
  dt <- dt[codon_pos > 0]
  dt[, .(psite_count = .N), by = .(transcript, codon_pos)]
}

# -----------------------------
# Annotation from GTF + auto ID mode across all BAMs
# -----------------------------
message("[Annot] Building annotation from GTF")
gtf <- rtracklayer::import(GTF)

tx_map <- as.data.frame(mcols(gtf)) %>%
  transmute(
    transcript = as.character(transcript_id),
    transcript_version = as.character(transcript_version),
    gene_id = as.character(gene_id),
    gene_name = as.character(gene_name)
  ) %>%
  filter(!is.na(transcript)) %>%
  distinct(transcript, .keep_all = TRUE)

gtf_dt <- as.data.table(gtf)
gtf_dt[, width_nt := width(gtf)]
gtf_dt[, `:=`(
  is_cds = toupper(type) == "CDS",
  is_utr5 = tolower(type) %in% c("five_prime_utr", "5utr", "five_prime_utr_region"),
  is_utr3 = tolower(type) %in% c("three_prime_utr", "3utr", "three_prime_utr_region")
)]
tx_len <- gtf_dt[!is.na(transcript_id), .(
  cds_length = sum(fifelse(is_cds, width_nt, 0L), na.rm = TRUE),
  utr5_length = sum(fifelse(is_utr5, width_nt, 0L), na.rm = TRUE),
  utr3_length = sum(fifelse(is_utr3, width_nt, 0L), na.rm = TRUE)
), by = .(transcript = as.character(transcript_id))]

annot2 <- tx_map %>%
  left_join(tx_len, by = "transcript") %>%
  mutate(transcript_versioned = if_else(
    !is.na(transcript_version) & transcript_version != "",
    paste0(transcript, ".", transcript_version),
    transcript
  ))

bam_files <- list.files(BAM_DIR, pattern = "\\.bam$", full.names = TRUE)
stopifnot(length(bam_files) > 0)

hdr_ids <- unique(unlist(lapply(bam_files, function(bf) names(scanBamHeader(bf)[[1]]$targets))))
annot_versioned <- annot2
annot_versioned$transcript <- annot_versioned$transcript_versioned
annot_unversioned <- annot2
annot_unversioned$transcript <- annot_unversioned$transcript
match_v <- mean(hdr_ids %in% annot_versioned$transcript)
match_u <- mean(hdr_ids %in% annot_unversioned$transcript)

message(sprintf("[Annot] Match rate (versioned): %.1f%%", 100 * match_v))
message(sprintf("[Annot] Match rate (unversioned): %.1f%%", 100 * match_u))

annotation_use <- if (match_u > match_v) annot_unversioned else annot_versioned
id_mode <- if (match_u > match_v) "UNversioned" else "VERSIONED"
message("[Annot] Using ", id_mode, " transcript IDs")

annot_rw <- annotation_use
if (!"l_utr5" %in% names(annot_rw) && "utr5_length" %in% names(annot_rw)) annot_rw$l_utr5 <- annot_rw$utr5_length
if (!"l_utr3" %in% names(annot_rw) && "utr3_length" %in% names(annot_rw)) annot_rw$l_utr3 <- annot_rw$utr3_length
if (!"l_cds" %in% names(annot_rw) && "cds_length" %in% names(annot_rw)) annot_rw$l_cds <- annot_rw$cds_length

tx2gene <- annotation_use %>%
  dplyr::select(transcript, gene_id, gene_name, cds_length, utr5_length, utr3_length)

# -----------------------------
# Metadata
# -----------------------------
meta <- data.frame(
  bam = bam_files,
  sample = tools::file_path_sans_ext(basename(bam_files)),
  fraction = vapply(bam_files, fraction_from_name, character(1)),
  pair_key = vapply(bam_files, sample_pair_key, character(1)),
  stringsAsFactors = FALSE
)
meta <- meta[meta$fraction %in% c("SSU", "RS", "DS"), , drop = FALSE]
stopifnot(nrow(meta) > 0)
message(sprintf("[Run] %d BAMs across SSU/RS/DS", nrow(meta)))

offset_dt <- fread(OFFCSV)

# -----------------------------
# MAIN LOOP: transcript-level counts + codon matrix
# -----------------------------
region_results <- list()
matrix_results <- list()
qa_results <- list()

for (i in seq_len(nrow(meta))) {
  bam <- meta$bam[i]
  sname <- meta$sample[i]
  frac <- meta$fraction[i]
  pair_key <- meta$pair_key[i]
  lenr <- length_range_for_fraction(frac)
  if (!length(lenr)) next
  message(sprintf("[%d/%d] %s (%s) lengths %d-%d", i, nrow(meta), sname, frac, min(lenr), max(lenr)))

  rl <- import_one_bam(bam, annot_rw, lenr)
  if (!nrow(rl[[1]])) next
  off_local <- offsets_for_sample(offset_dt, sname)
  rl <- apply_offsets(rl, off_local)
  if (!nrow(rl[[1]])) next

  tx_counts <- count_psites_by_regions(rl, annotation_use)
  region_qa <- attr(tx_counts, "region_qa")
  if (!is.null(region_qa)) {
    region_qa[, `:=`(sample = sname, fraction = frac, pair_key = pair_key)]
    qa_results[[length(qa_results) + 1]] <- region_qa
  }
  tx_counts$sample <- sname
  tx_counts$fraction <- frac
  tx_counts$pair_key <- pair_key
  region_results[[length(region_results) + 1]] <- tx_counts

  tx_matrix <- psite_codon_matrix_long(rl)
  if (nrow(tx_matrix)) {
    tx_matrix[, `:=`(sample = sname, fraction = frac, pair_key = pair_key)]
    matrix_results[[length(matrix_results) + 1]] <- tx_matrix
  }

  rm(rl, tx_counts, tx_matrix)
  gc()
}

counts_all <- rbindlist(region_results, fill = TRUE)
if (!nrow(counts_all)) stop("No transcript counts produced.")

# Apply biological transcript filters after read mapping to avoid losing mappable reads.
if (!"cds_length" %in% names(counts_all) || !"utr5_length" %in% names(counts_all)) {
  counts_all <- merge(
    counts_all,
    as.data.table(annotation_use[, c("transcript", "cds_length", "utr5_length")]),
    by = "transcript",
    all.x = TRUE
  )
}
counts_all <- counts_all[!is.na(cds_length) & cds_length > 0]
short_cds_rows <- counts_all[cds_length < MIN_CDS_FOR_CORE_METRICS, .N]
short_cds_transcripts <- counts_all[cds_length < MIN_CDS_FOR_CORE_METRICS, uniqueN(transcript)]
if (short_cds_rows > 0) {
  message(sprintf(
    "[Filter] Dropping %d rows across %d transcripts with CDS length < %d nt",
    short_cds_rows, short_cds_transcripts, MIN_CDS_FOR_CORE_METRICS
  ))
}
counts_all <- counts_all[cds_length >= MIN_CDS_FOR_CORE_METRICS]
if (!is.na(MIN_UTR5_FOR_METRICS) && MIN_UTR5_FOR_METRICS > 0) {
  counts_all <- counts_all[!is.na(utr5_length) & utr5_length >= MIN_UTR5_FOR_METRICS]
}

matrix_all <- rbindlist(matrix_results, fill = TRUE)
if (nrow(matrix_all)) {
  lib_sizes <- matrix_all[, .(lib_psites = sum(psite_count, na.rm = TRUE)), by = .(sample, fraction)]
  matrix_all <- merge(matrix_all, lib_sizes, by = c("sample", "fraction"), all.x = TRUE)
  matrix_all[, psite_cpm := (psite_count / pmax(lib_psites, 1)) * 1e6]
  fwrite(matrix_all, OUT_MATRIX)
  message("Wrote matrix: ", OUT_MATRIX)
} else {
  message("No codon-level matrix rows produced.")
}

# -----------------------------
# Transcript-level metrics
# -----------------------------
tx2gene_dt <- as.data.table(tx2gene[, c("transcript", "gene_id", "gene_name")])
counts_all <- merge(counts_all, tx2gene_dt, by = "transcript", all.x = TRUE)

rs_dt <- counts_all[fraction == "RS"]
ds_dt <- counts_all[fraction == "DS"]
ssu_dt <- counts_all[fraction == "SSU"]

pair_dups <- counts_all[, .(sample_n = uniqueN(sample)), by = .(fraction, pair_key, transcript)][sample_n > 1]
if (nrow(pair_dups) > 0) {
  fwrite(pair_dups, file.path(OUTDIR, "duplicate_fraction_pair_key_transcripts.csv"))
  stop(
    "Found non-1:1 fraction/pair_key/transcript mappings. ",
    "Wrote duplicate_fraction_pair_key_transcripts.csv; fix pair_key parsing or aggregate all fractions before computing metrics."
  )
}

lib_rs_core <- rs_dt[, .(rs_core_total = sum(n_core, na.rm = TRUE)), by = sample]
rs_dt <- merge(rs_dt, lib_rs_core, by = "sample", all.x = TRUE)
rs_dt[, rs_core_kb := n_core / pmax(core_len_nt / 1000, 1e-9)]
rs_dt[, `:=`(
  rs_core_cpm = n_core / pmax(rs_core_total / 1e6, 1e-9),
  rs_rate = rs_core_kb / pmax(rs_core_total / 1e6, 1e-9),
  init_index = n_start / pmax(n_cds, 1),
  stop_index = n_stop / pmax(n_cds, 1)
)]

ds_core <- ds_dt[, .(
  ds_core_kb = n_core / pmax(core_len_nt / 1000, 1e-9)
), by = .(pair_key, transcript)]
rs_dt <- merge(rs_dt, ds_core, by = c("pair_key", "transcript"), all.x = TRUE)
rs_dt[is.na(ds_core_kb), ds_core_kb := 0]
rs_dt[, collision_index := ds_core_kb / pmax(rs_core_kb, 1e-9)]

if (nrow(ssu_dt)) {
  ssu_utr5 <- ssu_dt[, .(
    ssu_utr5_kb = n_utr5 / pmax(utr5_len_nt / 1000, 1e-9)
  ), by = .(pair_key, transcript)]
  rs_dt <- merge(rs_dt, ssu_utr5, by = c("pair_key", "transcript"), all.x = TRUE)
  rs_dt[is.na(ssu_utr5_kb), ssu_utr5_kb := 0]
} else {
  rs_dt[, ssu_utr5_kb := 0]
}
rs_dt[, ssu_scanning_index := ssu_utr5_kb / pmax(rs_core_kb, 1e-9)]

rs_dt[, initiation_loading := ssu_scanning_index]
rs_dt[, initiation_clearance := 1 / pmax(init_index, eps)]
rs_dt[, initiation_rate_index := initiation_loading * initiation_clearance]
rs_dt[, elongation_throughput := rs_rate]
rs_dt[, elongation_efficiency := fifelse(
  is.na(ds_core_kb) | ds_core_kb == 0,
  NA_real_,
  rs_core_kb / ds_core_kb
)]
rs_dt[, elongation_rate_index := elongation_throughput * elongation_efficiency]
rs_dt[, total_translation_rate_proxy := initiation_rate_index * elongation_rate_index]

# Stabilized metrics to reduce ratio blow-up and discrete islands in downstream plots.
rs_dt[, `:=`(
  init_index_stable = (n_start + 1) / (n_cds + 3),
  collision_score = (ds_core_kb + 1e-3) / (rs_core_kb + 1e-3),
  scanning_score = (ssu_utr5_kb + 1e-3) / (rs_core_kb + 1e-3),
  elongation_score = (rs_core_kb + 1e-3) / (ds_core_kb + 1e-3)
)]
rs_dt[, initiation_clearance_stable := 1 / sqrt(pmax(init_index_stable, eps))]
rs_dt[, initiation_rate_index_stable := scanning_score * initiation_clearance_stable]
rs_dt[, elongation_rate_index_stable := rs_rate * elongation_score]
rs_dt[, total_translation_rate_proxy_stable := initiation_rate_index_stable * elongation_rate_index_stable]

# Reliability mask for stable metrics (keeps low-count rows but marks stable rates as NA).
rs_dt[, stable_mask := (n_cds >= MIN_CDS_COUNTS_STABLE & n_core >= MIN_CORE_COUNTS_STABLE)]
stable_cols <- c(
  "init_index_stable", "collision_score", "scanning_score",
  "elongation_score", "initiation_clearance_stable",
  "initiation_rate_index_stable", "elongation_rate_index_stable",
  "total_translation_rate_proxy_stable"
)
for (cc in stable_cols) {
  rs_dt[stable_mask == FALSE, (cc) := NA_real_]
}

metrics <- rs_dt[, .(
  sample, pair_key, transcript, gene_id, gene_name,
  n_utr5, n_start, n_core, n_stop, n_cds,
  rs_rate, init_index, stop_index, collision_index,
  rs_core_cpm, rs_core_kb, ds_core_kb, ssu_utr5_kb, ssu_scanning_index,
  stable_mask,
  initiation_clearance, initiation_rate_index, elongation_rate_index, total_translation_rate_proxy,
  init_index_stable, collision_score, scanning_score,
  elongation_score, initiation_clearance_stable, initiation_rate_index_stable,
  elongation_rate_index_stable, total_translation_rate_proxy_stable
)]
fwrite(metrics, OUT_METRICS)
message("Wrote metrics: ", OUT_METRICS)

qa_dt <- rbindlist(list(
  data.table(metric = "bam_count", value = nrow(meta)),
  data.table(metric = "short_cds_rows_dropped", value = short_cds_rows),
  data.table(metric = "short_cds_transcripts_dropped", value = short_cds_transcripts),
  data.table(metric = "min_cds_for_core_metrics_nt", value = MIN_CDS_FOR_CORE_METRICS),
  data.table(metric = "duplicate_fraction_pair_key_transcripts", value = nrow(pair_dups)),
  data.table(metric = "metric_rows_written", value = nrow(metrics)),
  data.table(metric = "stable_metric_rows", value = sum(!is.na(metrics$initiation_rate_index_stable))),
  data.table(metric = "total_proxy_formula", value = "initiation_rate_index * elongation_rate_index"),
  data.table(metric = "total_proxy_stable_formula", value = "initiation_rate_index_stable * elongation_rate_index_stable"),
  data.table(metric = "ribosome_efficiency_score_formula", value = "stable_mask TRUE only; log2((rs_rate + 1) / (baseline_cpm_line + 1))"),
  data.table(metric = "protein_output_score_formula", value = "stable_mask TRUE only; log2((baseline_cpm_line + 1) * (rs_rate + 1))"),
  data.table(metric = "initiation_clearance_note", value = "inverse start-proximal CDS ribosome density; stable metric uses sqrt damping")
), fill = TRUE)
if (length(qa_results) > 0) {
  sample_qa <- rbindlist(qa_results, fill = TRUE)
  fwrite(sample_qa, file.path(OUTDIR, "translation_index_sample_region_QA.csv"))
  qa_dt <- rbind(
    qa_dt,
    data.table(metric = "sample_region_QA_file", value = "translation_index_sample_region_QA.csv"),
    fill = TRUE
  )
}
fwrite(qa_dt, OUT_QA)
message("Wrote QA summary: ", OUT_QA)

# -----------------------------
# RNA baseline adjusted outputs
# -----------------------------
if (file.exists(BASELINE_FILE)) {
  base <- fread(BASELINE_FILE)
  req <- c("gene_id", "baseline_sensitive_cpm", "baseline_resistant_cpm")
  if (all(req %in% names(base))) {
    dt <- copy(metrics)
    dt[, gene_id_clean := clean_id(gene_id)]
    base[, gene_id_clean := clean_id(gene_id)]
    dt <- merge(
      dt,
      base[, .(gene_id_clean, baseline_sensitive_cpm, baseline_resistant_cpm)],
      by = "gene_id_clean",
      all.x = TRUE
    )
    dt[, cell_line := fifelse(
      grepl("^SU8R-", sample, ignore.case = TRUE) |
        grepl("^SU8-R-", sample, ignore.case = TRUE) |
        grepl("RES", sample, ignore.case = TRUE),
      "Resistant", "Sensitive"
    )]
    dt[, baseline_cpm_line := fifelse(cell_line == "Resistant", baseline_resistant_cpm, baseline_sensitive_cpm)]
    dt[, ribosome_efficiency_score := NA_real_]
    dt[!is.na(baseline_cpm_line) & stable_mask == TRUE, ribosome_efficiency_score := log2((rs_rate + 1) / (baseline_cpm_line + 1))]
    dt[, protein_output_score := NA_real_]
    dt[!is.na(baseline_cpm_line) & stable_mask == TRUE, protein_output_score := log2((baseline_cpm_line + 1) * (rs_rate + 1))]
    dt[, TE_base := pmax(ribosome_efficiency_score, eps)]
    dt[, initiation_rate_TEproxy := initiation_rate_index * TE_base]
    dt[, elongation_rate_TEproxy := elongation_rate_index * TE_base]
    dt[, total_translation_rate_TEproxy := total_translation_rate_proxy * TE_base]
    if ("initiation_rate_index_stable" %in% names(dt)) {
      dt[, initiation_rate_TEproxy_stable := initiation_rate_index_stable * TE_base]
    }
    if ("elongation_rate_index_stable" %in% names(dt)) {
      dt[, elongation_rate_TEproxy_stable := elongation_rate_index_stable * TE_base]
    }
    if ("total_translation_rate_proxy_stable" %in% names(dt)) {
      dt[, total_translation_rate_TEproxy_stable := total_translation_rate_proxy_stable * TE_base]
    }
    fwrite(dt, OUT_METRICS_RNA)
    message("Wrote RNA-adjusted metrics: ", OUT_METRICS_RNA)
  } else {
    message("Baseline file found but required columns missing; skipped RNA-adjusted output.")
  }
} else {
  message("Baseline RNA file not found; skipped RNA-adjusted output.")
}

message("Done. Outputs in: ", OUTDIR)
