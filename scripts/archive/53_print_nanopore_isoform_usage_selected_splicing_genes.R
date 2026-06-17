# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

rdata <- input_path("SUDHL.RData")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")

genes <- c("HNRNPA1", "SRSF2", "LSM4", "SNRNP70", "RBM3", "RBMX", "TRA2A", "SON")

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

exon_dt <- gtf_dt[feature == "exon"]
ex_attrs <- exon_dt$attributes
exon_dt[, `:=`(start_num = as.numeric(start), end_num = as.numeric(end))]
exon_map <- data.table(
  transcript_id = grab_attr(ex_attrs, "transcript_id"),
  transcript_version = grab_attr(ex_attrs, "transcript_version"),
  exon_len = exon_dt$end_num - exon_dt$start_num + 1
)
exon_map[, row_id := paste0(transcript_id, ".", transcript_version)]
tx_len <- exon_map[, .(exonic_length = sum(exon_len)), by = row_id]

tx_map <- merge(tx_map, tx_len, by = "row_id", all.x = TRUE)
tx_map <- tx_map[row_id %in% rownames(e$count_matrix)]

counts <- e$count_matrix[tx_map$row_id, samples, drop = FALSE]
counts[is.na(counts)] <- 0
count_dt <- as.data.table(as.table(counts))
setnames(count_dt, c("row_id", "Sample", "count"))
count_dt[, count := as.numeric(count)]
count_dt <- merge(count_dt, tx_map, by = "row_id", all.x = TRUE)
count_dt <- merge(count_dt, annotation[, .(Sample, Type, Replicate)], by = "Sample", all.x = TRUE)
count_dt <- count_dt[toupper(gene_name) %in% toupper(genes)]

cat("\nIsoform usage screen: summed Nanopore transcript counts by condition\n")
cat("Genes:", paste(genes, collapse = ", "), "\n\n")

summary_rows <- list()

for (g in genes) {
  d <- count_dt[toupper(gene_name) == toupper(g)]
  if (!nrow(d)) {
    cat("\n==============================\n")
    cat(g, "\n")
    cat("No transcripts found in count matrix.\n")
    next
  }

  iso <- d[, .(
    count = sum(count),
    exonic_length = suppressWarnings(max(exonic_length, na.rm = TRUE)),
    transcript_biotype = unique(transcript_biotype)[1]
  ), by = .(Type, row_id, transcript_name)]
  iso[, gene_total := sum(count), by = Type]
  iso[, proportion := fifelse(gene_total > 0, count / gene_total, NA_real_)]

  wide_count <- dcast(iso, row_id + transcript_name + transcript_biotype + exonic_length ~ Type, value.var = "count", fill = 0)
  if (!"Original" %in% names(wide_count)) wide_count[, Original := 0]
  if (!"Resistant" %in% names(wide_count)) wide_count[, Resistant := 0]
  wide_count[, total := Original + Resistant]
  wide_count <- wide_count[total > 0]
  setorder(wide_count, -total)

  wide_prop <- copy(wide_count)
  orig_total <- sum(wide_count$Original)
  res_total <- sum(wide_count$Resistant)
  wide_prop[, Original_prop := if (orig_total > 0) Original / orig_total else NA_real_]
  wide_prop[, Resistant_prop := if (res_total > 0) Resistant / res_total else NA_real_]
  wide_prop[, delta_prop := Resistant_prop - Original_prop]

  dominant_original <- wide_prop[which.max(Original_prop)]
  dominant_resistant <- wide_prop[which.max(Resistant_prop)]
  dominant_changed <- dominant_original$row_id != dominant_resistant$row_id

  mat <- as.matrix(wide_count[, .(Original, Resistant)])
  rownames(mat) <- wide_count$transcript_name
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  chi_p <- NA_real_
  fisher_p <- NA_real_
  if (nrow(mat) >= 2 && all(colSums(mat) > 0)) {
    chi_p <- suppressWarnings(chisq.test(mat)$p.value)
    if (nrow(mat) <= 8) {
      fisher_p <- tryCatch(fisher.test(mat, simulate.p.value = TRUE, B = 100000)$p.value, error = function(e) NA_real_)
    }
  }

  major_minor <- wide_prop[
    (Original_prop < 0.20 & Resistant_prop >= 0.20) |
      (Resistant_prop < 0.20 & Original_prop >= 0.20)
  ]

  nmd_shift <- wide_prop[grepl("nonsense_mediated_decay|NMD", transcript_biotype, ignore.case = TRUE)]
  nmd_orig <- if (nrow(nmd_shift)) sum(nmd_shift$Original) / orig_total else 0
  nmd_res <- if (nrow(nmd_shift)) sum(nmd_shift$Resistant) / res_total else 0

  shorter_flag <- dominant_resistant$exonic_length < dominant_original$exonic_length

  cat("\n==============================\n")
  cat(g, "\n")
  cat("Total counts: Original=", orig_total, "; Resistant=", res_total, "\n", sep = "")
  cat("Chi-square p=", signif(chi_p, 4), if (!is.na(fisher_p)) paste0("; Fisher simulated p=", signif(fisher_p, 4)) else "", "\n", sep = "")
  cat("Dominant Original: ", dominant_original$transcript_name, " (", dominant_original$Original, "/", orig_total,
      ", ", round(100 * dominant_original$Original_prop, 1), "%; biotype=", dominant_original$transcript_biotype,
      "; length=", dominant_original$exonic_length, ")\n", sep = "")
  cat("Dominant Resistant: ", dominant_resistant$transcript_name, " (", dominant_resistant$Resistant, "/", res_total,
      ", ", round(100 * dominant_resistant$Resistant_prop, 1), "%; biotype=", dominant_resistant$transcript_biotype,
      "; length=", dominant_resistant$exonic_length, ")\n", sep = "")
  cat("Dominant isoform changed: ", ifelse(dominant_changed, "YES", "NO"), "\n", sep = "")
  if (nmd_orig > 0 || nmd_res > 0) {
    cat("NMD-biotype proportion: Original=", round(100 * nmd_orig, 1), "%; Resistant=", round(100 * nmd_res, 1), "%\n", sep = "")
  }
  if (shorter_flag && dominant_changed) {
    cat("Length flag: resistant dominant isoform is shorter than original dominant isoform.\n")
  } else if (!shorter_flag && dominant_changed) {
    cat("Length flag: resistant dominant isoform is not shorter than original dominant isoform.\n")
  }

  cat("\nIsoform counts/proportions:\n")
  print(wide_prop[, .(
    transcript_name,
    row_id,
    transcript_biotype,
    exonic_length,
    Original,
    Original_prop = round(Original_prop, 3),
    Resistant,
    Resistant_prop = round(Resistant_prop, 3),
    delta_prop = round(delta_prop, 3)
  )])

  if (nrow(major_minor)) {
    cat("\nMinor-to-major or major-to-minor shifts at 20% threshold:\n")
    print(major_minor[, .(
      transcript_name,
      transcript_biotype,
      exonic_length,
      Original_prop = round(Original_prop, 3),
      Resistant_prop = round(Resistant_prop, 3),
      delta_prop = round(delta_prop, 3)
    )])
  } else {
    cat("\nNo isoform crossed the 20% minor/major threshold.\n")
  }

  summary_rows[[g]] <- data.table(
    gene = g,
    original_total = orig_total,
    resistant_total = res_total,
    chi_square_p = chi_p,
    fisher_sim_p = fisher_p,
    dominant_original = dominant_original$transcript_name,
    dominant_original_prop = dominant_original$Original_prop,
    dominant_resistant = dominant_resistant$transcript_name,
    dominant_resistant_prop = dominant_resistant$Resistant_prop,
    dominant_changed = dominant_changed,
    resistant_dominant_shorter = shorter_flag,
    nmd_original_prop = nmd_orig,
    nmd_resistant_prop = nmd_res,
    max_abs_delta_prop = max(abs(wide_prop$delta_prop), na.rm = TRUE)
  )
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
cat("\n==============================\n")
cat("Summary table\n")
print(summary_dt[, .(
  gene,
  original_total,
  resistant_total,
  chi_square_p = signif(chi_square_p, 4),
  dominant_original,
  dominant_original_pct = round(100 * dominant_original_prop, 1),
  dominant_resistant,
  dominant_resistant_pct = round(100 * dominant_resistant_prop, 1),
  dominant_changed,
  resistant_dominant_shorter,
  nmd_original_pct = round(100 * nmd_original_prop, 1),
  nmd_resistant_pct = round(100 * nmd_resistant_prop, 1),
  max_abs_delta_pct = round(100 * max_abs_delta_prop, 1)
)])
