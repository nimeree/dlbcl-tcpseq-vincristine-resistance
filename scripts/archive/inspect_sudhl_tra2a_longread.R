# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)

rdata <- input_path("SUDHL.RData")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")

e <- new.env()
load(rdata, envir = e)

gtf_lines <- readLines(gtf)
tra2a_tx_lines <- gtf_lines[
  grepl('\ttranscript\t', gtf_lines) &
    grepl('gene_id "ENSG00000164548"', gtf_lines)
]
attrs <- sub("^([^\t]*\t){8}", "", tra2a_tx_lines)

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

tx <- data.table(
  transcript_id = grab_attr(attrs, "transcript_id"),
  transcript_version = grab_attr(attrs, "transcript_version"),
  transcript_name = grab_attr(attrs, "transcript_name"),
  transcript_biotype = grab_attr(attrs, "transcript_biotype")
)
tx[, row_id := paste0(transcript_id, ".", transcript_version)]
tx[, in_matrix := row_id %in% rownames(e$count_matrix)]

cat("Objects in SUDHL.RData:\n")
print(data.table(
  object = ls(e),
  class = vapply(ls(e), function(x) paste(class(e[[x]]), collapse = "/"), character(1)),
  dim = vapply(ls(e), function(x) paste(dim(e[[x]]), collapse = "x"), character(1))
))

cat("\nAnnotation:\n")
print(e$Annotation)

cat("\nTRA2A annotated transcripts and matrix presence:\n")
print(tx)

present <- tx[in_matrix == TRUE]
if (nrow(present) > 0) {
  ann <- as.data.table(e$Annotation)

  summarize_matrix <- function(mat, value_name) {
    dt <- as.data.table(as.table(mat[present$row_id, , drop = FALSE]))
    setnames(dt, c("row_id", "Sample", value_name))
    dt <- merge(dt, ann, by = "Sample")
    dt[, .(
      mean = round(mean(get(value_name), na.rm = TRUE), 3),
      median = round(median(get(value_name), na.rm = TRUE), 3)
    ), by = .(row_id, Condition, Type)][order(row_id, Condition, Type)]
  }

  cat("\nMean transcript counts by cell line/type:\n")
  print(summarize_matrix(e$count_matrix, "count"))

  cat("\nMean fragment values by cell line/type:\n")
  print(summarize_matrix(e$frag_matrix, "frag"))

  cat("\nMean DTI values by cell line/type:\n")
  print(summarize_matrix(e$dti_matrix, "dti"))
}
