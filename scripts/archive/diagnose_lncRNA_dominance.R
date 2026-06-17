# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(rtracklayer)
  library(Rsubread)
})

gtf_file <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
outdir <- analysis_path("biotype_endRepair_subsample_pseudoreps_n20_10pct")

meta <- read_csv(file.path(outdir, "subsample_pseudoreplicate_bams_n20_10pct.csv"), show_col_types = FALSE)

fc <- featureCounts(
  files = meta$bam,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  allowMultiOverlap = FALSE,
  isPairedEnd = FALSE,
  nthreads = 4,
  verbose = FALSE
)

gtf <- import(gtf_file)
gene_rows <- gtf[mcols(gtf)$type == "gene"]
gene_map <- tibble(
  gene_id = sub("\\.\\d+$", "", as.character(mcols(gene_rows)$gene_id)),
  gene_name = as.character(mcols(gene_rows)$gene_name),
  gene_biotype = as.character(mcols(gene_rows)$gene_biotype)
) %>%
  distinct(gene_id, .keep_all = TRUE)

count_df <- as.data.frame(fc$counts, check.names = FALSE)
count_df$gene_id <- sub("\\.\\d+$", "", rownames(count_df))

long <- count_df %>%
  pivot_longer(-gene_id, names_to = "bam_name", values_to = "count") %>%
  left_join(meta, by = "bam_name") %>%
  left_join(gene_map, by = "gene_id")

ranked <- long %>%
  group_by(library_type, gene_id, gene_name, gene_biotype) %>%
  summarise(mean_count = mean(count), total_count = sum(count), .groups = "drop") %>%
  group_by(library_type) %>%
  mutate(lib_total = sum(mean_count), pct_of_library = 100 * mean_count / lib_total) %>%
  ungroup()

lnc <- ranked %>%
  filter(gene_biotype == "lncRNA") %>%
  group_by(library_type) %>%
  mutate(
    lnc_total = sum(mean_count),
    pct_of_lnc = 100 * mean_count / lnc_total
  ) %>%
  arrange(library_type, desc(mean_count)) %>%
  ungroup()

top_lnc <- lnc %>%
  group_by(library_type) %>%
  slice_head(n = 30) %>%
  ungroup()

lnc_concentration <- lnc %>%
  group_by(library_type) %>%
  summarise(
    n_lnc_genes_with_counts = sum(mean_count > 0),
    top1_gene = gene_name[which.max(mean_count)],
    top1_pct_of_lnc = max(pct_of_lnc),
    top1_pct_of_library = pct_of_library[which.max(mean_count)],
    top5_pct_of_lnc = sum(head(pct_of_lnc, 5)),
    top10_pct_of_lnc = sum(head(pct_of_lnc, 10)),
    .groups = "drop"
  )

write_csv(ranked, file.path(outdir, "diagnostic_gene_biotype_rankings_all.csv"))
write_csv(lnc, file.path(outdir, "diagnostic_lncRNA_gene_rankings.csv"))
write_csv(top_lnc, file.path(outdir, "diagnostic_lncRNA_top30_by_library.csv"))
write_csv(lnc_concentration, file.path(outdir, "diagnostic_lncRNA_concentration_summary.csv"))

print(lnc_concentration)
print(top_lnc, n = 90)
