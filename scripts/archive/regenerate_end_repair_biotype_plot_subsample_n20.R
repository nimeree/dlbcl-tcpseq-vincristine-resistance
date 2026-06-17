# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(Rsamtools)
  library(S4Vectors)
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(Rsubread)
})

set.seed(20260520)

n_pseudoreps <- 20L
subsample_fraction <- 0.10
force_resplit <- FALSE

bam_sources <- c(
  "Total RNA" = "/Volumes/T7/Initial_Test_Results/Bam/total.bam",
  "End-repaired" = "/Volumes/T7/Initial_Test_Results/Bam/Trimmed_end.bam",
  "Non end-repaired" = "/Volumes/T7/Initial_Test_Results/Bam/Trimmed_non.bam"
)

gtf_file <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
outdir <- analysis_path("biotype_endRepair_subsample_pseudoreps_n20_10pct")
pseudorep_dir <- file.path(outdir, "pseudo_replicate_bams")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(pseudorep_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(all(file.exists(bam_sources)), file.exists(gtf_file))

message("[1/5] Building gene biotype map from GTF")
gtf <- import(gtf_file)
gene_rows <- gtf[mcols(gtf)$type == "gene"]
gene_map <- tibble(
  gene_id = sub("\\.\\d+$", "", as.character(mcols(gene_rows)$gene_id)),
  gene_biotype = as.character(mcols(gene_rows)$gene_biotype)
) %>%
  distinct(gene_id, .keep_all = TRUE)

split_bam <- function(input_bam, label) {
  prefix <- gsub("[^A-Za-z0-9]+", "_", label)
  reps <- sprintf("ps%02d", seq_len(n_pseudoreps))
  output_bams <- file.path(pseudorep_dir, paste0(prefix, "_", reps, ".bam"))

  if (!force_resplit && all(file.exists(output_bams), file.exists(paste0(output_bams, ".bai")))) {
    message("[2/5] Reusing existing 10% pseudo-replicates for ", label)
    return(tibble(
      library_type = label,
      pseudorep = reps,
      bam = output_bams,
      bam_name = basename(output_bams)
    ))
  }

  message("[2/5] Creating ", n_pseudoreps, " x 10% pseudo-replicates for ", label)
  message("      reading qnames: ", input_bam)
  qnames <- scanBam(input_bam, param = ScanBamParam(what = "qname"))[[1]]$qname
  uq <- unique(qnames)
  sample_n <- max(1L, floor(length(uq) * subsample_fraction))
  rm(qnames)
  gc()

  if (any(file.exists(output_bams))) {
    unlink(c(output_bams, paste0(output_bams, ".bai")))
  }

  qname_sets_current <<- lapply(seq_len(n_pseudoreps), function(i) {
    sample(uq, size = sample_n, replace = FALSE)
  })
  rm(uq)
  gc()

  filters <- setNames(
    lapply(seq_len(n_pseudoreps), function(i) {
      substitute(qname %in% get("qname_sets_current", envir = .GlobalEnv)[[idx]], list(idx = i))
    }),
    reps
  )

  filterBam(
    input_bam,
    destination = output_bams,
    filter = filters,
    param = ScanBamParam(what = "qname"),
    indexDestination = TRUE
  )

  rm(qname_sets_current, envir = .GlobalEnv)
  gc()

  tibble(
    library_type = label,
    pseudorep = reps,
    bam = output_bams,
    bam_name = basename(output_bams)
  )
}

pseudorep_meta <- bind_rows(Map(split_bam, bam_sources, names(bam_sources)))
write_csv(pseudorep_meta, file.path(outdir, "subsample_pseudoreplicate_bams_n20_10pct.csv"))

message("[3/5] Counting assigned reads per gene with featureCounts")
fc <- featureCounts(
  files = pseudorep_meta$bam,
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

count_df <- as.data.frame(fc$counts, check.names = FALSE)
count_df$gene_id <- sub("\\.\\d+$", "", rownames(count_df))

long_counts <- count_df %>%
  pivot_longer(-gene_id, names_to = "bam_name", values_to = "assigned_count") %>%
  left_join(pseudorep_meta, by = "bam_name") %>%
  left_join(gene_map, by = "gene_id") %>%
  mutate(gene_biotype = if_else(is.na(gene_biotype) | gene_biotype == "", "unknown", gene_biotype))

biotype_counts <- long_counts %>%
  group_by(library_type, pseudorep, gene_biotype) %>%
  summarise(assigned_count = sum(assigned_count), .groups = "drop") %>%
  group_by(library_type, pseudorep) %>%
  mutate(total_assigned = sum(assigned_count), percent = 100 * assigned_count / total_assigned) %>%
  ungroup()

write_csv(biotype_counts, file.path(outdir, "RNA_biotype_endRepair_subsample_pseudorep_counts_long_n20_10pct.csv"))

message("[4/5] Summarising pseudo-replicates")
biotype_order <- c(
  "rRNA", "lncRNA", "protein_coding", "processed_pseudogene", "snoRNA",
  "misc_RNA", "snRNA", "Other", "unprocessed_pseudogene", "miRNA", "rRNA_pseudogene"
)

plot_counts <- biotype_counts %>%
  mutate(
    plot_biotype = if_else(gene_biotype %in% biotype_order, gene_biotype, "Other"),
    library_type = factor(library_type, levels = c("Total RNA", "Non end-repaired", "End-repaired"))
  ) %>%
  group_by(library_type, pseudorep, plot_biotype) %>%
  summarise(percent = sum(percent), assigned_count = sum(assigned_count), .groups = "drop")

summary_df <- plot_counts %>%
  group_by(library_type, plot_biotype) %>%
  summarise(
    mean_percent = mean(percent),
    sd_percent = sd(percent),
    n_pseudoreps = n(),
    mean_assigned_count = mean(assigned_count),
    .groups = "drop"
  ) %>%
  mutate(
    plot_biotype = factor(plot_biotype, levels = biotype_order),
    library_type = factor(library_type, levels = c("Total RNA", "Non end-repaired", "End-repaired"))
  )

write_csv(summary_df, file.path(outdir, "RNA_biotype_endRepair_subsample_pseudorep_summary_n20_10pct.csv"))

method_note <- tibble(
  n_pseudoreps = n_pseudoreps,
  subsample_fraction = subsample_fraction,
  error_bar = "SD",
  source_total = bam_sources[["Total RNA"]],
  source_end_repaired = bam_sources[["End-repaired"]],
  source_non_end_repaired = bam_sources[["Non end-repaired"]],
  gtf = gtf_file
)
write_csv(method_note, file.path(outdir, "RNA_biotype_endRepair_subsample_method_note.csv"))

message("[5/5] Plotting")
p <- ggplot(summary_df, aes(x = plot_biotype, y = mean_percent, fill = library_type)) +
  geom_col(position = position_dodge(width = 0.86), width = 0.72) +
  geom_errorbar(
    aes(ymin = pmax(mean_percent - sd_percent, 0), ymax = mean_percent + sd_percent),
    position = position_dodge(width = 0.86),
    width = 0.18,
    linewidth = 0.35
  ) +
  scale_fill_manual(values = c(
    "Total RNA" = "#F8766D",
    "Non end-repaired" = "#00BA38",
    "End-repaired" = "#619CFF"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "RNA biotype distribution in 80S-associated libraries (End-\nrepaired vs Non-end repaired) with NEB kit",
    x = "RNA biotype",
    y = "Percent of assigned counts (%)",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16, lineheight = 0.95),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(outdir, "RNA_biotype_endRepair_subsample_pseudoreps_n20_10pct_SD.png"), p, width = 11, height = 8, dpi = 300)
ggsave(file.path(outdir, "RNA_biotype_endRepair_subsample_pseudoreps_n20_10pct_SD.pdf"), p, width = 11, height = 8)

message("Done: ", outdir)
