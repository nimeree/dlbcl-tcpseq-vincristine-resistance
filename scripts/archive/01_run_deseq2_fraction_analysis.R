# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# Genome BAM DESeq2 per fraction (SSU/RS/DS) - LENGTH FILTERED
# RS/SSU: 15-33 nt, DS: 40-65 nt (from filtered BAMs)
#
# Outputs: <THESIS_ANALYSIS_DIR>
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
  library(Rsubread)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(rtracklayer)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

P_CUT <- 0.05
LFC_CUT <- 0.7
TOP_N_EACH <- 10

# -----------------------------
# Paths
# -----------------------------
BAM_DIR <- normalizePath(input_path("cDNA", "t2g_v3_lenFiltered"), winslash = "/")
GTF     <- normalizePath(input_path("Homo_sapiens.GRCh38.114.chr.gtf"), winslash = "/")

OUTBASE <- external_path("Thesis", "Analysis", "t2g_v3")
OUTDIR  <- file.path(OUTBASE, "DESeq2_genome_lenFiltered_lfc0.7")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(dir.exists(BAM_DIR), file.exists(GTF))

# -----------------------------
# Parse naming conventions
# -----------------------------
parse_sample <- function(fname) {
  b <- basename(fname)

  cell_line <- ifelse(grepl("^SU8R", b, ignore.case = TRUE), "Resistant", "Sensitive")
  treatment <- ifelse(grepl("-VIN-", b, ignore.case = TRUE), "Vin",
                      ifelse(grepl("-DMSO-", b, ignore.case = TRUE), "DMSO", NA))
  fraction  <- ifelse(grepl("-SSU_", b, ignore.case = TRUE), "SSU",
                      ifelse(grepl("-RS_", b, ignore.case = TRUE), "RS",
                             ifelse(grepl("-DS_", b, ignore.case = TRUE), "DS", NA)))
  rep <- str_match(b, "_Rep(\\d+)_")[,2]
  rep <- ifelse(is.na(rep), NA_integer_, as.integer(rep))

  sample_id <- tools::file_path_sans_ext(b)

  data.table(
    bam = fname,
    filename = b,
    sample = sample_id,
    cell_line = cell_line,
    treatment = treatment,
    fraction = fraction,
    rep = rep
  )
}

# -----------------------------
# Gene annotation map
# -----------------------------
message("[Annot] Loading GTF gene_id/gene_name map...")
gtf <- rtracklayer::import(GTF)
gtf_m <- as.data.frame(mcols(gtf))
gene_map <- unique(data.table(
  gene_id = as.character(gtf_m$gene_id),
  gene_name = as.character(gtf_m$gene_name)
))
gene_map <- gene_map[!is.na(gene_id)]
gene_map[, gene_id := sub("\\.\\d+$", "", gene_id)]
gene_map <- unique(gene_map, by = "gene_id")

message("[Annot] Loading org.Hs.eg.db gene descriptions...")
anno_keys <- unique(gene_map$gene_id)
org_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = anno_keys,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME", "GENETYPE")
)
org_map <- as.data.table(org_map)
setnames(org_map, "ENSEMBL", "gene_id")
first_non_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[1]
}
org_map <- org_map[, .(
  symbol = first_non_na(SYMBOL),
  gene_function = first_non_na(GENENAME),
  gene_type = first_non_na(GENETYPE)
), by = gene_id]

gene_map <- merge(gene_map, org_map, by = "gene_id", all.x = TRUE)
gene_map[, gene_name := fifelse(!is.na(gene_name) & gene_name != "", gene_name,
                                fifelse(!is.na(symbol) & symbol != "", symbol, gene_id))]
gene_map[, gene_function := fifelse(!is.na(gene_function) & gene_function != "", gene_function, NA_character_)]

# -----------------------------
# Find filtered BAMs
# -----------------------------
bam_files <- list.files(BAM_DIR, pattern = "\\.bam$", full.names = TRUE)
stopifnot(length(bam_files) > 0)

meta <- rbindlist(lapply(bam_files, parse_sample), fill = TRUE)
meta <- meta[!is.na(treatment) & !is.na(fraction) & !is.na(rep)]
stopifnot(nrow(meta) > 0)

message("[Info] Samples detected:")
print(meta[, .N, by = .(fraction, cell_line, treatment)][order(fraction, cell_line, treatment)])

# -----------------------------
# Significant gene export helper
# -----------------------------
write_significant_direction_csvs <- function(res_dt, out_prefix, out_dir,
                                             p_cut = P_CUT,
                                             lfc_cut = LFC_CUT) {
  df <- data.table::copy(as.data.table(res_dt))
  sig <- df[!is.na(pvalue) & !is.na(log2FoldChange) &
              pvalue < p_cut & abs(log2FoldChange) >= lfc_cut]

  sig_up <- sig[log2FoldChange >= lfc_cut][order(pvalue, -log2FoldChange)]
  sig_down <- sig[log2FoldChange <= -lfc_cut][order(pvalue, log2FoldChange)]

  fwrite(sig_up, file.path(out_dir, paste0(out_prefix, "_p05_lfc0.7_significant_up.csv")))
  fwrite(sig_down, file.path(out_dir, paste0(out_prefix, "_p05_lfc0.7_significant_down.csv")))
}

# -----------------------------
# Volcano helper
# -----------------------------
make_volcano <- function(res_dt, title, out_png,
                         p_cut = P_CUT,
                         lfc_cut = LFC_CUT,
                         top_n_each = TOP_N_EACH) {

  df <- data.table::copy(as.data.table(res_dt))
  df <- df[!is.na(pvalue) & !is.na(log2FoldChange)]

  df[, neglog10_p := -log10(pvalue)]
  df[, label := data.table::fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)]

  df[, direction := data.table::fifelse(
    pvalue < p_cut & log2FoldChange >=  lfc_cut, "Up",
    data.table::fifelse(pvalue < p_cut & log2FoldChange <= -lfc_cut, "Down", "NS")
  )]

  n_up <- nrow(df[direction == "Up"])
  n_down <- nrow(df[direction == "Down"])

  top_up <- df[direction == "Up"][order(pvalue)][1:min(top_n_each, .N)]
  top_dn <- df[direction == "Down"][order(pvalue)][1:min(top_n_each, .N)]
  lab <- data.table::rbindlist(list(top_up, top_dn), fill = TRUE)

  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10_p)) +
    geom_point(aes(color = direction), alpha = 0.75, size = 1.4) +
    scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
    geom_hline(yintercept = -log10(p_cut), linetype = "dotted", linewidth = 0.9) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", linewidth = 0.9) +
    ggrepel::geom_text_repel(
      data = lab,
      aes(label = label),
      size = 3.1,
      box.padding = 0.35,
      point.padding = 0.25,
      max.overlaps = Inf
    ) +
    annotate("label",
             x = Inf, y = Inf,
             hjust = 1.05, vjust = 1.3,
             label = paste0("Upregulated: ", n_up),
             color = "red", fill = "white", label.size = 0.25, size = 3.6) +
    annotate("label",
             x = -Inf, y = Inf,
             hjust = -0.05, vjust = 1.3,
             label = paste0("Downregulated: ", n_down),
             color = "blue", fill = "white", label.size = 0.25, size = 3.6) +
    labs(
      title = title,
      subtitle = paste0("p<", p_cut, " and |log2FC|>=", lfc_cut,
                        " | Up=", n_up, " Down=", n_down,
                        " | labels=top ", top_n_each, " up + top ", top_n_each, " down"),
      x = "log2 fold change",
      y = expression(-log[10](pvalue)),
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  ggsave(out_png, p, width = 6.8, height = 5.2, dpi = 160)
}

# -----------------------------
# PCA helper
# -----------------------------
save_pca <- function(dds, title, out_png) {
  vsd <- DESeq2::vst(dds, blind = TRUE)
  mat <- assay(vsd)
  pc <- prcomp(t(mat), scale. = FALSE)

  df <- as.data.frame(pc$x[, 1:2])
  df$sample <- rownames(df)
  cd <- as.data.frame(colData(dds))
  cd$sample <- rownames(cd)
  df <- merge(df, cd, by = "sample", all.x = TRUE)

  pv <- (pc$sdev^2) / sum(pc$sdev^2)
  pc1 <- round(pv[1] * 100, 1)
  pc2 <- round(pv[2] * 100, 1)

  p <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = sample), size = 2.7, max.overlaps = Inf) +
    labs(title = title,
         x = paste0("PC1 (", pc1, "%)"),
         y = paste0("PC2 (", pc2, "%)")) +
    theme_bw(base_size = 12)

  ggsave(out_png, p, width = 7.2, height = 5.6, dpi = 160)
}

# -----------------------------
# Run DESeq2 for a given subset
# -----------------------------
run_deseq_contrast <- function(meta_sub, design_formula, contrast_vec,
                               out_prefix, out_dir, gene_map, do_pca = TRUE) {

  fc <- Rsubread::featureCounts(
    files = meta_sub$bam,
    annot.ext = GTF,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    useMetaFeatures = TRUE,
    allowMultiOverlap = TRUE,
    nthreads = max(1, parallel::detectCores() - 1)
  )

  counts <- fc$counts
  rownames(counts) <- sub("\\.\\d+$", "", rownames(counts))
  colnames(counts) <- meta_sub$sample

  coldata <- as.data.frame(meta_sub[, .(sample, cell_line, treatment, fraction, rep)])
  rownames(coldata) <- meta_sub$sample

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design_formula)

  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)

  if (do_pca) {
    save_pca(dds,
             title = paste0("PCA (VST) - ", out_prefix),
             out_png = file.path(out_dir, paste0(out_prefix, "_PCA_vst.png")))
  }

  res <- DESeq2::results(dds, contrast = contrast_vec)
  res_dt <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")
  res_dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
  res_dt <- merge(res_dt, gene_map[, .(gene_id, gene_name, gene_function, gene_type)], by = "gene_id", all.x = TRUE)

  res_dt[, sig_p05_lfc0.7 := !is.na(pvalue) & pvalue < P_CUT & abs(log2FoldChange) >= LFC_CUT]
  res_dt[, direction := fifelse(!is.na(log2FoldChange) & log2FoldChange >= 0, "UP", "DOWN")]

  setcolorder(res_dt, c("gene_id","gene_name","gene_function","gene_type","log2FoldChange","lfcSE","stat","pvalue","padj","sig_p05_lfc0.7","direction"))

  out_csv <- file.path(out_dir, paste0(out_prefix, "_results_all.csv"))
  out_sig <- file.path(out_dir, paste0(out_prefix, "_p05_lfc0.7_sig.csv"))
  fwrite(res_dt, out_csv)
  fwrite(res_dt[sig_p05_lfc0.7 == TRUE], out_sig)
  write_significant_direction_csvs(res_dt, out_prefix, out_dir)

  out_png <- file.path(out_dir, paste0(out_prefix, "_volcano_pvalue_lfc0.7.png"))
  make_volcano(res_dt, title = out_prefix, out_png = out_png)

  invisible(list(dds = dds, res = res_dt))
}

# -----------------------------
# MAIN: run per fraction
# -----------------------------
for (frac in c("SSU","RS","DS")) {

  meta_f <- meta[fraction == frac]
  if (nrow(meta_f) == 0) next

  frac_dir <- file.path(OUTDIR, paste0("Fraction_", frac))
  dir.create(frac_dir, recursive = TRUE, showWarnings = FALSE)

  message("\n==============================")
  message("[Fraction] ", frac)
  message("==============================")

  message("[QC] PCA on all samples in fraction: ", frac)

  fc_all <- Rsubread::featureCounts(
    files = meta_f$bam,
    annot.ext = GTF,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    useMetaFeatures = TRUE,
    allowMultiOverlap = TRUE,
    nthreads = max(1, parallel::detectCores() - 1)
  )

  counts_all <- fc_all$counts
  rownames(counts_all) <- sub("\\.\\d+$", "", rownames(counts_all))
  colnames(counts_all) <- meta_f$sample

  cd_all <- as.data.frame(meta_f[, .(sample, cell_line, treatment, fraction, rep)])
  rownames(cd_all) <- meta_f$sample

  dds_all <- DESeq2::DESeqDataSetFromMatrix(countData = counts_all, colData = cd_all, design = ~ cell_line + treatment)
  keep_all <- rowSums(counts(dds_all) >= 10) >= 2
  dds_all <- dds_all[keep_all, ]
  dds_all <- DESeq2::DESeq(dds_all)

  save_pca(dds_all,
           title = paste0("PCA (VST) - ALL samples - ", frac),
           out_png = file.path(frac_dir, paste0("ALLsamples_", frac, "_PCA_vst.png")))

  meta_A <- meta_f[cell_line == "Sensitive" & treatment %in% c("DMSO","Vin")]
  if (nrow(meta_A) >= 4) {
    run_deseq_contrast(meta_A, ~ treatment, c("treatment","Vin","DMSO"),
                       paste0("Sensitive_Vin_vs_DMSO_", frac), frac_dir, gene_map)
  } else message("  [Skip] Sensitive Vin vs DMSO: not enough samples.")

  meta_B <- meta_f[cell_line == "Resistant" & treatment %in% c("DMSO","Vin")]
  if (nrow(meta_B) >= 4) {
    run_deseq_contrast(meta_B, ~ treatment, c("treatment","Vin","DMSO"),
                       paste0("Resistant_Vin_vs_DMSO_", frac), frac_dir, gene_map)
  } else message("  [Skip] Resistant Vin vs DMSO: not enough samples.")

  meta_C <- meta_f[treatment == "Vin" & cell_line %in% c("Sensitive","Resistant")]
  if (nrow(meta_C) >= 4) {
    run_deseq_contrast(meta_C, ~ cell_line, c("cell_line","Resistant","Sensitive"),
                       paste0("Vin_Resistant_vs_Sensitive_", frac), frac_dir, gene_map)
  } else message("  [Skip] Vin Resistant vs Sensitive: not enough samples.")
}

message("\nAll done. Outputs in: ", OUTDIR)
