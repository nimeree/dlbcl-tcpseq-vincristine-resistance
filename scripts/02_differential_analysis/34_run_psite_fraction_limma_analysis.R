# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(ggrepel)
})

base_dir <- analysis_path()
index_dir <- file.path(base_dir, "Translation_indexes_fixed")
out_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

psite_file <- file.path(index_dir, "transcript_psite_matrix_long_ALL_samples.csv")
metric_file <- file.path(index_dir, "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")

p_cut <- 0.05
lfc_cut <- 0.7
top_n_each <- 10

condition_from_sample <- function(sample) {
  fifelse(grepl("SU8R-DMSO|SU8-R-DMSO", sample), "Resistant_DMSO",
  fifelse(grepl("SU8R-Vin|SU8-R-Vin|SU8R-VIN|SU8-R-VIN", sample), "Resistant_Vin",
  fifelse(grepl("SU8-DMSO", sample), "Sensitive_DMSO",
  fifelse(grepl("SU8-Vin|SU8-VIN", sample), "Sensitive_Vin", NA_character_))))
}

cell_line_from_condition <- function(condition) {
  fifelse(grepl("^Resistant", condition), "Resistant", "Sensitive")
}

treatment_from_condition <- function(condition) {
  fifelse(grepl("Vin$", condition), "Vin", "DMSO")
}

replicate_from_sample <- function(sample) {
  out <- sub(".*_Rep([0-9]+).*", "\\1", sample)
  suppressWarnings(as.integer(fifelse(out == sample, NA_character_, out)))
}

make_volcano <- function(dt, title, out_png) {
  df <- copy(dt)
  df <- df[!is.na(P.Value) & !is.na(logFC)]
  df[, neglog10_p := -log10(P.Value)]
  df[, direction := fifelse(
    P.Value < p_cut & logFC >= lfc_cut, "Up",
    fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS")
  )]
  n_up <- df[direction == "Up", .N]
  n_down <- df[direction == "Down", .N]
  lab <- rbindlist(list(
    df[direction == "Up"][order(P.Value)][1:min(top_n_each, .N)],
    df[direction == "Down"][order(P.Value)][1:min(top_n_each, .N)]
  ), fill = TRUE)

  p <- ggplot(df, aes(logFC, neglog10_p)) +
    geom_point(aes(color = direction), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c(Up = "#D7191C", Down = "#2C7BB6", NS = "grey72")) +
    geom_hline(yintercept = -log10(p_cut), linetype = "dotted", linewidth = 0.7) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", linewidth = 0.7) +
    ggrepel::geom_text_repel(
      data = lab,
      aes(label = gene_name),
      size = 2.8,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25
    ) +
    annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.25,
             label = paste0("Up: ", n_up), color = "#D7191C", size = 3.2, fill = "white") +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.25,
             label = paste0("Down: ", n_down), color = "#2C7BB6", size = 3.2, fill = "white") +
    labs(
      title = title,
      subtitle = paste0("P-site offset gene counts; limma trend=TRUE; raw P < ", p_cut, ", |logFC| >= ", lfc_cut),
      x = "limma logFC",
      y = "-log10(P value)",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )

  ggsave(out_png, p, width = 7.2, height = 5.5, dpi = 300, bg = "white")
}

message("Reading transcript-to-gene map...")
tx_map <- unique(fread(metric_file, select = c("transcript", "gene_id_clean", "gene_name")))
tx_map <- tx_map[!is.na(gene_id_clean) & !is.na(gene_name) & gene_name != ""]

message("Reading P-site matrix and aggregating to gene/fraction/sample counts...")
psite <- fread(psite_file, select = c("sample", "fraction", "transcript", "psite_count"))
psite <- merge(psite, tx_map, by = "transcript", allow.cartesian = TRUE)
gene_counts_long <- psite[
  ,
  .(psite_count = sum(psite_count, na.rm = TRUE)),
  by = .(gene_id_clean, gene_name, sample, fraction)
]
rm(psite)

sample_meta <- unique(gene_counts_long[, .(sample, fraction)])
sample_meta[, condition := condition_from_sample(sample)]
sample_meta <- sample_meta[!is.na(condition)]
sample_meta[, cell_line := cell_line_from_condition(condition)]
sample_meta[, treatment := treatment_from_condition(condition)]
sample_meta[, replicate := replicate_from_sample(sample)]
sample_meta[, condition := factor(condition, levels = c("Sensitive_DMSO", "Sensitive_Vin", "Resistant_DMSO", "Resistant_Vin"))]

gene_counts_long <- merge(gene_counts_long, sample_meta[, .(sample, fraction, condition)], by = c("sample", "fraction"))

fwrite(gene_counts_long, file.path(out_dir, "psite_gene_counts_long_by_fraction_sample.csv"))
fwrite(sample_meta[order(fraction, condition, replicate)], file.path(out_dir, "psite_fraction_limma_sample_metadata.csv"))

run_fraction <- function(frac) {
  message("\n[Fraction] ", frac)
  frac_dir <- file.path(out_dir, paste0("Fraction_", frac))
  dir.create(frac_dir, recursive = TRUE, showWarnings = FALSE)

  d <- gene_counts_long[fraction == frac]
  meta <- unique(sample_meta[fraction == frac][order(condition, replicate)])
  sample_levels <- meta$sample

  mat_dt <- dcast(d, gene_id_clean + gene_name ~ sample, value.var = "psite_count", fill = 0)
  gene_info <- mat_dt[, .(gene_id_clean, gene_name)]
  mat <- as.matrix(mat_dt[, ..sample_levels])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(gene_info$gene_id_clean)

  keep <- rowSums(mat >= 10) >= 2
  mat <- mat[keep, , drop = FALSE]
  gene_info <- gene_info[keep]
  lib_sizes <- colSums(mat)
  log_cpm <- log2(t(t(mat + 0.5) / pmax(lib_sizes, 1)) * 1e6)

  group <- factor(meta$condition, levels = c("Sensitive_DMSO", "Sensitive_Vin", "Resistant_DMSO", "Resistant_Vin"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- lmFit(log_cpm, design)
  fit <- eBayes(fit, trend = TRUE)
  contrast_matrix <- makeContrasts(
    Sensitive_Vin_vs_DMSO = Sensitive_Vin - Sensitive_DMSO,
    Resistant_Vin_vs_DMSO = Resistant_Vin - Resistant_DMSO,
    Vin_Resistant_vs_Sensitive = Resistant_Vin - Sensitive_Vin,
    levels = design
  )
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE)

  fwrite(data.table(
    sample = colnames(mat),
    fraction = frac,
    lib_psites = as.numeric(lib_sizes),
    condition = as.character(group)
  ), file.path(frac_dir, paste0(frac, "_library_sizes_after_gene_filter.csv")))

  summary_rows <- list()
  for (contrast in colnames(contrast_matrix)) {
    tt <- as.data.table(topTable(fit2, coef = contrast, number = Inf, sort.by = "none"), keep.rownames = "row_key")
    tt[, row_index := seq_len(.N)]
    tt[, gene_id_clean := gene_info$gene_id_clean]
    tt[, gene_name := gene_info$gene_name]
    tt[, significant_rawP0.05_lfc0.7 := !is.na(P.Value) & P.Value < p_cut & abs(logFC) >= lfc_cut]
    tt[, direction := fifelse(
      P.Value < p_cut & logFC >= lfc_cut, "Up",
      fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS")
    )]
    setcolorder(tt, c("gene_id_clean", "gene_name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "significant_rawP0.05_lfc0.7", "direction"))
    tt <- tt[order(P.Value)]

    out_all <- file.path(frac_dir, paste0(contrast, "_", frac, "_psite_limma_all_genes.csv"))
    out_sig <- file.path(frac_dir, paste0(contrast, "_", frac, "_psite_limma_sig_rawP0.05_lfc0.7.csv"))
    out_up <- file.path(frac_dir, paste0(contrast, "_", frac, "_psite_limma_sig_up_rawP0.05_lfc0.7.csv"))
    out_down <- file.path(frac_dir, paste0(contrast, "_", frac, "_psite_limma_sig_down_rawP0.05_lfc0.7.csv"))
    fwrite(tt, out_all)
    fwrite(tt[significant_rawP0.05_lfc0.7 == TRUE], out_sig)
    fwrite(tt[direction == "Up"], out_up)
    fwrite(tt[direction == "Down"], out_down)

    make_volcano(
      tt,
      title = paste0(frac, " P-site limma: ", contrast),
      out_png = file.path(frac_dir, paste0(contrast, "_", frac, "_psite_limma_volcano.png"))
    )

    summary_rows[[contrast]] <- data.table(
      fraction = frac,
      contrast = contrast,
      genes_tested = nrow(tt),
      n_sig = tt[significant_rawP0.05_lfc0.7 == TRUE, .N],
      n_up = tt[direction == "Up", .N],
      n_down = tt[direction == "Down", .N],
      median_logFC = median(tt$logFC, na.rm = TRUE)
    )
  }

  summary <- rbindlist(summary_rows)
  fwrite(summary, file.path(frac_dir, paste0(frac, "_psite_limma_contrast_summary.csv")))
  summary
}

summary_all <- rbindlist(lapply(c("SSU", "RS", "DS"), run_fraction), fill = TRUE)
fwrite(summary_all, file.path(out_dir, "psite_fraction_limma_contrast_summary.csv"))

cat("\nP-site fraction limma complete.\n")
cat("Output directory:\n", out_dir, "\n", sep = "")
cat("\nSummary:\n")
print(summary_all[order(fraction, contrast)])
