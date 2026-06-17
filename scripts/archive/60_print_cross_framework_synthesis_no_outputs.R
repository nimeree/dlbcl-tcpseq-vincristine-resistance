# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# Print-only synthesis checks for p-site limma vs translation metric limma.
# No result tables or plots are written.

suppressPackageStartupMessages({
  library(data.table)
  library(gprofiler2)
})

BASE_DIR <- analysis_path()
PSITE_DIR <- file.path(BASE_DIR, "Psite_fraction_limma_lfc0.7_rawP0.05")
METRIC_DIR <- file.path(BASE_DIR, "Limma_translation_metrics_lfc0.7_rawP0.05")

P_CUT <- 0.05
LFC_CUT <- 0.7

pairs <- data.table(
  pair = c("Scanning_vs_SSU", "Ribosome_engagement_vs_RS", "Collision_vs_DS"),
  frac = c("SSU", "RS", "DS"),
  metric = c("scanning_score", "ribosome_efficiency_score", "collision_score")
)
contrasts <- c("Resistance_baseline", "Interaction")

read_psite <- function(frac, contrast) {
  f <- file.path(
    PSITE_DIR,
    paste0("Fraction_", frac),
    paste0(contrast, "_", frac, "_psite_limma_all_genes.csv")
  )
  d <- fread(f)
  if ("logFC" %in% names(d)) setnames(d, "logFC", "psite_logFC")
  if ("P.Value" %in% names(d)) setnames(d, "P.Value", "psite_P.Value")
  if ("adj.P.Val" %in% names(d)) setnames(d, "adj.P.Val", "psite_adj.P.Val")
  d[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  d[, psite_sig := psite_P.Value < P_CUT & abs(psite_logFC) >= LFC_CUT]
  d[, psite_dir := fifelse(psite_sig & psite_logFC > 0, "Up",
                    fifelse(psite_sig & psite_logFC < 0, "Down", "NS"))]
  unique(
    d[, .(gene_id_clean, gene_name, psite_logFC, psite_P.Value, psite_adj.P.Val, psite_sig, psite_dir)],
    by = "gene_id_clean"
  )
}

read_metric <- function(metric, contrast) {
  f <- file.path(METRIC_DIR, "Results", metric, paste0(contrast, "_limma_all_genes.csv"))
  d <- fread(f)
  if ("logFC" %in% names(d)) setnames(d, "logFC", "metric_logFC")
  if ("P.Value" %in% names(d)) setnames(d, "P.Value", "metric_P.Value")
  if ("adj.P.Val" %in% names(d)) setnames(d, "adj.P.Val", "metric_adj.P.Val")
  d[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  d[, metric_sig := metric_P.Value < P_CUT & abs(metric_logFC) >= LFC_CUT]
  d[, metric_dir := fifelse(metric_sig & metric_logFC > 0, "Up",
                     fifelse(metric_sig & metric_logFC < 0, "Down", "NS"))]
  unique(
    d[, .(gene_id_clean, gene_name, metric_logFC, metric_P.Value, metric_adj.P.Val, metric_sig, metric_dir)],
    by = "gene_id_clean"
  )
}

merge_pair <- function(frac, metric, contrast) {
  x <- merge(read_psite(frac, contrast), read_metric(metric, contrast),
             by = c("gene_id_clean", "gene_name"), all = TRUE)
  x[is.na(psite_sig), psite_sig := FALSE]
  x[is.na(metric_sig), metric_sig := FALSE]
  x[is.na(psite_dir), psite_dir := "NS"]
  x[is.na(metric_dir), metric_dir := "NS"]
  x[, same_dir := psite_sig & metric_sig & psite_dir == metric_dir]
  x[, opposite_dir := psite_sig & metric_sig & psite_dir != metric_dir]
  x[]
}

cat("\n1. CROSS-FRAMEWORK CONCORDANCE COUNTS\n")
conc <- rbindlist(lapply(contrasts, function(ct) {
  rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    x <- merge_pair(pairs$frac[i], pairs$metric[i], ct)
    data.table(
      contrast = ct,
      pair = pairs$pair[i],
      psite_only = sum(x$psite_sig & !x$metric_sig),
      metric_only = sum(x$metric_sig & !x$psite_sig),
      both_same_direction = sum(x$same_dir),
      both_opposite_direction = sum(x$opposite_dir),
      psite_sig_total = sum(x$psite_sig),
      metric_sig_total = sum(x$metric_sig)
    )
  }))
}))
print(conc)

same_genes <- function(frac, metric, contrast) {
  merge_pair(frac, metric, contrast)[
    same_dir == TRUE,
    .(gene_id_clean, gene_name, psite_dir, metric_dir)
  ]
}

robust_collision_ids <- intersect(
  same_genes("DS", "collision_score", "Resistance_baseline")$gene_id_clean,
  same_genes("DS", "collision_score", "Interaction")$gene_id_clean
)
robust_re_ids <- intersect(
  same_genes("RS", "ribosome_efficiency_score", "Resistance_baseline")$gene_id_clean,
  same_genes("RS", "ribosome_efficiency_score", "Interaction")$gene_id_clean
)

cat("\n2. ROBUST GENE SET SIZES\n")
cat("Robust collision genes:", length(robust_collision_ids), "\n")
cat("Robust ribosome engagement genes:", length(robust_re_ids), "\n")

get4 <- function(ids, frac, metric) {
  b <- merge_pair(frac, metric, "Resistance_baseline")[gene_id_clean %in% ids]
  i <- merge_pair(frac, metric, "Interaction")[gene_id_clean %in% ids]

  setnames(
    b,
    c("psite_logFC", "psite_P.Value", "psite_adj.P.Val",
      "metric_logFC", "metric_P.Value", "metric_adj.P.Val",
      "psite_dir", "metric_dir"),
    paste0(c("psite_logFC", "psite_P.Value", "psite_adj.P.Val",
             "metric_logFC", "metric_P.Value", "metric_adj.P.Val",
             "psite_dir", "metric_dir"), "_baseline")
  )
  setnames(
    i,
    c("psite_logFC", "psite_P.Value", "psite_adj.P.Val",
      "metric_logFC", "metric_P.Value", "metric_adj.P.Val",
      "psite_dir", "metric_dir"),
    paste0(c("psite_logFC", "psite_P.Value", "psite_adj.P.Val",
             "metric_logFC", "metric_P.Value", "metric_adj.P.Val",
             "psite_dir", "metric_dir"), "_interaction")
  )

  out <- merge(
    b[, .(
      gene_id_clean, gene_name,
      psite_logFC_baseline, psite_P.Value_baseline, psite_adj.P.Val_baseline,
      metric_logFC_baseline, metric_P.Value_baseline, metric_adj.P.Val_baseline,
      psite_dir_baseline, metric_dir_baseline
    )],
    i[, .(
      gene_id_clean, gene_name,
      psite_logFC_interaction, psite_P.Value_interaction, psite_adj.P.Val_interaction,
      metric_logFC_interaction, metric_P.Value_interaction, metric_adj.P.Val_interaction,
      psite_dir_interaction, metric_dir_interaction
    )],
    by = c("gene_id_clean", "gene_name"),
    all = TRUE
  )
  out[order(gene_name)]
}

fmt_num <- function(d) {
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  d[, (num) := lapply(.SD, function(z) signif(z, 3)), .SDcols = num]
  d
}

coll_det <- get4(robust_collision_ids, "DS", "collision_score")
re_det <- get4(robust_re_ids, "RS", "ribosome_efficiency_score")

cat("\n2a. ROBUST COLLISION GENE DETAILS\n")
print(fmt_num(copy(coll_det)))

cat("\n2b. ROBUST RIBOSOME ENGAGEMENT GENE DETAILS\n")
print(fmt_num(copy(re_det)))

heat_ids <- unique(c(robust_collision_ids, robust_re_ids))
all_gene_names <- unique(
  rbind(
    coll_det[, .(gene_id_clean, gene_name)],
    re_det[, .(gene_id_clean, gene_name)],
    fill = TRUE
  ),
  by = "gene_id_clean"
)
mat <- all_gene_names[gene_id_clean %in% heat_ids]

add_col <- function(dt, col, frac, metric, contrast, source) {
  x <- merge_pair(frac, metric, contrast)[, .(gene_id_clean, psite_logFC, metric_logFC)]
  valcol <- if (source == "psite") "psite_logFC" else "metric_logFC"
  x <- x[, .(gene_id_clean, val = get(valcol))]
  setnames(x, "val", col)
  merge(dt, x, by = "gene_id_clean", all.x = TRUE)
}

mat <- add_col(mat, "DS_psite_baseline", "DS", "collision_score", "Resistance_baseline", "psite")
mat <- add_col(mat, "Collision_score_baseline", "DS", "collision_score", "Resistance_baseline", "metric")
mat <- add_col(mat, "DS_psite_interaction", "DS", "collision_score", "Interaction", "psite")
mat <- add_col(mat, "Collision_score_interaction", "DS", "collision_score", "Interaction", "metric")
mat <- add_col(mat, "RS_psite_baseline", "RS", "ribosome_efficiency_score", "Resistance_baseline", "psite")
mat <- add_col(mat, "RE_score_baseline", "RS", "ribosome_efficiency_score", "Resistance_baseline", "metric")
mat <- add_col(mat, "RS_psite_interaction", "RS", "ribosome_efficiency_score", "Interaction", "psite")
mat <- add_col(mat, "RE_score_interaction", "RS", "ribosome_efficiency_score", "Interaction", "metric")
mat[, robust_set := fifelse(
  gene_id_clean %in% robust_collision_ids & gene_id_clean %in% robust_re_ids,
  "Collision+RE",
  fifelse(gene_id_clean %in% robust_collision_ids, "Collision", "Ribosome_engagement")
)]
setcolorder(
  mat,
  c("robust_set", "gene_id_clean", "gene_name",
    setdiff(names(mat), c("robust_set", "gene_id_clean", "gene_name")))
)

cat("\n3. HEATMAP LOGFC MATRIX VALUES (NO PLOT GENERATED)\n")
print(fmt_num(copy(mat[order(robust_set, gene_name)])))

run_ora <- function(ids, label) {
  cat("\n4. RELAXED ORA FOR ", label, " (", length(ids), " genes; min hit genes >= 3)\n", sep = "")
  if (length(ids) < 3) {
    cat("Too few genes.\n")
    return(invisible(NULL))
  }
  res <- tryCatch(
    gprofiler2::gost(
      query = ids,
      organism = "hsapiens",
      sources = c("REAC", "GO:BP", "GO:CC"),
      correction_method = "g_SCS",
      domain_scope = "annotated",
      user_threshold = 0.05,
      evcodes = TRUE
    ),
    error = function(e) e
  )
  if (inherits(res, "error") || is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
    cat("No g:SCS-significant terms.\n")
    return(invisible(NULL))
  }
  d <- as.data.table(res$result)
  d <- d[significant == TRUE & intersection_size >= 3 & term_size <= 500]
  if (nrow(d) == 0) {
    cat("No terms after relaxed filters.\n")
    return(invisible(NULL))
  }
  d <- d[order(p_value)][1:min(.N, 15),
                         .(source, term_name, p_value = signif(p_value, 3),
                           term_size, intersection_size, intersection)]
  print(d)
}

run_ora(robust_collision_ids, "ROBUST COLLISION")
run_ora(robust_re_ids, "ROBUST RIBOSOME ENGAGEMENT")

scan_special <- rbindlist(lapply(contrasts, function(ct) {
  scan <- read_metric("scanning_score", ct)[
    metric_sig == TRUE,
    .(gene_id_clean, gene_name, scan_logFC = metric_logFC, scan_dir = metric_dir)
  ]
  rs <- read_psite("RS", ct)[
    psite_sig == TRUE,
    .(gene_id_clean, rs_logFC = psite_logFC, rs_dir = psite_dir)
  ]
  x <- merge(scan, rs, by = "gene_id_clean", all.x = TRUE)
  x[, rs_sig := !is.na(rs_dir)]
  x[, rs_same_direction := rs_sig & scan_dir == rs_dir]
  x[, rs_opposite_direction := rs_sig & scan_dir != rs_dir]
  data.table(
    contrast = ct,
    scanning_sig_total = nrow(scan),
    scanning_up = sum(scan$scan_dir == "Up"),
    scanning_down = sum(scan$scan_dir == "Down"),
    also_RS_sig = sum(x$rs_sig),
    RS_same_direction = sum(x$rs_same_direction),
    RS_opposite_direction = sum(x$rs_opposite_direction),
    scan_up_RS_down = sum(x$scan_dir == "Up" & x$rs_dir == "Down", na.rm = TRUE),
    scan_down_RS_up = sum(x$scan_dir == "Down" & x$rs_dir == "Up", na.rm = TRUE),
    scan_up_RS_up = sum(x$scan_dir == "Up" & x$rs_dir == "Up", na.rm = TRUE),
    scan_down_RS_down = sum(x$scan_dir == "Down" & x$rs_dir == "Down", na.rm = TRUE)
  )
}))

cat("\n5. SCANNING SCORE SPECIAL CASE: OVERLAP WITH RS P-SITE SIGNIFICANCE\n")
print(scan_special)
