# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- analysis_path()
psite_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
pathway_dir <- file.path(psite_dir, "Pathway_gProfiler_Clean_REAC_GOBP_GOCC_true_interaction_lfc0.7_rawP0.05", "Tables")
out_dir <- file.path(psite_dir, "RS_DS_similarity_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_cut <- 0.05
lfc_cut <- 0.7

contrasts <- c("Sensitive_Vin_vs_DMSO", "Resistant_Vin_vs_DMSO", "Vin_Resistant_vs_Sensitive", "Interaction")

read_fraction <- function(frac, contrast) {
  f <- file.path(psite_dir, paste0("Fraction_", frac), paste0(contrast, "_", frac, "_psite_limma_all_genes.csv"))
  if (!file.exists(f)) stop("Missing file: ", f)
  d <- fread(f)
  d[, gene_key := sub("\\.\\d+$", "", gene_id_clean)]
  d[, direction := fifelse(P.Value < p_cut & logFC >= lfc_cut, "Up",
                    fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS"))]
  d[, sig := direction %in% c("Up", "Down")]
  d[, .(gene_key, gene_name, logFC, P.Value, adj.P.Val, direction, sig)]
}

gene_summaries <- list()
overlap_tables <- list()

for (ct in contrasts) {
  rs <- read_fraction("RS", ct)
  ds <- read_fraction("DS", ct)
  setnames(rs, c("logFC", "P.Value", "adj.P.Val", "direction", "sig"),
           c("RS_logFC", "RS_P.Value", "RS_adj.P.Val", "RS_direction", "RS_sig"))
  setnames(ds, c("gene_name", "logFC", "P.Value", "adj.P.Val", "direction", "sig"),
           c("gene_name_DS", "DS_logFC", "DS_P.Value", "DS_adj.P.Val", "DS_direction", "DS_sig"))
  joined <- merge(rs, ds, by = "gene_key", all = FALSE)
  joined[, gene_name := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_name_DS)]

  ok <- is.finite(joined$RS_logFC) & is.finite(joined$DS_logFC)
  pear <- suppressWarnings(cor.test(joined$RS_logFC[ok], joined$DS_logFC[ok], method = "pearson"))
  spear <- suppressWarnings(cor.test(joined$RS_logFC[ok], joined$DS_logFC[ok], method = "spearman", exact = FALSE))

  n_rs_sig <- joined[RS_sig == TRUE, .N]
  n_ds_sig <- joined[DS_sig == TRUE, .N]
  n_both <- joined[RS_sig == TRUE & DS_sig == TRUE, .N]
  n_same <- joined[RS_sig == TRUE & DS_sig == TRUE & RS_direction == DS_direction, .N]
  n_opp <- joined[RS_sig == TRUE & DS_sig == TRUE & RS_direction != DS_direction, .N]
  union_n <- joined[RS_sig == TRUE | DS_sig == TRUE, .N]

  gene_summaries[[ct]] <- data.table(
    contrast = ct,
    n_genes_compared = sum(ok),
    pearson_r = unname(pear$estimate),
    pearson_p = pear$p.value,
    spearman_rho = unname(spear$estimate),
    spearman_p = spear$p.value,
    RS_sig = n_rs_sig,
    DS_sig = n_ds_sig,
    both_sig = n_both,
    same_direction = n_same,
    opposite_direction = n_opp,
    jaccard_sig = ifelse(union_n == 0, NA_real_, n_both / union_n),
    RS_sig_recovered_by_DS_pct = ifelse(n_rs_sig == 0, NA_real_, 100 * n_both / n_rs_sig),
    DS_sig_recovered_by_RS_pct = ifelse(n_ds_sig == 0, NA_real_, 100 * n_both / n_ds_sig)
  )

  joined[, contrast := ct]
  overlap_tables[[ct]] <- joined[RS_sig == TRUE | DS_sig == TRUE, .(
    contrast, gene_name, gene_key,
    RS_logFC, RS_P.Value, RS_direction,
    DS_logFC, DS_P.Value, DS_direction,
    overlap_class = fifelse(RS_sig & DS_sig & RS_direction == DS_direction, "both_same_direction",
                     fifelse(RS_sig & DS_sig & RS_direction != DS_direction, "both_opposite_direction",
                     fifelse(RS_sig, "RS_only", "DS_only")))
  )]
}

gene_summary <- rbindlist(gene_summaries)
gene_overlap <- rbindlist(overlap_tables)
fwrite(gene_summary, file.path(out_dir, "RS_DS_gene_level_similarity_summary.csv"))
fwrite(gene_overlap, file.path(out_dir, "RS_DS_gene_level_sig_overlap_details.csv"))

scatter_dt <- rbindlist(lapply(contrasts, function(ct) {
  rs <- read_fraction("RS", ct)
  ds <- read_fraction("DS", ct)
  setnames(rs, "logFC", "RS_logFC")
  setnames(ds, c("gene_name", "logFC"), c("gene_name_DS", "DS_logFC"))
  j <- merge(rs[, .(gene_key, gene_name, RS_logFC, RS_direction = direction, RS_sig = sig)],
             ds[, .(gene_key, gene_name_DS, DS_logFC, DS_direction = direction, DS_sig = sig)],
             by = "gene_key")
  j[, contrast := ct]
  j[, class := fifelse(RS_sig & DS_sig & RS_direction == DS_direction, "Both sig same",
                fifelse(RS_sig & DS_sig & RS_direction != DS_direction, "Both sig opposite",
                fifelse(RS_sig, "RS only",
                fifelse(DS_sig, "DS only", "Neither"))))]
  j
}), fill = TRUE)

scatter_dt[, gene_name := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_name_DS)]
scatter_dt[, is_histone := grepl("^H(1|2A|2B|3|4)[A-Za-z0-9]+$", gene_name)]
histone_dt <- scatter_dt[is_histone == TRUE]
scatter_no_histone_dt <- scatter_dt[is_histone != TRUE]
histone_summary <- histone_dt[, .N, by = .(contrast, class)][order(contrast, class)]
fwrite(histone_dt, file.path(out_dir, "RS_DS_gene_level_logFC_concordance_histone_genes.csv"))
fwrite(histone_summary, file.path(out_dir, "RS_DS_gene_level_logFC_concordance_histone_summary.csv"))

base_scatter <- ggplot(scatter_no_histone_dt, aes(x = RS_logFC, y = DS_logFC)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_point(aes(color = class), alpha = 0.6, size = 0.8) +
  facet_wrap(~contrast, scales = "free") +
  scale_color_manual(values = c("Both sig same" = "#2F855A", "Both sig opposite" = "#C53030",
                                "RS only" = "#805AD5", "DS only" = "#DD6B20", "Neither" = "grey80")) +
  labs(
    title = "RS vs DS P-site limma count concordance",
    subtitle = "Gene-level logFC comparison; histone genes removed; significance = raw P < 0.05 and |logFC| >= 0.7",
    x = "RS limma logFC",
    y = "DS limma logFC",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

p <- base_scatter +
  geom_point(
    data = histone_dt,
    aes(fill = "Histone genes"),
    shape = 21,
    color = "black",
    stroke = 0.35,
    alpha = 0.95,
    size = 1.8
  ) +
  scale_fill_manual(values = c("Histone genes" = "#00A9E0")) +
  labs(
    title = "RS vs DS P-site limma count concordance",
    subtitle = "Gene-level logFC comparison; significance = raw P < 0.05 and |logFC| >= 0.7",
    x = "RS limma logFC",
    y = "DS limma logFC",
    color = NULL,
    fill = NULL
  )

ggsave(file.path(out_dir, "RS_DS_gene_level_logFC_concordance_scatter.png"), base_scatter, width = 10.5, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "RS_DS_gene_level_logFC_concordance_scatter_no_histones.png"), base_scatter, width = 10.5, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "RS_DS_gene_level_logFC_concordance_scatter_histone_highlight.png"), p, width = 10.5, height = 8, dpi = 300, bg = "white")

pathway_file <- file.path(pathway_dir, "gprofiler_clean_filtered_terms_termSizeLE500_hitGE5.csv")
path_summary <- data.table()
path_common <- data.table()

if (file.exists(pathway_file)) {
  pathways <- fread(pathway_file)
  compare_path <- function(dirn) {
    rs <- unique(pathways[fraction == "RS" & direction == dirn, .(source, term_id, term_name)])
    ds <- unique(pathways[fraction == "DS" & direction == dirn, .(source, term_id, term_name)])
    rsk <- rs[, .(source, term_id)]
    dsk <- ds[, .(source, term_id)]
    common <- merge(rs, ds[, .(source, term_id, DS_term_name = term_name)], by = c("source", "term_id"))
    union <- unique(rbind(rsk, dsk))
    list(
      summary = data.table(
        direction = dirn,
        RS_terms = nrow(rsk),
        DS_terms = nrow(dsk),
        shared_terms = nrow(common),
        jaccard = ifelse(nrow(union) == 0, NA_real_, nrow(common) / nrow(union)),
        RS_terms_recovered_by_DS_pct = ifelse(nrow(rsk) == 0, NA_real_, 100 * nrow(common) / nrow(rsk)),
        DS_terms_recovered_by_RS_pct = ifelse(nrow(dsk) == 0, NA_real_, 100 * nrow(common) / nrow(dsk))
      ),
      common = common[, .(direction = dirn, source, term_id, term_name)]
    )
  }
  tmp <- lapply(c("Up", "Down"), compare_path)
  path_summary <- rbindlist(lapply(tmp, `[[`, "summary"))
  path_common <- rbindlist(lapply(tmp, `[[`, "common"), fill = TRUE)
  fwrite(path_summary, file.path(out_dir, "RS_DS_pathway_level_similarity_summary.csv"))
  fwrite(path_common, file.path(out_dir, "RS_DS_pathway_level_shared_terms.csv"))
}

cat("\nRS vs DS gene-level similarity:\n")
print(gene_summary[, .(
  contrast, n_genes_compared,
  pearson_r = round(pearson_r, 3),
  spearman_rho = round(spearman_rho, 3),
  RS_sig, DS_sig, both_sig, same_direction, opposite_direction,
  jaccard_sig = round(jaccard_sig, 3),
  RS_recovered_by_DS_pct = round(RS_sig_recovered_by_DS_pct, 1),
  DS_recovered_by_RS_pct = round(DS_sig_recovered_by_RS_pct, 1)
)])

cat("\nRS vs DS pathway-level similarity for true interaction g:Profiler filtered terms:\n")
if (nrow(path_summary) == 0) {
  cat("No pathway comparison available.\n")
} else {
  print(path_summary[, .(
    direction, RS_terms, DS_terms, shared_terms,
    jaccard = round(jaccard, 3),
    RS_recovered_by_DS_pct = round(RS_terms_recovered_by_DS_pct, 1),
    DS_recovered_by_RS_pct = round(DS_terms_recovered_by_RS_pct, 1)
  )])
  cat("\nShared RS/DS pathway terms:\n")
  if (nrow(path_common) == 0) cat("None\n") else print(path_common)
}

cat("\nSaved outputs to:\n", out_dir, "\n", sep = "")
