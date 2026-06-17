# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

METRIC_FILE <- analysis_path("Translation_indexes_fixed", "transcript_translation_metrics_with_RNA_baseline_ALL_samples.csv")
OUT_DIR <- analysis_path("Translation_indexes_fixed", "Validation_Plots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

download_hrt <- function() {
  isc_zip <- tempfile(fileext = ".zip")
  download.file(
    "https://github.com/XiaoZhangryy/iSC.MEB/archive/refs/heads/main.zip",
    isc_zip,
    quiet = TRUE,
    mode = "wb"
  )
  isc_dir <- tempfile()
  dir.create(isc_dir)
  unzip(isc_zip, files = "iSC.MEB-main/data/Human_HK_genes.rda", exdir = isc_dir)
  load(file.path(isc_dir, "iSC.MEB-main/data/Human_HK_genes.rda"))
  hrt <- as.data.table(Human_HK_genes)
  hrt[, tx := sub("\\.\\d+$", "", Ensembl)]
  hrt
}

download_paxdb <- function() {
  pax_url <- "https://pax-db.org/downloads/6.0/datasets/paxdb-abundance-files-v6.0/9606/9606-WHOLE_ORGANISM-integrated.txt"
  pax_file <- tempfile(fileext = ".txt")
  download.file(pax_url, pax_file, quiet = TRUE, mode = "wb")
  pax <- fread(pax_file, skip = "#gene_name")
  setnames(pax, c("gene_name", "string_external_id", "paxdb_abundance"))
  pax[, paxdb_abundance := as.numeric(paxdb_abundance)]
  pax[!is.na(gene_name) & gene_name != "", .(
    paxdb_abundance = median(paxdb_abundance, na.rm = TRUE)
  ), by = gene_name]
}

save_plot <- function(plot, name, width, height) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), plot, width = width, height = height, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), plot, width = width, height = height)
}

message("Loading metrics")
cols <- c(
  "sample", "transcript", "gene_id_clean", "gene_name",
  "ribosome_efficiency_score", "protein_output_score", "rs_core_cpm", "rs_rate", "baseline_cpm_line"
)
dt <- fread(METRIC_FILE, select = cols)
dt <- dt[!is.na(gene_id_clean) & gene_id_clean != "" & !is.na(gene_name) & gene_name != ""]
for (cc in c("ribosome_efficiency_score", "protein_output_score", "rs_core_cpm", "rs_rate", "baseline_cpm_line")) {
  dt[, (cc) := as.numeric(get(cc))]
}
dt[, transcript_clean := sub("\\.\\d+$", "", transcript)]

message("Loading reference sets")
hrt <- download_hrt()
pax_gene <- download_paxdb()
hrt_tx <- unique(hrt$tx)
hrt_translation_tx <- unique(hrt[grepl("^(RPL|RPS|MRPL|MRPS|EEF|EIF)", Gene), tx])

hrt_gene_ids <- unique(dt[transcript_clean %in% hrt_tx, gene_id_clean])
hrt_translation_gene_ids <- unique(dt[transcript_clean %in% hrt_translation_tx, gene_id_clean])

message("Collapsing transcript rows to gene-sample and gene-level summaries")
gene_sample <- dt[, .(
  ribosome_efficiency_score = median(ribosome_efficiency_score, na.rm = TRUE),
  protein_output_score = median(protein_output_score, na.rm = TRUE),
  rs_core_cpm = median(rs_core_cpm, na.rm = TRUE),
  rs_rate = median(rs_rate, na.rm = TRUE),
  baseline_cpm_line = median(baseline_cpm_line, na.rm = TRUE),
  gene_name = names(sort(table(gene_name), decreasing = TRUE))[1]
), by = .(sample, gene_id_clean)]
gene_sample[, housekeeping_group := fifelse(
  gene_id_clean %in% hrt_translation_gene_ids,
  "Translation/ribosomal HK",
  fifelse(gene_id_clean %in% hrt_gene_ids, "HRT housekeeping", "Background")
)]
gene_sample[, housekeeping_binary := fifelse(gene_id_clean %in% hrt_gene_ids, "HRT housekeeping", "Background")]

gene_level <- gene_sample[, .(
  ribosome_efficiency_score = median(ribosome_efficiency_score, na.rm = TRUE),
  protein_output_score = median(protein_output_score, na.rm = TRUE),
  rs_core_cpm = median(rs_core_cpm, na.rm = TRUE),
  rs_rate = median(rs_rate, na.rm = TRUE),
  baseline_cpm_line = median(baseline_cpm_line, na.rm = TRUE),
  gene_name = names(sort(table(gene_name), decreasing = TRUE))[1],
  housekeeping_group = names(sort(table(housekeeping_group), decreasing = TRUE))[1],
  housekeeping_binary = names(sort(table(housekeeping_binary), decreasing = TRUE))[1]
), by = gene_id_clean]
gene_level <- merge(gene_level, pax_gene, by = "gene_name", all.x = TRUE)

metric_labels <- c(
  protein_output_score = "Protein output score",
  ribosome_efficiency_score = "Ribosome efficiency score",
  rs_core_cpm = "RS core CPM",
  rs_rate = "RS rate",
  baseline_cpm_line = "RNA baseline CPM"
)

message("Plotting housekeeping distributions")
hk_long <- melt(
  gene_sample,
  id.vars = c("sample", "gene_id_clean", "housekeeping_binary"),
  measure.vars = names(metric_labels),
  variable.name = "metric",
  value.name = "value"
)
hk_long <- hk_long[is.finite(value)]
hk_long[, metric := factor(metric, levels = names(metric_labels), labels = metric_labels)]
hk_long[, housekeeping_binary := factor(housekeeping_binary, levels = c("Background", "HRT housekeeping"))]

p_hk <- ggplot(hk_long, aes(x = housekeeping_binary, y = value, fill = housekeeping_binary)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.72, linewidth = 0.25) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.92, linewidth = 0.25) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 7, color = "black") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Background" = "#A7B0BA", "HRT housekeeping" = "#2C7A7B")) +
  labs(
    title = "Housekeeping Genes Show Higher Protein Output Score and Ribosome Signal",
    subtitle = "Gene-sample medians; HRT Atlas housekeeping set vs background",
    x = NULL,
    y = "Metric value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank()
  )
save_plot(p_hk, "housekeeping_metric_distributions", 14, 5.5)

message("Plotting per-sample housekeeping effects")
sample_summary <- gene_sample[, .(
  median_protein_output_score = median(protein_output_score, na.rm = TRUE),
  median_ribosome_efficiency_score = median(ribosome_efficiency_score, na.rm = TRUE),
  median_rs_core_cpm = median(rs_core_cpm, na.rm = TRUE)
), by = .(sample, housekeeping_binary)]
sample_long <- melt(
  sample_summary,
  id.vars = c("sample", "housekeeping_binary"),
  variable.name = "metric",
  value.name = "median_value"
)
sample_long[, metric := factor(
  metric,
  levels = c("median_protein_output_score", "median_ribosome_efficiency_score", "median_rs_core_cpm"),
  labels = c("Protein output score", "Ribosome efficiency score", "RS core CPM")
)]
sample_long[, housekeeping_binary := factor(housekeeping_binary, levels = c("Background", "HRT housekeeping"))]

p_sample <- ggplot(sample_long, aes(x = sample, y = median_value, color = housekeeping_binary, group = housekeeping_binary)) +
  geom_point(size = 2.2) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Background" = "#68717D", "HRT housekeeping" = "#2C7A7B")) +
  labs(
    title = "Housekeeping Validation Is Consistent Across Samples",
    subtitle = "Median gene-level metric per RS sample",
    x = NULL,
    y = "Median value",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
save_plot(p_sample, "housekeeping_per_sample_medians", 11, 8)

message("Plotting PaxDB scatter panels")
pax_long <- melt(
  gene_level[is.finite(paxdb_abundance)],
  id.vars = c("gene_id_clean", "gene_name", "paxdb_abundance", "housekeeping_binary"),
  measure.vars = names(metric_labels),
  variable.name = "metric",
  value.name = "value"
)
pax_long <- pax_long[is.finite(value)]
pax_long[, metric := factor(metric, levels = names(metric_labels), labels = metric_labels)]
pax_long[, log10_paxdb := log10(paxdb_abundance + 1)]

p_pax <- ggplot(pax_long, aes(x = log10_paxdb, y = value)) +
  geom_bin2d(bins = 70) +
  geom_smooth(method = "lm", se = FALSE, color = "#B8323B", linewidth = 0.65) +
  scale_fill_gradient(low = "#E8EEF3", high = "#234E52", name = "Genes") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(
    title = "External PaxDB Protein Abundance Tracks the Protein Output Score",
    subtitle = "Gene-level medians across samples; PaxDB whole-organism integrated abundance",
    x = "log10(PaxDB protein abundance + 1)",
    y = "Metric value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
save_plot(p_pax, "paxdb_metric_scatter_panels", 15, 5.5)

message("Plotting PaxDB quartile medians")
quart_dt <- gene_level[is.finite(paxdb_abundance) & is.finite(protein_output_score)]
q <- quantile(quart_dt$paxdb_abundance, probs = seq(0, 1, 0.25), na.rm = TRUE)
quart_dt[, paxdb_quartile := cut(paxdb_abundance, breaks = unique(q), include.lowest = TRUE, labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"))]
quart_summary <- quart_dt[, .(
  protein_output_score = median(protein_output_score, na.rm = TRUE),
  ribosome_efficiency_score = median(ribosome_efficiency_score, na.rm = TRUE),
  rs_core_cpm = median(rs_core_cpm, na.rm = TRUE),
  baseline_cpm_line = median(baseline_cpm_line, na.rm = TRUE)
), by = paxdb_quartile]
quart_long <- melt(quart_summary, id.vars = "paxdb_quartile", variable.name = "metric", value.name = "median_value")
quart_long[, metric := factor(
  metric,
  levels = c("protein_output_score", "ribosome_efficiency_score", "rs_core_cpm", "baseline_cpm_line"),
  labels = c("Protein output score", "Ribosome efficiency score", "RS core CPM", "RNA baseline CPM")
)]

p_quart <- ggplot(quart_long, aes(x = paxdb_quartile, y = median_value, group = metric, color = metric)) +
  geom_point(size = 2.6) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c(
    "Protein output score" = "#2C7A7B",
    "Ribosome efficiency score" = "#B8323B",
    "RS core CPM" = "#2B6CB0",
    "RNA baseline CPM" = "#805AD5"
  )) +
  labs(
    title = "Protein Output Score Increases Across PaxDB Protein-Abundance Quartiles",
    subtitle = "Gene-level medians across samples",
    x = "PaxDB abundance quartile",
    y = "Median metric value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank()
  )
save_plot(p_quart, "paxdb_quartile_metric_medians", 13, 5)

message("Plotting focused protein output score correlation")
translation_rate_dt <- gene_level[
  is.finite(paxdb_abundance) &
    is.finite(protein_output_score) &
    paxdb_abundance >= 0
]
translation_rate_dt[, log10_paxdb := log10(paxdb_abundance + 1)]
spearman_rho <- suppressWarnings(cor(
  translation_rate_dt$log10_paxdb,
  translation_rate_dt$protein_output_score,
  method = "spearman"
))
spearman_p <- suppressWarnings(cor.test(
  translation_rate_dt$log10_paxdb,
  translation_rate_dt$protein_output_score,
  method = "spearman",
  exact = FALSE
)$p.value)
ann <- sprintf(
  "Spearman rho = %.3f\np %s\nn = %s genes",
  spearman_rho,
  ifelse(spearman_p < 2.2e-16, "< 2.2e-16", paste0("= ", signif(spearman_p, 2))),
  format(nrow(translation_rate_dt), big.mark = ",")
)

p_translation_rate <- ggplot(
  translation_rate_dt,
  aes(x = log10_paxdb, y = protein_output_score)
) +
  geom_point(alpha = 0.22, size = 0.75, color = "#2C7A7B") +
  geom_smooth(method = "lm", se = TRUE, color = "#B8323B", fill = "#F1C9CE", linewidth = 0.8) +
  annotate(
    "label",
    x = quantile(translation_rate_dt$log10_paxdb, 0.03, na.rm = TRUE),
    y = quantile(translation_rate_dt$protein_output_score, 0.97, na.rm = TRUE),
    label = ann,
    hjust = 0,
    vjust = 1,
    size = 4.2,
    label.size = 0.25,
    fill = "white"
  ) +
  labs(
    title = "Protein Output Score Correlates with PaxDB Protein Abundance",
    subtitle = "Gene-level medians across samples",
    x = "log10(PaxDB protein abundance + 1)",
    y = "Protein output score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )
save_plot(p_translation_rate, "protein_output_score_paxdb_correlation_dotplot", 7.5, 5.7)

message("Plotting spread protein output score rank correlation")
translation_rate_dt[, paxdb_percentile := frank(paxdb_abundance, ties.method = "average") / .N]
spearman_rank_rho <- suppressWarnings(cor(
  translation_rate_dt$paxdb_percentile,
  translation_rate_dt$protein_output_score,
  method = "spearman"
))
spearman_rank_p <- suppressWarnings(cor.test(
  translation_rate_dt$paxdb_percentile,
  translation_rate_dt$protein_output_score,
  method = "spearman",
  exact = FALSE
)$p.value)
ann_rank <- sprintf(
  "Spearman rho = %.3f\np %s\nn = %s genes",
  spearman_rank_rho,
  ifelse(spearman_rank_p < 2.2e-16, "< 2.2e-16", paste0("= ", signif(spearman_rank_p, 2))),
  format(nrow(translation_rate_dt), big.mark = ",")
)

p_translation_rate_rank <- ggplot(
  translation_rate_dt,
  aes(x = paxdb_percentile, y = protein_output_score)
) +
  geom_point(alpha = 0.18, size = 0.65, color = "#2C7A7B", position = position_jitter(width = 0.0025, height = 0.03, seed = 7)) +
  geom_smooth(method = "loess", se = TRUE, color = "#B8323B", fill = "#F1C9CE", linewidth = 0.85, span = 0.65) +
  annotate(
    "label",
    x = 0.04,
    y = quantile(translation_rate_dt$protein_output_score, 0.97, na.rm = TRUE),
    label = ann_rank,
    hjust = 0,
    vjust = 1,
    size = 4.2,
    fill = "white"
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Protein Output Score vs PaxDB Abundance Rank",
    subtitle = "Gene-level medians; percentile x-axis spreads low-abundance genes",
    x = "PaxDB protein abundance percentile",
    y = "Protein output score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )
save_plot(p_translation_rate_rank, "protein_output_score_paxdb_percentile_dotplot", 8.4, 5.7)

message("Plotting HRT protein output score dotplot")
hrt_dot_dt <- gene_level[is.finite(protein_output_score)]
hrt_dot_dt[, hrt_plot_group := fifelse(
  housekeeping_group == "Translation/ribosomal HK",
  "Translation/ribosomal HK",
  fifelse(housekeeping_binary == "HRT housekeeping", "HRT housekeeping", "Background")
)]
hrt_dot_dt[, hrt_plot_group := factor(
  hrt_plot_group,
  levels = c("Background", "HRT housekeeping", "Translation/ribosomal HK")
)]
hrt_counts <- hrt_dot_dt[, .N, by = hrt_plot_group]
hrt_labels <- setNames(
  paste0(as.character(hrt_counts$hrt_plot_group), "\n", "n = ", format(hrt_counts$N, big.mark = ",")),
  as.character(hrt_counts$hrt_plot_group)
)
sig_label <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  paste0("p = ", formatC(p, format = "e", digits = 2))
}
wilcox_hk <- suppressWarnings(wilcox.test(
  hrt_dot_dt[hrt_plot_group == "HRT housekeeping", protein_output_score],
  hrt_dot_dt[hrt_plot_group == "Background", protein_output_score],
  alternative = "greater",
  exact = FALSE
)$p.value)
wilcox_translation_hk <- suppressWarnings(wilcox.test(
  hrt_dot_dt[hrt_plot_group == "Translation/ribosomal HK", protein_output_score],
  hrt_dot_dt[hrt_plot_group == "Background", protein_output_score],
  alternative = "greater",
  exact = FALSE
)$p.value)
sig_dt <- data.table(
  x = c(1, 1),
  xend = c(2, 3),
  y = c(
    quantile(hrt_dot_dt$protein_output_score, 0.965, na.rm = TRUE),
    quantile(hrt_dot_dt$protein_output_score, 0.995, na.rm = TRUE)
  ),
  label = c(
    paste0("Background vs HRT housekeeping\n", sig_label(wilcox_hk)),
    paste0("Background vs translation/ribosomal HK\n", sig_label(wilcox_translation_hk))
  )
)
sig_dt[, y_text := y + 0.45]

p_hrt_dot <- ggplot(hrt_dot_dt, aes(x = hrt_plot_group, y = protein_output_score, color = hrt_plot_group)) +
  geom_point(alpha = 0.18, size = 0.65, position = position_jitter(width = 0.22, height = 0.025, seed = 11)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.48, color = "black", linewidth = 0.45) +
  geom_segment(
    data = sig_dt,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_segment(
    data = sig_dt,
    aes(x = x, xend = x, y = y - 0.18, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_segment(
    data = sig_dt,
    aes(x = xend, xend = xend, y = y - 0.18, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_label(
    data = sig_dt,
    aes(x = (x + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    label.size = 0.2,
    fill = "white"
  ) +
  scale_x_discrete(labels = hrt_labels) +
  scale_color_manual(values = c(
    "Background" = "#7A838D",
    "HRT housekeeping" = "#2C7A7B",
    "Translation/ribosomal HK" = "#805AD5"
  )) +
  coord_cartesian(ylim = c(NA, max(sig_dt$y_text, na.rm = TRUE) + 0.75)) +
  labs(
    title = "HRT Housekeeping Genes Have Higher Protein Output Score",
    subtitle = "Gene-level medians across samples; Wilcoxon tests are one-sided vs background",
    x = NULL,
    y = "Protein output score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )
save_plot(p_hrt_dot, "protein_output_score_hrt_housekeeping_dotplot", 8.4, 6.2)

message("Plotting HRT ribosome efficiency score dotplot")
hrt_te_dot_dt <- gene_level[is.finite(ribosome_efficiency_score)]
hrt_te_dot_dt[, hrt_plot_group := fifelse(
  housekeeping_group == "Translation/ribosomal HK",
  "Translation/ribosomal HK",
  fifelse(housekeeping_binary == "HRT housekeeping", "HRT housekeeping", "Background")
)]
hrt_te_dot_dt[, hrt_plot_group := factor(
  hrt_plot_group,
  levels = c("Background", "HRT housekeeping", "Translation/ribosomal HK")
)]
hrt_te_counts <- hrt_te_dot_dt[, .N, by = hrt_plot_group]
hrt_te_labels <- setNames(
  paste0(as.character(hrt_te_counts$hrt_plot_group), "\n", "n = ", format(hrt_te_counts$N, big.mark = ",")),
  as.character(hrt_te_counts$hrt_plot_group)
)
wilcox_te_hk <- suppressWarnings(wilcox.test(
  hrt_te_dot_dt[hrt_plot_group == "HRT housekeeping", ribosome_efficiency_score],
  hrt_te_dot_dt[hrt_plot_group == "Background", ribosome_efficiency_score],
  alternative = "greater",
  exact = FALSE
)$p.value)
wilcox_te_translation_hk <- suppressWarnings(wilcox.test(
  hrt_te_dot_dt[hrt_plot_group == "Translation/ribosomal HK", ribosome_efficiency_score],
  hrt_te_dot_dt[hrt_plot_group == "Background", ribosome_efficiency_score],
  alternative = "greater",
  exact = FALSE
)$p.value)
sig_te_dt <- data.table(
  x = c(1, 1),
  xend = c(2, 3),
  y = c(
    quantile(hrt_te_dot_dt$ribosome_efficiency_score, 0.965, na.rm = TRUE),
    quantile(hrt_te_dot_dt$ribosome_efficiency_score, 0.995, na.rm = TRUE)
  ),
  label = c(
    paste0("Background vs HRT housekeeping\n", sig_label(wilcox_te_hk)),
    paste0("Background vs translation/ribosomal HK\n", sig_label(wilcox_te_translation_hk))
  )
)
sig_te_dt[, y_text := y + 0.45]

p_hrt_te_dot <- ggplot(hrt_te_dot_dt, aes(x = hrt_plot_group, y = ribosome_efficiency_score, color = hrt_plot_group)) +
  geom_point(alpha = 0.18, size = 0.65, position = position_jitter(width = 0.22, height = 0.025, seed = 11)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.48, color = "black", linewidth = 0.45) +
  geom_segment(
    data = sig_te_dt,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_segment(
    data = sig_te_dt,
    aes(x = x, xend = x, y = y - 0.18, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_segment(
    data = sig_te_dt,
    aes(x = xend, xend = xend, y = y - 0.18, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_label(
    data = sig_te_dt,
    aes(x = (x + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    label.size = 0.2,
    fill = "white"
  ) +
  scale_x_discrete(labels = hrt_te_labels) +
  scale_color_manual(values = c(
    "Background" = "#7A838D",
    "HRT housekeeping" = "#2C7A7B",
    "Translation/ribosomal HK" = "#805AD5"
  )) +
  coord_cartesian(ylim = c(NA, max(sig_te_dt$y_text, na.rm = TRUE) + 0.75)) +
  labs(
    title = "HRT Housekeeping Genes Have Higher Ribosome Efficiency Score",
    subtitle = "Gene-level medians across samples; Wilcoxon tests are one-sided vs background",
    x = NULL,
    y = "Ribosome efficiency score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )
save_plot(p_hrt_te_dot, "ribosome_efficiency_score_hrt_housekeeping_dotplot", 8.4, 6.2)

message("Wrote plots to: ", OUT_DIR)
