# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base <- analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05", "Results")
out_dir <- analysis_path("Limma_translation_metrics_lfc0.7_rawP0.05", "Metric_relationships")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_cut <- 0.05
lfc_cut <- 0.7
contrasts <- c("Resistance_baseline", "VCR_sensitive", "VCR_resistant", "Interaction")

read_metric <- function(metric, contrast) {
  d <- fread(file.path(base, metric, paste0(contrast, "_limma_all_genes.csv")))
  d <- d[!is.na(gene_name) & gene_name != ""]
  d[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  d[, dir := fifelse(
    P.Value < p_cut & logFC >= lfc_cut, "Up",
    fifelse(P.Value < p_cut & logFC <= -lfc_cut, "Down", "NS")
  )]
  d[, .(gene_id_clean, gene_name, logFC, P.Value, adj.P.Val, dir)]
}

all <- rbindlist(lapply(contrasts, function(ct) {
  s <- read_metric("scanning_score", ct)
  setnames(s, c("logFC", "P.Value", "adj.P.Val", "dir"), c("scanning_logFC", "scanning_P", "scanning_FDR", "scanning_dir"))
  r <- read_metric("ribosome_efficiency_score", ct)
  setnames(r, c("logFC", "P.Value", "adj.P.Val", "dir"), c("ribosome_logFC", "ribosome_P", "ribosome_FDR", "ribosome_dir"))
  w <- merge(s, r, by = c("gene_id_clean", "gene_name"))
  w[, contrast := ct]
  w
}), fill = TRUE)

cor_stats <- all[
  ,
  {
    ok <- is.finite(scanning_logFC) & is.finite(ribosome_logFC)
    sp <- cor.test(scanning_logFC[ok], ribosome_logFC[ok], method = "spearman", exact = FALSE)
    pe <- cor.test(scanning_logFC[ok], ribosome_logFC[ok], method = "pearson")
    .(
      n = sum(ok),
      spearman_rho = unname(sp$estimate),
      spearman_p = sp$p.value,
      pearson_r = unname(pe$estimate),
      pearson_p = pe$p.value
    )
  },
  by = contrast
]

direction_counts <- all[, .N, by = .(contrast, scanning_dir, ribosome_dir)][order(contrast, scanning_dir, ribosome_dir)]
scanning_up_ribosome <- all[scanning_dir == "Up", .N, by = .(contrast, ribosome_dir)][order(contrast, ribosome_dir)]
scanning_down_ribosome <- all[scanning_dir == "Down", .N, by = .(contrast, ribosome_dir)][order(contrast, ribosome_dir)]

fwrite(all, file.path(out_dir, "scanning_vs_ribosome_engagement_all_contrasts.csv"))
fwrite(cor_stats, file.path(out_dir, "scanning_vs_ribosome_engagement_correlations.csv"))
fwrite(direction_counts, file.path(out_dir, "scanning_vs_ribosome_engagement_direction_counts.csv"))

all[, highlight := fifelse(gene_name %in% c("TRA2A", "LPXN", "HNRNPD", "MAPKBP1", "SEC24C"), gene_name, "Other")]

p <- ggplot(all, aes(scanning_logFC, ribosome_logFC)) +
  geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", color = "grey50") +
  geom_point(aes(color = highlight != "Other"), alpha = 0.45, size = 1.1) +
  geom_text(
    data = all[highlight != "Other"],
    aes(label = gene_name),
    size = 2.7,
    vjust = -0.7,
    check_overlap = TRUE
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.55) +
  facet_wrap(~ contrast, scales = "free", ncol = 2) +
  scale_color_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "#9CA3AF"), guide = "none") +
  labs(
    title = "Scanning score versus ribosome engagement",
    subtitle = "limma logFC across genes; dotted lines mark |logFC| = 0.7",
    x = "Scanning score logFC",
    y = "Ribosome engagement logFC"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "scanning_vs_ribosome_engagement_scatter.png"), p, width = 9, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "scanning_vs_ribosome_engagement_scatter.pdf"), p, width = 9, height = 7, bg = "white")

cat("\nCorrelation stats:\n")
print(cor_stats)
cat("\nDirection overlap counts:\n")
print(direction_counts)
cat("\nScanning Up genes by ribosome engagement direction:\n")
print(scanning_up_ribosome)
cat("\nScanning Down genes by ribosome engagement direction:\n")
print(scanning_down_ribosome)
cat("\nOutput dir:\n")
cat(out_dir, "\n")
