# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)
library(Biostrings)
library(ggplot2)

base_out <- analysis_path("LongRead_TRA2A")
step3_dir <- file.path(base_out, "Step3_GenomeWide_isoform_usage_screen")
out_dir <- file.path(base_out, "Step3B_TRA2A_CDS_motif_enrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rdata <- input_path("SUDHL.RData")
cds_fasta <- input_path("Homo_sapiens.GRCh38.cds.all.fa.gz")
screen_file <- file.path(step3_dir, "SUDHL8_genomewide_isoform_usage_screen.csv")

motifs <- c(AGAA = "AGAA", GAAGAA = "GAAGAA", GAA = "GAA")

message("Loading Step 3 genome-wide screen...")
screen <- fread(screen_file)
eligible_genes <- screen$gene_id

message("Loading expressed transcript IDs...")
e <- new.env()
load(rdata, envir = e)
expressed_tx <- rownames(e$count_matrix)
expressed_tx_versionless <- sub("[.][0-9]+$", "", expressed_tx)

message("Reading CDS FASTA...")
cds <- readDNAStringSet(cds_fasta)
headers <- names(cds)
tx_id <- sub(" .*", "", headers)
tx_versionless <- sub("[.][0-9]+$", "", tx_id)
gene_id <- sub(".* gene:([^ .]+)([.][0-9]+)?.*", "\\1", headers)
gene_symbol <- sub(".* gene_symbol:([^ ]+).*", "\\1", headers)
gene_symbol[gene_symbol == headers] <- NA_character_

meta <- data.table(
  seq_index = seq_along(cds),
  tx_id = tx_id,
  tx_versionless = tx_versionless,
  gene_id = gene_id,
  gene_symbol = gene_symbol,
  width = width(cds)
)
meta <- meta[
  gene_id %in% eligible_genes &
    tx_versionless %in% expressed_tx_versionless &
    width >= 30
]

message("Counting motifs in ", nrow(meta), " expressed CDS sequences...")
count_one <- function(pattern, seqs) {
  vcountPattern(pattern, seqs, fixed = TRUE)
}

seq_subset <- cds[meta$seq_index]
for (nm in names(motifs)) {
  meta[, paste0("count_", nm) := count_one(motifs[[nm]], seq_subset)]
}

gene_motif <- meta[
  ,
  c(
    .(
      gene_symbol = gene_symbol[which.max(!is.na(gene_symbol))],
      n_cds_transcripts = .N,
      total_cds_bp = sum(width)
    ),
    lapply(.SD, sum)
  ),
  by = gene_id,
  .SDcols = patterns("^count_")
]
for (nm in names(motifs)) {
  count_col <- paste0("count_", nm)
  dens_col <- paste0("density_per_kb_", nm)
  present_col <- paste0("present_", nm)
  gene_motif[, (dens_col) := get(count_col) / total_cds_bp * 1000]
  gene_motif[, (present_col) := get(count_col) > 0]
}

gene_motif <- merge(gene_motif, screen, by = "gene_id", all.x = TRUE)

test_motif <- function(motif_name, switch_col) {
  dens_col <- paste0("density_per_kb_", motif_name)
  present_col <- paste0("present_", motif_name)
  top_quartile <- gene_motif[[dens_col]] >= quantile(gene_motif[[dens_col]], 0.75, na.rm = TRUE)

  wt <- wilcox.test(
    gene_motif[get(switch_col) == TRUE, get(dens_col)],
    gene_motif[get(switch_col) == FALSE, get(dens_col)],
    alternative = "greater"
  )
  ft_present <- fisher.test(table(switched = gene_motif[[switch_col]], motif_present = gene_motif[[present_col]]))
  ft_top <- fisher.test(table(switched = gene_motif[[switch_col]], motif_top_quartile = top_quartile))

  data.table(
    motif = motif_name,
    switch_definition = switch_col,
    switched_genes_with_cds = sum(gene_motif[[switch_col]], na.rm = TRUE),
    background_genes_with_cds = sum(!gene_motif[[switch_col]], na.rm = TRUE),
    switched_mean_density_per_kb = mean(gene_motif[get(switch_col) == TRUE, get(dens_col)], na.rm = TRUE),
    background_mean_density_per_kb = mean(gene_motif[get(switch_col) == FALSE, get(dens_col)], na.rm = TRUE),
    switched_median_density_per_kb = median(gene_motif[get(switch_col) == TRUE, get(dens_col)], na.rm = TRUE),
    background_median_density_per_kb = median(gene_motif[get(switch_col) == FALSE, get(dens_col)], na.rm = TRUE),
    wilcox_greater_pvalue = wt$p.value,
    motif_present_odds_ratio = unname(ft_present$estimate),
    motif_present_fisher_pvalue = ft_present$p.value,
    top_quartile_odds_ratio = unname(ft_top$estimate),
    top_quartile_fisher_pvalue = ft_top$p.value
  )
}

enrichment <- rbindlist(lapply(names(motifs), function(m) {
  rbindlist(list(test_motif(m, "switched_stringent"), test_motif(m, "switched_relaxed")))
}))
enrichment[, wilcox_padj := p.adjust(wilcox_greater_pvalue, method = "BH")]
enrichment[, top_quartile_fisher_padj := p.adjust(top_quartile_fisher_pvalue, method = "BH")]
enrichment[, motif_present_fisher_padj := p.adjust(motif_present_fisher_pvalue, method = "BH")]

message("Writing outputs...")
write.csv(meta, file.path(out_dir, "expressed_CDS_transcripts_motif_counts.csv"), row.names = FALSE)
write.csv(gene_motif, file.path(out_dir, "gene_level_TRA2A_CDS_motif_density_with_switch_status.csv"), row.names = FALSE)
write.csv(enrichment, file.path(out_dir, "TRA2A_CDS_motif_enrichment_switched_vs_background.csv"), row.names = FALSE)

plot_dt <- melt(
  gene_motif,
  id.vars = c("gene_id", "gene_name", "switched_stringent", "switched_relaxed", "TRA2A_eCLIP_target"),
  measure.vars = paste0("density_per_kb_", names(motifs)),
  variable.name = "motif",
  value.name = "density_per_kb"
)
plot_dt[, motif := sub("^density_per_kb_", "", motif)]
plot_dt[, stringent_class := fifelse(switched_stringent, "Stringent switched", "Background")]

p1 <- ggplot(plot_dt, aes(x = stringent_class, y = density_per_kb, fill = stringent_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.65) +
  geom_jitter(width = 0.18, alpha = 0.15, size = 0.5) +
  facet_wrap(~ motif, scales = "free_y") +
  scale_fill_manual(values = c("Background" = "grey70", "Stringent switched" = "#7B3294")) +
  labs(
    title = "TRA2A-like motif density in coding isoforms",
    subtitle = "SUDHL8 stringent isoform-switch genes vs expressed multi-isoform background",
    x = NULL,
    y = "Motif occurrences per kb CDS"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "TRA2A_CDS_motif_density_stringent_switches.png"), p1, width = 9, height = 5.5, dpi = 300)

readme <- c(
  "Step 3B-lite: TRA2A-like motif enrichment in expressed coding isoforms.",
  "This uses the locally available Ensembl GRCh38 CDS FASTA, not whole-genome intronic flanking sequence.",
  "Motifs tested: AGAA, GAAGAA, GAA on transcript-sense CDS sequences.",
  "Foreground: SUDHL8 isoform-switch genes from Step 3.",
  "Background: all Step 3 eligible expressed multi-isoform genes with expressed CDS sequences.",
  "Tests: one-sided Wilcoxon test for higher motif density in switched genes; Fisher tests for motif presence and top-quartile motif density.",
  "Interpretation: supportive CDS-level motif evidence only. The ideal TRA2A mechanistic test still requires whole-genome FASTA and differentially included/excluded exon flanking intronic sequence."
)
writeLines(readme, file.path(out_dir, "README_step3b_CDS_motif_method_and_limitations.txt"))

cat("\nGenes with expressed CDS motif data:", nrow(gene_motif), "\n")
cat("Stringent switched with CDS:", gene_motif[switched_stringent == TRUE, .N], "\n")
cat("Relaxed switched with CDS:", gene_motif[switched_relaxed == TRUE, .N], "\n")
cat("\nMotif enrichment results:\n")
print(enrichment)
cat("\nOutput directory:\n")
cat(out_dir, "\n")
