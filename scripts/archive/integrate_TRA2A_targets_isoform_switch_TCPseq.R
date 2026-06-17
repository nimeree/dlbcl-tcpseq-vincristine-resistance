# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)

base_dir <- analysis_path()
longread_dir <- file.path(base_dir, "LongRead_TRA2A")
step2_file <- file.path(longread_dir, "Step2_TRA2A_targets_isoform_usage", "TRA2A_target_genes_SUDHL8_isoform_usage_tests.csv")
step3_file <- file.path(longread_dir, "Step3_GenomeWide_isoform_usage_screen", "SUDHL8_genomewide_isoform_usage_screen.csv")
out_dir <- file.path(longread_dir, "Step4_IsoformSwitch_TCPseq_integration")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_fraction <- function(frac, contrast) {
  f <- file.path(base_dir, paste0("Fraction_", frac), paste0(contrast, "_", frac, "_results_all.csv"))
  d <- fread(f)
  d[, fraction := frac]
  d[, contrast := contrast]
  d[, tcp_log2FC := log2FoldChange]
  d[, tcp_pvalue := pvalue]
  d[, tcp_padj := padj]
  d[, tcp_sig_lfc0.7_fdr0.05 := padj < 0.05 & abs(log2FoldChange) >= 0.7]
  d[, .(gene_id, gene_name, fraction, contrast, tcp_log2FC, tcp_pvalue, tcp_padj, direction, tcp_sig_lfc0.7_fdr0.05, baseMean)]
}

contrasts <- c("Vin_Resistant_vs_Sensitive", "Resistant_Vin_vs_DMSO", "Sensitive_Vin_vs_DMSO")
tcp <- rbindlist(lapply(c("SSU", "RS", "DS"), function(fr) {
  rbindlist(lapply(contrasts, function(co) read_fraction(fr, co)))
}))

step2 <- fread(step2_file)
step3 <- fread(step3_file)

iso <- unique(rbindlist(list(
  step3[, .(
    gene_id, gene_name,
    source = "genomewide",
    switched_stringent,
    switched_relaxed,
    TRA2A_eCLIP_target,
    max_abs_proportion_shift,
    max_shift_isoform,
    max_shift_transcript_id,
    max_shift_original_prop,
    max_shift_resistant_prop,
    max_shift_resistant_minus_original,
    isoform_padj = padj,
    dominant_isoform_switch
  )],
  step2[, .(
    gene_id, gene_name,
    source = "TRA2A_eCLIP_target_screen",
    switched_stringent = padj < 0.05 & max_abs_proportion_shift >= 0.20,
    switched_relaxed = padj < 0.05 & max_abs_proportion_shift >= 0.10,
    TRA2A_eCLIP_target = TRUE,
    max_abs_proportion_shift,
    max_shift_isoform,
    max_shift_transcript_id,
    max_shift_original_prop,
    max_shift_resistant_prop,
    max_shift_resistant_minus_original,
    isoform_padj = padj,
    dominant_isoform_switch
  )]
), fill = TRUE), by = c("gene_id", "source"))

merged <- merge(iso, tcp, by = c("gene_id", "gene_name"), all.x = TRUE, allow.cartesian = TRUE)
merged[, tcp_raw_sig_lfc0.7 := tcp_pvalue < 0.05 & abs(tcp_log2FC) >= 0.7]
merged[, tcp_any_fdr_sig := tcp_padj < 0.05]
merged[, is_ribo_gene := grepl("^RP[SL]", gene_name)]

vin <- merged[contrast == "Vin_Resistant_vs_Sensitive"]
vin_sig <- vin[tcp_raw_sig_lfc0.7 == TRUE | tcp_any_fdr_sig == TRUE]
vin_stringent <- vin[switched_stringent == TRUE]

hsp90ab1 <- merged[gene_name == "HSP90AB1"]
rpl_rps <- merged[is_ribo_gene == TRUE & switched_stringent == TRUE]
rpl_rps_vin <- rpl_rps[contrast == "Vin_Resistant_vs_Sensitive"]
rpl_rps_vin_hits <- rpl_rps_vin[tcp_raw_sig_lfc0.7 == TRUE | tcp_any_fdr_sig == TRUE][
  order(gene_name, fraction)
]

summary_by_fraction <- vin_stringent[, .(
  switched_genes_with_tcp = uniqueN(gene_id),
  fdr_sig = uniqueN(gene_id[tcp_any_fdr_sig == TRUE]),
  rawP_lfc_sig = uniqueN(gene_id[tcp_raw_sig_lfc0.7 == TRUE]),
  up_rawP_lfc = uniqueN(gene_id[tcp_raw_sig_lfc0.7 == TRUE & tcp_log2FC > 0]),
  down_rawP_lfc = uniqueN(gene_id[tcp_raw_sig_lfc0.7 == TRUE & tcp_log2FC < 0])
), by = fraction][order(fraction)]

ribo_summary <- rpl_rps_vin[, .(
  rpl_rps_switched_genes = uniqueN(gene_id),
  fdr_sig = uniqueN(gene_id[tcp_any_fdr_sig == TRUE]),
  rawP_lfc_sig = uniqueN(gene_id[tcp_raw_sig_lfc0.7 == TRUE]),
  up_rawP_lfc = uniqueN(gene_id[tcp_raw_sig_lfc0.7 == TRUE & tcp_log2FC > 0]),
  down_rawP_lfc = uniqueN(gene_id[tcp_raw_sig_lfc0.7 == TRUE & tcp_log2FC < 0])
), by = fraction][order(fraction)]

top_integrated <- vin_stringent[
  tcp_raw_sig_lfc0.7 == TRUE | tcp_any_fdr_sig == TRUE
][
  order(tcp_padj, -max_abs_proportion_shift)
][
  ,
  .(gene_id, gene_name, TRA2A_eCLIP_target, is_ribo_gene, fraction, tcp_log2FC, tcp_pvalue, tcp_padj,
    direction, tcp_sig_lfc0.7_fdr0.05, tcp_raw_sig_lfc0.7, max_abs_proportion_shift,
    max_shift_isoform, max_shift_original_prop, max_shift_resistant_prop,
    isoform_padj, dominant_isoform_switch)
]

write.csv(merged, file.path(out_dir, "all_isoform_switch_TCPseq_fraction_integration.csv"), row.names = FALSE)
write.csv(hsp90ab1, file.path(out_dir, "HSP90AB1_isoform_switch_TCPseq_integration.csv"), row.names = FALSE)
write.csv(rpl_rps_vin, file.path(out_dir, "RPL_RPS_stringent_switches_VinResistantVsSensitive_TCPseq.csv"), row.names = FALSE)
write.csv(rpl_rps_vin_hits, file.path(out_dir, "RPL_RPS_stringent_switches_with_TCPseq_hits.csv"), row.names = FALSE)
write.csv(summary_by_fraction, file.path(out_dir, "stringent_isoform_switch_TCPseq_overlap_summary_by_fraction.csv"), row.names = FALSE)
write.csv(ribo_summary, file.path(out_dir, "RPL_RPS_isoform_switch_TCPseq_overlap_summary_by_fraction.csv"), row.names = FALSE)
write.csv(top_integrated, file.path(out_dir, "top_integrated_isoform_switch_TCPseq_hits.csv"), row.names = FALSE)

cat("\nHSP90AB1 integration:\n")
print(hsp90ab1[order(contrast, fraction)])

cat("\nStringent isoform switches with TCP-seq Vin Resistant vs Sensitive signal by fraction:\n")
print(summary_by_fraction)

cat("\nRPL/RPS stringent isoform switches with TCP-seq Vin Resistant vs Sensitive signal by fraction:\n")
print(ribo_summary)

cat("\nRPL/RPS Vin Resistant vs Sensitive TCP-seq hits:\n")
print(rpl_rps_vin_hits[, .(
  gene_name, fraction, tcp_log2FC, tcp_pvalue, tcp_padj, direction,
  tcp_sig_lfc0.7_fdr0.05, tcp_raw_sig_lfc0.7,
  max_abs_proportion_shift, max_shift_isoform,
  max_shift_original_prop, max_shift_resistant_prop,
  isoform_padj
)][order(fraction, tcp_padj)], nrows = 100)

cat("\nTop integrated hits:\n")
print(head(top_integrated, 50))

cat("\nOutput directory:\n")
cat(out_dir, "\n")
