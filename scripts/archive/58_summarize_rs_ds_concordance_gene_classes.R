# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- analysis_path()
psite_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
out_dir <- file.path(psite_dir, "RS_DS_similarity_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

overlap_file <- file.path(out_dir, "RS_DS_gene_level_sig_overlap_details.csv")
if (!file.exists(overlap_file)) stop("Missing overlap file: ", overlap_file)

contrast_files <- c(
  file.path(base_dir, "Fraction_DS", "Sensitive_Vin_vs_DMSO_DS_results_all.csv"),
  file.path(base_dir, "Fraction_DS", "Resistant_Vin_vs_DMSO_DS_results_all.csv"),
  file.path(base_dir, "Fraction_DS", "Vin_Resistant_vs_Sensitive_DS_results_all.csv"),
  file.path(base_dir, "Fraction_RS", "Sensitive_Vin_vs_DMSO_RS_results_all.csv"),
  file.path(base_dir, "Fraction_RS", "Resistant_Vin_vs_DMSO_RS_results_all.csv"),
  file.path(base_dir, "Fraction_RS", "Vin_Resistant_vs_Sensitive_RS_results_all.csv")
)

anno <- rbindlist(lapply(contrast_files[file.exists(contrast_files)], function(f) {
  d <- fread(f, select = c("gene_id", "gene_name", "gene_function", "gene_type"))
  d[, gene_key := sub("\\.\\d+$", "", gene_id)]
  d[, .(gene_key, gene_name, gene_function, gene_type)]
}), fill = TRUE)
anno <- unique(anno[!is.na(gene_key) & nzchar(gene_key)])
anno <- anno[, .(
  gene_name_anno = first(na.omit(gene_name)),
  gene_function = first(na.omit(gene_function)),
  gene_type = first(na.omit(gene_type))
), by = gene_key]

classify_gene <- function(gene_name, gene_function) {
  nm <- toupper(ifelse(is.na(gene_name), "", gene_name))
  fn <- tolower(ifelse(is.na(gene_function), "", gene_function))
  txt <- paste(nm, fn)

  fifelse(grepl("^H(1|2A|2B|3|4)[A-Z0-9]+$", nm), "Histone / chromatin core",
  fifelse(grepl("^(RPL|RPS)[0-9A-Z]|ribosomal protein|ribosome|translation initiation|translation elongation|^EIF|^EEF", txt), "Ribosome / translation",
  fifelse(grepl("splice|splicing|mrna processing|rna binding|ribonucleoprotein|snrnp|^SRSF|^HNRNP|^RBM|^DDX|^DHX|^SF3|^U2AF", txt), "RNA processing / splicing",
  fifelse(grepl("dna repair|dna replication|replication|cell cycle|mitotic|mitosis|chromatin|mcm|pcna|rfc|topoisomerase|recombination", txt), "DNA replication / repair / cell cycle",
  fifelse(grepl("^ZNF|zinc finger|transcription factor|transcription regulator|rna polymerase|mediator complex", txt), "Transcription / zinc-finger regulation",
  fifelse(grepl("mitochond|^MRPL|^MRPS|^NDUF|^COX|^ATP5|^UQCR", txt), "Mitochondrial / respiration",
  fifelse(grepl("proteasome|ubiquitin|chaperone|heat shock|^PSM|^HSP|protein folding", txt), "Proteostasis / ubiquitin / chaperone",
  fifelse(grepl("vesicle|transport|traffick|membrane|golgi|endosom|lysosom|receptor|solute carrier|^SLC", txt), "Membrane / trafficking / transport",
  fifelse(grepl("kinase|phosphatase|apoptosis|mapk|caspase|stress response|signaling", txt), "Signaling / stress / apoptosis",
  fifelse(grepl("metabolic|metabolism|dehydrogenase|synthase|synthetase|transferase|lipid|glycol|oxidase", txt), "Metabolism / enzymes",
          "Other / mixed"))))))))))
}

overlap <- fread(overlap_file)
overlap <- merge(overlap, anno, by = "gene_key", all.x = TRUE)
overlap[, gene_name_final := fifelse(!is.na(gene_name) & nzchar(gene_name), gene_name, gene_name_anno)]
overlap[, gene_class := classify_gene(gene_name_final, gene_function)]
overlap[, same_direction_sign := fifelse(
  overlap_class == "both_same_direction" & RS_direction == "Up" & DS_direction == "Up", "Both up",
  fifelse(overlap_class == "both_same_direction" & RS_direction == "Down" & DS_direction == "Down", "Both down", NA_character_)
)]

class_summary <- overlap[, .(
  n_genes = .N,
  example_genes = paste(head(unique(gene_name_final[!is.na(gene_name_final) & nzchar(gene_name_final)]), 12), collapse = "; ")
), by = .(contrast, overlap_class, gene_class)]
setorder(class_summary, contrast, overlap_class, -n_genes)

same_summary <- overlap[overlap_class == "both_same_direction", .(
  n_genes = .N,
  example_genes = paste(head(unique(gene_name_final[!is.na(gene_name_final) & nzchar(gene_name_final)]), 15), collapse = "; ")
), by = .(contrast, same_direction_sign, gene_class)]
setorder(same_summary, contrast, same_direction_sign, -n_genes)

overall_same <- overlap[overlap_class == "both_same_direction", .(
  n_genes = .N,
  n_contrasts = uniqueN(contrast),
  example_genes = paste(head(unique(gene_name_final[!is.na(gene_name_final) & nzchar(gene_name_final)]), 18), collapse = "; ")
), by = gene_class]
setorder(overall_same, -n_genes)

fwrite(overlap, file.path(out_dir, "RS_DS_gene_level_sig_overlap_details_with_gene_classes.csv"))
fwrite(class_summary, file.path(out_dir, "RS_DS_concordance_gene_class_summary_by_contrast.csv"))
fwrite(same_summary, file.path(out_dir, "RS_DS_same_direction_gene_class_summary_by_contrast_direction.csv"))
fwrite(overall_same, file.path(out_dir, "RS_DS_same_direction_gene_class_summary_overall.csv"))

cat("\nMajor classes among both-significant same-direction RS/DS genes:\n")
print(overall_same)

cat("\nSame-direction classes split by contrast and up/down:\n")
print(same_summary)

cat("\nSaved class summaries to:\n", out_dir, "\n", sep = "")
