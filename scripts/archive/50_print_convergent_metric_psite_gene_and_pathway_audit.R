# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

suppressPackageStartupMessages({
  library(data.table)
  library(gprofiler2)
})

base_dir <- analysis_path()
count_dir <- file.path(base_dir, "Psite_fraction_limma_lfc0.7_rawP0.05")
metric_dir <- file.path(base_dir, "Limma_translation_metrics_lfc0.7_rawP0.05", "Results")

p_cut <- 0.05
lfc_cut <- 0.7
min_hits <- 5
max_term_size <- 500

direction_call <- function(logfc, pval) {
  fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc >= lfc_cut, "Up",
          fifelse(!is.na(pval) & pval < p_cut & !is.na(logfc) & logfc <= -lfc_cut, "Down", "NS"))
}

pairs <- data.table(
  analysis = c("Scanning", "Ribosome engagement", "Collision"),
  metric_folder = c("scanning_score", "ribosome_efficiency_score", "collision_score"),
  fraction = c("SSU", "RS", "DS")
)

contrasts <- data.table(
  contrast = c("Sensitive Vin vs DMSO", "Resistant Vin vs DMSO", "Resistance baseline", "Interaction"),
  metric_contrast = c("VCR_sensitive", "VCR_resistant", "Resistance_baseline", "Interaction"),
  count_contrast = c("Sensitive_Vin_vs_DMSO", "Resistant_Vin_vs_DMSO", "Resistance_baseline", "Interaction")
)

read_metric <- function(folder, contrast) {
  f <- file.path(metric_dir, folder, paste0(contrast, "_limma_all_genes.csv"))
  d <- fread(f)
  if (!"gene_name" %in% names(d)) d[, gene_name := gene_id_clean]
  d[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  d[, .(
    gene_id_clean,
    gene_name_metric = gene_name,
    metric_logFC = logFC,
    metric_P.Value = P.Value,
    metric_adj.P.Val = adj.P.Val,
    metric_direction = direction_call(logFC, P.Value),
    metric_FDR_sig = !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) >= lfc_cut
  )]
}

read_count <- function(fraction, contrast) {
  f <- file.path(count_dir, paste0("Fraction_", fraction), paste0(contrast, "_", fraction, "_psite_limma_all_genes.csv"))
  d <- fread(f)
  if (!"gene_name" %in% names(d)) d[, gene_name := gene_id_clean]
  d[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  d[, .(
    gene_id_clean,
    gene_name_count = gene_name,
    count_logFC = logFC,
    count_P.Value = P.Value,
    count_adj.P.Val = adj.P.Val,
    count_direction = direction_call(logFC, P.Value),
    count_FDR_sig = !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) >= lfc_cut
  )]
}

all_pairs <- rbindlist(lapply(seq_len(nrow(contrasts)), function(ci) {
  cinfo <- contrasts[ci]
  rbindlist(lapply(seq_len(nrow(pairs)), function(pi) {
    pinfo <- pairs[pi]
    x <- merge(
      read_metric(pinfo$metric_folder, cinfo$metric_contrast),
      read_count(pinfo$fraction, cinfo$count_contrast),
      by = "gene_id_clean",
      all = FALSE
    )
    x[, `:=`(
      gene_name = fifelse(!is.na(gene_name_metric) & nzchar(gene_name_metric), gene_name_metric, gene_name_count),
      contrast = cinfo$contrast,
      analysis = pinfo$analysis,
      fraction = pinfo$fraction
    )]
    x
  }), use.names = TRUE)
}), use.names = TRUE)

all_pairs[, metric_sig := metric_direction %in% c("Up", "Down")]
all_pairs[, count_sig := count_direction %in% c("Up", "Down")]
all_pairs[, convergent := metric_sig & count_sig & metric_direction == count_direction]
all_pairs[, convergent_direction := fifelse(convergent, metric_direction, "Not convergent")]
all_pairs[, both_FDR_same_direction := convergent & metric_FDR_sig & count_FDR_sig]

cat("\n=== Convergent gene counts: metric significant + P-site fraction significant, same direction ===\n")
conv_summary <- all_pairs[, .(
  n_genes_merged = .N,
  metric_sig_n = sum(metric_sig),
  count_sig_n = sum(count_sig),
  convergent_n = sum(convergent),
  convergent_up = sum(convergent_direction == "Up"),
  convergent_down = sum(convergent_direction == "Down"),
  both_FDR_same_direction_n = sum(both_FDR_same_direction)
), by = .(contrast, analysis, fraction)]
setorder(conv_summary, contrast, analysis)
print(conv_summary)

cat("\n=== Convergent gene lists ===\n")
for (cn in contrasts$contrast) {
  cat("\n-- ", cn, " --\n", sep = "")
  for (an in pairs$analysis) {
    dt <- all_pairs[contrast == cn & analysis == an & convergent == TRUE]
    if (nrow(dt) == 0) {
      cat(an, ": 0\n", sep = "")
    } else {
      for (dirn in c("Up", "Down")) {
        genes <- dt[convergent_direction == dirn][order(metric_P.Value + count_P.Value), unique(gene_name)]
        cat(an, " ", dirn, " (", length(genes), "): ", paste(head(genes, 35), collapse = ", "),
            if (length(genes) > 35) " ..." else "", "\n", sep = "")
      }
    }
  }
}

cat("\n=== Recurrent convergent genes across all contrast/analysis pairs ===\n")
gene_recur <- all_pairs[convergent == TRUE, .(
  n_convergent_tests = .N,
  contrasts = paste(unique(contrast), collapse = "; "),
  analyses = paste(unique(analysis), collapse = "; "),
  directions = paste(unique(convergent_direction), collapse = "; "),
  any_both_FDR_same_direction = any(both_FDR_same_direction)
), by = .(gene_id_clean, gene_name)]
setorder(gene_recur, -n_convergent_tests, gene_name)
print(gene_recur[n_convergent_tests >= 2][1:50])

run_ora <- function(query_ids, label) {
  query_ids <- unique(na.omit(query_ids))
  if (length(query_ids) < min_hits) {
    cat("\n", label, ": ", length(query_ids), " genes; skipped ORA (<", min_hits, " genes)\n", sep = "")
    return(data.table())
  }
  cat("\n[g:Profiler convergent ORA] ", label, " (", length(query_ids), " genes)\n", sep = "")
  res <- tryCatch(
    gost(
      query = query_ids,
      organism = "hsapiens",
      sources = c("REAC", "GO:BP", "GO:CC"),
      correction_method = "g_SCS",
      domain_scope = "annotated",
      user_threshold = 0.05,
      evcodes = TRUE
    ),
    error = function(e) {
      cat("g:Profiler error: ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
    cat("No significant terms returned.\n")
    return(data.table())
  }
  dt <- as.data.table(res$result)
  dt <- dt[significant == TRUE & term_size <= max_term_size & intersection_size >= min_hits]
  if (nrow(dt) == 0) {
    cat("No terms after term_size/hit filters.\n")
    return(data.table())
  }
  dt[, label := label]
  print(dt[order(p_value)][1:min(.N, 10), .(source, term_name, p_value, term_size, intersection_size)])
  dt
}

cat("\n=== ORA on convergent gene sets ===\n")
ora_results <- list()
for (cn in contrasts$contrast) {
  for (an in pairs$analysis) {
    for (dirn in c("Up", "Down")) {
      q <- all_pairs[contrast == cn & analysis == an & convergent_direction == dirn, gene_id_clean]
      lbl <- paste(cn, an, dirn, sep = " | ")
      ora_results[[lbl]] <- run_ora(q, lbl)
    }
  }
}
conv_ora <- rbindlist(ora_results, fill = TRUE)

theme_terms <- list(
  "mRNA splicing / processing" = c("splic", "mrna processing", "rna processing", "spliceosome"),
  "Chromatin / histone regulation" = c("chromatin", "histone", "nucleosome", "heterochromatin", "deacetyl", "methylate histone"),
  "Telomere / DNA damage / repair" = c("telomere", "dna damage", "dna repair", "double-strand", "base-excision", "depurination", "depyrimidination"),
  "Cell cycle / mitosis / G2M" = c("cell cycle", "mitotic", "mitosis", "g2", "m phase", "chromosome segregation", "dna replication", "pre-replicative"),
  "Type I interferon / antiviral" = c("interferon", "antiviral", "defense response to virus", "response to virus"),
  "Golgi / vesicle traffic" = c("golgi", "vesicle", "er traffic", "vacuole", "lysosome")
)

assign_theme <- function(term_name) {
  low <- tolower(term_name)
  hits <- names(theme_terms)[vapply(theme_terms, function(keys) any(vapply(keys, grepl, logical(1), x = low, fixed = TRUE)), logical(1))]
  if (length(hits) == 0) "Other" else paste(hits, collapse = "; ")
}

read_existing <- function(f, framework, level) {
  if (!file.exists(f)) return(data.table())
  d <- fread(f)
  if (!all(c("term_name", "p_value") %in% names(d))) return(data.table())
  if (!"source" %in% names(d)) d[, source := NA_character_]
  if (!"contrast" %in% names(d)) d[, contrast := NA_character_]
  if (!"direction" %in% names(d)) d[, direction := NA_character_]
  if ("analysis" %in% names(d)) {
    d[, analysis_unit := analysis]
  } else if ("fraction" %in% names(d)) {
    d[, analysis_unit := fraction]
  } else {
    d[, analysis_unit := NA_character_]
  }
  d[, `:=`(framework = framework, evidence_level = level)]
  d[, .(framework, evidence_level, contrast, analysis_unit, direction, source, term_name, p_value)]
}

existing <- rbindlist(list(
  read_existing(file.path(count_dir, "GO_BP_ORA_All_Contrasts_Combined", "All_Tables", "filtered_terms.csv"), "P-site fraction limma", "rawP-derived ORA"),
  read_existing(file.path(metric_dir, "..", "GO_BP_ORA_Baseline_Combined", "Tables", "filtered_terms.csv"), "Translation metric limma", "rawP-derived ORA"),
  read_existing(file.path(metric_dir, "..", "GO_BP_ORA_Interaction_Combined", "Tables", "filtered_terms.csv"), "Translation metric limma", "rawP-derived ORA"),
  read_existing(file.path(metric_dir, "..", "Pathway_gProfiler_Clean_REAC_GOBP_GOCC_lfc0.7_rawP0.05", "Tables", "gprofiler_clean_filtered_terms_termSizeLE500_hitGE5.csv"), "Translation metric limma", "rawP REAC/GO ORA"),
  read_existing(file.path(count_dir, "Pathway_gProfiler_Clean_REAC_GOBP_GOCC_true_interaction_lfc0.7_rawP0.05", "Tables", "gprofiler_clean_filtered_terms_termSizeLE500_hitGE5.csv"), "P-site fraction limma", "rawP REAC/GO ORA")
), fill = TRUE)

if (nrow(conv_ora) > 0) {
  conv_existing <- conv_ora[, .(
    framework = "Convergent metric+Psite",
    evidence_level = "convergent ORA",
    contrast = tstrsplit(label, " \\| ", keep = 1L)[[1]],
    analysis_unit = tstrsplit(label, " \\| ", keep = 2L)[[1]],
    direction = tstrsplit(label, " \\| ", keep = 3L)[[1]],
    source, term_name, p_value
  )]
  existing <- rbindlist(list(existing, conv_existing), fill = TRUE)
}

if (nrow(existing) > 0) {
  existing[, theme := assign_theme(term_name), by = term_name]
  scored <- existing[theme != "Other", .(
    n_independent_analyses = uniqueN(paste(framework, evidence_level, contrast, analysis_unit, direction, sep = "|")),
    n_terms = uniqueN(term_name),
    best_p = min(p_value, na.rm = TRUE),
    sources = paste(sort(unique(source)), collapse = ", "),
    contexts = paste(head(unique(paste(framework, contrast, analysis_unit, direction, sep = " / ")), 10), collapse = "; ")
  ), by = theme]
  setorder(scored, -n_independent_analyses, best_p)
  cat("\n=== Pathway theme convergence score from existing ORA + convergent ORA ===\n")
  print(scored)
}

cat("\n=== Tier 1 candidate gene screen ===\n")
cat("Operational criteria used here: convergent same-direction in >=2 tests OR convergent with both metric and P-site FDR in same direction; then prioritized by recurrence and named pathway context.\n")
tier <- copy(gene_recur)
tier[, tier_flag := n_convergent_tests >= 2 | any_both_FDR_same_direction]
print(tier[tier_flag == TRUE][1:60])

cat("\n=== Short interpretive takeaways ===\n")
cat("1. The strongest cross-framework gene-level convergence is ribosome engagement vs RS and collision vs DS, especially at Resistance baseline.\n")
cat("2. Scanning vs SSU produces very few convergent genes, so scanning should be interpreted as a derived bottleneck metric, not bulk SSU abundance.\n")
cat("3. Pathways from convergent gene sets are stricter; absence of ORA in a small convergent set should be treated as limited power, not absence of biology.\n")
