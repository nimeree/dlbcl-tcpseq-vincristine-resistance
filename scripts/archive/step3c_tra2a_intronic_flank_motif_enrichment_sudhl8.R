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
out_dir <- file.path(base_out, "Step3C_TRA2A_intronic_flank_motif_enrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rdata <- input_path("SUDHL.RData")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
genome_fasta <- input_path("Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
screen_file <- file.path(step3_dir, "SUDHL8_genomewide_isoform_usage_screen.csv")

flank_bp <- 250L
motifs <- c(AGAA = "AGAA", GAAGAA = "GAAGAA", GAA = "GAA")

standard_chr <- function(x) {
  x <- sub("^chr", "", x)
  fifelse(x == "M", "MT", x)
}

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

message("Loading Step 3 screen and SUDHL8 expressed transcripts...")
screen <- fread(screen_file)
e <- new.env()
load(rdata, envir = e)
annotation <- as.data.table(e$Annotation)[Condition == "SUDHL8"]
samples <- annotation$Sample
counts <- e$count_matrix[, samples, drop = FALSE]
counts[is.na(counts)] <- 0
expressed_tx <- rownames(counts)[rowSums(counts) > 0]

message("Parsing GTF transcript/exon annotation...")
gtf_dt <- fread(
  gtf,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
)

tx_dt <- gtf_dt[feature == "transcript"]
tx_attrs <- tx_dt$attributes
tx_map <- data.table(
  chr = standard_chr(tx_dt$chr),
  strand = tx_dt$strand,
  gene_id = sub("[.][0-9]+$", "", grab_attr(tx_attrs, "gene_id")),
  gene_name = grab_attr(tx_attrs, "gene_name"),
  transcript_id = grab_attr(tx_attrs, "transcript_id"),
  transcript_version = grab_attr(tx_attrs, "transcript_version"),
  transcript_name = grab_attr(tx_attrs, "transcript_name")
)
tx_map[, row_id := paste0(transcript_id, ".", transcript_version)]
tx_map <- tx_map[row_id %in% expressed_tx]

exon_dt <- gtf_dt[feature == "exon"]
exon_attrs <- exon_dt$attributes
exons <- data.table(
  chr = standard_chr(exon_dt$chr),
  start = as.integer(exon_dt$start),
  end = as.integer(exon_dt$end),
  strand = exon_dt$strand,
  gene_id = sub("[.][0-9]+$", "", grab_attr(exon_attrs, "gene_id")),
  gene_name = grab_attr(exon_attrs, "gene_name"),
  transcript_id = grab_attr(exon_attrs, "transcript_id"),
  transcript_version = grab_attr(exon_attrs, "transcript_version"),
  exon_number = suppressWarnings(as.integer(grab_attr(exon_attrs, "exon_number")))
)
exons[, row_id := paste0(transcript_id, ".", transcript_version)]
exons <- exons[row_id %in% tx_map$row_id]
exons <- merge(exons, tx_map[, .(row_id, transcript_name)], by = "row_id", all.x = TRUE)

exon_counts <- exons[, .(n_exons = .N), by = row_id]
exons <- merge(exons, exon_counts, by = "row_id", all.x = TRUE)
exons[, is_internal_exon := exon_number > 1 & exon_number < n_exons]
exons <- exons[is_internal_exon == TRUE]

message("Loading genome FASTA...")
genome <- readDNAStringSet(genome_fasta)
names(genome) <- sub(" .*", "", names(genome))
seq_lengths <- width(genome)
names(seq_lengths) <- names(genome)

message("Defining foreground and background exons...")
stringent_genes <- screen[switched_stringent == TRUE, gene_id]
background_genes <- screen[switched_stringent == FALSE, gene_id]

foreground <- exons[gene_id %in% stringent_genes]
foreground[, set := "Stringent switched"]

background <- exons[gene_id %in% background_genes]
background[, set := "Background"]

# Keep the background manageable and reproducible while preserving broad coverage.
set.seed(114)
max_background_exons <- min(nrow(background), max(50000L, nrow(foreground) * 5L))
if (nrow(background) > max_background_exons) {
  background <- background[sample(.N, max_background_exons)]
}

flank_exons <- rbindlist(list(foreground, background), use.names = TRUE, fill = TRUE)
flank_exons <- flank_exons[chr %in% names(genome)]
flank_exons[, seq_length := seq_lengths[chr]]

make_flanks <- function(dt) {
  left <- copy(dt)
  left[, flank_side := "acceptor_upstream"]
  left[, flank_start := pmax(1L, start - flank_bp)]
  left[, flank_end := start - 1L]

  right <- copy(dt)
  right[, flank_side := "donor_downstream"]
  right[, flank_start := end + 1L]
  right[, flank_end := pmin(seq_length, end + flank_bp)]

  out <- rbindlist(list(left, right), use.names = TRUE, fill = TRUE)
  out <- out[flank_start <= flank_end & flank_end <= seq_length]
  out
}

flanks <- make_flanks(flank_exons)
flanks[, flank_width := flank_end - flank_start + 1L]
flanks <- flanks[flank_width >= 50]

message("Extracting ", nrow(flanks), " intronic flank sequences...")
extract_one <- function(chr, start, end, strand) {
  s <- subseq(genome[[chr]], start = start, end = end)
  if (strand == "-") s <- reverseComplement(s)
  as.character(s)
}
flanks[, sequence := mapply(extract_one, chr, flank_start, flank_end, strand, USE.NAMES = FALSE)]
flanks[, sequence := toupper(sequence)]

message("Counting TRA2A-like motifs...")
for (nm in names(motifs)) {
  pat <- motifs[[nm]]
  count_col <- paste0("count_", nm)
  present_col <- paste0("present_", nm)
  flanks[, (count_col) := vcountPattern(pat, DNAStringSet(sequence), fixed = TRUE)]
  flanks[, (present_col) := get(count_col) > 0]
}

gene_summary <- flanks[
  ,
  c(
    .(
      gene_name = gene_name[1],
      set = set[1],
      n_flanks = .N,
      total_flank_bp = sum(flank_width)
    ),
    lapply(.SD, sum)
  ),
  by = gene_id,
  .SDcols = patterns("^count_")
]
for (nm in names(motifs)) {
  gene_summary[, paste0("density_per_kb_", nm) := get(paste0("count_", nm)) / total_flank_bp * 1000]
  gene_summary[, paste0("present_", nm) := get(paste0("count_", nm)) > 0]
}

test_motif <- function(motif_name) {
  dens_col <- paste0("density_per_kb_", motif_name)
  present_col <- paste0("present_", motif_name)
  top_quartile <- gene_summary[[dens_col]] >= quantile(gene_summary[set == "Background", get(dens_col)], 0.75, na.rm = TRUE)

  wt <- wilcox.test(
    gene_summary[set == "Stringent switched", get(dens_col)],
    gene_summary[set == "Background", get(dens_col)],
    alternative = "greater"
  )
  present_tab <- table(
    set = factor(gene_summary$set, levels = c("Background", "Stringent switched")),
    motif_present = factor(gene_summary[[present_col]], levels = c(FALSE, TRUE))
  )
  top_tab <- table(
    set = factor(gene_summary$set, levels = c("Background", "Stringent switched")),
    motif_top_background_quartile = factor(top_quartile, levels = c(FALSE, TRUE))
  )
  ft_present <- fisher.test(present_tab)
  ft_top <- fisher.test(top_tab)

  data.table(
    motif = motif_name,
    switched_genes = gene_summary[set == "Stringent switched", .N],
    background_genes = gene_summary[set == "Background", .N],
    switched_mean_density_per_kb = mean(gene_summary[set == "Stringent switched", get(dens_col)], na.rm = TRUE),
    background_mean_density_per_kb = mean(gene_summary[set == "Background", get(dens_col)], na.rm = TRUE),
    switched_median_density_per_kb = median(gene_summary[set == "Stringent switched", get(dens_col)], na.rm = TRUE),
    background_median_density_per_kb = median(gene_summary[set == "Background", get(dens_col)], na.rm = TRUE),
    wilcox_greater_pvalue = wt$p.value,
    motif_present_odds_ratio = unname(ft_present$estimate),
    motif_present_fisher_pvalue = ft_present$p.value,
    top_quartile_odds_ratio = unname(ft_top$estimate),
    top_quartile_fisher_pvalue = ft_top$p.value
  )
}

enrichment <- rbindlist(lapply(names(motifs), test_motif))
enrichment[, wilcox_padj := p.adjust(wilcox_greater_pvalue, method = "BH")]
enrichment[, motif_present_fisher_padj := p.adjust(motif_present_fisher_pvalue, method = "BH")]
enrichment[, top_quartile_fisher_padj := p.adjust(top_quartile_fisher_pvalue, method = "BH")]

write.csv(flanks[, !"sequence"], file.path(out_dir, "intronic_flanks_motif_counts_no_sequence.csv"), row.names = FALSE)
write.csv(gene_summary, file.path(out_dir, "gene_level_intronic_flank_TRA2A_motif_density.csv"), row.names = FALSE)
write.csv(enrichment, file.path(out_dir, "TRA2A_intronic_flank_motif_enrichment.csv"), row.names = FALSE)

plot_dt <- melt(
  gene_summary,
  id.vars = c("gene_id", "gene_name", "set"),
  measure.vars = paste0("density_per_kb_", names(motifs)),
  variable.name = "motif",
  value.name = "density_per_kb"
)
plot_dt[, motif := sub("^density_per_kb_", "", motif)]
plot_dt[, set := factor(set, levels = c("Background", "Stringent switched"))]

p <- ggplot(plot_dt, aes(x = set, y = density_per_kb, fill = set)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.78, width = 0.65) +
  geom_jitter(width = 0.18, alpha = 0.12, size = 0.45) +
  facet_wrap(~ motif, scales = "free_y") +
  scale_fill_manual(values = c("Background" = "grey72", "Stringent switched" = "#7B3294")) +
  labs(
    title = "TRA2A-like motif density in intronic flanks",
    subtitle = paste0(flank_bp, " bp flanks around internal exons; SUDHL8 switched genes vs non-switched background"),
    x = NULL,
    y = "Motif occurrences per kb intronic flank"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "TRA2A_intronic_flank_motif_density_stringent_switches.png"), p, width = 9, height = 5.5, dpi = 300)

readme <- c(
  "Step 3C: TRA2A-like motif enrichment in intronic flanking sequence.",
  paste0("Genome FASTA: ", genome_fasta),
  paste0("GTF: ", gtf),
  "Foreground: genes classified as stringent SUDHL8 isoform switches in Step 3.",
  "Background: non-switched eligible expressed multi-isoform genes from Step 3.",
  paste0("Sequences: ", flank_bp, " bp upstream and downstream of internal exons from expressed transcripts; minus-strand sequences reverse-complemented to transcript orientation."),
  "Motifs tested: AGAA, GAAGAA, GAA.",
  "Tests: one-sided Wilcoxon test for higher motif density in switched genes; Fisher tests for motif presence and high motif density relative to background top quartile.",
  "Caveat: this approximates regulated regions using internal exon flanks from switched genes. A full exon-level DIU model would refine the foreground to specific included/excluded exons."
)
writeLines(readme, file.path(out_dir, "README_step3c_intronic_flank_motif_method.txt"))

cat("\nForeground genes with flanks:", gene_summary[set == "Stringent switched", .N], "\n")
cat("Background genes with flanks:", gene_summary[set == "Background", .N], "\n")
cat("Foreground flanks:", flanks[set == "Stringent switched", .N], "\n")
cat("Background flanks:", flanks[set == "Background", .N], "\n")
cat("\nIntronic flank motif enrichment:\n")
print(enrichment)
cat("\nOutput directory:\n")
cat(out_dir, "\n")
