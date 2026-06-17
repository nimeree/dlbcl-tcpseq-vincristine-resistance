# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)
library(ggplot2)
library(scales)

gene_query <- "LPXN"

rdata <- input_path("SUDHL.RData")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
out_dir <- analysis_path("LongRead_LPXN", "Step1_LPXN_isoform_usage")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

e <- new.env()
load(rdata, envir = e)

annotation <- as.data.table(e$Annotation)[Condition == "SUDHL8"]
samples <- annotation$Sample

gtf_dt <- fread(
  gtf,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
)

tx_dt <- gtf_dt[feature == "transcript"]
attrs <- tx_dt$attributes
tx <- data.table(
  gene_id = grab_attr(attrs, "gene_id"),
  gene_name = grab_attr(attrs, "gene_name"),
  gene_biotype = grab_attr(attrs, "gene_biotype"),
  transcript_id = grab_attr(attrs, "transcript_id"),
  transcript_version = grab_attr(attrs, "transcript_version"),
  transcript_name = grab_attr(attrs, "transcript_name"),
  transcript_biotype = grab_attr(attrs, "transcript_biotype")
)
tx[, row_id := paste0(transcript_id, ".", transcript_version)]
tx[, in_matrix := row_id %in% rownames(e$count_matrix)]

gene_tx <- tx[toupper(gene_name) == toupper(gene_query)]
if (!nrow(gene_tx)) stop("No annotated transcripts found for ", gene_query)

present_tx <- gene_tx[in_matrix == TRUE]
if (!nrow(present_tx)) stop("No count-matrix transcripts found for ", gene_query)

counts <- e$count_matrix[present_tx$row_id, samples, drop = FALSE]
counts[is.na(counts)] <- 0

count_dt <- as.data.table(as.table(counts))
setnames(count_dt, c("row_id", "Sample", "count"))
count_dt[, count := as.numeric(count)]
count_dt <- merge(count_dt, present_tx, by = "row_id", all.x = TRUE)
count_dt <- merge(count_dt, annotation, by = "Sample", all.x = TRUE)

sample_totals <- count_dt[, .(gene_total_count = sum(count)), by = .(Sample, Condition, Type, Replicate)]
usage_dt <- merge(count_dt, sample_totals, by = c("Sample", "Condition", "Type", "Replicate"))
usage_dt[, proportion := fifelse(gene_total_count > 0, count / gene_total_count, NA_real_)]

group_usage <- usage_dt[, .(
  mean_count = mean(count),
  median_count = median(count),
  total_count = sum(count),
  mean_proportion = mean(proportion, na.rm = TRUE),
  median_proportion = median(proportion, na.rm = TRUE)
), by = .(Condition, Type, row_id, transcript_name, transcript_biotype)]

wide_usage <- dcast(
  group_usage,
  Condition + row_id + transcript_name + transcript_biotype ~ Type,
  value.var = "mean_proportion"
)
if (!"Original" %in% names(wide_usage)) wide_usage[, Original := NA_real_]
if (!"Resistant" %in% names(wide_usage)) wide_usage[, Resistant := NA_real_]
wide_usage[, resistant_minus_original := Resistant - Original]
wide_usage[, fold_change_prop := fifelse(Original > 0, Resistant / Original, NA_real_)]

total_summary <- sample_totals[, .(
  n = .N,
  counts = paste(gene_total_count, collapse = ","),
  mean_total_count = mean(gene_total_count),
  median_total_count = median(gene_total_count),
  sd_total_count = sd(gene_total_count)
), by = .(Condition, Type)]

total_tests <- data.table(
  test = c("Wilcoxon rank-sum", "Welch t-test"),
  p_value = c(
    tryCatch(wilcox.test(gene_total_count ~ Type, data = sample_totals, exact = FALSE)$p.value, error = function(err) NA_real_),
    tryCatch(t.test(gene_total_count ~ Type, data = sample_totals)$p.value, error = function(err) NA_real_)
  )
)

usage_mat <- xtabs(count ~ Type + transcript_name, data = usage_dt)
usage_mat <- usage_mat[, colSums(usage_mat) > 0, drop = FALSE]
usage_tests <- data.table(test = character(), p_value = numeric())
if (nrow(usage_mat) >= 2 && ncol(usage_mat) >= 2) {
  usage_tests <- rbind(
    usage_tests,
    data.table(test = "Chi-square aggregated isoform usage", p_value = suppressWarnings(chisq.test(usage_mat)$p.value))
  )
  set.seed(1)
  usage_tests <- rbind(
    usage_tests,
    data.table(test = "Monte Carlo chi-square aggregated isoform usage", p_value = chisq.test(usage_mat, simulate.p.value = TRUE, B = 100000)$p.value)
  )
}

write.csv(gene_tx, file.path(out_dir, "LPXN_all_annotated_transcripts_matrix_presence.csv"), row.names = FALSE)
write.csv(usage_dt, file.path(out_dir, "LPXN_isoform_usage_per_sample.csv"), row.names = FALSE)
write.csv(group_usage, file.path(out_dir, "LPXN_isoform_usage_group_summary.csv"), row.names = FALSE)
write.csv(wide_usage[order(Condition, -abs(resistant_minus_original))], file.path(out_dir, "LPXN_isoform_usage_resistant_vs_original_shift.csv"), row.names = FALSE)
write.csv(sample_totals, file.path(out_dir, "LPXN_total_counts_per_sample.csv"), row.names = FALSE)
write.csv(total_summary, file.path(out_dir, "LPXN_total_counts_group_summary.csv"), row.names = FALSE)
write.csv(total_tests, file.path(out_dir, "LPXN_total_counts_tests.csv"), row.names = FALSE)
usage_mat_dt <- as.data.table(as.matrix(usage_mat), keep.rownames = "Type")
write.csv(usage_mat_dt, file.path(out_dir, "LPXN_isoform_usage_contingency_matrix.csv"), row.names = FALSE)
write.csv(usage_tests, file.path(out_dir, "LPXN_isoform_usage_tests.csv"), row.names = FALSE)

plot_usage_dt <- usage_dt[
  ,
  .(proportion = sum(count) / unique(gene_total_count)),
  by = .(Sample, Condition, Type, Replicate, transcript_name, transcript_biotype)
]
plot_usage_dt[, Type := factor(Type, levels = c("Original", "Resistant"))]
plot_usage_dt[, sample_label := factor(Sample, levels = annotation$Sample)]

usage_plot <- ggplot(plot_usage_dt, aes(x = sample_label, y = proportion, fill = transcript_name)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.15) +
  facet_grid(. ~ Type, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "LPXN isoform usage in SUDHL8 long-read data",
    subtitle = "Isoform count divided by total LPXN transcript count per sample",
    x = NULL,
    y = "Within-gene isoform proportion",
    fill = "LPXN isoform"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

total_plot <- ggplot(sample_totals, aes(x = Type, y = gene_total_count, fill = Type)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.4) +
  scale_fill_manual(values = c(Original = "#4C78A8", Resistant = "#D84A3A")) +
  labs(
    title = "Total LPXN long-read transcript counts in SUDHL8",
    subtitle = "Sum of all LPXN isoform counts per sample",
    x = NULL,
    y = "Total LPXN count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "LPXN_isoform_usage_stacked_proportions.png"), usage_plot, width = 13, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "LPXN_total_RNA_counts.png"), total_plot, width = 7.5, height = 5, dpi = 300)

cat("\nLPXN transcript presence:\n")
print(gene_tx)

cat("\nPer-sample total LPXN long-read counts:\n")
print(sample_totals[order(Type, Replicate)])

cat("\nTotal LPXN RNA count summary:\n")
print(total_summary)
print(total_tests)

cat("\nLargest resistant-vs-original isoform usage shifts:\n")
print(wide_usage[order(Condition, -abs(resistant_minus_original))])

cat("\nAggregated isoform counts by group:\n")
print(usage_mat_dt)

cat("\nIsoform usage tests:\n")
print(usage_tests)

cat("\nOutput directory:\n")
cat(out_dir, "\n")
