# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

rdata <- input_path("SUDHL.RData")
gtf <- input_path("Homo_sapiens.GRCh38.114.chr.gtf")
out_dir <- analysis_path("LongRead_TRA2A", "Step1_TRA2A_isoform_usage")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

e <- new.env()
load(rdata, envir = e)

gtf_lines <- readLines(gtf)
tra2a_tx_lines <- gtf_lines[
  grepl("\ttranscript\t", gtf_lines) &
    grepl('gene_id "ENSG00000164548"', gtf_lines)
]
attrs <- sub("^([^\t]*\t){8}", "", tra2a_tx_lines)

grab_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[out == x] <- NA_character_
  out
}

tx <- data.table(
  transcript_id = grab_attr(attrs, "transcript_id"),
  transcript_version = grab_attr(attrs, "transcript_version"),
  transcript_name = grab_attr(attrs, "transcript_name"),
  transcript_biotype = grab_attr(attrs, "transcript_biotype")
)
tx[, row_id := paste0(transcript_id, ".", transcript_version)]
tx[, in_matrix := row_id %in% rownames(e$count_matrix)]

present_tx <- tx[in_matrix == TRUE]
counts <- e$count_matrix[present_tx$row_id, , drop = FALSE]
counts[is.na(counts)] <- 0

count_dt <- as.data.table(as.table(counts))
setnames(count_dt, c("row_id", "Sample", "count"))
count_dt <- merge(count_dt, present_tx, by = "row_id", all.x = TRUE)
count_dt <- merge(count_dt, as.data.table(e$Annotation), by = "Sample", all.x = TRUE)
count_dt[, count := as.numeric(count)]
count_dt <- count_dt[Condition == "SUDHL8"]

sample_totals <- count_dt[, .(TRA2A_total_count = sum(count)), by = .(Sample, Condition, Type, Replicate)]
usage_dt <- merge(count_dt, sample_totals, by = c("Sample", "Condition", "Type", "Replicate"))
usage_dt[, proportion := fifelse(TRA2A_total_count > 0, count / TRA2A_total_count, NA_real_)]

group_usage <- usage_dt[, .(
  mean_count = mean(count),
  median_count = median(count),
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
  mean_total_count = mean(TRA2A_total_count),
  median_total_count = median(TRA2A_total_count),
  sd_total_count = sd(TRA2A_total_count)
), by = .(Condition, Type)]
total_wide <- dcast(total_summary, Condition ~ Type, value.var = "mean_total_count")
if (!"Original" %in% names(total_wide)) total_wide[, Original := NA_real_]
if (!"Resistant" %in% names(total_wide)) total_wide[, Resistant := NA_real_]
total_wide[, resistant_minus_original := Resistant - Original]
total_wide[, resistant_over_original := Resistant / Original]

write.csv(tx, file.path(out_dir, "TRA2A_all_annotated_transcripts_matrix_presence.csv"), row.names = FALSE)
write.csv(usage_dt, file.path(out_dir, "TRA2A_isoform_usage_per_sample.csv"), row.names = FALSE)
write.csv(group_usage, file.path(out_dir, "TRA2A_isoform_usage_group_summary.csv"), row.names = FALSE)
write.csv(wide_usage[order(Condition, -abs(resistant_minus_original))], file.path(out_dir, "TRA2A_isoform_usage_resistant_vs_original_shift.csv"), row.names = FALSE)
write.csv(sample_totals, file.path(out_dir, "TRA2A_total_counts_per_sample.csv"), row.names = FALSE)
write.csv(total_wide, file.path(out_dir, "TRA2A_total_counts_resistant_vs_original_summary.csv"), row.names = FALSE)

plot_usage_dt <- usage_dt[
  ,
  .(proportion = sum(count) / unique(TRA2A_total_count)),
  by = .(Sample, Condition, Type, Replicate, transcript_name, transcript_biotype)
]
plot_usage_dt[, Type := factor(Type, levels = c("Original", "Resistant"))]
plot_usage_dt[, sample_label := factor(Sample, levels = e$Annotation$Sample[e$Annotation$Condition == "SUDHL8"])]

usage_plot <- ggplot(plot_usage_dt, aes(x = sample_label, y = proportion, fill = transcript_name)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.15) +
  facet_grid(. ~ Type, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "TRA2A isoform usage in SUDHL8 long-read data",
    subtitle = "Isoform count divided by total TRA2A transcript count per sample",
    x = NULL,
    y = "Within-gene isoform proportion",
    fill = "TRA2A isoform"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

total_plot <- ggplot(sample_totals, aes(x = Type, y = TRA2A_total_count, fill = Type)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
  geom_point(aes(group = Replicate), position = position_jitter(width = 0.08, height = 0), size = 2.4) +
  scale_fill_manual(values = c(Original = "#4C78A8", Resistant = "#D84A3A")) +
  labs(
    title = "Total TRA2A long-read transcript counts in SUDHL8",
    subtitle = "Sum of all TRA2A isoform counts per sample",
    x = NULL,
    y = "Total TRA2A count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "TRA2A_isoform_usage_stacked_proportions.png"), usage_plot, width = 13, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "TRA2A_total_RNA_counts.png"), total_plot, width = 7.5, height = 5, dpi = 300)

cat("\nTRA2A transcript presence:\n")
print(tx)

cat("\nLargest resistant-vs-original isoform usage shifts:\n")
print(wide_usage[order(Condition, -abs(resistant_minus_original))])

cat("\nTotal TRA2A RNA count summary:\n")
print(total_wide)

cat("\nOutput directory:\n")
cat(out_dir, "\n")
