# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# Regenerate combined GO:BP ORA plots with Word-friendly pathway labels.
# This uses the already-saved terms_shown_in_plot.csv files, so ORA results
# are unchanged; only figure typography/layout is updated.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(stringr)
})

BASE_DIR <- analysis_path()

table_files <- sort(system(
  paste(
    "find",
    shQuote(BASE_DIR),
    "-path '*GO_BP_ORA*Combined*' -type f -name '*terms_shown_in_plot.csv'"
  ),
  intern = TRUE
))

table_files <- table_files[!grepl("/All_Tables/", table_files, fixed = TRUE)]

if (length(table_files) == 0) {
  stop("No terms_shown_in_plot.csv files found.")
}

read_terms <- function(f) {
  if (!file.exists(f) || file.info(f)$size == 0) return(data.table())
  dt <- tryCatch(fread(f), error = function(e) data.table())
  dt
}

infer_root <- function(table_file) {
  sub("/Tables/[^/]+$", "", table_file)
}

infer_key <- function(table_file) {
  sub("_terms_shown_in_plot\\.csv$", "", basename(table_file))
}

find_plot_base <- function(table_file) {
  root <- infer_root(table_file)
  key <- infer_key(table_file)
  plot_dir <- file.path(root, "Plots")
  if (!dir.exists(plot_dir)) return(NA_character_)
  candidates <- list.files(plot_dir, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
  candidates <- candidates[startsWith(basename(candidates), paste0(key, "_GO_BP_ORA"))]
  candidates <- candidates[!grepl("_large_labels_preview\\.png$", candidates)]
  if (length(candidates) == 0) return(NA_character_)
  sub("\\.png$", "", candidates[[1]])
}

nice_title_from_data <- function(dt, table_file) {
  if (nrow(dt) > 0) {
    if ("analysis_label" %in% names(dt) && !is.na(dt$analysis_label[1]) && dt$analysis_label[1] != "") {
      return(paste(dt$analysis_label[1], "GO:BP ORA"))
    }
    if (all(c("contrast_label", "fraction") %in% names(dt)) &&
        !is.na(dt$contrast_label[1]) && dt$contrast_label[1] != "") {
      return(paste(dt$contrast_label[1], dt$fraction[1], "GO:BP ORA"))
    }
  }

  key <- infer_key(table_file)
  root <- infer_root(table_file)
  context <- basename(root)
  paste(gsub("_", " ", key), context, "GO:BP ORA")
}

subtitle_from_data <- function(dt) {
  if (nrow(dt) > 0 && all(c("user_threshold", "max_term_size_filter", "min_intersection_filter") %in% names(dt))) {
    return(paste0(
      "GO Biological Process only; Up + Down combined; g:SCS < ", dt$user_threshold[1],
      "; term size <= ", dt$max_term_size_filter[1],
      "; hit genes >= ", dt$min_intersection_filter[1]
    ))
  }
  "GO Biological Process only; Up + Down combined; g:SCS < 0.05; term size <= 500; hit genes >= 5"
}

plot_large <- function(dt, title, subtitle, out_base) {
  if (nrow(dt) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No GO Biological Process terms after filters", size = 7) +
      labs(title = title) +
      theme_void(base_size = 16) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20))
    ggsave(paste0(out_base, ".png"), p, width = 13.5, height = 7.2, dpi = 300, bg = "white")
    ggsave(paste0(out_base, ".pdf"), p, width = 13.5, height = 7.2, bg = "white")
    return(invisible(NULL))
  }

  if (!"neglog10_gscs" %in% names(dt)) dt[, neglog10_gscs := -log10(p_value)]
  dt <- dt[order(p_value, -intersection_size)]
  dt[, term_label := stringr::str_wrap(term_name, width = 34)]

  p <- ggplot(dt, aes(
    x = neglog10_gscs,
    y = forcats::fct_reorder(term_label, neglog10_gscs)
  )) +
    geom_segment(aes(
      x = 0, xend = neglog10_gscs,
      y = forcats::fct_reorder(term_label, neglog10_gscs),
      yend = forcats::fct_reorder(term_label, neglog10_gscs)
    ), linewidth = 0.9, color = "grey76") +
    geom_point(aes(size = intersection_size, color = mean_hit_log2FC), alpha = 0.95) +
    scale_color_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
    scale_size_continuous(range = c(5, 10)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "-log10(g:SCS corrected p-value)",
      y = NULL,
      size = "Hit genes",
      color = "Mean hit\nlogFC"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 13),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(color = "black", size = 15, lineheight = 0.95),
      panel.grid.major.y = element_blank(),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      plot.margin = margin(t = 10, r = 18, b = 10, l = 16)
    )

  height <- max(7.2, min(13.5, 3.0 + 0.58 * nrow(dt)))
  ggsave(paste0(out_base, ".png"), p, width = 13.5, height = height, dpi = 300, limitsize = FALSE, bg = "white")
  ggsave(paste0(out_base, ".pdf"), p, width = 13.5, height = height, limitsize = FALSE, bg = "white")
}

updated <- data.table()
skipped <- data.table()

for (f in table_files) {
  out_base <- find_plot_base(f)
  if (is.na(out_base)) {
    skipped <- rbind(skipped, data.table(table_file = f, reason = "No matching PNG plot found"))
    next
  }

  dt <- read_terms(f)
  plot_large(dt, nice_title_from_data(dt, f), subtitle_from_data(dt), out_base)
  updated <- rbind(updated, data.table(table_file = f, plot_base = out_base, rows = nrow(dt)))
}

cat("Updated combined ORA plots:", nrow(updated), "\n")
if (nrow(skipped) > 0) {
  cat("Skipped:", nrow(skipped), "\n")
  print(skipped)
}
cat("\nUpdated plot bases:\n")
print(updated)
