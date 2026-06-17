# Load portable path helpers when run from the repository root or scripts subfolders.
.local_config_candidates <- file.path(c(".", "..", "../.."), "config", "paths.R")
.local_config <- .local_config_candidates[file.exists(.local_config_candidates)][1]
if (!is.na(.local_config)) source(.local_config)
rm(.local_config, .local_config_candidates)

# ============================================================
# UpSet plots for P-site fraction-specific limma count outputs
# - 3-set UpSet plots: SSU, RS, DS
# - Separate Up and Down plots for each contrast
# - Includes true interaction:
#   (Resistant Vin - Resistant DMSO) - (Sensitive Vin - Sensitive DMSO)
# - Threshold already encoded in input files: raw P < 0.05 and |logFC| >= 0.7
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
})

BASE_DIR <- analysis_path()
IN_DIR <- file.path(BASE_DIR, "Psite_fraction_limma_lfc0.7_rawP0.05")
OUT_DIR <- file.path(IN_DIR, "UpSet_Plots_lfc0.7_rawP0.05")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

P_CUT <- 0.05
LFC_CUT <- 0.7
FRACTIONS <- c("SSU", "RS", "DS")

COMPARISONS <- c(
  "Resistance_baseline",
  "Sensitive_Vin_vs_DMSO",
  "Resistant_Vin_vs_DMSO",
  "Vin_Resistant_vs_Sensitive",
  "Interaction"
)

COMPARISON_TITLES <- c(
  Resistance_baseline = "Baseline resistance",
  Sensitive_Vin_vs_DMSO = "Sensitive Vin vs DMSO",
  Resistant_Vin_vs_DMSO = "Resistant Vin vs DMSO",
  Vin_Resistant_vs_Sensitive = "Vin Resistant vs Sensitive",
  Interaction = "True interaction"
)

COMPARISON_SUBTITLES <- c(
  Resistance_baseline = "P-site fraction limma: Resistant DMSO - Sensitive DMSO",
  Sensitive_Vin_vs_DMSO = "P-site fraction limma: Sensitive Vin - Sensitive DMSO",
  Resistant_Vin_vs_DMSO = "P-site fraction limma: Resistant Vin - Resistant DMSO",
  Vin_Resistant_vs_Sensitive = "P-site fraction limma: Resistant Vin - Sensitive Vin",
  Interaction = "P-site fraction limma: (Resistant Vin - Resistant DMSO) - (Sensitive Vin - Sensitive DMSO)"
)

read_gene_set <- function(comparison, fraction, direction) {
  f <- file.path(
    IN_DIR,
    paste0("Fraction_", fraction),
    paste0("Interaction_", fraction, "_psite_limma_sig_", tolower(direction), "_rawP0.05_lfc0.7.csv")
  )
  if (comparison != "Interaction") {
    f <- file.path(
      IN_DIR,
      paste0("Fraction_", fraction),
      paste0(comparison, "_", fraction, "_psite_limma_sig_", tolower(direction), "_rawP0.05_lfc0.7.csv")
    )
  }
  if (!file.exists(f)) stop("Missing input file: ", f)

  dt <- fread(f)
  if (!"gene_id_clean" %in% names(dt)) stop("gene_id_clean missing from: ", f)
  if (!"gene_name" %in% names(dt)) dt[, gene_name := NA_character_]
  dt[, gene_id_clean := sub("\\.\\d+$", "", gene_id_clean)]
  unique(dt[!is.na(gene_id_clean) & gene_id_clean != "", .(gene_id = gene_id_clean, gene_name)])
}

intersection_table <- function(sets) {
  all_genes <- unique(rbindlist(sets, fill = TRUE))
  all_genes <- unique(all_genes[order(gene_id)], by = "gene_id")
  mat <- copy(all_genes)

  for (fraction in FRACTIONS) {
    mat[, (fraction) := gene_id %in% sets[[fraction]]$gene_id]
  }

  mat[, intersection := fifelse(
    SSU & !RS & !DS, "SSU only",
    fifelse(!SSU & RS & !DS, "RS only",
      fifelse(!SSU & !RS & DS, "DS only",
        fifelse(SSU & RS & !DS, "SSU+RS",
          fifelse(SSU & !RS & DS, "SSU+DS",
            fifelse(!SSU & RS & DS, "RS+DS", "SSU+RS+DS")
          )
        )
      )
    )
  )]

  counts <- mat[, .(count = .N), by = .(intersection, SSU, RS, DS)]
  order_levels <- c("SSU+RS+DS", "SSU+RS", "SSU+DS", "RS+DS", "SSU only", "RS only", "DS only")
  counts[, intersection := factor(intersection, levels = order_levels)]
  setorder(counts, intersection)

  list(membership = mat, counts = counts)
}

make_upset_plot <- function(sets, title, subtitle, direction) {
  tbl <- intersection_table(sets)
  counts <- copy(tbl$counts)
  color <- if (direction == "Up") "#B8323B" else "#2B6CB0"

  set_sizes <- data.table(
    fraction = factor(FRACTIONS, levels = rev(FRACTIONS)),
    count = as.integer(vapply(FRACTIONS, function(x) nrow(sets[[x]]), integer(1)))
  )

  point_dt <- melt(
    counts[, .(intersection, SSU, RS, DS)],
    id.vars = "intersection",
    variable.name = "fraction",
    value.name = "present"
  )
  point_dt[, fraction := factor(fraction, levels = rev(FRACTIONS))]
  point_dt[, y_pos := as.integer(fraction)]

  line_dt <- point_dt[present == TRUE, .(
    y_min = min(y_pos),
    y_max = max(y_pos)
  ), by = intersection]

  top <- ggplot(counts, aes(x = intersection, y = count)) +
    geom_col(width = 0.68, fill = color, alpha = 0.92) +
    geom_text(aes(label = count), vjust = -0.35, size = 4.2, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    coord_cartesian(clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Intersection genes"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
      plot.subtitle = element_text(size = 11.5, hjust = 0.5, color = "grey30"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 14, 0, 14)
    )

  matrix <- ggplot(point_dt, aes(x = intersection, y = y_pos)) +
    geom_segment(
      data = line_dt[y_min != y_max],
      aes(x = intersection, xend = intersection, y = y_min, yend = y_max),
      inherit.aes = FALSE,
      linewidth = 1.1,
      color = color,
      alpha = 0.8
    ) +
    geom_point(
      aes(fill = present),
      shape = 21,
      size = 5.2,
      stroke = 0.55,
      color = "grey25"
    ) +
    scale_fill_manual(values = c(`TRUE` = color, `FALSE` = "grey88"), guide = "none") +
    scale_y_continuous(
      breaks = seq_along(rev(FRACTIONS)),
      labels = rev(FRACTIONS),
      expand = expansion(add = 0.32)
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(size = 12, color = "black", angle = 35, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black", face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(0, 14, 6, 14)
    )

  sizes <- ggplot(set_sizes, aes(x = count, y = fraction)) +
    geom_col(width = 0.55, fill = "grey35") +
    geom_text(aes(label = count), hjust = -0.2, size = 4) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    coord_cartesian(clip = "off") +
    labs(x = "Set size", y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0, 14, 6, 0)
    )

  top_row <- top + plot_spacer() + plot_layout(widths = c(4.8, 1.25))
  bottom <- matrix + sizes + plot_layout(widths = c(4.8, 1.25))
  top_row / bottom + plot_layout(heights = c(3.2, 1.65))
}

summary_rows <- list()

for (comparison in COMPARISONS) {
  for (direction in c("Up", "Down")) {
    sets <- setNames(
      lapply(FRACTIONS, function(frac) read_gene_set(comparison, frac, direction)),
      FRACTIONS
    )
    tbl <- intersection_table(sets)
    counts <- tbl$counts
    membership <- tbl$membership

    title <- paste(COMPARISON_TITLES[[comparison]], direction, "genes")
    subtitle <- paste0(
      COMPARISON_SUBTITLES[[comparison]],
      "; raw P < ", P_CUT, " and ",
      ifelse(direction == "Up", "logFC >= ", "logFC <= -"),
      LFC_CUT
    )

    p <- make_upset_plot(sets, title, subtitle, direction)

    out_base <- paste0(comparison, "_", direction, "_SSU_RS_DS_psite_limma_upset_lfc0.7_rawP0.05")
    ggsave(file.path(OUT_DIR, paste0(out_base, ".png")), p, width = 10.5, height = 7.4, dpi = 300, bg = "white")
    ggsave(file.path(OUT_DIR, paste0(out_base, ".pdf")), p, width = 10.5, height = 7.4, bg = "white")

    wb <- createWorkbook()
    for (intersection_name in as.character(counts$intersection)) {
      addWorksheet(wb, intersection_name)
      genes <- membership[intersection == intersection_name, .(gene_id, gene_name)]
      writeDataTable(wb, intersection_name, genes[order(gene_name, gene_id)], tableStyle = "TableStyleMedium2")
      freezePane(wb, intersection_name, firstRow = TRUE)
      setColWidths(wb, intersection_name, cols = 1:2, widths = "auto")
    }
    saveWorkbook(
      wb,
      file.path(OUT_DIR, paste0(out_base, "_gene_lists.xlsx")),
      overwrite = TRUE
    )

    summary_rows[[paste(comparison, direction, sep = "_")]] <- data.table(
      comparison = comparison,
      direction = direction,
      total_SSU = nrow(sets$SSU),
      total_RS = nrow(sets$RS),
      total_DS = nrow(sets$DS),
      intersection = as.character(counts$intersection),
      count = counts$count
    )
  }
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, file.path(OUT_DIR, "psite_limma_upset_intersection_counts_summary.csv"))

cat("\nDone. P-site limma UpSet plots and exact-intersection gene lists saved to:\n", OUT_DIR, "\n\n", sep = "")
print(summary_dt)
