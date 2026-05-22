# =============================================================================
# R/generate_figure_s1_pdf.R
# Render Figure S1 (sensitivity forest plot) as vector PDF
# =============================================================================
# Reads output/sensitivity/sensitivity_summary.csv and produces
# output/sensitivity/sensitivity_forest.pdf using cairo_pdf (vector).
# Plot specification mirrors generate_forest_plot() in sensitivity_analysis.R.
# =============================================================================

if (file.exists("/analysis/R")) setwd("/analysis")

summary_path <- file.path("output", "sensitivity", "sensitivity_summary.csv")
pdf_path     <- file.path("output", "sensitivity", "sensitivity_forest.pdf")

if (!file.exists(summary_path)) {
  stop("Missing input: ", summary_path)
}

summary_df <- read.csv(summary_path, stringsAsFactors = FALSE)

generate_forest_pdf <- function(summary_df, output_file) {
  plot_data <- summary_df[!is.na(summary_df$HapB3_OR), ]
  if (nrow(plot_data) == 0) stop("No valid HapB3_OR rows in summary CSV.")

  is_loo     <- grepl("^loo_", plot_data$scenario)
  is_primary <- plot_data$scenario == "primary"

  plot_data$order <- ifelse(is_primary, 1, ifelse(!is_loo, 2, 3))
  plot_data <- plot_data[order(plot_data$order, plot_data$scenario), ]

  is_loo     <- grepl("^loo_", plot_data$scenario)
  is_primary <- plot_data$scenario == "primary"

  plot_data$label <- ifelse(
    is_loo,
    gsub("^loo_", "LOO: ", plot_data$scenario),
    plot_data$description
  )

  is_nonconverged <- !is.na(plot_data$converged) & plot_data$converged == FALSE

  plot_data$color <- ifelse(
    is_primary, "blue",
    ifelse(is_nonconverged, "darkorange",
    ifelse(!is.na(plot_data$conclusion) & plot_data$conclusion == "Materially different",
           "red", "gray40"))
  )

  plot_data$pch <- ifelse(is_nonconverged, 17, ifelse(is_loo, 1, 19))

  primary_or <- plot_data$HapB3_OR[plot_data$scenario == "primary"]
  if (length(primary_or) == 0) primary_or <- 2.005

  n <- nrow(plot_data)

  cairo_pdf(output_file, width = 8, height = max(6, 0.4 * n + 2))
  on.exit(dev.off(), add = TRUE)

  par(mar = c(5, 15, 3, 2))

  x_max <- max(8, max(plot_data$HapB3_CrI_high, na.rm = TRUE) * 1.1)
  x_range <- c(0.5, x_max)

  plot(NA, xlim = x_range, ylim = c(0.5, n + 0.5),
       xlab = "Odds Ratio (log scale)", ylab = "",
       yaxt = "n", log = "x",
       main = "Sensitivity Analysis: HapB3 OR vs WT")

  axis(2, at = n:1, labels = plot_data$label, las = 1, cex.axis = 0.7)

  abline(v = 1.0, lty = 2, col = "gray60")
  abline(v = primary_or, lty = 3, col = "blue")

  for (i in seq_len(n)) {
    y_pos <- n - i + 1
    segments(
      x0 = plot_data$HapB3_CrI_low[i],
      x1 = plot_data$HapB3_CrI_high[i],
      y0 = y_pos, y1 = y_pos,
      col = plot_data$color[i], lwd = 1.5
    )
    points(plot_data$HapB3_OR[i], y_pos,
           pch = plot_data$pch[i], col = plot_data$color[i],
           bg = plot_data$color[i], cex = 1.2)
  }

  # Build legend dynamically; include only colour categories that appear in the
  # plotted data so we don't advertise a category (e.g. "Materially different")
  # that has zero scenarios.
  has_nonconverged <- any(is_nonconverged)
  has_material <- any(!is.na(plot_data$conclusion) &
                      plot_data$conclusion == "Materially different")

  legend_labels <- c("Primary", "Alternative scenario", "LOO study")
  legend_pch    <- c(19, 19, 1)
  legend_col    <- c("blue", "gray40", "gray40")
  legend_lty    <- c(NA, NA, NA)
  if (has_nonconverged) {
    legend_labels <- c(legend_labels, "Non-converged")
    legend_pch    <- c(legend_pch, 17)
    legend_col    <- c(legend_col, "darkorange")
    legend_lty    <- c(legend_lty, NA)
  }
  if (has_material) {
    legend_labels <- c(legend_labels, "Materially different")
    legend_pch    <- c(legend_pch, 19)
    legend_col    <- c(legend_col, "red")
    legend_lty    <- c(legend_lty, NA)
  }
  legend_labels <- c(legend_labels, "OR = 1.0", "Primary HapB3 OR")
  legend_pch    <- c(legend_pch, NA, NA)
  legend_col    <- c(legend_col, "gray60", "blue")
  legend_lty    <- c(legend_lty, 2, 3)

  legend("topright",
         legend = legend_labels,
         pch = legend_pch,
         col = legend_col,
         lty = legend_lty,
         cex = 0.7, bg = "white")
}

generate_forest_pdf(summary_df, pdf_path)
cat(sprintf("Saved: %s\n", pdf_path))
