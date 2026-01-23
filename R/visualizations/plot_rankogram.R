# =============================================================================
# plot_rankogram.R
# Generate Figure 6: Rank-O-Gram with SUCRA values
# =============================================================================
# This script creates a combined visualization showing:
#   Panel A: Cumulative ranking probability curves
#   Panel B: SUCRA (Surface Under the Cumulative Ranking) bar chart
# =============================================================================

# Source constants for colors, labels, and theme (if not already loaded)
if (!exists("VARIANT_COLORS")) {
  source("R/visualizations/constants.R")
}

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# -----------------------------------------------------------------------------
#' Generate Rank-O-Gram with SUCRA Values
#'
#' Creates a combined figure showing cumulative ranking probabilities and
#' SUCRA values for treatment ranking in the DPYD NMA.
#'
#' @param result_rds_path Path to the pcnetmeta result RDS file
#' @param output_dir Directory for saving output files
#' @return A ggplot object containing the combined figure
#' @export
# -----------------------------------------------------------------------------

plot_rankogram <- function(result_rds_path = "output/wt_unified/wt_unified_result.rds",
                           output_dir = "output/wt_unified") {

  cat("=== Generating Figure 6: Rank-O-Gram with SUCRA ===\n\n")

  # ---------------------------------------------------------------------------
  # 1. Load and extract rank probabilities
  # ---------------------------------------------------------------------------

  cat("Loading result file...\n")
  if (!file.exists(result_rds_path)) {
    stop("Result file not found: ", result_rds_path)
  }

  result <- readRDS(result_rds_path)

  # Extract TrtRankProb matrix
  if (!"TrtRankProb" %in% names(result)) {
    stop("TrtRankProb not found in result object")
  }

  rank_probs_raw <- result$TrtRankProb
  cat("  Rank probability matrix dimensions:", dim(rank_probs_raw), "\n")

  # Convert to numeric matrix
  ntrt <- nrow(rank_probs_raw)
  rank_probs <- matrix(as.numeric(rank_probs_raw), nrow = ntrt, ncol = ntrt)

  # Get treatment names - try multiple possible file patterns
  base_dir <- dirname(result_rds_path)
  possible_paths <- c(
    file.path(base_dir, "wt_unified_treatment_names.rds"),
    file.path(base_dir, "treatment_names.rds")
  )

  trt_names <- NULL
  for (path in possible_paths) {
    if (file.exists(path)) {
      trt_names <- readRDS(path)
      break
    }
  }

  if (is.null(trt_names)) {
    # Default treatment names if no file found
    trt_names <- c("WT", "HapB3", "2846hetho", "2Ahetho", "13hetho")
    cat("  Warning: Using default treatment names\n")
  }

  rownames(rank_probs) <- trt_names
  colnames(rank_probs) <- paste0("Rank", 1:ntrt)

  cat("  Treatments:", paste(trt_names, collapse = ", "), "\n\n")

  # ---------------------------------------------------------------------------
  # 2. Calculate SUCRA for each treatment
  # ---------------------------------------------------------------------------

  cat("Calculating SUCRA values...\n")

  # SUCRA = (1/(n-1)) * sum(cumulative_prob[1:(n-1)])
  # where cumulative_prob[k] = sum of rank_prob from rank 1 to k

  # Calculate cumulative probabilities
  cumulative_probs <- t(apply(rank_probs, 1, cumsum))
  colnames(cumulative_probs) <- paste0("Rank", 1:ntrt)

  # Calculate SUCRA for each treatment
  sucra <- apply(cumulative_probs[, 1:(ntrt - 1), drop = FALSE], 1, sum) / (ntrt - 1)

  # Create SUCRA data frame
  sucra_df <- data.frame(
    Treatment = names(sucra),
    SUCRA = sucra,
    SUCRA_pct = sucra * 100,
    stringsAsFactors = FALSE
  )

  # Apply treatment labels
  sucra_df$Label <- TREATMENT_LABELS[sucra_df$Treatment]

  cat("  SUCRA values:\n")
  for (i in 1:nrow(sucra_df)) {
    cat(sprintf("    %s: %.1f%%\n", sucra_df$Label[i], sucra_df$SUCRA_pct[i]))
  }
  cat("\n")

  # ---------------------------------------------------------------------------
  # 3. Reshape cumulative probabilities for ggplot
  # ---------------------------------------------------------------------------

  cumulative_df <- as.data.frame(cumulative_probs)
  cumulative_df$Treatment <- rownames(cumulative_df)

  cumulative_long <- cumulative_df %>%
    pivot_longer(
      cols = starts_with("Rank"),
      names_to = "Rank",
      values_to = "CumulativeProbability"
    ) %>%
    mutate(
      Rank = as.integer(gsub("Rank", "", Rank)),
      Label = TREATMENT_LABELS[Treatment]
    )

  # Set factor order for consistent plotting
  cumulative_long$Treatment <- factor(
    cumulative_long$Treatment,
    levels = TREATMENT_ORDER
  )
  cumulative_long$Label <- factor(
    cumulative_long$Label,
    levels = TREATMENT_LABELS[TREATMENT_ORDER]
  )

  # ---------------------------------------------------------------------------
  # 4. Create Panel A: Cumulative ranking curves
  # ---------------------------------------------------------------------------

  cat("Creating Panel A: Cumulative ranking curves...\n")

  p_cumulative <- ggplot(cumulative_long,
                         aes(x = Rank, y = CumulativeProbability,
                             color = Treatment, group = Treatment)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(
      values = VARIANT_COLORS,
      labels = TREATMENT_LABELS,
      name = "Treatment"
    ) +
    scale_x_continuous(breaks = 1:ntrt) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = "A. Cumulative Ranking Probability",
      x = "Rank (1 = lowest toxicity)",
      y = "Cumulative Probability"
    ) +
    theme_publication() +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  # ---------------------------------------------------------------------------
  # 5. Create Panel B: SUCRA bar chart
  # ---------------------------------------------------------------------------

  cat("Creating Panel B: SUCRA bar chart...\n")

  # Order treatments by SUCRA value (highest first for coord_flip)
  sucra_df$Treatment <- factor(
    sucra_df$Treatment,
    levels = sucra_df$Treatment[order(sucra_df$SUCRA)]
  )
  sucra_df$Label <- factor(
    sucra_df$Label,
    levels = TREATMENT_LABELS[as.character(sucra_df$Treatment[order(sucra_df$SUCRA)])]
  )

  p_sucra <- ggplot(sucra_df,
                    aes(x = Treatment, y = SUCRA_pct, fill = Treatment)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = sprintf("%.0f%%", SUCRA_pct)),
      hjust = -0.2,
      size = 4
    ) +
    scale_fill_manual(
      values = VARIANT_COLORS,
      guide = "none"
    ) +
    scale_x_discrete(labels = TREATMENT_LABELS) +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
    coord_flip() +
    labs(
      title = "B. SUCRA Values",
      x = NULL,
      y = "SUCRA (%)"
    ) +
    theme_publication() +
    theme(
      panel.grid.major.y = element_blank()
    )

  # ---------------------------------------------------------------------------
  # 6. Combine panels using patchwork
  # ---------------------------------------------------------------------------

  cat("Combining panels...\n")

  combined_plot <- p_cumulative + p_sucra +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(
      title = "Treatment Ranking Analysis",
      subtitle = "Rank 1 = lowest toxicity risk (best outcome)",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )

  # ---------------------------------------------------------------------------
  # 7. Save rank probabilities and SUCRA values to CSV
  # ---------------------------------------------------------------------------

  cat("Saving data files...\n")

  # Rank probabilities CSV
  rank_probs_df <- as.data.frame(rank_probs)
  rank_probs_df$Treatment <- rownames(rank_probs_df)
  rank_probs_df <- rank_probs_df[, c("Treatment", colnames(rank_probs))]

  rank_probs_file <- file.path(output_dir, "wt_unified_rank_probabilities.csv")
  write.csv(rank_probs_df, rank_probs_file, row.names = FALSE)
  cat("  Saved:", rank_probs_file, "\n")

  # SUCRA CSV
  sucra_output <- data.frame(
    Treatment = as.character(sucra_df$Treatment),
    Label = as.character(sucra_df$Label),
    SUCRA = sucra_df$SUCRA,
    SUCRA_pct = sucra_df$SUCRA_pct,
    stringsAsFactors = FALSE
  )
  sucra_file <- file.path(output_dir, "wt_unified_sucra.csv")
  write.csv(sucra_output, sucra_file, row.names = FALSE)
  cat("  Saved:", sucra_file, "\n")

  # ---------------------------------------------------------------------------
  # 8. Save figure
  # ---------------------------------------------------------------------------

  # Create figures directory if not exists
  figures_dir <- file.path(output_dir, "figures")
  if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
    cat("  Created directory:", figures_dir, "\n")
  }

  # Save combined plot
  figure_file <- file.path(figures_dir, "figure6_rankogram.png")
  ggsave(
    figure_file,
    plot = combined_plot,
    width = FIG_DIMS$rankogram["width"],
    height = FIG_DIMS$rankogram["height"],
    dpi = FIG_DPI
  )
  cat("  Saved:", figure_file, "\n")

  # ---------------------------------------------------------------------------
  # 9. Print confirmation and return
  # ---------------------------------------------------------------------------

  cat("\n=== Figure 6 Generation Complete ===\n")
  cat("Output files:\n")
  cat("  - ", rank_probs_file, " (rank probabilities)\n", sep = "")
  cat("  - ", sucra_file, " (SUCRA values)\n", sep = "")
  cat("  - ", figure_file, " (combined figure)\n", sep = "")

  return(combined_plot)
}

# -----------------------------------------------------------------------------
# Run if executed directly
# -----------------------------------------------------------------------------

if (!interactive() && sys.nframe() == 0) {
  plot_rankogram()
}
