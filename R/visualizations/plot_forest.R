# =============================================================================
# Figure 5: Forest Plot of Odds Ratios vs Wild-Type
# =============================================================================
# Generates a publication-ready forest plot showing odds ratios and 95%
# credible intervals for each DPYD variant compared to wild-type reference.
# =============================================================================

# Source shared constants (colors, theme, etc.) if not already loaded
if (!exists("VARIANT_COLORS")) {
  source("R/visualizations/constants.R")
}

# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)

# -----------------------------------------------------------------------------
#' Generate Forest Plot of Odds Ratios vs Wild-Type
#'
#' Creates a horizontal forest plot displaying point estimates and 95% credible
#' intervals for each DPYD variant's odds ratio compared to wild-type.
#'
#' @param or_csv_path Path to CSV file containing OR data.
#'   Default: "output/wt_unified/wt_unified_odds_ratios.csv"
#' @param output_dir Directory where figure will be saved.
#'   Default: "output/wt_unified"
#'
#' @return A ggplot2 object (invisibly)
#'
#' @examples
#' plot_forest()
#' plot_forest(or_csv_path = "output/custom/odds_ratios.csv",
#'             output_dir = "output/custom")
# -----------------------------------------------------------------------------

plot_forest <- function(or_csv_path = "output/wt_unified/wt_unified_odds_ratios.csv",
                        output_dir = "output/wt_unified") {

  cat("\n=== Generating Figure 5: Forest Plot ===\n\n")

  # ---------------------------------------------------------------------------
  # 1. Load and parse the OR data
  # ---------------------------------------------------------------------------

  cat("Loading OR data from:", or_csv_path, "\n")
  or_raw <- read.csv(or_csv_path, stringsAsFactors = FALSE)

  # Filter out WT (which has "--" as OR_vs_Reference)
  or_filtered <- or_raw %>%
    filter(OR_vs_Reference != "--")

  cat("  Variants found:", nrow(or_filtered), "\n")

  # ---------------------------------------------------------------------------
  # 2. Parse the OR_vs_Reference string to extract OR, lower, upper
  # ---------------------------------------------------------------------------
  # Format: "2.0050 (1.2940, 3.1860)"
  # Regex explanation:
  #   ^([0-9.]+)         - OR value (first number at start)
  #   \s*\(\s*           - opening paren with optional whitespace

  #   ([0-9.]+)          - lower bound
  #   \s*,\s*            - comma with optional whitespace
  #   ([0-9.]+)          - upper bound
  #   \s*\)$             - closing paren at end

  parsed_data <- or_filtered %>%
    mutate(
      # Extract the OR (first number before the parenthesis)
      or = as.numeric(str_extract(OR_vs_Reference, "^[0-9.]+")),
      # Extract lower bound (number after opening paren)
      lower = as.numeric(str_match(OR_vs_Reference, "\\(([0-9.]+),")[, 2]),
      # Extract upper bound (number before closing paren)
      upper = as.numeric(str_match(OR_vs_Reference, ",\\s*([0-9.]+)\\)")[, 2])
    ) %>%
    rename(variant = Treatment) %>%
    select(variant, or, lower, upper)

  cat("  Parsed OR data:\n")
  print(parsed_data)

  # ---------------------------------------------------------------------------
  # 3. Order variants by OR magnitude (ascending)
  # ---------------------------------------------------------------------------

  parsed_data <- parsed_data %>%
    arrange(or) %>%
    mutate(variant = factor(variant, levels = variant))

  # ---------------------------------------------------------------------------
  # 4. Create formatted label for OR (95% CrI)
  # ---------------------------------------------------------------------------

  parsed_data <- parsed_data %>%
    mutate(
      or_label = sprintf("%.2f (%.2f, %.2f)", or, lower, upper)
    )

  # ---------------------------------------------------------------------------
  # 5. Create the forest plot
  # ---------------------------------------------------------------------------

  cat("\nCreating forest plot...\n")

  # Calculate x position for text labels (right side of plot)
  label_x <- 30

  p <- ggplot(parsed_data, aes(x = or, y = variant, color = variant)) +
    # Vertical reference line at OR = 1 (null effect)
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +

    # Confidence interval error bars
    geom_errorbarh(
      aes(xmin = lower, xmax = upper),
      height = 0.2,
      linewidth = 0.8
    ) +

    # Point estimates
    geom_point(size = 3.5) +

    # OR labels on right side
    geom_text(
      aes(x = label_x, label = or_label),
      hjust = 0,
      size = 3.5,
      color = "black",
      fontface = "plain"
    ) +

    # Log scale for x-axis
    scale_x_log10(
      breaks = c(0.5, 1, 2, 5, 10, 20),
      limits = c(0.5, 35)
    ) +

    # Use variant colors from constants
    scale_color_manual(values = VARIANT_COLORS, guide = "none") +

    # Labels
    labs(
      title = "Odds Ratios for Severe Toxicity vs Wild-Type",
      x = "Odds Ratio (log scale)",
      y = NULL,
      caption = "OR (95% CrI)"
    ) +

    # Apply publication theme
    theme_publication() +

    # Additional theme customizations for forest plot
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(face = "bold", size = 11),
      plot.caption = element_text(hjust = 1, face = "italic", size = 9)
    )

  # ---------------------------------------------------------------------------
  # 6. Create figures directory and save
  # ---------------------------------------------------------------------------

  figures_dir <- file.path(output_dir, "figures")
  if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
    cat("  Created directory:", figures_dir, "\n")
  }

  output_path <- file.path(figures_dir, "figure5_forest.png")

  ggsave(
    filename = output_path,
    plot = p,
    width = FIG_DIMS$forest["width"],
    height = FIG_DIMS$forest["height"],
    dpi = FIG_DPI,
    bg = "white"
  )

  cat("\n=== Figure 5 Complete ===\n")
  cat("Saved to:", output_path, "\n")

  # Return the plot object invisibly
  invisible(p)
}
