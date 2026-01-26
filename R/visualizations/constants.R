# =============================================================================
# DPYD NMA Visualization Constants
# =============================================================================
# This file defines shared constants for all visualization scripts in the
# DPYD pharmacogenetics network meta-analysis project.
# =============================================================================

library(ggplot2)

# -----------------------------------------------------------------------------
# Color Palette for DPYD Variants
# -----------------------------------------------------------------------------
# Colors are chosen to reflect functional classification:
#   - Green: Normal Function (Wild-type)
#   - Blue/Orange: Decreased Function (HapB3, c.2846A>T)
#   - Red: No Function (*2A, *13)

VARIANT_COLORS <- c(
  "WT"        = "#2E7D32",
  "HapB3"     = "#1976D2",
  "2846hetho" = "#F57C00",
  "2Ahetho"   = "#C62828",
  "13hetho"   = "#8E24AA"
)

# -----------------------------------------------------------------------------
# Treatment Labels for Display
# -----------------------------------------------------------------------------
# Maps internal data names to human-readable labels for figures and tables

TREATMENT_LABELS <- c(
  "WT"        = "Wild-type",
  "HapB3"     = "HapB3",
  "2846hetho" = "c.2846A>T",
  "2Ahetho"   = "*2A",
  "13hetho"   = "*13"
)

# -----------------------------------------------------------------------------
# Treatment Display Order
# -----------------------------------------------------------------------------
# Preferred order for displaying treatments in figures (reference first,
# then by functional classification: Decreased Function, No Function)

TREATMENT_ORDER <- c("WT", "HapB3", "2846hetho", "2Ahetho", "13hetho")

# -----------------------------------------------------------------------------
# Figure Dimensions (in inches)
# -----------------------------------------------------------------------------

FIG_DIMS <- list(
  network   = c(width = 6, height = 6),
  forest    = c(width = 10, height = 5),
  rankogram = c(width = 10, height = 5)
)

# -----------------------------------------------------------------------------
# Figure Resolution
# -----------------------------------------------------------------------------

FIG_DPI <- 300

# -----------------------------------------------------------------------------
# Publication-Ready ggplot2 Theme
# -----------------------------------------------------------------------------
#' Create a publication-ready ggplot2 theme
#'
#' @param base_size Base font size (default: 12)
#' @return A ggplot2 theme object
#' @export
#'
#' @examples
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   theme_publication()

theme_publication <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      # Title styling
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = rel(1.2)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = rel(1.0)
      ),

      # Legend styling
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),

      # Axis styling
      axis.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),

      # Grid styling
      panel.grid.minor = element_blank(),

      # Strip styling (for faceted plots)
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA),

      # Plot margins
      plot.margin = margin(15, 15, 15, 15)
    )
}
