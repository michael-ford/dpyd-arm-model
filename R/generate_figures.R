# =============================================================================
# R/generate_figures.R
# Generate all publication figures for DPYD NMA
# =============================================================================

cat("\n")
cat("###############################################################\n")
cat("### FIGURE GENERATION                                       ###\n")
cat("###############################################################\n\n")

# -----------------------------------------------------------------------------
# Source individual plot functions
# -----------------------------------------------------------------------------

source("R/visualizations/constants.R")
source("R/visualizations/plot_network.R")
source("R/visualizations/plot_forest.R")
source("R/visualizations/plot_rankogram.R")

# -----------------------------------------------------------------------------
# Load required data
# -----------------------------------------------------------------------------

cat("Loading data...\n")
nma_data <- readRDS("output/wt_unified/wt_unified_nma_data.rds")
treatment_names <- readRDS("output/wt_unified/wt_unified_treatment_names.rds")
cat("  Loaded nma_data:", nrow(nma_data), "arms\n")
cat("  Treatments:", paste(treatment_names, collapse = ", "), "\n\n")

output_dir <- "output/wt_unified"

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------

fig_dir <- file.path(output_dir, "figures")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  cat("Created output directory:", fig_dir, "\n\n")
}

# -----------------------------------------------------------------------------
# Generate Figure 2: Network Diagram
# -----------------------------------------------------------------------------

cat("Generating Figure 2 (Network Diagram)...\n")
tryCatch({
  plot_network(nma_data, treatment_names, output_dir)
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
})

# -----------------------------------------------------------------------------
# Generate Figure 5: Forest Plot
# -----------------------------------------------------------------------------

cat("\nGenerating Figure 5 (Forest Plot)...\n")
tryCatch({
  plot_forest(output_dir = output_dir)
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
})

# -----------------------------------------------------------------------------
# Generate Figure 6: Rank-O-Gram
# -----------------------------------------------------------------------------

cat("\nGenerating Figure 6 (Rank-O-Gram)...\n")
tryCatch({
  plot_rankogram(output_dir = output_dir)
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
})

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

cat("\n")
cat("=== Figure Generation Complete ===\n")
cat("Output directory:", fig_dir, "\n")
cat("Generated files:\n")

generated_files <- list.files(fig_dir, pattern = "\\.png$")
if (length(generated_files) > 0) {
  for (f in generated_files) {
    file_path <- file.path(fig_dir, f)
    file_size <- file.info(file_path)$size / 1024
    cat(sprintf("  - %s (%.1f KB)\n", f, file_size))
  }
} else {
  cat("  (No PNG files found - check for errors above)\n")
}

cat("\n")
