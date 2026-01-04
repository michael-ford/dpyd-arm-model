# =============================================================================
# run_analysis.R
# Main entry point for DPYD Arm-Based NMA
# =============================================================================

# Set working directory (for Docker)
if (file.exists("/analysis/R")) {
  setwd("/analysis")
}

# Create output directory if needed
if (!dir.exists("output")) {
  dir.create("output")
}

# -----------------------------------------------------------------------------
# Setup logging - capture all console output to timestamped log file
# -----------------------------------------------------------------------------

log_filename <- sprintf("output/analysis_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S"))

# Open connection for logging (captures both stdout and messages)
log_con <- file(log_filename, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message", append = TRUE)

# Ensure sinks are closed on exit (even if error occurs)
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  # Print final message to console (after sinks closed)
  message(sprintf("Log saved to: %s", log_filename))
}, add = TRUE)

cat("=============================================================\n")
cat("DPYD Arm-Based Network Meta-Analysis\n")
cat("Resolving Sparsity-Driven Inconsistency via AB Parameterization\n")
cat("=============================================================\n\n")

cat(sprintf("Analysis started: %s\n", Sys.time()))
cat(sprintf("Log file: %s\n\n", log_filename))

# -----------------------------------------------------------------------------
# Step 1: Prepare Data
# -----------------------------------------------------------------------------

cat("\n>>> Step 1: Preparing data...\n")
cat("-----------------------------------------------------------\n")
source("R/01_prepare_data.R")

# -----------------------------------------------------------------------------
# Step 2: Run NMA Model
# -----------------------------------------------------------------------------

cat("\n>>> Step 2: Running arm-based NMA model...\n")
cat("-----------------------------------------------------------\n")
source("R/02_run_nma.R")

# -----------------------------------------------------------------------------
# Step 3: Extract and Compare Results
# -----------------------------------------------------------------------------

cat("\n>>> Step 3: Extracting results...\n")
cat("-----------------------------------------------------------\n")
source("R/03_results.R")

# -----------------------------------------------------------------------------
# Step 4: Generate Visualizations
# -----------------------------------------------------------------------------

cat("\n>>> Step 4: Generating visualizations...\n")
cat("-----------------------------------------------------------\n")
source("R/04_visualizations.R")

# -----------------------------------------------------------------------------
# Complete
# -----------------------------------------------------------------------------

cat("\n")
cat("=============================================================\n")
cat("Analysis Complete\n")
cat("=============================================================\n")
cat("\nOutput files:\n")
cat("  - output/nma_data.rds        (prepared data)\n")
cat("  - output/nma_result.rds      (full model result)\n")
cat("  - output/results_summary.rds (formatted results)\n")
cat("\nVisualizations:\n")
cat("  - output/network_plot.png          (network structure)\n")
cat("  - output/contrast_plot.png         (ORs vs WT_Clean)\n")
cat("  - output/absolute_effects_plot.png (toxicity rates)\n")
cat("  - output/ranking_plot.png          (treatment rankings)\n")
cat("\nConvergence diagnostics:\n")
cat("  - Check PSRF*.txt files for Gelman-Rubin statistics\n")
cat("  - Check LOR*.png files for trace plots\n")
