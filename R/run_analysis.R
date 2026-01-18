# =============================================================================
# run_analysis.R
# Main entry point for DPYD Arm-Based NMA (WT Unified Model)
# =============================================================================

# Set working directory (for Docker)
if (file.exists("/analysis/R")) {
  setwd("/analysis")
}

# Create output directories
for (dir in c("output", "output/wt_unified")) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
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
cat("WT Unified Model (het_cor)\n")
cat("=============================================================\n\n")

cat(sprintf("Analysis started: %s\n", Sys.time()))
cat(sprintf("Log file: %s\n\n", log_filename))

# -----------------------------------------------------------------------------
# WT Unified Model (WT_Clean and WT_Biased merged into WT)
# -----------------------------------------------------------------------------

cat("\n")
cat("###############################################################\n")
cat("### WT UNIFIED MODEL                                        ###\n")
cat("###############################################################\n")
source("R/model_wt_unified/run.R")

# -----------------------------------------------------------------------------
# Pairwise Probability Analysis
# -----------------------------------------------------------------------------

source("R/pairwise_analysis.R")

# -----------------------------------------------------------------------------
# Complete
# -----------------------------------------------------------------------------

cat("\n")
cat("=============================================================\n")
cat("Analysis Complete\n")
cat("=============================================================\n")
cat("\nWT Unified outputs:\n")
cat("  - output/wt_unified/wt_unified_result.rds\n")
cat("  - output/wt_unified/wt_unified_publication_summary.csv\n")
cat("  - output/wt_unified/wt_unified_model_summary.csv\n")
cat("  - output/wt_unified/wt_unified_absolute_risks.csv\n")
cat("  - output/wt_unified/wt_unified_odds_ratios.csv\n")
cat("\nPairwise analysis:\n")
cat("  - output/pairwise_probability_results.csv\n")
cat("\nMCMC samples (for custom analyses):\n")
cat("  - output/wt_unified/wt_unified_mcmc_samples.rds\n")
cat("\nConvergence diagnostics:\n")
cat("  - output/wt_unified/wt_unified_convergence.rds for PSRF values\n")
cat("  - output/wt_unified/wt_unified_ConvergenceDiagnostic.txt\n")
