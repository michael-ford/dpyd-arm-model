# =============================================================================
# run_analysis.R
# Main entry point for DPYD Arm-Based NMA (Multi-Model)
# =============================================================================

# Set working directory (for Docker)
if (file.exists("/analysis/R")) {
  setwd("/analysis")
}

# Create output directories
for (dir in c("output", "output/wt_binary", "output/wt_unified")) {
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
cat("Multi-Model Comparison: WT Binary vs WT Unified\n")
cat("=============================================================\n\n")

cat(sprintf("Analysis started: %s\n", Sys.time()))
cat(sprintf("Log file: %s\n\n", log_filename))

# -----------------------------------------------------------------------------
# Model 1: WT Binary (WT_Clean and WT_Biased as separate nodes)
# -----------------------------------------------------------------------------

cat("\n")
cat("###############################################################\n")
cat("### MODEL 1: WT BINARY                                      ###\n")
cat("###############################################################\n")
source("R/model_wt_binary/run.R")

# -----------------------------------------------------------------------------
# Model 2: WT Unified (WT_Clean and WT_Biased merged into WT)
# -----------------------------------------------------------------------------

cat("\n")
cat("###############################################################\n")
cat("### MODEL 2: WT UNIFIED                                     ###\n")
cat("###############################################################\n")
source("R/model_wt_unified/run.R")

# -----------------------------------------------------------------------------
# Compare Models
# -----------------------------------------------------------------------------

cat("\n")
cat("###############################################################\n")
cat("### MODEL COMPARISON                                        ###\n")
cat("###############################################################\n")
source("R/compare_models.R")

# -----------------------------------------------------------------------------
# Pairwise Probability Analysis
# -----------------------------------------------------------------------------

source("R/pairwise_analysis.R")

# -----------------------------------------------------------------------------
# Generate Publication Figures
# -----------------------------------------------------------------------------

source("R/generate_figures.R")

# -----------------------------------------------------------------------------
# Generate Supplementary Tables
# -----------------------------------------------------------------------------

source("R/generate_tables.R")

# -----------------------------------------------------------------------------
# Complete
# -----------------------------------------------------------------------------

cat("\n")
cat("=============================================================\n")
cat("Multi-Model Analysis Complete\n")
cat("=============================================================\n")
cat("\nModel 1 (WT Binary) outputs:\n")
cat("  - output/wt_binary/wt_binary_result.rds\n")
cat("  - output/wt_binary/wt_binary_publication_summary.csv\n")
cat("  - output/wt_binary/wt_binary_model_summary.csv\n")
cat("\nModel 2 (WT Unified) outputs:\n")
cat("  - output/wt_unified/wt_unified_result.rds\n")
cat("  - output/wt_unified/wt_unified_publication_summary.csv\n")
cat("  - output/wt_unified/wt_unified_model_summary.csv\n")
cat("\nComparison outputs:\n")
cat("  - output/model_comparison_dic.csv\n")
cat("  - output/pairwise_probability_results.csv\n")
cat("\nMCMC samples (for custom analyses):\n")
cat("  - output/wt_binary/wt_binary_mcmc_samples.rds\n")
cat("  - output/wt_unified/wt_unified_mcmc_samples.rds\n")
cat("\nConvergence diagnostics:\n")
cat("  - Check output/*/wt_*_convergence.rds for PSRF values\n")
cat("  - Check output/*/ for trace plots\n")
cat("\nPublication figures:\n")
cat("  - output/wt_unified/figures/figure2_network.png\n")
cat("  - output/wt_unified/figures/figure5_forest.png\n")
cat("  - output/wt_unified/figures/figure6_rankogram.png\n")
cat("\nSupplementary tables:\n")
cat("  - output/wt_unified/tables/tableS5_absolute_risks.csv\n")
cat("  - output/wt_unified/tables/tableS6_diagnostics.csv\n")
cat("  - output/wt_unified/tables/tableS7_heterogeneity.csv\n")
cat("  - output/wt_unified/tables/table2_pairwise_or.csv\n")
