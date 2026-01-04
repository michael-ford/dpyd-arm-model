# =============================================================================
# R/model_wt_binary/run.R
# Entry point for WT Binary model (WT_Clean and WT_Biased as separate nodes)
# =============================================================================

cat("=============================================================================\n")
cat("DPYD Arm-Based NMA: WT Binary Model\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Setup paths and source dependencies
# -----------------------------------------------------------------------------

# Use fixed paths relative to project root (works when sourced from run_analysis.R)
script_dir <- "R/model_wt_binary"
project_root <- "."

# Source configuration
source(file.path(script_dir, "config.R"))

# Source common utilities
source(file.path(project_root, "R/common/data_utils.R"))
source(file.path(project_root, "R/common/model_runner.R"))
source(file.path(project_root, "R/common/result_utils.R"))
source(file.path(project_root, "R/common/visualization.R"))

cat("Model:", MODEL_NAME, "\n")
cat("Description:", MODEL_DESCRIPTION, "\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

# -----------------------------------------------------------------------------
# Step 1: Load and Prepare Data
# -----------------------------------------------------------------------------

cat("=== Step 1: Data Preparation ===\n\n")

raw_data <- load_raw_data(DATA_FILE)
nma_data <- transform_to_nma_format(raw_data, TREATMENT_MAP, aggregate_wt = FALSE)
validate_nma_data(nma_data)
report_data_quality(nma_data, TREATMENT_NAMES)
print_treatment_summary(nma_data, TREATMENT_NAMES)
report_sparse_data(nma_data)

prepared <- save_prepared_data(nma_data, TREATMENT_NAMES, OUTPUT_DIR, MODEL_NAME)

# -----------------------------------------------------------------------------
# Step 2: Run NMA Model
# -----------------------------------------------------------------------------

cat("\n=== Step 2: Run NMA Model ===\n\n")

model_output <- run_nma_model(
  nma_data = prepared$nma_data,
  treatment_names = TREATMENT_NAMES,
  mcmc_config = get_default_mcmc_config(),
  model_name = MODEL_NAME,
  output_dir = OUTPUT_DIR,
  seed = 12345
)

result <- model_output$result
convergence_status <- model_output$convergence_status

# -----------------------------------------------------------------------------
# Step 3: Extract and Export Results
# -----------------------------------------------------------------------------

cat("\n=== Step 3: Extract Results ===\n\n")

results_summary <- extract_results(result, TREATMENT_NAMES, MODEL_NAME, OUTPUT_DIR)
export_csv_results(result, TREATMENT_NAMES, MODEL_NAME, OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Step 4: Generate Visualizations
# -----------------------------------------------------------------------------

cat("\n=== Step 4: Generate Visualizations ===\n\n")

generate_visualizations(result, prepared$nma_data, TREATMENT_NAMES, MODEL_NAME, OUTPUT_DIR)
organize_trace_plots(MODEL_NAME, OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Step 5: Publication Summary and Key Comparisons
# -----------------------------------------------------------------------------

cat("\n=== Step 5: Publication Summary ===\n\n")

generate_publication_summary(result, TREATMENT_NAMES, convergence_status, MODEL_NAME, OUTPUT_DIR)
analyze_key_comparison(result, TREATMENT_NAMES, KEY_COMPARISON_TREATMENT)

# -----------------------------------------------------------------------------
# Complete
# -----------------------------------------------------------------------------

cat("\n=============================================================================\n")
cat("WT Binary Model Analysis Complete\n")
cat("=============================================================================\n")
cat("\nOutput files saved to:", OUTPUT_DIR, "\n")
cat("Key files:\n")
cat("  -", file.path(OUTPUT_DIR, paste0(MODEL_NAME, "_result.rds")), "\n")
cat("  -", file.path(OUTPUT_DIR, paste0(MODEL_NAME, "_publication_summary.csv")), "\n")
cat("  -", file.path(OUTPUT_DIR, paste0(MODEL_NAME, "_model_summary.csv")), "\n")
