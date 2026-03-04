# =============================================================================
# R/generate_tables.R
# Generate all supplementary tables for DPYD NMA
# =============================================================================

cat("\n")
cat("###############################################################\n")
cat("### TABLE GENERATION                                        ###\n")
cat("###############################################################\n\n")

library(dplyr)

output_dir <- "output/wt_unified"
table_dir <- file.path(output_dir, "tables")

if (!dir.exists(table_dir)) {
  dir.create(table_dir, recursive = TRUE)
  cat("Created output directory:", table_dir, "\n\n")
}

# -----------------------------------------------------------------------------
# Table S5: Absolute Toxicity Probability Estimates
# -----------------------------------------------------------------------------
# Data format in absolute_risks.csv:
#   Treatment, Toxicity_Rate
#   WT, "0.2681 (0.2158, 0.3319)"
# -----------------------------------------------------------------------------

cat("Generating Table S5 (Absolute Risks)...\n")

tryCatch({
  # Load source data
  abs_risks <- read.csv(file.path(output_dir, "wt_unified_absolute_risks.csv"),
                        stringsAsFactors = FALSE)
  nma_data <- readRDS(file.path(output_dir, "wt_unified_nma_data.rds"))
  treatment_names <- readRDS(file.path(output_dir, "wt_unified_treatment_names.rds"))

  # Parse the Toxicity_Rate string to extract median and CrI
  # Format: "0.2681 (0.2158, 0.3319)"
  parse_risk <- function(risk_str) {
    # Remove any whitespace around parentheses
    risk_str <- gsub("\\s*\\(\\s*", " (", risk_str)
    risk_str <- gsub("\\s*,\\s*", ", ", risk_str)
    risk_str <- gsub("\\s*\\)\\s*", ")", risk_str)

    # Extract median (before parenthesis)
    median_val <- as.numeric(sub("\\s*\\(.*", "", risk_str))

    # Extract CrI bounds from parentheses
    cri_match <- regmatches(risk_str, regexec("\\(([0-9.]+),\\s*([0-9.]+)\\)", risk_str))[[1]]
    lower_val <- as.numeric(cri_match[2])
    upper_val <- as.numeric(cri_match[3])

    return(c(median = median_val, lower = lower_val, upper = upper_val))
  }

  # Parse all risk values
  parsed_risks <- t(sapply(abs_risks$Toxicity_Rate, parse_risk))
  abs_risks$Risk_Median <- parsed_risks[, "median"]
  abs_risks$Risk_Lower <- parsed_risks[, "lower"]
  abs_risks$Risk_Upper <- parsed_risks[, "upper"]

  # Aggregate patient counts per treatment
  patient_counts <- nma_data %>%
    mutate(Treatment = treatment_names[t.id]) %>%
    group_by(Treatment) %>%
    summarise(
      N_Studies = n(),
      N_Patients = sum(n),
      Events = sum(r),
      .groups = "drop"
    )

  # Merge with risk estimates and format as percentages
  table_s5 <- patient_counts %>%
    left_join(abs_risks, by = "Treatment") %>%
    mutate(
      Estimated_Risk_Pct = round(Risk_Median * 100, 1),
      CrI_Lower_Pct = round(Risk_Lower * 100, 1),
      CrI_Upper_Pct = round(Risk_Upper * 100, 1)
    ) %>%
    select(
      Variant = Treatment,
      N_Studies,
      N_Patients,
      Events,
      Estimated_Risk_Pct,
      CrI_Lower_Pct,
      CrI_Upper_Pct
    )

  write.csv(table_s5, file.path(table_dir, "tableS5_absolute_risks.csv"), row.names = FALSE)
  cat("  Saved: tables/tableS5_absolute_risks.csv\n")

}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
})

# -----------------------------------------------------------------------------
# Table S6: Model Diagnostics
# -----------------------------------------------------------------------------
# Data format in model_summary.csv:
#   Metric, Value
#   Model, wt_unified
#   DIC, 151.49
#   pD, 59.45
#   ...
# -----------------------------------------------------------------------------

cat("\nGenerating Table S6 (Model Diagnostics)...\n")

tryCatch({
  model_summary <- read.csv(file.path(output_dir, "wt_unified_model_summary.csv"),
                            stringsAsFactors = FALSE)

  # Extract values from the key-value format
  get_value <- function(metric) {
    row <- model_summary[model_summary$Metric == metric, ]
    if (nrow(row) > 0) return(row$Value[1])
    return(NA)
  }

  dic_val <- get_value("DIC")
  pd_val <- get_value("pD")
  dbar_val <- get_value("D.bar")
  n_arms_val <- get_value("n_arms")
  totresdev_per_arm_val <- get_value("totresdev_per_arm")
  max_psrf_val <- get_value("Max_PSRF")

  # Format as key-value table with descriptions
  table_s6 <- data.frame(
    Metric = c(
      "DIC",
      "pD (effective parameters)",
      "Posterior mean deviance (D-bar)",
      "Total data points (arms)",
      "D-bar residual / n ratio",
      "Max Gelman-Rubin PSRF",
      "MCMC chains",
      "Iterations (post burn-in)",
      "Thinning interval"
    ),
    Value = c(
      dic_val,
      pd_val,
      dbar_val,
      n_arms_val,
      totresdev_per_arm_val,
      max_psrf_val,
      "3",
      "75000",
      "10"
    ),
    stringsAsFactors = FALSE
  )

  write.csv(table_s6, file.path(table_dir, "tableS6_diagnostics.csv"), row.names = FALSE)
  cat("  Saved: tables/tableS6_diagnostics.csv\n")

}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
})

# -----------------------------------------------------------------------------
# Table S7: Heterogeneity Parameters
# -----------------------------------------------------------------------------
# Source: /Users/mikeford/dpyd-arm-model-2/heterogeneity_params.csv
# Format: Parameter, Type, Posterior_Mean, CI_Lower, CI_Upper, Interpretation
# -----------------------------------------------------------------------------

cat("\nGenerating Table S7 (Heterogeneity)...\n")

tryCatch({
  # Search for heterogeneity params in multiple possible locations
  possible_paths <- c(
    "data/heterogeneity_params.csv",                           # Docker location
    file.path(output_dir, "wt_unified_heterogeneity_params.csv"), # Output directory
    "../heterogeneity_params.csv"                              # Parent directory (local)
  )

  het_found <- FALSE
  for (het_path in possible_paths) {
    if (file.exists(het_path)) {
      het_params <- read.csv(het_path, stringsAsFactors = FALSE)
      write.csv(het_params, file.path(table_dir, "tableS7_heterogeneity.csv"), row.names = FALSE)
      cat("  Saved: tables/tableS7_heterogeneity.csv\n")
      cat("  Source:", het_path, "\n")
      het_found <- TRUE
      break
    }
  }

  if (!het_found) {
    cat("  WARNING: heterogeneity_params.csv not found at expected locations\n")
    for (p in possible_paths) {
      cat("           Searched:", p, "\n")
    }
  }

}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
})

# -----------------------------------------------------------------------------
# Table 2 / Supplement Table S6: Full Pairwise OR Matrix
# -----------------------------------------------------------------------------
# Computes the full 5x5 pairwise OR matrix directly from MCMC posterior
# samples (LOR[i,j] columns). Convention: row = reference, col = comparator.
# -----------------------------------------------------------------------------

cat("\nGenerating Table 2 (Full Pairwise OR Matrix from posterior)...\n")

tryCatch({
  # Load MCMC posterior samples for full pairwise computation
  mcmc_file <- file.path(output_dir, "wt_unified_mcmc_samples.rds")
  treatment_names <- readRDS(file.path(output_dir, "wt_unified_treatment_names.rds"))

  if (!file.exists(mcmc_file)) {
    stop("MCMC samples not found: ", mcmc_file,
         "\nRun analysis with mcmc.samples=TRUE to generate.")
  }

  samples <- readRDS(mcmc_file)
  mat <- do.call(rbind, samples)
  n_trt <- length(treatment_names)

  # Build full pairwise OR matrix from posterior LOR samples
  # pcnetmeta convention: LOR[i,j] where i < j = log(OR of treatment i vs treatment j)
  # Cell [row, col] = OR of col vs row (row is reference, col is comparator)
  or_matrix <- matrix("Ref", nrow = n_trt, ncol = n_trt,
                       dimnames = list(treatment_names, treatment_names))

  for (i in 1:n_trt) {
    for (j in 1:n_trt) {
      if (i == j) next

      lo <- min(i, j)
      hi <- max(i, j)
      col_name <- sprintf("LOR[%d,%d]", lo, hi)
      lor_samples <- mat[, col_name]

      # LOR[lo,hi] = log(OR of lo vs hi)
      # Cell [row=i, col=j] = OR of j vs i (j as comparator, i as reference)
      # If i < j (i=lo, j=hi): OR j/i = exp(-LOR[i,j])
      # If i > j (j=lo, i=hi): OR j/i = exp(LOR[j,i])
      if (i < j) {
        or_samples <- exp(-lor_samples)
      } else {
        or_samples <- exp(lor_samples)
      }

      med <- median(or_samples)
      ci <- quantile(or_samples, c(0.025, 0.975))
      or_matrix[i, j] <- sprintf("%.2f (%.2f-%.2f)", med, ci[1], ci[2])
    }
  }

  # Save as CSV
  or_df <- as.data.frame(or_matrix)
  or_df <- cbind(Treatment = rownames(or_df), or_df)
  write.csv(or_df, file.path(table_dir, "table2_pairwise_or.csv"), row.names = FALSE)
  cat("  Saved: tables/table2_pairwise_or.csv\n")
  cat("  Full 5x5 pairwise OR matrix computed from", nrow(mat), "posterior samples.\n")

}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
})

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

cat("\n")
cat("=== Table Generation Complete ===\n")
cat("Output directory:", table_dir, "\n")
cat("Generated files:\n")

generated_files <- list.files(table_dir, pattern = "\\.csv$")
if (length(generated_files) > 0) {
  for (f in generated_files) {
    file_path <- file.path(table_dir, f)
    n_rows <- nrow(read.csv(file_path))
    cat(sprintf("  - %s (%d rows)\n", f, n_rows))
  }
} else {
  cat("  (No CSV files found - check for errors above)\n")
}

cat("\n")
