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
# Table 2: Pairwise OR Matrix
# -----------------------------------------------------------------------------
# Data format in odds_ratios.csv:
#   Treatment, OR_vs_Reference
#   WT, "--"
#   HapB3, "2.0050 (1.2940, 3.1860)"
#
# NOTE: This file only contains ORs vs WT (reference).
# Full pairwise matrix would require posterior samples to calculate.
# For now, we create a table showing ORs vs WT reference.
# -----------------------------------------------------------------------------

cat("\nGenerating Table 2 (Pairwise OR Matrix)...\n")

tryCatch({
  or_data <- read.csv(file.path(output_dir, "wt_unified_odds_ratios.csv"),
                      stringsAsFactors = FALSE)

  # Parse the OR_vs_Reference string to extract OR and CrI
  # Format: "2.0050 (1.2940, 3.1860)" or "--" for reference
  parse_or <- function(or_str) {
    if (or_str == "--" || is.na(or_str)) {
      return(c(OR = 1.0, Lower = NA, Upper = NA))
    }

    # Extract OR (before parenthesis)
    or_val <- as.numeric(sub("\\s*\\(.*", "", or_str))

    # Extract CrI bounds from parentheses
    cri_match <- regmatches(or_str, regexec("\\(([0-9.]+),\\s*([0-9.]+)\\)", or_str))[[1]]
    lower_val <- as.numeric(cri_match[2])
    upper_val <- as.numeric(cri_match[3])

    return(c(OR = or_val, Lower = lower_val, Upper = upper_val))
  }

  # Parse all OR values
  parsed_ors <- t(sapply(or_data$OR_vs_Reference, parse_or))
  or_data$OR <- parsed_ors[, "OR"]
  or_data$OR_Lower <- parsed_ors[, "Lower"]
  or_data$OR_Upper <- parsed_ors[, "Upper"]

  # Create the pairwise matrix (ORs vs WT reference)
  # Since we only have ORs vs WT, this will be a column format table
  treatments <- c("WT", "HapB3", "2846hetho", "2Ahetho", "13hetho")
  n_trt <- length(treatments)

  # Initialize matrix with em-dash
  or_matrix <- matrix("---", nrow = n_trt, ncol = n_trt,
                      dimnames = list(treatments, treatments))

  # Fill diagonal with reference indicator
  diag(or_matrix) <- "Ref"

  # Fill ORs vs WT (first row = WT as row treatment, comparing to column treatments)
  # Convention: OR in row i, col j = OR of j vs i
  for (i in 1:nrow(or_data)) {
    trt <- or_data$Treatment[i]
    if (trt != "WT" && trt %in% treatments) {
      or_val <- or_data$OR[i]
      lower_val <- or_data$OR_Lower[i]
      upper_val <- or_data$OR_Upper[i]

      if (!is.na(or_val)) {
        # OR of variant vs WT
        formatted <- sprintf("%.2f (%.2f-%.2f)", or_val, lower_val, upper_val)
        or_matrix["WT", trt] <- formatted

        # Reciprocal: OR of WT vs variant
        recip_or <- 1 / or_val
        recip_upper <- 1 / lower_val
        recip_lower <- 1 / upper_val
        recip_formatted <- sprintf("%.2f (%.2f-%.2f)", recip_or, recip_lower, recip_upper)
        or_matrix[trt, "WT"] <- recip_formatted
      }
    }
  }

  # NOTE: Pairwise ORs between non-WT variants would require posterior samples
  # to calculate properly. Mark as needing posterior calculation.

  # Save as CSV
  or_df <- as.data.frame(or_matrix)
  or_df <- cbind(Treatment = rownames(or_df), or_df)
  write.csv(or_df, file.path(table_dir, "table2_pairwise_or.csv"), row.names = FALSE)
  cat("  Saved: tables/table2_pairwise_or.csv\n")
  cat("  NOTE: Only ORs vs WT reference are available from current output files.\n")
  cat("        Full pairwise ORs between variants require posterior samples.\n")

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
