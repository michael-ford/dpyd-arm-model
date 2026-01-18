# =============================================================================
# R/pairwise_analysis.R
# Calculate pairwise probabilities from MCMC samples (WT Unified Model)
# =============================================================================

cat("\n")
cat("###############################################################\n")
cat("### PAIRWISE PROBABILITY ANALYSIS                           ###\n")
cat("###############################################################\n\n")

# Source the model runner for the calculate_pairwise_probability function
source("R/common/model_runner.R")

# -----------------------------------------------------------------------------
# Load MCMC samples
# -----------------------------------------------------------------------------

cat("Loading MCMC samples...\n")

unified_samples_file <- "output/wt_unified/wt_unified_mcmc_samples.rds"

results <- list()

# -----------------------------------------------------------------------------
# Helper function for index-based pairwise comparison
# -----------------------------------------------------------------------------

calculate_pairwise_by_index <- function(mcmc_samples, idx1, idx2, ref_idx,
                                         trt1_name, trt2_name, ref_name) {
  # Combine chains
  samples_matrix <- do.call(rbind, mcmc_samples)

  cat(sprintf("\n=== Pairwise Probability: %s vs %s ===\n", trt1_name, trt2_name))
  cat(sprintf("Reference: %s (index %d)\n", ref_name, ref_idx))
  cat(sprintf("Treatment indices: %s=%d, %s=%d\n", trt1_name, idx1, trt2_name, idx2))

  # LOR columns: LOR[i,j] where i < j (pcnetmeta convention)
  get_lor_col <- function(i, j) {
    if (i < j) {
      sprintf("LOR[%d,%d]", i, j)
    } else {
      sprintf("LOR[%d,%d]", j, i)
    }
  }

  # Get LOR vs reference for each treatment
  col1 <- get_lor_col(ref_idx, idx1)
  col2 <- get_lor_col(ref_idx, idx2)

  cat(sprintf("Using columns: %s, %s\n", col1, col2))

  # Extract samples (handle sign based on index ordering)
  if (ref_idx < idx1) {
    samples1 <- samples_matrix[, col1]  # LOR[ref, trt1] = log(OR_trt1/OR_ref)
  } else {
    samples1 <- -samples_matrix[, col1]  # Negate if reversed
  }

  if (ref_idx < idx2) {
    samples2 <- samples_matrix[, col2]
  } else {
    samples2 <- -samples_matrix[, col2]
  }

  # Difference: LOR_trt1 - LOR_trt2 = log(OR_trt1/OR_trt2)
  samples_diff <- samples1 - samples2

  # Calculate probability
  p_trt1_greater <- mean(samples_diff > 0)

  # Summary statistics
  or_ratio_median <- median(exp(samples_diff))
  or_ratio_ci <- quantile(exp(samples_diff), c(0.025, 0.975))

  cat(sprintf("\nResults:\n"))
  cat(sprintf("  P(%s OR > %s OR) = %.3f\n", trt1_name, trt2_name, p_trt1_greater))
  cat(sprintf("  P(%s OR < %s OR) = %.3f\n", trt1_name, trt2_name, 1 - p_trt1_greater))
  cat(sprintf("\n  Median OR ratio (%s/%s): %.2f\n", trt1_name, trt2_name, or_ratio_median))
  cat(sprintf("  95%% CrI: [%.2f, %.2f]\n", or_ratio_ci[1], or_ratio_ci[2]))

  # Interpretation
  cat("\nInterpretation:\n")
  if (p_trt1_greater > 0.6) {
    cat(sprintf("  %s likely has higher OR (more toxic) than %s\n", trt1_name, trt2_name))
  } else if (p_trt1_greater < 0.4) {
    cat(sprintf("  %s likely has lower OR (less toxic) than %s\n", trt1_name, trt2_name))
  } else {
    cat("  Ranking is essentially a coin flip - no meaningful difference\n")
  }

  return(list(
    p_trt1_greater = p_trt1_greater,
    p_trt2_greater = 1 - p_trt1_greater,
    or_ratio_median = or_ratio_median,
    or_ratio_ci = or_ratio_ci,
    n_samples = length(samples_diff)
  ))
}

# -----------------------------------------------------------------------------
# WT Unified Model Analysis
# Treatment order: 1=WT, 2=HapB3, 3=2846hetho, 4=2Ahetho, 5=13hetho
# -----------------------------------------------------------------------------

if (file.exists(unified_samples_file)) {
  cat("\n=== WT Unified Model ===\n")
  unified_samples <- readRDS(unified_samples_file)

  # Check structure
  cat("Sample structure:\n")
  cat("  Chains:", length(unified_samples), "\n")
  cat("  Iterations per chain:", nrow(unified_samples[[1]]), "\n")
  cat("  Total samples:", length(unified_samples) * nrow(unified_samples[[1]]), "\n")

  # List available parameters
  params <- colnames(unified_samples[[1]])
  lor_params <- params[grepl("^LOR", params)]
  sigma_params <- params[grepl("^Sigma|^sigma|^sd\\.", params)]
  cat("  LOR parameters:", length(lor_params), "\n")
  cat("  Sigma/variance parameters:", length(sigma_params), "\n")
  if (length(sigma_params) > 0) {
    cat("    ", paste(sigma_params, collapse = ", "), "\n")
  }
  cat("\n")

  # Calculate pairwise probability: P(OR_2846hetho > OR_13hetho)
  # Indices: 3=2846hetho, 5=13hetho, 1=WT (reference)
  tryCatch({
    results$unified <- calculate_pairwise_by_index(
      mcmc_samples = unified_samples,
      idx1 = 3,  # 2846hetho
      idx2 = 5,  # 13hetho
      ref_idx = 1,  # WT
      trt1_name = "2846hetho",
      trt2_name = "13hetho",
      ref_name = "WT"
    )
  }, error = function(e) {
    cat("Error in unified model analysis:", e$message, "\n")
  })
} else {
  cat("WT Unified samples not found. Run analysis with mcmc.samples=TRUE first.\n")
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

cat("\n")
cat("===============================================================\n")
cat("PAIRWISE PROBABILITY SUMMARY\n")
cat("===============================================================\n\n")

cat("Question: Is the ranking difference between 2846hetho and 13hetho\n")
cat("          statistically meaningful, or just noise?\n\n")

if (length(results) > 0) {
  summary_df <- data.frame(
    Model = character(),
    P_2846_greater = numeric(),
    P_13_greater = numeric(),
    OR_ratio_median = numeric(),
    OR_ratio_CI_lower = numeric(),
    OR_ratio_CI_upper = numeric(),
    Interpretation = character(),
    stringsAsFactors = FALSE
  )

  for (model_name in names(results)) {
    r <- results[[model_name]]
    interp <- if (r$p_trt1_greater > 0.6) {
      "2846hetho likely more toxic"
    } else if (r$p_trt1_greater < 0.4) {
      "13hetho likely more toxic"
    } else {
      "No meaningful difference (coin flip)"
    }

    summary_df <- rbind(summary_df, data.frame(
      Model = model_name,
      P_2846_greater = round(r$p_trt1_greater, 3),
      P_13_greater = round(r$p_trt2_greater, 3),
      OR_ratio_median = round(r$or_ratio_median, 3),
      OR_ratio_CI_lower = round(r$or_ratio_ci[1], 3),
      OR_ratio_CI_upper = round(r$or_ratio_ci[2], 3),
      Interpretation = interp,
      stringsAsFactors = FALSE
    ))
  }

  print(summary_df)

  # Save results
  write.csv(summary_df, "output/pairwise_probability_results.csv", row.names = FALSE)
  saveRDS(results, "output/pairwise_probability_full.rds")

  cat("\nResults saved to:\n")
  cat("  - output/pairwise_probability_results.csv\n")
  cat("  - output/pairwise_probability_full.rds\n")
} else {
  cat("No results available. Ensure models were run with mcmc.samples=TRUE\n")
}

cat("\n")
