# =============================================================================
# 03_results.R
# Extract and format NMA results, compare with contrast-based model
# =============================================================================

library(dplyr)

# -----------------------------------------------------------------------------
# Load results
# -----------------------------------------------------------------------------

result <- readRDS("output/nma_result.rds")
model_used <- readRDS("output/model_used.rds")
treatment_names <- readRDS("output/treatment_names.rds")

cat("=== Arm-Based NMA Results ===\n")
cat("Model:", model_used, "\n\n")

# -----------------------------------------------------------------------------
# Absolute Risk (Event probability per treatment)
# -----------------------------------------------------------------------------

cat("=== Absolute Risk (Toxicity Rate) ===\n")
cat("Interpretation: Probability of toxicity event per treatment\n\n")

if (!is.null(result$AbsoluteRisk)) {
  print(result$AbsoluteRisk$Median_CI)
}

# -----------------------------------------------------------------------------
# Odds Ratios vs Reference (WT_Clean)
# This is the KEY comparison - especially for 13hetho
# -----------------------------------------------------------------------------

cat("\n=== Odds Ratios vs WT_Clean (Reference) ===\n")
cat("Interpretation: OR > 1 means higher toxicity than WT_Clean\n\n")

if (!is.null(result$OddsRatio)) {
  # Full matrix
  or_matrix <- result$OddsRatio$Median_CI

  # Extract ORs vs WT_Clean (column 1)
  cat("Treatment vs WT_Clean:\n")
  for (i in 1:length(treatment_names)) {
    cat(sprintf("  %s: %s\n", treatment_names[i], or_matrix[i, 1]))
  }
}

# -----------------------------------------------------------------------------
# Log Odds Ratios (for statistical significance)
# -----------------------------------------------------------------------------

cat("\n=== Log Odds Ratios vs WT_Clean ===\n")
cat("Interpretation: 95% CrI excluding 0 = statistically significant\n\n")

if (!is.null(result$LogOddsRatio)) {
  lor_matrix <- result$LogOddsRatio$Median_CI

  for (i in 1:length(treatment_names)) {
    cat(sprintf("  %s: %s\n", treatment_names[i], lor_matrix[i, 1]))
  }
}

# -----------------------------------------------------------------------------
# Treatment Rankings
# -----------------------------------------------------------------------------

cat("\n=== Treatment Rankings ===\n")
cat("Interpretation: Probability of each rank (Rank 1 = lowest toxicity)\n\n")

if (!is.null(result$TrtRankProb)) {
  print(result$TrtRankProb)
}

# -----------------------------------------------------------------------------
# KEY COMPARISON: 13hetho in AB vs CB model
# From project-background.md:
#   - CB Direct Estimate: OR 1.14 (0.2 - 4.7)
#   - CB Indirect Estimate: OR 390
# Success: AB model OR >> 1.0, closer to indirect estimate
# -----------------------------------------------------------------------------

cat("\n")
cat("=====================================================\n")
cat("=== KEY COMPARISON: 13hetho (Arm-Based vs Contrast-Based) ===\n")
cat("=====================================================\n\n")

# Previous contrast-based results (from project-background.md)
cb_results <- data.frame(
  source = c("CB Direct", "CB Indirect"),
  OR = c(1.14, 390),
  CI_low = c(0.2, NA),
  CI_high = c(4.7, NA)
)

cat("Previous Contrast-Based Model:\n")
cat("  Direct Estimate:   OR = 1.14 (95% CI: 0.2 - 4.7)\n")
cat("  Indirect Estimate: OR ~ 390\n")
cat("  Node-split p-value: 0.000175 (INCONSISTENT)\n\n")

# Extract AB result for 13hetho (dynamic lookup)
if (!is.null(result$OddsRatio)) {
  idx_13 <- which(treatment_names == "13hetho")
  if (length(idx_13) == 0) {
    stop("Treatment '13hetho' not found in treatment_names")
  }
  ab_or_13 <- result$OddsRatio$Median_CI[idx_13, 1]
  cat("New Arm-Based Model:\n")
  cat(sprintf("  13hetho vs WT_Clean: %s\n", ab_or_13))
  cat("\n")

  # Parse the OR value and credible interval for statistical assessment
  # Format is typically "X.XX (Y.YY, Z.ZZ)" where X.XX is median, Y.YY is lower CrI, Z.ZZ is upper CrI
  # Robust regex to capture all three values
  or_pattern <- "([0-9.]+)\\s*\\(\\s*([0-9.]+)\\s*,\\s*([0-9.]+)\\s*\\)"
  or_match <- regmatches(ab_or_13, regexec(or_pattern, ab_or_13))

  if (length(or_match[[1]]) >= 4) {
    ab_or_median <- as.numeric(or_match[[1]][2])
    ab_or_lower <- as.numeric(or_match[[1]][3])
    ab_or_upper <- as.numeric(or_match[[1]][4])

    cat("Parsed Values:\n")
    cat(sprintf("  OR Median: %.4f\n", ab_or_median))
    cat(sprintf("  95%% CrI: (%.4f, %.4f)\n", ab_or_lower, ab_or_upper))
    cat("\n")

    cat("Statistical Assessment:\n")

    # Primary criterion: 95% CrI lower bound > 1.0 indicates statistical significance
    is_significant <- ab_or_lower > 1.0

    if (is_significant) {
      cat("  STATISTICALLY SIGNIFICANT: 95% CrI lower bound > 1.0\n")
      cat("  The arm-based model shows 13hetho has significantly higher toxicity than WT_Clean.\n")
      cat("  Interpretation: We can be 95% confident that OR > 1.0.\n")
    } else if (ab_or_median > 1.0 && ab_or_lower < 1.0) {
      cat("  NOT STATISTICALLY SIGNIFICANT: 95% CrI includes 1.0\n")
      cat(sprintf("  Point estimate (OR = %.2f) suggests higher toxicity, but uncertainty is high.\n", ab_or_median))
      cat("  Interpretation: Cannot rule out no difference from WT_Clean at 95% level.\n")
    } else {
      cat("  NO EVIDENCE OF INCREASED TOXICITY: OR median <= 1.0\n")
      cat("  The arm-based model does not show higher toxicity for 13hetho.\n")
    }

    # Also calculate probability that OR > 1 from the posterior if LOR is available
    if (!is.null(result$LogOddsRatio)) {
      lor_13 <- result$LogOddsRatio$Median_CI[idx_13, 1]
      lor_match <- regmatches(lor_13, regexec(or_pattern, lor_13))

      if (length(lor_match[[1]]) >= 4) {
        lor_median <- as.numeric(lor_match[[1]][2])
        lor_lower <- as.numeric(lor_match[[1]][3])
        lor_upper <- as.numeric(lor_match[[1]][4])

        # Approximate probability OR > 1 (LOR > 0) using normal approximation
        # SD approx = (upper - lower) / (2 * 1.96)
        lor_sd <- (lor_upper - lor_lower) / (2 * 1.96)
        if (lor_sd > 0) {
          prob_or_gt_1 <- pnorm(lor_median / lor_sd)
          cat(sprintf("\n  Probability OR > 1.0: %.1f%%\n", prob_or_gt_1 * 100))
        }
      }
    }

    cat("\n")
    cat("Comparison with Contrast-Based Model:\n")
    cat(sprintf("  CB Direct: OR = 1.14 | AB Model: OR = %.2f\n", ab_or_median))
    if (ab_or_median > 1.14) {
      cat("  The arm-based model gives a higher OR estimate than CB direct.\n")
    }
  } else {
    cat("  WARNING: Could not parse OR format. Raw value: ", ab_or_13, "\n")
  }
}

# -----------------------------------------------------------------------------
# Model diagnostics summary
# -----------------------------------------------------------------------------

cat("\n=== Model Diagnostics ===\n")

# Helper function to safely extract DIC values
get_dic_value <- function(dic_obj, name) {
  if (is.null(dic_obj)) return(NA)
  if (is.list(dic_obj) && name %in% names(dic_obj)) return(dic_obj[[name]])
  if (is.numeric(dic_obj) && name %in% names(dic_obj)) return(dic_obj[name])
  if (is.numeric(dic_obj) && length(dic_obj) >= 1 && name == "DIC") return(dic_obj[1])
  return(NA)
}

if (!is.null(result$DIC)) {
  dic_val <- get_dic_value(result$DIC, "DIC")
  pd_val <- get_dic_value(result$DIC, "pD")
  if (!is.na(dic_val)) cat(sprintf("DIC: %.2f\n", dic_val))
  if (!is.na(pd_val)) cat(sprintf("pD (effective parameters): %.2f\n", pd_val))
}

cat("\nConvergence: Check PSRF values in output files.\n")
cat("Target: All PSRF < 1.05\n")

# -----------------------------------------------------------------------------
# Save formatted results
# -----------------------------------------------------------------------------

results_summary <- list(
  model = model_used,
  absolute_risk = if (!is.null(result$AbsoluteRisk)) result$AbsoluteRisk$Median_CI else NULL,
  odds_ratios = if (!is.null(result$OddsRatio)) result$OddsRatio$Median_CI else NULL,
  log_odds_ratios = if (!is.null(result$LogOddsRatio)) result$LogOddsRatio$Median_CI else NULL,
  rank_probs = if (!is.null(result$TrtRankProb)) result$TrtRankProb else NULL,
  dic = if (!is.null(result$DIC)) result$DIC else NULL,
  treatment_names = treatment_names
)

saveRDS(results_summary, "output/results_summary.rds")

# -----------------------------------------------------------------------------
# Export CSV files for reproducibility
# -----------------------------------------------------------------------------

cat("\n=== Exporting CSV Files ===\n")

# 1. Odds Ratios vs Reference (WT_Clean)
if (!is.null(result$OddsRatio)) {
  or_matrix <- result$OddsRatio$Median_CI
  or_df <- data.frame(
    Treatment = treatment_names,
    OR_vs_WT_Clean = or_matrix[, 1],
    stringsAsFactors = FALSE
  )
  write.csv(or_df, "output/odds_ratios_vs_reference.csv", row.names = FALSE)
  cat("  Saved: output/odds_ratios_vs_reference.csv\n")
}

# 2. Absolute Risks (Toxicity Rate)
if (!is.null(result$AbsoluteRisk)) {
  ar_matrix <- result$AbsoluteRisk$Median_CI
  ar_df <- data.frame(
    Treatment = treatment_names,
    Toxicity_Rate = ar_matrix[, 1],
    stringsAsFactors = FALSE
  )
  write.csv(ar_df, "output/absolute_risks.csv", row.names = FALSE)
  cat("  Saved: output/absolute_risks.csv\n")
}

# 3. Treatment Rankings
if (!is.null(result$TrtRankProb)) {
  tryCatch({
    rank_mat <- result$TrtRankProb
    if (is.matrix(rank_mat) || is.array(rank_mat)) {
      rank_df <- as.data.frame(rank_mat)
      if (!is.null(rownames(rank_mat))) {
        rank_df <- cbind(Treatment = rownames(rank_mat), rank_df)
      }
      rownames(rank_df) <- NULL
      write.csv(rank_df, "output/treatment_rankings.csv", row.names = FALSE)
      cat("  Saved: output/treatment_rankings.csv\n")
    } else {
      cat("  Warning: TrtRankProb has unexpected structure, skipping CSV export\n")
    }
  }, error = function(e) {
    cat("  Warning: Could not export treatment rankings:", e$message, "\n")
  })
}

# 4. Model Summary
dic_val <- get_dic_value(result$DIC, "DIC")
pd_val <- get_dic_value(result$DIC, "pD")
model_summary <- data.frame(
  Metric = c("Model", "DIC", "pD", "Convergence_Status"),
  Value = c(
    model_used,
    if (!is.na(dic_val)) sprintf("%.2f", dic_val) else "NA",
    if (!is.na(pd_val)) sprintf("%.2f", pd_val) else "NA",
    "Check PSRF files for Gelman-Rubin < 1.05"
  ),
  stringsAsFactors = FALSE
)
write.csv(model_summary, "output/model_summary.csv", row.names = FALSE)
cat("  Saved: output/model_summary.csv\n")

# -----------------------------------------------------------------------------
# Publication-Ready Summary Table
# -----------------------------------------------------------------------------

cat("\n")
cat("=====================================================\n")
cat("=== PUBLICATION-READY SUMMARY TABLE ===\n")
cat("=====================================================\n\n")

if (!is.null(result$OddsRatio)) {
  tryCatch({
  # Parse OR values and CrI to determine significance
  or_matrix <- result$OddsRatio$Median_CI

  # Safely extract rank probabilities
  rank_probs <- NULL
  if (!is.null(result$TrtRankProb)) {
    rp <- result$TrtRankProb
    if (is.matrix(rp) && ncol(rp) >= 1) {
      rank_probs <- rp
    } else if (is.data.frame(rp) && ncol(rp) >= 1) {
      rank_probs <- as.matrix(rp)
    }
  }

  # Helper function to parse OR string and check if CrI excludes 1
  parse_or_significance <- function(or_string) {
    # Format is typically "X.XX (Y.YY, Z.ZZ)"
    match <- regmatches(or_string, regexec("([0-9.]+)\\s*\\(([0-9.]+),\\s*([0-9.]+)\\)", or_string))
    if (length(match[[1]]) >= 4) {
      ci_low <- as.numeric(match[[1]][3])
      ci_high <- as.numeric(match[[1]][4])
      # Significant if CrI excludes 1
      return(ci_low > 1 || ci_high < 1)
    }
    return(NA)
  }

  # Build summary table
  summary_table <- data.frame(
    Treatment = treatment_names,
    OR_vs_WT_Clean_95CrI = or_matrix[, 1],
    stringsAsFactors = FALSE
  )

  # Add significance column
  summary_table$Significant <- sapply(summary_table$OR_vs_WT_Clean_95CrI, function(x) {
    sig <- parse_or_significance(x)
    if (is.na(sig)) return("--")
    if (sig) return("Yes")
    return("No")
  })

  # Add P(Best) - probability of being rank 1 (lowest toxicity)
  if (!is.null(rank_probs) && is.numeric(rank_probs[, 1])) {
    summary_table$P_Best <- sprintf("%.1f%%", as.numeric(rank_probs[, 1]) * 100)
  } else {
    summary_table$P_Best <- rep("N/A", nrow(summary_table))
  }

  # Print formatted table
  cat(sprintf("%-20s %-25s %-12s %-10s\n", "Treatment", "OR vs WT_Clean (95% CrI)", "Significant", "P(Best)"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in 1:nrow(summary_table)) {
    cat(sprintf("%-20s %-25s %-12s %-10s\n",
                summary_table$Treatment[i],
                summary_table$OR_vs_WT_Clean_95CrI[i],
                summary_table$Significant[i],
                summary_table$P_Best[i]))
  }

  cat("\n")
  cat("Note: Significant = 95% Credible Interval excludes 1.0\n")
  cat("      P(Best) = Probability of having lowest toxicity (Rank 1)\n")

  # Save publication table as CSV
  write.csv(summary_table, "output/publication_summary.csv", row.names = FALSE)
  cat("\nSaved: output/publication_summary.csv\n")
  }, error = function(e) {
    cat("Warning: Could not generate publication summary:", e$message, "\n")
  })
}

cat("\n=== Results Complete ===\n")
cat("Saved: output/results_summary.rds\n")
