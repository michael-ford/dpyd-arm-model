# =============================================================================
# R/common/result_utils.R
# Shared result extraction and formatting utilities
# =============================================================================

library(dplyr)

# -----------------------------------------------------------------------------
# Extract and print results summary
# -----------------------------------------------------------------------------

extract_results <- function(result, treatment_names, model_name, output_dir = "output") {
  cat("=== Arm-Based NMA Results ===\n")
  cat("Model:", model_name, "\n\n")

  # Absolute Risk
  cat("=== Absolute Risk (Toxicity Rate) ===\n")
  cat("Interpretation: Probability of toxicity event per treatment\n\n")
  if (!is.null(result$AbsoluteRisk)) {
    print(result$AbsoluteRisk$Median_CI)
  }

  # Odds Ratios vs Reference
  reference_name <- treatment_names[1]
  cat("\n=== Odds Ratios vs", reference_name, "(Reference) ===\n")
  cat("Interpretation: OR > 1 means higher toxicity than", reference_name, "\n\n")

  if (!is.null(result$OddsRatio)) {
    or_matrix <- result$OddsRatio$Median_CI
    cat("Treatment vs", reference_name, ":\n")
    for (i in 1:length(treatment_names)) {
      cat(sprintf("  %s: %s\n", treatment_names[i], or_matrix[i, 1]))
    }
  }

  # Log Odds Ratios
  cat("\n=== Log Odds Ratios vs", reference_name, "===\n")
  cat("Interpretation: 95% CrI excluding 0 = statistically significant\n\n")

  if (!is.null(result$LogOddsRatio)) {
    lor_matrix <- result$LogOddsRatio$Median_CI
    for (i in 1:length(treatment_names)) {
      cat(sprintf("  %s: %s\n", treatment_names[i], lor_matrix[i, 1]))
    }
  }

  # Treatment Rankings
  cat("\n=== Treatment Rankings ===\n")
  cat("Interpretation: Probability of each rank (Rank 1 = lowest toxicity)\n\n")
  if (!is.null(result$TrtRankProb)) {
    print(result$TrtRankProb)
  }

  # Build results summary
  results_summary <- list(
    model_name = model_name,
    treatment_names = treatment_names,
    absolute_risk = if (!is.null(result$AbsoluteRisk)) result$AbsoluteRisk$Median_CI else NULL,
    odds_ratios = if (!is.null(result$OddsRatio)) result$OddsRatio$Median_CI else NULL,
    log_odds_ratios = if (!is.null(result$LogOddsRatio)) result$LogOddsRatio$Median_CI else NULL,
    rank_probs = if (!is.null(result$TrtRankProb)) result$TrtRankProb else NULL,
    dic = if (!is.null(result$DIC)) result$DIC else NULL
  )

  # Save results summary
  saveRDS(results_summary, file.path(output_dir, paste0(model_name, "_results_summary.rds")))

  return(results_summary)
}

# -----------------------------------------------------------------------------
# Export CSV files
# -----------------------------------------------------------------------------

export_csv_results <- function(result, treatment_names, model_name, output_dir = "output") {
  cat("\n=== Exporting CSV Files ===\n")

  reference_name <- treatment_names[1]

  # Odds Ratios vs Reference
  if (!is.null(result$OddsRatio)) {
    or_matrix <- result$OddsRatio$Median_CI
    or_df <- data.frame(
      Treatment = treatment_names,
      OR_vs_Reference = or_matrix[, 1],
      stringsAsFactors = FALSE
    )
    filename <- file.path(output_dir, paste0(model_name, "_odds_ratios.csv"))
    write.csv(or_df, filename, row.names = FALSE)
    cat("  Saved:", filename, "\n")
  }

  # Absolute Risks
  if (!is.null(result$AbsoluteRisk)) {
    ar_matrix <- result$AbsoluteRisk$Median_CI
    ar_df <- data.frame(
      Treatment = treatment_names,
      Toxicity_Rate = ar_matrix[, 1],
      stringsAsFactors = FALSE
    )
    filename <- file.path(output_dir, paste0(model_name, "_absolute_risks.csv"))
    write.csv(ar_df, filename, row.names = FALSE)
    cat("  Saved:", filename, "\n")
  }

  # Treatment Rankings
  if (!is.null(result$TrtRankProb)) {
    tryCatch({
      rank_mat <- result$TrtRankProb
      if (is.matrix(rank_mat) || is.array(rank_mat)) {
        rank_df <- as.data.frame(rank_mat)
        if (!is.null(rownames(rank_mat))) {
          rank_df <- cbind(Treatment = rownames(rank_mat), rank_df)
        }
        rownames(rank_df) <- NULL
        filename <- file.path(output_dir, paste0(model_name, "_rankings.csv"))
        write.csv(rank_df, filename, row.names = FALSE)
        cat("  Saved:", filename, "\n")
      }
    }, error = function(e) {
      cat("  Warning: Could not export treatment rankings:", e$message, "\n")
    })
  }
}

# -----------------------------------------------------------------------------
# Generate publication-ready summary table
# -----------------------------------------------------------------------------

generate_publication_summary <- function(result, treatment_names, convergence_status,
                                          model_name, output_dir = "output") {
  cat("\n")
  cat("=====================================================\n")
  cat("=== PUBLICATION-READY SUMMARY TABLE ===\n")
  cat("=====================================================\n\n")

  reference_name <- treatment_names[1]

  if (is.null(result$OddsRatio)) {
    cat("No odds ratio results available.\n")
    return(NULL)
  }

  tryCatch({
    or_matrix <- result$OddsRatio$Median_CI

    # Helper to parse OR and check significance
    parse_or_significance <- function(or_string) {
      match <- regmatches(or_string, regexec("([0-9.]+)\\s*\\(([0-9.]+),\\s*([0-9.]+)\\)", or_string))
      if (length(match[[1]]) >= 4) {
        ci_low <- as.numeric(match[[1]][3])
        ci_high <- as.numeric(match[[1]][4])
        return(ci_low > 1 || ci_high < 1)
      }
      return(NA)
    }

    # Build summary table
    summary_table <- data.frame(
      Treatment = treatment_names,
      OR_vs_Reference_95CrI = or_matrix[, 1],
      stringsAsFactors = FALSE
    )

    # Add significance column
    summary_table$Significant <- sapply(summary_table$OR_vs_Reference_95CrI, function(x) {
      sig <- parse_or_significance(x)
      if (is.na(sig)) return("--")
      if (sig) return("Yes")
      return("No")
    })

    # Add P(Best) if available
    if (!is.null(result$TrtRankProb)) {
      rp <- result$TrtRankProb
      if (is.matrix(rp) && ncol(rp) >= 1 && is.numeric(rp[, 1])) {
        summary_table$P_Best <- sprintf("%.1f%%", as.numeric(rp[, 1]) * 100)
      } else {
        summary_table$P_Best <- rep("N/A", nrow(summary_table))
      }
    } else {
      summary_table$P_Best <- rep("N/A", nrow(summary_table))
    }

    # Print formatted table
    cat(sprintf("%-20s %-25s %-12s %-10s\n",
                "Treatment", paste0("OR vs ", reference_name, " (95% CrI)"), "Significant", "P(Best)"))
    cat(paste(rep("-", 70), collapse = ""), "\n")

    for (i in 1:nrow(summary_table)) {
      cat(sprintf("%-20s %-25s %-12s %-10s\n",
                  summary_table$Treatment[i],
                  summary_table$OR_vs_Reference_95CrI[i],
                  summary_table$Significant[i],
                  summary_table$P_Best[i]))
    }

    cat("\n")
    cat("Note: Significant = 95% Credible Interval excludes 1.0\n")
    cat("      P(Best) = Probability of having lowest toxicity (Rank 1)\n")

    # Add model diagnostics
    cat("\nModel Diagnostics:\n")
    if (!is.na(convergence_status$dic_value)) {
      cat(sprintf("  DIC: %.2f\n", convergence_status$dic_value))
    }
    if (!is.na(convergence_status$pd_value)) {
      cat(sprintf("  pD: %.2f\n", convergence_status$pd_value))
    }
    if (!is.na(convergence_status$max_psrf)) {
      cat(sprintf("  Max PSRF: %.4f (threshold: 1.05)\n", convergence_status$max_psrf))
    }
    if (!is.na(convergence_status$converged)) {
      cat(sprintf("  Convergence: %s\n", ifelse(convergence_status$converged, "YES", "NO")))
    }

    # Save publication summary
    filename <- file.path(output_dir, paste0(model_name, "_publication_summary.csv"))
    write.csv(summary_table, filename, row.names = FALSE)
    cat("\nSaved:", filename, "\n")

    # Save model summary with fit statistics
    dres_value <- if (!is.null(convergence_status$dres_value) && !is.na(convergence_status$dres_value)) {
      sprintf("%.2f", convergence_status$dres_value)
    } else { "NA" }

    n_arms_value <- if (!is.null(convergence_status$n_arms) && !is.na(convergence_status$n_arms)) {
      as.character(convergence_status$n_arms)
    } else { "NA" }

    dres_ratio <- if (!is.null(convergence_status$dres_value) && !is.na(convergence_status$dres_value) &&
                      !is.null(convergence_status$n_arms) && !is.na(convergence_status$n_arms)) {
      sprintf("%.2f", convergence_status$dres_value / convergence_status$n_arms)
    } else { "NA" }

    # D.bar = deviance at posterior mean D(θ̄)
    # totresdev = posterior mean of residual deviance (from MCMC samples)
    # For absolute fit: totresdev/n_arms should be ≈ 1
    dbar_value <- if (!is.null(convergence_status$dbar_value) && !is.na(convergence_status$dbar_value)) {
      sprintf("%.2f", convergence_status$dbar_value)
    } else { "NA" }

    model_summary <- data.frame(
      Metric = c("Model", "DIC", "pD", "D.bar", "totresdev", "n_arms", "totresdev_per_arm", "Max_PSRF", "Converged"),
      Value = c(
        model_name,
        ifelse(!is.na(convergence_status$dic_value), sprintf("%.2f", convergence_status$dic_value), "NA"),
        ifelse(!is.na(convergence_status$pd_value), sprintf("%.2f", convergence_status$pd_value), "NA"),
        dbar_value,
        dres_value,
        n_arms_value,
        dres_ratio,
        ifelse(!is.na(convergence_status$max_psrf), sprintf("%.4f", convergence_status$max_psrf), "NA"),
        ifelse(!is.na(convergence_status$converged), as.character(convergence_status$converged), "NA")
      ),
      stringsAsFactors = FALSE
    )
    filename <- file.path(output_dir, paste0(model_name, "_model_summary.csv"))
    write.csv(model_summary, filename, row.names = FALSE)
    cat("Saved:", filename, "\n")

    return(summary_table)
  }, error = function(e) {
    cat("Warning: Could not generate publication summary:", e$message, "\n")
    return(NULL)
  })
}

# -----------------------------------------------------------------------------
# Analyze key comparison (e.g., 13hetho)
# -----------------------------------------------------------------------------

analyze_key_comparison <- function(result, treatment_names, target_treatment = "13hetho",
                                    reference_treatment = NULL) {
  if (is.null(reference_treatment)) {
    reference_treatment <- treatment_names[1]
  }

  cat("\n")
  cat("=====================================================\n")
  cat("=== KEY COMPARISON:", target_treatment, "vs", reference_treatment, "===\n")
  cat("=====================================================\n\n")

  idx_target <- which(treatment_names == target_treatment)
  if (length(idx_target) == 0) {
    cat("Treatment '", target_treatment, "' not found in treatment_names\n", sep = "")
    return(NULL)
  }

  if (!is.null(result$OddsRatio)) {
    or_string <- result$OddsRatio$Median_CI[idx_target, 1]
    cat(target_treatment, "vs", reference_treatment, ":", or_string, "\n\n")

    # Parse OR value
    or_pattern <- "([0-9.]+)\\s*\\(\\s*([0-9.]+)\\s*,\\s*([0-9.]+)\\s*\\)"
    or_match <- regmatches(or_string, regexec(or_pattern, or_string))

    if (length(or_match[[1]]) >= 4) {
      or_median <- as.numeric(or_match[[1]][2])
      or_lower <- as.numeric(or_match[[1]][3])
      or_upper <- as.numeric(or_match[[1]][4])

      cat("Parsed Values:\n")
      cat(sprintf("  OR Median: %.4f\n", or_median))
      cat(sprintf("  95%% CrI: (%.4f, %.4f)\n", or_lower, or_upper))

      cat("\nStatistical Assessment:\n")
      if (or_lower > 1.0) {
        cat("  STATISTICALLY SIGNIFICANT: 95% CrI lower bound > 1.0\n")
        cat("  Interpretation: ", target_treatment, " has significantly higher toxicity than ", reference_treatment, "\n", sep = "")
      } else if (or_median > 1.0 && or_lower < 1.0) {
        cat("  NOT STATISTICALLY SIGNIFICANT: 95% CrI includes 1.0\n")
        cat(sprintf("  Point estimate (OR = %.2f) suggests higher toxicity, but uncertainty is high.\n", or_median))
      } else {
        cat("  NO EVIDENCE OF INCREASED TOXICITY: OR median <= 1.0\n")
      }

      return(list(
        treatment = target_treatment,
        or_median = or_median,
        or_lower = or_lower,
        or_upper = or_upper,
        significant = or_lower > 1.0
      ))
    }
  }

  return(NULL)
}
