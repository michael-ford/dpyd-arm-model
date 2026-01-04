# =============================================================================
# R/compare_models.R
# Compare results across NMA model variants
# =============================================================================

library(dplyr)

cat("=============================================================================\n")
cat("DPYD Arm-Based NMA: Model Comparison\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

models <- list(
  wt_binary = list(
    name = "WT Binary",
    output_dir = "output/wt_binary",
    prefix = "wt_binary"
  ),
  wt_unified = list(
    name = "WT Unified",
    output_dir = "output/wt_unified",
    prefix = "wt_unified"
  )
)

# Treatments present in both models (for comparison)
COMMON_TREATMENTS <- c("HapB3", "2846hetho", "2Ahetho", "13hetho")
KEY_TREATMENT <- "13hetho"

# -----------------------------------------------------------------------------
# Load model results
# -----------------------------------------------------------------------------

load_model_results <- function(model_config) {
  result_file <- file.path(model_config$output_dir, paste0(model_config$prefix, "_result.rds"))
  convergence_file <- file.path(model_config$output_dir, paste0(model_config$prefix, "_convergence.rds"))
  treatments_file <- file.path(model_config$output_dir, paste0(model_config$prefix, "_treatment_names.rds"))

  if (!file.exists(result_file)) {
    cat("  Warning: Result file not found:", result_file, "\n")
    return(NULL)
  }

  list(
    result = readRDS(result_file),
    convergence = if (file.exists(convergence_file)) readRDS(convergence_file) else NULL,
    treatments = if (file.exists(treatments_file)) readRDS(treatments_file) else NULL
  )
}

cat("Loading model results...\n\n")

model_data <- lapply(models, function(m) {
  cat("  Loading", m$name, "from", m$output_dir, "...\n")
  load_model_results(m)
})

# Check which models loaded successfully
loaded_models <- names(model_data)[!sapply(model_data, is.null)]
cat("\nLoaded models:", paste(loaded_models, collapse = ", "), "\n")

if (length(loaded_models) < 2) {
  cat("\nNeed at least 2 models to compare. Run both models first.\n")
  cat("Usage:\n")
  cat("  Rscript R/model_wt_binary/run.R\n")
  cat("  Rscript R/model_wt_unified/run.R\n")
  cat("  Rscript R/compare_models.R\n")
  quit(status = 1)
}

# -----------------------------------------------------------------------------
# Compare Model Fit (DIC)
# -----------------------------------------------------------------------------

cat("\n")
cat("=====================================================\n")
cat("=== MODEL FIT COMPARISON (DIC) ===\n")
cat("=====================================================\n\n")

dic_comparison <- data.frame(
  Model = character(),
  DIC = numeric(),
  pD = numeric(),
  Converged = character(),
  Max_PSRF = numeric(),
  stringsAsFactors = FALSE
)

for (model_name in loaded_models) {
  m <- models[[model_name]]
  d <- model_data[[model_name]]

  dic_val <- NA
  pd_val <- NA
  converged <- "Unknown"
  max_psrf <- NA

  if (!is.null(d$convergence)) {
    dic_val <- d$convergence$dic_value
    pd_val <- d$convergence$pd_value
    converged <- ifelse(!is.na(d$convergence$converged),
                        ifelse(d$convergence$converged, "Yes", "No"),
                        "Unknown")
    max_psrf <- d$convergence$max_psrf
  }

  dic_comparison <- rbind(dic_comparison, data.frame(
    Model = m$name,
    DIC = dic_val,
    pD = pd_val,
    Converged = converged,
    Max_PSRF = max_psrf,
    stringsAsFactors = FALSE
  ))
}

print(dic_comparison)

# Identify best model by DIC
if (all(!is.na(dic_comparison$DIC))) {
  best_idx <- which.min(dic_comparison$DIC)
  dic_diff <- abs(diff(dic_comparison$DIC))
  cat("\n")
  cat("Best model by DIC:", dic_comparison$Model[best_idx], "\n")
  cat("DIC difference:", round(dic_diff, 2), "\n")
  if (dic_diff < 2) {
    cat("Interpretation: Difference < 2 suggests models are essentially equivalent\n")
  } else if (dic_diff < 5) {
    cat("Interpretation: Difference 2-5 suggests weak preference for lower DIC model\n")
  } else {
    cat("Interpretation: Difference > 5 suggests meaningful improvement with lower DIC model\n")
  }
}

# -----------------------------------------------------------------------------
# Compare Odds Ratios for Key Treatment (13hetho)
# -----------------------------------------------------------------------------

cat("\n")
cat("=====================================================\n")
cat("=== KEY TREATMENT COMPARISON:", KEY_TREATMENT, "===\n")
cat("=====================================================\n\n")

or_comparison <- data.frame(
  Model = character(),
  Reference = character(),
  OR_Median_95CrI = character(),
  stringsAsFactors = FALSE
)

for (model_name in loaded_models) {
  m <- models[[model_name]]
  d <- model_data[[model_name]]

  if (is.null(d$treatments)) next

  ref_name <- d$treatments[1]
  key_idx <- which(d$treatments == KEY_TREATMENT)

  if (length(key_idx) == 0) {
    cat("  ", KEY_TREATMENT, "not found in", m$name, "\n")
    next
  }

  or_string <- "N/A"
  if (!is.null(d$result$OddsRatio)) {
    or_matrix <- d$result$OddsRatio$Median_CI
    or_string <- or_matrix[key_idx, 1]
  }

  or_comparison <- rbind(or_comparison, data.frame(
    Model = m$name,
    Reference = ref_name,
    OR_Median_95CrI = or_string,
    stringsAsFactors = FALSE
  ))
}

cat(sprintf("%-15s %-15s %-30s\n", "Model", "Reference", paste0("OR (", KEY_TREATMENT, " vs Ref)")))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in seq_len(nrow(or_comparison))) {
  cat(sprintf("%-15s %-15s %-30s\n",
              or_comparison$Model[i],
              or_comparison$Reference[i],
              or_comparison$OR_Median_95CrI[i]))
}

cat("\nNote: WT Binary uses WT_Clean as reference\n")
cat("      WT Unified merges WT_Clean + WT_Biased into single WT reference\n")

# -----------------------------------------------------------------------------
# Compare All Common Treatment ORs
# -----------------------------------------------------------------------------

cat("\n")
cat("=====================================================\n")
cat("=== ALL TREATMENT ODDS RATIOS ===\n")
cat("=====================================================\n\n")

for (trt in COMMON_TREATMENTS) {
  cat(trt, "vs Reference:\n")
  for (model_name in loaded_models) {
    m <- models[[model_name]]
    d <- model_data[[model_name]]

    if (is.null(d$treatments)) next

    trt_idx <- which(d$treatments == trt)
    if (length(trt_idx) == 0) next

    or_string <- "N/A"
    if (!is.null(d$result$OddsRatio)) {
      or_matrix <- d$result$OddsRatio$Median_CI
      or_string <- or_matrix[trt_idx, 1]
    }

    cat(sprintf("  %-12s: %s\n", m$name, or_string))
  }
  cat("\n")
}

# -----------------------------------------------------------------------------
# Compare Treatment Rankings (Probability of Best)
# -----------------------------------------------------------------------------

cat("=====================================================\n")
cat("=== TREATMENT RANKINGS (P(Best) = Lowest Toxicity) ===\n")
cat("=====================================================\n\n")

for (model_name in loaded_models) {
  m <- models[[model_name]]
  d <- model_data[[model_name]]

  cat(m$name, ":\n")

  tryCatch({
    if (!is.null(d$result$TrtRankProb)) {
      rp <- d$result$TrtRankProb
      if (is.matrix(rp) && ncol(rp) >= 1) {
        p_best <- as.numeric(rp[, 1])
        if (all(!is.na(p_best))) {
          order_idx <- order(p_best, decreasing = TRUE)
          for (i in order_idx) {
            trt_name <- d$treatments[i]
            cat(sprintf("  %-12s: %.1f%% probability of being best (lowest toxicity)\n",
                        trt_name, p_best[i] * 100))
          }
        } else {
          cat("  Rankings contain non-numeric values\n")
        }
      } else {
        cat("  Rankings format not recognized\n")
      }
    } else {
      cat("  Rankings not available\n")
    }
  }, error = function(e) {
    cat("  Error extracting rankings:", e$message, "\n")
  })
  cat("\n")
}

# -----------------------------------------------------------------------------
# Save Comparison Results
# -----------------------------------------------------------------------------

cat("=====================================================\n")
cat("=== SAVING COMPARISON RESULTS ===\n")
cat("=====================================================\n\n")

output_dir <- "output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# DIC comparison
write.csv(dic_comparison, file.path(output_dir, "model_comparison_dic.csv"), row.names = FALSE)
cat("Saved:", file.path(output_dir, "model_comparison_dic.csv"), "\n")

# OR comparison for key treatment
write.csv(or_comparison, file.path(output_dir, "model_comparison_13hetho.csv"), row.names = FALSE)
cat("Saved:", file.path(output_dir, "model_comparison_13hetho.csv"), "\n")

cat("\n=== Model Comparison Complete ===\n")
