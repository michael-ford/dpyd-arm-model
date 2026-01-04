# =============================================================================
# 02_run_nma.R
# Run arm-based network meta-analysis using pcnetmeta
# =============================================================================

library(pcnetmeta)
library(rjags)

# -----------------------------------------------------------------------------
# Load prepared data
# -----------------------------------------------------------------------------

nma_data <- readRDS("output/nma_data.rds")
treatment_names <- readRDS("output/treatment_names.rds")

cat("=== Running Arm-Based NMA ===\n")
cat("Data: ", nrow(nma_data), " arms\n")
cat("Treatments: ", paste(treatment_names, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# Model Configuration
# Based on project-background.md recommendations
# -----------------------------------------------------------------------------

# MCMC settings (conservative for sparse data)
mcmc_config <- list(
  n.adapt   = 10000,   # Adaptation iterations
  n.iter    = 100000,  # Total iterations
  n.burnin  = 50000,   # Burn-in (discarded)
  n.chains  = 3        # Number of chains for convergence diagnostics
)

cat("MCMC Configuration:\n")
cat("  Adaptation:", mcmc_config$n.adapt, "\n")
cat("  Iterations:", mcmc_config$n.iter, "\n")
cat("  Burn-in:", mcmc_config$n.burnin, "\n")
cat("  Chains:", mcmc_config$n.chains, "\n\n")

# -----------------------------------------------------------------------------
# Run Arm-Based Model (het_cor)
# Model 4 from Hong et al. - Random Effects with Correlation
# -----------------------------------------------------------------------------

cat("Running het_cor model (heterogeneous variance with correlation)...\n")
cat("This may take several minutes.\n\n")

set.seed(12345)  # Reproducibility

result_het_cor <- tryCatch({
  nma.ab.bin(
    s.id = s.id,
    t.id = t.id,
    event.n = r,
    total.n = n,
    data = nma_data,
    trtname = treatment_names,

    # Model specification
    model = "het_cor",           # Heterogeneous correlation (most flexible)
    link = "probit",             # Probit link (recommended for sparse data)
    param = c("AR", "OR", "LOR", "RD", "rank.prob"),

    # Priors (inverse-Wishart for het_cor, weakly informative)
    prior.type = "invwishart",

    # MCMC settings
    n.adapt = mcmc_config$n.adapt,
    n.iter = mcmc_config$n.iter,
    n.burnin = mcmc_config$n.burnin,
    n.chains = mcmc_config$n.chains,

    # Diagnostics
    conv.diag = TRUE,            # Generate PSRF diagnostics
    trace = c("LOR"),            # Trace plots for log odds ratios
    dic = TRUE,                  # Calculate DIC for model comparison
    postdens = FALSE,            # Skip posterior density plots

    # Output
    higher.better = FALSE,       # Lower toxicity is better
    digits = 4
  )
}, error = function(e) {
  cat("ERROR in het_cor model:", e$message, "\n")
  return(NULL)
})

# -----------------------------------------------------------------------------
# Fallback: hom_eqcor if het_cor fails
# -----------------------------------------------------------------------------

if (is.null(result_het_cor)) {
  cat("\nhet_cor failed. Trying hom_eqcor (simpler model)...\n\n")

  result_hom_eqcor <- tryCatch({
    nma.ab.bin(
      s.id = s.id,
      t.id = t.id,
      event.n = r,
      total.n = n,
      data = nma_data,
      trtname = treatment_names,

      model = "hom_eqcor",         # Homogeneous equal correlation
      link = "probit",
      param = c("AR", "OR", "LOR", "RD", "rank.prob"),

      prior.type = "invgamma",     # For hom_eqcor
      a = 0.001,                   # Weakly informative
      b = 0.001,

      n.adapt = mcmc_config$n.adapt,
      n.iter = mcmc_config$n.iter,
      n.burnin = mcmc_config$n.burnin,
      n.chains = mcmc_config$n.chains,

      conv.diag = TRUE,
      trace = c("LOR"),
      dic = TRUE,
      postdens = FALSE,

      higher.better = FALSE,
      digits = 4
    )
  }, error = function(e) {
    cat("ERROR in hom_eqcor model:", e$message, "\n")
    return(NULL)
  })

  result <- result_hom_eqcor
  model_used <- "hom_eqcor"
} else {
  result <- result_het_cor
  model_used <- "het_cor"
}

# -----------------------------------------------------------------------------
# Check convergence with automated PSRF validation
# -----------------------------------------------------------------------------

if (!is.null(result)) {
  cat("\n=== Model Completed ===\n")
  cat("Model used:", model_used, "\n")

  # Initialize convergence status
  convergence_status <- list(
    converged = NA,
    max_psrf = NA,
    psrf_threshold = 1.05,
    psrf_file = NA,
    warning_message = NULL
  )

  # PSRF diagnostics are saved to file by pcnetmeta
  # Find, read, and validate PSRF values
  # Check both current directory and output directory
  psrf_files <- list.files(pattern = "PSRF.*\\.txt$")
  if (length(psrf_files) == 0) {
    psrf_files <- list.files(pattern = ".*[Cc]onvergence.*\\.txt$")
  }
  if (length(psrf_files) == 0) {
    psrf_files <- list.files(".", pattern = ".*\\.txt$", full.names = FALSE)
    psrf_files <- psrf_files[grepl("PSRF|[Cc]onvergence|[Dd]iagnostic", psrf_files)]
  }

  if (length(psrf_files) > 0) {
    psrf_file <- psrf_files[1]
    convergence_status$psrf_file <- psrf_file
    cat("\nPSRF diagnostics file:", psrf_file, "\n")

    # Read and parse PSRF file
    tryCatch({
      psrf_lines <- readLines(psrf_file)

      # Extract numeric PSRF values from the file
      # PSRF files typically have format: "parameter_name    value"
      # Look for lines with numeric values
      psrf_values <- c()
      psrf_params <- c()

      for (line in psrf_lines) {
        # Skip empty lines and headers
        if (nchar(trimws(line)) == 0) next
        if (grepl("^#|^Parameter|^---", line)) next

        # Try to extract parameter name and PSRF value
        parts <- strsplit(trimws(line), "\\s+")[[1]]
        if (length(parts) >= 2) {
          # Last numeric value is typically the PSRF
          numeric_vals <- suppressWarnings(as.numeric(parts))
          valid_nums <- numeric_vals[!is.na(numeric_vals)]
          if (length(valid_nums) > 0) {
            psrf_val <- valid_nums[length(valid_nums)]
            if (psrf_val > 0 && psrf_val < 100) {  # Sanity check
              psrf_values <- c(psrf_values, psrf_val)
              psrf_params <- c(psrf_params, parts[1])
            }
          }
        }
      }

      if (length(psrf_values) > 0) {
        max_psrf <- max(psrf_values)
        convergence_status$max_psrf <- max_psrf

        cat(sprintf("\nPSRF Validation Results:\n"))
        cat(sprintf("  Total parameters checked: %d\n", length(psrf_values)))
        cat(sprintf("  Maximum PSRF: %.4f\n", max_psrf))
        cat(sprintf("  Threshold: %.2f\n", convergence_status$psrf_threshold))

        # Check convergence
        if (max_psrf <= convergence_status$psrf_threshold) {
          convergence_status$converged <- TRUE
          cat("  Status: CONVERGED (all PSRF <= 1.05)\n")
        } else {
          convergence_status$converged <- FALSE
          # Find which parameters failed
          failed_idx <- which(psrf_values > convergence_status$psrf_threshold)
          failed_params <- psrf_params[failed_idx]
          failed_vals <- psrf_values[failed_idx]

          warning_msg <- sprintf(
            "CONVERGENCE WARNING: %d parameter(s) have PSRF > 1.05",
            length(failed_idx)
          )
          convergence_status$warning_message <- warning_msg

          cat(sprintf("  Status: NOT CONVERGED\n"))
          cat(sprintf("  WARNING: %d parameter(s) exceed threshold:\n", length(failed_idx)))

          # Show up to 5 worst offenders
          order_idx <- order(failed_vals, decreasing = TRUE)
          show_n <- min(5, length(order_idx))
          for (i in 1:show_n) {
            idx <- order_idx[i]
            cat(sprintf("    - %s: %.4f\n", failed_params[idx], failed_vals[idx]))
          }
          if (length(failed_idx) > 5) {
            cat(sprintf("    ... and %d more\n", length(failed_idx) - 5))
          }

          cat("\n  Recommendation: Increase n.iter or n.burnin and re-run.\n")
          warning(warning_msg)
        }
      } else {
        cat("  Could not parse PSRF values from file.\n")
        convergence_status$warning_message <- "Could not parse PSRF file"
      }
    }, error = function(e) {
      cat("  Error reading PSRF file:", e$message, "\n")
      convergence_status$warning_message <- paste("Error reading PSRF file:", e$message)
    })
  } else {
    cat("\nNo PSRF file found. Convergence diagnostics unavailable.\n")
    convergence_status$warning_message <- "No PSRF file generated"
  }

  # DIC - handle different result structures
  if (!is.null(result$DIC)) {
    tryCatch({
      # pcnetmeta stores DIC as a named numeric vector or list
      if (is.list(result$DIC)) {
        cat("\nDIC:", result$DIC$DIC, "\n")
        cat("pD (effective parameters):", result$DIC$pD, "\n")
      } else if (is.numeric(result$DIC)) {
        # DIC might be stored as named vector: c(DIC=..., pD=...)
        if ("DIC" %in% names(result$DIC)) {
          cat("\nDIC:", result$DIC["DIC"], "\n")
          cat("pD (effective parameters):", result$DIC["pD"], "\n")
        } else {
          cat("\nDIC (raw):", result$DIC, "\n")
        }
      } else {
        cat("\nDIC structure:", class(result$DIC), "\n")
        print(result$DIC)
      }
    }, error = function(e) {
      cat("\nCould not extract DIC:", e$message, "\n")
      cat("DIC object class:", class(result$DIC), "\n")
    })
  }

  # Save result with convergence status
  saveRDS(result, "output/nma_result.rds")
  saveRDS(model_used, "output/model_used.rds")
  saveRDS(convergence_status, "output/convergence_status.rds")

  cat("\nSaved: output/nma_result.rds\n")
  cat("Saved: output/convergence_status.rds\n")
} else {
  stop("Both models failed. Check data and JAGS installation.")
}
