# =============================================================================
# R/common/model_runner.R
# Shared NMA model execution utilities
# =============================================================================

library(pcnetmeta)
library(rjags)

# -----------------------------------------------------------------------------
# Default MCMC configuration (improved for convergence)
# -----------------------------------------------------------------------------

get_default_mcmc_config <- function() {
  list(
    n.adapt   = 10000,   # Adaptation iterations
    n.iter    = 150000,  # Total iterations (increased)
    n.burnin  = 75000,   # Burn-in (increased for better convergence)
    n.chains  = 3,       # Number of chains
    n.thin    = 10       # Thinning interval (reduces autocorrelation)
  )
}

# -----------------------------------------------------------------------------
# Run arm-based NMA model
# -----------------------------------------------------------------------------

run_nma_model <- function(nma_data, treatment_names, mcmc_config = NULL,
                          model_name = "nma", output_dir = "output",
                          seed = 12345) {

  if (is.null(mcmc_config)) {
    mcmc_config <- get_default_mcmc_config()
  }

  cat("=== Running Arm-Based NMA ===\n")
  cat("Model:", model_name, "\n")
  cat("Data:", nrow(nma_data), "arms\n")
  cat("Treatments:", paste(treatment_names, collapse = ", "), "\n\n")

  cat("MCMC Configuration:\n")
  cat("  Adaptation:", mcmc_config$n.adapt, "\n")
  cat("  Iterations:", mcmc_config$n.iter, "\n")
  cat("  Burn-in:", mcmc_config$n.burnin, "\n")
  cat("  Thinning:", mcmc_config$n.thin, "\n")
  cat("  Chains:", mcmc_config$n.chains, "\n")
  effective_samples <- (mcmc_config$n.iter - mcmc_config$n.burnin) / mcmc_config$n.thin
  cat("  Effective samples per chain:", effective_samples, "\n\n")

  set.seed(seed)

  # Try het_cor model first (most flexible)
  cat("Running het_cor model (heterogeneous variance with correlation)...\n")
  cat("This may take several minutes.\n\n")

  result_het_cor <- tryCatch({
    nma.ab.bin(
      s.id = s.id,
      t.id = t.id,
      event.n = r,
      total.n = n,
      data = nma_data,
      trtname = treatment_names,

      model = "het_cor",
      link = "probit",
      param = c("AR", "OR", "LOR", "RD", "rank.prob"),
      prior.type = "invwishart",

      n.adapt = mcmc_config$n.adapt,
      n.iter = mcmc_config$n.iter,
      n.burnin = mcmc_config$n.burnin,
      n.chains = mcmc_config$n.chains,
      n.thin = mcmc_config$n.thin,

      conv.diag = TRUE,
      trace = c("LOR"),
      dic = TRUE,
      postdens = FALSE,

      higher.better = FALSE,
      digits = 4
    )
  }, error = function(e) {
    cat("ERROR in het_cor model:", e$message, "\n")
    return(NULL)
  })

  # Fallback to hom_eqcor if het_cor fails
  if (is.null(result_het_cor)) {
    cat("\nhet_cor failed. Trying hom_eqcor (simpler model)...\n\n")

    result <- tryCatch({
      nma.ab.bin(
        s.id = s.id,
        t.id = t.id,
        event.n = r,
        total.n = n,
        data = nma_data,
        trtname = treatment_names,

        model = "hom_eqcor",
        link = "probit",
        param = c("AR", "OR", "LOR", "RD", "rank.prob"),
        prior.type = "invgamma",
        a = 0.001,
        b = 0.001,

        n.adapt = mcmc_config$n.adapt,
        n.iter = mcmc_config$n.iter,
        n.burnin = mcmc_config$n.burnin,
        n.chains = mcmc_config$n.chains,
        n.thin = mcmc_config$n.thin,

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

    model_type <- "hom_eqcor"
  } else {
    result <- result_het_cor
    model_type <- "het_cor"
  }

  if (is.null(result)) {
    stop("Both models failed. Check data and JAGS installation.")
  }

  cat("\n=== Model Completed ===\n")
  cat("Model type:", model_type, "\n")

  # Check convergence and extract DIC
  convergence_status <- check_convergence(result, output_dir, model_name)

  # Save results
  save_model_results(result, model_type, convergence_status,
                     treatment_names, model_name, output_dir)

  return(list(
    result = result,
    model_type = model_type,
    convergence_status = convergence_status,
    treatment_names = treatment_names
  ))
}

# -----------------------------------------------------------------------------
# Check convergence via PSRF
# -----------------------------------------------------------------------------

check_convergence <- function(result, output_dir, model_name) {
  convergence_status <- list(
    converged = NA,
    max_psrf = NA,
    psrf_threshold = 1.05,
    psrf_file = NA,
    dic_value = NA,
    pd_value = NA,
    dbar_value = NA,
    warning_message = NULL
  )

  # Find PSRF file
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

    tryCatch({
      psrf_lines <- readLines(psrf_file)
      psrf_values <- c()
      psrf_params <- c()

      for (line in psrf_lines) {
        if (nchar(trimws(line)) == 0) next
        if (grepl("^#|^Parameter|^---", line)) next

        parts <- strsplit(trimws(line), "\\s+")[[1]]
        if (length(parts) >= 2) {
          numeric_vals <- suppressWarnings(as.numeric(parts))
          valid_nums <- numeric_vals[!is.na(numeric_vals)]
          if (length(valid_nums) > 0) {
            psrf_val <- valid_nums[length(valid_nums)]
            if (psrf_val > 0 && psrf_val < 100) {
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

        if (max_psrf <= convergence_status$psrf_threshold) {
          convergence_status$converged <- TRUE
          cat("  Status: CONVERGED (all PSRF <= 1.05)\n")
        } else {
          convergence_status$converged <- FALSE
          failed_idx <- which(psrf_values > convergence_status$psrf_threshold)
          cat(sprintf("  Status: NOT CONVERGED (%d parameters exceed threshold)\n",
                      length(failed_idx)))
          cat("  Recommendation: Increase n.iter or n.burnin and re-run.\n")
          convergence_status$warning_message <- sprintf(
            "CONVERGENCE WARNING: %d parameter(s) have PSRF > 1.05", length(failed_idx))
        }
      }
    }, error = function(e) {
      cat("  Error reading PSRF file:", e$message, "\n")
      convergence_status$warning_message <- paste("Error reading PSRF:", e$message)
    })
  } else {
    cat("\nNo PSRF file found. Convergence diagnostics unavailable.\n")
    convergence_status$warning_message <- "No PSRF file generated"
  }

  # Extract DIC - pcnetmeta returns numeric vector: c(DIC, pD, Dbar)
  if (!is.null(result$DIC)) {
    tryCatch({
      dic_vec <- result$DIC
      if (is.numeric(dic_vec) && length(dic_vec) >= 2) {
        convergence_status$dic_value <- dic_vec[1]
        convergence_status$pd_value <- dic_vec[2]
        if (length(dic_vec) >= 3) {
          convergence_status$dbar_value <- dic_vec[3]
        }
        cat("\nDIC:", round(convergence_status$dic_value, 2), "\n")
        cat("pD (effective parameters):", round(convergence_status$pd_value, 2), "\n")
        if (!is.na(convergence_status$dbar_value)) {
          cat("Dbar (mean deviance):", round(convergence_status$dbar_value, 2), "\n")
        }
      } else if (is.list(dic_vec)) {
        convergence_status$dic_value <- dic_vec$DIC
        convergence_status$pd_value <- dic_vec$pD
        cat("\nDIC:", round(convergence_status$dic_value, 2), "\n")
        cat("pD (effective parameters):", round(convergence_status$pd_value, 2), "\n")
      } else {
        cat("\nDIC (raw):", dic_vec, "\n")
      }
    }, error = function(e) {
      cat("\nCould not extract DIC:", e$message, "\n")
    })
  }

  return(convergence_status)
}

# -----------------------------------------------------------------------------
# Save model results
# -----------------------------------------------------------------------------

save_model_results <- function(result, model_type, convergence_status,
                               treatment_names, model_name, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  saveRDS(result, file.path(output_dir, paste0(model_name, "_result.rds")))
  saveRDS(model_type, file.path(output_dir, paste0(model_name, "_model_type.rds")))
  saveRDS(convergence_status, file.path(output_dir, paste0(model_name, "_convergence.rds")))
  saveRDS(treatment_names, file.path(output_dir, paste0(model_name, "_treatments.rds")))

  cat("\nSaved:", file.path(output_dir, paste0(model_name, "_result.rds")), "\n")
  cat("Saved:", file.path(output_dir, paste0(model_name, "_convergence.rds")), "\n")
}
