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
                          seed = 12345, save_samples = TRUE) {

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
      mcmc.samples = save_samples,  # Return raw MCMC samples for pairwise comparisons

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
        mcmc.samples = save_samples,  # Return raw MCMC samples

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
  n_arms <- nrow(nma_data)
  convergence_status <- check_convergence(result, output_dir, model_name, n_arms = n_arms)

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
# Check convergence via PSRF and extract model fit statistics
# -----------------------------------------------------------------------------

check_convergence <- function(result, output_dir, model_name, n_arms = NA) {
  convergence_status <- list(
    converged = NA,
    max_psrf = NA,
    psrf_threshold = 1.05,
    psrf_file = NA,
    dic_value = NA,
    pd_value = NA,
    dbar_value = NA,
    dres_value = NA,      # Posterior mean total residual deviance
    n_arms = n_arms,      # Number of data points (study arms)
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

  # Extract DIC - pcnetmeta returns matrix with rownames: D.bar, pD, DIC
  # D.bar = deviance at posterior mean D(θ̄)
  # pD = effective number of parameters
  # DIC = D.bar + pD
  if (!is.null(result$DIC)) {
    tryCatch({
      dic_vec <- result$DIC
      if (is.matrix(dic_vec) || is.array(dic_vec)) {
        # pcnetmeta returns a matrix with rownames
        convergence_status$dbar_value <- dic_vec["D.bar", 1]  # D(θ̄)
        convergence_status$pd_value <- dic_vec["pD", 1]
        convergence_status$dic_value <- dic_vec["DIC", 1]
      } else if (is.numeric(dic_vec) && length(dic_vec) >= 3) {
        # Fallback: assume order is D.bar, pD, DIC
        convergence_status$dbar_value <- dic_vec[1]
        convergence_status$pd_value <- dic_vec[2]
        convergence_status$dic_value <- dic_vec[3]
      }
      cat("\nD.bar (deviance at posterior mean):", round(convergence_status$dbar_value, 2), "\n")
      cat("pD (effective parameters):", round(convergence_status$pd_value, 2), "\n")
      cat("DIC:", round(convergence_status$dic_value, 2), "\n")
    }, error = function(e) {
      cat("\nCould not extract DIC:", e$message, "\n")
    })
  }

  # Extract posterior mean of total residual deviance from MCMC samples
  # totresdev = D̄res = posterior mean of residual deviance (NOT the same as D.bar)
  # D.bar = D(θ̄) = deviance evaluated at posterior mean of parameters
  # For absolute fit: totresdev should ≈ number of unconstrained data points (arms)
  if (!is.null(result$mcmc.samples)) {
    tryCatch({
      samples_matrix <- do.call(rbind, result$mcmc.samples)
      if ("totresdev" %in% colnames(samples_matrix)) {
        convergence_status$dres_value <- mean(samples_matrix[, "totresdev"])
        cat("totresdev (posterior mean):", round(convergence_status$dres_value, 2), "\n")
        if (!is.na(n_arms)) {
          ratio <- convergence_status$dres_value / n_arms
          cat(sprintf("Dres/n_arms ratio: %.2f (should be ~1 for good fit)\n", ratio))
        }
      }
    }, error = function(e) {
      cat("Could not extract totresdev from MCMC samples:", e$message, "\n")
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

  # Save MCMC samples separately if they exist
  if (!is.null(result$mcmc.samples)) {
    saveRDS(result$mcmc.samples, file.path(output_dir, paste0(model_name, "_mcmc_samples.rds")))
    cat("Saved:", file.path(output_dir, paste0(model_name, "_mcmc_samples.rds")), "\n")
  }

  cat("\nSaved:", file.path(output_dir, paste0(model_name, "_result.rds")), "\n")
  cat("Saved:", file.path(output_dir, paste0(model_name, "_convergence.rds")), "\n")
}

# -----------------------------------------------------------------------------
# Pairwise probability comparison
# -----------------------------------------------------------------------------

#' Calculate probability that one treatment's OR exceeds another's
#'
#' @param mcmc_samples mcmc.list object from nma.ab.bin with mcmc.samples=TRUE
#' @param trt1 Name of first treatment (e.g., "2846hetho")
#' @param trt2 Name of second treatment (e.g., "13hetho")
#' @param reference Reference treatment for OR calculation (e.g., "WT_Clean")
#' @return List with probability and sample statistics
calculate_pairwise_probability <- function(mcmc_samples, trt1, trt2, reference = NULL) {

  if (is.null(mcmc_samples)) {
    stop("No MCMC samples available. Re-run model with mcmc.samples=TRUE")
  }

  # Combine chains into single matrix
  samples_matrix <- do.call(rbind, mcmc_samples)
  param_names <- colnames(samples_matrix)

  cat("\n=== Pairwise Probability Calculation ===\n")
  cat(sprintf("Comparing: %s vs %s\n", trt1, trt2))

  # Find the LOR columns for each treatment vs reference
  # pcnetmeta names them like "LOR[trt1,ref]" or "LOR.trt1.ref"
  find_lor_column <- function(trt, ref, params) {
    # Try different naming patterns
    patterns <- c(
      sprintf("LOR\\[%s,%s\\]", trt, ref),
      sprintf("LOR\\[%s, %s\\]", trt, ref),
      sprintf("LOR\\.%s\\.%s", trt, ref),
      sprintf("LOR_%s_%s", trt, ref)
    )

    for (pat in patterns) {
      matches <- grep(pat, params, value = TRUE)
      if (length(matches) > 0) return(matches[1])
    }

    # Try reversed order
    patterns_rev <- c(
      sprintf("LOR\\[%s,%s\\]", ref, trt),
      sprintf("LOR\\[%s, %s\\]", ref, trt),
      sprintf("LOR\\.%s\\.%s", ref, trt),
      sprintf("LOR_%s_%s", ref, trt)
    )

    for (pat in patterns_rev) {
      matches <- grep(pat, params, value = TRUE)
      if (length(matches) > 0) return(paste0("-", matches[1]))  # Negate
    }

    return(NULL)
  }

  # If no reference specified, find direct comparison between trt1 and trt2
  if (is.null(reference)) {
    # Try to find direct LOR comparison
    col1 <- find_lor_column(trt1, trt2, param_names)
    if (!is.null(col1)) {
      if (startsWith(col1, "-")) {
        samples_diff <- -samples_matrix[, substring(col1, 2)]
      } else {
        samples_diff <- samples_matrix[, col1]
      }
    } else {
      stop(sprintf("Cannot find LOR comparison between %s and %s", trt1, trt2))
    }
  } else {
    # Find LOR vs reference for each treatment
    col1 <- find_lor_column(trt1, reference, param_names)
    col2 <- find_lor_column(trt2, reference, param_names)

    if (is.null(col1) || is.null(col2)) {
      cat("Available parameters:\n")
      cat(paste(head(param_names, 30), collapse = "\n"), "\n")
      stop(sprintf("Cannot find LOR columns for %s or %s vs %s", trt1, trt2, reference))
    }

    # Extract samples
    if (startsWith(col1, "-")) {
      samples1 <- -samples_matrix[, substring(col1, 2)]
    } else {
      samples1 <- samples_matrix[, col1]
    }

    if (startsWith(col2, "-")) {
      samples2 <- -samples_matrix[, substring(col2, 2)]
    } else {
      samples2 <- samples_matrix[, col2]
    }

    samples_diff <- samples1 - samples2
  }

  # Calculate probability
  p_trt1_greater <- mean(samples_diff > 0)

  # Summary statistics
  or_diff_median <- median(exp(samples_diff))
  or_diff_ci <- quantile(exp(samples_diff), c(0.025, 0.975))

  cat(sprintf("\nResults:\n"))
  cat(sprintf("  P(%s OR > %s OR) = %.3f\n", trt1, trt2, p_trt1_greater))
  cat(sprintf("  P(%s OR < %s OR) = %.3f\n", trt1, trt2, 1 - p_trt1_greater))
  cat(sprintf("\n  Median OR ratio (%s/%s): %.2f\n", trt1, trt2, or_diff_median))
  cat(sprintf("  95%% CrI: [%.2f, %.2f]\n", or_diff_ci[1], or_diff_ci[2]))

  # Interpretation
  cat("\nInterpretation:\n")
  if (p_trt1_greater > 0.6) {
    cat(sprintf("  %s likely has higher OR (more toxic) than %s\n", trt1, trt2))
  } else if (p_trt1_greater < 0.4) {
    cat(sprintf("  %s likely has lower OR (less toxic) than %s\n", trt1, trt2))
  } else {
    cat("  Ranking is essentially a coin flip - no meaningful difference\n")
  }

  return(list(
    p_trt1_greater = p_trt1_greater,
    p_trt2_greater = 1 - p_trt1_greater,
    or_ratio_median = or_diff_median,
    or_ratio_ci = or_diff_ci,
    n_samples = length(samples_diff)
  ))
}
