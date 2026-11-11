# =============================================================================
# Fixed design-by-treatment inconsistency diagnostic
# =============================================================================

.dbt_require_functions <- function(functions) {
  missing <- functions[!vapply(functions, exists, logical(1L), mode = "function")]
  if (length(missing)) {
    stop(
      "Source R/network_geometry.R and R/fixed_study_model.R first; missing: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
}

.dbt_rank <- function(x, tolerance = 1e-8) {
  if (!length(x) || !nrow(x) || !ncol(x)) return(0L)
  as.integer(qr(x, tol = tolerance)$rank)
}

.dbt_validate_prior_sd <- function(gamma_prior_sd) {
  if (length(gamma_prior_sd) != 1L || is.na(gamma_prior_sd) ||
      !is.numeric(gamma_prior_sd) || !gamma_prior_sd %in% c(1, 2)) {
    stop("gamma_prior_sd must be 2 for the main diagnostic or 1 for the bounded check.", call. = FALSE)
  }
  as.numeric(gamma_prior_sd)
}

# Conservative reporting policy fixed before seeing the inconsistency results.
# An 80% posterior/prior SD ratio is 0.8^2 = 0.64 on the variance scale. The
# prior comparison additionally limits median movement to half of the pooled
# posterior SD so conclusions cannot depend materially on the chosen bound.
design_by_treatment_gate_policy <- function() {
  list(
    posterior_to_prior_sd_ratio = 0.80,
    posterior_to_prior_variance_ratio = 0.64,
    maximum_standardized_median_shift = 0.50,
    information_rationale = paste(
      "A reported inconsistency direction must have posterior SD no more than",
      "80% of its prior SD; squaring gives the maximum variance ratio 0.64."
    ),
    stability_rationale = paste(
      "The SD=2 and SD=1 fits must agree in global and interaction classifications,",
      "and each interaction median may move by no more than 0.5 pooled posterior SD."
    )
  )
}

.dbt_sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) stop("Model file is unavailable for fingerprinting: ", path, call. = FALSE)
  if (exists("sha256_file", mode = "function")) return(sha256_file(path))
  executable <- Sys.which("sha256sum")
  if (nzchar(executable)) {
    output <- system2(executable, shQuote(path), stdout = TRUE, stderr = TRUE)
  } else {
    executable <- Sys.which("shasum")
    if (!nzchar(executable)) stop("No SHA-256 utility is available for model fingerprinting.", call. = FALSE)
    output <- system2(executable, c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  }
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 model fingerprinting failed.", call. = FALSE)
  fingerprint <- sub(" .*", "", output[1L])
  if (!grepl("^[0-9a-f]{64}$", fingerprint)) stop("Invalid SHA-256 model fingerprint.", call. = FALSE)
  fingerprint
}

# Construct a deterministic orthonormal basis for the complement of col(A).
# Projected coordinate vectors are considered in row order. This removes the
# arbitrary sign and rotation choices returned by a generic null-space SVD.
design_by_treatment_basis <- function(geometry, tolerance = 1e-8) {
  if (!is.list(geometry) || !all(c("A", "m", "q", "rank_A", "treatments") %in% names(geometry))) {
    stop("geometry must be an analyze_network_geometry() result.", call. = FALSE)
  }
  A <- as.matrix(geometry$A)
  m <- as.integer(geometry$m)
  q <- as.integer(geometry$q)
  expected_rank <- length(geometry$treatments) - 1L
  if (!identical(dim(A), c(m, expected_rank))) {
    stop("Geometry A dimensions do not agree with m and the treatment count.", call. = FALSE)
  }
  if (geometry$rank_A != expected_rank || .dbt_rank(A, tolerance) != expected_rank) {
    stop("Geometry A must have rank T - 1.", call. = FALSE)
  }
  if (q <= 0L || q != m - expected_rank) {
    stop("The network has no estimable design-by-treatment inconsistency degrees of freedom.", call. = FALSE)
  }

  gram <- crossprod(A)
  if (any(!is.finite(gram)) || .dbt_rank(gram, tolerance) != ncol(A)) {
    stop("Geometry A has a singular or non-finite Gram matrix.", call. = FALSE)
  }
  projector <- diag(m) - A %*% solve(gram, t(A))
  H <- matrix(0, nrow = m, ncol = 0L)
  for (coordinate in seq_len(m)) {
    candidate <- projector[, coordinate]
    if (ncol(H)) candidate <- candidate - H %*% as.vector(crossprod(H, candidate))
    norm <- sqrt(sum(candidate^2))
    if (is.finite(norm) && norm > tolerance) {
      candidate <- candidate / norm
      first <- which(abs(candidate) > tolerance)[1L]
      if (candidate[first] < 0) candidate <- -candidate
      H <- cbind(H, candidate)
      if (ncol(H) == q) break
    }
  }
  if (ncol(H) != q) stop("Could not construct the full inconsistency basis H.", call. = FALSE)
  colnames(H) <- sprintf("gamma[%d]", seq_len(q))
  rownames(H) <- rownames(A)

  orthonormal_error <- max(abs(crossprod(H) - diag(q)))
  complement_error <- max(abs(crossprod(H, A)))
  if (!is.finite(orthonormal_error) || orthonormal_error > tolerance * 10 ||
      !is.finite(complement_error) || complement_error > tolerance * 10) {
    stop("Constructed H failed the orthonormal-complement check.", call. = FALSE)
  }
  H
}

.dbt_study_design_map <- function(arms, geometry) {
  studies <- sort(unique(as.character(arms$Study)), method = "radix")
  keys <- vapply(studies, function(study) {
    paste(sort(unique(arms$treatment[arms$Study == study])), collapse = "-")
  }, character(1L))
  design_id <- geometry$designs$design_id[match(keys, geometry$designs$design)]
  if (anyNA(design_id)) stop("At least one study has no unique geometry design.", call. = FALSE)
  data.frame(Study = studies, design = keys, design_id = design_id, stringsAsFactors = FALSE)
}

# Build the JAGS design rows from the final, prepared arm input. Baseline arms
# use row 1 of H_jags (all zero); every other arm maps to one A/H contrast row.
prepare_design_by_treatment_data <- function(prepared, gamma_prior_sd = 2, tolerance = 1e-8) {
  .dbt_require_functions(c("analyze_network_geometry"))
  gamma_prior_sd <- .dbt_validate_prior_sd(gamma_prior_sd)
  if (!is.list(prepared) || !all(c("arms", "jags_data", "treatment_names", "network_connected") %in% names(prepared))) {
    stop("prepared must be a prepare_fixed_study_data() result.", call. = FALSE)
  }
  arms <- prepared$arms
  required <- c("Study", "treatment", "R", "N")
  if (!is.data.frame(arms) || !all(required %in% names(arms)) || !nrow(arms)) {
    stop("Prepared data are missing required arm fields.", call. = FALSE)
  }
  if (length(prepared$dropped_single_arm_studies)) {
    stop("The inconsistency diagnostic does not accept dropped or retained single-arm studies.", call. = FALSE)
  }
  if (any(table(arms$Study) < 2L)) {
    stop("The inconsistency diagnostic does not accept single-arm studies.", call. = FALSE)
  }
  if (!isTRUE(prepared$network_connected)) stop("The treatment network must be connected.", call. = FALSE)
  if (anyDuplicated(arms[c("Study", "treatment")])) {
    stop("Prepared arms contain aliased duplicate study-treatment rows.", call. = FALSE)
  }
  wt_id <- match("WT", prepared$treatment_names)
  if (!is.na(wt_id) && wt_id != 1L) {
    stop("WT must be treatment ID 1 so the diagnostic matches the primary model reference.", call. = FALSE)
  }

  geometry_input <- data.frame(
    Study = as.character(arms$Study),
    t.id = as.integer(arms$treatment),
    stringsAsFactors = FALSE
  )
  geometry <- analyze_network_geometry(geometry_input)
  treatment_count <- length(geometry$treatments)
  expected_rank <- treatment_count - 1L
  if (!identical(sort(geometry$treatments), seq_len(prepared$jags_data$Ntrt))) {
    stop("Geometry treatment IDs must cover the primary model nodes without aliases or gaps.", call. = FALSE)
  }
  if (!geometry$connectedness$connected || geometry$rank_A != expected_rank ||
      geometry$q <= 0L || geometry$q != geometry$q_from_rank ||
      !isTRUE(geometry$arm_level_rank_check$passes) ||
      geometry$q != geometry$arm_level_rank_check$inconsistency_rank_increase) {
    stop("Network geometry failed connectivity, rank(A), q, or arm-rank preflight.", call. = FALSE)
  }
  expected_baselines <- vapply(geometry$designs$treatments, function(x) {
    if (1L %in% x) 1L else min(x)
  }, integer(1L))
  if (!identical(as.integer(geometry$designs$baseline), expected_baselines)) {
    stop("Design baselines must be WT when present and otherwise the lowest treatment ID.", call. = FALSE)
  }

  H <- design_by_treatment_basis(geometry, tolerance)
  study_design <- .dbt_study_design_map(arms, geometry)
  design_for_arm <- study_design$design_id[match(arms$Study, study_design$Study)]
  baseline_for_arm <- geometry$designs$baseline[match(design_for_arm, geometry$designs$design_id)]
  contrast_row <- integer(nrow(arms))
  for (arm in seq_len(nrow(arms))) {
    if (arms$treatment[arm] == baseline_for_arm[arm]) next
    candidates <- which(
      geometry$A_mapping$design_id == design_for_arm[arm] &
        geometry$A_mapping$treatment == arms$treatment[arm]
    )
    if (length(candidates) != 1L) {
      stop("A non-baseline arm does not map to exactly one design contrast; aliases are not allowed.", call. = FALSE)
    }
    contrast_row[arm] <- candidates
  }
  if (any(contrast_row == 0L & arms$treatment != baseline_for_arm) ||
      any(contrast_row > 0L & arms$treatment == baseline_for_arm)) {
    stop("Baseline and non-baseline arm mappings are inconsistent.", call. = FALSE)
  }

  H_jags <- rbind(rep(0, geometry$q), H)
  interaction_row <- contrast_row + 1L
  checks <- list(
    connected = geometry$connectedness$connected,
    no_single_arm_studies = all(table(arms$Study) >= 2L),
    rank_A_is_T_minus_1 = geometry$rank_A == expected_rank,
    q_positive = geometry$q > 0L,
    q_agrees_with_A_rank = geometry$q == geometry$q_from_rank,
    q_agrees_with_arm_rank = geometry$q == geometry$arm_level_rank_check$inconsistency_rank_increase,
    H_orthonormal = max(abs(crossprod(H) - diag(geometry$q))) <= tolerance * 10,
    H_complements_A = max(abs(crossprod(H, geometry$A))) <= tolerance * 10,
    no_aliases = !anyDuplicated(geometry$A_mapping[c("design_id", "treatment")]) &&
      identical(as.integer(geometry$A_mapping$row), seq_len(geometry$m)) &&
      .dbt_rank(cbind(geometry$A, H), tolerance) == geometry$m &&
      all(contrast_row >= 0L & contrast_row <= geometry$m),
    deterministic_design_baselines = identical(as.integer(geometry$designs$baseline), expected_baselines)
  )
  if (!all(unlist(checks, use.names = FALSE))) {
    stop("Design-by-treatment preflight failed.", call. = FALSE)
  }

  jags_data <- c(prepared$jags_data, list(
    Ninconsistency = as.integer(geometry$q),
    Ncontrast = as.integer(geometry$m),
    H = unname(H_jags),
    interaction_row = as.integer(interaction_row),
    gamma_prior_sd = gamma_prior_sd
  ))
  arm_mapping <- data.frame(
    arm = arms$arm,
    Study = arms$Study,
    treatment = arms$treatment,
    design_id = design_for_arm,
    design_baseline = baseline_for_arm,
    contrast_row = contrast_row,
    jags_interaction_row = interaction_row,
    stringsAsFactors = FALSE
  )
  list(
    jags_data = jags_data,
    geometry = geometry,
    H = H,
    H_jags = H_jags,
    arm_mapping = arm_mapping,
    interaction_mapping = geometry$A_mapping,
    gamma_prior_sd = gamma_prior_sd,
    checks = checks,
    passed = TRUE
  )
}

design_by_treatment_preflight <- prepare_design_by_treatment_data

.dbt_parameter_diagnostics <- function(samples, parameters, psrf_threshold, ess_threshold) {
  if (!inherits(samples, "mcmc.list")) stop("samples must be a coda mcmc.list.", call. = FALSE)
  if (!requireNamespace("coda", quietly = TRUE)) stop("Package 'coda' is required.", call. = FALSE)
  missing <- setdiff(parameters, colnames(samples[[1L]]))
  if (length(missing)) stop("Missing inconsistency diagnostic parameter(s): ", paste(missing, collapse = ", "), call. = FALSE)
  rows <- lapply(parameters, function(parameter) {
    chains <- lapply(samples, function(chain) {
      coda::mcmc(as.matrix(chain)[, parameter, drop = FALSE],
                 start = stats::start(chain), thin = coda::thin(chain))
    })
    one <- do.call(coda::mcmc.list, chains)
    psrf <- tryCatch(
      coda::gelman.diag(one, confidence = 0.95, transform = FALSE,
                        autoburnin = FALSE, multivariate = FALSE)$psrf[1L, ],
      error = function(e) c("Point est." = NA_real_, "Upper C.I." = NA_real_)
    )
    chain_ess <- vapply(chains, function(chain) {
      tryCatch(as.numeric(coda::effectiveSize(chain)[1L]), error = function(e) NA_real_)
    }, numeric(1L))
    data.frame(
      parameter = parameter,
      parameter_type = if (startsWith(parameter, "gamma[")) "gamma" else "design_interaction",
      psrf_point = unname(psrf["Point est."]),
      psrf_upper_95 = unname(psrf["Upper C.I."]),
      min_chain_effective_size = if (all(is.na(chain_ess))) NA_real_ else min(chain_ess, na.rm = TRUE),
      psrf_pass = is.finite(psrf["Upper C.I."]) && psrf["Upper C.I."] <= psrf_threshold,
      ess_pass = all(is.finite(chain_ess)) && all(chain_ess >= ess_threshold),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

design_by_treatment_diagnostics <- function(samples, q, n_interactions,
                                            psrf_threshold = 1.05,
                                            ess_threshold = 1000) {
  parameters <- c(sprintf("gamma[%d]", seq_len(q)),
                  sprintf("design_interaction[%d]", seq_len(n_interactions)))
  table <- .dbt_parameter_diagnostics(samples, parameters, psrf_threshold, ess_threshold)
  table$passed <- table$psrf_pass & table$ess_pass
  list(
    table = table,
    converged = all(table$passed),
    max_psrf_upper_95 = if (all(is.na(table$psrf_upper_95))) NA_real_ else max(table$psrf_upper_95, na.rm = TRUE),
    min_chain_effective_size = if (all(is.na(table$min_chain_effective_size))) NA_real_ else min(table$min_chain_effective_size, na.rm = TRUE),
    psrf_threshold = psrf_threshold,
    ess_threshold = ess_threshold
  )
}

.dbt_sample_matrix <- function(x) {
  if (is.list(x) && !inherits(x, "mcmc.list") && !is.null(x$samples)) x <- x$samples
  as.matrix(x)
}

.dbt_normalize_scalar_gamma_name <- function(samples) {
  if (!inherits(samples, "mcmc.list") || !"gamma" %in% colnames(samples[[1L]])) return(samples)
  chains <- lapply(samples, function(chain) {
    values <- as.matrix(chain)
    colnames(values)[colnames(values) == "gamma"] <- "gamma[1]"
    coda::mcmc(values, start = stats::start(chain), thin = coda::thin(chain))
  })
  do.call(coda::mcmc.list, chains)
}

.dbt_parameter_columns <- function(draws, prefix) {
  columns <- grep(paste0("^", prefix, "\\[[0-9]+\\]$"), colnames(draws), value = TRUE)
  if (!length(columns)) stop("No ", prefix, " draws were found.", call. = FALSE)
  index <- as.integer(sub(paste0("^", prefix, "\\[([0-9]+)\\]$"), "\\1", columns))
  columns[order(index)]
}

design_by_treatment_information_gate <- function(samples, gamma_prior_sd,
                                                 variance_ratio_threshold = design_by_treatment_gate_policy()$posterior_to_prior_variance_ratio,
                                                 tolerance = 1e-10) {
  gamma_prior_sd <- .dbt_validate_prior_sd(gamma_prior_sd)
  if (!identical(as.numeric(variance_ratio_threshold), 0.64)) {
    stop("The prespecified information-gate variance ratio must remain 0.64.", call. = FALSE)
  }
  draws <- .dbt_sample_matrix(samples)
  columns <- .dbt_parameter_columns(draws, "gamma")
  values <- draws[, columns, drop = FALSE]
  covariance <- stats::cov(values)
  if (length(columns) == 1L) covariance <- matrix(covariance, 1L, 1L, dimnames = list(columns, columns))
  finite <- all(is.finite(covariance))
  eigenvalues <- if (finite) eigen(covariance, symmetric = TRUE, only.values = TRUE)$values else rep(NA_real_, length(columns))
  nonsingular <- finite && all(eigenvalues > tolerance * max(1, max(eigenvalues)))
  maximum_ratio <- if (finite) max(eigenvalues) / gamma_prior_sd^2 else Inf
  passed <- nonsingular && is.finite(maximum_ratio) && maximum_ratio <= variance_ratio_threshold
  list(
    passed = passed,
    reliably_estimable = passed,
    classification = if (passed) "reliably_estimable" else "not_reliably_estimable",
    covariance = covariance,
    eigenvalues = eigenvalues,
    maximum_posterior_to_prior_variance_ratio = maximum_ratio,
    max_variance_ratio = maximum_ratio,
    variance_ratio_threshold = variance_ratio_threshold,
    posterior_sd_fraction_threshold = sqrt(variance_ratio_threshold),
    rationale = design_by_treatment_gate_policy()$information_rationale,
    prior_variance = gamma_prior_sd^2,
    finite_covariance = finite,
    nonsingular_covariance = nonsingular
  )
}

design_by_treatment_omnibus <- function(samples, probability = 0.95,
                                       tolerance = 1e-10) {
  draws <- .dbt_sample_matrix(samples)
  interaction_columns <- .dbt_parameter_columns(draws, "design_interaction")
  interactions <- draws[, interaction_columns, drop = FALSE]
  center <- apply(interactions, 2L, stats::median)
  scale <- apply(interactions, 2L, stats::sd)
  active <- is.finite(scale) & scale > tolerance
  structurally_zero <- !active & apply(abs(interactions), 2L, max) <= tolerance
  valid <- all(is.finite(interactions)) && all(active | structurally_zero) && any(active)
  if (valid) {
    standardized <- sweep(interactions[, active, drop = FALSE], 2L, center[active], "-")
    standardized <- sweep(standardized, 2L, scale[active], "/")
    maximum_deviation <- apply(abs(standardized), 1L, max)
    simultaneous_critical_value <- unname(stats::quantile(
      maximum_deviation, probability, names = FALSE, type = 8
    ))
  } else {
    simultaneous_critical_value <- NA_real_
  }
  lower <- center - simultaneous_critical_value * scale
  upper <- center + simultaneous_critical_value * scale
  lower[structurally_zero] <- 0
  upper[structurally_zero] <- 0
  excludes_zero <- valid & (lower > 0 | upper < 0)

  # Retain the covariance-ellipsoid calculation descriptively, but do not use
  # its normal approximation as the reporting classification.
  gamma_columns <- .dbt_parameter_columns(draws, "gamma")
  gamma <- draws[, gamma_columns, drop = FALSE]
  gamma_center <- colMeans(gamma)
  gamma_covariance <- stats::cov(gamma)
  if (length(gamma_columns) == 1L) gamma_covariance <- matrix(gamma_covariance, 1L, 1L)
  eigenvalues <- if (all(is.finite(gamma_covariance))) {
    eigen(gamma_covariance, symmetric = TRUE, only.values = TRUE)$values
  } else rep(NA_real_, length(gamma_columns))
  ellipsoid_valid <- all(is.finite(gamma_center)) && all(is.finite(eigenvalues)) &&
    all(eigenvalues > tolerance * max(1, max(eigenvalues)))
  distance <- if (ellipsoid_valid) {
    as.numeric(t(gamma_center) %*% solve(gamma_covariance, gamma_center))
  } else NA_real_
  ellipsoid_threshold <- stats::qchisq(probability, df = length(gamma_columns))
  list(
    valid = valid,
    method = "empirical simultaneous posterior band over reported design interactions",
    probability = probability,
    dimensions = length(interaction_columns),
    center = center,
    posterior_sd = scale,
    simultaneous_critical_value = simultaneous_critical_value,
    simultaneous_lower = lower,
    simultaneous_upper = upper,
    interaction_excludes_zero = excludes_zero,
    all_interactions_include_zero = valid && !any(excludes_zero),
    classification = if (!valid) {
      "inconclusive"
    } else if (any(excludes_zero)) {
      "evidence_of_inconsistency"
    } else {
      "no_evidence_of_inconsistency"
    },
    descriptive_ellipsoid = list(
      valid = ellipsoid_valid,
      squared_mahalanobis_at_zero = distance,
      threshold = ellipsoid_threshold,
      zero_inside_ellipsoid = ellipsoid_valid && distance <= ellipsoid_threshold,
      classification = if (!ellipsoid_valid) "inconclusive" else if (distance > ellipsoid_threshold) "zero_excluded" else "zero_included",
      reporting_basis = FALSE
    )
  )
}

summarize_design_interactions <- function(samples, probability = 0.95) {
  draws <- .dbt_sample_matrix(samples)
  columns <- .dbt_parameter_columns(draws, "design_interaction")
  alpha <- (1 - probability) / 2
  do.call(rbind, lapply(seq_along(columns), function(i) {
    values <- draws[, columns[i]]
    interval <- stats::quantile(values, c(0.5, alpha, 1 - alpha), names = FALSE)
    classification <- if (interval[2L] > 0) "positive" else if (interval[3L] < 0) "negative" else "zero_compatible"
    data.frame(
      interaction = columns[i],
      median = interval[1L],
      lower = interval[2L],
      upper = interval[3L],
      sd = stats::sd(values),
      zero_sign_classification = classification,
      stringsAsFactors = FALSE
    )
  }))
}

.dbt_fit_converged <- function(x) {
  if (is.list(x) && !inherits(x, "mcmc.list") && !is.null(x$diagnostics$converged)) {
    return(isTRUE(x$diagnostics$converged))
  }
  TRUE
}

.dbt_fit_prior_sd <- function(x, fallback) {
  if (is.list(x) && !inherits(x, "mcmc.list") && !is.null(x$metadata$gamma_prior_sd)) {
    return(as.numeric(x$metadata$gamma_prior_sd))
  }
  fallback
}

.dbt_comparison_contract_checks <- function(main, bounded) {
  main_contract <- if (is.list(main) && !inherits(main, "mcmc.list")) main$comparison_contract else NULL
  bounded_contract <- if (is.list(bounded) && !inherits(bounded, "mcmc.list")) bounded$comparison_contract else NULL
  required <- c(
    "analysis_scenario", "prepared_arms", "treatment_names", "treatment_map",
    "A", "H", "arm_mapping", "interaction_mapping", "config", "model_fingerprint"
  )
  complete <- is.list(main_contract) && is.list(bounded_contract) &&
    all(required %in% names(main_contract)) && all(required %in% names(bounded_contract))
  if (!complete) {
    return(c(list(contract_complete = FALSE),
             setNames(as.list(rep(FALSE, length(required))), paste0("same_", required))))
  }
  checks <- lapply(required, function(field) identical(main_contract[[field]], bounded_contract[[field]]))
  names(checks) <- paste0("same_", required)
  c(list(contract_complete = TRUE), checks)
}

# Compare the main SD=2 diagnostic with the bounded SD=1 check. A reportable
# conclusion requires convergence, information, identical omnibus and local
# classifications, and small posterior-median movement under the tighter prior.
compare_design_by_treatment_priors <- function(main, bounded,
                                               main_prior_sd = 2,
                                               bounded_prior_sd = 1,
                                               maximum_standardized_shift = design_by_treatment_gate_policy()$maximum_standardized_median_shift) {
  if (!identical(as.numeric(maximum_standardized_shift), 0.5)) {
    stop("The prespecified stability limit must remain 0.5 pooled posterior SD.", call. = FALSE)
  }
  main_sd <- .dbt_fit_prior_sd(main, main_prior_sd)
  bounded_sd <- .dbt_fit_prior_sd(bounded, bounded_prior_sd)
  if (main_sd != 2 || bounded_sd != 1) {
    stop("Prior stability requires the main SD=2 fit and bounded SD=1 fit in that order.", call. = FALSE)
  }
  main_info <- design_by_treatment_information_gate(main, main_sd)
  bounded_info <- design_by_treatment_information_gate(bounded, bounded_sd)
  contract_checks <- .dbt_comparison_contract_checks(main, bounded)
  main_omnibus <- design_by_treatment_omnibus(main)
  bounded_omnibus <- design_by_treatment_omnibus(bounded)
  main_interactions <- summarize_design_interactions(main)
  bounded_interactions <- summarize_design_interactions(bounded)
  if (!identical(main_interactions$interaction, bounded_interactions$interaction)) {
    stop("Main and bounded fits report different design interactions.", call. = FALSE)
  }
  pooled_sd <- sqrt((main_interactions$sd^2 + bounded_interactions$sd^2) / 2)
  shifts <- abs(main_interactions$median - bounded_interactions$median) / pooled_sd
  shifts[pooled_sd == 0 & main_interactions$median == bounded_interactions$median] <- 0
  interaction_comparison <- data.frame(
    interaction = main_interactions$interaction,
    main_classification = main_interactions$zero_sign_classification,
    bounded_classification = bounded_interactions$zero_sign_classification,
    main_median = main_interactions$median,
    bounded_median = bounded_interactions$median,
    pooled_posterior_sd = pooled_sd,
    standardized_median_shift = shifts,
    same_classification = main_interactions$zero_sign_classification == bounded_interactions$zero_sign_classification,
    shift_pass = is.finite(shifts) & shifts <= maximum_standardized_shift,
    stringsAsFactors = FALSE
  )
  gates <- list(
    main_converged = .dbt_fit_converged(main),
    bounded_converged = .dbt_fit_converged(bounded),
    main_information = main_info$passed,
    bounded_information = bounded_info$passed,
    paired_run_contract = all(unlist(contract_checks, use.names = FALSE)),
    same_omnibus_classification = identical(main_omnibus$classification, bounded_omnibus$classification) &&
      main_omnibus$classification != "inconclusive",
    same_interaction_classifications = all(interaction_comparison$same_classification),
    median_shifts_within_limit = all(interaction_comparison$shift_pass)
  )
  passed <- all(unlist(gates, use.names = FALSE))
  conclusion <- if (!passed) {
    "inconclusive"
  } else if (main_omnibus$classification == "evidence_of_inconsistency") {
    "evidence_of_inconsistency"
  } else {
    "no_evidence_of_inconsistency"
  }
  reporting_statement <- if (!passed) {
    "Inconsistency could not be reliably estimated."
  } else if (conclusion == "evidence_of_inconsistency") {
    "Evidence of inconsistency under the stated model and priors."
  } else {
    "No evidence of inconsistency under the stated model and priors."
  }
  list(
    passed = passed,
    stable = passed,
    reliably_estimable = passed,
    conclusion = conclusion,
    reporting_statement = reporting_statement,
    gates = gates,
    contract_checks = contract_checks,
    main_information = main_info,
    bounded_information = bounded_info,
    main_omnibus = main_omnibus,
    bounded_omnibus = bounded_omnibus,
    interaction_comparison = interaction_comparison,
    maximum_standardized_shift = maximum_standardized_shift,
    stability_rationale = design_by_treatment_gate_policy()$stability_rationale,
    gate_policy = design_by_treatment_gate_policy()
  )
}

.dbt_combined_diagnostic_table <- function(primary, inconsistency) {
  primary_table <- primary$table
  primary_table$diagnostic_group <- "primary_reported_parameter"
  primary_table$min_chain_effective_size <- NA_real_
  inconsistency_table <- inconsistency$table
  inconsistency_table$diagnostic_group <- "inconsistency_parameter"
  inconsistency_table$effective_size <- NA_real_
  columns <- c(
    "diagnostic_group", "parameter", "parameter_type", "psrf_point",
    "psrf_upper_95", "effective_size", "min_chain_effective_size",
    "psrf_pass", "ess_pass", "passed"
  )
  rbind(primary_table[, columns], inconsistency_table[, columns])
}

fit_design_by_treatment_model <- function(
  scenario_name,
  raw_rows,
  treatment_map,
  treatment_names,
  analysis_scenario = NULL,
  gamma_prior_sd = 2,
  model_file = file.path("R", "model_fixed_study", "design_by_treatment.jags"),
  config = fixed_study_production_config(),
  chain_seeds = NULL,
  master_seed = 12345L,
  scheduled_run_names = scenario_name,
  preflight = FALSE,
  test_only = FALSE,
  quiet = TRUE,
  output_controls = NULL
) {
  .dbt_require_functions(c(
    "fixed_study_production_config", "validate_fixed_study_config",
    "fixed_study_chain_seeds", "prepare_fixed_study_data",
    "augment_fixed_study_samples", "summarize_fixed_study_posterior",
    "fixed_study_diagnostics", "extract_fixed_study_p_draws",
    "fixed_study_residual_deviance", "fixed_study_descriptive_dic"
  ))
  if (length(scenario_name) != 1L || is.na(scenario_name) || !nzchar(scenario_name)) {
    stop("scenario_name must be one non-empty string.", call. = FALSE)
  }
  if (is.null(analysis_scenario)) analysis_scenario <- sub("_sd[12]$", "", scenario_name)
  if (length(analysis_scenario) != 1L || is.na(analysis_scenario) || !nzchar(analysis_scenario)) {
    stop("analysis_scenario must be one non-empty paired-run scenario identifier.", call. = FALSE)
  }
  gamma_prior_sd <- .dbt_validate_prior_sd(gamma_prior_sd)
  config <- validate_fixed_study_config(config, test_only = test_only)
  seeds <- fixed_study_chain_seeds(scenario_name, chain_seeds, scheduled_run_names, master_seed)
  # Reject rather than silently alter the final diagnostic network.
  prepared <- prepare_fixed_study_data(raw_rows, treatment_map, treatment_names, single_arm = "reject")
  design <- prepare_design_by_treatment_data(prepared, gamma_prior_sd)
  if (!file.exists(model_file)) stop("Design-by-treatment JAGS model is unavailable: ", model_file, call. = FALSE)
  normalized_model_file <- normalizePath(model_file, mustWork = TRUE)
  model_fingerprint <- .dbt_sha256_file(normalized_model_file)
  comparison_contract <- list(
    analysis_scenario = as.character(analysis_scenario),
    prepared_arms = prepared$arms,
    treatment_names = prepared$treatment_names,
    treatment_map = prepared$treatment_map,
    A = design$geometry$A,
    H = design$H,
    arm_mapping = design$arm_mapping,
    interaction_mapping = design$interaction_mapping,
    config = config,
    model_fingerprint = model_fingerprint
  )

  controls <- list(progress_bar = if (isTRUE(quiet)) "none" else "text", messages = !isTRUE(quiet))
  if (!is.null(output_controls)) {
    if (!is.list(output_controls) || !all(names(output_controls) %in% names(controls))) {
      stop("output_controls may contain only progress_bar and messages.", call. = FALSE)
    }
    controls[names(output_controls)] <- output_controls
  }
  if (!controls$progress_bar %in% c("none", "text", "gui") ||
      length(controls$messages) != 1L || !is.logical(controls$messages) || is.na(controls$messages)) {
    stop("Invalid output_controls.", call. = FALSE)
  }
  preflight_result <- list(
    scenario_name = scenario_name,
    analysis_scenario = analysis_scenario,
    model_file = normalized_model_file,
    model_fingerprint = model_fingerprint,
    config = config,
    chain_seeds = seeds,
    prepared_data = prepared,
    design = design,
    comparison_contract = comparison_contract,
    test_only = isTRUE(test_only),
    preflight = TRUE
  )
  if (isTRUE(preflight)) return(preflight_result)
  if (!requireNamespace("rjags", quietly = TRUE) || !requireNamespace("coda", quietly = TRUE)) {
    stop("Packages 'rjags' and 'coda' are required to fit the design-by-treatment model.", call. = FALSE)
  }

  if (isTRUE(controls$messages)) message("Fitting fixed design-by-treatment model for ", scenario_name, ".")
  inits <- lapply(seeds, function(seed) {
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)
  })
  model <- rjags::jags.model(
    file = model_file,
    data = design$jags_data,
    inits = inits,
    n.chains = config$n.chains,
    n.adapt = config$n.adapt,
    quiet = isTRUE(quiet)
  )
  stats::update(model, n.iter = config$n.burnin, progress.bar = controls$progress_bar)
  monitored <- c("gamma", "design_interaction", "d", "tau_arm", "tau_contrast", "p")
  raw_samples <- rjags::coda.samples(
    model = model,
    variable.names = monitored,
    n.iter = config$n.iter,
    thin = config$n.thin,
    progress.bar = controls$progress_bar
  )
  raw_samples <- .dbt_normalize_scalar_gamma_name(raw_samples)
  samples <- augment_fixed_study_samples(raw_samples, length(treatment_names))
  posterior <- summarize_fixed_study_posterior(samples, treatment_names)
  inconsistency_diagnostics <- design_by_treatment_diagnostics(samples, design$geometry$q, design$geometry$m)
  primary_diagnostics <- fixed_study_diagnostics(samples, length(treatment_names))
  diagnostics <- list(
    primary = primary_diagnostics,
    inconsistency = inconsistency_diagnostics,
    converged = primary_diagnostics$converged && inconsistency_diagnostics$converged,
    psrf_threshold = 1.05,
    ess_threshold = 1000
  )
  diagnostic_table <- .dbt_combined_diagnostic_table(primary_diagnostics, inconsistency_diagnostics)
  diagnostics$table <- diagnostic_table
  diagnostics$max_psrf_upper_95 <- max(
    primary_diagnostics$max_psrf_upper_95,
    inconsistency_diagnostics$max_psrf_upper_95,
    na.rm = TRUE
  )
  diagnostics$min_effective_size <- min(
    primary_diagnostics$min_effective_size,
    inconsistency_diagnostics$min_chain_effective_size,
    na.rm = TRUE
  )
  information <- design_by_treatment_information_gate(samples, gamma_prior_sd)
  omnibus <- design_by_treatment_omnibus(samples)
  interactions <- cbind(design$interaction_mapping, summarize_design_interactions(samples))
  p_draws <- extract_fixed_study_p_draws(samples, prepared$jags_data$Narm)
  residual <- fixed_study_residual_deviance(p_draws, prepared$arms$R, prepared$arms$N)
  dic <- fixed_study_descriptive_dic(p_draws, prepared$arms$R, prepared$arms$N)
  initial_gate <- diagnostics$converged && information$passed && omnibus$valid
  warnings <- character()
  if (!primary_diagnostics$converged) {
    warnings <- c(warnings, "Primary reported d, pairwise, or tau parameters failed PSRF/ESS thresholds.")
  }
  if (!inconsistency_diagnostics$converged) {
    warnings <- c(warnings, "Gamma or reported design interactions failed PSRF/chain-wise ESS thresholds.")
  }
  if (!information$passed) warnings <- c(warnings, "Inconsistency directions failed the information gate.")
  if (!omnibus$valid) warnings <- c(warnings, "The empirical simultaneous posterior-band check was invalid.")

  list(
    scenario_name = scenario_name,
    analysis_scenario = analysis_scenario,
    model_type = "fixed_design_by_treatment_arm_heterogeneity",
    samples = samples,
    monitored_samples = raw_samples,
    posterior = posterior,
    treatment_effects = posterior$treatment_effects,
    pairwise_effects = posterior$pairwise_effects,
    heterogeneity = posterior$heterogeneity,
    design_interactions = interactions,
    diagnostics = diagnostics,
    primary_diagnostics = primary_diagnostics,
    inconsistency_diagnostics = inconsistency_diagnostics,
    diagnostic_table = diagnostic_table,
    information_gate = information,
    omnibus = omnibus,
    gate_policy = design_by_treatment_gate_policy(),
    reporting_status = if (initial_gate) "awaiting_bounded_prior_stability_check" else "inconclusive",
    reporting_statement = if (initial_gate) {
      "Initial gates passed; a bounded-prior stability check is still required."
    } else {
      "Inconsistency could not be reliably estimated."
    },
    residual_deviance = residual,
    residual_deviance_arms = residual$arm_summary,
    residual_deviance_summary = residual$total_summary,
    dic = dic,
    warnings = warnings,
    prepared_data = prepared,
    design = design,
    comparison_contract = comparison_contract,
    metadata = list(
      scenario_name = scenario_name,
      analysis_scenario = analysis_scenario,
      model_type = "fixed_design_by_treatment_arm_heterogeneity",
      gamma_prior_sd = gamma_prior_sd,
      model_fingerprint = model_fingerprint,
      n_inconsistency_df = design$geometry$q,
      n_design_contrasts = design$geometry$m,
      n_adapt = config$n.adapt,
      n_burnin = config$n.burnin,
      n_iter = config$n.iter,
      n_thin = config$n.thin,
      n_chains = config$n.chains,
      chain_seeds = seeds,
      monitored_parameters = paste(monitored, collapse = "|"),
      production_mcmc_contract = identical(config, fixed_study_production_config()),
      max_psrf_upper_95 = diagnostics$max_psrf_upper_95,
      min_effective_size = diagnostics$min_effective_size,
      converged = diagnostics$converged,
      information_gate_passed = information$passed,
      information_variance_ratio_threshold = information$variance_ratio_threshold,
      information_posterior_sd_fraction_threshold = information$posterior_sd_fraction_threshold,
      information_gate_rationale = information$rationale,
      stability_maximum_standardized_median_shift = design_by_treatment_gate_policy()$maximum_standardized_median_shift,
      stability_gate_rationale = design_by_treatment_gate_policy()$stability_rationale,
      residual_deviance_mean = residual$total_summary$mean,
      residual_deviance_per_arm = residual$total_summary$per_arm_mean,
      dic = dic$DIC,
      dic_descriptive_only = TRUE,
      no_fallback_or_extension = TRUE
    )
  )
}
