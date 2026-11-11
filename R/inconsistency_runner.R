# =============================================================================
# Immutable runner for the paired fixed-study design-by-treatment diagnostic
# =============================================================================

inconsistency_runs <- function() c("inconsistency_sd2", "inconsistency_sd1")

inconsistency_sd <- function(scenario_name) {
  values <- c(inconsistency_sd2 = 2, inconsistency_sd1 = 1)
  if (length(scenario_name) != 1L || is.na(scenario_name) ||
      !scenario_name %in% names(values)) {
    stop("Only inconsistency_sd2 or inconsistency_sd1 may be run by this runner.", call. = FALSE)
  }
  unname(values[[scenario_name]])
}

inconsistency_runtime_hashes <- function(project_root) {
  core <- core_runtime_hashes(project_root)
  sha <- core_required("sha256_file")
  paths <- c(
    inconsistency_runner = file.path(project_root, "R", "inconsistency_runner.R"),
    design_by_treatment = file.path(project_root, "R", "design_by_treatment.R"),
    design_by_treatment_jags = file.path(project_root, "R", "model_fixed_study", "design_by_treatment.jags")
  )
  rbind(core, data.frame(
    item = names(paths), path = unname(paths),
    sha256 = vapply(paths, sha, character(1)), stringsAsFactors = FALSE
  ))
}

inconsistency_seed_schedule <- function(output_root, scenario_name) {
  schedule_path <- file.path(output_root, "core_seed_schedule.csv")
  schedule <- read.csv(schedule_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("scenario_name", "chain", "seed")
  if (!all(required %in% names(schedule))) stop("Committed seed schedule is incomplete.", call. = FALSE)
  runs <- inconsistency_runs()
  if (!all(runs %in% schedule$scenario_name)) {
    stop("Committed seed schedule does not contain both inconsistency runs.", call. = FALSE)
  }
  rows <- schedule[as.character(schedule$scenario_name) == scenario_name, required, drop = FALSE]
  rows <- rows[order(rows$chain), , drop = FALSE]
  if (!identical(as.integer(rows$chain), 1:4) || anyNA(rows$seed) ||
      anyDuplicated(rows$seed) || any(as.integer(rows$seed) <= 0L)) {
    stop("Committed seed schedule must provide four unique positive chains for the selected inconsistency run.", call. = FALSE)
  }
  rows$chain <- as.integer(rows$chain)
  rows$seed <- as.integer(rows$seed)
  rows
}

inconsistency_verify <- function(project_root, output_root, scenario_name,
                                               input_manifest_sha256, checkpoint_commit,
                                               check_git = TRUE) {
  gamma_prior_sd <- inconsistency_sd(scenario_name)
  # The core verifier remains the single authority for the committed checkpoint,
  # primary raw input, plan, model, geometry, and protected artifact ledger.
  verified <- core_verify_checkpoint(
    project_root, output_root, "primary", input_manifest_sha256,
    checkpoint_commit, check_git
  )
  prepared <- core_assert_raw_fit_contract(verified)
  if (length(prepared$dropped_single_arm_studies)) {
    stop("The manifested fixed-study primary input has dropped single-arm studies and cannot support the diagnostic.", call. = FALSE)
  }
  seeds <- inconsistency_seed_schedule(output_root, scenario_name)
  list(
    verified = verified,
    prepared = prepared,
    seeds = seeds,
    gamma_prior_sd = gamma_prior_sd,
    runtime_hashes = inconsistency_runtime_hashes(project_root)
  )
}

inconsistency_next_attempt <- function(run_root, scenario_name, rerun_reason = "") {
  core_next_attempt(run_root, scenario_name, rerun_reason)
}

.inconsistency_write_csv <- function(value, path) {
  if (is.null(value)) value <- data.frame()
  write.csv(value, path, row.names = FALSE)
}

inconsistency_write_warnings <- function(runtime_warnings, gate_warnings,
                                                       attempt_path) {
  .inconsistency_write_csv(data.frame(
    warning = unique(as.character(runtime_warnings)), stringsAsFactors = FALSE
  ), file.path(attempt_path, "r_jags_warnings.csv"))
  .inconsistency_write_csv(data.frame(
    warning = unique(as.character(gate_warnings)), stringsAsFactors = FALSE
  ), file.path(attempt_path, "gate_warnings.csv"))
  .inconsistency_write_csv(data.frame(
    source = c(rep("R_or_JAGS", length(runtime_warnings)), rep("gate", length(gate_warnings))),
    warning = c(as.character(runtime_warnings), as.character(gate_warnings)),
    stringsAsFactors = FALSE
  ), file.path(attempt_path, "warnings.csv"))
}

inconsistency_write_artifacts <- function(result, attempt_path,
                                                        runtime_warnings) {
  .inconsistency_write_csv(result$treatment_effects, file.path(attempt_path, "posterior_treatment_effects.csv"))
  .inconsistency_write_csv(result$pairwise_effects, file.path(attempt_path, "posterior_pairwise_effects.csv"))
  .inconsistency_write_csv(result$heterogeneity, file.path(attempt_path, "posterior_heterogeneity.csv"))
  .inconsistency_write_csv(result$design_interactions, file.path(attempt_path, "design_interactions.csv"))
  .inconsistency_write_csv(result$diagnostic_table, file.path(attempt_path, "diagnostics.csv"))
  saveRDS(result$diagnostics, file.path(attempt_path, "diagnostics.rds"))
  saveRDS(result$primary_diagnostics, file.path(attempt_path, "primary_diagnostics.rds"))
  saveRDS(result$inconsistency_diagnostics, file.path(attempt_path, "inconsistency_diagnostics.rds"))
  saveRDS(result$information_gate, file.path(attempt_path, "information_gate.rds"))
  saveRDS(result$omnibus, file.path(attempt_path, "omnibus.rds"))
  .inconsistency_write_csv(data.frame(
    passed = result$information_gate$passed,
    classification = result$information_gate$classification,
    maximum_posterior_to_prior_variance_ratio = result$information_gate$maximum_posterior_to_prior_variance_ratio,
    variance_ratio_threshold = result$information_gate$variance_ratio_threshold,
    stringsAsFactors = FALSE
  ), file.path(attempt_path, "information_gate.csv"))
  .inconsistency_write_csv(data.frame(
    valid = result$omnibus$valid, classification = result$omnibus$classification,
    method = result$omnibus$method, probability = result$omnibus$probability,
    simultaneous_critical_value = result$omnibus$simultaneous_critical_value,
    stringsAsFactors = FALSE
  ), file.path(attempt_path, "omnibus.csv"))
  .inconsistency_write_csv(result$residual_deviance_arms, file.path(attempt_path, "residual_deviance_arms.csv"))
  .inconsistency_write_csv(result$residual_deviance_summary, file.path(attempt_path, "residual_deviance_summary.csv"))
  .inconsistency_write_csv(as.data.frame(result$dic, stringsAsFactors = FALSE), file.path(attempt_path, "descriptive_dic.csv"))
  inconsistency_write_warnings(runtime_warnings, result$warnings, attempt_path)
  saveRDS(result$samples, file.path(attempt_path, "design_by_treatment_mcmc_samples.rds"))
  invisible(TRUE)
}

inconsistency_metadata <- function(attempt, scenario_name, verified, check,
                                                analysis_commit, checkpoint_commit, docker_image,
                                                input_manifest_sha256, config, rerun_reason,
                                                status = "started", started_utc,
                                                completed_utc = NA_character_, posterior_draws = NA_integer_,
                                                combined_converged = NA, information_passed = NA,
                                                omnibus_valid = NA, code_hashes_unchanged = NA,
                                                protected_artifacts_unchanged = NA,
                                                deviation = "none", result = NULL) {
  data.frame(
    scenario = scenario_name, analysis_scenario = "primary",
    definition = sprintf("Fixed design-by-treatment diagnostic with gamma prior SD=%d", check$gamma_prior_sd),
    gamma_prior_sd = check$gamma_prior_sd, attempt = attempt$attempt,
    attempt_path = attempt$relative_path, status = status, started_utc = started_utc,
    completed_utc = completed_utc, analysis_commit = analysis_commit,
    checkpoint_commit = checkpoint_commit,
    checkpoint_input_builder_commit = verified$metadata$analysis_commit,
    docker_image = docker_image, input_manifest_sha256 = input_manifest_sha256,
    raw_input_sha256 = verified$manifest$raw_input_sha256,
    nma_input_sha256 = verified$manifest$nma_input_sha256,
    n_adapt = config$n.adapt, n_burnin = config$n.burnin, n_iter = config$n.iter,
    n_thin = config$n.thin, n_chains = config$n.chains,
    expected_posterior_draws = as.integer(config$n.iter / config$n.thin * config$n.chains),
    posterior_draws = posterior_draws,
    chain_seed_1 = check$seeds$seed[1], chain_seed_2 = check$seeds$seed[2],
    chain_seed_3 = check$seeds$seed[3], chain_seed_4 = check$seeds$seed[4],
    rerun_reason = rerun_reason, test_only = isTRUE(result$test_only),
    production_mcmc_contract = identical(config, fixed_study_production_config()),
    combined_converged = combined_converged, information_gate_passed = information_passed,
    empirical_omnibus_valid = omnibus_valid, code_hashes_unchanged = code_hashes_unchanged,
    protected_artifacts_unchanged = protected_artifacts_unchanged,
    max_psrf_upper_95 = if (is.null(result)) NA_real_ else result$diagnostics$max_psrf_upper_95,
    min_effective_size = if (is.null(result)) NA_real_ else result$diagnostics$min_effective_size,
    residual_deviance = if (is.null(result)) NA_real_ else result$residual_deviance_summary$mean,
    residual_deviance_per_arm = if (is.null(result)) NA_real_ else result$residual_deviance_summary$per_arm_mean,
    descriptive_dic = if (is.null(result)) NA_real_ else result$dic$DIC,
    no_fallback_or_extension = TRUE, paired_prior_stability_is_not_single_run_gate = TRUE,
    deviation = deviation, stringsAsFactors = FALSE
  )
}

run_inconsistency <- function(
  project_root, scenario_name, input_manifest_sha256, analysis_commit, checkpoint_commit,
  docker_image, output_root = file.path(project_root, "results", "analysis", "core"),
  rerun_reason = "", preflight = FALSE, config = fixed_study_production_config(),
  test_only = FALSE, fit_function = fit_design_by_treatment_model, check_git = TRUE
) {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  output_root <- normalizePath(output_root, mustWork = TRUE)
  if (!grepl("^[0-9a-f]{40}$", analysis_commit) || !grepl("^[0-9a-f]{40}$", checkpoint_commit)) {
    stop("ANALYSIS_COMMIT and CHECKPOINT_COMMIT must be full commit hashes.", call. = FALSE)
  }
  config <- validate_fixed_study_config(config, test_only = test_only)
  if (!isTRUE(test_only) && !identical(config, fixed_study_production_config())) {
    stop("Only the locked four-chain production MCMC contract may be run.", call. = FALSE)
  }
  check <- inconsistency_verify(
    project_root, output_root, scenario_name, input_manifest_sha256, checkpoint_commit, check_git
  )
  if (isTRUE(preflight)) return(invisible(c(check, list(config = config))))

  attempt <- inconsistency_next_attempt(publication_run_root(output_root), scenario_name, rerun_reason)
  dir.create(attempt$path, recursive = TRUE, showWarnings = FALSE)
  attempt$relative_path <- publication_relative_attempt(attempt$path, output_root)
  started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  metadata <- inconsistency_metadata(
    attempt, scenario_name, check$verified, check, analysis_commit, checkpoint_commit,
    docker_image, input_manifest_sha256, config, rerun_reason, started_utc = started
  )
  publication_write_metadata(metadata, file.path(attempt$path, "run_metadata.csv"))
  write.csv(check$seeds, file.path(attempt$path, "seed_schedule_used.csv"), row.names = FALSE)
  write.csv(check$runtime_hashes, file.path(attempt$path, "runtime_code_hashes_before.csv"), row.names = FALSE)
  write.csv(check$verified$protected_before, file.path(attempt$path, "protected_artifacts_before.csv"), row.names = FALSE)
  write.csv(check$prepared$arms, file.path(attempt$path, "fitted_analysis_arms.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(attempt$path, "sessionInfo.txt"))
  cat(paste(started, scenario_name, sprintf("attempt_%03d", attempt$attempt), analysis_commit), "\n",
      file = file.path(dirname(publication_run_root(output_root)), "MCMC_STARTED"), sep = "", append = TRUE)

  runtime_warnings <- character()
  result <- tryCatch(withCallingHandlers(
    fit_function(
      scenario_name = scenario_name, analysis_scenario = "primary",
      raw_rows = read.csv(check$verified$raw_path, stringsAsFactors = FALSE, check.names = FALSE),
      treatment_map = core_treatment_map(), treatment_names = core_treatment_names(),
      gamma_prior_sd = check$gamma_prior_sd,
      model_file = file.path(project_root, "R", "model_fixed_study", "design_by_treatment.jags"),
      config = config, chain_seeds = check$seeds, master_seed = 12345L,
      scheduled_run_names = inconsistency_runs(), preflight = FALSE,
      test_only = test_only, quiet = TRUE
    ), warning = function(w) {
      runtime_warnings <<- c(runtime_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
  completed <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  failed <- inherits(result, "error")
  failure <- if (failed) conditionMessage(result) else NA_character_
  if (failed) inconsistency_write_warnings(runtime_warnings, character(), attempt$path)
  if (!failed) {
    write_error <- tryCatch({
      inconsistency_write_artifacts(result, attempt$path, runtime_warnings)
      NULL
    }, error = identity)
    if (inherits(write_error, "error")) { failed <- TRUE; failure <- conditionMessage(write_error) }
  }
  protected_after <- tryCatch(core_collect_protected(), error = identity)
  protected_ok <- FALSE
  if (inherits(protected_after, "error")) {
    failed <- TRUE; failure <- conditionMessage(protected_after)
  } else {
    protected_ok <- all(protected_after$matches_expected)
    write.csv(protected_after, file.path(attempt$path, "protected_artifacts_after.csv"), row.names = FALSE)
    if (!protected_ok) { failed <- TRUE; failure <- "Protected artifact checksum mismatch after run." }
  }
  hashes_after <- tryCatch(inconsistency_runtime_hashes(project_root), error = identity)
  hashes_ok <- FALSE
  if (inherits(hashes_after, "error")) {
    failed <- TRUE; failure <- conditionMessage(hashes_after)
  } else {
    hashes_ok <- identical(check$runtime_hashes$item, hashes_after$item) &&
      identical(check$runtime_hashes$sha256, hashes_after$sha256)
    write.csv(hashes_after, file.path(attempt$path, "runtime_code_hashes_after.csv"), row.names = FALSE)
    if (!hashes_ok) { failed <- TRUE; failure <- "Core or design-by-treatment code changed during MCMC." }
  }
  draws <- NA_integer_; combined_converged <- FALSE; information_passed <- FALSE; omnibus_valid <- FALSE
  if (!failed) {
    draws <- nrow(as.matrix(result$samples))
    combined_converged <- isTRUE(result$diagnostics$converged)
    information_passed <- isTRUE(result$information_gate$passed)
    omnibus_valid <- isTRUE(result$omnibus$valid)
    if (!combined_converged) failure <- "Combined primary and inconsistency diagnostics did not meet the PSRF/ESS contract."
    else if (!information_passed) failure <- "The design-by-treatment information gate did not pass."
    else if (!omnibus_valid) failure <- "The empirical inconsistency omnibus result was invalid."
    else if (!identical(draws, as.integer(config$n.iter / config$n.thin * config$n.chains))) failure <- "Retained posterior draw count did not match the locked contract."
  }
  status <- if (failed) "failed" else if (is.na(failure)) "complete" else "nonconverged"
  metadata <- inconsistency_metadata(
    attempt, scenario_name, check$verified, check, analysis_commit, checkpoint_commit,
    docker_image, input_manifest_sha256, config, rerun_reason, status, started, completed,
    draws, combined_converged, information_passed, omnibus_valid, hashes_ok, protected_ok,
    if (is.na(failure)) "none" else failure, if (failed) NULL else result
  )
  publication_write_metadata(metadata, file.path(attempt$path, "run_metadata.csv"))
  marker <- if (identical(status, "complete")) "RUN_COMPLETE" else if (identical(status, "nonconverged")) "RUN_NONCONVERGED" else "RUN_FAILED"
  writeLines(if (is.na(failure)) completed else failure, file.path(attempt$path, marker))
  core_write_artifact_index(attempt$path)
  if (!identical(status, "complete")) {
    stop("Fixed-study inconsistency diagnostic did not meet the execution contract: ", scenario_name, ". ", failure, call. = FALSE)
  }
  invisible(list(attempt = attempt, metadata = metadata, result = result))
}
