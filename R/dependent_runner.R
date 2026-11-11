# =============================================================================
# Reference runner for manuscript-facing fixed-study dependent fixed-study fits
# =============================================================================

dependent_pinned_docker_image <- function() {
  "runtime-recorded"
}

dependent_plan_file <- function(project_root) {
  core_plan_file(project_root)
}

dependent_checkpoint_files <- function() c(
  "VALIDATION_REPORT.md", "checkpoint_metadata.csv", "fixed_study_seed_schedule.csv",
  "protected_artifact_checksums.csv", "scenario_input_manifest.csv",
  "scenario_specification.csv", "scenario_row_deltas.csv", "validation_checks.csv"
)

dependent_required <- function(name) {
  if (!exists(name, mode = "function", inherits = TRUE)) {
    stop("Dependent runner dependency is unavailable: ", name, call. = FALSE)
  }
  get(name, mode = "function", inherits = TRUE)
}

dependent_runtime_hashes <- function(project_root) {
  sha <- dependent_required("sha256_file")
  paths <- c(
    public_runner = file.path(project_root, "R", "run_scenario.R"),
    public_utils = file.path(project_root, "R", "publication_utils.R"),
    data_utils = file.path(project_root, "R", "common", "data_utils.R"),
    input_manifest = dependent_plan_file(project_root),
    fixed_model = file.path(project_root, "R", "fixed_study_model.R"),
    fixed_model_jags = file.path(project_root, "R", "model_fixed_study", "model.jags"),
    scenario_builder = file.path(project_root, "R", "dependent_scenarios.R"),
    checkpoint_builder = file.path(project_root, "R", "dependent_checkpoint.R"),
    core_checkpoint = file.path(project_root, "R", "core_checkpoint.R"),
    core_runner = file.path(project_root, "R", "core_runner.R"),
    dependent_runner = file.path(project_root, "R", "dependent_runner.R")
  )
  if (any(!file.exists(paths))) {
    stop("Dependent runner code artifact(s) missing: ",
         paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
  }
  data.frame(
    item = names(paths), path = unname(paths), sha256 = vapply(paths, sha, character(1)),
    stringsAsFactors = FALSE
  )
}

dependent_verify_checkpoint_commit <- function(project_root, output_root, checkpoint_commit, check_git = TRUE) {
  publication_verify_inputs(project_root)
  invisible(TRUE)
}

dependent_treatment_contract <- function(manifest_row) {
  split_wt <- identical(as.character(manifest_row$analysis_family), "genotyping")
  map <- if (split_wt) {
    dependent_split_wt_map()
  } else {
    dependent_treatment_map()
  }
  names <- if (split_wt) {
    dependent_split_wt_names()
  } else {
    dependent_treatment_names()
  }
  expected_nodes <- if (split_wt) 6L else 5L
  manifested_names <- strsplit(as.character(manifest_row$treatment_names), "|", fixed = TRUE)[[1]]
  if (length(names) != expected_nodes || !identical(unname(as.character(names)), manifested_names) ||
      !setequal(as.integer(map), seq_len(expected_nodes))) {
    stop("Manifest treatment map/names do not satisfy the dependent fixed-study contract.", call. = FALSE)
  }
  list(map = map, names = names, split_wt = split_wt, expected_nodes = expected_nodes)
}

dependent_seed_schedule <- function(output_root, scenario_name) {
  schedule_path <- file.path(output_root, "fixed_study_seed_schedule.csv")
  schedule <- read.csv(schedule_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("scenario_name", "chain", "seed")
  if (!all(required %in% names(schedule)) || anyDuplicated(schedule$seed)) {
    stop("Committed dependent seed schedule is incomplete or non-unique.", call. = FALSE)
  }
  rows <- schedule[as.character(schedule$scenario_name) == scenario_name, required, drop = FALSE]
  rows <- rows[order(rows$chain), , drop = FALSE]
  if (!identical(as.integer(rows$chain), 1:4) || anyNA(rows$seed) ||
      anyDuplicated(rows$seed) || any(as.integer(rows$seed) <= 0L)) {
    stop("Committed dependent seed schedule must provide four unique positive chains for the selected scenario.", call. = FALSE)
  }
  rows$chain <- as.integer(rows$chain)
  rows$seed <- as.integer(rows$seed)
  rows
}

dependent_assert_raw_fit_contract <- function(verified) {
  raw <- read.csv(verified$raw_path, stringsAsFactors = FALSE, check.names = FALSE)
  contract <- verified$treatment_contract
  prepared <- prepare_fixed_study_data(raw, contract$map, contract$names, single_arm = "drop")
  expected <- c(
    fixed_model_studies = prepared$jags_data$Nstudy,
    fixed_model_arms = prepared$jags_data$Narm,
    fixed_model_n = sum(prepared$arms$N),
    fixed_model_r = sum(prepared$arms$R),
    dropped_single_arm_studies = paste(prepared$dropped_single_arm_studies, collapse = "|")
  )
  actual <- c(
    fixed_model_studies = verified$manifest$fixed_model_studies,
    fixed_model_arms = verified$manifest$fixed_model_arms,
    fixed_model_n = verified$manifest$fixed_model_n,
    fixed_model_r = verified$manifest$fixed_model_r,
    dropped_single_arm_studies = ifelse(is.na(verified$manifest$dropped_single_arm_studies), "",
                                        verified$manifest$dropped_single_arm_studies)
  )
  if (!identical(as.character(expected), as.character(actual))) {
    stop("Dependent raw input no longer produces the manifested fixed-study analysis input.", call. = FALSE)
  }
  nma <- read.csv(verified$nma_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(nma), c("Study", "s.id", "t.id", "r", "n"))) {
    stop("Committed dependent transformed input has an unexpected arm schema.", call. = FALSE)
  }
  recreated <- transform_to_nma_format(
    prepared$arms[, c("Study", "T", "N", "R"), drop = FALSE],
    contract$map, aggregate_wt = !contract$split_wt
  )
  recreated <- recreated[, c("Study", "s.id", "t.id", "r", "n"), drop = FALSE]
  comparable <- function(value) data.frame(
    Study = as.character(value$Study), s.id = as.integer(value[["s.id"]]),
    t.id = as.integer(value[["t.id"]]), r = as.integer(value$r), n = as.integer(value$n),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!identical(comparable(recreated), comparable(nma))) {
    stop("Prepared dependent fixed-study arms do not match the manifested transformed input row-for-row.", call. = FALSE)
  }
  if (!setequal(unique(as.integer(nma[["t.id"]])), seq_len(contract$expected_nodes))) {
    stop("Committed dependent transformed input does not contain every required treatment node.", call. = FALSE)
  }
  prepared
}

dependent_verify_checkpoint <- function(
  project_root, output_root, scenario_name, input_manifest_sha256,
  checkpoint_commit, parent_checkpoint_commit, docker_image, check_git = TRUE
) {
  sha <- dependent_required("sha256_file")
  missing <- dependent_checkpoint_files()[
    !file.exists(file.path(output_root, dependent_checkpoint_files()))
  ]
  if (length(missing)) {
    stop("Immutable dependent checkpoint artifact(s) missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  metadata <- core_read_one_row(file.path(output_root, "checkpoint_metadata.csv"),
                                         "Dependent checkpoint metadata")
  required <- c("analysis_commit", "docker_image", "input_manifest_sha256",
                "source_parent_manifest_sha256", "source_parent_primary_raw_sha256",
                "scenario_builder_sha256", "checkpoint_code_sha256")
  if (!all(required %in% names(metadata))) stop("Dependent checkpoint metadata is incomplete.", call. = FALSE)
  if (!identical(as.character(metadata$input_manifest_sha256), as.character(input_manifest_sha256))) {
    stop("Dependent input-manifest checksum differs from the locked contract.", call. = FALSE)
  }
  runtime_hashes <- dependent_runtime_hashes(project_root)
  observed <- setNames(runtime_hashes$sha256, runtime_hashes$item)
  expected <- c(
    input_manifest = metadata$input_manifest_sha256,
    scenario_builder = metadata$scenario_builder_sha256,
    checkpoint_builder = metadata$checkpoint_code_sha256
  )
  bad <- names(expected)[observed[names(expected)] != unname(expected)]
  if (length(bad)) stop("Dependent input-manifest/checkpoint code hash mismatch: ", paste(bad, collapse = ", "), call. = FALSE)

  manifest <- read.csv(file.path(output_root, "scenario_input_manifest.csv"), stringsAsFactors = FALSE,
                       check.names = FALSE)
  specification <- read.csv(file.path(output_root, "scenario_specification.csv"), stringsAsFactors = FALSE,
                            check.names = FALSE)
  checks <- read.csv(file.path(output_root, "validation_checks.csv"), stringsAsFactors = FALSE,
                     check.names = FALSE)
  if (!nrow(manifest) || anyDuplicated(manifest$scenario) || anyDuplicated(specification$scenario) ||
      !setequal(manifest$scenario, specification$scenario) || !scenario_name %in% manifest$scenario) {
    stop("Scenario is not an approved, uniquely manifested dependent scenario: ", scenario_name, call. = FALSE)
  }
  row <- manifest[manifest$scenario == scenario_name, , drop = FALSE]
  spec_row <- specification[specification$scenario == scenario_name, , drop = FALSE]
  if (nrow(row) != 1L || nrow(spec_row) != 1L ||
      !identical(as.character(row$run_status), "planned") ||
      !identical(as.character(spec_row$run_status), "planned")) {
    stop("Only a planned dependent scenario may be run; aliases and infeasible records are not runnable.", call. = FALSE)
  }
  if (!isTRUE(row$validation_passed) || !all(checks$passed)) {
    stop("Dependent checkpoint validation does not clear this scenario.", call. = FALSE)
  }
  full_schedule <- read.csv(file.path(output_root, "fixed_study_seed_schedule.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
  planned_scenarios <- as.character(manifest$scenario[manifest$run_status == "planned"])
  schedule_counts <- table(as.character(full_schedule$scenario_name))
  if (!setequal(names(schedule_counts), planned_scenarios) ||
      any(schedule_counts[planned_scenarios] != 4L)) {
    stop("Committed dependent seed schedule does not exactly cover the planned scenarios.", call. = FALSE)
  }
  raw_path <- publication_input_path(output_root, row$raw_input_path)
  nma_path <- publication_input_path(output_root, row$nma_input_path)
  if (!file.exists(raw_path) || !file.exists(nma_path) ||
      !identical(sha(raw_path), as.character(row$raw_input_sha256)) ||
      !identical(sha(nma_path), as.character(row$nma_input_sha256))) {
    stop("Dependent raw or transformed input checksum mismatch for ", scenario_name, call. = FALSE)
  }
  stored_protected <- read.csv(file.path(output_root, "protected_artifact_checksums.csv"), stringsAsFactors = FALSE)
  protected_before <- core_collect_protected()
  if (nrow(stored_protected) != nrow(protected_before) || !all(stored_protected$matches_expected) ||
      !identical(as.character(stored_protected$artifact), as.character(protected_before$artifact)) ||
      !identical(as.character(stored_protected$expected_sha256), as.character(protected_before$expected_sha256))) {
    stop("Committed dependent protected-artifact ledger is invalid.", call. = FALSE)
  }
  core_assert_protected(protected_before, "dependent preflight")

  parent_root <- dependent_parent_root(project_root)
  parent <- core_verify_checkpoint(
    project_root, parent_root, "primary", input_manifest_sha256,
    parent_checkpoint_commit, check_git
  )
  if (!identical(as.character(parent$manifest$raw_input_sha256), as.character(metadata$source_parent_primary_raw_sha256)) ||
      !identical(sha(parent$raw_path), as.character(metadata$source_parent_primary_raw_sha256)) ||
      !identical(sha(file.path(parent_root, "scenario_input_manifest.csv")),
                 as.character(metadata$source_parent_manifest_sha256))) {
    stop("Dependent checkpoint parent-primary checksum or provenance does not match the reference core checkpoint.", call. = FALSE)
  }
  dependent_verify_checkpoint_commit(project_root, output_root, checkpoint_commit, check_git)
  contract <- dependent_treatment_contract(row)
  verified <- list(
    manifest = row, specification = spec_row, metadata = metadata, raw_path = raw_path,
    nma_path = nma_path, protected_before = protected_before, runtime_hashes = runtime_hashes,
    seeds = dependent_seed_schedule(output_root, scenario_name),
    treatment_contract = contract, parent = parent
  )
  verified$prepared <- dependent_assert_raw_fit_contract(verified)
  verified
}

dependent_write_result_artifacts <- function(result, runtime_warnings, attempt_path) {
  result$warnings <- unique(c(as.character(result$warnings), as.character(runtime_warnings)))
  core_write_result_artifacts(result, attempt_path)
  write.csv(data.frame(warning = unique(as.character(runtime_warnings)), stringsAsFactors = FALSE),
            file.path(attempt_path, "runtime_warnings.csv"), row.names = FALSE)
  invisible(TRUE)
}

dependent_metadata <- function(attempt, verified, scenario_name, analysis_commit,
                                            checkpoint_commit, parent_checkpoint_commit, docker_image,
                                            input_manifest_sha256, config, rerun_reason, status,
                                            started_utc, completed_utc = NA_character_,
                                            posterior_draws = NA_integer_, converged = NA,
                                            code_hashes_unchanged = NA,
                                            protected_artifacts_unchanged = NA,
                                            deviation = "none", result = NULL, test_only = FALSE) {
  data.frame(
    scenario = scenario_name, analysis_family = verified$manifest$analysis_family,
    definition = verified$manifest$definition, attempt = attempt$attempt,
    attempt_path = attempt$relative_path, status = status, started_utc = started_utc,
    completed_utc = completed_utc, analysis_commit = analysis_commit,
    checkpoint_commit = checkpoint_commit, parent_checkpoint_commit = parent_checkpoint_commit,
    checkpoint_input_builder_commit = verified$metadata$analysis_commit,
    parent_checkpoint_input_builder_commit = verified$parent$metadata$analysis_commit,
    docker_image = docker_image, input_manifest_sha256 = input_manifest_sha256,
    parent_primary_raw_sha256 = verified$metadata$source_parent_primary_raw_sha256,
    raw_input_sha256 = verified$manifest$raw_input_sha256,
    nma_input_sha256 = verified$manifest$nma_input_sha256,
    treatment_names = verified$manifest$treatment_names,
    n_treatments = verified$treatment_contract$expected_nodes,
    n_adapt = config$n.adapt, n_burnin = config$n.burnin, n_iter = config$n.iter,
    n_thin = config$n.thin, n_chains = config$n.chains,
    expected_posterior_draws = 100000L, posterior_draws = posterior_draws,
    chain_seed_1 = verified$seeds$seed[1], chain_seed_2 = verified$seeds$seed[2],
    chain_seed_3 = verified$seeds$seed[3], chain_seed_4 = verified$seeds$seed[4],
    rerun_reason = rerun_reason, test_only = isTRUE(test_only),
    production_mcmc_contract = identical(config, fixed_study_production_config()),
    converged = converged, code_hashes_unchanged = code_hashes_unchanged,
    protected_artifacts_unchanged = protected_artifacts_unchanged,
    max_psrf_upper_95 = if (is.null(result)) NA_real_ else result$diagnostics$max_psrf_upper_95,
    min_effective_size = if (is.null(result)) NA_real_ else result$diagnostics$min_effective_size,
    residual_deviance = if (is.null(result)) NA_real_ else result$residual_deviance_summary$mean,
    residual_deviance_per_arm = if (is.null(result)) NA_real_ else result$residual_deviance_summary$per_arm_mean,
    descriptive_dic = if (is.null(result)) NA_real_ else result$dic$DIC,
    no_fallback_or_extension = TRUE, deviation = deviation, stringsAsFactors = FALSE
  )
}

run_dependent_scenario <- function(
  project_root, scenario_name, input_manifest_sha256, analysis_commit, checkpoint_commit,
  parent_checkpoint_commit, docker_image,
  output_root = dependent_output_root(project_root), rerun_reason = "",
  preflight = FALSE, config = fixed_study_production_config(), test_only = FALSE,
  fit_function = fit_fixed_study_model, check_git = TRUE
) {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  output_root <- normalizePath(output_root, mustWork = TRUE)
  commits <- c(analysis_commit, checkpoint_commit, parent_checkpoint_commit)
  if (any(!grepl("^[0-9a-f]{40}$", commits))) {
    stop("ANALYSIS_COMMIT, CHECKPOINT_COMMIT, and PARENT_CHECKPOINT_COMMIT must be full commit hashes.", call. = FALSE)
  }
  config <- validate_fixed_study_config(config, test_only = test_only)
  if (!isTRUE(test_only) && !identical(config, fixed_study_production_config())) {
    stop("Only the locked four-chain production MCMC contract may be run.", call. = FALSE)
  }
  verified <- dependent_verify_checkpoint(
    project_root, output_root, scenario_name, input_manifest_sha256, checkpoint_commit,
    parent_checkpoint_commit, docker_image, check_git
  )
  if (isTRUE(preflight)) return(invisible(c(verified, list(config = config))))

  attempt <- core_next_attempt(publication_run_root(output_root), scenario_name, rerun_reason)
  dir.create(attempt$path, recursive = TRUE, showWarnings = FALSE)
  attempt$relative_path <- publication_relative_attempt(attempt$path, output_root)
  started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  meta <- dependent_metadata(
    attempt, verified, scenario_name, analysis_commit, checkpoint_commit, parent_checkpoint_commit,
    docker_image, input_manifest_sha256, config, rerun_reason, "started", started, test_only = test_only
  )
  publication_write_metadata(meta, file.path(attempt$path, "run_metadata.csv"))
  write.csv(verified$seeds, file.path(attempt$path, "seed_schedule_used.csv"), row.names = FALSE)
  write.csv(verified$runtime_hashes, file.path(attempt$path, "runtime_code_hashes_before.csv"), row.names = FALSE)
  write.csv(verified$protected_before, file.path(attempt$path, "protected_artifacts_before.csv"), row.names = FALSE)
  write.csv(verified$prepared$arms, file.path(attempt$path, "fitted_analysis_arms.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(attempt$path, "sessionInfo.txt"))
  cat(paste(started, scenario_name, sprintf("attempt_%03d", attempt$attempt), analysis_commit), "\n",
      file = file.path(dirname(publication_run_root(output_root)), "MCMC_STARTED"), sep = "", append = TRUE)

  runtime_warnings <- character()
  result <- tryCatch(withCallingHandlers(
    fit_function(
      scenario_name = scenario_name,
      raw_rows = read.csv(verified$raw_path, stringsAsFactors = FALSE, check.names = FALSE),
      treatment_map = verified$treatment_contract$map,
      treatment_names = verified$treatment_contract$names,
      model_file = file.path(project_root, "R", "model_fixed_study", "model.jags"),
      config = config, chain_seeds = verified$seeds, single_arm = "drop", master_seed = 12345L,
      scheduled_run_names = verified$seeds$scenario_name, preflight = FALSE,
      test_only = test_only, quiet = TRUE
    ), warning = function(w) {
      runtime_warnings <<- c(runtime_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
  completed <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  failed <- inherits(result, "error")
  failure <- if (failed) conditionMessage(result) else NA_character_
  if (!failed) {
    write_error <- tryCatch({
      dependent_write_result_artifacts(result, runtime_warnings, attempt$path)
      NULL
    }, error = identity)
    if (inherits(write_error, "error")) { failed <- TRUE; failure <- conditionMessage(write_error) }
  } else {
    write.csv(data.frame(warning = unique(runtime_warnings), stringsAsFactors = FALSE),
              file.path(attempt$path, "runtime_warnings.csv"), row.names = FALSE)
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
  hashes_after <- tryCatch(dependent_runtime_hashes(project_root), error = identity)
  hashes_ok <- FALSE
  if (inherits(hashes_after, "error")) {
    failed <- TRUE; failure <- conditionMessage(hashes_after)
  } else {
    hashes_ok <- identical(verified$runtime_hashes$item, hashes_after$item) &&
      identical(verified$runtime_hashes$sha256, hashes_after$sha256)
    write.csv(hashes_after, file.path(attempt$path, "runtime_code_hashes_after.csv"), row.names = FALSE)
    if (!hashes_ok) { failed <- TRUE; failure <- "Plan/model/runner code changed during MCMC." }
  }
  draws <- NA_integer_
  converged <- FALSE
  if (!failed) {
    draws <- nrow(as.matrix(result$samples))
    converged <- isTRUE(result$diagnostics$converged)
    if (!converged) failure <- "MCMC completed but the PSRF/ESS convergence contract failed."
    else if (!identical(config, fixed_study_production_config())) failure <- "A test-only MCMC configuration cannot be marked complete."
    else if (!identical(draws, 100000L)) failure <- "Retained posterior draw count did not equal the locked 100000-draw contract."
  }
  status <- if (failed) "failed" else if (is.na(failure)) "complete" else "nonconverged"
  meta <- dependent_metadata(
    attempt, verified, scenario_name, analysis_commit, checkpoint_commit, parent_checkpoint_commit,
    docker_image, input_manifest_sha256, config, rerun_reason, status, started, completed, draws,
    converged, hashes_ok, protected_ok, if (is.na(failure)) "none" else failure,
    if (failed) NULL else result, test_only
  )
  publication_write_metadata(meta, file.path(attempt$path, "run_metadata.csv"))
  marker <- if (identical(status, "complete")) "RUN_COMPLETE" else if (identical(status, "nonconverged")) "RUN_NONCONVERGED" else "RUN_FAILED"
  writeLines(if (is.na(failure)) completed else failure, file.path(attempt$path, marker))
  core_write_artifact_index(attempt$path)
  if (!identical(status, "complete")) {
    stop("Fixed-study dependent scenario did not meet the execution contract: ", scenario_name,
         ". ", failure, call. = FALSE)
  }
  invisible(list(attempt = attempt, metadata = meta, result = result))
}
