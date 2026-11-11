# Reference runner for the fixed-study-effect MCMC analyses.

core_plan_file <- function(project_root) {
  file.path(project_root, "data", "input_checksums.csv")
}

core_checkpoint_files <- function() c(
  "VALIDATION_REPORT.md", "checkpoint_metadata.csv", "core_seed_schedule.csv",
  "protected_artifact_checksums.csv", "scenario_input_manifest.csv",
  "scenario_specification.csv", "validation_checks.csv"
)

core_required <- function(name) {
  if (!exists(name, mode = "function", inherits = TRUE)) stop("Runner dependency is unavailable: ", name, call. = FALSE)
  get(name, mode = "function", inherits = TRUE)
}

core_read_one_row <- function(path, label) {
  if (!file.exists(path)) stop(label, " is unavailable: ", path, call. = FALSE)
  value <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(value) != 1L) stop(label, " must contain exactly one row.", call. = FALSE)
  value
}

core_runtime_hashes <- function(project_root) {
  sha <- core_required("sha256_file")
  paths <- c(
    public_runner = file.path(project_root, "R", "run_scenario.R"),
    public_utils = file.path(project_root, "R", "publication_utils.R"),
    data_utils = file.path(project_root, "R", "common", "data_utils.R"),
    input_manifest = core_plan_file(project_root),
    fixed_model = file.path(project_root, "R", "fixed_study_model.R"),
    fixed_model_jags = file.path(project_root, "R", "model_fixed_study", "model.jags"),
    scenario_builder = file.path(project_root, "R", "core_scenarios.R"),
    checkpoint_builder = file.path(project_root, "R", "core_checkpoint.R"),
    network_geometry = file.path(project_root, "R", "network_geometry.R"),
    core_runner = file.path(project_root, "R", "core_runner.R")
  )
  data.frame(item = names(paths), path = unname(paths), sha256 = vapply(paths, sha, character(1)), stringsAsFactors = FALSE)
}

core_verify_checkpoint_commit <- function(project_root, output_root, checkpoint_commit, check_git = TRUE) {
  publication_verify_inputs(project_root)
  invisible(TRUE)
}

core_collect_protected <- function() {
  checks <- core_required("checksum_protected_artifacts")(core_required("core_protected_artifacts")())
  if (nrow(checks) != 6L) stop("Expected exactly six protected artifacts.", call. = FALSE)
  checks
}

core_assert_protected <- function(checks, phase) {
  if (!all(checks$matches_expected)) {
    stop("Protected artifact checksum mismatch at ", phase, ": ", paste(checks$artifact[!checks$matches_expected], collapse = ", "), call. = FALSE)
  }
  invisible(checks)
}

core_verify_checkpoint <- function(project_root, output_root, scenario_name, input_manifest_sha256,
                                           checkpoint_commit = NA_character_, check_git = TRUE) {
  sha <- core_required("sha256_file")
  missing <- core_checkpoint_files()[!file.exists(file.path(output_root, core_checkpoint_files()))]
  if (length(missing)) stop("Immutable checkpoint artifact(s) missing: ", paste(missing, collapse = ", "), call. = FALSE)
  metadata <- core_read_one_row(file.path(output_root, "checkpoint_metadata.csv"), "Checkpoint metadata")
  required <- c("analysis_commit", "docker_image", "input_manifest_sha256", "scenario_builder_sha256", "checkpoint_code_sha256", "fixed_model_code_sha256", "fixed_model_jags_sha256", "geometry_code_sha256")
  if (!all(required %in% names(metadata))) stop("Checkpoint metadata is incomplete.", call. = FALSE)
  if (!identical(as.character(metadata$input_manifest_sha256), as.character(input_manifest_sha256))) stop("INPUT_MANIFEST_SHA256 differs from checkpoint.", call. = FALSE)
  hashes <- core_runtime_hashes(project_root)
  expected <- c(input_manifest = metadata$input_manifest_sha256, fixed_model = metadata$fixed_model_code_sha256,
                fixed_model_jags = metadata$fixed_model_jags_sha256, scenario_builder = metadata$scenario_builder_sha256,
                checkpoint_builder = metadata$checkpoint_code_sha256, network_geometry = metadata$geometry_code_sha256)
  observed <- setNames(hashes$sha256, hashes$item)
  bad <- names(expected)[observed[names(expected)] != unname(expected)]
  if (length(bad)) stop("Current input-manifest/model/checkpoint code hash mismatch: ", paste(bad, collapse = ", "), call. = FALSE)
  manifest <- read.csv(file.path(output_root, "scenario_input_manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  specification <- read.csv(file.path(output_root, "scenario_specification.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  checks <- read.csv(file.path(output_root, "validation_checks.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(manifest) || anyDuplicated(manifest$scenario) || !setequal(manifest$scenario, specification$scenario) || !scenario_name %in% manifest$scenario) {
    stop("Scenario is not an approved, uniquely manifested fixed-study scenario: ", scenario_name, call. = FALSE)
  }
  row <- manifest[manifest$scenario == scenario_name, , drop = FALSE]
  if (nrow(row) != 1L || !isTRUE(row$validation_passed) || !all(checks$passed)) stop("Checkpoint validation does not clear this scenario.", call. = FALSE)
  raw_path <- publication_input_path(output_root, row$raw_input_path)
  nma_path <- publication_input_path(output_root, row$nma_input_path)
  if (!file.exists(raw_path) || !file.exists(nma_path) || !identical(sha(raw_path), as.character(row$raw_input_sha256)) || !identical(sha(nma_path), as.character(row$nma_input_sha256))) {
    stop("Checkpoint input checksum mismatch for ", scenario_name, call. = FALSE)
  }
  stored <- read.csv(file.path(output_root, "protected_artifact_checksums.csv"), stringsAsFactors = FALSE)
  protected <- core_collect_protected()
  if (nrow(stored) != 6L || !all(stored$matches_expected) || !identical(as.character(stored$artifact), as.character(protected$artifact)) || !identical(as.character(stored$expected_sha256), as.character(protected$expected_sha256))) {
    stop("Committed protected-artifact ledger is invalid.", call. = FALSE)
  }
  core_assert_protected(protected, "preflight")
  schedule <- read.csv(file.path(output_root, "core_seed_schedule.csv"), stringsAsFactors = FALSE)
  schedule <- schedule[schedule$scenario_name == scenario_name, , drop = FALSE]
  schedule <- schedule[order(schedule$chain), , drop = FALSE]
  if (!identical(as.integer(schedule$chain), 1:4) || anyDuplicated(schedule$seed)) stop("Seed schedule does not provide four unique chains.", call. = FALSE)
  core_verify_checkpoint_commit(project_root, output_root, checkpoint_commit, check_git)
  list(manifest = row, raw_path = raw_path, nma_path = nma_path, metadata = metadata, seed_schedule = schedule, protected_before = protected, runtime_hashes = hashes)
}

core_assert_raw_fit_contract <- function(verified) {
  raw <- read.csv(verified$raw_path, stringsAsFactors = FALSE, check.names = FALSE)
  prepared <- core_required("prepare_fixed_study_data")(raw, core_required("core_treatment_map")(), core_required("core_treatment_names")(), single_arm = "drop")
  expected <- c(fixed_model_studies = prepared$jags_data$Nstudy, fixed_model_arms = prepared$jags_data$Narm,
                fixed_model_n = sum(prepared$arms$N), fixed_model_r = sum(prepared$arms$R),
                dropped_single_arm_studies = paste(prepared$dropped_single_arm_studies, collapse = "|"))
  actual <- c(fixed_model_studies = verified$manifest$fixed_model_studies, fixed_model_arms = verified$manifest$fixed_model_arms,
              fixed_model_n = verified$manifest$fixed_model_n, fixed_model_r = verified$manifest$fixed_model_r,
              dropped_single_arm_studies = ifelse(is.na(verified$manifest$dropped_single_arm_studies), "", verified$manifest$dropped_single_arm_studies))
  if (!identical(as.character(expected), as.character(actual))) stop("Raw input no longer produces the manifested fixed-study analysis input.", call. = FALSE)
  nma <- read.csv(verified$nma_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_nma_fields <- c("s.id", "t.id", "r", "n", "Study", "T")
  if (!identical(names(nma), required_nma_fields)) {
    stop("Committed transformed input has an unexpected fixed-study arm schema.", call. = FALSE)
  }
  prepared_as_nma <- data.frame(
    s.id = as.integer(prepared$arms$study),
    t.id = as.integer(prepared$arms$treatment),
    r = as.integer(prepared$arms$R),
    n = as.integer(prepared$arms$N),
    Study = as.character(prepared$arms$Study),
    T = as.character(prepared$arms$T),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  manifested_nma <- data.frame(
    s.id = as.integer(nma[["s.id"]]),
    t.id = as.integer(nma[["t.id"]]),
    r = as.integer(nma$r),
    n = as.integer(nma$n),
    Study = as.character(nma$Study),
    T = as.character(nma$T),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!identical(prepared_as_nma, manifested_nma)) {
    stop("Prepared fixed-study arms do not match the manifested transformed input row-for-row.", call. = FALSE)
  }
  prepared
}

core_next_attempt <- function(run_root, scenario_name, rerun_reason = "") {
  scenario_root <- file.path(run_root, scenario_name)
  existing <- if (dir.exists(scenario_root)) list.dirs(scenario_root, recursive = FALSE, full.names = FALSE) else character()
  numbers <- suppressWarnings(as.integer(sub("^attempt_", "", existing[grepl("^attempt_[0-9]{3}$", existing)])))
  next_number <- if (length(numbers)) max(numbers) + 1L else 1L
  if (next_number > 1L && !nzchar(trimws(rerun_reason))) stop("RERUN_REASON is required before creating a subsequent immutable attempt.", call. = FALSE)
  path <- file.path(scenario_root, sprintf("attempt_%03d", next_number))
  if (file.exists(path)) stop("Refusing to overwrite an existing attempt: ", path, call. = FALSE)
  list(attempt = next_number, path = path, rerun_reason = rerun_reason)
}

core_write_artifact_index <- function(attempt_path) {
  sha <- core_required("sha256_file")
  files <- list.files(attempt_path, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  files <- files[basename(files) != "artifact_checksums.csv"]
  base <- normalizePath(attempt_path, winslash = "/", mustWork = TRUE)
  table <- data.frame(artifact = substring(files, nchar(base) + 2L), size_bytes = file.info(files)$size,
                      sha256 = vapply(files, sha, character(1)), stringsAsFactors = FALSE)
  write.csv(table, file.path(attempt_path, "artifact_checksums.csv"), row.names = FALSE)
  table
}

core_write_result_artifacts <- function(result, attempt_path) {
  write.csv(result$treatment_effects, file.path(attempt_path, "posterior_treatment_effects.csv"), row.names = FALSE)
  write.csv(result$pairwise_effects, file.path(attempt_path, "posterior_pairwise_effects.csv"), row.names = FALSE)
  write.csv(result$heterogeneity, file.path(attempt_path, "posterior_heterogeneity.csv"), row.names = FALSE)
  write.csv(result$diagnostic_table, file.path(attempt_path, "diagnostics.csv"), row.names = FALSE)
  write.csv(result$residual_deviance_arms, file.path(attempt_path, "residual_deviance_arms.csv"), row.names = FALSE)
  write.csv(result$residual_deviance_summary, file.path(attempt_path, "residual_deviance_summary.csv"), row.names = FALSE)
  write.csv(as.data.frame(result$dic, stringsAsFactors = FALSE), file.path(attempt_path, "descriptive_dic.csv"), row.names = FALSE)
  write.csv(data.frame(warning = as.character(result$warnings), stringsAsFactors = FALSE), file.path(attempt_path, "warnings.csv"), row.names = FALSE)
  saveRDS(result$samples, file.path(attempt_path, "fixed_study_mcmc_samples.rds"))
  invisible(TRUE)
}

core_run_scenario <- function(project_root, scenario_name, input_manifest_sha256, analysis_commit, checkpoint_commit, docker_image,
                                      output_root = file.path(project_root, "results", "analysis", "core"), rerun_reason = "",
                                      preflight = FALSE, config = fixed_study_production_config(), test_only = FALSE,
                                      fit_function = fit_fixed_study_model, check_git = TRUE) {
  project_root <- normalizePath(project_root, mustWork = TRUE); output_root <- normalizePath(output_root, mustWork = TRUE)
  if (!grepl("^[0-9a-f]{40}$", analysis_commit) || !grepl("^[0-9a-f]{40}$", checkpoint_commit)) stop("ANALYSIS_COMMIT and CHECKPOINT_COMMIT must be full commit hashes.", call. = FALSE)
  verified <- core_verify_checkpoint(project_root, output_root, scenario_name, input_manifest_sha256, checkpoint_commit, check_git)
  prepared <- core_assert_raw_fit_contract(verified)
  config <- validate_fixed_study_config(config, test_only = test_only)
  if (!isTRUE(test_only) && !identical(config, fixed_study_production_config())) stop("Only the locked four-chain production MCMC contract may be run.", call. = FALSE)
  if (isTRUE(preflight)) return(invisible(list(verified = verified, prepared = prepared, config = config)))
  attempt <- core_next_attempt(publication_run_root(output_root), scenario_name, rerun_reason); dir.create(attempt$path, recursive = TRUE, showWarnings = FALSE)
  started <- format(Sys.time(), tz = "UTC", usetz = TRUE); relative <- publication_relative_attempt(attempt$path, output_root)
  meta <- data.frame(scenario = scenario_name, definition = verified$manifest$definition, attempt = attempt$attempt, attempt_path = relative,
    status = "started", started_utc = started, completed_utc = NA_character_, analysis_commit = analysis_commit, checkpoint_commit = checkpoint_commit,
    checkpoint_input_builder_commit = verified$metadata$analysis_commit, docker_image = docker_image, input_manifest_sha256 = input_manifest_sha256,
    raw_input_sha256 = verified$manifest$raw_input_sha256, nma_input_sha256 = verified$manifest$nma_input_sha256,
    n_adapt = config$n.adapt, n_burnin = config$n.burnin, n_iter = config$n.iter, n_thin = config$n.thin, n_chains = config$n.chains,
    expected_posterior_draws = as.integer(config$n.iter / config$n.thin * config$n.chains),
    chain_seed_1 = verified$seed_schedule$seed[1], chain_seed_2 = verified$seed_schedule$seed[2], chain_seed_3 = verified$seed_schedule$seed[3], chain_seed_4 = verified$seed_schedule$seed[4],
    single_arm_policy = "drop", dropped_single_arm_studies = paste(prepared$dropped_single_arm_studies, collapse = "|"), rerun_reason = rerun_reason,
    test_only = isTRUE(test_only), no_fallback_or_extension = TRUE, stringsAsFactors = FALSE)
  publication_write_metadata(meta, file.path(attempt$path, "run_metadata.csv"))
  write.csv(verified$seed_schedule, file.path(attempt$path, "seed_schedule_used.csv"), row.names = FALSE)
  write.csv(verified$runtime_hashes, file.path(attempt$path, "runtime_code_hashes_before.csv"), row.names = FALSE)
  write.csv(verified$protected_before, file.path(attempt$path, "protected_artifacts_before.csv"), row.names = FALSE)
  write.csv(prepared$arms, file.path(attempt$path, "fitted_analysis_arms.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(attempt$path, "sessionInfo.txt"))
  cat(paste(started, scenario_name, sprintf("attempt_%03d", attempt$attempt), analysis_commit), "\n",
      file = file.path(dirname(publication_run_root(output_root)), "MCMC_STARTED"), sep = "", append = TRUE)
  warnings <- character()
  result <- tryCatch(withCallingHandlers(fit_function(scenario_name = scenario_name, raw_rows = read.csv(verified$raw_path, stringsAsFactors = FALSE, check.names = FALSE),
    treatment_map = core_treatment_map(), treatment_names = core_treatment_names(), model_file = file.path(project_root, "R", "model_fixed_study", "model.jags"),
    config = config, chain_seeds = verified$seed_schedule, single_arm = "drop", master_seed = 12345L, scheduled_run_names = verified$seed_schedule$scenario_name,
    preflight = FALSE, test_only = test_only, quiet = TRUE), warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
  completed <- format(Sys.time(), tz = "UTC", usetz = TRUE); failed <- inherits(result, "error"); failure <- if (failed) conditionMessage(result) else NA_character_
  if (!failed) { result$warnings <- unique(c(result$warnings, warnings)); write_error <- tryCatch({ core_write_result_artifacts(result, attempt$path); NULL }, error = identity); if (inherits(write_error, "error")) { failed <- TRUE; failure <- conditionMessage(write_error) } }
  protected_after <- tryCatch(core_collect_protected(), error = identity)
  if (inherits(protected_after, "error")) { failed <- TRUE; failure <- conditionMessage(protected_after) } else { write.csv(protected_after, file.path(attempt$path, "protected_artifacts_after.csv"), row.names = FALSE); if (!all(protected_after$matches_expected)) { failed <- TRUE; failure <- "Protected artifact checksum mismatch after run." } }
  hashes_after <- tryCatch(core_runtime_hashes(project_root), error = identity)
  if (inherits(hashes_after, "error")) { failed <- TRUE; failure <- conditionMessage(hashes_after) } else { write.csv(hashes_after, file.path(attempt$path, "runtime_code_hashes_after.csv"), row.names = FALSE); if (!identical(verified$runtime_hashes$sha256, hashes_after$sha256)) { failed <- TRUE; failure <- "Plan/model/runner code changed during MCMC." } }
  draws <- NA_integer_; converged <- FALSE
  if (!failed) { draws <- nrow(as.matrix(result$samples)); converged <- isTRUE(result$diagnostics$converged) && identical(draws, as.integer(config$n.iter / config$n.thin * config$n.chains)); if (!converged) failure <- "MCMC completed but the PSRF/ESS or retained-draw convergence contract failed." }
  meta$status <- if (failed) "failed" else if (converged) "complete" else "nonconverged"; meta$completed_utc <- completed; meta$posterior_draws <- draws; meta$converged <- converged; meta$deviation <- if (is.na(failure)) "none" else failure
  if (!inherits(result, "error")) { meta$max_psrf_upper_95 <- result$diagnostics$max_psrf_upper_95; meta$min_effective_size <- result$diagnostics$min_effective_size; meta$residual_deviance <- result$residual_deviance_summary$mean; meta$residual_deviance_per_arm <- result$residual_deviance_summary$per_arm_mean; meta$descriptive_dic <- result$dic$DIC }
  publication_write_metadata(meta, file.path(attempt$path, "run_metadata.csv"))
  marker <- if (identical(meta$status, "complete")) "RUN_COMPLETE" else if (identical(meta$status, "nonconverged")) "RUN_NONCONVERGED" else "RUN_FAILED"; writeLines(if (is.na(failure)) completed else failure, file.path(attempt$path, marker))
  core_write_artifact_index(attempt$path)
  if (!identical(meta$status, "complete")) stop("Fixed-study scenario did not meet the execution contract: ", scenario_name, ". ", failure, call. = FALSE)
  invisible(list(attempt = attempt, metadata = meta, result = result))
}
