#!/usr/bin/env Rscript

# Internal checkpoint helpers; not a public rerun aggregation command.
# Verify and aggregate the bounded fixed-study dependent analyses. The
# inventory is policy-agnostic: the LD family may contain either four runs
# (one explicitly selected second target) or six runs (both tied targets).

dependent_results_read <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

dependent_results_required <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

dependent_results_one_row <- function(path, label) {
  data <- dependent_results_read(path, label)
  if (nrow(data) != 1L) stop(label, " must contain exactly one row.", call. = FALSE)
  data
}

dependent_results_bool <- function(value) {
  !is.na(value) & as.logical(value)
}

dependent_validate_inventory <- function(manifest) {
  dependent_results_required(
    manifest,
    c("scenario", "analysis_family", "run_status", "definition", "alias_of", "treatment_names"),
    "Dependent manifest"
  )
  if (anyDuplicated(manifest$scenario) || !all(manifest$run_status %in% c("planned", "alias", "infeasible"))) {
    stop("Dependent manifest contains duplicate scenarios or an unsupported run status.", call. = FALSE)
  }
  counts <- table(manifest$analysis_family, manifest$run_status)
  count <- function(family, status) {
    if (!family %in% rownames(counts) || !status %in% colnames(counts)) return(0L)
    unname(counts[family, status])
  }
  valid <- count("genotyping", "planned") == 2L &&
    count("sparse", "planned") == 1L && count("sparse", "infeasible") == 1L &&
    count("ld", "planned") %in% c(4L, 6L) &&
    count("loo", "planned") == 20L && count("loo", "alias") == 3L &&
    sum(manifest$run_status == "alias") == 3L && sum(manifest$run_status == "infeasible") == 1L
  if (!valid) stop("Dependent manifest does not match the bounded approved family inventory.", call. = FALSE)
  invisible(TRUE)
}

dependent_attempt_inventory <- function(output_root, manifest) {
  planned <- manifest[manifest$run_status == "planned", , drop = FALSE]
  metadata_paths <- sort(list.files(
    file.path(output_root, "runs"), pattern = "^run_metadata\\.csv$",
    recursive = TRUE, full.names = TRUE
  ))
  if (!length(metadata_paths)) stop("No dependent run metadata were found.", call. = FALSE)
  attempts <- do.call(rbind, lapply(metadata_paths, function(path) {
    row <- dependent_results_one_row(path, "Dependent run metadata")
    row$metadata_path <- substring(path, nchar(output_root) + 2L)
    row
  }))
  rownames(attempts) <- NULL
  dependent_results_required(
    attempts,
    c("scenario", "attempt", "attempt_path", "status", "analysis_family", "analysis_commit",
      "checkpoint_commit", "parent_checkpoint_commit", "docker_image", "input_manifest_sha256",
      "raw_input_sha256", "nma_input_sha256", "n_adapt", "n_burnin", "n_iter", "n_thin",
      "n_chains", "expected_posterior_draws", "posterior_draws", "test_only",
      "production_mcmc_contract", "converged", "code_hashes_unchanged",
      "protected_artifacts_unchanged", "max_psrf_upper_95", "min_effective_size",
      "no_fallback_or_extension", "deviation"),
    "Dependent run metadata inventory"
  )
  if (anyDuplicated(attempts[c("scenario", "attempt")]) ||
      !all(attempts$scenario %in% planned$scenario) ||
      !setequal(unique(attempts$scenario), planned$scenario)) {
    stop("Dependent attempt inventory does not exactly cover the planned manifest scenarios.", call. = FALSE)
  }
  selected_rows <- lapply(planned$scenario, function(scenario) {
    rows <- attempts[attempts$scenario == scenario & attempts$status == "complete", , drop = FALSE]
    if (nrow(rows) != 1L) {
      stop("Expected exactly one complete dependent attempt for ", scenario, "; found ", nrow(rows), ".", call. = FALSE)
    }
    rows
  })
  selected <- do.call(rbind, selected_rows)
  rownames(selected) <- NULL
  selected <- selected[match(planned$scenario, selected$scenario), , drop = FALSE]
  list(attempts = attempts, selected = selected, planned = planned)
}

dependent_verify_selected_contract <- function(selected, planned, checkpoint, seed_schedule) {
  mapped <- planned[match(selected$scenario, planned$scenario), , drop = FALSE]
  if (anyNA(mapped$scenario)) stop("Selected dependent attempts do not map to the manifest.", call. = FALSE)
  expected_image <- "runtime-recorded"
  checks <- c(
    family_mapping = identical(as.character(selected$analysis_family), as.character(mapped$analysis_family)),
    raw_input_hash = identical(as.character(selected$raw_input_sha256), as.character(mapped$raw_input_sha256)),
    nma_input_hash = identical(as.character(selected$nma_input_sha256), as.character(mapped$nma_input_sha256)),
    software_recorded = all(nzchar(selected$docker_image)),
    plan_hash = all(selected$input_manifest_sha256 == checkpoint$input_manifest_sha256),
    adaptation = all(selected$n_adapt == 50000L), burnin = all(selected$n_burnin == 125000L),
    iterations = all(selected$n_iter == 250000L), thinning = all(selected$n_thin == 10L),
    chains = all(selected$n_chains == 4L),
    expected_draws = all(selected$expected_posterior_draws == 100000L),
    retained_draws = all(selected$posterior_draws == 100000L),
    production_only = all(!dependent_results_bool(selected$test_only)),
    production_contract = all(dependent_results_bool(selected$production_mcmc_contract)),
    converged = all(dependent_results_bool(selected$converged)),
    code_unchanged = all(dependent_results_bool(selected$code_hashes_unchanged)),
    protected_unchanged = all(dependent_results_bool(selected$protected_artifacts_unchanged)),
    no_fallback = all(dependent_results_bool(selected$no_fallback_or_extension)),
    psrf = all(selected$max_psrf_upper_95 <= 1.05), ess = all(selected$min_effective_size >= 1000),
    no_deviation = all(selected$deviation == "none")
  )
  if (!all(checks)) {
    stop("Selected dependent attempts violate contract check(s): ",
         paste(names(checks)[!checks], collapse = ", "), ".", call. = FALSE)
  }
  expected_seed_counts <- table(seed_schedule$scenario_name)
  if (!setequal(names(expected_seed_counts), planned$scenario) || any(expected_seed_counts != 4L)) {
    stop("Dependent seed schedule does not exactly cover every planned scenario with four chains.", call. = FALSE)
  }
  for (index in seq_len(nrow(selected))) {
    scenario <- selected$scenario[[index]]
    expected <- seed_schedule[seed_schedule$scenario_name == scenario, , drop = FALSE]
    expected <- expected[order(expected$chain), , drop = FALSE]
    observed <- as.integer(selected[index, paste0("chain_seed_", 1:4), drop = TRUE])
    if (!identical(as.integer(expected$chain), 1:4) || !identical(observed, as.integer(expected$seed))) {
      stop("Selected attempt seed mismatch for ", scenario, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

dependent_verify_attempt_artifacts <- function(output_root, selected, checkpoint_protected) {
  rows <- list()
  counter <- 0L
  for (index in seq_len(nrow(selected))) {
    attempt_dir <- file.path(output_root, selected$attempt_path[[index]])
    if (!file.exists(file.path(attempt_dir, "RUN_COMPLETE"))) {
      stop("Missing RUN_COMPLETE marker for ", selected$scenario[[index]], ".", call. = FALSE)
    }
    runtime_before <- dependent_results_read(
      file.path(attempt_dir, "runtime_code_hashes_before.csv"), "runtime hashes before"
    )
    runtime_after <- dependent_results_read(
      file.path(attempt_dir, "runtime_code_hashes_after.csv"), "runtime hashes after"
    )
    if (!identical(runtime_before[c("item", "sha256")], runtime_after[c("item", "sha256")])) {
      stop("Runtime code hashes changed within ", selected$scenario[[index]], ".", call. = FALSE)
    }
    protected_after <- dependent_results_read(
      file.path(attempt_dir, "protected_artifacts_after.csv"), "protected artifacts after"
    )
    if (!all(dependent_results_bool(protected_after$matches_expected)) ||
        !identical(as.character(protected_after$artifact), as.character(checkpoint_protected$artifact)) ||
        !identical(as.character(protected_after$expected_sha256), as.character(checkpoint_protected$expected_sha256))) {
      stop("Protected-artifact ledger mismatch in ", selected$scenario[[index]], ".", call. = FALSE)
    }
    artifact_index <- dependent_results_read(
      file.path(attempt_dir, "artifact_checksums.csv"), "attempt artifact index"
    )
    dependent_results_required(artifact_index, c("artifact", "size_bytes", "sha256"), "Attempt artifact index")
    for (row_index in seq_len(nrow(artifact_index))) {
      artifact_path <- file.path(attempt_dir, artifact_index$artifact[[row_index]])
      observed_hash <- if (file.exists(artifact_path)) sha256_file(artifact_path) else NA_character_
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        scenario = selected$scenario[[index]], attempt = selected$attempt[[index]],
        artifact = artifact_index$artifact[[row_index]],
        relative_path = substring(artifact_path, nchar(output_root) + 2L),
        size_bytes = artifact_index$size_bytes[[row_index]],
        recorded_sha256 = artifact_index$sha256[[row_index]], observed_sha256 = observed_hash,
        exists = file.exists(artifact_path),
        size_matches = file.exists(artifact_path) && file.info(artifact_path)$size == artifact_index$size_bytes[[row_index]],
        checksum_matches = !is.na(observed_hash) && identical(observed_hash, artifact_index$sha256[[row_index]]),
        stringsAsFactors = FALSE
      )
    }
  }
  artifacts <- do.call(rbind, rows)
  if (!all(artifacts$exists & artifacts$size_matches & artifacts$checksum_matches)) {
    stop("One or more dependent attempt artifacts failed checksum verification.", call. = FALSE)
  }
  artifacts
}

dependent_key_rows <- function(pairwise, family, scenario) {
  rows <- pairwise[pairwise$scenario == scenario & pairwise$numerator == "HapB3", , drop = FALSE]
  wanted <- if (identical(family, "genotyping")) c("WT_Clean", "WT_Biased") else "WT"
  rows <- rows[rows$denominator %in% wanted, , drop = FALSE]
  if (nrow(rows) != length(wanted) || !setequal(rows$denominator, wanted)) {
    stop("Missing key HapB3 comparison(s) for ", scenario, ".", call. = FALSE)
  }
  rows$comparison <- paste(rows$numerator, "vs", rows$denominator)
  rows
}

