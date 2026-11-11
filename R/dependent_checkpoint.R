# =============================================================================
# Fixed-study dependent input checkpoint (no MCMC)
# =============================================================================

dependent_protected_artifacts <- function() {
  core_protected_artifacts()
}

dependent_expected_delta <- function(primary, scenario, spec_row, loo_map) {
  scenario_name <- as.character(spec_row$scenario)
  delta <- scenario_row_delta(primary, scenario, scenario_name)
  expected <- switch(
    scenario_name,
    genotyping_wt_split_loose = empty_row_delta(),
    genotyping_wt_split_strict = {
      wanted <- primary[
        primary$Study %in% c("Froehlich_2015", "Jennings_2013") & primary$T == "WT_Clean",
        , drop = FALSE
      ]
      if (nrow(wanted) != 2L) stop("Strict WT delta source rows are unavailable.", call. = FALSE)
      strict <- primary
      strict$T[strict$Study %in% wanted$Study & strict$T == "WT_Clean"] <- "WT_Biased"
      scenario_row_delta(primary, strict, scenario_name)
    },
    sparse_exclude_n2 = scenario_row_delta(
      primary, dependent_sparse_raw(primary, 2L), scenario_name
    ),
    {
      if (grepl("^ld_(best|worst)_", scenario_name)) {
        expected_delta <- delta
        if (nrow(expected_delta) != 1L || expected_delta$change_type != "modified" ||
            expected_delta$T != "HapB3_1129or1236" || expected_delta$scenario_N != expected_delta$primary_N - 1L ||
            !(expected_delta$scenario_R %in% c(expected_delta$primary_R, expected_delta$primary_R - 1L))) {
          stop("LD scenario must change exactly one HapB3 arm by one carrier.", call. = FALSE)
        }
        expected_delta
      } else if (grepl("^loo_exclude_", scenario_name)) {
        group <- sub("^loo_exclude_", "", scenario_name)
        members <- loo_map$Study[loo_map$loo_group == group]
        expected <- primary[!primary$Study %in% members, , drop = FALSE]
        scenario_row_delta(primary, expected, scenario_name)
      } else {
        stop("No exact row-delta contract is registered for ", scenario_name, call. = FALSE)
      }
    }
  )
  if (!identical(delta, expected)) stop("Row delta contract failed for ", scenario_name, call. = FALSE)
  delta
}

dependent_alias_row <- function(spec_row, parent_metadata) {
  alias_manifest <- parent_metadata$parent_manifest
  core <- alias_manifest[alias_manifest$scenario == spec_row$alias_of, , drop = FALSE]
  if (nrow(core) != 1L || !isTRUE(core$validation_passed) ||
      !file.exists(publication_input_path(parent_metadata$parent_root, core$raw_input_path)) ||
      !identical(sha256_file(publication_input_path(parent_metadata$parent_root, core$raw_input_path)), core$raw_input_sha256)) {
    stop("Completed-core alias does not resolve to a validated parent input: ", spec_row$alias_of, call. = FALSE)
  }
  data.frame(
    scenario = spec_row$scenario,
    analysis_family = spec_row$analysis_family,
    run_status = "alias",
    definition = spec_row$definition,
    target_study = spec_row$target_study,
    alias_of = spec_row$alias_of,
    infeasible_reason = spec_row$infeasible_reason,
    treatment_names = spec_row$treatment_names,
    studies = core$studies, arms = core$arms, total_n = core$total_n, total_r = core$total_r,
    fixed_model_studies = core$fixed_model_studies,
    fixed_model_arms = core$fixed_model_arms,
    fixed_model_n = core$fixed_model_n,
    fixed_model_r = core$fixed_model_r,
    dropped_single_arm_studies = core$dropped_single_arm_studies,
    raw_input_path = NA_character_, raw_input_sha256 = NA_character_,
    nma_input_path = NA_character_, nma_input_sha256 = NA_character_,
    alias_parent_raw_input_path = core$raw_input_path,
    alias_parent_raw_input_sha256 = core$raw_input_sha256,
    validation_passed = TRUE,
    stringsAsFactors = FALSE
  )
}

dependent_infeasible_row <- function(spec_row) {
  data.frame(
    scenario = spec_row$scenario,
    analysis_family = spec_row$analysis_family,
    run_status = "infeasible",
    definition = spec_row$definition,
    target_study = spec_row$target_study,
    alias_of = spec_row$alias_of,
    infeasible_reason = spec_row$infeasible_reason,
    treatment_names = spec_row$treatment_names,
    studies = NA_integer_, arms = NA_integer_, total_n = NA_integer_, total_r = NA_integer_,
    fixed_model_studies = NA_integer_, fixed_model_arms = NA_integer_,
    fixed_model_n = NA_integer_, fixed_model_r = NA_integer_,
    dropped_single_arm_studies = NA_character_,
    raw_input_path = NA_character_, raw_input_sha256 = NA_character_,
    nma_input_path = NA_character_, nma_input_sha256 = NA_character_,
    alias_parent_raw_input_path = NA_character_, alias_parent_raw_input_sha256 = NA_character_,
    validation_passed = TRUE,
    stringsAsFactors = FALSE
  )
}

dependent_manifest_row <- function(spec_row, raw, nma, prepared, output_root) {
  scenario <- spec_row$scenario
  raw_relative <- file.path("inputs", paste0(scenario, "_raw.csv"))
  nma_relative <- file.path("inputs", paste0(scenario, "_nma.csv"))
  raw_path <- file.path(output_root, raw_relative)
  nma_path <- file.path(output_root, nma_relative)
  write.csv(raw, raw_path, row.names = FALSE)
  write.csv(nma, nma_path, row.names = FALSE)
  totals <- dependent_totals(raw)
  data.frame(
    scenario = scenario,
    analysis_family = spec_row$analysis_family,
    run_status = "planned",
    definition = spec_row$definition,
    target_study = spec_row$target_study,
    alias_of = NA_character_,
    infeasible_reason = NA_character_,
    treatment_names = spec_row$treatment_names,
    studies = unname(totals["studies"]), arms = unname(totals["arms"]),
    total_n = unname(totals["N"]), total_r = unname(totals["R"]),
    fixed_model_studies = prepared$jags_data$Nstudy,
    fixed_model_arms = prepared$jags_data$Narm,
    fixed_model_n = sum(prepared$arms$N), fixed_model_r = sum(prepared$arms$R),
    dropped_single_arm_studies = paste(prepared$dropped_single_arm_studies, collapse = "|"),
    raw_input_path = raw_relative, raw_input_sha256 = sha256_file(raw_path),
    nma_input_path = nma_relative, nma_input_sha256 = sha256_file(nma_path),
    alias_parent_raw_input_path = NA_character_, alias_parent_raw_input_sha256 = NA_character_,
    validation_passed = TRUE,
    stringsAsFactors = FALSE
  )
}

write_dependent_validation_report <- function(output_root, manifest, checks, protected, parent_metadata) {
  status <- if (all(checks$passed) && all(protected$matches_expected)) "PASS" else "FAIL"
  lines <- c(
    "# Fixed-study dependent pre-MCMC validation report", "",
    paste0("Status: **", status, "**"), "",
    "Dependent transforms use only the committed parent primary input; the workbook is not read.",
    paste0("Parent raw SHA-256: `", parent_metadata$primary_raw_sha256, "`"),
    paste0("Parent manifest SHA-256: `", parent_metadata$parent_manifest_sha256, "`"), "",
    paste0("Runnable new fits: ", sum(manifest$run_status == "planned"), "; ",
           "completed-core aliases: ", sum(manifest$run_status == "alias"), "; infeasible records: ",
           sum(manifest$run_status == "infeasible"), "."), "",
    "MCMC is not run or authorized by this generator. The sparse N <= 5 record has no input files because *13 is absent. Cevik, Morel, and Onesti leave-one-source-out records are aliases to their completed core scenarios, not new fits."
  )
  writeLines(lines, file.path(output_root, "VALIDATION_REPORT.md"))
}

generate_dependent_checkpoint <- function(
  project_root = getwd(),
  output_root = dependent_output_root(project_root),
  metadata,
  ld_second_target_policy = NULL,
  protected_artifacts = dependent_protected_artifacts()
) {
  required_metadata <- c("analysis_commit", "docker_image", "input_manifest_sha256")
  if (length(setdiff(required_metadata, names(metadata)))) {
    stop("Missing checkpoint metadata: ", paste(setdiff(required_metadata, names(metadata)), collapse = ", "), call. = FALSE)
  }
  if (file.exists(file.path(output_root, "MCMC_STARTED"))) {
    stop("MCMC has started; refusing to regenerate the fixed-study dependent checkpoint.", call. = FALSE)
  }
  parent_metadata <- dependent_parent_metadata(project_root)
  primary <- dependent_primary(project_root)
  specialized <- build_dependent_scenarios(primary, ld_second_target_policy)
  dir.create(file.path(output_root, "inputs"), recursive = TRUE, showWarnings = FALSE)

  protected_before <- checksum_protected_artifacts(protected_artifacts)
  if (!all(protected_before$matches_expected)) stop("Protected artifact checksum mismatch before generation.", call. = FALSE)
  checks <- list()
  add_check <- function(scenario, check, passed, detail) {
    checks[[length(checks) + 1L]] <<- validation_row(scenario, check, passed, detail)
  }
  add_check("all", "source_parent_checksum",
    identical(sha256_file(parent_metadata$primary_raw_path), parent_metadata$primary_raw_sha256),
    "Committed parent primary raw input matches its parent manifest")
  add_check("all", "loo_inventory",
    nrow(specialized$spec[specialized$spec$analysis_family == "loo", ]) == 23L &&
      sum(specialized$spec$analysis_family == "loo" & specialized$spec$run_status == "planned") == 20L,
    "23 grouped LOO records; Lee and Rosmarin grouped; three core aliases")

  manifest_rows <- list()
  delta_rows <- list()
  for (index in seq_len(nrow(specialized$spec))) {
    spec_row <- specialized$spec[index, , drop = FALSE]
    scenario <- spec_row$scenario
    if (spec_row$run_status == "infeasible") {
      manifest_rows[[index]] <- dependent_infeasible_row(spec_row)
      add_check(scenario, "infeasible_no_input", !scenario %in% names(specialized$raw), "*13 is missing; no runnable input is emitted")
      next
    }
    if (spec_row$run_status == "alias") {
      manifest_rows[[index]] <- dependent_alias_row(spec_row, parent_metadata)
      add_check(scenario, "completed_core_alias", !scenario %in% names(specialized$raw), paste("Alias of", spec_row$alias_of))
      next
    }
    raw <- specialized$raw[[scenario]]
    split_wt <- identical(spec_row$analysis_family, "genotyping")
    map <- if (split_wt) dependent_split_wt_map() else dependent_treatment_map()
    names <- if (split_wt) dependent_split_wt_names() else dependent_treatment_names()
    prepared <- prepare_fixed_study_data(raw, map, names, single_arm = "drop")
    # The emitted model input is exactly the identifiable fixed-study subset,
    # not merely a transformation of the source-provenance rows.
    nma <- transform_to_nma_format(
      prepared$arms[, c("Study", "T", "N", "R"), drop = FALSE],
      map, aggregate_wt = !split_wt
    )
    validate_nma_data(nma)
    nma <- nma[, c("Study", "s.id", "t.id", "r", "n")]
    delta <- dependent_expected_delta(primary, raw, spec_row, specialized$loo_map)
    delta_rows[[length(delta_rows) + 1L]] <- delta
    manifest_rows[[index]] <- dependent_manifest_row(spec_row, raw, nma, prepared, output_root)
    add_check(scenario, "exact_row_delta", TRUE, paste(nrow(delta), "locked row deltas"))
    add_check(scenario, "exact_totals", isTRUE(all.equal(
      as.numeric(dependent_totals(raw)),
      c(length(unique(raw$Study)), nrow(raw), sum(raw$N), sum(raw$R))
    )), "Manifest totals are recomputed from emitted raw input")
    add_check(scenario, "treatment_mapping", setequal(sort(unique(nma$t.id)), seq_along(names)), paste("Nodes", paste(seq_along(names), collapse = ",")))
    add_check(scenario, "fixed_study_identifiable", prepared$network_connected && prepared$jags_data$Nstudy > 0L,
      paste("Connected after dropping:", paste(prepared$dropped_single_arm_studies, collapse = "|")))
    add_check(scenario, "input_checksums", nchar(manifest_rows[[index]]$raw_input_sha256) == 64L && nchar(manifest_rows[[index]]$nma_input_sha256) == 64L,
      "Raw and NMA inputs are checksummed")
  }
  manifest <- do.call(rbind, manifest_rows)
  deltas <- if (length(delta_rows)) do.call(rbind, delta_rows) else empty_row_delta()
  runnable <- manifest$scenario[manifest$run_status == "planned"]
  seed_schedule <- fixed_study_seed_schedule(runnable, master_seed = 12345L)
  add_check("all", "seed_schedule", nrow(seed_schedule) == length(runnable) * 4L && !anyDuplicated(seed_schedule$seed),
    "Four deterministic fixed-study chain seeds for runnable new fits only")
  protected_after <- checksum_protected_artifacts(protected_artifacts)
  protected <- data.frame(
    artifact = protected_before$artifact, path = protected_before$path, kind = protected_before$kind,
    expected_sha256 = protected_before$expected_sha256,
    sha256_before = protected_before$observed_sha256, sha256_after = protected_after$observed_sha256,
    matches_expected = protected_before$matches_expected & protected_after$matches_expected &
      protected_before$observed_sha256 == protected_after$observed_sha256,
    stringsAsFactors = FALSE
  )
  add_check("all", "protected_artifacts", all(protected$matches_expected), "Protected artifacts remain byte-identical")
  checks <- do.call(rbind, checks)
  if (!all(checks$passed) || !all(protected$matches_expected)) stop("Fixed-study dependent checkpoint failed; inputs are not cleared.", call. = FALSE)

  metadata_row <- data.frame(
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    analysis_commit = metadata$analysis_commit, docker_image = metadata$docker_image,
    input_manifest_sha256 = metadata$input_manifest_sha256,
    ld_second_target_policy = ld_second_target_policy,
    source_parent_manifest_sha256 = parent_metadata$parent_manifest_sha256,
    source_parent_primary_raw_sha256 = parent_metadata$primary_raw_sha256,
    scenario_builder_sha256 = sha256_file("R/dependent_scenarios.R"),
    checkpoint_code_sha256 = sha256_file("R/dependent_checkpoint.R"),
    stringsAsFactors = FALSE
  )
  write.csv(specialized$spec, file.path(output_root, "scenario_specification.csv"), row.names = FALSE, na = "")
  write.csv(manifest, file.path(output_root, "scenario_input_manifest.csv"), row.names = FALSE, na = "")
  write.csv(deltas, file.path(output_root, "scenario_row_deltas.csv"), row.names = FALSE, na = "")
  write.csv(specialized$loo_map, file.path(output_root, "leave_one_source_out_map.csv"), row.names = FALSE, na = "")
  write.csv(specialized$ld_targets, file.path(output_root, "ld_target_policy.csv"), row.names = FALSE)
  write.csv(seed_schedule, file.path(output_root, "fixed_study_seed_schedule.csv"), row.names = FALSE)
  write.csv(checks, file.path(output_root, "validation_checks.csv"), row.names = FALSE)
  write.csv(protected, file.path(output_root, "protected_artifact_checksums.csv"), row.names = FALSE)
  write.csv(metadata_row, file.path(output_root, "checkpoint_metadata.csv"), row.names = FALSE)
  write_dependent_validation_report(output_root, manifest, checks, protected, parent_metadata)
  list(manifest = manifest, checks = checks, protected = protected, deltas = deltas, seed_schedule = seed_schedule)
}
