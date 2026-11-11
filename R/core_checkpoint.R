# =============================================================================
# Fixed-study immutable input checkpoint and pre-MCMC validation
# =============================================================================

core_protected_artifacts <- function() {
  manifest <- read.csv("data/scientific_checksums.csv", stringsAsFactors = FALSE)
  data.frame(artifact = manifest$path, path = manifest$path, kind = "file",
             expected_sha256 = manifest$sha256, stringsAsFactors = FALSE)
}

core_expected_delta_counts <- function() {
  c(
    primary = 0L,
    sensitivity_exclude_cevik = 2L,
    sensitivity_include_iachetta_kleibl = 6L,
    sensitivity_include_boisdron_overlap = 3L,
    sensitivity_exclude_onesti = 2L,
    sensitivity_exclude_genotype_informed_carriers = 8L,
    sensitivity_exclude_pk_informed_morel = 4L,
    sensitivity_exclude_all_treatment_moderated = 12L,
    sensitivity_wigle_early = 4L,
    sensitivity_cremolini_first_four = 2L,
    sensitivity_compound_lower_function = 3L,
    sensitivity_compound_reciprocal = 3L
  )
}

core_expected_single_arm_drops <- function(scenario_name) {
  if (scenario_name %in% c(
    "sensitivity_exclude_genotype_informed_carriers",
    "sensitivity_exclude_all_treatment_moderated"
  )) {
    return(c("Bozina_2021", "Wigle_2021_prospective"))
  }
  character()
}

core_flatten_designs <- function(designs) {
  data.frame(
    design_id = designs$design_id,
    design = designs$design,
    n_treatments = designs$n_treatments,
    baseline = designs$baseline,
    n_studies = designs$n_studies,
    treatments = vapply(designs$treatments, paste, character(1L), collapse = "|"),
    studies = vapply(designs$studies, paste, character(1L), collapse = "|"),
    stringsAsFactors = FALSE
  )
}

core_flatten_edges <- function(edges) {
  data.frame(
    treatment1 = edges$treatment1,
    treatment2 = edges$treatment2,
    n_studies = edges$n_studies,
    studies = vapply(edges$studies, paste, character(1L), collapse = "|"),
    stringsAsFactors = FALSE
  )
}

core_validation_report <- function(output_root, manifest, checks, protected, geometry) {
  status <- if (all(checks$passed) && all(protected$matches_expected)) "PASS" else "FAIL"
  scenario_lines <- vapply(seq_len(nrow(manifest)), function(index) {
    paste0(
      "| `", manifest$scenario[index], "` | ", manifest$studies[index], " | ",
      manifest$arms[index], " | ", manifest$total_n[index], " | ",
      manifest$total_r[index], " | `", manifest$raw_input_sha256[index], "` |"
    )
  }, character(1L))
  protected_lines <- vapply(seq_len(nrow(protected)), function(index) {
    paste0(
      "| `", protected$artifact[index], "` | `", protected$sha256_before[index],
      "` | `", protected$sha256_after[index], "` | ", protected$matches_expected[index], " |"
    )
  }, character(1L))
  lines <- c(
    "# Fixed-study pre-MCMC validation report", "",
    paste0("Status: **", status, "**"), "",
    paste0("Scenarios: ", nrow(manifest), "; checks: ", nrow(checks),
           "; passed: ", sum(checks$passed), "."), "",
    "## Primary network geometry", "",
    paste0("- Treatments: `", paste(geometry$treatments, collapse = ","), "`."),
    paste0("- Connected: `", geometry$connectedness$connected, "`."),
    paste0("- Unique designs: `", nrow(geometry$designs), "`."),
    paste0("- Direct comparisons: `", nrow(geometry$direct_edges), "`."),
    paste0("- Enumerated elementary closed loops: `", nrow(geometry$loop_table), "`."),
    paste0("- Design contrasts `m`: `", geometry$m, "`; rank(A): `", geometry$rank_A,
           "`; inconsistency df `q`: `", geometry$q, "`."),
    paste0("- Independent arm-level rank increase: `",
           geometry$arm_level_rank_check$inconsistency_rank_increase,
           "` (expected `", geometry$q, "`)."), "",
    "## Scenario inputs", "",
    "| Scenario | Source strata | Source arms | Source N | Source events | Raw SHA-256 |",
    "|---|---:|---:|---:|---:|---|", scenario_lines, "",
    "## Protected artifacts", "",
    "| Artifact | Before | After | Match |", "|---|---|---|---|", protected_lines, "",
    "## Validation coverage", "",
    "The checkpoint verifies all twelve approved source transformations, exact scenario totals and row deltas, source-specific primary rows, treatment mappings, unique study/treatment arms, integer binomial validity, fixed-study identifiability, presence of all five nodes, network connectivity, designs, direct comparisons, closed loops, two independent inconsistency-df calculations, immutable input checksums, deterministic seed scheduling, and protected-artifact stability.",
    "",
    "The genotype-informed and combined treatment-moderation source scenarios retain the prospective Wigle and Bozina standard-care comparators for provenance. Once their carrier arms are removed, those comparator rows are single-arm strata and therefore cannot identify a relative effect under fixed study effects. The transformed model inputs mechanically omit those two orphan strata; the manifest records both source-level and fitted-input totals.",
    "",
    "MCMC is not authorized unless this report and every machine-readable check pass."
  )
  writeLines(lines, file.path(output_root, "VALIDATION_REPORT.md"))
}

generate_core_checkpoint <- function(
  raw_data,
  output_root,
  metadata,
  protected_artifacts = core_protected_artifacts()
) {
  required_metadata <- c("analysis_commit", "docker_image", "input_manifest_sha256")
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata)) {
    stop("Missing checkpoint metadata: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  input_dir <- file.path(output_root, "inputs")
  geometry_dir <- output_root
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(geometry_dir, recursive = TRUE, showWarnings = FALSE)

  protected_before <- checksum_protected_artifacts(protected_artifacts)
  if (!all(protected_before$matches_expected)) {
    stop("Protected artifact checksum mismatch before checkpoint generation.", call. = FALSE)
  }

  spec <- core_scenario_spec()
  primary <- build_core_scenario(raw_data, "primary")
  expected_delta_counts <- core_expected_delta_counts()
  check_rows <- list()
  manifest_rows <- list()
  delta_rows <- list()
  geometry_rows <- list()

  add_check <- function(scenario, check, passed, detail) {
    check_rows[[length(check_rows) + 1L]] <<- validation_row(scenario, check, passed, detail)
  }
  add_check("all", "scenario_count", nrow(spec) == 12L, paste(nrow(spec), "approved scenarios"))

  for (index in seq_len(nrow(spec))) {
    scenario_name <- spec$scenario[index]
    scenario <- build_core_scenario(raw_data, scenario_name)
    validate_core_scenario(scenario, scenario_name)
    prepared <- prepare_fixed_study_data(
      scenario,
      core_treatment_map(),
      core_treatment_names(),
      single_arm = "drop"
    )
    expected_single_arm_drops <- core_expected_single_arm_drops(scenario_name)
    nma <- transform_to_nma_format(
      prepared$arms[, c("Study", "T", "N", "R")],
      core_treatment_map(),
      aggregate_wt = TRUE
    )
    validate_nma_data(nma)
    geometry <- analyze_network_geometry(nma[, c("Study", "t.id")])

    raw_relative <- file.path("inputs", paste0(scenario_name, "_raw.csv"))
    nma_relative <- file.path("inputs", paste0(scenario_name, "_nma.csv"))
    raw_path <- file.path(output_root, raw_relative)
    nma_path <- file.path(output_root, nma_relative)
    write.csv(scenario, raw_path, row.names = FALSE)
    write.csv(nma, nma_path, row.names = FALSE)

    totals <- core_totals(scenario)
    nodes <- core_node_totals(scenario)
    manifest_rows[[index]] <- data.frame(
      scenario = scenario_name,
      definition = spec$definition[index],
      studies = unname(totals["studies"]),
      arms = unname(totals["arms"]),
      total_n = unname(totals["N"]),
      total_r = unname(totals["R"]),
      hapb3_n = unname(nodes$HapB3["N"]),
      hapb3_r = unname(nodes$HapB3["R"]),
      fixed_model_studies = prepared$jags_data$Nstudy,
      fixed_model_arms = prepared$jags_data$Narm,
      fixed_model_n = sum(prepared$arms$N),
      fixed_model_r = sum(prepared$arms$R),
      dropped_single_arm_studies = paste(prepared$dropped_single_arm_studies, collapse = "|"),
      n_designs = nrow(geometry$designs),
      n_direct_comparisons = nrow(geometry$direct_edges),
      n_closed_loops = nrow(geometry$loop_table),
      inconsistency_df = geometry$q,
      raw_input_path = raw_relative,
      raw_input_sha256 = sha256_file(raw_path),
      nma_input_path = nma_relative,
      nma_input_sha256 = sha256_file(nma_path),
      validation_passed = TRUE,
      stringsAsFactors = FALSE
    )
    delta <- scenario_row_delta(primary, scenario, scenario_name)
    delta_rows[[index]] <- delta
    geometry_rows[[index]] <- data.frame(
      scenario = scenario_name,
      connected = geometry$connectedness$connected,
      treatments = paste(geometry$treatments, collapse = "|"),
      designs = nrow(geometry$designs),
      direct_comparisons = nrow(geometry$direct_edges),
      closed_loops = nrow(geometry$loop_table),
      m = geometry$m,
      rank_A = geometry$rank_A,
      q = geometry$q,
      q_from_rank = geometry$q_from_rank,
      q_from_arm_rank = geometry$arm_level_rank_check$inconsistency_rank_increase,
      stringsAsFactors = FALSE
    )

    add_check(scenario_name, "unique_study_treatment", !anyDuplicated(scenario[c("Study", "T")]), "No duplicate source arms")
    add_check(scenario_name, "binomial_valid", all(scenario$N > 0 & scenario$R >= 0 & scenario$R <= scenario$N & scenario$N %% 1 == 0 & scenario$R %% 1 == 0), "Integer N > 0 and 0 <= R <= N")
    add_check(
      scenario_name,
      "fixed_study_identifiable",
      identical(prepared$dropped_single_arm_studies, expected_single_arm_drops) && prepared$network_connected,
      if (length(expected_single_arm_drops)) {
        paste("Mechanically dropped orphan comparator strata:", paste(expected_single_arm_drops, collapse = ", "))
      } else {
        "No single-arm strata; the five-node network is connected"
      }
    )
    add_check(scenario_name, "all_nodes_present", setequal(sort(unique(prepared$arms$treatment)), 1:5), "All five modeled treatment nodes are observed")
    add_check(scenario_name, "geometry_rank", geometry$rank_A == 4L && geometry$q == geometry$q_from_rank && geometry$arm_level_rank_check$passes, paste("m/rank/q", geometry$m, geometry$rank_A, geometry$q, sep = "/"))
    add_check(scenario_name, "row_delta_contract", nrow(delta) == expected_delta_counts[[scenario_name]], paste(nrow(delta), "changed rows"))
    add_check(scenario_name, "input_checksums", nchar(sha256_file(raw_path)) == 64L && nchar(sha256_file(nma_path)) == 64L, "Raw and transformed inputs checksummed")
  }

  manifest <- do.call(rbind, manifest_rows)
  deltas <- do.call(rbind, delta_rows)
  if (!nrow(deltas)) deltas <- empty_row_delta()
  geometry_summary <- do.call(rbind, geometry_rows)
  primary_nma <- read.csv(file.path(output_root, manifest$nma_input_path[manifest$scenario == "primary"]), stringsAsFactors = FALSE)
  primary_geometry <- analyze_network_geometry(primary_nma[, c("Study", "t.id")])

  seed_run_ids <- c(spec$scenario, "inconsistency_sd2", "inconsistency_sd1")
  seed_schedule <- fixed_study_seed_schedule(seed_run_ids, master_seed = 12345L)
  add_check("all", "seed_schedule", nrow(seed_schedule) == length(seed_run_ids) * 4L && !anyDuplicated(seed_schedule$seed), "Four unique deterministic chain seeds for every core and inconsistency run")

  protected_after <- checksum_protected_artifacts(protected_artifacts)
  protected <- data.frame(
    artifact = protected_before$artifact,
    path = protected_before$path,
    kind = protected_before$kind,
    expected_sha256 = protected_before$expected_sha256,
    sha256_before = protected_before$observed_sha256,
    sha256_after = protected_after$observed_sha256,
    matches_expected = protected_before$matches_expected & protected_after$matches_expected &
      protected_before$observed_sha256 == protected_after$observed_sha256,
    stringsAsFactors = FALSE
  )
  add_check("all", "protected_artifacts", all(protected$matches_expected), "All six protected scientific artifacts remained byte-identical")
  checks <- do.call(rbind, check_rows)
  if (!all(checks$passed) || !all(protected$matches_expected)) {
    stop("Fixed-study pre-MCMC validation failed; inputs are not cleared.", call. = FALSE)
  }

  metadata_row <- data.frame(
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    analysis_commit = metadata$analysis_commit,
    docker_image = metadata$docker_image,
    input_manifest_sha256 = metadata$input_manifest_sha256,
    scenario_builder_sha256 = sha256_file("R/core_scenarios.R"),
    checkpoint_code_sha256 = sha256_file("R/core_checkpoint.R"),
    fixed_model_code_sha256 = sha256_file("R/fixed_study_model.R"),
    fixed_model_jags_sha256 = sha256_file("R/model_fixed_study/model.jags"),
    geometry_code_sha256 = sha256_file("R/network_geometry.R"),
    stringsAsFactors = FALSE
  )

  write.csv(spec, file.path(output_root, "scenario_specification.csv"), row.names = FALSE)
  write.csv(manifest, file.path(output_root, "scenario_input_manifest.csv"), row.names = FALSE)
  write.csv(deltas, file.path(output_root, "scenario_row_deltas.csv"), row.names = FALSE, na = "")
  write.csv(checks, file.path(output_root, "validation_checks.csv"), row.names = FALSE)
  write.csv(protected, file.path(output_root, "protected_artifact_checksums.csv"), row.names = FALSE)
  write.csv(metadata_row, file.path(output_root, "checkpoint_metadata.csv"), row.names = FALSE)
  write.csv(seed_schedule, file.path(output_root, "core_seed_schedule.csv"), row.names = FALSE)
  write.csv(geometry_summary, file.path(geometry_dir, "scenario_geometry.csv"), row.names = FALSE)
  write.csv(core_flatten_designs(primary_geometry$designs), file.path(geometry_dir, "primary_designs.csv"), row.names = FALSE)
  write.csv(core_flatten_edges(primary_geometry$direct_edges), file.path(geometry_dir, "primary_direct_comparisons.csv"), row.names = FALSE)
  write.csv(primary_geometry$loop_table, file.path(geometry_dir, "primary_closed_loops.csv"), row.names = FALSE)
  write.csv(primary_geometry$A_mapping, file.path(geometry_dir, "primary_A_mapping.csv"), row.names = FALSE)
  write.csv(cbind(contrast = rownames(primary_geometry$A), as.data.frame(primary_geometry$A)), file.path(geometry_dir, "primary_A_matrix.csv"), row.names = FALSE)
  core_validation_report(output_root, manifest, checks, protected, primary_geometry)

  list(
    manifest = read.csv(file.path(output_root, "scenario_input_manifest.csv"), stringsAsFactors = FALSE),
    checks = read.csv(file.path(output_root, "validation_checks.csv"), stringsAsFactors = FALSE),
    protected = protected,
    deltas = deltas,
    geometry = primary_geometry,
    seed_schedule = seed_schedule
  )
}
