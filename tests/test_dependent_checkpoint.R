source("R/common/data_utils.R")
source("R/arm_rows.R")
source("R/scenario_helpers.R")
source("R/checksums.R")
source("R/core_scenarios.R")
source("R/core_checkpoint.R")
source("R/fixed_study_model.R")
source("R/dependent_scenarios.R")
source("R/dependent_checkpoint.R")

output_root <- tempfile("fixed-study-dependent-")
dir.create(output_root, recursive = TRUE)
on.exit(unlink(output_root, recursive = TRUE), add = TRUE)
metadata <- list(
  analysis_commit = "test-commit",
  docker_image = "test-image",
  input_manifest_sha256 = "0000000000000000000000000000000000000000000000000000000000000000"
)

unresolved <- tryCatch({
  generate_dependent_checkpoint(getwd(), output_root, metadata)
  NULL
}, error = identity)
stopifnot(inherits(unresolved, "error"), grepl("explicit LD second-target policy", conditionMessage(unresolved)))

result <- generate_dependent_checkpoint(
  project_root = getwd(), output_root = output_root, metadata = metadata,
  ld_second_target_policy = "wigle_prospective"
)
required <- c(
  "scenario_specification.csv", "scenario_input_manifest.csv", "scenario_row_deltas.csv",
  "leave_one_source_out_map.csv", "ld_target_policy.csv", "fixed_study_seed_schedule.csv",
  "validation_checks.csv", "protected_artifact_checksums.csv", "checkpoint_metadata.csv",
  "VALIDATION_REPORT.md"
)
stopifnot(all(file.exists(file.path(output_root, required))))

manifest <- read.csv(file.path(output_root, "scenario_input_manifest.csv"), stringsAsFactors = FALSE)
checks <- read.csv(file.path(output_root, "validation_checks.csv"), stringsAsFactors = FALSE)
protected <- read.csv(file.path(output_root, "protected_artifact_checksums.csv"), stringsAsFactors = FALSE)
seeds <- read.csv(file.path(output_root, "fixed_study_seed_schedule.csv"), stringsAsFactors = FALSE)
metadata_row <- read.csv(file.path(output_root, "checkpoint_metadata.csv"), stringsAsFactors = FALSE)
deltas <- read.csv(file.path(output_root, "scenario_row_deltas.csv"), stringsAsFactors = FALSE)

stopifnot(
  nrow(manifest) == 31L,
  sum(manifest$run_status == "planned") == 27L,
  sum(manifest$run_status == "alias") == 3L,
  sum(manifest$run_status == "infeasible") == 1L,
  all(manifest$validation_passed),
  all(file.exists(file.path(output_root, manifest$raw_input_path[manifest$run_status == "planned"]))),
  all(nchar(manifest$raw_input_sha256[manifest$run_status == "planned"]) == 64L),
  all(nchar(manifest$nma_input_sha256[manifest$run_status == "planned"]) == 64L),
  nrow(seeds) == 108L,
  !anyDuplicated(seeds$seed),
  identical(sort(unique(seeds$chain)), 1:4),
  !any(seeds$scenario_name %in% manifest$scenario[manifest$run_status != "planned"]),
  all(checks$passed),
  all(protected$matches_expected),
  nchar(metadata_row$source_parent_primary_raw_sha256) == 64L
)

strict_delta <- deltas[deltas$scenario == "genotyping_wt_split_strict", , drop = FALSE]
stopifnot(
  nrow(strict_delta) == 4L,
  setequal(strict_delta$Study, c("Froehlich_2015", "Jennings_2013")),
  all(strict_delta$T %in% c("WT_Clean", "WT_Biased"))
)
aliases <- manifest[manifest$run_status == "alias", , drop = FALSE]
stopifnot(
  identical(sort(aliases$alias_of), sort(c("sensitivity_exclude_cevik", "sensitivity_exclude_onesti", "sensitivity_exclude_pk_informed_morel"))),
  all(!nzchar(aliases$raw_input_path)),
  all(nchar(aliases$alias_parent_raw_input_sha256) == 64L)
)
infeasible <- manifest[manifest$run_status == "infeasible", , drop = FALSE]
stopifnot(
  nrow(infeasible) == 1L,
  identical(infeasible$scenario, "sparse_exclude_n5_missing_13"),
  !nzchar(infeasible$raw_input_path)
)
stopifnot(all(result$checks$passed), !file.exists(file.path(output_root, "MCMC_STARTED")))

cat("All fixed-study dependent checkpoint tests passed.\n")
