source("R/common/data_utils.R")
source("R/arm_rows.R")
source("R/scenario_helpers.R")
source("R/checksums.R")
source("R/core_scenarios.R")
source("R/fixed_study_model.R")
source("R/network_geometry.R")
source("R/core_checkpoint.R")

output_root <- tempfile("fixed-study-checkpoint-")
dir.create(output_root, recursive = TRUE)
on.exit(unlink(output_root, recursive = TRUE), add = TRUE)
metadata <- list(
  analysis_commit = "test-commit",
  docker_image = "test-image",
  input_manifest_sha256 = "0000000000000000000000000000000000000000000000000000000000000000"
)
result <- generate_core_checkpoint(read.csv("data/core/primary_raw.csv", stringsAsFactors = FALSE), output_root, metadata)

required <- c(
  "scenario_specification.csv", "scenario_input_manifest.csv", "scenario_row_deltas.csv",
  "validation_checks.csv",
  "protected_artifact_checksums.csv", "checkpoint_metadata.csv", "core_seed_schedule.csv",
  "VALIDATION_REPORT.md", "scenario_geometry.csv",
  "primary_designs.csv", "primary_direct_comparisons.csv",
  "primary_closed_loops.csv", "primary_A_mapping.csv",
  "primary_A_matrix.csv"
)
stopifnot(all(file.exists(file.path(output_root, required))))

manifest <- read.csv(file.path(output_root, "scenario_input_manifest.csv"), stringsAsFactors = FALSE)
checks <- read.csv(file.path(output_root, "validation_checks.csv"), stringsAsFactors = FALSE)
protected <- read.csv(file.path(output_root, "protected_artifact_checksums.csv"), stringsAsFactors = FALSE)
seeds <- read.csv(file.path(output_root, "core_seed_schedule.csv"), stringsAsFactors = FALSE)
geometry <- read.csv(file.path(output_root, "scenario_geometry.csv"), stringsAsFactors = FALSE)

stopifnot(
  identical(manifest$scenario, core_scenario_spec()$scenario),
  nrow(manifest) == 12L,
  all(manifest$validation_passed),
  all(nchar(manifest$raw_input_sha256) == 64L),
  all(nchar(manifest$nma_input_sha256) == 64L),
  all(file.exists(file.path(output_root, manifest$raw_input_path))),
  all(file.exists(file.path(output_root, manifest$nma_input_path)))
)
primary <- manifest[manifest$scenario == "primary", , drop = FALSE]
stopifnot(
  primary$studies == 25L, primary$arms == 79L, primary$total_n == 12486L,
  primary$total_r == 3821L, primary$hapb3_n == 311L, primary$hapb3_r == 133L,
  primary$fixed_model_studies == 25L, primary$fixed_model_arms == 79L,
  primary$fixed_model_n == 12486L, primary$fixed_model_r == 3821L,
  identical(primary$dropped_single_arm_studies, "")
)
genotype <- manifest[manifest$scenario == "sensitivity_exclude_genotype_informed_carriers", , drop = FALSE]
combined <- manifest[manifest$scenario == "sensitivity_exclude_all_treatment_moderated", , drop = FALSE]
stopifnot(
  genotype$fixed_model_studies == 23L, genotype$fixed_model_arms == 69L,
  genotype$fixed_model_n == 10592L, genotype$fixed_model_r == 3111L,
  identical(genotype$dropped_single_arm_studies, "Bozina_2021|Wigle_2021_prospective"),
  combined$fixed_model_studies == 22L, combined$fixed_model_arms == 65L,
  combined$fixed_model_n == 10106L, combined$fixed_model_r == 3068L,
  identical(combined$dropped_single_arm_studies, "Bozina_2021|Wigle_2021_prospective")
)
stopifnot(nrow(checks) > 0L, all(checks$passed))
stopifnot(all(c(
  "scenario_count", "unique_study_treatment", "binomial_valid", "fixed_study_identifiable",
  "all_nodes_present", "geometry_rank", "row_delta_contract", "input_checksums",
  "seed_schedule", "protected_artifacts"
) %in% checks$check))
stopifnot(nrow(protected) == 6L, all(protected$matches_expected), all(protected$sha256_before == protected$sha256_after))
stopifnot(nrow(seeds) == 56L, !anyDuplicated(seeds$seed), identical(sort(unique(seeds$chain)), 1:4))
stopifnot(nrow(geometry) == 12L, all(geometry$connected), all(geometry$rank_A == 4L))
stopifnot(all(geometry$q == geometry$q_from_rank), all(geometry$q == geometry$q_from_arm_rank))
stopifnot(result$geometry$connectedness$connected, result$geometry$q > 0L)

cat("All fixed-study checkpoint tests passed.\n")
