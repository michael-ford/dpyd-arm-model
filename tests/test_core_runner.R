source("R/common/data_utils.R")
source("R/checksums.R")
source("R/core_scenarios.R")
source("R/fixed_study_model.R")
source("R/core_checkpoint.R")
source("R/core_runner.R")

project_root <- normalizePath(getwd(), mustWork = TRUE)
checkpoint <- file.path(project_root, "results", "analysis", "core")
fixture_parent <- tempfile("fixed-study-runner-")
dir.create(fixture_parent)
fixture <- file.path(fixture_parent, basename(checkpoint))
dir.create(fixture)
checkpoint_entries <- list.files(checkpoint, full.names = TRUE, all.files = TRUE, no.. = TRUE)
checkpoint_entries <- checkpoint_entries[!basename(checkpoint_entries) %in% c("runs", "MCMC_STARTED")]
stopifnot(length(checkpoint_entries) > 0L, all(file.copy(checkpoint_entries, fixture, recursive = TRUE)))
publication_copy_checkpoint_inputs(fixture)
on.exit(unlink(fixture_parent, recursive = TRUE), add = TRUE)
metadata <- read.csv(file.path(fixture, "checkpoint_metadata.csv"), stringsAsFactors = FALSE)

before_files <- list.files(fixture, recursive = TRUE, full.names = TRUE)
before_hashes <- vapply(before_files, sha256_file, character(1))
preflight <- core_run_scenario(
  project_root = project_root, output_root = fixture, scenario_name = "primary",
  input_manifest_sha256 = metadata$input_manifest_sha256,
  analysis_commit = "0000000000000000000000000000000000000000",
  checkpoint_commit = "0000000000000000000000000000000000000000",
  docker_image = metadata$docker_image, preflight = TRUE, check_git = FALSE
)
stopifnot(!dir.exists(file.path(fixture, "runs")), !file.exists(file.path(fixture, "MCMC_STARTED")))
after_files <- list.files(fixture, recursive = TRUE, full.names = TRUE)
stopifnot(identical(before_hashes, vapply(after_files, sha256_file, character(1))))
stopifnot(identical(preflight$prepared$single_arm_policy, "drop"))

# This deliberately tiny test-only fit exercises the runner's attempt capture
# without permitting a reduced configuration through the production entrypoint.
tiny_config <- list(n.adapt = 20L, n.burnin = 20L, n.iter = 40L, n.thin = 1L, n.chains = 4L)
tiny_error <- tryCatch({
  core_run_scenario(
    project_root = project_root, output_root = fixture, scenario_name = "primary",
    input_manifest_sha256 = metadata$input_manifest_sha256,
    analysis_commit = "0000000000000000000000000000000000000000",
    checkpoint_commit = "0000000000000000000000000000000000000000",
    docker_image = metadata$docker_image, config = tiny_config, test_only = TRUE,
    check_git = FALSE
  )
  NULL
}, error = identity)
attempt_one <- file.path(fixture, "runs", "primary", "attempt_001")
stopifnot(
  inherits(tiny_error, "error"), dir.exists(attempt_one),
  file.exists(file.path(fixture, "MCMC_STARTED")),
  file.exists(file.path(attempt_one, "run_metadata.csv")),
  file.exists(file.path(attempt_one, "fixed_study_mcmc_samples.rds")),
  file.exists(file.path(attempt_one, "artifact_checksums.csv")),
  file.exists(file.path(attempt_one, "RUN_NONCONVERGED"))
)

# Equal aggregate totals are insufficient: every prepared source arm must match
# the checkpoint's transformed input in deterministic study/treatment order.
mutated_raw <- read.csv(preflight$verified$raw_path, stringsAsFactors = FALSE)
slack_rows <- which(mutated_raw$N > mutated_raw$R)
stopifnot(length(slack_rows) >= 2L)
mutated_raw$N[slack_rows[1]] <- mutated_raw$N[slack_rows[1]] + 1L
mutated_raw$N[slack_rows[2]] <- mutated_raw$N[slack_rows[2]] - 1L
stopifnot(sum(mutated_raw$N) == sum(read.csv(preflight$verified$raw_path, stringsAsFactors = FALSE)$N))
write.csv(mutated_raw, preflight$verified$raw_path, row.names = FALSE)
row_mutation_error <- tryCatch({ core_assert_raw_fit_contract(preflight$verified); NULL }, error = identity)
stopifnot(
  inherits(row_mutation_error, "error"),
  grepl("row-for-row", conditionMessage(row_mutation_error), fixed = TRUE)
)

manifest <- read.csv(file.path(fixture, "scenario_input_manifest.csv"), stringsAsFactors = FALSE)
raw_path <- publication_input_path(fixture, manifest$raw_input_path[manifest$scenario == "sensitivity_exclude_genotype_informed_carriers"])
prepared <- prepare_fixed_study_data(read.csv(raw_path, stringsAsFactors = FALSE), core_treatment_map(), core_treatment_names(), single_arm = "drop")
stopifnot(identical(prepared$dropped_single_arm_studies, c("Bozina_2021", "Wigle_2021_prospective")))

rerun_error <- tryCatch({ core_next_attempt(file.path(fixture, "runs"), "primary"); NULL }, error = identity)
stopifnot(inherits(rerun_error, "error"), grepl("RERUN_REASON", conditionMessage(rerun_error)))

cat("Fixed-study runner preflight tests passed.\n")
