# Reference dependent-runner checks: copied-checkpoint preflight and tiny attempt.

source("R/common/data_utils.R")
source("R/arm_rows.R")
source("R/scenario_helpers.R")
source("R/checksums.R")
source("R/core_scenarios.R")
source("R/fixed_study_model.R")
source("R/core_checkpoint.R")
source("R/core_runner.R")
source("R/dependent_scenarios.R")
source("R/dependent_checkpoint.R")
source("R/dependent_runner.R")

project_root <- normalizePath(getwd(), mustWork = TRUE)
parent_root <- dependent_parent_root(project_root)
parent_metadata <- read.csv(file.path(parent_root, "checkpoint_metadata.csv"), stringsAsFactors = FALSE)
commits <- list(
  analysis = "0000000000000000000000000000000000000000",
  dependent_checkpoint = "0000000000000000000000000000000000000000",
  parent_checkpoint = "0000000000000000000000000000000000000000"
)

# Generate once, then preflight a copied checkpoint so the preflight cannot
# accidentally depend on the staging path or write to the reference fixture.
staging_parent <- tempfile("fixed-study-dependent-staging-")
fixture_parent <- tempfile("fixed-study-dependent-runner-")
dir.create(staging_parent)
dir.create(fixture_parent)
on.exit(unlink(c(staging_parent, fixture_parent), recursive = TRUE), add = TRUE)
staging <- file.path(staging_parent, "dependent_checkpoint_fixture")
checkpoint_result <- generate_dependent_checkpoint(
  project_root = project_root, output_root = staging,
  metadata = list(
    analysis_commit = commits$analysis,
    docker_image = dependent_pinned_docker_image(),
    input_manifest_sha256 = parent_metadata$input_manifest_sha256
  ),
  ld_second_target_policy = "wigle_prospective"
)
stopifnot(file.copy(staging, fixture_parent, recursive = TRUE))
fixture <- file.path(fixture_parent, basename(staging))

before_files <- list.files(fixture, recursive = TRUE, full.names = TRUE)
before_hashes <- vapply(before_files, sha256_file, character(1))
preflight <- run_dependent_scenario(
  project_root = project_root, output_root = fixture,
  scenario_name = "genotyping_wt_split_loose",
  input_manifest_sha256 = parent_metadata$input_manifest_sha256,
  analysis_commit = commits$analysis, checkpoint_commit = commits$dependent_checkpoint,
  parent_checkpoint_commit = commits$parent_checkpoint,
  docker_image = dependent_pinned_docker_image(),
  preflight = TRUE, check_git = FALSE
)
after_files <- list.files(fixture, recursive = TRUE, full.names = TRUE)
stopifnot(
  identical(before_hashes, vapply(after_files, sha256_file, character(1))),
  !dir.exists(file.path(fixture, "runs")), !file.exists(file.path(fixture, "MCMC_STARTED")),
  preflight$treatment_contract$expected_nodes == 6L,
  identical(preflight$treatment_contract$names, dependent_split_wt_names()),
  identical(as.integer(preflight$seeds$chain), 1:4),
  preflight$prepared$jags_data$Ntrt == 6L
)

# The copied-checkpoint preflight above exercised the real protected-artifact
# verifier. Reuse that immutable ledger for the remaining tiny checks so they
# do not repeatedly hash unrelated historical MCMC directories.
core_collect_protected <- function() preflight$protected_before

unified_preflight <- run_dependent_scenario(
  project_root = project_root, output_root = fixture,
  scenario_name = "sparse_exclude_n2",
  input_manifest_sha256 = parent_metadata$input_manifest_sha256,
  analysis_commit = commits$analysis, checkpoint_commit = commits$dependent_checkpoint,
  parent_checkpoint_commit = commits$parent_checkpoint,
  docker_image = dependent_pinned_docker_image(),
  preflight = TRUE, check_git = FALSE
)
stopifnot(
  unified_preflight$treatment_contract$expected_nodes == 5L,
  identical(unified_preflight$treatment_contract$names, dependent_treatment_names()),
  unified_preflight$prepared$jags_data$Ntrt == 5L
)

alias_error <- tryCatch({
  run_dependent_scenario(
    project_root = project_root, output_root = fixture, scenario_name = "loo_exclude_Cevik_2022",
    input_manifest_sha256 = parent_metadata$input_manifest_sha256,
    analysis_commit = commits$analysis, checkpoint_commit = commits$dependent_checkpoint,
    parent_checkpoint_commit = commits$parent_checkpoint,
    docker_image = dependent_pinned_docker_image(), preflight = TRUE, check_git = FALSE
  )
  NULL
}, error = identity)
infeasible_error <- tryCatch({
  run_dependent_scenario(
    project_root = project_root, output_root = fixture, scenario_name = "sparse_exclude_n5_missing_13",
    input_manifest_sha256 = parent_metadata$input_manifest_sha256,
    analysis_commit = commits$analysis, checkpoint_commit = commits$dependent_checkpoint,
    parent_checkpoint_commit = commits$parent_checkpoint,
    docker_image = dependent_pinned_docker_image(), preflight = TRUE, check_git = FALSE
  )
  NULL
}, error = identity)
stopifnot(
  inherits(alias_error, "error"), inherits(infeasible_error, "error"),
  grepl("Only a planned dependent scenario", conditionMessage(alias_error), fixed = TRUE),
  grepl("Only a planned dependent scenario", conditionMessage(infeasible_error), fixed = TRUE)
)

# The fake is only the MCMC boundary. It captures a real warning and supplies
# real runner-shaped artifacts without starting a JAGS chain.
tiny_fit <- function(scenario_name, treatment_names, config, chain_seeds, ...) {
  stopifnot(
    identical(scenario_name, "genotyping_wt_split_loose"),
    identical(treatment_names, dependent_split_wt_names()),
    identical(as.integer(chain_seeds$chain), 1:4), config$n.chains == 4L
  )
  warning("test-only dependent JAGS warning")
  draws <- matrix(0, nrow = config$n.iter * config$n.chains / config$n.thin, ncol = 1L,
                  dimnames = list(NULL, "d[2]"))
  list(
    samples = draws,
    treatment_effects = data.frame(treatment = "WT_Clean", median = 1),
    pairwise_effects = data.frame(comparison = "WT_Biased vs WT_Clean", median = 1),
    heterogeneity = data.frame(parameter = "tau_arm", median = 0),
    diagnostic_table = data.frame(parameter = "d[2]", passed = FALSE),
    diagnostics = list(converged = FALSE, max_psrf_upper_95 = 1.2, min_effective_size = 10),
    residual_deviance_arms = data.frame(arm = 1L, mean = 1),
    residual_deviance_summary = data.frame(mean = 1, per_arm_mean = 1),
    dic = list(DIC = 2, descriptive_only = TRUE),
    warnings = "PSRF/ESS convergence contract failed."
  )
}
tiny_config <- list(n.adapt = 2L, n.burnin = 2L, n.iter = 10L, n.thin = 1L, n.chains = 4L)
tiny_error <- tryCatch({
  run_dependent_scenario(
    project_root = project_root, output_root = fixture,
    scenario_name = "genotyping_wt_split_loose",
    input_manifest_sha256 = parent_metadata$input_manifest_sha256,
    analysis_commit = commits$analysis, checkpoint_commit = commits$dependent_checkpoint,
    parent_checkpoint_commit = commits$parent_checkpoint,
    docker_image = dependent_pinned_docker_image(),
    config = tiny_config, test_only = TRUE, fit_function = tiny_fit, check_git = FALSE
  )
  NULL
}, error = identity)
attempt <- file.path(fixture, "runs", "genotyping_wt_split_loose", "attempt_001")
attempt_metadata <- read.csv(file.path(attempt, "run_metadata.csv"), stringsAsFactors = FALSE)
stopifnot(
  inherits(tiny_error, "error"), dir.exists(attempt),
  file.exists(file.path(attempt, "RUN_NONCONVERGED")),
  file.exists(file.path(attempt, "fixed_study_mcmc_samples.rds")),
  file.exists(file.path(attempt, "posterior_treatment_effects.csv")),
  file.exists(file.path(attempt, "posterior_pairwise_effects.csv")),
  file.exists(file.path(attempt, "posterior_heterogeneity.csv")),
  file.exists(file.path(attempt, "diagnostics.csv")),
  file.exists(file.path(attempt, "residual_deviance_arms.csv")),
  file.exists(file.path(attempt, "descriptive_dic.csv")),
  file.exists(file.path(attempt, "runtime_warnings.csv")),
  file.exists(file.path(attempt, "runtime_code_hashes_before.csv")),
  file.exists(file.path(attempt, "runtime_code_hashes_after.csv")),
  file.exists(file.path(attempt, "protected_artifacts_before.csv")),
  file.exists(file.path(attempt, "protected_artifacts_after.csv")),
  file.exists(file.path(attempt, "sessionInfo.txt")),
  file.exists(file.path(attempt, "artifact_checksums.csv")),
  grepl("test-only dependent JAGS warning",
        paste(read.csv(file.path(attempt, "runtime_warnings.csv"))$warning, collapse = " "), fixed = TRUE),
  identical(attempt_metadata$status, "nonconverged"),
  attempt_metadata$posterior_draws == 40L,
  attempt_metadata$expected_posterior_draws == 100000L,
  !attempt_metadata$converged,
  attempt_metadata$n_treatments == 6L
)

rerun_error <- tryCatch({
  core_next_attempt(file.path(fixture, "runs"), "genotyping_wt_split_loose")
  NULL
}, error = identity)
stopifnot(inherits(rerun_error, "error"), grepl("RERUN_REASON", conditionMessage(rerun_error)))

cat("Fixed-study dependent runner preflight and immutable-attempt tests passed.\n")
