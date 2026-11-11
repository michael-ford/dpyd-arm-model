source("R/common/data_utils.R")
source("R/checksums.R")
source("R/core_scenarios.R")
source("R/fixed_study_model.R")
source("R/network_geometry.R")
source("R/design_by_treatment.R")
source("R/core_checkpoint.R")
source("R/core_runner.R")
source("R/inconsistency_runner.R")

project_root <- normalizePath(getwd(), mustWork = TRUE)
checkpoint <- file.path(project_root, "results", "analysis", "core")
fixture_parent <- tempfile("fixed-study-inconsistency-runner-")
dir.create(fixture_parent)
fixture <- file.path(fixture_parent, basename(checkpoint))
dir.create(fixture)
checkpoint_entries <- list.files(checkpoint, full.names = TRUE, all.files = TRUE, no.. = TRUE)
checkpoint_entries <- checkpoint_entries[!basename(checkpoint_entries) %in% c("runs", "MCMC_STARTED")]
stopifnot(length(checkpoint_entries) > 0L, all(file.copy(checkpoint_entries, fixture, recursive = TRUE)))
publication_copy_checkpoint_inputs(fixture)
metadata <- read.csv(file.path(fixture, "checkpoint_metadata.csv"), stringsAsFactors = FALSE)
on.exit(unlink(fixture_parent, recursive = TRUE), add = TRUE)

before_files <- list.files(fixture, recursive = TRUE, full.names = TRUE)
before_hashes <- vapply(before_files, sha256_file, character(1))
preflight_sd2 <- run_inconsistency(
  project_root = project_root, output_root = fixture, scenario_name = "inconsistency_sd2",
  input_manifest_sha256 = metadata$input_manifest_sha256,
  analysis_commit = "0000000000000000000000000000000000000000",
  checkpoint_commit = "0000000000000000000000000000000000000000",
  docker_image = metadata$docker_image, preflight = TRUE, check_git = FALSE
)
sd1_seeds <- inconsistency_seed_schedule(fixture, "inconsistency_sd1")
after_files <- list.files(fixture, recursive = TRUE, full.names = TRUE)
stopifnot(
  identical(before_hashes, vapply(after_files, sha256_file, character(1))),
  !dir.exists(file.path(fixture, "runs", "inconsistency_sd2")),
  preflight_sd2$gamma_prior_sd == 2, inconsistency_sd("inconsistency_sd1") == 1,
  identical(preflight_sd2$seeds$scenario_name, rep("inconsistency_sd2", 4L)),
  identical(sd1_seeds$scenario_name, rep("inconsistency_sd1", 4L)),
  identical(as.integer(preflight_sd2$seeds$chain), 1:4),
  identical(as.integer(sd1_seeds$chain), 1:4)
)

# The copied-checkpoint preflight above exercised the real protected-artifact
# verifier. Reuse its verified, immutable ledger for the tiny attempt below so
# this test does not repeatedly hash unrelated historical MCMC directories.
core_collect_protected <- function() preflight_sd2$verified$protected_before

# The fake is only the MCMC boundary.  It emits a real R warning so attempt
# capture is exercised without starting a production JAGS chain.
tiny_fit <- function(scenario_name, gamma_prior_sd, config, chain_seeds, ...) {
  stopifnot(identical(scenario_name, "inconsistency_sd2"), gamma_prior_sd == 2,
            identical(as.integer(chain_seeds$chain), 1:4), config$n.chains == 4L)
  warning("test-only simulated JAGS warning")
  draws <- matrix(0, nrow = config$n.iter * config$n.chains / config$n.thin, ncol = 1L,
                  dimnames = list(NULL, "gamma[1]"))
  list(
    samples = draws,
    treatment_effects = data.frame(treatment = "WT", median = 1),
    pairwise_effects = data.frame(comparison = "WT vs WT", median = 1),
    heterogeneity = data.frame(parameter = "tau_arm", median = 0),
    design_interactions = data.frame(interaction = "design_interaction[1]", median = 0),
    diagnostic_table = data.frame(parameter = "gamma[1]", passed = FALSE),
    diagnostics = list(converged = FALSE, max_psrf_upper_95 = 1.2, min_effective_size = 10),
    primary_diagnostics = list(converged = TRUE), inconsistency_diagnostics = list(converged = FALSE),
    information_gate = list(passed = TRUE, classification = "reliably_estimable",
                            maximum_posterior_to_prior_variance_ratio = 0.5,
                            variance_ratio_threshold = 0.64),
    omnibus = list(valid = TRUE, classification = "no_evidence_of_inconsistency",
                   method = "test-only empirical band", probability = 0.95,
                   simultaneous_critical_value = 2),
    residual_deviance_arms = data.frame(arm = 1L, mean = 1),
    residual_deviance_summary = data.frame(mean = 1, per_arm_mean = 1),
    dic = list(DIC = 2, descriptive_only = TRUE),
    warnings = "Gamma diagnostic did not converge.", test_only = TRUE
  )
}
tiny_config <- list(n.adapt = 2L, n.burnin = 2L, n.iter = 10L, n.thin = 1L, n.chains = 4L)
tiny_error <- tryCatch({
  run_inconsistency(
    project_root = project_root, output_root = fixture, scenario_name = "inconsistency_sd2",
    input_manifest_sha256 = metadata$input_manifest_sha256,
    analysis_commit = "0000000000000000000000000000000000000000",
    checkpoint_commit = "0000000000000000000000000000000000000000",
    docker_image = metadata$docker_image, config = tiny_config, test_only = TRUE,
    fit_function = tiny_fit, check_git = FALSE
  )
  NULL
}, error = identity)
attempt_one <- file.path(fixture, "runs", "inconsistency_sd2", "attempt_001")
attempt_metadata <- read.csv(file.path(attempt_one, "run_metadata.csv"), stringsAsFactors = FALSE)
stopifnot(
  inherits(tiny_error, "error"), dir.exists(attempt_one),
  file.exists(file.path(attempt_one, "RUN_NONCONVERGED")),
  file.exists(file.path(attempt_one, "design_by_treatment_mcmc_samples.rds")),
  file.exists(file.path(attempt_one, "diagnostics.rds")),
  file.exists(file.path(attempt_one, "information_gate.rds")),
  file.exists(file.path(attempt_one, "omnibus.rds")),
  file.exists(file.path(attempt_one, "artifact_checksums.csv")),
  grepl("test-only simulated JAGS warning", paste(read.csv(file.path(attempt_one, "r_jags_warnings.csv"))$warning, collapse = " "), fixed = TRUE),
  grepl("Gamma diagnostic did not converge", paste(read.csv(file.path(attempt_one, "gate_warnings.csv"))$warning, collapse = " "), fixed = TRUE),
  identical(attempt_metadata$status, "nonconverged"),
  !attempt_metadata$combined_converged,
  attempt_metadata$information_gate_passed,
  attempt_metadata$empirical_omnibus_valid,
  attempt_metadata$posterior_draws == 40L,
  attempt_metadata$gamma_prior_sd == 2
)

rerun_error <- tryCatch({
  inconsistency_next_attempt(file.path(fixture, "runs"), "inconsistency_sd2")
  NULL
}, error = identity)
stopifnot(inherits(rerun_error, "error"), grepl("RERUN_REASON", conditionMessage(rerun_error)))

cat("Fixed-study inconsistency runner preflight and immutable-attempt tests passed.\n")
