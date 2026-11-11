source("R/common/data_utils.R")
source("R/checksums.R")
source("R/summarize_dependent.R")

scenario_row <- function(scenario, family, status = "planned", alias_of = NA_character_) {
  data.frame(
    scenario = scenario, analysis_family = family, run_status = status,
    definition = scenario, alias_of = alias_of, infeasible_reason = NA_character_,
    treatment_names = if (family == "genotyping") {
      "WT_Clean|WT_Biased|HapB3|2846hetho|2Ahetho|13hetho"
    } else {
      "WT|HapB3|2846hetho|2Ahetho|13hetho"
    },
    raw_input_sha256 = if (status == "planned") paste0(sprintf("%02x", nchar(scenario)), strrep("a", 62)) else NA_character_,
    nma_input_sha256 = if (status == "planned") paste0(sprintf("%02x", nchar(scenario) + 1L), strrep("b", 62)) else NA_character_,
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, c(
  list(
    scenario_row("genotyping_wt_split_loose", "genotyping"),
    scenario_row("genotyping_wt_split_strict", "genotyping"),
    scenario_row("sparse_exclude_n2", "sparse"),
    scenario_row("sparse_exclude_n5_missing_13", "sparse", "infeasible")
  ),
  lapply(paste0("ld_", 1:4), scenario_row, family = "ld"),
  lapply(paste0("loo_", 1:20), scenario_row, family = "loo"),
  list(
    scenario_row("loo_alias_cevik", "loo", "alias", "sensitivity_exclude_cevik"),
    scenario_row("loo_alias_morel", "loo", "alias", "sensitivity_exclude_pk_informed_morel"),
    scenario_row("loo_alias_onesti", "loo", "alias", "sensitivity_exclude_onesti")
  )
))
stopifnot(isTRUE(dependent_validate_inventory(manifest)))
manifest_both <- rbind(manifest, scenario_row("ld_5", "ld"), scenario_row("ld_6", "ld"))
stopifnot(isTRUE(dependent_validate_inventory(manifest_both)))
invalid <- rbind(manifest, scenario_row("ld_5", "ld"))
invalid_error <- tryCatch({ dependent_validate_inventory(invalid); NULL }, error = identity)
stopifnot(inherits(invalid_error, "error"))

fixture <- tempfile("fixed-study-dependent-results-")
dir.create(fixture, recursive = TRUE)
on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
planned <- manifest[manifest$run_status == "planned", , drop = FALSE]
seed_rows <- list()
for (index in seq_len(nrow(planned))) {
  scenario <- planned$scenario[[index]]
  attempt_path <- file.path("runs", scenario, "attempt_001")
  dir.create(file.path(fixture, attempt_path), recursive = TRUE)
  seeds <- as.integer(index * 100L + 1:4)
  seed_rows[[index]] <- data.frame(scenario_name = scenario, chain = 1:4, seed = seeds)
  metadata <- data.frame(
    scenario = scenario, attempt = 1L, attempt_path = attempt_path, status = "complete",
    analysis_family = planned$analysis_family[[index]],
    analysis_commit = strrep("1", 40), checkpoint_commit = strrep("2", 40),
    parent_checkpoint_commit = strrep("3", 40),
    docker_image = "runtime-recorded",
    input_manifest_sha256 = strrep("4a", 32),
    raw_input_sha256 = planned$raw_input_sha256[[index]],
    nma_input_sha256 = planned$nma_input_sha256[[index]],
    n_adapt = 50000L, n_burnin = 125000L, n_iter = 250000L, n_thin = 10L, n_chains = 4L,
    expected_posterior_draws = 100000L, posterior_draws = 100000L,
    chain_seed_1 = seeds[[1]], chain_seed_2 = seeds[[2]], chain_seed_3 = seeds[[3]], chain_seed_4 = seeds[[4]],
    test_only = FALSE, production_mcmc_contract = TRUE, converged = TRUE,
    code_hashes_unchanged = TRUE, protected_artifacts_unchanged = TRUE,
    max_psrf_upper_95 = 1.01, min_effective_size = 5000,
    no_fallback_or_extension = TRUE, deviation = "none",
    stringsAsFactors = FALSE
  )
  write.csv(metadata, file.path(fixture, attempt_path, "run_metadata.csv"), row.names = FALSE)
}
seed_schedule <- do.call(rbind, seed_rows)
inventory <- dependent_attempt_inventory(fixture, manifest)
stopifnot(nrow(inventory$selected) == nrow(planned), nrow(inventory$attempts) == nrow(planned))
checkpoint <- data.frame(input_manifest_sha256 = strrep("4a", 32), stringsAsFactors = FALSE)
stopifnot(isTRUE(dependent_verify_selected_contract(
  inventory$selected, planned, checkpoint, seed_schedule
)))
bad_selected <- inventory$selected
bad_selected$min_effective_size[[1]] <- 999
contract_error <- tryCatch({
  dependent_verify_selected_contract(bad_selected, planned, checkpoint, seed_schedule)
  NULL
}, error = identity)
stopifnot(inherits(contract_error, "error"))

pairwise <- data.frame(
  numerator = c("HapB3", "HapB3", "HapB3"), denominator = c("WT", "WT_Clean", "WT_Biased"),
  or_median = c(1.7, 1.8, 1.6), or_lower_95_cri = c(1.0, 0.8, 1.1),
  or_upper_95_cri = c(2.9, 4.0, 2.7), scenario = c("ordinary", "split", "split"),
  stringsAsFactors = FALSE
)
ordinary <- dependent_key_rows(pairwise, "loo", "ordinary")
split <- dependent_key_rows(pairwise, "genotyping", "split")
stopifnot(nrow(ordinary) == 1L, ordinary$comparison == "HapB3 vs WT", nrow(split) == 2L)

artifact_root <- tempfile("fixed-study-dependent-artifact-")
dir.create(file.path(artifact_root, "runs", "one", "attempt_001"), recursive = TRUE)
on.exit(unlink(artifact_root, recursive = TRUE), add = TRUE)
attempt_dir <- file.path(artifact_root, "runs", "one", "attempt_001")
writeLines("complete", file.path(attempt_dir, "RUN_COMPLETE"))
runtime <- data.frame(item = "model", path = "model", sha256 = strrep("5", 64), stringsAsFactors = FALSE)
write.csv(runtime, file.path(attempt_dir, "runtime_code_hashes_before.csv"), row.names = FALSE)
write.csv(runtime, file.path(attempt_dir, "runtime_code_hashes_after.csv"), row.names = FALSE)
protected <- data.frame(
  artifact = "original", path = "original", kind = "file", expected_sha256 = strrep("6a", 32),
  observed_sha256 = strrep("6a", 32), matches_expected = TRUE, stringsAsFactors = FALSE
)
write.csv(protected, file.path(attempt_dir, "protected_artifacts_after.csv"), row.names = FALSE)
writeLines("payload", file.path(attempt_dir, "payload.txt"))
artifact_index <- data.frame(
  artifact = "payload.txt", size_bytes = file.info(file.path(attempt_dir, "payload.txt"))$size,
  sha256 = sha256_file(file.path(attempt_dir, "payload.txt")), stringsAsFactors = FALSE
)
write.csv(artifact_index, file.path(attempt_dir, "artifact_checksums.csv"), row.names = FALSE)
selected_one <- data.frame(scenario = "one", attempt = 1L, attempt_path = "runs/one/attempt_001")
verified <- dependent_verify_attempt_artifacts(artifact_root, selected_one, protected)
stopifnot(nrow(verified) == 1L, verified$checksum_matches)

cat("Fixed-study dependent result summarizer tests passed.\n")
