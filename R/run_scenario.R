#!/usr/bin/env Rscript
# The public CLI accepts scenario selection only; the production contract is locked.
source("R/common/data_utils.R")
source("R/checksums.R")
source("R/arm_rows.R")
source("R/scenario_helpers.R")
source("R/core_scenarios.R")
source("R/fixed_study_model.R")
source("R/network_geometry.R")
source("R/design_by_treatment.R")
source("R/core_checkpoint.R")
source("R/core_runner.R")
source("R/dependent_scenarios.R")
source("R/dependent_checkpoint.R")
source("R/dependent_runner.R")
source("R/inconsistency_runner.R")

publication_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  root <- normalizePath(getwd(), mustWork = TRUE)
  core <- read.csv("results/analysis/core/scenario_input_manifest.csv", stringsAsFactors = FALSE)
  dependent <- read.csv("results/analysis/dependent/scenario_input_manifest.csv", stringsAsFactors = FALSE)
  all <- c(core$scenario, dependent$scenario[dependent$run_status == "planned"], inconsistency_runs())
  if (identical(args, "--list")) { cat(paste(all, collapse = "\n"), "\n"); return(invisible(all)) }
  preflight <- "--preflight" %in% args
  args <- args[args != "--preflight"]
  if (length(args) != 1L) stop("Usage: Rscript R/run_scenario.R <scenario|--all|--list> [--preflight]", call. = FALSE)
  scenarios <- if (args == "--all") all else args
  publication_verify_inputs(root)
  commit <- publication_git_commit()
  manifest_hash <- sha256_file("data/input_checksums.csv")
  for (scenario in scenarios) {
    row <- dependent[dependent$scenario == scenario, , drop = FALSE]
    if (nrow(row) && row$run_status == "infeasible") stop(row$infeasible_reason, call. = FALSE)
    selected <- if (nrow(row) && row$run_status == "alias") row$alias_of else scenario
    if (!selected %in% all) stop("Unknown scenario: ", scenario, call. = FALSE)
    common <- list(project_root = root, scenario_name = selected, input_manifest_sha256 = manifest_hash,
                   analysis_commit = commit, checkpoint_commit = commit, docker_image = "runtime-recorded",
                   rerun_reason = Sys.getenv("RERUN_REASON", ""), preflight = preflight)
    if (selected %in% core$scenario) {
      value <- do.call(core_run_scenario, common)
    } else if (selected %in% inconsistency_runs()) {
      value <- tryCatch(do.call(run_inconsistency, common), error = function(e) {
        # An information-gate failure is a diagnostic result, not permission to change the model.
        dirs <- sort(list.dirs(file.path("runs", selected), recursive = FALSE, full.names = TRUE))
        if (!length(dirs)) stop(e)
        latest <- tail(dirs, 1L)
        path <- file.path(latest, "run_metadata.csv")
        if (!file.exists(path)) stop(e)
        meta <- read.csv(path, stringsAsFactors = FALSE)
        if (!isTRUE(meta$combined_converged) || !identical(meta$status, "nonconverged") ||
            isTRUE(meta$information_gate_passed) || !file.exists(file.path(latest, "RUN_NONCONVERGED"))) stop(e)
        message(selected, ": MCMC converged; inconsistency is not reliably estimable (information gate failed).")
        NULL
      })
    } else {
      common$parent_checkpoint_commit <- commit
      value <- do.call(run_dependent_scenario, common)
    }
    if (preflight) cat(selected, ": input, row, seed, and production-contract checks passed.\n", sep = "")
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) publication_main()
