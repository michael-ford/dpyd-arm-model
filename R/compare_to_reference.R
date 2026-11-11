#!/usr/bin/env Rscript
source("R/checksums.R")

compare_to_reference <- function(scenario, attempt_path = NULL) {
  publication_verify_inputs()
  dependent <- read.csv("results/analysis/dependent/scenario_input_manifest.csv", stringsAsFactors = FALSE)
  row <- dependent[dependent$scenario == scenario, , drop = FALSE]
  execution <- if (nrow(row) && row$run_status == "alias") row$alias_of else scenario
  if (nrow(row) && row$run_status == "infeasible") stop(row$infeasible_reason, call. = FALSE)
  if (is.null(attempt_path)) {
    paths <- sort(list.dirs(file.path("runs", execution), recursive = FALSE, full.names = TRUE))
    paths <- paths[grepl("^attempt_[0-9]{3}$", basename(paths))]
    if (!length(paths)) stop("No rerun found for ", execution, call. = FALSE)
    attempt_path <- tail(paths, 1L)
  }
  read <- function(name) read.csv(file.path(attempt_path, name), stringsAsFactors = FALSE, check.names = FALSE)
  meta <- read("run_metadata.csv")
  diagnostic_only <- execution %in% c("inconsistency_sd1", "inconsistency_sd2")
  if (nrow(meta) != 1L || meta$scenario != execution || isTRUE(meta$test_only) ||
      (!diagnostic_only && (meta$status != "complete" || !isTRUE(meta$converged))) ||
      (diagnostic_only && !isTRUE(meta$combined_converged))) {
    stop("Comparison requires a completed production fit with passed MCMC gates.", call. = FALSE)
  }
  checks <- read("artifact_checksums.csv")
  actual <- vapply(file.path(attempt_path, checks$artifact), sha256_file, character(1))
  if (!identical(unname(actual), checks$sha256)) stop("Rerun artifact checksum mismatch.", call. = FALSE)
  differences <- list()
  add <- function(label, actual, reference) {
    if (!length(actual) || length(actual) != length(reference) || any(!is.finite(c(actual, reference)))) {
      stop("Missing, unmatched, or non-finite values for ", label, call. = FALSE)
    }
    differences[[length(differences) + 1L]] <<- data.frame(
      quantity = label, compared_values = length(actual), max_absolute_difference = max(abs(actual - reference)))
  }
  cols <- c("or_median", "or_lower_95_cri", "or_upper_95_cri")
  core <- read.csv("results/analysis/core/all_treatment_posterior_summary.csv", stringsAsFactors = FALSE)
  ref <- core[core$scenario == execution, , drop = FALSE]
  if (nrow(ref)) {
    rerun <- read("posterior_treatment_effects.csv")
    if (anyDuplicated(rerun$treatment) || !setequal(ref$treatment, rerun$treatment)) stop("Treatment nodes differ.", call. = FALSE)
    rerun <- rerun[match(ref$treatment, rerun$treatment), , drop = FALSE]
    for (col in cols) add(col, rerun[[col]], ref[[col]])
    tau <- read.csv("results/analysis/core/heterogeneity_fit_summary.csv", stringsAsFactors = FALSE)
    tau <- tau[tau$scenario == execution, , drop = FALSE]
    got <- read("posterior_heterogeneity.csv")
    got <- got[match(tau$parameter, got$parameter), , drop = FALSE]
    for (i in seq_along(c("median", "lower_95_cri", "upper_95_cri"))) {
      col <- c("median", "lower_95_cri", "upper_95_cri")[[i]]
      add(paste0("tau_", col), got[[col]], tau[[paste0("heterogeneity_", col)]])
    }
    extrema <- read.csv("results/analysis/core/diagnostic_extrema.csv", stringsAsFactors = FALSE)
  } else if (!diagnostic_only) {
    ref <- read.csv("results/analysis/dependent/dependent_key_comparison_summary.csv", stringsAsFactors = FALSE)
    ref <- ref[ref$scenario == scenario, , drop = FALSE]
    if (!nrow(ref)) stop("No reference posterior summary for ", scenario, call. = FALSE)
    got <- read("posterior_pairwise_effects.csv")
    keys <- function(x) paste(x$numerator, x$denominator, sep = "|")
    matched <- match(keys(ref), keys(got))
    if (anyNA(matched)) stop("Missing reference pairwise contrasts.", call. = FALSE)
    for (col in cols) add(col, got[[col]][matched], ref[[col]])
    cat("Reference coverage: dependent OR/CrI includes the committed HapB3 contrasts only.\n",
        "Tau comparison unavailable: V4 contains no dependent tau summary.\n", sep = "")
    extrema <- read.csv("results/analysis/dependent/dependent_selected_run_index.csv", stringsAsFactors = FALSE)
  } else {
    cat("Reference coverage: inconsistency diagnostics only; OR/CrI and tau summaries were not packaged.\n")
    extrema <- read.csv("results/analysis/inconsistency/inconsistency_run_gate_summary.csv", stringsAsFactors = FALSE)
  }
  extrema <- extrema[extrema$scenario == execution, , drop = FALSE]
  if (nrow(extrema) != 1L) stop("Reference diagnostic extrema are missing.", call. = FALSE)
  got <- read("diagnostics.csv")
  add("PSRF_upper_95_max", max(got$psrf_upper_95), extrema$max_psrf_upper_95)
  add("ESS_min", min(got$effective_size), extrema$min_effective_size)
  result <- do.call(rbind, differences)
  cat("Scenario:", scenario, "\nRerun:", attempt_path, "\n")
  print(result, row.names = FALSE, digits = 12)
  cat("PSRF/ESS comparisons use extrema because per-parameter reference diagnostics were not packaged.\n")
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || length(args) > 2L) stop("Usage: Rscript R/compare_to_reference.R <scenario> [attempt-directory]", call. = FALSE)
  compare_to_reference(args[[1]], if (length(args) == 2L) args[[2]] else NULL)
}
