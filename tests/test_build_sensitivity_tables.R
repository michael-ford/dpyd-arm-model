source("R/checksums.R")
source("R/build_sensitivity_tables.R")

root <- tempfile("fixed-study-final-tables-")
core_root <- file.path(root, "core")
dependent_root <- file.path(root, "dependent")
output_dir <- file.path(root, "tables")
dir.create(file.path(core_root, "tables"), recursive = TRUE)
dir.create(core_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dependent_root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

core_aliases <- c("sensitivity_exclude_cevik", "sensitivity_exclude_pk_informed_morel", "sensitivity_exclude_onesti")
core_scenarios <- c("primary", core_aliases, paste0("core_", 1:8))
core <- data.frame(
  scenario = core_scenarios, definition = core_scenarios, studies = 25L, arms = 79L,
  participants = 12486L, events = 3821L,
  HapB3_OR = seq(1.65, 2.05, length.out = 12), lower_95_CrI = seq(0.94, 1.10, length.out = 12),
  upper_95_CrI = seq(2.8, 3.2, length.out = 12), stringsAsFactors = FALSE
)
core_diagnostics <- data.frame(
  scenario = core_scenarios, max_psrf_upper_95 = 1.01, min_effective_size = 5000,
  stringsAsFactors = FALSE
)
write.csv(core, file.path(core_root, "tables", "table_core_hapb3_sensitivities.csv"), row.names = FALSE)
write.csv(core_diagnostics, file.path(core_root, "diagnostic_extrema.csv"), row.names = FALSE)

manifest_rows <- list()
key_rows <- list(data.frame(
  scenario = "primary", analysis_family = "primary", comparison = "HapB3 vs WT",
  or_median = 1.73, or_lower_95_cri = 1.0004, or_upper_95_cri = 2.91,
  stringsAsFactors = FALSE
))
selected_rows <- list()
alias_rows <- list()
add_scenario <- function(scenario, family, status = "planned", target = NA_character_, alias_of = NA_character_, comparisons = "HapB3 vs WT") {
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
    scenario = scenario, analysis_family = family, run_status = status, definition = scenario,
    target_study = target, alias_of = alias_of, studies = 24L, arms = 75L, total_n = 12000L, total_r = 3700L,
    fixed_model_studies = 24L, fixed_model_arms = 75L, fixed_model_n = 12000L, fixed_model_r = 3700L,
    stringsAsFactors = FALSE
  )
  for (comparison in comparisons) {
    key_rows[[length(key_rows) + 1L]] <<- data.frame(
      scenario = scenario, analysis_family = family, comparison = comparison,
      or_median = 1.75, or_lower_95_cri = 0.98, or_upper_95_cri = 3.05,
      stringsAsFactors = FALSE
    )
  }
  if (status == "planned") {
    selected_rows[[length(selected_rows) + 1L]] <<- data.frame(
      scenario = scenario, max_psrf_upper_95 = 1.01, min_effective_size = 4000,
      stringsAsFactors = FALSE
    )
  } else if (status == "alias") {
    alias_rows[[length(alias_rows) + 1L]] <<- data.frame(
      scenario = scenario, analysis_family = family, definition = scenario, attempt = NA_integer_,
      attempt_path = NA_character_, status = "alias", execution_scenario = alias_of,
      stringsAsFactors = FALSE
    )
  }
}
add_scenario("genotyping_wt_split_loose", "genotyping", comparisons = c("HapB3 vs WT_Clean", "HapB3 vs WT_Biased"))
add_scenario("genotyping_wt_split_strict", "genotyping", comparisons = c("HapB3 vs WT_Clean", "HapB3 vs WT_Biased"))
add_scenario("sparse_exclude_n2", "sparse")
manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
  scenario = "sparse_exclude_n5_missing_13", analysis_family = "sparse", run_status = "infeasible",
  definition = "infeasible", target_study = NA_character_, alias_of = NA_character_,
  studies = NA, arms = NA, total_n = NA, total_r = NA,
  fixed_model_studies = NA, fixed_model_arms = NA, fixed_model_n = NA, fixed_model_r = NA,
  stringsAsFactors = FALSE
)
for (index in 1:6) add_scenario(paste0("ld_", index), "ld", target = paste0("LD_", ceiling(index / 2)))
for (index in 1:23) {
  status <- if (index <= 3L) "alias" else "planned"
  add_scenario(
    paste0("loo_", index), "loo", status = status, target = paste0("Source_", index),
    alias_of = if (status == "alias") core_aliases[[index]] else NA_character_
  )
}
manifest <- do.call(rbind, manifest_rows)
key <- do.call(rbind, key_rows)
selected <- do.call(rbind, selected_rows)
aliases <- do.call(rbind, alias_rows)
infeasible <- data.frame(
  scenario = "sparse_exclude_n5_missing_13", analysis_family = "sparse",
  definition = "infeasible", infeasible_reason = "*13 absent", stringsAsFactors = FALSE
)
write.csv(manifest, file.path(dependent_root, "scenario_input_manifest.csv"), row.names = FALSE)
write.csv(key, file.path(dependent_root, "dependent_key_comparison_summary.csv"), row.names = FALSE)
write.csv(selected, file.path(dependent_root, "dependent_selected_run_index.csv"), row.names = FALSE)
write.csv(aliases, file.path(dependent_root, "dependent_alias_index.csv"), row.names = FALSE)
write.csv(infeasible, file.path(dependent_root, "dependent_infeasible_records.csv"), row.names = FALSE)

result <- build_sensitivity_tables(
  project_root = root, core_root = core_root, dependent_root = dependent_root, output_dir = output_dir
)
stopifnot(
  nrow(result$core) == 12L, nrow(result$wt_split) == 4L, nrow(result$ld) == 6L,
  nrow(result$sparse) == 1L, nrow(result$loo) == 23L, nrow(result$infeasible) == 1L,
  nrow(result$all) == 46L,
  !any(result$dependent$execution_status == "planned"),
  all(result$dependent$execution_status[result$dependent$execution_status != "alias"] == "complete"),
  all(file.exists(file.path(output_dir, c(
    "table_s2_all_final_sensitivities.csv", "table_s8_wt_split.csv", "table_s10_ld_bounds.csv",
    "table_supported_sparse_data.csv", "table_grouped_source_influence.csv",
    "table_infeasible_sensitivity_records.csv", "final_table_input_checksums.csv",
    "SENSITIVITY_TABLES.md"
  ))))
)

cat("Fixed-study final table tests passed.\n")
