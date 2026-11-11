source("R/checksums.R")
source("R/render_sensitivity_figures.R")

root <- tempfile("fixed-study-dependent-figures-")
core_root <- file.path(root, "core")
dependent_root <- file.path(root, "dependent")
output_dir <- file.path(root, "figures")
dir.create(core_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dependent_root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

core_scenarios <- names(core_sensitivity_labels)
core <- data.frame(
  scenario = core_scenarios, definition = unname(core_sensitivity_labels),
  or_median = seq(1.65, 2.09, length.out = 12),
  or_lower_95_cri = seq(0.92, 1.11, length.out = 12),
  or_upper_95_cri = seq(2.8, 3.3, length.out = 12),
  stringsAsFactors = FALSE
)
core_index <- data.frame(
  scenario = core_scenarios, status = "complete", converged = TRUE,
  all_diagnostics_passed = TRUE, artifact_checksums_match = TRUE,
  stringsAsFactors = FALSE
)
write.csv(core, file.path(core_root, "hapb3_sensitivity_summary.csv"), row.names = FALSE)
write.csv(core_index, file.path(core_root, "core_run_index.csv"), row.names = FALSE)

manifest_rows <- list()
key_rows <- list(data.frame(
  scenario = "primary", analysis_family = "primary", numerator = "HapB3",
  denominator = "WT", comparison = "HapB3 vs WT", or_median = 1.73,
  or_lower_95_cri = 1.0004, or_upper_95_cri = 2.91, stringsAsFactors = FALSE
))
add_manifest <- function(scenario, family, status = "planned", target = NA_character_) {
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
    scenario = scenario, analysis_family = family, run_status = status,
    definition = scenario, target_study = target, stringsAsFactors = FALSE
  )
}
add_key <- function(scenario, family, denominator = "WT", value = 1.7) {
  key_rows[[length(key_rows) + 1L]] <<- data.frame(
    scenario = scenario, analysis_family = family, numerator = "HapB3", denominator = denominator,
    comparison = paste("HapB3 vs", denominator), or_median = value,
    or_lower_95_cri = value * 0.55, or_upper_95_cri = value * 1.75,
    stringsAsFactors = FALSE
  )
}
for (classification in c("loose", "strict")) {
  scenario <- paste0("genotyping_wt_split_", classification)
  add_manifest(scenario, "genotyping")
  add_key(scenario, "genotyping", "WT_Clean", 1.8)
  add_key(scenario, "genotyping", "WT_Biased", 1.7)
}
add_manifest("sparse_exclude_n2", "sparse")
add_key("sparse_exclude_n2", "sparse", value = 1.68)
add_manifest("sparse_exclude_n5_missing_13", "sparse", "infeasible")
for (index in 1:6) {
  scenario <- paste0(if (index %% 2L) "ld_best_" else "ld_worst_", index)
  add_manifest(scenario, "ld", target = paste0("LD_target_", ceiling(index / 2)))
  add_key(scenario, "ld", value = 1.65 + index / 50)
}
for (index in 1:23) {
  scenario <- paste0("loo_", index)
  status <- if (index <= 3L) "alias" else "planned"
  add_manifest(scenario, "loo", status, paste0("Source_", index))
  add_key(scenario, "loo", value = 1.55 + index / 30)
}
manifest <- do.call(rbind, manifest_rows)
key <- do.call(rbind, key_rows)
write.csv(manifest, file.path(dependent_root, "scenario_input_manifest.csv"), row.names = FALSE)
write.csv(key, file.path(dependent_root, "dependent_key_comparison_summary.csv"), row.names = FALSE)

built <- build_dependent_figure_data(core_root, dependent_root)
stopifnot(
  nrow(built$data) == 46L,
  sum(built$data$analysis_family == "core") == 12L,
  sum(built$data$analysis_family == "genotyping") == 4L,
  sum(built$data$analysis_family == "ld") == 6L,
  sum(built$data$analysis_family == "sparse") == 1L,
  sum(built$data$analysis_family == "loo") == 23L,
  identical(built$primary_or, 1.73)
)
rendered <- render_sensitivity_figures(
  project_root = root, core_root = core_root, dependent_root = dependent_root, output_dir = output_dir
)
stopifnot(
  file.exists(rendered$pdf), file.info(rendered$pdf)$size > 1000,
  length(rendered$png) == 4L, all(file.exists(rendered$png)), all(file.info(rendered$png)$size > 1000),
  file.exists(file.path(output_dir, "all_sensitivity_figure_data.csv")),
  file.exists(file.path(output_dir, "dependent_figure_input_checksums.csv")),
  file.exists(file.path(output_dir, "README.md"))
)

cat("Fixed-study dependent figure tests passed.\n")
