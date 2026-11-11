#!/usr/bin/env Rscript

# Build reviewer-ready final sensitivity tables from the verified core and
# dependent summaries. This script never reads posterior chains or runs MCMC.

sensitivity_table_read <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

sensitivity_table_required <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
}

sensitivity_table_md <- function(data, columns) {
  data <- data[, columns, drop = FALSE]
  data[] <- lapply(data, function(value) {
    value <- as.character(value); value[is.na(value)] <- ""; gsub("|", "\\\\|", value, fixed = TRUE)
  })
  c(
    paste0("| ", paste(names(data), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(data)), collapse = " | "), " |"),
    apply(data, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  )
}

build_sensitivity_tables <- function(
  project_root = getwd(),
  core_root = file.path(project_root, "results", "analysis", "core"),
  dependent_root = file.path(project_root, "results", "analysis", "dependent"),
  output_dir = file.path(project_root, "runs", "rebuilt", "tables", "sensitivity")
) {
  core <- sensitivity_table_read(
    if (identical(core_root, file.path(project_root, "results", "analysis", "core"))) file.path(project_root, "results", "tables", "core", "table_core_hapb3_sensitivities.csv") else file.path(core_root, "tables", "table_core_hapb3_sensitivities.csv"), "core sensitivity table"
  )
  core_diagnostics <- sensitivity_table_read(
    file.path(core_root, "diagnostic_extrema.csv"), "core diagnostics"
  )
  key <- sensitivity_table_read(
    file.path(dependent_root, "dependent_key_comparison_summary.csv"), "dependent key comparisons"
  )
  manifest <- sensitivity_table_read(
    file.path(dependent_root, "scenario_input_manifest.csv"), "dependent manifest"
  )
  selected <- sensitivity_table_read(
    file.path(dependent_root, "dependent_selected_run_index.csv"), "dependent selected-run index"
  )
  aliases <- sensitivity_table_read(
    file.path(dependent_root, "dependent_alias_index.csv"), "dependent alias index"
  )
  infeasible <- sensitivity_table_read(
    file.path(dependent_root, "dependent_infeasible_records.csv"), "dependent infeasible records"
  )
  sensitivity_table_required(
    core, c("scenario", "definition", "studies", "arms", "participants", "events", "HapB3_OR", "lower_95_CrI", "upper_95_CrI"),
    "Core sensitivity table"
  )
  sensitivity_table_required(
    key, c("scenario", "analysis_family", "comparison", "or_median", "or_lower_95_cri", "or_upper_95_cri"),
    "Dependent key comparisons"
  )
  sensitivity_table_required(
    manifest, c("scenario", "analysis_family", "run_status", "definition", "target_study", "alias_of", "studies", "arms", "total_n", "total_r"),
    "Dependent manifest"
  )
  sensitivity_table_required(
    selected, c("scenario", "max_psrf_upper_95", "min_effective_size"), "Dependent selected-run index"
  )
  if (nrow(core) != 12L || nrow(infeasible) != 1L || nrow(aliases) != 3L) {
    stop("Final table inputs do not match the bounded core/dependent inventory.", call. = FALSE)
  }

  core_diag <- core_diagnostics[match(core$scenario, core_diagnostics$scenario), , drop = FALSE]
  if (anyNA(core_diag$scenario)) stop("Core sensitivity diagnostics are incomplete.", call. = FALSE)
  core_table <- data.frame(
    family = "core", scenario = core$scenario, definition = core$definition,
    comparison = "HapB3 vs WT", studies = core$studies, arms = core$arms,
    participants = core$participants, events = core$events,
    OR = core$HapB3_OR, lower_95_CrI = core$lower_95_CrI, upper_95_CrI = core$upper_95_CrI,
    max_PSRF_upper_95 = core_diag$max_psrf_upper_95, min_ESS = core_diag$min_effective_size,
    execution_status = "complete", execution_source = core$scenario,
    stringsAsFactors = FALSE
  )

  dependent_key <- key[key$analysis_family != "primary", , drop = FALSE]
  mapped <- manifest[match(dependent_key$scenario, manifest$scenario), , drop = FALSE]
  if (anyNA(mapped$scenario) || any(mapped$run_status == "infeasible")) {
    stop("Dependent key rows do not map to runnable or aliased manifest entries.", call. = FALSE)
  }
  planned_diag <- selected[match(dependent_key$scenario, selected$scenario), , drop = FALSE]
  core_alias_diag <- core_diagnostics[match(mapped$alias_of, core_diagnostics$scenario), , drop = FALSE]
  is_alias <- mapped$run_status == "alias"
  max_psrf <- ifelse(is_alias, core_alias_diag$max_psrf_upper_95, planned_diag$max_psrf_upper_95)
  min_ess <- ifelse(is_alias, core_alias_diag$min_effective_size, planned_diag$min_effective_size)
  if (anyNA(max_psrf) || anyNA(min_ess)) stop("Dependent scenario diagnostics are incomplete.", call. = FALSE)
  dependent_table <- data.frame(
    family = dependent_key$analysis_family, scenario = dependent_key$scenario,
    definition = mapped$definition, comparison = dependent_key$comparison,
    studies = mapped$fixed_model_studies, arms = mapped$fixed_model_arms,
    participants = mapped$fixed_model_n, events = mapped$fixed_model_r,
    OR = dependent_key$or_median, lower_95_CrI = dependent_key$or_lower_95_cri,
    upper_95_CrI = dependent_key$or_upper_95_cri,
    max_PSRF_upper_95 = max_psrf, min_ESS = min_ess,
    # The manifest records whether a scenario was planned before execution.
    # Presence in the selected-run index means the production fit completed;
    # do not carry the pre-run "planned" label into manuscript-facing tables.
    execution_status = ifelse(is_alias, "alias", "complete"),
    execution_source = ifelse(is_alias, mapped$alias_of, mapped$scenario),
    stringsAsFactors = FALSE
  )
  counts <- table(dependent_table$family)
  if (counts[["genotyping"]] != 4L || counts[["sparse"]] != 1L ||
      !counts[["ld"]] %in% c(4L, 6L) || counts[["loo"]] != 23L) {
    stop("Final dependent table family counts are invalid.", call. = FALSE)
  }
  all_sensitivities <- rbind(core_table, dependent_table)
  wt_split <- dependent_table[dependent_table$family == "genotyping", , drop = FALSE]
  ld <- dependent_table[dependent_table$family == "ld", , drop = FALSE]
  sparse <- dependent_table[dependent_table$family == "sparse", , drop = FALSE]
  loo <- dependent_table[dependent_table$family == "loo", , drop = FALSE]
  loo$excluded_source <- manifest$target_study[match(loo$scenario, manifest$scenario)]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(all_sensitivities, file.path(output_dir, "table_s2_all_final_sensitivities.csv"), row.names = FALSE, na = "")
  write.csv(wt_split, file.path(output_dir, "table_s8_wt_split.csv"), row.names = FALSE, na = "")
  write.csv(ld, file.path(output_dir, "table_s10_ld_bounds.csv"), row.names = FALSE, na = "")
  write.csv(sparse, file.path(output_dir, "table_supported_sparse_data.csv"), row.names = FALSE, na = "")
  write.csv(loo, file.path(output_dir, "table_grouped_source_influence.csv"), row.names = FALSE, na = "")
  write.csv(infeasible, file.path(output_dir, "table_infeasible_sensitivity_records.csv"), row.names = FALSE, na = "")

  input_paths <- c(
    core_table = if (identical(core_root, file.path(project_root, "results", "analysis", "core"))) file.path(project_root, "results", "tables", "core", "table_core_hapb3_sensitivities.csv") else file.path(core_root, "tables", "table_core_hapb3_sensitivities.csv"),
    core_diagnostics = file.path(core_root, "diagnostic_extrema.csv"),
    dependent_key = file.path(dependent_root, "dependent_key_comparison_summary.csv"),
    dependent_manifest = file.path(dependent_root, "scenario_input_manifest.csv"),
    selected_runs = file.path(dependent_root, "dependent_selected_run_index.csv"),
    aliases = file.path(dependent_root, "dependent_alias_index.csv"),
    infeasible = file.path(dependent_root, "dependent_infeasible_records.csv")
  )
  checksums <- data.frame(
    input = names(input_paths), path = substring(normalizePath(unname(input_paths), winslash = "/", mustWork = TRUE),
                     nchar(normalizePath(project_root, winslash = "/", mustWork = TRUE)) + 2L), sha256 = vapply(input_paths, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(checksums, file.path(output_dir, "final_table_input_checksums.csv"), row.names = FALSE)
  report <- c(
    "# Fixed-study final sensitivity tables", "",
    paste0("Core rows: ", nrow(core_table), "; WT-split comparisons: ", nrow(wt_split),
           "; LD bounds: ", nrow(ld), "; supported sparse-data rows: ", nrow(sparse),
           "; grouped source exclusions: ", nrow(loo), "."),
    "The disconnected N<=5 sparse scenario remains an infeasible record and is not represented as a fitted result.", "",
    "## WT-split", "",
    sensitivity_table_md(wt_split, c("scenario", "comparison", "OR", "lower_95_CrI", "upper_95_CrI")), "",
    "## LD bounds", "",
    sensitivity_table_md(ld, c("scenario", "OR", "lower_95_CrI", "upper_95_CrI")), "",
    "## Supported sparse-data result", "",
    sensitivity_table_md(sparse, c("scenario", "OR", "lower_95_CrI", "upper_95_CrI")), "",
    "## Grouped source-influence envelope", "",
    paste0("Median OR range: ", sprintf("%.4f", min(loo$OR)), "–", sprintf("%.4f", max(loo$OR)),
           "; minimum lower 95% CrI: ", sprintf("%.4f", min(loo$lower_95_CrI)), "."), "",
    "All rows are linked to a verified new attempt or an explicitly named completed-core alias."
  )
  writeLines(report, file.path(output_dir, "SENSITIVITY_TABLES.md"))
  invisible(list(
    all = all_sensitivities, core = core_table, dependent = dependent_table,
    wt_split = wt_split, ld = ld, sparse = sparse, loo = loo,
    infeasible = infeasible, checksums = checksums
  ))
}

if (sys.nframe() == 0L) {
  if (file.exists("/analysis/R")) setwd("/analysis")
  source("R/checksums.R")
  result <- build_sensitivity_tables()
  cat("Built fixed-study final tables with", nrow(result$all), "reported comparison rows.\n")
}
