#!/usr/bin/env Rscript

# Render the complete manuscript-facing sensitivity set after the verified
# dependent results exist. The four panels keep scientifically different
# families separate rather than placing unlike WT references on one axis.

dependent_figure_read <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

dependent_figure_required <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
}

core_sensitivity_labels <- c(
  primary = "Primary",
  sensitivity_exclude_cevik = "Exclude Cevik",
  sensitivity_include_iachetta_kleibl = "Include Iachetta and Kleibl",
  sensitivity_include_boisdron_overlap = "Include Boisdron-Celle overlap",
  sensitivity_exclude_onesti = "Exclude Onesti",
  sensitivity_exclude_genotype_informed_carriers = "Exclude genotype-informed carriers",
  sensitivity_exclude_pk_informed_morel = "Exclude PK-informed Morel",
  sensitivity_exclude_all_treatment_moderated = "Exclude all treatment-moderated",
  sensitivity_wigle_early = "Wigle early window",
  sensitivity_cremolini_first_four = "Cremolini first four cycles",
  sensitivity_compound_lower_function = "Compound lower-function allocation",
  sensitivity_compound_reciprocal = "Compound reciprocal allocation"
)

build_dependent_figure_data <- function(
  core_root,
  dependent_root
) {
  core <- dependent_figure_read(
    file.path(core_root, "hapb3_sensitivity_summary.csv"), "core sensitivity summary"
  )
  core_index <- dependent_figure_read(
    file.path(core_root, "core_run_index.csv"), "core run index"
  )
  key <- dependent_figure_read(
    file.path(dependent_root, "dependent_key_comparison_summary.csv"), "dependent key comparisons"
  )
  manifest <- dependent_figure_read(
    file.path(dependent_root, "scenario_input_manifest.csv"), "dependent manifest"
  )
  dependent_figure_required(
    core, c("scenario", "definition", "or_median", "or_lower_95_cri", "or_upper_95_cri"), "Core sensitivity summary"
  )
  dependent_figure_required(
    core_index, c("scenario", "status", "converged", "all_diagnostics_passed", "artifact_checksums_match"), "Core run index"
  )
  dependent_figure_required(
    key, c("scenario", "analysis_family", "numerator", "denominator", "comparison", "or_median", "or_lower_95_cri", "or_upper_95_cri"),
    "Dependent key comparisons"
  )
  dependent_figure_required(
    manifest, c("scenario", "analysis_family", "run_status", "definition", "target_study"), "Dependent manifest"
  )
  if (nrow(core) != 12L || nrow(core_index) != 12L ||
      !all(core_index$status == "complete" & core_index$converged &
           core_index$all_diagnostics_passed & core_index$artifact_checksums_match) ||
      !setequal(core$scenario, names(core_sensitivity_labels))) {
    stop("Core sensitivity figure inputs do not represent twelve verified scenarios.", call. = FALSE)
  }
  core$analysis_family <- "core"
  core$comparison <- "HapB3 vs WT"
  core$label <- unname(core_sensitivity_labels[core$scenario])
  core <- core[, c("scenario", "analysis_family", "comparison", "label", "or_median", "or_lower_95_cri", "or_upper_95_cri")]

  primary <- key[key$scenario == "primary" & key$analysis_family == "primary", , drop = FALSE]
  if (nrow(primary) != 1L || primary$comparison != "HapB3 vs WT") {
    stop("Dependent key summary must contain one verified primary comparison.", call. = FALSE)
  }
  dependent <- key[key$analysis_family != "primary", , drop = FALSE]
  mapped <- manifest[match(dependent$scenario, manifest$scenario), , drop = FALSE]
  if (anyNA(mapped$scenario) || any(mapped$run_status == "infeasible")) {
    stop("Dependent key comparisons do not map exactly to runnable or aliased manifest rows.", call. = FALSE)
  }
  dependent$definition <- mapped$definition
  dependent$target_study <- mapped$target_study
  clean_target <- function(value) gsub("_", " ", value, fixed = TRUE)
  dependent$label <- vapply(seq_len(nrow(dependent)), function(index) {
    scenario <- dependent$scenario[[index]]
    family <- dependent$analysis_family[[index]]
    target <- dependent$target_study[[index]]
    if (family == "genotyping") {
      classification <- if (grepl("strict", scenario, fixed = TRUE)) "Strict" else "Loose"
      paste0(classification, ": ", dependent$comparison[[index]])
    } else if (family == "ld") {
      allocation <- if (grepl("ld_best_", scenario, fixed = TRUE)) "best case" else "worst case"
      paste0(clean_target(target), ": ", allocation)
    } else if (family == "sparse") {
      "Exclude carrier arms with N <= 2"
    } else if (family == "loo") {
      paste0("Exclude ", clean_target(target))
    } else {
      stop("Unsupported dependent figure family: ", family, call. = FALSE)
    }
  }, character(1))
  dependent <- dependent[, c("scenario", "analysis_family", "comparison", "label", "or_median", "or_lower_95_cri", "or_upper_95_cri")]
  counts <- table(dependent$analysis_family)
  expected <- counts[c("genotyping", "sparse", "ld", "loo")]
  if (anyNA(expected) || expected[[1L]] != 4L || expected[[2L]] != 1L ||
      !expected[[3L]] %in% c(4L, 6L) || expected[[4L]] != 23L) {
    stop("Dependent figure families do not match the bounded manifest inventory.", call. = FALSE)
  }
  data <- rbind(core, dependent)
  if (any(!is.finite(data$or_median) | !is.finite(data$or_lower_95_cri) |
          !is.finite(data$or_upper_95_cri) | data$or_lower_95_cri <= 0 |
          data$or_lower_95_cri > data$or_median | data$or_median > data$or_upper_95_cri)) {
    stop("Sensitivity figure contains invalid odds-ratio intervals.", call. = FALSE)
  }
  list(data = data, primary_or = primary$or_median[[1L]], manifest = manifest)
}

plot_core_sensitivity_panel <- function(rows, title, primary_or, point_color = "#1976D2") {
  rows <- rows[rev(seq_len(nrow(rows))), , drop = FALSE]
  y <- seq_len(nrow(rows))
  limits <- range(c(rows$or_lower_95_cri, rows$or_upper_95_cri, 1, primary_or), finite = TRUE)
  limits <- exp(log(limits) + c(-0.08, 0.08) * diff(log(limits)))
  old <- par(mar = c(4.6, max(9, min(21, max(nchar(rows$label)) * 0.42)), 3.2, 2.2))
  on.exit(par(old), add = TRUE)
  plot(
    rows$or_median, y, log = "x", xlim = limits, ylim = c(0.5, nrow(rows) + 0.5),
    yaxt = "n", ylab = "", xlab = "Odds ratio", pch = 19, col = point_color,
    main = title, bty = "n"
  )
  grid(nx = NA, ny = NULL, col = "#ECEFF1", lty = 1)
  segments(rows$or_lower_95_cri, y, rows$or_upper_95_cri, y, col = "#455A64", lwd = 1.6)
  points(rows$or_median, y, pch = 19, col = point_color)
  axis(2, at = y, labels = rows$label, las = 1, tick = FALSE, cex.axis = 0.82)
  abline(v = 1, lty = 2, col = "#6D6D6D")
  abline(v = primary_or, lty = 3, col = "#1976D2")
  box(bty = "l", col = "#607D8B")
}

render_sensitivity_figures <- function(
  project_root = getwd(),
  core_root = file.path(project_root, "results", "analysis", "core"),
  dependent_root = file.path(project_root, "results", "analysis", "dependent"),
  output_dir = file.path(project_root, "runs", "rebuilt", "figures", "sensitivity")
) {
  built <- build_dependent_figure_data(core_root, dependent_root)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  families <- list(
    core = list(title = "Primary and core sensitivity analyses", file = "core"),
    genotyping = list(title = "WT comparator-classification sensitivities", file = "wt_split"),
    ld_sparse = list(title = "Tag/causal-SNP and sparse-data sensitivities", file = "ld_sparse"),
    loo = list(title = "Grouped leave-one-source-out sensitivities", file = "grouped_source_influence")
  )
  panel_rows <- list(
    core = built$data[built$data$analysis_family == "core", , drop = FALSE],
    genotyping = built$data[built$data$analysis_family == "genotyping", , drop = FALSE],
    ld_sparse = built$data[built$data$analysis_family %in% c("ld", "sparse"), , drop = FALSE],
    loo = built$data[built$data$analysis_family == "loo", , drop = FALSE]
  )
  pdf_path <- file.path(output_dir, "all_sensitivity_panels.pdf")
  pdf(pdf_path, width = 12, height = 8.5, onefile = TRUE)
  for (name in names(families)) {
    plot_core_sensitivity_panel(panel_rows[[name]], families[[name]]$title, built$primary_or)
  }
  dev.off()
  png_paths <- character()
  for (name in names(families)) {
    path <- file.path(output_dir, paste0("", families[[name]]$file, "_sensitivity.png"))
    height <- max(1200, 280 + nrow(panel_rows[[name]]) * 62)
    png(path, width = 2200, height = height, res = 180)
    plot_core_sensitivity_panel(panel_rows[[name]], families[[name]]$title, built$primary_or)
    dev.off()
    png_paths <- c(png_paths, path)
  }
  data_path <- file.path(output_dir, "all_sensitivity_figure_data.csv")
  write.csv(built$data, data_path, row.names = FALSE)
  inputs <- c(
    core_results = file.path(core_root, "hapb3_sensitivity_summary.csv"),
    core_index = file.path(core_root, "core_run_index.csv"),
    dependent_results = file.path(dependent_root, "dependent_key_comparison_summary.csv"),
    dependent_manifest = file.path(dependent_root, "scenario_input_manifest.csv")
  )
  checksums <- data.frame(
    input = names(inputs), path = substring(normalizePath(unname(inputs), winslash = "/", mustWork = TRUE),
                     nchar(normalizePath(project_root, winslash = "/", mustWork = TRUE)) + 2L), sha256 = vapply(inputs, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(checksums, file.path(output_dir, "dependent_figure_input_checksums.csv"), row.names = FALSE)
  writeLines(c(
    "# Fixed-study complete sensitivity figures", "",
    "The multi-page PDF and four PNG panels separate core, WT-split, LD/sparse, and grouped-source-influence analyses.",
    "Dashed vertical line: OR=1. Dotted vertical line: fixed-study primary HapB3-versus-WT median.",
    "WT-split panels use WT_Clean or WT_Biased as labelled and must not be interpreted as the same comparison as the unified-WT panels.",
    "All plotted values come from checksum-verified core and dependent result summaries."
  ), file.path(output_dir, "README.md"))
  invisible(list(pdf = pdf_path, png = png_paths, data = built$data, checksums = checksums))
}

if (sys.nframe() == 0L) {
  if (file.exists("/analysis/R")) setwd("/analysis")
  source("R/checksums.R")
  result <- render_sensitivity_figures()
  cat("Rendered complete fixed-study sensitivity figures with", nrow(result$data), "rows.\n")
}
