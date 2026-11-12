#!/usr/bin/env Rscript

# Render manuscript-facing figures from the reference fixed-study core reports.
# This renderer intentionally does not load pcnetmeta output, rankings, or
# absolute-risk calculations.  Every displayed estimate and edge is taken from
# a committed report or the committed primary network-geometry export.

core_figure_root <- function(project_root = getwd()) {
  file.path(project_root, "results", "analysis", "core")
}

core_figure_input_paths <- function(project_root = getwd()) {
  root <- core_figure_root(project_root)
  c(
    primary_posterior = file.path(root, "all_treatment_posterior_summary.csv"),
    hapb3_sensitivity = file.path(root, "hapb3_sensitivity_summary.csv"),
    core_run_index = file.path(root, "core_run_index.csv"),
    direct_comparisons = file.path(root, "primary_direct_comparisons.csv")
  )
}

core_figure_output_dir <- function(project_root = getwd()) {
  file.path(project_root, "runs", "rebuilt", "figures", "core")
}

core_figure_required_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
}

core_figure_read_csv <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

core_figure_sha256 <- function(path) {
  output <- suppressWarnings(system2("sha256sum", path, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L || length(output) != 1L) {
    stop("Could not calculate SHA-256 for ", path, call. = FALSE)
  }
  sub("[[:space:]].*$", "", output)
}

core_treatment_ids <- c(WT = 1L, HapB3 = 2L, `2846hetho` = 3L, `2Ahetho` = 4L, `13hetho` = 5L)
core_treatment_labels <- c(
  WT = "Wild-type", HapB3 = "HapB3", `2846hetho` = "c.2846A>T", `2Ahetho` = "*2A", `13hetho` = "*13"
)
core_treatment_colors <- c(
  WT = "#2E7D32", HapB3 = "#1976D2", `2846hetho` = "#F57C00", `2Ahetho` = "#C62828", `13hetho` = "#8E24AA"
)

core_scenario_labels <- c(
  primary = "Primary",
  sensitivity_exclude_cevik = "Exclude Cevik",
  sensitivity_include_iachetta_kleibl = "Include Iachetta and Kleibl",
  sensitivity_include_boisdron_overlap = "Include Boisdron-Celle overlap",
  sensitivity_exclude_onesti = "Exclude Onesti",
  sensitivity_exclude_genotype_informed_carriers = "Exclude genotype-informed carriers",
  sensitivity_exclude_pk_informed_morel = "Exclude PK-informed Morel",
  sensitivity_exclude_all_treatment_moderated = "Exclude all treatment-moderated",
  sensitivity_wigle_early = "Wigle: early toxicity window",
  sensitivity_cremolini_first_four = "Cremolini: first four cycles",
  sensitivity_compound_lower_function = "Compound: lower-function allocation",
  sensitivity_compound_reciprocal = "Compound: reciprocal allocation"
)

core_expected_edge_keys <- function() {
  apply(utils::combn(1:5, 2L), 2L, paste, collapse = "-")
}

build_core_figure_data <- function(project_root = getwd()) {
  paths <- core_figure_input_paths(project_root)
  primary <- core_figure_read_csv(paths[["primary_posterior"]], "primary posterior report")
  sensitivity <- core_figure_read_csv(paths[["hapb3_sensitivity"]], "HapB3 sensitivity report")
  run_index <- core_figure_read_csv(paths[["core_run_index"]], "core-run index")
  edges <- core_figure_read_csv(paths[["direct_comparisons"]], "primary direct-comparisons geometry")

  core_figure_required_columns(
    primary, c("scenario", "treatment_id", "treatment", "reference", "or_median", "or_lower_95_cri", "or_upper_95_cri"),
    "Primary posterior report"
  )
  core_figure_required_columns(
    sensitivity, c("scenario", "definition", "treatment", "reference", "or_median", "or_lower_95_cri", "or_upper_95_cri"),
    "HapB3 sensitivity report"
  )
  core_figure_required_columns(edges, c("treatment1", "treatment2", "n_studies", "studies"), "Primary direct-comparisons geometry")
  core_figure_required_columns(run_index, c("scenario", "status", "converged", "all_diagnostics_passed", "artifact_checksums_match"), "Core-run index")

  if (nrow(run_index) != 12L || !setequal(run_index$scenario, names(core_scenario_labels)) ||
      !all(run_index$status == "complete") || !all(run_index$converged) ||
      !all(run_index$all_diagnostics_passed) || !all(run_index$artifact_checksums_match)) {
    stop("The committed fixed-study core-run report does not attest 12 complete, converged, checksummed scenarios.", call. = FALSE)
  }

  primary <- primary[primary$scenario == "primary" & primary$reference == "WT", , drop = FALSE]
  primary <- primary[order(primary$treatment_id), , drop = FALSE]
  if (nrow(primary) != 5L || !identical(as.integer(primary$treatment_id), unname(core_treatment_ids)) ||
      !identical(as.character(primary$treatment), names(core_treatment_ids)) ||
      any(!is.finite(as.matrix(primary[c("or_median", "or_lower_95_cri", "or_upper_95_cri")]))) ||
      any(primary$or_lower_95_cri > primary$or_median | primary$or_median > primary$or_upper_95_cri)) {
    stop("Primary posterior report does not contain one valid WT-referenced row for each of the five treatment nodes.", call. = FALSE)
  }
  primary$display_label <- unname(core_treatment_labels[primary$treatment])
  primary$estimate_label <- ifelse(
    primary$treatment == "WT", "Reference (OR 1.00)",
    sprintf("%.2f (%.2f to %.2f)", primary$or_median, primary$or_lower_95_cri, primary$or_upper_95_cri)
  )
  primary$plot_order <- seq_len(nrow(primary))

  expected_scenarios <- names(core_scenario_labels)
  sensitivity <- sensitivity[match(expected_scenarios, sensitivity$scenario), , drop = FALSE]
  if (nrow(sensitivity) != 12L || anyNA(sensitivity$scenario) ||
      !identical(as.character(sensitivity$scenario), expected_scenarios) ||
      !all(sensitivity$treatment == "HapB3") || !all(sensitivity$reference == "WT") ||
      any(!is.finite(as.matrix(sensitivity[c("or_median", "or_lower_95_cri", "or_upper_95_cri")]))) ||
      any(sensitivity$or_lower_95_cri > sensitivity$or_median | sensitivity$or_median > sensitivity$or_upper_95_cri)) {
    stop("HapB3 sensitivity report does not contain the exact 12 valid HapB3-versus-WT scenarios.", call. = FALSE)
  }
  sensitivity$scenario_label <- unname(core_scenario_labels[sensitivity$scenario])
  sensitivity$estimate_label <- sprintf("%.2f (%.2f to %.2f)", sensitivity$or_median, sensitivity$or_lower_95_cri, sensitivity$or_upper_95_cri)
  sensitivity$plot_order <- seq_len(nrow(sensitivity))

  edges$treatment1 <- as.integer(edges$treatment1)
  edges$treatment2 <- as.integer(edges$treatment2)
  edge_key <- paste(pmin(edges$treatment1, edges$treatment2), pmax(edges$treatment1, edges$treatment2), sep = "-")
  if (nrow(edges) != 10L || anyDuplicated(edge_key) || !identical(sort(edge_key), core_expected_edge_keys()) ||
      any(!is.finite(edges$n_studies)) || any(edges$n_studies < 1L)) {
    stop("Primary network geometry must provide each of the 10 observed five-node direct comparisons exactly once.", call. = FALSE)
  }
  edges$edge_key <- edge_key
  edges$comparison_label <- paste(
    core_treatment_labels[names(core_treatment_ids)[edges$treatment1]],
    "vs",
    core_treatment_labels[names(core_treatment_ids)[edges$treatment2]]
  )

  node_ids <- seq_along(core_treatment_ids)
  angles <- seq(pi / 2, pi / 2 + 2 * pi, length.out = length(node_ids) + 1L)[-length(node_ids) - 1L]
  nodes <- data.frame(
    treatment_id = node_ids,
    treatment = names(core_treatment_ids),
    display_label = unname(core_treatment_labels),
    color = unname(core_treatment_colors),
    x = cos(angles), y = sin(angles),
    direct_evidence_study_count = vapply(node_ids, function(node) sum(edges$n_studies[edges$treatment1 == node | edges$treatment2 == node]), numeric(1L)),
    stringsAsFactors = FALSE
  )
  edge_coordinates <- merge(edges, nodes[, c("treatment_id", "x", "y")], by.x = "treatment1", by.y = "treatment_id", all.x = TRUE)
  names(edge_coordinates)[names(edge_coordinates) %in% c("x", "y")] <- c("x1", "y1")
  edge_coordinates <- merge(edge_coordinates, nodes[, c("treatment_id", "x", "y")], by.x = "treatment2", by.y = "treatment_id", all.x = TRUE)
  names(edge_coordinates)[names(edge_coordinates) %in% c("x", "y")] <- c("x2", "y2")

  list(primary_forest = primary, sensitivity_forest = sensitivity, network_nodes = nodes, network_edges = edge_coordinates, input_paths = paths)
}

core_figure_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(hjust = 0),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 26, 10, 10)
    )
}

core_primary_forest_plot <- function(data) {
  display_levels <- rev(data$display_label)
  data$display_label <- factor(data$display_label, levels = display_levels)
  ggplot2::ggplot(data, ggplot2::aes(x = or_median, y = display_label, color = treatment)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey45", linewidth = 0.45) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = or_lower_95_cri, xmax = or_upper_95_cri), orientation = "y", width = 0.18, linewidth = 0.75) +
    ggplot2::geom_point(size = 3.1) +
    ggplot2::geom_text(ggplot2::aes(x = 28, label = estimate_label), inherit.aes = TRUE, hjust = 0, colour = "black", size = 3.3) +
    ggplot2::scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50), limits = c(0.5, 60)) +
    ggplot2::scale_color_manual(values = core_treatment_colors, guide = "none") +
    ggplot2::labs(
      title = "Primary odds ratios for severe toxicity",
      subtitle = "Fixed-study-effect model; 95% credible intervals versus wild type",
      x = "Odds ratio (log scale)"
    ) +
    core_figure_theme()
}

core_network_plot <- function(nodes, edges) {
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = n_studies),
      colour = "grey45", lineend = "round"
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, size = direct_evidence_study_count, colour = treatment),
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = 1.24 * x, y = 1.24 * y, label = display_label),
      fontface = "bold", size = 4.1
    ) +
    ggplot2::scale_colour_manual(values = core_treatment_colors, guide = "none") +
    ggplot2::scale_size_continuous(range = c(7, 13), guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(0.45, 1.75), name = "Contributing strata") +
    ggplot2::coord_equal(xlim = c(-1.55, 1.55), ylim = c(-1.45, 1.45), clip = "off") +
    ggplot2::labs(
      title = "Observed study-co-occurrence network",
      subtitle = "Edges connect variants observed together within at least one study stratum\nWidth = number of contributing strata"
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(hjust = 0),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(10, 20, 10, 20)
    )
}

core_sensitivity_plot <- function(data) {
  display_levels <- rev(data$scenario_label)
  data$scenario_label <- factor(data$scenario_label, levels = display_levels)
  ggplot2::ggplot(data, ggplot2::aes(x = or_median, y = scenario_label)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey45", linewidth = 0.45) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = or_lower_95_cri, xmax = or_upper_95_cri), orientation = "y", width = 0.16, linewidth = 0.7, colour = "#1976D2") +
    ggplot2::geom_point(size = 2.6, colour = "#1976D2") +
    ggplot2::geom_text(ggplot2::aes(x = 3.6, label = estimate_label), hjust = 0, size = 3.05) +
    ggplot2::scale_x_log10(breaks = c(0.75, 1, 1.5, 2, 3, 4), limits = c(0.75, 5.5)) +
    ggplot2::labs(
      title = "HapB3 sensitivity analyses",
      subtitle = "Fixed-study-effect model; 95% credible intervals versus wild type",
      x = "Odds ratio (log scale)"
    ) +
    core_figure_theme(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8.8))
}

core_write_figure_data <- function(figure_data, path) {
  primary <- transform(figure_data$primary_forest,
    figure = "primary_forest", scenario = "primary", scenario_label = NA_character_,
    treatment_id = as.integer(treatment_id), treatment = as.character(treatment),
    node_x = NA_real_, node_y = NA_real_, direct_evidence_study_count = NA_real_,
    edge_treatment1 = NA_integer_, edge_treatment2 = NA_integer_, edge_n_studies = NA_real_, edge_studies = NA_character_
  )
  sensitivity <- transform(figure_data$sensitivity_forest,
    figure = "hapb3_sensitivity_forest", display_label = scenario_label,
    treatment_id = 2L, node_x = NA_real_, node_y = NA_real_, direct_evidence_study_count = NA_real_,
    edge_treatment1 = NA_integer_, edge_treatment2 = NA_integer_, edge_n_studies = NA_real_, edge_studies = NA_character_
  )
  nodes <- transform(figure_data$network_nodes,
    figure = "network_node", scenario = "primary", scenario_label = NA_character_,
    reference = "WT", or_median = NA_real_, or_lower_95_cri = NA_real_, or_upper_95_cri = NA_real_, estimate_label = NA_character_,
    plot_order = treatment_id, node_x = x, node_y = y,
    edge_treatment1 = NA_integer_, edge_treatment2 = NA_integer_, edge_n_studies = NA_real_, edge_studies = NA_character_
  )
  edges <- transform(figure_data$network_edges,
    figure = "network_edge", scenario = "primary", scenario_label = NA_character_, display_label = comparison_label,
    treatment_id = NA_integer_, treatment = NA_character_, reference = "WT", or_median = NA_real_, or_lower_95_cri = NA_real_, or_upper_95_cri = NA_real_,
    estimate_label = NA_character_, plot_order = NA_integer_, node_x = NA_real_, node_y = NA_real_, direct_evidence_study_count = NA_real_,
    edge_treatment1 = treatment1, edge_treatment2 = treatment2, edge_n_studies = n_studies, edge_studies = studies
  )
  columns <- c("figure", "scenario", "scenario_label", "treatment_id", "treatment", "reference", "display_label", "or_median", "or_lower_95_cri", "or_upper_95_cri", "estimate_label", "plot_order", "node_x", "node_y", "direct_evidence_study_count", "edge_treatment1", "edge_treatment2", "edge_n_studies", "edge_studies")
  output <- do.call(rbind, lapply(list(primary, sensitivity, nodes, edges), function(value) {
    missing <- setdiff(columns, names(value))
    value[missing] <- NA
    value[, columns, drop = FALSE]
  }))
  write.csv(output, path, row.names = FALSE, na = "")
}

render_core_figures <- function(project_root = getwd()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required to render fixed-study figures.", call. = FALSE)
  figure_data <- build_core_figure_data(project_root)
  output_dir <- core_figure_output_dir(project_root)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  plots <- list(
    primary_forest = core_primary_forest_plot(figure_data$primary_forest),
    primary_network = core_network_plot(figure_data$network_nodes, figure_data$network_edges),
    hapb3_sensitivity_forest = core_sensitivity_plot(figure_data$sensitivity_forest)
  )
  dimensions <- list(
    primary_forest = c(10, 4.6),
    primary_network = c(7.2, 6.6),
    hapb3_sensitivity_forest = c(10.5, 7.8)
  )
  for (name in names(plots)) {
    dims <- dimensions[[name]]
    ggplot2::ggsave(file.path(output_dir, paste0(name, ".png")), plots[[name]], width = dims[1L], height = dims[2L], dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(output_dir, paste0(name, ".pdf")), plots[[name]], width = dims[1L], height = dims[2L], device = grDevices::cairo_pdf, bg = "white")
  }

  core_write_figure_data(figure_data, file.path(output_dir, "core_figure_data.csv"))
  input_paths <- figure_data$input_paths
  checksums <- data.frame(
    input = names(input_paths),
    relative_path = sub(paste0("^", normalizePath(project_root, winslash = "/", mustWork = TRUE), "/"), "", normalizePath(input_paths, winslash = "/", mustWork = TRUE)),
    sha256 = vapply(input_paths, core_figure_sha256, character(1L)),
    stringsAsFactors = FALSE
  )
  write.csv(checksums, file.path(output_dir, "input_checksums.csv"), row.names = FALSE)
  writeLines(c(
    "# Fixed-study core figures", "",
    "Rendered from the reference fixed-study core reports and primary network geometry only.",
    "No pcnetmeta refit, ranking, or absolute-risk output is used.", "",
    "- `primary_forest.{png,pdf}` — all five modeled nodes versus wild-type.",
    "- `primary_network.{png,pdf}` — all 10 observed study-co-occurrence edges.",
    "- `hapb3_sensitivity_forest.{png,pdf}` — 12 prespecified HapB3 scenarios.",
    "- `core_figure_data.csv` — machine-readable plotted data, including derived node coordinates and direct edges.",
    "- `input_checksums.csv` — SHA-256 checksums of the four committed input artifacts.", "",
    "Renderer: `R/render_core_figures.R`; execution environment: pinned Docker image `runtime-recorded`."
  ), file.path(output_dir, "README.md"))
  invisible(list(plots = plots, figure_data = figure_data, output_dir = output_dir))
}

if (sys.nframe() == 0L) {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
  rendered <- render_core_figures(project_root)
  cat("Rendered fixed-study core figures to", rendered$output_dir, "\n")
}
