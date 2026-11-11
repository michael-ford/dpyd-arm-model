source("R/render_core_figures.R")

project_root <- normalizePath(getwd(), mustWork = TRUE)
figure_data <- build_core_figure_data(project_root)

stopifnot(
  identical(figure_data$primary_forest$treatment, names(core_treatment_ids)),
  identical(figure_data$primary_forest$display_label, unname(core_treatment_labels)),
  nrow(figure_data$primary_forest) == 5L,
  all(figure_data$primary_forest$reference == "WT")
)

edge_keys <- paste(figure_data$network_edges$treatment1, figure_data$network_edges$treatment2, sep = "-")
expected_edge_labels <- c(
  "Wild-type vs HapB3", "Wild-type vs c.2846A>T", "Wild-type vs *2A", "Wild-type vs *13",
  "HapB3 vs c.2846A>T", "HapB3 vs *2A", "HapB3 vs *13", "c.2846A>T vs *2A",
  "c.2846A>T vs *13", "*2A vs *13"
)
stopifnot(
  nrow(figure_data$network_nodes) == 5L,
  identical(figure_data$network_nodes$treatment, names(core_treatment_ids)),
  nrow(figure_data$network_edges) == 10L,
  identical(sort(edge_keys), core_expected_edge_keys()),
  setequal(figure_data$network_edges$comparison_label, expected_edge_labels),
  all(figure_data$network_edges$n_studies >= 1L)
)

stopifnot(
  nrow(figure_data$sensitivity_forest) == 12L,
  identical(figure_data$sensitivity_forest$scenario, names(core_scenario_labels)),
  identical(figure_data$sensitivity_forest$scenario_label, unname(core_scenario_labels)),
  all(figure_data$sensitivity_forest$treatment == "HapB3"),
  all(figure_data$sensitivity_forest$reference == "WT")
)

cat("Fixed-study core figure data checks passed.\n")
