# =============================================================================
# Fixed-study dependent scenarios derived from the reference primary input only
# =============================================================================

dependent_output_root <- function(project_root = getwd()) {
  file.path(project_root, "results", "analysis", "dependent")
}

dependent_parent_root <- function(project_root = getwd()) {
  file.path(project_root, "results", "analysis", "core")
}

dependent_treatment_map <- function() core_treatment_map()
dependent_treatment_names <- function() core_treatment_names()

dependent_split_wt_map <- function() {
  c(
    "WT_Clean" = 1L,
    "WT_Biased" = 2L,
    "HapB3_1129or1236" = 3L,
    "2846hetho" = 4L,
    "2Ahetho" = 5L,
    "13hetho" = 6L
  )
}

dependent_split_wt_names <- function() {
  c("WT_Clean", "WT_Biased", "HapB3", "2846hetho", "2Ahetho", "13hetho")
}

dependent_primary <- function(project_root = getwd()) {
  parent_root <- dependent_parent_root(project_root)
  manifest_path <- file.path(parent_root, "scenario_input_manifest.csv")
  if (!file.exists(manifest_path)) {
    stop("Committed fixed-study parent manifest is unavailable: ", manifest_path, call. = FALSE)
  }
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  row <- manifest[manifest$scenario == "primary", , drop = FALSE]
  if (nrow(row) != 1L || !isTRUE(row$validation_passed)) {
    stop("Parent manifest must contain one validated primary row.", call. = FALSE)
  }
  input_path <- publication_input_path(parent_root, row$raw_input_path)
  if (!file.exists(input_path) || !identical(sha256_file(input_path), row$raw_input_sha256)) {
    stop("Committed fixed-study primary raw input does not match its manifest checksum.", call. = FALSE)
  }
  primary <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Study", "T", "N", "R")
  if (!all(required %in% names(primary))) stop("Parent primary input has an invalid schema.", call. = FALSE)
  primary
}

dependent_parent_metadata <- function(project_root = getwd()) {
  parent_root <- dependent_parent_root(project_root)
  manifest_path <- file.path(parent_root, "scenario_input_manifest.csv")
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  row <- manifest[manifest$scenario == "primary", , drop = FALSE]
  if (nrow(row) != 1L) stop("Parent primary row is not uniquely manifested.", call. = FALSE)
  list(
    parent_root = parent_root,
    parent_manifest_path = manifest_path,
    parent_manifest_sha256 = sha256_file(manifest_path),
    primary_raw_path = publication_input_path(parent_root, row$raw_input_path),
    primary_raw_sha256 = as.character(row$raw_input_sha256),
    parent_manifest = manifest,
    parent_primary = row
  )
}

dependent_totals <- function(data) core_totals(data)

dependent_sparse_raw <- function(primary, threshold) {
  filtered <- primary[primary$N > threshold, , drop = FALSE]
  variant_studies <- unique(filtered$Study[!filtered$T %in% c("WT_Clean", "WT_Biased")])
  filtered <- filtered[filtered$Study %in% variant_studies, , drop = FALSE]
  rownames(filtered) <- NULL
  filtered
}

dependent_loo_map <- function(primary) {
  studies <- sort(unique(as.character(primary$Study)), method = "radix")
  groups <- studies
  groups[studies %in% c("Lee_2016", "Lee_N0147_RareVariants")] <- "Lee_N0147"
  groups[studies %in% c("Rosmarin_2015_arm_A", "Rosmarin_2015_arm_B")] <- "Rosmarin_2015"
  aliases <- c(
    Cevik_2022 = "sensitivity_exclude_cevik",
    Morel_2006 = "sensitivity_exclude_pk_informed_morel",
    Onesti_2017 = "sensitivity_exclude_onesti"
  )
  data.frame(
    Study = studies,
    loo_group = groups,
    alias_of = unname(aliases[groups]),
    stringsAsFactors = FALSE
  )
}

dependent_ld_second_target_policies <- function() {
  c("wigle_prospective", "rosmarin_arm_b", "both")
}

dependent_ld_targets <- function(policy = NULL) {
  if (length(policy) != 1L || is.na(policy) || !policy %in% dependent_ld_second_target_policies()) {
    stop(
      "An explicit LD second-target policy is required: ",
      paste(dependent_ld_second_target_policies(), collapse = ", "), ".",
      call. = FALSE
    )
  }
  tied <- switch(
    policy,
    wigle_prospective = "Wigle_2021_prospective",
    rosmarin_arm_b = "Rosmarin_2015_arm_B",
    both = c("Wigle_2021_prospective", "Rosmarin_2015_arm_B")
  )
  data.frame(
    Study = c("Bozina_2021", tied),
    expected_n = c(22L, ifelse(tied == "Wigle_2021_prospective", 18L, 18L)),
    expected_r = c(10L, ifelse(tied == "Wigle_2021_prospective", 3L, 12L)),
    target_label = c("bozina", ifelse(tied == "Wigle_2021_prospective", "wigle_prospective", "rosmarin_arm_b")),
    stringsAsFactors = FALSE
  )
}

dependent_spec_row <- function(
  scenario, analysis_family, definition, run_status = "planned", target_study = NA_character_,
  alias_of = NA_character_, infeasible_reason = NA_character_, treatment_names = NULL
) {
  if (is.null(treatment_names)) treatment_names <- paste(dependent_treatment_names(), collapse = "|")
  data.frame(
    scenario = scenario,
    analysis_family = analysis_family,
    definition = definition,
    run_status = run_status,
    target_study = target_study,
    alias_of = alias_of,
    infeasible_reason = infeasible_reason,
    treatment_names = treatment_names,
    stringsAsFactors = FALSE
  )
}

dependent_transform_unified <- function(raw) {
  nma <- transform_to_nma_format(raw, dependent_treatment_map(), aggregate_wt = TRUE)
  validate_nma_data(nma)
  nma[, c("Study", "s.id", "t.id", "r", "n")]
}

build_dependent_scenarios <- function(primary, ld_second_target_policy = NULL) {
  # Do not provide a default: the second LD target is tied and requires an
  # explicit prospective policy decision before any inputs can be generated.
  ld_targets <- dependent_ld_targets(ld_second_target_policy)
  if (anyDuplicated(primary[c("Study", "T")])) stop("Parent primary has duplicate source arms.", call. = FALSE)
  if (!identical(as.numeric(dependent_totals(primary)), c(25, 79, 12486, 3821))) {
    stop("Parent primary totals do not match the committed fixed-study manifest.", call. = FALSE)
  }

  raw <- list()
  nma <- list()
  spec_rows <- list()
  add_planned <- function(scenario, family, definition, data, target = NA_character_, split_wt = FALSE) {
    transformed <- if (split_wt) {
      out <- transform_to_nma_format(data, dependent_split_wt_map(), aggregate_wt = FALSE)
      validate_nma_data(out)
      out[, c("Study", "s.id", "t.id", "r", "n")]
    } else dependent_transform_unified(data)
    spec_rows[[length(spec_rows) + 1L]] <<- dependent_spec_row(
      scenario, family, definition, target_study = target,
      treatment_names = if (split_wt) paste(dependent_split_wt_names(), collapse = "|") else NULL
    )
    raw[[scenario]] <<- data
    nma[[scenario]] <<- transformed
  }

  loose <- primary
  strict <- primary
  strict_rows <- strict$Study %in% c("Froehlich_2015", "Jennings_2013") & strict$T == "WT_Clean"
  if (sum(strict_rows) != 2L) {
    stop("Strict WT split requires exactly the Froehlich 469/61 and Jennings 238/37 rows.", call. = FALSE)
  }
  strict_source <- strict[strict_rows, c("Study", "N", "R"), drop = FALSE]
  if (!identical(sort(strict_source$N), c(238L, 469L)) || !identical(sort(strict_source$R), c(37L, 61L))) {
    stop("Strict WT split source values are not the locked 469/61 and 238/37 rows.", call. = FALSE)
  }
  strict$T[strict_rows] <- "WT_Biased"
  add_planned("genotyping_wt_split_loose", "genotyping", "Separate committed WT_Clean and WT_Biased source labels", loose, split_wt = TRUE)
  add_planned("genotyping_wt_split_strict", "genotyping", "Reclassify Froehlich 469/61 and Jennings 238/37 WT_Clean rows as WT_Biased", strict, split_wt = TRUE)

  sparse2 <- dependent_sparse_raw(primary, 2L)
  add_planned("sparse_exclude_n2", "sparse", "Exclude arms with N <= 2 and then drop WT-only strata", sparse2)
  sparse5 <- dependent_sparse_raw(primary, 5L)
  if ("13hetho" %in% sparse5$T) stop("N <= 5 sparse input unexpectedly retains *13.", call. = FALSE)
  spec_rows[[length(spec_rows) + 1L]] <- dependent_spec_row(
    "sparse_exclude_n5_missing_13", "sparse",
    "Exclude arms with N <= 5; *13 is absent and the fixed-study network is infeasible",
    run_status = "infeasible", infeasible_reason = "The *13 treatment node is absent after N <= 5 filtering"
  )

  for (index in seq_len(nrow(ld_targets))) {
    target <- ld_targets$Study[index]
    row_index <- which(primary$Study == target & primary$T == "HapB3_1129or1236")
    if (length(row_index) != 1L || primary$N[row_index] != ld_targets$expected_n[index] || primary$R[row_index] != ld_targets$expected_r[index]) {
      stop("LD target does not match the committed primary HapB3 arm: ", target, call. = FALSE)
    }
    for (allocation in c("best", "worst")) {
      adjusted <- primary
      adjusted$N[row_index] <- adjusted$N[row_index] - 1L
      if (allocation == "worst") adjusted$R[row_index] <- adjusted$R[row_index] - 1L
      add_planned(
        paste("ld", allocation, ld_targets$target_label[index], sep = "_"), "ld",
        paste0("Remove one bounded discordant carrier from ", target, if (allocation == "best") " as a non-event" else " as an event"),
        adjusted, target
      )
    }
  }

  loo_map <- dependent_loo_map(primary)
  for (group in sort(unique(loo_map$loo_group), method = "radix")) {
    members <- loo_map$Study[loo_map$loo_group == group]
    aliases <- unique(na.omit(loo_map$alias_of[loo_map$loo_group == group]))
    scenario <- paste0("loo_exclude_", group)
    if (length(aliases)) {
      if (length(aliases) != 1L) stop("LOO alias group is ambiguous: ", group, call. = FALSE)
      spec_rows[[length(spec_rows) + 1L]] <- dependent_spec_row(
        scenario, "loo", paste0("Leave out source group ", group, "; completed core alias"),
        run_status = "alias", target_study = group, alias_of = aliases
      )
    } else {
      adjusted <- primary[!primary$Study %in% members, , drop = FALSE]
      rownames(adjusted) <- NULL
      add_planned(scenario, "loo", paste0("Leave out source group ", group, ": ", paste(members, collapse = "; ")), adjusted, group)
    }
  }

  spec <- do.call(rbind, spec_rows)
  rownames(spec) <- NULL
  loo_spec <- spec[spec$analysis_family == "loo", , drop = FALSE]
  if (nrow(loo_spec) != 23L || sum(loo_spec$run_status == "planned") != 20L || sum(loo_spec$run_status == "alias") != 3L) {
    stop("LOO inventory must be 23 groups: 20 new fits and three completed-core aliases.", call. = FALSE)
  }
  if (anyDuplicated(spec$scenario) || !identical(sort(names(raw)), sort(spec$scenario[spec$run_status == "planned"]))) {
    stop("Dependent scenario specification and generated input lists differ.", call. = FALSE)
  }
  list(spec = spec, raw = raw, nma = nma, loo_map = loo_map, ld_targets = ld_targets)
}
