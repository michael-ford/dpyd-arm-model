# =============================================================================
# Fixed-study-effect scenario inputs
# =============================================================================

core_treatment_map <- function() {
  c(
    "WT_Clean" = 1,
    "WT_Biased" = 1,
    "HapB3_1129or1236" = 2,
    "2846hetho" = 3,
    "2Ahetho" = 4,
    "13hetho" = 5
  )
}

core_treatment_names <- function() {
  c("WT", "HapB3", "2846hetho", "2Ahetho", "13hetho")
}

core_scenario_spec <- function() {
  data.frame(
    scenario = c(
      "primary",
      "sensitivity_exclude_cevik",
      "sensitivity_include_iachetta_kleibl",
      "sensitivity_include_boisdron_overlap",
      "sensitivity_exclude_onesti",
      "sensitivity_exclude_genotype_informed_carriers",
      "sensitivity_exclude_pk_informed_morel",
      "sensitivity_exclude_all_treatment_moderated",
      "sensitivity_wigle_early",
      "sensitivity_cremolini_first_four",
      "sensitivity_compound_lower_function",
      "sensitivity_compound_reciprocal"
    ),
    expected_studies = c(25, 24, 27, 26, 24, 25, 24, 24, 25, 25, 25, 25),
    expected_arms = c(79, 77, 85, 82, 77, 71, 75, 67, 79, 79, 79, 79),
    expected_n = c(12486, 12383, 12982, 12658, 12360, 12395, 12000, 11909, 12486, 12486, 12489, 12489),
    expected_r = c(3821, 3795, 4047, 3836, 3793, 3778, 3778, 3735, 3682, 3771, 3824, 3824),
    definition = c(
      "Primary: fixed-study-effect eligible strata and mutually exclusive pure-variant arms",
      "Remove the complete approved Cevik HapB3 comparison",
      "Restore source-reconciled outcome-sampled Iachetta and Kleibl rows",
      "Restore Boisdron-Celle as a pure-arm overlap sensitivity alongside Morel",
      "Remove the smaller probably overlapping Onesti stratum",
      "Remove prospective guided Wigle and guideline-guided Bozina carrier arms only",
      "Remove the complete PK-informed Morel stratum; Boisdron-Celle is already absent",
      "Apply genotype-informed carrier-arm and PK-informed Morel exclusions together",
      "Use the source-supported early toxicity window for prospective Wigle only",
      "Use the source-supported first-four-cycle toxicity window for Cremolini only",
      "Allocate each modeled compound once under the lower-function rule",
      "Allocate each modeled compound once under the reciprocal rule"
    ),
    stringsAsFactors = FALSE
  )
}

core_totals <- function(data) {
  c(
    studies = length(unique(data$Study)),
    arms = nrow(data),
    N = sum(data$N),
    R = sum(data$R)
  )
}

core_node_totals <- function(data) {
  labels <- c("WT_Clean", "WT_Biased", "HapB3_1129or1236", "2846hetho", "2Ahetho", "13hetho")
  node <- c("WT", "WT", "HapB3", "2846", "2A", "13")
  out <- lapply(unique(node), function(name) {
    rows <- data[data$T %in% labels[node == name], , drop = FALSE]
    c(arms = nrow(rows), N = sum(rows$N), R = sum(rows$R))
  })
  names(out) <- unique(node)
  out
}

core_assert_exact_rows <- function(data, study, expected) {
  actual <- normalize_arm_rows(data[data$Study == study, , drop = FALSE])
  expected <- normalize_arm_rows(expected)
  if (!identical(actual, expected)) {
    stop("Study ", study, " does not match its locked fixed-study rows.", call. = FALSE)
  }
  invisible(TRUE)
}

core_replace_arm <- function(data, study, treatment, expected_n, expected_r, n, r) {
  index <- which(data$Study == study & data$T == treatment)
  if (length(index) != 1L || data$N[index] != expected_n || data$R[index] != expected_r) {
    stop("Cannot apply locked fixed-study change to ", study, "/", treatment, ".", call. = FALSE)
  }
  data$N[index] <- n
  data$R[index] <- r
  data
}

core_primary_replacements <- function() {
  list(
    Etienne_Grimaldi_2017 = make_arm_rows(
      "Etienne_Grimaldi_2017",
      c("13hetho", "HapB3_1129or1236", "2846hetho", "2Ahetho", "WT_Biased"),
      c(1, 4, 3, 3, 231), c(1, 0, 2, 2, 25), c(1, 1, 1, 1, 2)
    ),
    Froehlich_2015 = make_arm_rows(
      "Froehlich_2015",
      c("HapB3_1129or1236", "13hetho", "2846hetho", "2Ahetho", "WT_Clean"),
      c(22, 1, 3, 4, 469), c(10, 0, 0, 0, 61)
    ),
    Largillier_2006 = make_arm_rows(
      "Largillier_2006", c("2Ahetho", "WT_Biased"), c(1, 83), c(1, 14), c(1, 2)
    ),
    Loganayagam_2013 = make_arm_rows(
      "Loganayagam_2013",
      c("13hetho", "HapB3_1129or1236", "2846hetho", "2Ahetho", "WT_Clean"),
      c(1, 15, 3, 3, 391), c(1, 6, 3, 3, 74)
    ),
    Toffoli_2015 = make_arm_rows(
      "Toffoli_2015",
      c("13hetho", "2846hetho", "2Ahetho", "WT_Biased"),
      c(1, 5, 11, 585), c(0, 3, 7, 84)
    ),
    Rosmarin_2015_arm_A = make_arm_rows(
      "Rosmarin_2015_arm_A", c("HapB3_1129or1236", "WT_Biased"), c(17, 422), c(8, 126)
    ),
    Rosmarin_2015_arm_B = make_arm_rows(
      "Rosmarin_2015_arm_B", c("HapB3_1129or1236", "WT_Biased"), c(18, 483), c(12, 171)
    )
  )
}

core_boisdron_rows <- function() {
  make_arm_rows(
    "Boisdron_Celle_2007", c("2846hetho", "2Ahetho", "WT_Biased"),
    c(7, 2, 163), c(5, 1, 9), c(1, 1, 2)
  )
}

build_core_primary <- function(raw_data) {
  primary <- normalize_arm_rows(raw_data)
  committed <- normalize_arm_rows(read.csv("data/core/primary_raw.csv", stringsAsFactors = FALSE))
  if (!identical(primary, committed)) stop("Primary rows differ from the committed mutually exclusive arms.", call. = FALSE)
  validate_core_scenario(raw_data, "primary")
  as.data.frame(raw_data, stringsAsFactors = FALSE)
}

core_apply_genotype_informed_carrier_exclusion <- function(data) {
  remove <-
    (data$Study == "Wigle_2021_prospective" & data$T != "WT_Clean") |
    (data$Study == "Bozina_2021" & data$T != "WT_Biased")
  data[!remove, , drop = FALSE]
}

core_apply_wigle_early_endpoint <- function(data) {
  data <- core_replace_arm(data, "Wigle_2021_prospective", "WT_Clean", 1347, 418, 1347, 284)
  data <- core_replace_arm(data, "Wigle_2021_prospective", "HapB3_1129or1236", 18, 3, 18, 1)
  data <- core_replace_arm(data, "Wigle_2021_prospective", "2Ahetho", 9, 3, 9, 2)
  data <- core_replace_arm(data, "Wigle_2021_prospective", "2846hetho", 19, 5, 19, 3)
  core_replace_arm(data, "Wigle_2021_prospective", "13hetho", 1, 0, 1, 0)
}

core_apply_cremolini_first_four_endpoint <- function(data) {
  data <- core_replace_arm(data, "Cremolini_2018", "WT_Biased", 429, 215, 429, 166)
  data <- core_replace_arm(data, "Cremolini_2018", "2Ahetho", 5, 4, 5, 3)
  core_replace_arm(data, "Cremolini_2018", "2846hetho", 5, 4, 5, 4)
}

core_apply_compound_allocation <- function(data, allocation) {
  if (identical(allocation, "lower_function")) {
    data <- core_replace_arm(data, "Froehlich_2015", "13hetho", 1, 0, 2, 1)
    data <- core_replace_arm(data, "Lee_N0147_RareVariants", "2Ahetho", 24, 21, 25, 22)
    return(core_replace_arm(data, "Toffoli_2015", "2Ahetho", 11, 7, 12, 8))
  }
  if (identical(allocation, "reciprocal")) {
    data <- core_replace_arm(data, "Froehlich_2015", "HapB3_1129or1236", 22, 10, 23, 11)
    data <- core_replace_arm(data, "Lee_N0147_RareVariants", "2846hetho", 26, 21, 27, 22)
    return(core_replace_arm(data, "Toffoli_2015", "13hetho", 1, 0, 2, 1))
  }
  stop("Unknown compound allocation: ", allocation, call. = FALSE)
}

build_core_scenario <- function(raw_data, scenario_name) {
  spec <- core_scenario_spec()
  if (!scenario_name %in% spec$scenario) {
    stop("Unknown fixed-study scenario: ", scenario_name, call. = FALSE)
  }

  scenario <- build_core_primary(raw_data)
  if (scenario_name == "sensitivity_exclude_cevik") {
    scenario <- scenario[scenario$Study != "Cevik_2022", , drop = FALSE]
  } else if (scenario_name == "sensitivity_include_iachetta_kleibl") {
    scenario <- bind_replacement(scenario, outcome_sampled_reinclusion_rows())
  } else if (scenario_name == "sensitivity_include_boisdron_overlap") {
    scenario <- bind_replacement(scenario, core_boisdron_rows())
  } else if (scenario_name == "sensitivity_exclude_onesti") {
    scenario <- scenario[scenario$Study != "Onesti_2017", , drop = FALSE]
  } else if (scenario_name == "sensitivity_exclude_genotype_informed_carriers") {
    scenario <- core_apply_genotype_informed_carrier_exclusion(scenario)
  } else if (scenario_name == "sensitivity_exclude_pk_informed_morel") {
    scenario <- scenario[scenario$Study != "Morel_2006", , drop = FALSE]
  } else if (scenario_name == "sensitivity_exclude_all_treatment_moderated") {
    scenario <- core_apply_genotype_informed_carrier_exclusion(scenario)
    scenario <- scenario[scenario$Study != "Morel_2006", , drop = FALSE]
  } else if (scenario_name == "sensitivity_wigle_early") {
    scenario <- core_apply_wigle_early_endpoint(scenario)
  } else if (scenario_name == "sensitivity_cremolini_first_four") {
    scenario <- core_apply_cremolini_first_four_endpoint(scenario)
  } else if (scenario_name == "sensitivity_compound_lower_function") {
    scenario <- core_apply_compound_allocation(scenario, "lower_function")
  } else if (scenario_name == "sensitivity_compound_reciprocal") {
    scenario <- core_apply_compound_allocation(scenario, "reciprocal")
  }

  rownames(scenario) <- NULL
  validate_core_scenario(scenario, scenario_name)
  scenario
}

validate_core_scenario <- function(data, scenario_name = NULL) {
  validate_verified_scenario(data)

  allowed <- names(core_treatment_map())
  unmapped <- setdiff(unique(data$T), allowed)
  if (length(unmapped) > 0L) {
    stop("Scenario has unmapped treatment labels: ", paste(unmapped, collapse = ", "), call. = FALSE)
  }
  if ("Wigle_2021_retrospective" %in% data$Study) {
    stop("Retrospective single-arm Wigle is ineligible under fixed study effects.", call. = FALSE)
  }
  if (all(c("Rosmarin_2015_arm_A", "Rosmarin_2015_arm_B") %in% data$Study)) {
    for (study in c("Rosmarin_2015_arm_A", "Rosmarin_2015_arm_B")) {
      core_assert_exact_rows(data, study, core_primary_replacements()[[study]])
    }
  }
  if ("Morel_2006" %in% data$Study) {
    core_assert_exact_rows(data, "Morel_2006", verified_replacements()$Morel_2006)
  }
  if ("Boisdron_Celle_2007" %in% data$Study) {
    if (!identical(scenario_name, "sensitivity_include_boisdron_overlap")) {
      stop("Boisdron-Celle is present outside its approved overlap sensitivity.", call. = FALSE)
    }
    core_assert_exact_rows(data, "Boisdron_Celle_2007", core_boisdron_rows())
  }
  if (any(c("Iachetta_2019", "Kleibl_2009") %in% data$Study) &&
      !identical(scenario_name, "sensitivity_include_iachetta_kleibl")) {
    stop("Outcome-sampled cohorts are present outside their approved inclusive sensitivity.", call. = FALSE)
  }

  spec <- core_scenario_spec()
  if (!is.null(scenario_name) && scenario_name %in% spec$scenario) {
    expected <- spec[spec$scenario == scenario_name, , drop = FALSE]
    actual <- as.numeric(core_totals(data))
    target <- as.numeric(expected[1, c("expected_studies", "expected_arms", "expected_n", "expected_r")])
    if (!identical(actual, target)) {
      stop(
        "Scenario ", scenario_name, " totals do not match the locked manifest. Expected ",
        paste(target, collapse = "/"), "; found ", paste(actual, collapse = "/"), ".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}
