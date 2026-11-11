source("R/common/data_utils.R")
source("R/arm_rows.R")
source("R/scenario_helpers.R")
source("R/checksums.R")
source("R/core_scenarios.R")
source("R/dependent_scenarios.R")

primary <- dependent_primary()
parent <- dependent_parent_metadata()
stopifnot(
  identical(sha256_file(parent$primary_raw_path), parent$primary_raw_sha256),
  identical(as.numeric(dependent_totals(primary)), c(25, 79, 12486, 3821))
)

unresolved <- tryCatch({ build_dependent_scenarios(primary); NULL }, error = identity)
stopifnot(inherits(unresolved, "error"), grepl("explicit LD second-target policy", conditionMessage(unresolved)))

scenarios <- build_dependent_scenarios(primary, "wigle_prospective")
spec <- scenarios$spec
stopifnot(
  nrow(spec) == 31L,
  sum(spec$run_status == "planned") == 27L,
  sum(spec$run_status == "alias") == 3L,
  sum(spec$run_status == "infeasible") == 1L
)

loose <- scenarios$raw[["genotyping_wt_split_loose"]]
strict <- scenarios$raw[["genotyping_wt_split_strict"]]
strict_rows <- strict[strict$Study %in% c("Froehlich_2015", "Jennings_2013") & grepl("^WT_", strict$T), , drop = FALSE]
stopifnot(
  sum(loose$T == "WT_Clean") == 7L,
  sum(strict$T == "WT_Clean") == 5L,
  identical(sort(strict_rows$T), c("WT_Biased", "WT_Biased")),
  identical(sort(strict_rows$N), c(238L, 469L)),
  identical(sort(strict_rows$R), c(37L, 61L)),
  nrow(scenario_row_delta(primary, strict, "strict")) == 4L
)

sparse2 <- scenarios$raw[["sparse_exclude_n2"]]
stopifnot(
  all(sparse2$N > 2L),
  !"sparse_exclude_n5_missing_13" %in% names(scenarios$raw),
  spec$run_status[spec$scenario == "sparse_exclude_n5_missing_13"] == "infeasible",
  grepl("\\*13", spec$infeasible_reason[spec$scenario == "sparse_exclude_n5_missing_13"])
)

loo <- scenarios$loo_map
stopifnot(
  length(unique(loo$loo_group)) == 23L,
  setequal(loo$Study[loo$loo_group == "Lee_N0147"], c("Lee_2016", "Lee_N0147_RareVariants")),
  setequal(loo$Study[loo$loo_group == "Rosmarin_2015"], c("Rosmarin_2015_arm_A", "Rosmarin_2015_arm_B")),
  sum(spec$analysis_family == "loo" & spec$run_status == "planned") == 20L,
  identical(spec$alias_of[spec$scenario == "loo_exclude_Cevik_2022"], "sensitivity_exclude_cevik"),
  identical(spec$alias_of[spec$scenario == "loo_exclude_Morel_2006"], "sensitivity_exclude_pk_informed_morel"),
  identical(spec$alias_of[spec$scenario == "loo_exclude_Onesti_2017"], "sensitivity_exclude_onesti")
)

bozina_best <- scenarios$raw[["ld_best_bozina"]]
bozina_worst <- scenarios$raw[["ld_worst_bozina"]]
row_value <- function(data, study, treatment) {
  row <- data[data$Study == study & data$T == treatment, c("N", "R"), drop = FALSE]
  stopifnot(nrow(row) == 1L)
  as.integer(row[1, ])
}
stopifnot(
  identical(row_value(bozina_best, "Bozina_2021", "HapB3_1129or1236"), c(21L, 10L)),
  identical(row_value(bozina_worst, "Bozina_2021", "HapB3_1129or1236"), c(21L, 9L)),
  all(c("ld_best_wigle_prospective", "ld_worst_wigle_prospective") %in% names(scenarios$raw))
)

both <- build_dependent_scenarios(primary, "both")
stopifnot(
  sum(both$spec$run_status == "planned") == 29L,
  all(c("ld_best_rosmarin_arm_b", "ld_worst_rosmarin_arm_b") %in% names(both$raw))
)

cat("All fixed-study dependent scenario tests passed.\n")
