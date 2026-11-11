source("R/network_geometry.R")

expect_error <- function(expression, pattern) {
  error <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = function(e) e
  )
  stopifnot(inherits(error, "error"))
  stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
}

arms <- function(study, treatment) {
  data.frame(Study = study, t.id = as.integer(treatment), stringsAsFactors = FALSE)
}

# Repeated AB evidence remains one unique design and has no inconsistency df.
repeated_ab <- rbind(
  arms(c("AB_1", "AB_1", "AB_2", "AB_2"), c(1, 2, 1, 2))
)
geometry_ab <- analyze_network_geometry(repeated_ab)
stopifnot(
  geometry_ab$connectedness$connected,
  geometry_ab$baseline == 1L,
  nrow(geometry_ab$designs) == 1L,
  geometry_ab$designs$n_studies == 2L,
  geometry_ab$m == 1L,
  geometry_ab$rank_A == 1L,
  geometry_ab$q == 0L,
  geometry_ab$q_from_rank == 0L,
  geometry_ab$arm_level_rank_check$treatment_rank_check$passes,
  geometry_ab$arm_level_rank_check$inconsistency_rank_increase == 0L,
  geometry_ab$arm_level_rank_check$passes,
  geometry_ab$direct_edges$n_studies == 2L,
  identical(geometry_ab$direct_edges$studies[[1L]], c("AB_1", "AB_2"))
)

# Three two-arm designs create one independent inconsistency degree of freedom.
triangle <- rbind(
  arms(c("AB", "AB", "AC", "AC", "BC", "BC"), c(1, 2, 1, 3, 2, 3))
)
geometry_triangle <- analyze_network_geometry(triangle)
stopifnot(
  geometry_triangle$m == 3L,
  geometry_triangle$rank_A == 2L,
  geometry_triangle$q == 1L,
  geometry_triangle$q_from_rank == 1L,
  geometry_triangle$arm_level_rank_check$treatment_rank_increase == 2L,
  geometry_triangle$arm_level_rank_check$treatment_rank_check$passes,
  geometry_triangle$arm_level_rank_check$rank_consistent == 5L,
  geometry_triangle$arm_level_rank_check$rank_saturated == 6L,
  geometry_triangle$arm_level_rank_check$inconsistency_rank_increase == 1L,
  geometry_triangle$arm_level_rank_check$passes,
  nrow(geometry_triangle$direct_edges) == 3L,
  length(geometry_triangle$loops) == 1L,
  identical(geometry_triangle$loops[[1L]], c(1L, 2L, 3L)),
  identical(geometry_triangle$loop_table$cycle, "1-2-3-1")
)

# A single multi-arm design contains graph cycles but no design-by-treatment df.
single_abc <- arms(rep("ABC", 3L), c(1, 2, 3))
geometry_abc <- analyze_network_geometry(single_abc)
stopifnot(
  geometry_abc$m == 2L,
  geometry_abc$rank_A == 2L,
  geometry_abc$q == 0L,
  geometry_abc$arm_level_rank_check$inconsistency_rank_increase == 0L,
  geometry_abc$arm_level_rank_check$passes,
  length(geometry_abc$loops) == 1L
)

# Adding an already-observed design only adds studies, never degrees of freedom.
triangle_with_repeat <- rbind(triangle, arms(c("AB_repeat", "AB_repeat"), c(1, 2)))
geometry_repeat <- analyze_network_geometry(triangle_with_repeat)
stopifnot(
  geometry_repeat$m == geometry_triangle$m,
  geometry_repeat$q == geometry_triangle$q,
  geometry_repeat$arm_level_rank_check$inconsistency_rank_increase == 1L,
  geometry_repeat$arm_level_rank_check$passes,
  geometry_repeat$designs$n_studies[geometry_repeat$designs$design == "1-2"] == 2L
)

# In a connected network without WT node 1, the lowest treatment is the baseline.
no_wt <- rbind(
  arms(c("23", "23", "24", "24"), c(2, 3, 2, 4))
)
geometry_no_wt <- analyze_network_geometry(no_wt)
stopifnot(
  geometry_no_wt$connectedness$connected,
  geometry_no_wt$baseline == 2L,
  all(geometry_no_wt$designs$baseline == 2L),
  geometry_no_wt$rank_A == 2L,
  geometry_no_wt$q == 0L,
  identical(colnames(geometry_no_wt$A), c("treatment_3", "treatment_4"))
)

# Components are inspectable, while the full analysis refuses a disconnected input.
disconnected <- rbind(
  arms(c("AB", "AB", "CD", "CD"), c(1, 2, 3, 4))
)
components <- network_connectedness(disconnected)
stopifnot(!components$connected, components$n_components == 2L)
expect_error(analyze_network_geometry(disconnected), "Network is disconnected")
expect_error(analyze_network_geometry(arms("single", 1)), "single-arm")
expect_error(analyze_network_geometry(arms(c("duplicate", "duplicate"), c(1, 1))), "duplicate study-treatment")
expect_error(analyze_network_geometry(arms(c("missing", "missing"), c(1, NA))), "missing or non-numeric treatment node")

# All public results are invariant to arm and study ordering.
set.seed(20260825)
shuffled_triangle <- triangle[sample(seq_len(nrow(triangle))), , drop = FALSE]
geometry_shuffled <- analyze_network_geometry(shuffled_triangle)
stopifnot(
  identical(geometry_shuffled$designs, geometry_triangle$designs),
  identical(geometry_shuffled$direct_edges, geometry_triangle$direct_edges),
  identical(geometry_shuffled$loops, geometry_triangle$loops),
  identical(geometry_shuffled$A, geometry_triangle$A),
  identical(geometry_shuffled$A_mapping, geometry_triangle$A_mapping),
  identical(geometry_shuffled$arm_level_rank_check, geometry_triangle$arm_level_rank_check),
  identical(geometry_shuffled$q, geometry_triangle$q)
)

cat("All network geometry tests passed.\n")
