source("R/fixed_study_model.R")

production <- fixed_study_production_config()
stopifnot(identical(production, list(
  n.adapt = 50000L,
  n.burnin = 125000L,
  n.iter = 250000L,
  n.thin = 10L,
  n.chains = 4L
)))
stopifnot(identical(fixed_study_mcmc_config(), production))
mutated <- production
mutated$n.iter <- 1L
stopifnot(fixed_study_production_config()$n.iter == 250000L)
locked_error <- tryCatch({ validate_fixed_study_config(mutated); NULL }, error = identity)
stopifnot(inherits(locked_error, "error"), grepl("test_only", conditionMessage(locked_error)))

runs_a <- fixed_study_seed_schedule(c("zeta", "alpha", "middle"))
runs_b <- fixed_study_seed_schedule(c("middle", "zeta", "alpha"))
stopifnot(identical(runs_a, runs_b))
stopifnot(nrow(runs_a) == 12L, !anyDuplicated(runs_a$seed), all(runs_a$seed > 0))
stopifnot(identical(unique(runs_a$scenario_name), sort(unique(runs_a$scenario_name))))
stopifnot(identical(
  fixed_study_chain_seeds("middle", runs_a),
  runs_a$seed[runs_a$scenario_name == "middle"]
))

raw <- data.frame(
  Study = c("S1", "S1", "S2", "S2", "S3"),
  T = c("WT", "HapB3", "HapB3", "Rare", "WT"),
  N = c(20L, 10L, 12L, 8L, 5L),
  R = c(2L, 3L, 4L, 2L, 0L),
  stringsAsFactors = FALSE
)
map <- c(WT = 1L, HapB3 = 2L, Rare = 3L)
trt_names <- c("WT", "HapB3", "Rare")
prepared <- prepare_fixed_study_data(raw, map, trt_names, single_arm = "drop")
stopifnot(
  identical(prepared$dropped_single_arm_studies, "S3"),
  prepared$jags_data$Narm == 4L,
  prepared$jags_data$Nstudy == 2L,
  prepared$network_connected,
  all(c("HapB3", "Rare") %in% prepared$arms$T[prepared$arms$Study == "S2"]),
  !any(prepared$arms$Study == "S3")
)
reject_error <- tryCatch({
  prepare_fixed_study_data(raw, map, trt_names, single_arm = "reject")
  NULL
}, error = identity)
stopifnot(inherits(reject_error, "error"), grepl("Single-arm", conditionMessage(reject_error)))

duplicate_alias <- rbind(raw[1:2, ], data.frame(Study = "S1", T = "WT alias", N = 10L, R = 1L))
duplicate_map_error <- tryCatch({
  prepare_fixed_study_data(duplicate_alias, c(map, "WT alias" = 1L), trt_names)
  NULL
}, error = identity)
stopifnot(inherits(duplicate_map_error, "error"), grepl("duplicate Study/treatment", conditionMessage(duplicate_map_error)))

missing_node_error <- tryCatch({
  prepare_fixed_study_data(raw[1:2, ], map, trt_names, single_arm = "reject")
  NULL
}, error = identity)
stopifnot(
  inherits(missing_node_error, "error"),
  grepl("Treatment node.*absent", conditionMessage(missing_node_error)),
  grepl("Rare", conditionMessage(missing_node_error))
)

disconnected_raw <- data.frame(
  Study = c("Component A", "Component A", "Component B", "Component B"),
  T = c("A", "B", "C", "D"),
  N = rep(10L, 4L),
  R = c(1L, 2L, 3L, 4L),
  stringsAsFactors = FALSE
)
disconnected_error <- tryCatch({
  prepare_fixed_study_data(
    disconnected_raw,
    treatment_map = c(A = 1L, B = 2L, C = 3L, D = 4L),
    treatment_names = c("A", "B", "C", "D"),
    single_arm = "reject"
  )
  NULL
}, error = identity)
stopifnot(
  inherits(disconnected_error, "error"),
  grepl("network is disconnected", conditionMessage(disconnected_error))
)

# Saturated-binomial residual deviance is finite and exactly zero at matching
# zero-event/all-event probability boundaries.
boundary_p <- rbind(c(0, 1, 0.5), c(0.1, 0.9, 0.25))
boundary_dev <- fixed_study_residual_deviance(
  boundary_p, events = c(0, 10, 2), participants = c(10, 10, 8)
)
stopifnot(
  all(is.finite(boundary_dev$draws_by_arm)),
  boundary_dev$draws_by_arm[1, 1] < 1e-12,
  boundary_dev$draws_by_arm[1, 2] < 1e-12,
  identical(boundary_dev$total_summary$per_arm_mean,
            boundary_dev$total_summary$mean / 3)
)

if (!requireNamespace("coda", quietly = TRUE)) {
  stop("Pinned test environment must provide coda.")
}
set.seed(80625)
synthetic_chains <- lapply(seq_len(4L), function(chain) {
  n <- 2500L
  d2 <- stats::rnorm(n, 0.4, 0.2)
  d3 <- stats::rnorm(n, 0.8, 0.25)
  tau <- abs(stats::rnorm(n, 0.15, 0.03))
  coda::mcmc(cbind(
    "d[1]" = 0,
    "d[2]" = d2,
    "d[3]" = d3,
    "p[1]" = stats::plogis(stats::rnorm(n, -2, 0.2)),
    "p[2]" = stats::plogis(stats::rnorm(n, -1.5, 0.2)),
    tau_arm = tau,
    tau_contrast = sqrt(2) * tau
  ))
})
synthetic <- do.call(coda::mcmc.list, synthetic_chains)
augmented <- augment_fixed_study_samples(synthetic, 3L)
stopifnot(all(c(
  "OR[1]", "OR[2]", "OR[3]",
  "log_contrast[2,1]", "log_contrast[3,1]", "log_contrast[3,2]",
  "pairwise_OR[3,2]"
) %in% colnames(augmented[[1]])))
stopifnot(all(as.matrix(augmented[[1]])[, "OR[1]"] == 1))

summary <- summarize_fixed_study_posterior(augmented, trt_names)
stopifnot(
  nrow(summary$treatment_effects) == 3L,
  nrow(summary$pairwise_effects) == choose(3L, 2L),
  nrow(summary$heterogeneity) == 2L,
  summary$treatment_effects$or_median[1] == 1
)
diagnostics <- fixed_study_diagnostics(augmented, 3L)
stopifnot(
  all(c("psrf_point", "psrf_upper_95", "effective_size", "passed") %in% names(diagnostics$table)),
  all(is.finite(diagnostics$table$psrf_point)),
  all(diagnostics$table$parameter_type %in% c("treatment_effect", "heterogeneity")),
  diagnostics$converged,
  diagnostics$max_psrf_upper_95 <= 1.05,
  diagnostics$min_effective_size >= 1000
)

preflight <- fit_fixed_study_model(
  scenario_name = "mapping_preflight",
  raw_rows = raw,
  treatment_map = map,
  treatment_names = trt_names,
  chain_seeds = 101:104,
  single_arm = "drop",
  preflight = TRUE
)
stopifnot(
  preflight$preflight,
  identical(preflight$config, fixed_study_production_config()),
  identical(preflight$chain_seeds, 101:104),
  preflight$prepared_data$jags_data$Narm == 4L
)

# A deliberately tiny test-only fit verifies the JAGS data/inits/model wiring.
# It is not an inferential run and is never accepted as production config.
if (requireNamespace("rjags", quietly = TRUE)) {
  tiny <- list(n.adapt = 100L, n.burnin = 100L, n.iter = 200L, n.thin = 1L, n.chains = 4L)
  smoke <- suppressWarnings(fit_fixed_study_model(
    scenario_name = "tiny_jags_smoke",
    raw_rows = raw[1:4, ],
    treatment_map = map,
    treatment_names = trt_names,
    config = tiny,
    chain_seeds = c(1101L, 1102L, 1103L, 1104L),
    single_arm = "reject",
    test_only = TRUE,
    quiet = TRUE
  ))
  stopifnot(
    inherits(smoke$samples, "mcmc.list"),
    nrow(as.matrix(smoke$samples)) == 800L,
    nrow(smoke$treatment_effects) == 3L,
    nrow(smoke$pairwise_effects) == 3L,
    smoke$metadata$n_chains == 4L,
    smoke$metadata$no_fallback_or_extension,
    isTRUE(smoke$metadata$dic_descriptive_only),
    is.finite(smoke$metadata$residual_deviance_mean)
  )
}

cat("All fixed-study model engine tests passed.\n")
