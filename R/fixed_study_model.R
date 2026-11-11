# =============================================================================
# Fixed-study-effect arm-based binomial-logit model engine
# =============================================================================

fixed_study_production_config <- function() {
  # Return a fresh value on every call so callers cannot mutate the production
  # contract for a later run.
  list(
    n.adapt = 50000L,
    n.burnin = 125000L,
    n.iter = 250000L,
    n.thin = 10L,
    n.chains = 4L
  )
}

# Short public alias for callers that do not need to distinguish configuration
# construction from production-contract validation.
fixed_study_mcmc_config <- fixed_study_production_config

validate_fixed_study_config <- function(config, test_only = FALSE) {
  required <- c("n.adapt", "n.burnin", "n.iter", "n.thin", "n.chains")
  if (!is.list(config) || !setequal(names(config), required)) {
    stop("MCMC config must contain exactly: ", paste(required, collapse = ", "), call. = FALSE)
  }

  for (field in required) {
    value <- config[[field]]
    if (length(value) != 1L || is.na(value) || !is.numeric(value) ||
        value <= 0 || value != as.integer(value)) {
      stop("MCMC config value must be a positive integer: ", field, call. = FALSE)
    }
    config[[field]] <- as.integer(value)
  }
  config <- config[required]
  if (config$n.chains != 4L) {
    stop("The fixed-study model requires exactly four chains.", call. = FALSE)
  }
  if (!isTRUE(test_only) && !identical(config, fixed_study_production_config())) {
    stop(
      "Non-production MCMC settings require test_only = TRUE; production settings are locked.",
      call. = FALSE
    )
  }
  config
}

fixed_study_seed_schedule <- function(run_names, master_seed = 12345L) {
  if (!length(run_names) || anyNA(run_names) || any(!nzchar(as.character(run_names)))) {
    stop("run_names must contain at least one non-empty scenario name.", call. = FALSE)
  }
  if (length(master_seed) != 1L || is.na(master_seed) || !is.numeric(master_seed) ||
      master_seed <= 0 || master_seed > .Machine$integer.max ||
      master_seed != as.integer(master_seed)) {
    stop("master_seed must be one positive R integer.", call. = FALSE)
  }

  sorted_runs <- sort(unique(as.character(run_names)), method = "radix")
  draw_count <- length(sorted_runs) * 4L
  if (draw_count > .Machine$integer.max) stop("Too many runs for a unique seed schedule.", call. = FALSE)

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(master_seed))
  seeds <- sample.int(.Machine$integer.max, draw_count, replace = FALSE)
  schedule <- expand.grid(
    chain = seq_len(4L),
    scenario_name = sorted_runs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  schedule <- schedule[order(match(schedule$scenario_name, sorted_runs), schedule$chain), , drop = FALSE]
  schedule$seed <- as.integer(seeds)
  schedule <- schedule[, c("scenario_name", "chain", "seed")]
  rownames(schedule) <- NULL
  schedule
}

fixed_study_chain_seeds <- function(scenario_name, chain_seeds = NULL,
                                    scheduled_run_names = scenario_name,
                                    master_seed = 12345L) {
  if (is.null(chain_seeds)) {
    schedule <- fixed_study_seed_schedule(scheduled_run_names, master_seed)
    if (!scenario_name %in% schedule$scenario_name) {
      stop("scenario_name is absent from scheduled_run_names.", call. = FALSE)
    }
    chain_seeds <- schedule$seed[schedule$scenario_name == scenario_name]
  } else if (is.data.frame(chain_seeds)) {
    required <- c("scenario_name", "chain", "seed")
    if (!all(required %in% names(chain_seeds))) {
      stop("Seed schedule must contain scenario_name, chain, and seed.", call. = FALSE)
    }
    rows <- chain_seeds[as.character(chain_seeds$scenario_name) == scenario_name, , drop = FALSE]
    rows <- rows[order(rows$chain), , drop = FALSE]
    if (!identical(as.integer(rows$chain), 1:4)) {
      stop("Seed schedule must contain chains 1 through 4 exactly once for the scenario.", call. = FALSE)
    }
    chain_seeds <- rows$seed
  }

  if (length(chain_seeds) != 4L || anyNA(chain_seeds) || !is.numeric(chain_seeds) ||
      any(chain_seeds <= 0) || any(chain_seeds > .Machine$integer.max) ||
      any(chain_seeds != as.integer(chain_seeds)) || anyDuplicated(chain_seeds)) {
    stop("chain_seeds must be four unique positive R integers.", call. = FALSE)
  }
  as.integer(chain_seeds)
}

fixed_study_network_connected <- function(treatment_ids, study_ids) {
  nodes <- sort(unique(treatment_ids))
  if (length(nodes) <= 1L) return(TRUE)
  neighbours <- setNames(vector("list", length(nodes)), as.character(nodes))
  for (study in unique(study_ids)) {
    members <- unique(treatment_ids[study_ids == study])
    for (member in members) {
      neighbours[[as.character(member)]] <- unique(c(
        neighbours[[as.character(member)]], members[members != member]
      ))
    }
  }
  visited <- nodes[1]
  frontier <- nodes[1]
  while (length(frontier)) {
    current <- frontier[1]
    frontier <- frontier[-1]
    new_nodes <- setdiff(neighbours[[as.character(current)]], visited)
    visited <- c(visited, new_nodes)
    frontier <- c(frontier, new_nodes)
  }
  setequal(visited, nodes)
}

prepare_fixed_study_data <- function(raw_rows, treatment_map, treatment_names,
                                     single_arm = c("drop", "reject")) {
  single_arm <- match.arg(single_arm)
  required <- c("Study", "T", "N", "R")
  if (!is.data.frame(raw_rows) || !all(required %in% names(raw_rows))) {
    stop("Audited rows must be a data frame containing Study, T, N, and R.", call. = FALSE)
  }
  raw <- raw_rows[, required, drop = FALSE]
  if (!nrow(raw)) stop("Audited rows are empty.", call. = FALSE)
  if (anyNA(raw) || any(!nzchar(trimws(as.character(raw$Study)))) ||
      any(!nzchar(trimws(as.character(raw$T))))) {
    stop("Study, T, N, and R must be complete and study/treatment labels non-empty.", call. = FALSE)
  }
  if (anyDuplicated(raw[c("Study", "T")])) {
    stop("Each Study/T arm must occur exactly once in audited rows.", call. = FALSE)
  }

  if (is.null(names(treatment_map)) || any(!nzchar(names(treatment_map))) ||
      anyDuplicated(names(treatment_map))) {
    stop("treatment_map must be a named vector with unique source labels.", call. = FALSE)
  }
  if (!is.character(treatment_names) || !length(treatment_names) || anyNA(treatment_names) ||
      any(!nzchar(treatment_names)) || anyDuplicated(treatment_names)) {
    stop("treatment_names must be a unique, non-empty character vector.", call. = FALSE)
  }
  map_values <- suppressWarnings(as.numeric(treatment_map))
  if (anyNA(map_values) || any(map_values != as.integer(map_values)) ||
      !setequal(unique(as.integer(map_values)), seq_along(treatment_names))) {
    stop("treatment_map IDs must cover treatment_names as consecutive IDs starting at 1.", call. = FALSE)
  }
  names(map_values) <- names(treatment_map)
  unmapped <- setdiff(unique(as.character(raw$T)), names(map_values))
  if (length(unmapped)) {
    stop("Unmapped treatment label(s): ", paste(sort(unmapped), collapse = ", "), call. = FALSE)
  }

  numeric_n <- suppressWarnings(as.numeric(as.character(raw$N)))
  numeric_r <- suppressWarnings(as.numeric(as.character(raw$R)))
  if (anyNA(numeric_n) || anyNA(numeric_r) || any(numeric_n != as.integer(numeric_n)) ||
      any(numeric_r != as.integer(numeric_r)) || any(numeric_n <= 0) ||
      any(numeric_r < 0) || any(numeric_r > numeric_n)) {
    stop("N and R must be integer counts satisfying N > 0 and 0 <= R <= N.", call. = FALSE)
  }

  raw$Study <- as.character(raw$Study)
  raw$T <- as.character(raw$T)
  raw$N <- as.integer(numeric_n)
  raw$R <- as.integer(numeric_r)
  raw$treatment <- as.integer(map_values[raw$T])
  if (anyDuplicated(raw[c("Study", "treatment")])) {
    stop("Treatment mapping creates duplicate Study/treatment arms.", call. = FALSE)
  }

  arm_counts <- table(raw$Study)
  single_studies <- sort(names(arm_counts)[arm_counts == 1L], method = "radix")
  if (length(single_studies) && single_arm == "reject") {
    stop(
      "Single-arm studies are not identifiable with fixed study effects: ",
      paste(single_studies, collapse = ", "),
      call. = FALSE
    )
  }
  kept <- if (single_arm == "drop") !raw$Study %in% single_studies else rep(TRUE, nrow(raw))
  analysis <- raw[kept, , drop = FALSE]
  if (!nrow(analysis)) stop("No multi-arm studies remain after applying single-arm policy.", call. = FALSE)

  # Stable IDs make the JAGS data invariant to input row order. A study need not
  # contain treatment 1: connected multi-arm non-reference studies are retained.
  analysis <- analysis[order(analysis$Study, analysis$treatment, method = "radix"), , drop = FALSE]
  study_levels <- sort(unique(analysis$Study), method = "radix")
  analysis$study <- match(analysis$Study, study_levels)
  analysis$arm <- seq_len(nrow(analysis))
  rownames(analysis) <- NULL

  missing_treatment_ids <- setdiff(seq_along(treatment_names), unique(analysis$treatment))
  if (length(missing_treatment_ids)) {
    stop(
      "Treatment node(s) absent after single-arm handling: ",
      paste(
        sprintf("%d (%s)", missing_treatment_ids, treatment_names[missing_treatment_ids]),
        collapse = ", "
      ),
      ". Prior-only treatment effects are not permitted.",
      call. = FALSE
    )
  }
  network_connected <- fixed_study_network_connected(analysis$treatment, analysis$study)
  if (!network_connected) {
    stop("Treatment network is disconnected after single-arm handling.", call. = FALSE)
  }

  jags_data <- list(
    Narm = nrow(analysis),
    Nstudy = length(study_levels),
    Ntrt = length(treatment_names),
    r = analysis$R,
    n = analysis$N,
    study = as.integer(analysis$study),
    treatment = as.integer(analysis$treatment)
  )
  list(
    arms = analysis,
    jags_data = jags_data,
    treatment_names = treatment_names,
    treatment_map = treatment_map,
    study_levels = study_levels,
    single_arm_policy = single_arm,
    dropped_single_arm_studies = if (single_arm == "drop") single_studies else character(),
    network_connected = network_connected
  )
}

fixed_study_pair_names <- function(n_treatments) {
  if (n_treatments < 2L) return(character())
  unlist(lapply(2:n_treatments, function(i) {
    sprintf("log_contrast[%d,%d]", i, seq_len(i - 1L))
  }), use.names = FALSE)
}

augment_fixed_study_samples <- function(samples, n_treatments) {
  if (!inherits(samples, "mcmc.list")) stop("samples must be a coda mcmc.list.", call. = FALSE)
  if (!requireNamespace("coda", quietly = TRUE)) stop("Package 'coda' is required.", call. = FALSE)
  d_names <- sprintf("d[%d]", seq_len(n_treatments))
  missing <- setdiff(d_names, colnames(samples[[1]]))
  if (length(missing)) stop("Missing monitored treatment effect(s): ", paste(missing, collapse = ", "), call. = FALSE)

  augmented <- lapply(samples, function(chain) {
    values <- as.matrix(chain)
    additions <- matrix(numeric(nrow(values) * 0L), nrow = nrow(values), ncol = 0L)
    for (i in seq_len(n_treatments)) {
      additions <- cbind(additions, exp(values[, d_names[i]]))
      colnames(additions)[ncol(additions)] <- sprintf("OR[%d]", i)
    }
    if (n_treatments >= 2L) {
      for (i in 2:n_treatments) {
        for (j in seq_len(i - 1L)) {
          contrast <- values[, d_names[i]] - values[, d_names[j]]
          additions <- cbind(additions, contrast, exp(contrast))
          colnames(additions)[(ncol(additions) - 1L):ncol(additions)] <- c(
            sprintf("log_contrast[%d,%d]", i, j),
            sprintf("pairwise_OR[%d,%d]", i, j)
          )
        }
      }
    }
    coda::mcmc(
      cbind(values, additions),
      start = stats::start(chain),
      thin = coda::thin(chain)
    )
  })
  do.call(coda::mcmc.list, augmented)
}

fixed_study_interval <- function(values) {
  quantiles <- stats::quantile(values, c(0.5, 0.025, 0.975), names = FALSE)
  c(median = quantiles[1], lower_95_cri = quantiles[2], upper_95_cri = quantiles[3])
}

summarize_fixed_study_posterior <- function(samples, treatment_names) {
  matrix_samples <- as.matrix(samples)
  n_treatments <- length(treatment_names)
  treatment_rows <- lapply(seq_len(n_treatments), function(i) {
    d_interval <- fixed_study_interval(matrix_samples[, sprintf("d[%d]", i)])
    or_interval <- fixed_study_interval(matrix_samples[, sprintf("OR[%d]", i)])
    data.frame(
      treatment_id = i,
      treatment = treatment_names[i],
      reference = treatment_names[1],
      log_or_median = unname(d_interval["median"]),
      log_or_lower_95_cri = unname(d_interval["lower_95_cri"]),
      log_or_upper_95_cri = unname(d_interval["upper_95_cri"]),
      or_median = unname(or_interval["median"]),
      or_lower_95_cri = unname(or_interval["lower_95_cri"]),
      or_upper_95_cri = unname(or_interval["upper_95_cri"]),
      stringsAsFactors = FALSE
    )
  })

  pair_rows <- list()
  if (n_treatments >= 2L) {
    for (i in 2:n_treatments) {
      for (j in seq_len(i - 1L)) {
        log_values <- matrix_samples[, sprintf("log_contrast[%d,%d]", i, j)]
        or_values <- matrix_samples[, sprintf("pairwise_OR[%d,%d]", i, j)]
        log_interval <- fixed_study_interval(log_values)
        or_interval <- fixed_study_interval(or_values)
        pair_rows[[length(pair_rows) + 1L]] <- data.frame(
          numerator_id = i,
          numerator = treatment_names[i],
          denominator_id = j,
          denominator = treatment_names[j],
          log_or_median = unname(log_interval["median"]),
          log_or_lower_95_cri = unname(log_interval["lower_95_cri"]),
          log_or_upper_95_cri = unname(log_interval["upper_95_cri"]),
          or_median = unname(or_interval["median"]),
          or_lower_95_cri = unname(or_interval["lower_95_cri"]),
          or_upper_95_cri = unname(or_interval["upper_95_cri"]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  pairwise <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame()

  heterogeneity <- do.call(rbind, lapply(c("tau_arm", "tau_contrast"), function(parameter) {
    interval <- fixed_study_interval(matrix_samples[, parameter])
    data.frame(
      parameter = parameter,
      median = unname(interval["median"]),
      lower_95_cri = unname(interval["lower_95_cri"]),
      upper_95_cri = unname(interval["upper_95_cri"]),
      stringsAsFactors = FALSE
    )
  }))
  list(
    treatment_effects = do.call(rbind, treatment_rows),
    pairwise_effects = pairwise,
    heterogeneity = heterogeneity
  )
}

fixed_study_diagnostics <- function(samples, n_treatments,
                                    psrf_threshold = 1.05,
                                    ess_threshold = 1000) {
  if (!inherits(samples, "mcmc.list")) stop("samples must be a coda mcmc.list.", call. = FALSE)
  if (!requireNamespace("coda", quietly = TRUE)) stop("Package 'coda' is required.", call. = FALSE)
  nonreference_effects <- if (n_treatments >= 2L) {
    sprintf("d[%d]", 2:n_treatments)
  } else {
    character()
  }
  parameters <- unique(c(
    nonreference_effects,
    fixed_study_pair_names(n_treatments),
    "tau_arm", "tau_contrast"
  ))
  missing <- setdiff(parameters, colnames(samples[[1]]))
  if (length(missing)) stop("Missing diagnostic parameter(s): ", paste(missing, collapse = ", "), call. = FALSE)

  rows <- lapply(parameters, function(parameter) {
    one_parameter <- do.call(coda::mcmc.list, lapply(samples, function(chain) {
      coda::mcmc(as.matrix(chain)[, parameter, drop = FALSE],
                 start = stats::start(chain), thin = coda::thin(chain))
    }))
    psrf <- tryCatch(
      coda::gelman.diag(
        one_parameter, confidence = 0.95, transform = FALSE,
        autoburnin = FALSE, multivariate = FALSE
      )$psrf[1, c("Point est.", "Upper C.I.")],
      error = function(e) c("Point est." = NA_real_, "Upper C.I." = NA_real_)
    )
    ess <- tryCatch(as.numeric(coda::effectiveSize(one_parameter)[1]), error = function(e) NA_real_)
    data.frame(
      parameter = parameter,
      parameter_type = if (grepl("^tau_", parameter)) "heterogeneity" else "treatment_effect",
      psrf_point = unname(psrf["Point est."]),
      psrf_upper_95 = unname(psrf["Upper C.I."]),
      effective_size = ess,
      psrf_pass = is.finite(psrf["Upper C.I."]) && psrf["Upper C.I."] <= psrf_threshold,
      ess_pass = is.finite(ess) && ess >= ess_threshold,
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  table$passed <- table$psrf_pass & table$ess_pass
  list(
    table = table,
    converged = all(table$passed),
    max_psrf_upper_95 = if (all(is.na(table$psrf_upper_95))) NA_real_ else max(table$psrf_upper_95, na.rm = TRUE),
    min_effective_size = if (all(is.na(table$effective_size))) NA_real_ else min(table$effective_size, na.rm = TRUE),
    psrf_threshold = psrf_threshold,
    ess_threshold = ess_threshold
  )
}

extract_fixed_study_p_draws <- function(samples, n_arms) {
  matrix_samples <- as.matrix(samples)
  expected <- sprintf("p[%d]", seq_len(n_arms))
  missing <- setdiff(expected, colnames(matrix_samples))
  if (length(missing)) stop("Missing monitored arm probability(s): ", paste(missing, collapse = ", "), call. = FALSE)
  matrix_samples[, expected, drop = FALSE]
}

fixed_study_residual_deviance <- function(p_draws, events, participants) {
  p_draws <- as.matrix(p_draws)
  events <- as.numeric(events)
  participants <- as.numeric(participants)
  if (ncol(p_draws) != length(events) || length(events) != length(participants)) {
    stop("p_draws columns must align one-for-one with events and participants.", call. = FALSE)
  }
  if (anyNA(p_draws) || any(!is.finite(p_draws)) || any(p_draws < 0 | p_draws > 1) ||
      anyNA(events) || anyNA(participants) || any(events < 0) ||
      any(participants <= 0) || any(events > participants)) {
    stop("Invalid probabilities or binomial counts for residual deviance.", call. = FALSE)
  }

  # Observed zero-event and all-event arms have valid saturated contributions:
  # 0 * log(0 / expected) is defined as zero. Clipping only protects log(0)
  # for numerical probability boundaries and does not continuity-correct counts.
  epsilon <- .Machine$double.eps
  safe_p <- pmin(pmax(p_draws, epsilon), 1 - epsilon)
  result <- matrix(0, nrow = nrow(safe_p), ncol = ncol(safe_p))
  for (arm in seq_along(events)) {
    r <- events[arm]
    n <- participants[arm]
    event_term <- if (r == 0) 0 else r * log(r / (n * safe_p[, arm]))
    non_event_term <- if (r == n) 0 else (n - r) * log((n - r) / (n * (1 - safe_p[, arm])))
    result[, arm] <- 2 * (event_term + non_event_term)
  }
  colnames(result) <- sprintf("arm[%d]", seq_len(ncol(result)))
  total <- rowSums(result)
  arm_summary <- data.frame(
    arm = seq_len(ncol(result)),
    mean = colMeans(result),
    median = apply(result, 2, stats::median),
    lower_95_cri = apply(result, 2, stats::quantile, probs = 0.025, names = FALSE),
    upper_95_cri = apply(result, 2, stats::quantile, probs = 0.975, names = FALSE),
    stringsAsFactors = FALSE
  )
  total_interval <- fixed_study_interval(total)
  list(
    draws_by_arm = result,
    total_draws = total,
    arm_summary = arm_summary,
    total_summary = data.frame(
      mean = mean(total),
      median = unname(total_interval["median"]),
      lower_95_cri = unname(total_interval["lower_95_cri"]),
      upper_95_cri = unname(total_interval["upper_95_cri"]),
      per_arm_mean = mean(total) / ncol(result),
      stringsAsFactors = FALSE
    )
  )
}

fixed_study_descriptive_dic <- function(p_draws, events, participants) {
  p_draws <- as.matrix(p_draws)
  epsilon <- .Machine$double.eps
  safe_p <- pmin(pmax(p_draws, epsilon), 1 - epsilon)
  deviance_draws <- -2 * rowSums(vapply(seq_along(events), function(arm) {
    stats::dbinom(events[arm], size = participants[arm], prob = safe_p[, arm], log = TRUE)
  }, numeric(nrow(safe_p))))
  dbar <- mean(deviance_draws)
  p_hat <- colMeans(safe_p)
  dhat <- -2 * sum(stats::dbinom(events, size = participants, prob = p_hat, log = TRUE))
  p_d <- dbar - dhat
  list(
    Dbar = dbar,
    Dhat = dhat,
    pD = p_d,
    DIC = dbar + p_d,
    descriptive_only = TRUE,
    method = "posterior mean deviance minus deviance at posterior mean arm probabilities; descriptive only"
  )
}

fit_fixed_study_model <- function(
  scenario_name,
  raw_rows,
  treatment_map,
  treatment_names,
  model_file = file.path("R", "model_fixed_study", "model.jags"),
  config = fixed_study_production_config(),
  chain_seeds = NULL,
  single_arm = c("drop", "reject"),
  master_seed = 12345L,
  scheduled_run_names = scenario_name,
  preflight = FALSE,
  test_only = FALSE,
  quiet = TRUE,
  output_controls = NULL
) {
  if (length(scenario_name) != 1L || is.na(scenario_name) || !nzchar(scenario_name)) {
    stop("scenario_name must be one non-empty string.", call. = FALSE)
  }
  single_arm <- match.arg(single_arm)
  config <- validate_fixed_study_config(config, test_only = test_only)
  seeds <- fixed_study_chain_seeds(
    scenario_name, chain_seeds, scheduled_run_names, master_seed
  )
  prepared <- prepare_fixed_study_data(raw_rows, treatment_map, treatment_names, single_arm)
  if (!file.exists(model_file)) stop("Fixed-study JAGS model is unavailable: ", model_file, call. = FALSE)

  controls <- list(
    progress_bar = if (isTRUE(quiet)) "none" else "text",
    messages = !isTRUE(quiet)
  )
  if (!is.null(output_controls)) {
    if (!is.list(output_controls) || !all(names(output_controls) %in% names(controls))) {
      stop("output_controls may contain only progress_bar and messages.", call. = FALSE)
    }
    controls[names(output_controls)] <- output_controls
  }
  if (!controls$progress_bar %in% c("none", "text", "gui") ||
      length(controls$messages) != 1L || !is.logical(controls$messages) ||
      is.na(controls$messages)) {
    stop("Invalid output_controls.", call. = FALSE)
  }

  preflight_result <- list(
    scenario_name = scenario_name,
    model_file = normalizePath(model_file, mustWork = TRUE),
    config = config,
    chain_seeds = seeds,
    prepared_data = prepared,
    test_only = isTRUE(test_only),
    preflight = TRUE
  )
  if (isTRUE(preflight)) return(preflight_result)
  if (!requireNamespace("rjags", quietly = TRUE) || !requireNamespace("coda", quietly = TRUE)) {
    stop("Packages 'rjags' and 'coda' are required to fit the fixed-study model.", call. = FALSE)
  }

  if (isTRUE(controls$messages)) {
    message("Fitting fixed-study model for ", scenario_name, " with four prespecified chains.")
  }
  inits <- lapply(seeds, function(seed) {
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)
  })
  model <- rjags::jags.model(
    file = model_file,
    data = prepared$jags_data,
    inits = inits,
    n.chains = config$n.chains,
    n.adapt = config$n.adapt,
    quiet = isTRUE(quiet)
  )
  stats::update(model, n.iter = config$n.burnin, progress.bar = controls$progress_bar)
  monitored <- c("d", "tau_arm", "tau_contrast", "p")
  raw_samples <- rjags::coda.samples(
    model = model,
    variable.names = monitored,
    n.iter = config$n.iter,
    thin = config$n.thin,
    progress.bar = controls$progress_bar
  )

  samples <- augment_fixed_study_samples(raw_samples, length(treatment_names))
  posterior <- summarize_fixed_study_posterior(samples, treatment_names)
  diagnostics <- fixed_study_diagnostics(samples, length(treatment_names))
  p_draws <- extract_fixed_study_p_draws(samples, prepared$jags_data$Narm)
  residual <- fixed_study_residual_deviance(
    p_draws, prepared$arms$R, prepared$arms$N
  )
  dic <- fixed_study_descriptive_dic(p_draws, prepared$arms$R, prepared$arms$N)
  retained_draws <- as.integer(config$n.iter / config$n.thin * config$n.chains)
  warnings <- character()
  if (length(prepared$dropped_single_arm_studies)) {
    warnings <- sprintf(
      "Dropped non-identifiable single-arm study strata: %s",
      paste(prepared$dropped_single_arm_studies, collapse = ", ")
    )
  }

  list(
    scenario_name = scenario_name,
    model_type = "fixed_study_arm_heterogeneity",
    samples = samples,
    monitored_samples = raw_samples,
    posterior = posterior,
    treatment_effects = posterior$treatment_effects,
    pairwise_effects = posterior$pairwise_effects,
    heterogeneity = posterior$heterogeneity,
    diagnostics = diagnostics,
    diagnostic_table = diagnostics$table,
    residual_deviance = residual,
    residual_deviance_arms = residual$arm_summary,
    residual_deviance_summary = residual$total_summary,
    dic = dic,
    warnings = warnings,
    prepared_data = prepared,
    metadata = list(
      scenario_name = scenario_name,
      model_type = "fixed_study_arm_heterogeneity",
      n_arms = prepared$jags_data$Narm,
      n_studies = prepared$jags_data$Nstudy,
      n_treatments = prepared$jags_data$Ntrt,
      n_adapt = config$n.adapt,
      n_burnin = config$n.burnin,
      n_iter = config$n.iter,
      n_thin = config$n.thin,
      n_chains = config$n.chains,
      retained_draws = retained_draws,
      master_seed = as.integer(master_seed),
      chain_seeds = seeds,
      chain_seed_1 = seeds[1],
      chain_seed_2 = seeds[2],
      chain_seed_3 = seeds[3],
      chain_seed_4 = seeds[4],
      monitored_parameters = paste(monitored, collapse = "|"),
      production_mcmc_contract = identical(config, fixed_study_production_config()),
      single_arm_policy = single_arm,
      dropped_single_arm_studies = paste(prepared$dropped_single_arm_studies, collapse = "|"),
      network_connected = prepared$network_connected,
      max_psrf_upper_95 = diagnostics$max_psrf_upper_95,
      min_effective_size = diagnostics$min_effective_size,
      converged = diagnostics$converged,
      residual_deviance_mean = residual$total_summary$mean,
      residual_deviance_per_arm = residual$total_summary$per_arm_mean,
      dic = dic$DIC,
      dic_descriptive_only = TRUE,
      no_fallback_or_extension = TRUE
    )
  )
}
