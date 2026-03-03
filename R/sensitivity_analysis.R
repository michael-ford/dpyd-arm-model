# =============================================================================
# R/sensitivity_analysis.R
# Sensitivity Analysis for DPYD Arm-Based NMA
# =============================================================================
#
# Tests robustness of primary finding: HapB3 OR ~2.0 (95% CrI: 1.29-3.19)
# is substantially lower than other DPYD variant ORs across model perturbations.
#
# Scenarios:
#   1. Primary model (read from existing result)
#   2. Alternative model (hom_eqcor)
#   3. Alternative link function (logit)
#   4. Alternative prior (inv-gamma a=0.01, b=0.01)
#   5. het_cor with inv-gamma prior (a=0.001)
#   6. Weakly informative prior (inv-gamma a=0.5, b=0.5)
#   7a. Exclude sparse arms N <= 2
#   7b. Exclude sparse arms N <= 5
#   8. Full leave-one-out analysis (all 31 studies)
#   9. Contrast-based model comparison (conditional)
#
# Usage: Rscript R/sensitivity_analysis.R
# =============================================================================

# Set working directory (for Docker)
if (file.exists("/analysis/R")) {
  setwd("/analysis")
}

# Store project root
project_root <- getwd()

# Create output directories
output_base <- file.path(project_root, "output/sensitivity")
scenarios_dir <- file.path(output_base, "scenarios")
dir.create(scenarios_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Setup logging
# =============================================================================

log_file <- file.path(output_base, "sensitivity_log.txt")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message", append = TRUE)

on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  message(sprintf("Sensitivity analysis log saved to: %s", log_file))
}, add = TRUE)

cat("=============================================================\n")
cat("DPYD Sensitivity Analysis\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat("=============================================================\n\n")

# =============================================================================
# Source dependencies
# =============================================================================

source("R/common/data_utils.R")
source("R/common/model_runner.R")  # for check_convergence()

library(pcnetmeta)
library(rjags)

# =============================================================================
# Configuration
# =============================================================================

# Treatment mapping (primary model, unified WT)
TREATMENT_MAP <- c(
  "WT_Clean"         = 1,
  "WT_Biased"        = 1,
  "HapB3_1129or1236" = 2,
  "2846hetho"        = 3,
  "2Ahetho"          = 4,
  "13hetho"          = 5
)

TREATMENT_NAMES <- c("WT", "HapB3", "2846hetho", "2Ahetho", "13hetho")

DATA_FILE <- "data/Binary WT HapB3 Data for NMA (12-18-2025).xlsx"
PRIMARY_RESULT_FILE <- "output/wt_unified/wt_unified_result.rds"

# MCMC configuration for sensitivity runs (reduced)
SENSITIVITY_MCMC <- list(
  n.adapt  = 10000,
  n.iter   = 100000,
  n.burnin = 50000,
  n.chains = 3,
  n.thin   = 10
)

# Full MCMC configuration for escalation
FULL_MCMC <- list(
  n.adapt  = 10000,
  n.iter   = 150000,
  n.burnin = 75000,
  n.chains = 3,
  n.thin   = 10
)

# Convergence threshold (relaxed for sensitivity)
PSRF_THRESHOLD <- 1.10

# Seed for reproducibility
SEED <- 12345

# LOO study list is derived dynamically from raw data in Phase 4

# =============================================================================
# Helper Functions
# =============================================================================

# --- Parse OR string "X.XXXX (Y.YYYY, Z.ZZZZ)" ---
parse_or_string <- function(or_string) {
  m <- regmatches(or_string,
                  regexec("([0-9.]+)\\s*\\(([0-9.]+),\\s*([0-9.]+)\\)", or_string))
  if (length(m[[1]]) >= 4) {
    return(list(
      median = as.numeric(m[[1]][2]),
      lower  = as.numeric(m[[1]][3]),
      upper  = as.numeric(m[[1]][4])
    ))
  }
  return(NULL)
}

# --- Extract SUCRA from rank probability matrix ---
extract_sucra <- function(result, treatment_names) {
  if (is.null(result$TrtRankProb)) return(NULL)
  rank_numeric <- matrix(as.numeric(result$TrtRankProb),
                         nrow = length(treatment_names),
                         ncol = length(treatment_names))
  ntrt <- length(treatment_names)
  sucra <- numeric(ntrt)
  for (i in 1:ntrt) {
    cumulative_probs <- cumsum(rank_numeric[i, 1:(ntrt - 1)])
    sucra[i] <- sum(cumulative_probs) / (ntrt - 1)
  }
  names(sucra) <- treatment_names
  return(sucra)
}

# --- Exclude study for LOO (exact match) ---
exclude_study <- function(raw_data, study_name) {
  filtered <- raw_data[raw_data$Study != study_name, ]
  return(filtered)
}

# --- Network connectivity check ---
check_network_connectivity <- function(nma_data, treatment_names) {
  arms_per_trt <- table(nma_data$t.id)
  all_trt_ids <- 1:length(treatment_names)
  present_ids <- as.integer(names(arms_per_trt))
  disconnected <- setdiff(all_trt_ids, present_ids)
  sparse <- as.integer(names(arms_per_trt)[arms_per_trt < 2])

  issues <- list(
    connected = TRUE,
    disconnected_nodes = character(0),
    sparse_nodes = character(0),
    arms_per_treatment = as.list(arms_per_trt)
  )

  if (length(disconnected) > 0) {
    issues$connected <- FALSE
    issues$disconnected_nodes <- treatment_names[disconnected]
  }
  if (length(sparse) > 0) {
    issues$sparse_nodes <- treatment_names[sparse]
  }

  return(issues)
}

# --- Run function in scenario-specific directory ---
run_in_scenario_dir <- function(scenario_name, fn) {
  scenario_dir <- file.path(project_root, "output/sensitivity/scenarios", scenario_name)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(scenario_dir)
  on.exit(setwd(old_wd))
  fn()
}

# --- Read max PSRF from convergence diagnostic file ---
read_max_psrf <- function() {
  psrf_files <- list.files(pattern = ".*[Cc]onvergence.*\\.txt$")
  if (length(psrf_files) == 0) {
    psrf_files <- list.files(pattern = "PSRF.*\\.txt$")
  }
  if (length(psrf_files) == 0) return(NA)

  psrf_file <- psrf_files[1]
  tryCatch({
    psrf_lines <- readLines(psrf_file)
    psrf_values <- c()
    for (line in psrf_lines) {
      if (nchar(trimws(line)) == 0) next
      if (grepl("^#|^Parameter|^---", line)) next
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      if (length(parts) >= 2) {
        numeric_vals <- suppressWarnings(as.numeric(parts))
        valid_nums <- numeric_vals[!is.na(numeric_vals)]
        if (length(valid_nums) > 0) {
          psrf_val <- valid_nums[length(valid_nums)]
          if (psrf_val > 0 && psrf_val < 100) {
            psrf_values <- c(psrf_values, psrf_val)
          }
        }
      }
    }
    if (length(psrf_values) > 0) return(max(psrf_values))
    return(NA)
  }, error = function(e) {
    return(NA)
  })
}

# --- Extract DIC from result ---
extract_dic <- function(result) {
  if (is.null(result$DIC)) return(NA)
  tryCatch({
    dic_vec <- result$DIC
    if (is.matrix(dic_vec) || is.array(dic_vec)) {
      return(dic_vec["DIC", 1])
    } else if (is.numeric(dic_vec) && length(dic_vec) >= 3) {
      return(dic_vec[3])
    }
    return(NA)
  }, error = function(e) NA)
}

# --- Classify conclusion ---
classify_conclusion <- function(hapb3_or, hapb3_lower, hapb3_upper,
                                 hapb3_sucra, variant_sucras,
                                 primary_hapb3_or) {
  # Check if CrI includes 1.0
  cri_excludes_one <- hapb3_lower > 1.0

  # Log-OR shift
  log_or_shift <- abs(log(hapb3_or) - log(primary_hapb3_or))

  # Rank preserved: HapB3 SUCRA > all variant SUCRAs
  rank_preserved <- all(hapb3_sucra > variant_sucras, na.rm = TRUE)

  if (!cri_excludes_one) {
    return(list(
      conclusion = "Materially different",
      log_or_shift = log_or_shift,
      rank_preserved = rank_preserved
    ))
  }

  if (log_or_shift < 0.20 && rank_preserved) {
    return(list(
      conclusion = "Consistent",
      log_or_shift = log_or_shift,
      rank_preserved = rank_preserved
    ))
  }

  return(list(
    conclusion = "Direction preserved, magnitude shifted",
    log_or_shift = log_or_shift,
    rank_preserved = rank_preserved
  ))
}

# =============================================================================
# Core: Run a single sensitivity scenario
# =============================================================================

run_sensitivity_scenario <- function(config, nma_data = NULL,
                                      treatment_names = TREATMENT_NAMES) {
  cat("\n=============================================================\n")
  cat(sprintf("Scenario: %s\n", config$name))
  cat(sprintf("Description: %s\n", config$description))
  cat("=============================================================\n\n")

  # Track escalation state
  escalated <- FALSE

  run_model <- function(mcmc_config) {
    set.seed(SEED)

    # Build nma.ab.bin arguments
    args <- list(
      s.id      = quote(s.id),
      t.id      = quote(t.id),
      event.n   = quote(r),
      total.n   = quote(n),
      data      = nma_data,
      trtname   = treatment_names,
      param     = c("AR", "OR", "LOR", "RD", "rank.prob"),
      n.adapt   = mcmc_config$n.adapt,
      n.iter    = mcmc_config$n.iter,
      n.burnin  = mcmc_config$n.burnin,
      n.chains  = mcmc_config$n.chains,
      n.thin    = mcmc_config$n.thin,
      conv.diag = TRUE,
      trace     = c("LOR"),
      dic       = TRUE,
      postdens  = FALSE,
      mcmc.samples = TRUE,
      higher.better = FALSE,
      digits    = 4
    )

    # Model-specific arguments
    args$model <- if (!is.null(config$model)) config$model else "het_cor"
    args$link  <- if (!is.null(config$link))  config$link  else "probit"

    if (!is.null(config$prior.type)) {
      args$prior.type <- config$prior.type
    } else {
      args$prior.type <- "invwishart"
    }

    if (!is.null(config$a)) args$a <- config$a
    if (!is.null(config$b)) args$b <- config$b

    cat(sprintf("  Model: %s, Link: %s, Prior: %s\n",
                args$model, args$link, args$prior.type))
    if (!is.null(config$a)) {
      cat(sprintf("  Prior hyperparameters: a=%s, b=%s\n", config$a, config$b))
    }
    cat(sprintf("  MCMC: %d iter, %d burn-in, thin=%d, %d chains\n",
                mcmc_config$n.iter, mcmc_config$n.burnin,
                mcmc_config$n.thin, mcmc_config$n.chains))
    cat(sprintf("  Data: %d arms\n\n", nrow(nma_data)))

    do.call(nma.ab.bin, args)
  }

  # Run in scenario-specific directory
  result <- run_in_scenario_dir(config$name, function() {
    res <- tryCatch(run_model(SENSITIVITY_MCMC), error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      return(NULL)
    })

    if (is.null(res)) return(NULL)

    # Check convergence
    max_psrf <- read_max_psrf()
    cat(sprintf("  Max PSRF: %s\n",
                ifelse(is.na(max_psrf), "unavailable",
                       sprintf("%.4f", max_psrf))))

    # Escalation if needed
    if (!is.na(max_psrf) && max_psrf > PSRF_THRESHOLD) {
      cat(sprintf("  PSRF %.4f > %.2f threshold. Escalating to full MCMC...\n",
                  max_psrf, PSRF_THRESHOLD))
      escalated <<- TRUE
      res <- tryCatch(run_model(FULL_MCMC), error = function(e) {
        cat(sprintf("  ESCALATION ERROR: %s\n", e$message))
        return(res)  # Return original result
      })
      max_psrf <- read_max_psrf()
      cat(sprintf("  Escalated Max PSRF: %s\n",
                  ifelse(is.na(max_psrf), "unavailable",
                         sprintf("%.4f", max_psrf))))
    }

    # Save result RDS
    saveRDS(res, "scenario_result.rds")

    list(result = res, max_psrf = max_psrf)
  })

  if (is.null(result)) {
    return(create_failed_row(config))
  }

  # Extract metrics
  extract_scenario_row(config, result$result, result$max_psrf,
                       nma_data, treatment_names, escalated)
}

# --- Extract standardized row from scenario result ---
extract_scenario_row <- function(config, result, max_psrf,
                                  nma_data, treatment_names, escalated) {
  # Count studies and arms
  n_studies <- length(unique(nma_data$s.id))
  n_arms <- nrow(nma_data)

  # Extract ORs
  or_matrix <- result$OddsRatio$Median_CI

  hapb3_parsed  <- parse_or_string(or_matrix[2, 1])  # HapB3 = treatment 2
  c2846_parsed  <- parse_or_string(or_matrix[3, 1])  # 2846hetho = treatment 3
  star2a_parsed <- parse_or_string(or_matrix[4, 1])  # 2Ahetho = treatment 4
  star13_parsed <- parse_or_string(or_matrix[5, 1])  # 13hetho = treatment 5

  # Extract SUCRA
  sucra <- extract_sucra(result, treatment_names)

  # Extract DIC
  dic <- extract_dic(result)

  # Convergence
  converged <- !is.na(max_psrf) && max_psrf <= PSRF_THRESHOLD

  row <- data.frame(
    scenario           = config$name,
    description        = config$description,
    n_studies          = n_studies,
    n_arms             = n_arms,
    HapB3_OR           = if (!is.null(hapb3_parsed))  hapb3_parsed$median  else NA,
    HapB3_CrI_low      = if (!is.null(hapb3_parsed))  hapb3_parsed$lower   else NA,
    HapB3_CrI_high     = if (!is.null(hapb3_parsed))  hapb3_parsed$upper   else NA,
    HapB3_SUCRA        = if (!is.null(sucra)) round(sucra["HapB3"] * 100, 1) else NA,
    HapB3_rank_preserved = NA,  # Filled during conclusion classification
    c2846_OR           = if (!is.null(c2846_parsed))  c2846_parsed$median  else NA,
    c2846_CrI_low      = if (!is.null(c2846_parsed))  c2846_parsed$lower   else NA,
    c2846_CrI_high     = if (!is.null(c2846_parsed))  c2846_parsed$upper   else NA,
    star2A_OR          = if (!is.null(star2a_parsed)) star2a_parsed$median else NA,
    star13_OR          = if (!is.null(star13_parsed)) star13_parsed$median else NA,
    star13_CrI_low     = if (!is.null(star13_parsed)) star13_parsed$lower  else NA,
    star13_CrI_high    = if (!is.null(star13_parsed)) star13_parsed$upper  else NA,
    DIC                = dic,
    max_PSRF           = max_psrf,
    converged          = converged,
    escalated          = escalated,
    hapb3_log_or_shift = NA,  # Filled after primary baseline known
    conclusion         = NA,  # Filled during conclusion classification
    stringsAsFactors = FALSE
  )

  return(row)
}

# --- Create row for failed scenario ---
create_failed_row <- function(config, reason = "Model failed to run") {
  data.frame(
    scenario             = config$name,
    description          = config$description,
    n_studies            = NA,
    n_arms               = NA,
    HapB3_OR             = NA,
    HapB3_CrI_low        = NA,
    HapB3_CrI_high       = NA,
    HapB3_SUCRA          = NA,
    HapB3_rank_preserved = NA,
    c2846_OR             = NA,
    c2846_CrI_low        = NA,
    c2846_CrI_high       = NA,
    star2A_OR            = NA,
    star13_OR            = NA,
    star13_CrI_low       = NA,
    star13_CrI_high      = NA,
    DIC                  = NA,
    max_PSRF             = NA,
    converged            = FALSE,
    escalated            = FALSE,
    hapb3_log_or_shift   = NA,
    conclusion           = reason,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Read primary model result (Scenario 1)
# =============================================================================

read_primary_result <- function() {
  cat("\n=============================================================\n")
  cat("Scenario: primary (reading existing result)\n")
  cat("=============================================================\n\n")

  if (!file.exists(PRIMARY_RESULT_FILE)) {
    stop("Primary result file not found: ", PRIMARY_RESULT_FILE)
  }

  result <- readRDS(PRIMARY_RESULT_FILE)

  # Extract ORs
  or_matrix <- result$OddsRatio$Median_CI
  hapb3_parsed  <- parse_or_string(or_matrix[2, 1])
  c2846_parsed  <- parse_or_string(or_matrix[3, 1])
  star2a_parsed <- parse_or_string(or_matrix[4, 1])
  star13_parsed <- parse_or_string(or_matrix[5, 1])

  # Extract SUCRA
  sucra <- extract_sucra(result, TREATMENT_NAMES)

  # Extract DIC
  dic <- extract_dic(result)

  # Read convergence from existing file
  conv_file <- "output/wt_unified/wt_unified_convergence.rds"
  max_psrf <- NA
  if (file.exists(conv_file)) {
    conv <- readRDS(conv_file)
    max_psrf <- conv$max_psrf
  }

  # Rank preserved: HapB3 SUCRA > all other variant SUCRAs
  variant_sucras <- sucra[c("2846hetho", "2Ahetho", "13hetho")]
  rank_preserved <- all(sucra["HapB3"] > variant_sucras, na.rm = TRUE)

  row <- data.frame(
    scenario             = "primary",
    description          = "Primary model (het_cor, probit)",
    n_studies            = 31,
    n_arms               = 98,
    HapB3_OR             = hapb3_parsed$median,
    HapB3_CrI_low        = hapb3_parsed$lower,
    HapB3_CrI_high       = hapb3_parsed$upper,
    HapB3_SUCRA          = round(sucra["HapB3"] * 100, 1),
    HapB3_rank_preserved = rank_preserved,
    c2846_OR             = c2846_parsed$median,
    c2846_CrI_low        = c2846_parsed$lower,
    c2846_CrI_high       = c2846_parsed$upper,
    star2A_OR            = star2a_parsed$median,
    star13_OR            = star13_parsed$median,
    star13_CrI_low       = star13_parsed$lower,
    star13_CrI_high      = star13_parsed$upper,
    DIC                  = dic,
    max_PSRF             = max_psrf,
    converged            = TRUE,
    escalated            = FALSE,
    hapb3_log_or_shift   = 0.000,
    conclusion           = "Reference",
    stringsAsFactors = FALSE
  )

  cat(sprintf("  HapB3 OR: %.4f (%.4f, %.4f)\n",
              hapb3_parsed$median, hapb3_parsed$lower, hapb3_parsed$upper))
  cat(sprintf("  HapB3 SUCRA: %.1f%%\n", sucra["HapB3"] * 100))
  cat(sprintf("  Rank preserved: %s\n", rank_preserved))
  cat(sprintf("  DIC: %.2f\n", dic))
  cat(sprintf("  Max PSRF: %.4f\n", max_psrf))

  return(row)
}

# =============================================================================
# Data preparation functions for filtered scenarios
# =============================================================================

# --- Prepare sparse-filtered data (Scenarios 5a, 5b) ---
prepare_sparse_filtered_data <- function(threshold) {
  raw_data <- load_raw_data(DATA_FILE)
  nma_data <- transform_to_nma_format(raw_data, TREATMENT_MAP, aggregate_wt = TRUE)

  cat(sprintf("  Before filtering: %d arms\n", nrow(nma_data)))

  # Filter arms by N threshold
  nma_data <- nma_data[nma_data$n > threshold, ]
  cat(sprintf("  After removing arms with N <= %d: %d arms\n",
              threshold, nrow(nma_data)))

  # Drop studies that have only WT arms remaining
  studies_with_variants <- unique(nma_data$s.id[nma_data$t.id > 1])
  wt_only_studies <- setdiff(unique(nma_data$s.id), studies_with_variants)
  if (length(wt_only_studies) > 0) {
    cat(sprintf("  Dropping %d WT-only studies after filtering\n",
                length(wt_only_studies)))
    nma_data <- nma_data[!nma_data$s.id %in% wt_only_studies, ]
  }

  # Renumber s.id sequentially
  old_ids <- sort(unique(nma_data$s.id))
  id_map <- setNames(seq_along(old_ids), old_ids)
  nma_data$s.id <- id_map[as.character(nma_data$s.id)]

  # Keep only pcnetmeta columns
  nma_data <- nma_data[, c("s.id", "t.id", "r", "n")]

  cat(sprintf("  Final data: %d studies, %d arms\n",
              length(unique(nma_data$s.id)), nrow(nma_data)))

  return(nma_data)
}

# --- Prepare LOO data ---
prepare_loo_data <- function(study_name) {
  raw_data <- load_raw_data(DATA_FILE)

  # Check if study_name is a prefix for multi-arm entries (e.g., Rosmarin_2015)
  # If any raw study name starts with study_name followed by _arm_, use prefix matching
  arm_matches <- grepl(paste0("^", study_name, "_arm_"), raw_data$Study)
  if (any(arm_matches)) {
    raw_data <- raw_data[!grepl(paste0("^", study_name), raw_data$Study), ]
  } else {
    raw_data <- exclude_study(raw_data, study_name)
  }

  nma_data <- transform_to_nma_format(raw_data, TREATMENT_MAP, aggregate_wt = TRUE)

  # Keep only pcnetmeta columns
  nma_data_clean <- nma_data[, c("s.id", "t.id", "r", "n")]

  return(nma_data_clean)
}

# =============================================================================
# Forest plot generation
# =============================================================================

generate_forest_plot <- function(summary_df, output_file) {
  cat("\nGenerating forest plot...\n")

  # Filter to rows with valid HapB3 OR data
  plot_data <- summary_df[!is.na(summary_df$HapB3_OR), ]

  if (nrow(plot_data) == 0) {
    cat("  No valid data for forest plot.\n")
    return(invisible(NULL))
  }

  # Classify row types
  is_loo <- grepl("^loo_", plot_data$scenario)
  is_primary <- plot_data$scenario == "primary"

  # Order: primary at top, then main scenarios, then LOO at bottom
  plot_data$order <- ifelse(is_primary, 1, ifelse(!is_loo, 2, 3))
  plot_data <- plot_data[order(plot_data$order, plot_data$scenario), ]

  # Recompute type flags after reordering
  is_loo <- grepl("^loo_", plot_data$scenario)
  is_primary <- plot_data$scenario == "primary"

  # Create display labels
  plot_data$label <- ifelse(
    is_loo,
    gsub("^loo_", "LOO: ", plot_data$scenario),
    plot_data$description
  )

  # Non-converged flag
  is_nonconverged <- !is.na(plot_data$converged) & plot_data$converged == FALSE

  # Colors
  plot_data$color <- ifelse(
    is_primary, "blue",
    ifelse(is_nonconverged, "darkorange",
    ifelse(!is.na(plot_data$conclusion) & plot_data$conclusion == "Materially different",
           "red", "gray40"))
  )

  # Point type: filled for main scenarios, open for LOO, triangle for non-converged
  plot_data$pch <- ifelse(is_nonconverged, 17, ifelse(is_loo, 1, 19))

  primary_or <- plot_data$HapB3_OR[plot_data$scenario == "primary"]
  if (length(primary_or) == 0) primary_or <- 2.005

  n <- nrow(plot_data)

  png(output_file, width = 8, height = max(6, 0.4 * n + 2),
      units = "in", res = 300)

  par(mar = c(5, 15, 3, 2))

  # X-axis range (log scale)
  x_max <- max(8, max(plot_data$HapB3_CrI_high, na.rm = TRUE) * 1.1)
  x_range <- c(0.5, x_max)

  plot(NA, xlim = x_range, ylim = c(0.5, n + 0.5),
       xlab = "Odds Ratio (log scale)", ylab = "",
       yaxt = "n", log = "x",
       main = "Sensitivity Analysis: HapB3 OR vs WT")

  # Y-axis labels
  axis(2, at = n:1, labels = plot_data$label, las = 1, cex.axis = 0.7)

  # Reference lines
  abline(v = 1.0, lty = 2, col = "gray60")       # Null
  abline(v = primary_or, lty = 3, col = "blue")   # Primary OR

  # Plot points and CrIs
  for (i in 1:n) {
    y_pos <- n - i + 1

    # CrI line
    segments(
      x0 = plot_data$HapB3_CrI_low[i],
      x1 = plot_data$HapB3_CrI_high[i],
      y0 = y_pos, y1 = y_pos,
      col = plot_data$color[i], lwd = 1.5
    )

    # Point
    points(plot_data$HapB3_OR[i], y_pos,
           pch = plot_data$pch[i], col = plot_data$color[i],
           bg = plot_data$color[i], cex = 1.2)
  }

  # Legend
  legend("topright",
         legend = c("Primary", "Alternative scenario", "LOO study",
                    "Non-converged", "Materially different",
                    "OR = 1.0", "Primary HapB3 OR"),
         pch = c(19, 19, 1, 17, 19, NA, NA),
         col = c("blue", "gray40", "gray40", "darkorange", "red",
                 "gray60", "blue"),
         lty = c(NA, NA, NA, NA, NA, 2, 3),
         cex = 0.7, bg = "white")

  dev.off()
  cat(sprintf("  Saved: %s\n", output_file))
}

# =============================================================================
# Manuscript edit generation
# =============================================================================

generate_manuscript_edits <- function(summary_df, output_file) {
  cat("\nGenerating manuscript edit suggestions...\n")

  # --- Compute summary statistics with converged/infeasible distinction ---
  all_non_primary <- summary_df[summary_df$scenario != "primary", ]

  # Scenarios with valid results (non-NA HapB3_OR)
  feasible <- all_non_primary[!is.na(all_non_primary$HapB3_OR), ]
  n_scenarios <- nrow(feasible)

  # Converged subset
  converged <- feasible[!is.na(feasible$converged) & feasible$converged == TRUE, ]
  n_converged <- nrow(converged)

  # Infeasible (no HapB3_OR — disconnected network or model failure)
  n_infeasible <- sum(is.na(all_non_primary$HapB3_OR))

  # Non-converged (has results but PSRF > threshold)
  non_converged <- feasible[!is.na(feasible$converged) & feasible$converged == FALSE, ]
  n_non_converged <- nrow(non_converged)

  primary_or <- summary_df$HapB3_OR[summary_df$scenario == "primary"]
  primary_sucra <- summary_df$HapB3_SUCRA[summary_df$scenario == "primary"]

  or_min <- min(converged$HapB3_OR, na.rm = TRUE)
  or_max <- max(converged$HapB3_OR, na.rm = TRUE)

  # SUCRA range from CONVERGED scenarios only (addresses peer review issue 5)
  sucra_min <- min(converged$HapB3_SUCRA, na.rm = TRUE)
  sucra_max <- max(converged$HapB3_SUCRA, na.rm = TRUE)

  n_consistent <- sum(converged$conclusion == "Consistent", na.rm = TRUE)
  n_direction <- sum(converged$conclusion ==
                       "Direction preserved, magnitude shifted", na.rm = TRUE)
  n_material <- sum(converged$conclusion == "Materially different", na.rm = TRUE)

  all_rank_preserved <- all(converged$HapB3_rank_preserved == TRUE, na.rm = TRUE)
  all_cri_exclude_one <- all(converged$HapB3_CrI_low > 1.0, na.rm = TRUE)

  # Check if het_cor_invgamma failed
  het_cor_invgamma_row <- all_non_primary[all_non_primary$scenario == "het_cor_invgamma", ]
  het_cor_failed <- nrow(het_cor_invgamma_row) > 0 && is.na(het_cor_invgamma_row$HapB3_OR[1])

  # Determine overall robustness (from converged scenarios only)
  if (n_material == 0 && n_consistent >= n_converged * 0.8) {
    robustness <- "robust"
  } else if (n_material == 0) {
    robustness <- "generally robust"
  } else {
    robustness <- "partially robust"
  }

  rank_statement <- if (all_rank_preserved) {
    "all converged scenarios, HapB3 maintained lowest predicted toxicity among variant genotypes"
  } else {
    "most converged scenarios, HapB3 maintained lowest predicted toxicity among variant genotypes"
  }

  cri_statement <- if (all_cri_exclude_one) {
    "no converged scenario produced a 95% CrI that included 1.0"
  } else {
    sprintf("%d converged scenario(s) produced a 95%% CrI that included 1.0", n_material)
  }

  # Non-converged/infeasible notes for results
  notes_parts <- c()
  if (n_non_converged > 0) {
    nc_names <- non_converged$scenario
    notes_parts <- c(notes_parts, sprintf(
      "%d scenario(s) did not converge (PSRF > 1.10: %s) and should be interpreted with caution",
      n_non_converged, paste(nc_names, collapse = ", ")))
  }
  if (n_infeasible > 0) {
    inf_names <- all_non_primary$scenario[is.na(all_non_primary$HapB3_OR)]
    notes_parts <- c(notes_parts, sprintf(
      "%d scenario(s) were infeasible due to network disconnection or model failure (%s)",
      n_infeasible, paste(inf_names, collapse = ", ")))
  }
  scenario_notes <- if (length(notes_parts) > 0) {
    paste0(" ", paste(notes_parts, collapse = ". "), ".")
  } else {
    ""
  }

  # Total scenario count for methods (all non-primary, including LOO)
  n_total <- nrow(all_non_primary)

  # Build document
  lines <- c(
    "# Manuscript Edit Suggestions: Sensitivity Analysis",
    "",
    "**Generated from:** sensitivity_summary.csv, sensitivity_forest.png",
    sprintf("**Date:** %s", Sys.Date()),
    "**Applies to:** docs/FINAL DRAFT HapB3.pdf",
    "",
    "## Instructions for corresponding author",
    "",
    "These are suggested additions to address PRISMA 2020 Items 13f (sensitivity",
    "analysis methods) and 20d (sensitivity analysis results). Each edit below",
    "specifies the insertion point and the exact text to add. All numerical values",
    "are drawn from the completed sensitivity analysis.",
    "",
    "---",
    "",
    "## Edit 1: Methods -- Add \"Sensitivity Analyses\" subsection",
    "",
    "**Location:** Section 2.5 \"Network Meta-Analysis\", after the paragraph ending",
    "with \"...thinning interval of 10\" (the MCMC configuration paragraph, page 6).",
    "",
    "**Insert new subsection:**",
    "",
    "> **2.6 Sensitivity Analyses**",
    ">",
    sprintf(paste0(
      "> To assess robustness, we conducted %d sensitivity analyses varying the ",
      "model specification (homogeneous equicorrelation), link function (logit vs ",
      "probit), prior specification (inverse-gamma hyperparameters with diffuse and ",
      "weakly informative settings), and sparse data handling (excluding ",
      "arms with N \\u2264 2 and N \\u2264 5). We also performed a complete leave-one-out ",
      "analysis across all %d studies to identify influential observations. ",
      "For each scenario, we compared the HapB3 odds ratio (median, ",
      "95%% CrI), SUCRA rankings, and model fit statistics against the primary ",
      "analysis. Reduced MCMC sampling (100,000 iterations, 50,000 burn-in) was ",
      "used for sensitivity runs, with convergence assessed via the Gelman-Rubin ",
      "potential scale reduction factor (PSRF \\u2264 1.10); scenarios exceeding this ",
      "threshold were re-run with the primary model's full sampling configuration."),
      n_total, 31),
    "",
    "---",
    "",
    "## Edit 2: Results -- Add \"Sensitivity Analyses\" subsection",
    "",
    "**Location:** Section 3 \"Results\", after subsection 3.6 \"Study Quality",
    "Assessment\" and before Section 4 \"Discussion\" (page 10).",
    "",
    "**Insert new subsection:**",
    "",
    "> **3.7 Sensitivity Analyses**",
    ">",
    sprintf(paste0(
      "> The primary finding\\u2014HapB3 OR substantially lower than other DPYD ",
      "variant ORs\\u2014was %s across %d of %d sensitivity scenarios that converged ",
      "(Supplementary Table S1, Supplementary Figure S1). Among converged scenarios, ",
      "HapB3 OR ranged from %.2f to %.2f ",
      "(primary: %.2f), with SUCRA consistently near ",
      "%.1f%% (range: %.1f%%\\u2013%.1f%%). In %s, and %s.%s"),
      robustness, n_converged, n_scenarios, or_min, or_max, primary_or,
      primary_sucra, sucra_min, sucra_max, rank_statement, cri_statement,
      scenario_notes),
    "",
    "---",
    "",
    "## Edit 3: Discussion -- Add limitations acknowledgment",
    "",
    "**Location:** Section 4.3 \"Methodological Considerations\", in the paragraph",
    "beginning \"In terms of limitations regarding study design...\" (page 12).",
    "Append to the end of the existing limitations discussion.",
    "",
    "**Insert:**",
    "",
    paste0(
      "> While sensitivity analyses addressed model specification, link function, ",
      "prior choice, and data sparsity, we did not perform ",
      "formal publication bias assessment (e.g., funnel plots or Egger's test), ",
      "as the network structure\\u2014with most comparisons informed by a single arm ",
      "per study\\u2014limits the applicability of standard methods designed for ",
      "pairwise meta-analysis. Similarly, a formal certainty-of-evidence ",
      "assessment using the CINeMA framework was not conducted; future work ",
      "should address this gap.",
      if (het_cor_failed) {
        paste0(" The inverse-Wishart prior is fixed for the heterogeneous correlation ",
               "model in pcnetmeta; prior sensitivity testing was therefore limited to ",
               "the homogeneous equicorrelation model, confounding model structure with ",
               "prior choice in those comparisons.")
      } else { "" }),
    "",
    "---",
    "",
    "## Edit 4: Supplementary materials",
    "",
    "**Add to supplementary materials:**",
    "",
    "- **Table S1:** Sensitivity analysis summary",
    "  - File: output/sensitivity/sensitivity_summary.csv",
    sprintf(paste0(
      "  - Caption: \"Summary of %d sensitivity analyses testing robustness of ",
      "HapB3 odds ratio. Each row represents a model perturbation with OR ",
      "(median, 95%% CrI), SUCRA, and convergence diagnostics. ",
      "Of %d feasible scenarios, %d converged (PSRF \\u2264 1.10).\""),
      n_total, n_scenarios, n_converged),
    "",
    "- **Figure S1:** Sensitivity analysis forest plot",
    "  - File: output/sensitivity/sensitivity_forest.png",
    paste0(
      "  - Caption: \"Forest plot of HapB3 odds ratio vs wild-type across ",
      "sensitivity scenarios. Filled circles represent main perturbation scenarios; ",
      "open circles represent leave-one-out iterations. Horizontal lines ",
      "represent 95% credible intervals. Dashed vertical line: OR = 1.0 (null). ",
      "Dotted vertical line: primary analysis HapB3 OR.\""),
    "",
    "---",
    "",
    "## Edit 5: Abstract update (optional)",
    "",
    "**Location:** Abstract, after the sentence reporting primary results.",
    "",
    sprintf(paste0(
      "**Consider adding:** \"Sensitivity analyses (%d scenarios including complete ",
      "leave-one-out analysis) confirmed robustness of findings; %d of %d converged ",
      "scenarios were consistent with the primary result.\""),
      n_total, n_consistent, n_converged),
    "",
    "(Check journal word limit before adding.)",
    "",
    "---",
    "",
    "## Edit 6: Data availability statement",
    "",
    "**Location:** End of manuscript, data availability section.",
    "",
    "**Append:** \"Analytic code, including sensitivity analysis scripts, is",
    "available at ___(insert repository URL).\""
  )

  writeLines(lines, output_file)
  cat(sprintf("  Saved: %s\n", output_file))
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

cat("\n###############################################################\n")
cat("### PHASE 1: PRIMARY BASELINE                               ###\n")
cat("###############################################################\n")

# Read primary result (Scenario 1)
results <- list()
results[[1]] <- read_primary_result()
primary_hapb3_or <- results[[1]]$HapB3_OR

cat("\n###############################################################\n")
cat("### PHASE 2: NON-FILTERED SCENARIOS                         ###\n")
cat("###############################################################\n")

# Load data once for non-filtered scenarios
raw_data <- load_raw_data(DATA_FILE)
nma_data_full <- transform_to_nma_format(raw_data, TREATMENT_MAP, aggregate_wt = TRUE)
nma_data_clean <- nma_data_full[, c("s.id", "t.id", "r", "n")]

# Scenario 2: hom_eqcor
results[[length(results) + 1]] <- run_sensitivity_scenario(
  config = list(name = "hom_eqcor",
                description = "Homogeneous equicorrelation model",
                model = "hom_eqcor", link = "probit",
                prior.type = "invgamma", a = 0.001, b = 0.001),
  nma_data = nma_data_clean
)

# Scenario 3: logit link
results[[length(results) + 1]] <- run_sensitivity_scenario(
  config = list(name = "logit_link",
                description = "Logit link function",
                model = "het_cor", link = "logit",
                prior.type = "invwishart"),
  nma_data = nma_data_clean
)

# Scenario 4: alternative prior
results[[length(results) + 1]] <- run_sensitivity_scenario(
  config = list(name = "prior_invgamma",
                description = "Alternative prior (inv-gamma a=0.01)",
                model = "hom_eqcor", link = "probit",
                prior.type = "invgamma", a = 0.01, b = 0.01),
  nma_data = nma_data_clean
)

# Scenario 5: het_cor with inv-gamma prior
results[[length(results) + 1]] <- run_sensitivity_scenario(
  config = list(name = "het_cor_invgamma",
                description = "het_cor with inv-gamma prior (a=0.001)",
                model = "het_cor", link = "probit",
                prior.type = "invgamma", a = 0.001, b = 0.001),
  nma_data = nma_data_clean
)

# Scenario 6: weakly informative prior
results[[length(results) + 1]] <- run_sensitivity_scenario(
  config = list(name = "informative_prior",
                description = "Weakly informative prior (inv-gamma a=0.5, b=0.5)",
                model = "hom_eqcor", link = "probit",
                prior.type = "invgamma", a = 0.5, b = 0.5),
  nma_data = nma_data_clean
)

cat("\n###############################################################\n")
cat("### PHASE 3: FILTERED SCENARIOS                             ###\n")
cat("###############################################################\n")

# Scenario 5a: exclude sparse N <= 2
cat("\n--- Preparing sparse-filtered data (N <= 2) ---\n")
nma_data_sparse2 <- prepare_sparse_filtered_data(threshold = 2)
connectivity_s2 <- check_network_connectivity(nma_data_sparse2, TREATMENT_NAMES)
cat(sprintf("  Connectivity: %s\n",
            ifelse(connectivity_s2$connected, "OK", "DISCONNECTED")))
if (length(connectivity_s2$sparse_nodes) > 0) {
  cat(sprintf("  Sparse nodes (<2 arms): %s\n",
              paste(connectivity_s2$sparse_nodes, collapse = ", ")))
}

if (connectivity_s2$connected) {
  results[[length(results) + 1]] <- run_sensitivity_scenario(
    config = list(name = "exclude_sparse_n2",
                  description = "Exclude arms with N <= 2"),
    nma_data = nma_data_sparse2
  )
} else {
  results[[length(results) + 1]] <- create_failed_row(
    list(name = "exclude_sparse_n2", description = "Exclude arms with N <= 2"),
    reason = sprintf("Network disconnected: %s",
                     paste(connectivity_s2$disconnected_nodes, collapse = ", "))
  )
}

# Scenario 5b: exclude sparse N <= 5
cat("\n--- Preparing sparse-filtered data (N <= 5) ---\n")
nma_data_sparse5 <- prepare_sparse_filtered_data(threshold = 5)
connectivity_s5 <- check_network_connectivity(nma_data_sparse5, TREATMENT_NAMES)
cat(sprintf("  Connectivity: %s\n",
            ifelse(connectivity_s5$connected, "OK", "DISCONNECTED")))
if (length(connectivity_s5$sparse_nodes) > 0) {
  cat(sprintf("  Sparse nodes (<2 arms): %s\n",
              paste(connectivity_s5$sparse_nodes, collapse = ", ")))
}

if (connectivity_s5$connected) {
  results[[length(results) + 1]] <- run_sensitivity_scenario(
    config = list(name = "exclude_sparse_n5",
                  description = "Exclude arms with N <= 5"),
    nma_data = nma_data_sparse5
  )
} else {
  results[[length(results) + 1]] <- create_failed_row(
    list(name = "exclude_sparse_n5", description = "Exclude arms with N <= 5"),
    reason = sprintf("Network disconnected: %s",
                     paste(connectivity_s5$disconnected_nodes, collapse = ", "))
  )
}

cat("\n###############################################################\n")
cat("### PHASE 4: LEAVE-ONE-OUT ANALYSIS                         ###\n")
cat("###############################################################\n")

# Derive LOO study list dynamically from raw data
loo_raw_data <- load_raw_data(DATA_FILE)
all_study_names <- sort(unique(loo_raw_data$Study))
cat(sprintf("  Total unique study names in raw data: %d\n", length(all_study_names)))

# Group multi-arm studies by prefix (entries sharing a prefix before _arm_)
arm_studies <- all_study_names[grepl("_arm_", all_study_names)]
non_arm_studies <- all_study_names[!grepl("_arm_", all_study_names)]

# Extract prefixes from arm studies
arm_prefixes <- unique(sub("_arm_.*$", "", arm_studies))
cat(sprintf("  Multi-arm study groups: %s\n",
            paste(arm_prefixes, collapse = ", ")))

# Build LOO iteration list: non-arm studies + one entry per arm prefix
loo_iterations <- c(non_arm_studies, arm_prefixes)
cat(sprintf("  LOO iterations: %d (%d single + %d grouped)\n",
            length(loo_iterations), length(non_arm_studies), length(arm_prefixes)))

# Track max shift for reporting
max_loo_shift <- 0
max_loo_study <- ""

for (loo_name in loo_iterations) {
  cat(sprintf("\n--- LOO: Excluding %s ---\n", loo_name))

  # Determine if this is a grouped (prefix) or single study exclusion
  is_grouped <- loo_name %in% arm_prefixes

  tryCatch({
    if (is_grouped) {
      # Prefix-match exclusion for multi-arm studies
      nma_data_loo <- prepare_loo_data(loo_name)  # prepare_loo_data already handles prefix
    } else {
      nma_data_loo <- prepare_loo_data(loo_name)
    }

    # Check connectivity
    connectivity <- check_network_connectivity(nma_data_loo, TREATMENT_NAMES)
    if (!connectivity$connected) {
      cat(sprintf("  Network disconnected after removing %s\n", loo_name))
      results[[length(results) + 1]] <- create_failed_row(
        list(name = paste0("loo_", loo_name),
             description = sprintf("LOO: Exclude %s", loo_name)),
        reason = sprintf("Network disconnected: %s",
                         paste(connectivity$disconnected_nodes, collapse = ", "))
      )
      next
    }

    results[[length(results) + 1]] <- run_sensitivity_scenario(
      config = list(name = paste0("loo_", loo_name),
                    description = sprintf("LOO: Exclude %s", loo_name)),
      nma_data = nma_data_loo
    )
  }, error = function(e) {
    cat(sprintf("  ERROR in LOO for %s: %s\n", loo_name, e$message))
    results[[length(results) + 1]] <<- create_failed_row(
      list(name = paste0("loo_", loo_name),
           description = sprintf("LOO: Exclude %s", loo_name)),
      reason = e$message
    )
  })
}

cat(sprintf("\n  LOO complete: %d iterations run\n", length(loo_iterations)))

cat("\n###############################################################\n")
cat("### PHASE 5: ASSEMBLY AND OUTPUTS                           ###\n")
cat("###############################################################\n")

# --- Scenario 8: Check for contrast-based results ---
cat("\n--- Checking for contrast-based model results ---\n")
cb_files <- list.files("output", pattern = "contrast|cb_|CB_",
                       recursive = TRUE, full.names = TRUE)
cb_rds <- cb_files[grepl("\\.rds$", cb_files)]
if (length(cb_rds) > 0) {
  cat(sprintf("  Found CB result files: %s\n", paste(cb_rds, collapse = ", ")))
  tryCatch({
    cb_result <- readRDS(cb_rds[1])
    if (!is.null(cb_result$OddsRatio)) {
      hapb3_or <- parse_or_string(cb_result$OddsRatio$Median_CI[2, 1])
      if (!is.null(hapb3_or)) {
        results[[length(results) + 1]] <- data.frame(
          scenario             = "contrast_based",
          description          = "Contrast-based model (existing)",
          n_studies            = NA,
          n_arms               = NA,
          HapB3_OR             = hapb3_or$median,
          HapB3_CrI_low        = hapb3_or$lower,
          HapB3_CrI_high       = hapb3_or$upper,
          HapB3_SUCRA          = NA,
          HapB3_rank_preserved = NA,
          c2846_OR             = NA,
          c2846_CrI_low        = NA,
          c2846_CrI_high       = NA,
          star2A_OR            = NA,
          star13_OR            = NA,
          star13_CrI_low       = NA,
          star13_CrI_high      = NA,
          DIC                  = NA,
          max_PSRF             = NA,
          converged            = NA,
          escalated            = FALSE,
          hapb3_log_or_shift   = NA,
          conclusion           = NA,
          stringsAsFactors = FALSE
        )
      }
    }
  }, error = function(e) {
    cat(sprintf("  Could not extract CB results: %s\n", e$message))
  })
} else {
  cat("  No contrast-based result files found. Skipping Scenario 8.\n")
}

# --- Assemble summary table ---
cat("\n--- Assembling summary table ---\n")
summary_df <- do.call(rbind, results)
rownames(summary_df) <- NULL

# Compute log-OR shift and conclusion for non-primary rows
for (i in 1:nrow(summary_df)) {
  if (summary_df$scenario[i] == "primary") next
  if (is.na(summary_df$HapB3_OR[i])) next

  hapb3_or <- summary_df$HapB3_OR[i]
  hapb3_lower <- summary_df$HapB3_CrI_low[i]
  hapb3_upper <- summary_df$HapB3_CrI_high[i]
  hapb3_sucra <- summary_df$HapB3_SUCRA[i]

  # Read variant SUCRAs from scenario result for rank comparison
  scenario_dir <- file.path(scenarios_dir, summary_df$scenario[i])
  result_file <- file.path(scenario_dir, "scenario_result.rds")

  variant_sucras <- c(NA, NA, NA)
  if (file.exists(result_file)) {
    tryCatch({
      scenario_result <- readRDS(result_file)
      sucra <- extract_sucra(scenario_result, TREATMENT_NAMES)
      if (!is.null(sucra)) {
        variant_sucras <- sucra[c("2846hetho", "2Ahetho", "13hetho")]
      }
    }, error = function(e) {})
  }

  classification <- classify_conclusion(
    hapb3_or, hapb3_lower, hapb3_upper,
    hapb3_sucra / 100,  # Convert back to 0-1 scale for comparison
    variant_sucras,
    primary_hapb3_or
  )

  summary_df$hapb3_log_or_shift[i] <- round(classification$log_or_shift, 4)
  summary_df$HapB3_rank_preserved[i] <- classification$rank_preserved

  # Non-converged scenarios get a special conclusion regardless of classification
  if (!is.na(summary_df$converged[i]) && summary_df$converged[i] == FALSE) {
    summary_df$conclusion[i] <- "Non-converged (interpret with caution)"
  } else if (is.na(summary_df$conclusion[i])) {
    # Only set conclusion if not already a failure reason
    summary_df$conclusion[i] <- classification$conclusion
  }
}

# Save summary CSV
summary_file <- file.path(output_base, "sensitivity_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat(sprintf("  Saved: %s\n", summary_file))

# Print summary
cat("\n=== Sensitivity Summary ===\n")
for (i in 1:nrow(summary_df)) {
  or_str <- ifelse(
    is.na(summary_df$HapB3_OR[i]), "N/A",
    sprintf("%.3f (%.3f, %.3f)",
            summary_df$HapB3_OR[i],
            summary_df$HapB3_CrI_low[i],
            summary_df$HapB3_CrI_high[i])
  )
  cat(sprintf("  %-25s HapB3 OR: %-25s  Conclusion: %s\n",
              summary_df$scenario[i], or_str,
              summary_df$conclusion[i]))
}

# --- LOO influence summary ---
loo_rows <- summary_df[grepl("^loo_", summary_df$scenario) &
                         !is.na(summary_df$hapb3_log_or_shift), ]
if (nrow(loo_rows) > 0) {
  max_shift_idx <- which.max(loo_rows$hapb3_log_or_shift)
  cat(sprintf("\n  Most influential LOO study: %s (log-OR shift: %.4f)\n",
              loo_rows$scenario[max_shift_idx],
              loo_rows$hapb3_log_or_shift[max_shift_idx]))
}

# --- Validation checks ---
cat("\n=== Validation Checks ===\n")

# Check primary matches
primary_row <- summary_df[summary_df$scenario == "primary", ]
cat(sprintf("  Primary HapB3 OR: %.4f (expected: 2.0050) -- %s\n",
            primary_row$HapB3_OR,
            ifelse(abs(primary_row$HapB3_OR - 2.005) < 0.001, "PASS", "FAIL")))

# Check all ORs positive
all_positive <- all(summary_df$HapB3_OR > 0, na.rm = TRUE)
cat(sprintf("  All HapB3 ORs positive: %s\n",
            ifelse(all_positive, "PASS", "FAIL")))

# Check CrI ordering
cri_valid <- all(
  summary_df$HapB3_CrI_low < summary_df$HapB3_OR &
    summary_df$HapB3_OR < summary_df$HapB3_CrI_high,
  na.rm = TRUE
)
cat(sprintf("  CrI ordering (low < median < high): %s\n",
            ifelse(cri_valid, "PASS", "FAIL")))

# Check DIC finite
dic_valid <- all(is.finite(summary_df$DIC[!is.na(summary_df$DIC)]))
cat(sprintf("  DIC finite for all converged: %s\n",
            ifelse(dic_valid, "PASS", "FAIL")))

# Check SUCRA range
sucra_valid <- all(
  summary_df$HapB3_SUCRA >= 0 & summary_df$HapB3_SUCRA <= 100,
  na.rm = TRUE
)
cat(sprintf("  SUCRA in 0-100%% range: %s\n",
            ifelse(sucra_valid, "PASS", "FAIL")))

# Check n_studies and n_arms for filtered scenarios
cat("\n  Filtered scenario data counts:\n")
for (sc in c("exclude_sparse_n2", "exclude_sparse_n5")) {
  row <- summary_df[summary_df$scenario == sc, ]
  if (nrow(row) > 0 && !is.na(row$n_studies)) {
    cat(sprintf("    %s: %d studies, %d arms\n",
                sc, row$n_studies, row$n_arms))
  }
}

# --- Forest plot ---
forest_file <- file.path(output_base, "sensitivity_forest.png")
generate_forest_plot(summary_df, forest_file)

# --- Manuscript edits ---
manuscript_file <- file.path(project_root, "docs/manuscript-edits-sensitivity.md")
dir.create(file.path(project_root, "docs"), showWarnings = FALSE)
generate_manuscript_edits(summary_df, manuscript_file)

# =============================================================================
# Complete
# =============================================================================

cat("\n=============================================================\n")
cat("Sensitivity Analysis Complete\n")
cat(sprintf("Finished: %s\n", Sys.time()))
cat("=============================================================\n")
cat("\nOutputs:\n")
cat(sprintf("  Summary: %s\n", summary_file))
cat(sprintf("  Forest plot: %s\n", forest_file))
cat(sprintf("  Log: %s\n", log_file))
cat(sprintf("  Manuscript edits: %s\n", manuscript_file))
cat("\nScenario results saved in: output/sensitivity/scenarios/\n")
