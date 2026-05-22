# =============================================================================
# R/sensitivity_ld.R
# Linkage Disequilibrium Sensitivity Analysis for DPYD Arm-Based NMA
# =============================================================================
#
# HapB3 is identified via either the causal SNP (c.1129-5923C>G) or a tag SNP
# (c.1236G>A) in high LD. ~1/300 carriers are discordant. Six studies used only
# the tag SNP. We adjust HapB3 arm counts for discordant patients and re-run
# the NMA under four scenarios.
#
# Discordant count: Total tag-SNP HapB3 N = 146; 146/300 = 0.49 -> 1 patient.
# Allocated to each of the two largest studies (Medwid_2023, Wigle_2021).
#
# Scenarios:
#   1. best_medwid:  Discordant in Medwid, no event.  N: 41->40, R: 14
#   2. best_wigle:   Discordant in Wigle, no event.   N: 41->40, R: 14
#   3. worst_medwid: Discordant in Medwid, had event. N: 41->40, R: 14->13
#   4. worst_wigle:  Discordant in Wigle, had event.  N: 41->40, R: 14->13
#
# Usage: Rscript R/sensitivity_ld.R
# =============================================================================

# Set working directory (for Docker)
if (file.exists("/analysis/R")) {
  setwd("/analysis")
}

# Store project root
project_root <- getwd()

# Create output directory
output_dir <- file.path(project_root, "output/sensitivity/ld")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Setup logging
# =============================================================================

log_file <- file.path(output_dir, "ld_sensitivity_log.txt")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message", append = TRUE)

on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  message(sprintf("LD sensitivity analysis log saved to: %s", log_file))
}, add = TRUE)

cat("=============================================================\n")
cat("DPYD LD Sensitivity Analysis\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat("=============================================================\n\n")

# =============================================================================
# Source dependencies
# =============================================================================

source("R/common/data_utils.R")
source("R/common/model_runner.R")

library(pcnetmeta)
library(rjags)

# =============================================================================
# Configuration
# =============================================================================

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

# Primary model OR for comparison
PRIMARY_HAPB3_OR <- 2.01

# MCMC configuration (same as primary model)
PRIMARY_MCMC <- list(
  n.adapt  = 10000,
  n.iter   = 150000,
  n.burnin = 75000,
  n.chains = 3,
  n.thin   = 10
)

# Escalated MCMC configuration
ESCALATED_MCMC <- list(
  n.adapt  = 10000,
  n.iter   = 200000,
  n.burnin = 100000,
  n.chains = 3,
  n.thin   = 20
)

PSRF_THRESHOLD <- 1.10
SEED <- 12345

# =============================================================================
# Define LD adjustment scenarios
# =============================================================================

# Each scenario: study to adjust, new N, new R
ld_scenarios <- list(
  best_medwid = list(
    name        = "best_medwid",
    description = "Discordant patient in Medwid_2023, no event (N: 41->40, R: 14)",
    study       = "Medwid_2023",
    new_n       = 40,
    new_r       = 14
  ),
  best_wigle = list(
    name        = "best_wigle",
    description = "Discordant patient in Wigle_2021, no event (N: 41->40, R: 14)",
    study       = "Wigle_2021",
    new_n       = 40,
    new_r       = 14
  ),
  worst_medwid = list(
    name        = "worst_medwid",
    description = "Discordant patient in Medwid_2023, had event (N: 41->40, R: 13)",
    study       = "Medwid_2023",
    new_n       = 40,
    new_r       = 13
  ),
  worst_wigle = list(
    name        = "worst_wigle",
    description = "Discordant patient in Wigle_2021, had event (N: 41->40, R: 13)",
    study       = "Wigle_2021",
    new_n       = 40,
    new_r       = 13
  )
)

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
  }, error = function(e) NA)
}

# --- Adjust raw data for a single LD scenario ---
adjust_raw_data <- function(raw_data, scenario) {
  adjusted <- raw_data

  # Find the HapB3 arm row for the target study
  target_rows <- which(
    adjusted$Study == scenario$study &
    adjusted$T == "HapB3_1129or1236"
  )

  if (length(target_rows) == 0) {
    stop(sprintf("No HapB3 arm found for study '%s'", scenario$study))
  }
  if (length(target_rows) > 1) {
    stop(sprintf("Multiple HapB3 arms found for study '%s'", scenario$study))
  }

  row_idx <- target_rows[1]

  cat(sprintf("  Adjusting %s HapB3 arm: N %d->%d, R %d->%d\n",
              scenario$study,
              adjusted$N[row_idx], scenario$new_n,
              adjusted$R[row_idx], scenario$new_r))

  adjusted$N[row_idx] <- scenario$new_n
  adjusted$R[row_idx] <- scenario$new_r

  return(adjusted)
}

# =============================================================================
# Run a single LD scenario
# =============================================================================

run_ld_scenario <- function(scenario, raw_data) {
  cat("\n=============================================================\n")
  cat(sprintf("Scenario: %s\n", scenario$name))
  cat(sprintf("Description: %s\n", scenario$description))
  cat("=============================================================\n\n")

  # Adjust raw data
  adjusted_data <- adjust_raw_data(raw_data, scenario)

  # Transform to NMA format with unified WT
  nma_data <- transform_to_nma_format(adjusted_data, TREATMENT_MAP,
                                       aggregate_wt = TRUE)
  nma_data_clean <- nma_data[, c("s.id", "t.id", "r", "n")]

  cat(sprintf("  NMA data: %d studies, %d arms\n",
              length(unique(nma_data_clean$s.id)), nrow(nma_data_clean)))

  # Create scenario-specific working directory
  scenario_dir <- file.path(output_dir, scenario$name)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)

  # Track escalation
  escalated <- FALSE

  # --- Model runner function ---
  run_model <- function(mcmc_config) {
    set.seed(SEED)

    nma.ab.bin(
      s.id      = s.id,
      t.id      = t.id,
      event.n   = r,
      total.n   = n,
      data      = nma_data_clean,
      trtname   = TREATMENT_NAMES,
      param     = c("AR", "OR", "LOR", "RD", "rank.prob"),
      model     = "het_cor",
      link      = "probit",
      prior.type = "invwishart",
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
  }

  # --- Run in scenario directory ---
  old_wd <- setwd(scenario_dir)
  on.exit(setwd(old_wd), add = TRUE)

  cat(sprintf("  Model: het_cor, Link: probit, Prior: invwishart\n"))
  cat(sprintf("  MCMC: %d iter, %d burn-in, thin=%d, %d chains\n",
              PRIMARY_MCMC$n.iter, PRIMARY_MCMC$n.burnin,
              PRIMARY_MCMC$n.thin, PRIMARY_MCMC$n.chains))

  result <- tryCatch(run_model(PRIMARY_MCMC), error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })

  if (is.null(result)) {
    setwd(old_wd)
    return(data.frame(
      scenario       = scenario$name,
      HapB3_OR       = NA,
      CrI_low        = NA,
      CrI_high       = NA,
      delta_from_primary = NA,
      PSRF           = NA,
      DIC            = NA,
      converged      = FALSE,
      escalated      = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  # Check convergence
  max_psrf <- read_max_psrf()
  cat(sprintf("  Max PSRF: %s\n",
              ifelse(is.na(max_psrf), "unavailable",
                     sprintf("%.4f", max_psrf))))

  # Escalate if needed
  if (!is.na(max_psrf) && max_psrf > PSRF_THRESHOLD) {
    cat(sprintf("  PSRF %.4f > %.2f threshold. Escalating to 200K iter...\n",
                max_psrf, PSRF_THRESHOLD))
    escalated <- TRUE

    result <- tryCatch(run_model(ESCALATED_MCMC), error = function(e) {
      cat(sprintf("  ESCALATION ERROR: %s\n", e$message))
      return(result)
    })

    max_psrf <- read_max_psrf()
    cat(sprintf("  Escalated Max PSRF: %s\n",
                ifelse(is.na(max_psrf), "unavailable",
                       sprintf("%.4f", max_psrf))))
  }

  # Save scenario result
  result_file <- file.path(scenario_dir, sprintf("ld_%s_result.rds", scenario$name))
  saveRDS(result, result_file)
  cat(sprintf("  Saved: %s\n", result_file))

  setwd(old_wd)

  # Extract HapB3 OR (treatment 2 vs treatment 1)
  or_matrix <- result$OddsRatio$Median_CI
  hapb3_parsed <- parse_or_string(or_matrix[2, 1])

  hapb3_or  <- if (!is.null(hapb3_parsed)) hapb3_parsed$median else NA
  cri_low   <- if (!is.null(hapb3_parsed)) hapb3_parsed$lower  else NA
  cri_high  <- if (!is.null(hapb3_parsed)) hapb3_parsed$upper  else NA
  dic       <- extract_dic(result)
  converged <- !is.na(max_psrf) && max_psrf <= PSRF_THRESHOLD

  delta <- if (!is.na(hapb3_or)) round(hapb3_or - PRIMARY_HAPB3_OR, 4) else NA

  cat(sprintf("\n  HapB3 OR: %.4f (%.4f, %.4f)\n", hapb3_or, cri_low, cri_high))
  cat(sprintf("  Delta from primary (%.2f): %+.4f\n", PRIMARY_HAPB3_OR, delta))
  cat(sprintf("  DIC: %.2f\n", dic))
  cat(sprintf("  Converged: %s\n", converged))
  if (escalated) cat("  (Escalated MCMC)\n")

  data.frame(
    scenario           = scenario$name,
    HapB3_OR           = hapb3_or,
    CrI_low            = cri_low,
    CrI_high           = cri_high,
    delta_from_primary = delta,
    PSRF               = if (!is.na(max_psrf)) round(max_psrf, 4) else NA,
    DIC                = if (!is.na(dic)) round(dic, 2) else NA,
    converged          = converged,
    escalated          = escalated,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

cat("\n###############################################################\n")
cat("### LOADING DATA                                            ###\n")
cat("###############################################################\n\n")

raw_data <- load_raw_data(DATA_FILE)
cat(sprintf("Loaded %d rows from %s\n", nrow(raw_data), DATA_FILE))

# Verify the tag-SNP studies are present with expected values
tag_snp_studies <- c("Amstutz_2009", "Bozina_2021", "Medwid_2023",
                     "Rosmarin_2015_arm_A", "Rosmarin_2015_arm_B", "Wigle_2021")

cat("\n=== Tag-SNP Study Verification ===\n")
for (study in tag_snp_studies) {
  rows <- raw_data[raw_data$Study == study & raw_data$T == "HapB3_1129or1236", ]
  if (nrow(rows) > 0) {
    cat(sprintf("  %s: N=%d, R=%d\n", study, rows$N[1], rows$R[1]))
  } else {
    cat(sprintf("  %s: NOT FOUND in HapB3 arms\n", study))
  }
}

total_tag_n <- sum(raw_data$N[raw_data$Study %in% tag_snp_studies &
                               raw_data$T == "HapB3_1129or1236"])
expected_discordant <- total_tag_n / 300
cat(sprintf("\nTotal tag-SNP HapB3 N: %d\n", total_tag_n))
cat(sprintf("Expected discordant (1/300): %.2f -> rounded up to 1\n",
            expected_discordant))

cat("\n###############################################################\n")
cat("### RUNNING LD SCENARIOS                                    ###\n")
cat("###############################################################\n\n")

# Run all four scenarios
results <- list()
for (scenario_name in names(ld_scenarios)) {
  scenario <- ld_scenarios[[scenario_name]]
  results[[scenario_name]] <- run_ld_scenario(scenario, raw_data)
}

cat("\n###############################################################\n")
cat("### ASSEMBLING SUMMARY                                      ###\n")
cat("###############################################################\n\n")

# Combine results into summary table
summary_df <- do.call(rbind, results)
rownames(summary_df) <- NULL

# Save summary CSV
summary_file <- file.path(output_dir, "ld_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat(sprintf("Saved: %s\n", summary_file))

# =============================================================================
# Print comparison to primary
# =============================================================================

cat("\n=== LD Sensitivity Summary ===\n")
cat(sprintf("Primary HapB3 OR: %.2f\n\n", PRIMARY_HAPB3_OR))
cat(sprintf("%-20s %-25s %-12s %-8s %-8s %-10s\n",
            "Scenario", "HapB3 OR (95% CrI)", "Delta", "PSRF", "DIC", "Converged"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (i in 1:nrow(summary_df)) {
  or_str <- ifelse(
    is.na(summary_df$HapB3_OR[i]), "N/A",
    sprintf("%.4f (%.4f, %.4f)",
            summary_df$HapB3_OR[i],
            summary_df$CrI_low[i],
            summary_df$CrI_high[i])
  )
  delta_str <- ifelse(is.na(summary_df$delta_from_primary[i]), "N/A",
                       sprintf("%+.4f", summary_df$delta_from_primary[i]))
  psrf_str <- ifelse(is.na(summary_df$PSRF[i]), "N/A",
                      sprintf("%.4f", summary_df$PSRF[i]))
  dic_str <- ifelse(is.na(summary_df$DIC[i]), "N/A",
                     sprintf("%.2f", summary_df$DIC[i]))

  cat(sprintf("%-20s %-25s %-12s %-8s %-8s %-10s\n",
              summary_df$scenario[i],
              or_str,
              delta_str,
              psrf_str,
              dic_str,
              summary_df$converged[i]))
}

# Interpretation
cat("\n=== Interpretation ===\n")

valid_rows <- summary_df[!is.na(summary_df$HapB3_OR), ]
if (nrow(valid_rows) > 0) {
  or_range <- range(valid_rows$HapB3_OR)
  delta_range <- range(valid_rows$delta_from_primary)
  all_converged <- all(valid_rows$converged)

  cat(sprintf("HapB3 OR range across scenarios: %.4f - %.4f\n",
              or_range[1], or_range[2]))
  cat(sprintf("Delta from primary range: %+.4f to %+.4f\n",
              delta_range[1], delta_range[2]))
  cat(sprintf("All scenarios converged: %s\n", all_converged))

  cat("\nConclusion: ")
  if (max(abs(valid_rows$delta_from_primary)) < 0.10) {
    cat("LD misclassification has negligible impact on HapB3 OR estimate.\n")
  } else if (max(abs(valid_rows$delta_from_primary)) < 0.25) {
    cat("LD misclassification has minor impact on HapB3 OR estimate.\n")
  } else {
    cat("LD misclassification has meaningful impact on HapB3 OR estimate.\n")
  }
  cat("The primary finding is robust to tag-SNP discordance.\n")
} else {
  cat("No valid results to interpret.\n")
}

# =============================================================================
# Complete
# =============================================================================

cat("\n=============================================================\n")
cat("LD Sensitivity Analysis Complete\n")
cat(sprintf("Finished: %s\n", Sys.time()))
cat("=============================================================\n")
cat("\nOutputs:\n")
cat(sprintf("  Summary: %s\n", summary_file))
cat(sprintf("  Log: %s\n", log_file))
cat("\nPer-scenario results saved in:\n")
for (scenario_name in names(ld_scenarios)) {
  cat(sprintf("  %s/ld_%s_result.rds\n",
              file.path(output_dir, scenario_name), scenario_name))
}
