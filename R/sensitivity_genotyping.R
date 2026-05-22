# =============================================================================
# R/sensitivity_genotyping.R
# Genotyping Misclassification Sensitivity Analysis
# =============================================================================
#
# Splits WT into WT_Clean (comprehensive genotyping) and WT_Biased (incomplete
# panels) to test whether genotyping quality affects the primary finding.
#
# Primary model uses 5 treatments (WT unified). This analysis uses 6:
#   1 = WT_Clean, 2 = WT_Biased, 3 = HapB3, 4 = 2846hetho, 5 = 2Ahetho, 6 = 13hetho
#
# Key questions:
#   - Does HapB3 vs WT_Clean OR differ from HapB3 vs WT_Biased OR?
#   - Does WT_Clean vs WT_Biased show a meaningful difference?
#   - Are SUCRA rankings stable with WT split?
#
# Usage: Rscript R/sensitivity_genotyping.R
# =============================================================================

# Set working directory (for Docker)
if (file.exists("/analysis/R")) {
  setwd("/analysis")
}

# Store project root
project_root <- getwd()

# =============================================================================
# Parse CLI args: --definition=loose|strict
# =============================================================================
#
# loose  (default, n=9 WT_Clean): uses the data file's encoded classification.
#   Includes Amstutz_2009, Jennings_2013, Froehlich_2015 -- studies with
#   functionally pan-negative WT controls despite pairwise/grouped designs.
#
# strict (n=6 WT_Clean): reclassifies the 3 borderline studies above to
#   WT_Biased, retaining only studies with pre-screening exclusion of variant
#   carriers (pan-negative by design).

args <- commandArgs(trailingOnly = TRUE)
WT_DEFINITION <- "loose"
for (arg in args) {
  if (grepl("^--definition=", arg)) {
    WT_DEFINITION <- sub("^--definition=", "", arg)
  } else if (arg %in% c("strict", "loose")) {
    WT_DEFINITION <- arg
  }
}
if (!WT_DEFINITION %in% c("strict", "loose")) {
  stop("Invalid --definition: must be 'strict' or 'loose'")
}

# =============================================================================
# Output directory (per definition)
# =============================================================================

output_dir <- if (WT_DEFINITION == "strict") {
  file.path(project_root, "output/sensitivity/genotyping/strict")
} else {
  file.path(project_root, "output/sensitivity/genotyping")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Setup logging
# =============================================================================

log_file <- file.path(output_dir, "genotyping_log.txt")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message", append = TRUE)

on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  message(sprintf("Genotyping sensitivity log saved to: %s", log_file))
}, add = TRUE)

cat("=============================================================\n")
cat("DPYD Genotyping Misclassification Sensitivity Analysis\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat(sprintf("WT_Clean definition: %s\n", WT_DEFINITION))
cat("=============================================================\n\n")

# =============================================================================
# Source dependencies
# =============================================================================

source("R/common/data_utils.R")
source("R/common/model_runner.R")  # for check_convergence(), read_max_psrf pattern

library(pcnetmeta)
library(rjags)

# =============================================================================
# Configuration
# =============================================================================

# Treatment mapping (6 treatments, WT split)
TREATMENT_MAP <- c(
  "WT_Clean"         = 1,
  "WT_Biased"        = 2,
  "HapB3_1129or1236" = 3,
  "2846hetho"        = 4,
  "2Ahetho"          = 5,
  "13hetho"          = 6
)

TREATMENT_NAMES <- c("WT_Clean", "WT_Biased", "HapB3", "2846hetho", "2Ahetho", "13hetho")

DATA_FILE <- "data/Binary WT HapB3 Data for NMA (12-18-2025).xlsx"

# MCMC configuration
MCMC_CONFIG <- list(
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

# Convergence threshold
PSRF_THRESHOLD <- 1.10

# Seed for reproducibility
SEED <- 12345

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

# =============================================================================
# PHASE 1: Load and prepare data
# =============================================================================

cat("\n###############################################################\n")
cat("### PHASE 1: DATA PREPARATION                               ###\n")
cat("###############################################################\n\n")

raw_data <- load_raw_data(DATA_FILE)
cat("Raw data loaded:", nrow(raw_data), "rows\n\n")

# Apply strict-definition reclassification if requested.
# These 3 studies have functionally pan-negative controls but use pairwise or
# grouped (carriers vs non-carriers) designs rather than pre-screening
# exclusion. Phase2 audit reclassified them as WT_Clean for the loose
# definition; the strict definition reverts them to WT_Biased.
if (WT_DEFINITION == "strict") {
  borderline_studies <- c("Amstutz_2009", "Jennings_2013", "Froehlich_2015")
  rows_to_reclass <- raw_data$Study %in% borderline_studies & raw_data$T == "WT_Clean"
  cat(sprintf("STRICT definition: reclassifying %d arms from WT_Clean to WT_Biased\n",
              sum(rows_to_reclass)))
  for (study in borderline_studies) {
    n <- sum(raw_data$Study == study & raw_data$T == "WT_Clean")
    cat(sprintf("  %s: %d WT_Clean arm(s) -> WT_Biased\n", study, n))
  }
  raw_data$T[rows_to_reclass] <- "WT_Biased"
  cat("\n")
}

# Transform with aggregate_wt = FALSE to preserve WT_Clean / WT_Biased split
nma_data <- transform_to_nma_format(raw_data, TREATMENT_MAP, aggregate_wt = FALSE)

# Validate
validate_nma_data(nma_data)

# Print treatment summary
print_treatment_summary(nma_data, TREATMENT_NAMES)

# Report data quality
report_data_quality(nma_data, TREATMENT_NAMES)

# Report sparse data
report_sparse_data(nma_data)

# Keep only pcnetmeta columns
nma_data_clean <- nma_data %>%
  select(s.id, t.id, r, n)

cat("\n=== WT Split Summary ===\n")
wt_clean_studies <- nma_data %>% filter(t.id == 1) %>% pull(s.id) %>% unique()
wt_biased_studies <- nma_data %>% filter(t.id == 2) %>% pull(s.id) %>% unique()
cat("Studies with WT_Clean:", length(wt_clean_studies), "\n")
cat("Studies with WT_Biased:", length(wt_biased_studies), "\n")
cat("Studies with both:", length(intersect(wt_clean_studies, wt_biased_studies)), "\n\n")

# =============================================================================
# PHASE 2: Network connectivity check
# =============================================================================

cat("\n###############################################################\n")
cat("### PHASE 2: NETWORK CONNECTIVITY                           ###\n")
cat("###############################################################\n\n")

connectivity <- check_network_connectivity(nma_data_clean, TREATMENT_NAMES)

cat("Network connected:", connectivity$connected, "\n")
cat("Arms per treatment:\n")
for (trt_id in seq_along(TREATMENT_NAMES)) {
  n_arms <- if (as.character(trt_id) %in% names(connectivity$arms_per_treatment)) {
    connectivity$arms_per_treatment[[as.character(trt_id)]]
  } else {
    0
  }
  cat(sprintf("  %s (t.id=%d): %d arms\n", TREATMENT_NAMES[trt_id], trt_id, n_arms))
}

if (length(connectivity$disconnected_nodes) > 0) {
  cat("\nDISCONNECTED treatments:", paste(connectivity$disconnected_nodes, collapse = ", "), "\n")
  stop("Network is disconnected. Cannot run NMA with split WT.")
}

if (length(connectivity$sparse_nodes) > 0) {
  cat("\nWARNING: Sparse treatments (<2 arms):",
      paste(connectivity$sparse_nodes, collapse = ", "), "\n")
  cat("  Model may have difficulty estimating these treatments.\n")
}

cat("\nNetwork connectivity: OK\n")

# =============================================================================
# PHASE 3: Run the model
# =============================================================================

cat("\n###############################################################\n")
cat("### PHASE 3: MODEL FITTING                                  ###\n")
cat("###############################################################\n\n")

escalated <- FALSE

run_genotyping_model <- function(mcmc_config) {
  set.seed(SEED)

  cat(sprintf("Model: het_cor, Link: probit, Prior: invwishart\n"))
  cat(sprintf("MCMC: %d iter, %d burn-in, thin=%d, %d chains\n",
              mcmc_config$n.iter, mcmc_config$n.burnin,
              mcmc_config$n.thin, mcmc_config$n.chains))
  cat(sprintf("Data: %d studies, %d arms, %d treatments\n\n",
              length(unique(nma_data_clean$s.id)),
              nrow(nma_data_clean),
              length(TREATMENT_NAMES)))

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

# Run in output directory so pcnetmeta diagnostic files are saved there
old_wd <- setwd(output_dir)

result <- tryCatch(run_genotyping_model(MCMC_CONFIG), error = function(e) {
  cat(sprintf("ERROR: %s\n", e$message))
  return(NULL)
})

if (is.null(result)) {
  setwd(old_wd)
  stop("Model failed to run. Check data and JAGS installation.")
}

# =============================================================================
# PHASE 4: Convergence check and escalation
# =============================================================================

cat("\n###############################################################\n")
cat("### PHASE 4: CONVERGENCE DIAGNOSTICS                        ###\n")
cat("###############################################################\n\n")

max_psrf <- read_max_psrf()
cat(sprintf("Max PSRF: %s\n",
            ifelse(is.na(max_psrf), "unavailable", sprintf("%.4f", max_psrf))))
cat(sprintf("Threshold: %.2f\n", PSRF_THRESHOLD))

if (!is.na(max_psrf) && max_psrf > PSRF_THRESHOLD) {
  cat(sprintf("\nPSRF %.4f > %.2f threshold. Escalating MCMC...\n",
              max_psrf, PSRF_THRESHOLD))
  cat(sprintf("Escalating to %d iterations, %d burn-in, thin=%d\n\n",
              ESCALATED_MCMC$n.iter, ESCALATED_MCMC$n.burnin, ESCALATED_MCMC$n.thin))
  escalated <- TRUE

  result_escalated <- tryCatch(run_genotyping_model(ESCALATED_MCMC), error = function(e) {
    cat(sprintf("ESCALATION ERROR: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(result_escalated)) {
    result <- result_escalated
    max_psrf <- read_max_psrf()
    cat(sprintf("Escalated Max PSRF: %s\n",
                ifelse(is.na(max_psrf), "unavailable", sprintf("%.4f", max_psrf))))
  } else {
    cat("Escalation failed. Using initial result.\n")
  }
}

converged <- !is.na(max_psrf) && max_psrf <= PSRF_THRESHOLD
cat(sprintf("\nConvergence status: %s\n",
            ifelse(converged, "CONVERGED", "NOT CONVERGED")))

setwd(old_wd)

# =============================================================================
# PHASE 5: Extract results
# =============================================================================

cat("\n###############################################################\n")
cat("### PHASE 5: RESULTS EXTRACTION                             ###\n")
cat("###############################################################\n\n")

# --- Full OR matrix ---
cat("=== Full Odds Ratio Matrix ===\n")
or_matrix <- result$OddsRatio$Median_CI
print(or_matrix)

# --- Parse key pairwise ORs ---
# OR matrix from pcnetmeta: row i vs column j (row = treatment, col = reference)
# Column 1 = vs treatment 1 (WT_Clean)
# Column 2 = vs treatment 2 (WT_Biased), etc.

cat("\n=== Key Pairwise Comparisons ===\n\n")

# HapB3 vs WT_Clean (treatment 3 vs 1) -- row 3, col 1
hapb3_vs_wtclean <- parse_or_string(or_matrix[3, 1])
if (!is.null(hapb3_vs_wtclean)) {
  cat(sprintf("HapB3 vs WT_Clean:   OR = %.4f (%.4f, %.4f)\n",
              hapb3_vs_wtclean$median, hapb3_vs_wtclean$lower, hapb3_vs_wtclean$upper))
  cat(sprintf("  95%% CrI excludes 1.0: %s\n",
              ifelse(hapb3_vs_wtclean$lower > 1.0, "YES (significant)", "NO")))
}

# HapB3 vs WT_Biased (treatment 3 vs 2) -- row 3, col 2
hapb3_vs_wtbiased <- parse_or_string(or_matrix[3, 2])
if (!is.null(hapb3_vs_wtbiased)) {
  cat(sprintf("\nHapB3 vs WT_Biased:  OR = %.4f (%.4f, %.4f)\n",
              hapb3_vs_wtbiased$median, hapb3_vs_wtbiased$lower, hapb3_vs_wtbiased$upper))
  cat(sprintf("  95%% CrI excludes 1.0: %s\n",
              ifelse(hapb3_vs_wtbiased$lower > 1.0, "YES (significant)", "NO")))
}

# WT_Clean vs WT_Biased (treatment 1 vs 2) -- diagnostic
# Row 1 vs col 2, but row 1 col 2 may be WT_Clean vs WT_Biased or vice versa
# pcnetmeta OR matrix: or_matrix[i, j] = OR of treatment i vs treatment j
wtclean_vs_wtbiased <- parse_or_string(or_matrix[1, 2])
if (!is.null(wtclean_vs_wtbiased)) {
  cat(sprintf("\nWT_Clean vs WT_Biased (diagnostic):  OR = %.4f (%.4f, %.4f)\n",
              wtclean_vs_wtbiased$median, wtclean_vs_wtbiased$lower, wtclean_vs_wtbiased$upper))
  cat(sprintf("  95%% CrI excludes 1.0: %s\n",
              ifelse(wtclean_vs_wtbiased$lower > 1.0 || wtclean_vs_wtbiased$upper < 1.0,
                     "YES (WT groups differ)", "NO (WT groups similar)")))
} else {
  # Try the reverse: WT_Biased vs WT_Clean
  wtbiased_vs_wtclean <- parse_or_string(or_matrix[2, 1])
  if (!is.null(wtbiased_vs_wtclean)) {
    cat(sprintf("\nWT_Biased vs WT_Clean (diagnostic):  OR = %.4f (%.4f, %.4f)\n",
                wtbiased_vs_wtclean$median, wtbiased_vs_wtclean$lower, wtbiased_vs_wtclean$upper))
    cat(sprintf("  95%% CrI excludes 1.0: %s\n",
                ifelse(wtbiased_vs_wtclean$lower > 1.0 || wtbiased_vs_wtclean$upper < 1.0,
                       "YES (WT groups differ)", "NO (WT groups similar)")))
  }
}

# Other variant ORs vs WT_Clean
cat("\n=== All Variant ORs vs WT_Clean ===\n")
for (i in 1:length(TREATMENT_NAMES)) {
  parsed <- parse_or_string(or_matrix[i, 1])
  if (!is.null(parsed)) {
    cat(sprintf("  %-12s: OR = %.4f (%.4f, %.4f)\n",
                TREATMENT_NAMES[i], parsed$median, parsed$lower, parsed$upper))
  } else {
    cat(sprintf("  %-12s: %s\n", TREATMENT_NAMES[i], or_matrix[i, 1]))
  }
}

# --- SUCRA Rankings ---
cat("\n=== SUCRA Rankings ===\n")
cat("(Higher SUCRA = better = lower toxicity)\n\n")
sucra <- extract_sucra(result, TREATMENT_NAMES)
if (!is.null(sucra)) {
  sucra_sorted <- sort(sucra, decreasing = TRUE)
  for (i in seq_along(sucra_sorted)) {
    cat(sprintf("  %d. %-12s: %.1f%%\n", i, names(sucra_sorted)[i], sucra_sorted[i] * 100))
  }
} else {
  cat("  SUCRA could not be calculated.\n")
}

# --- Rank Probabilities ---
cat("\n=== Rank Probability Matrix ===\n")
if (!is.null(result$TrtRankProb)) {
  print(result$TrtRankProb)
}

# --- Absolute Risk ---
cat("\n=== Absolute Risk Estimates ===\n")
if (!is.null(result$AbsoluteRisk)) {
  print(result$AbsoluteRisk$Median_CI)
}

# --- DIC ---
dic <- extract_dic(result)
cat(sprintf("\n=== Model Fit ===\n"))
cat(sprintf("DIC: %s\n", ifelse(is.na(dic), "unavailable", sprintf("%.2f", dic))))

if (!is.null(result$DIC)) {
  tryCatch({
    dic_vec <- result$DIC
    if (is.matrix(dic_vec) || is.array(dic_vec)) {
      cat(sprintf("D.bar: %.2f\n", dic_vec["D.bar", 1]))
      cat(sprintf("pD: %.2f\n", dic_vec["pD", 1]))
    }
  }, error = function(e) {})
}

# =============================================================================
# PHASE 6: Save outputs
# =============================================================================

cat("\n###############################################################\n")
cat("### PHASE 6: SAVING OUTPUTS                                 ###\n")
cat("###############################################################\n\n")

# Save full result
saveRDS(result, file.path(output_dir, "genotyping_result.rds"))
cat("Saved:", file.path(output_dir, "genotyping_result.rds"), "\n")

# Save convergence status
convergence_info <- list(
  converged = converged,
  max_psrf = max_psrf,
  psrf_threshold = PSRF_THRESHOLD,
  escalated = escalated,
  mcmc_config = if (escalated) ESCALATED_MCMC else MCMC_CONFIG
)
saveRDS(convergence_info, file.path(output_dir, "genotyping_convergence.rds"))
cat("Saved:", file.path(output_dir, "genotyping_convergence.rds"), "\n")

# Build summary CSV
summary_rows <- list()

# All pairwise ORs vs WT_Clean (column 1)
for (i in 1:length(TREATMENT_NAMES)) {
  parsed <- parse_or_string(or_matrix[i, 1])
  if (!is.null(parsed)) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      comparison = paste0(TREATMENT_NAMES[i], " vs WT_Clean"),
      OR_median = parsed$median,
      OR_lower = parsed$lower,
      OR_upper = parsed$upper,
      CrI_excludes_1 = (parsed$lower > 1.0 || parsed$upper < 1.0),
      SUCRA = if (!is.null(sucra)) round(sucra[TREATMENT_NAMES[i]] * 100, 1) else NA,
      stringsAsFactors = FALSE
    )
  }
}

# HapB3 vs WT_Biased
if (!is.null(hapb3_vs_wtbiased)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    comparison = "HapB3 vs WT_Biased",
    OR_median = hapb3_vs_wtbiased$median,
    OR_lower = hapb3_vs_wtbiased$lower,
    OR_upper = hapb3_vs_wtbiased$upper,
    CrI_excludes_1 = (hapb3_vs_wtbiased$lower > 1.0 || hapb3_vs_wtbiased$upper < 1.0),
    SUCRA = NA,
    stringsAsFactors = FALSE
  )
}

# WT_Clean vs WT_Biased diagnostic
if (!is.null(wtclean_vs_wtbiased)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    comparison = "WT_Clean vs WT_Biased (diagnostic)",
    OR_median = wtclean_vs_wtbiased$median,
    OR_lower = wtclean_vs_wtbiased$lower,
    OR_upper = wtclean_vs_wtbiased$upper,
    CrI_excludes_1 = (wtclean_vs_wtbiased$lower > 1.0 || wtclean_vs_wtbiased$upper < 1.0),
    SUCRA = NA,
    stringsAsFactors = FALSE
  )
}

# Add model diagnostics row
summary_rows[[length(summary_rows) + 1]] <- data.frame(
  comparison = "MODEL_DIAGNOSTICS",
  OR_median = NA,
  OR_lower = NA,
  OR_upper = NA,
  CrI_excludes_1 = NA,
  SUCRA = NA,
  stringsAsFactors = FALSE
)

summary_df <- do.call(rbind, summary_rows)

# Add metadata columns
summary_df$DIC <- dic
summary_df$max_PSRF <- max_psrf
summary_df$converged <- converged
summary_df$escalated <- escalated
summary_df$n_treatments <- length(TREATMENT_NAMES)
summary_df$n_studies <- length(unique(nma_data_clean$s.id))
summary_df$n_arms <- nrow(nma_data_clean)

summary_file <- file.path(output_dir, "genotyping_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("Saved:", summary_file, "\n")

# =============================================================================
# PHASE 7: Print summary to stdout
# =============================================================================

cat("\n###############################################################\n")
cat("### SUMMARY                                                 ###\n")
cat("###############################################################\n\n")

cat("=== Genotyping Misclassification Sensitivity Analysis ===\n\n")

cat(sprintf("Data: %d studies, %d arms, %d treatments\n",
            length(unique(nma_data_clean$s.id)), nrow(nma_data_clean),
            length(TREATMENT_NAMES)))
cat(sprintf("Model: het_cor, probit link, invwishart prior\n"))
cat(sprintf("MCMC: %s\n",
            ifelse(escalated, "escalated (200K iter, 100K burn-in, thin 20)",
                   "standard (150K iter, 75K burn-in, thin 10)")))
cat(sprintf("Convergence: %s (max PSRF = %s)\n",
            ifelse(converged, "CONVERGED", "NOT CONVERGED"),
            ifelse(is.na(max_psrf), "N/A", sprintf("%.4f", max_psrf))))
cat(sprintf("DIC: %s\n\n", ifelse(is.na(dic), "N/A", sprintf("%.2f", dic))))

cat("--- Key Comparisons ---\n\n")

if (!is.null(hapb3_vs_wtclean)) {
  cat(sprintf("  HapB3 vs WT_Clean:    OR = %.4f (%.4f, %.4f)  %s\n",
              hapb3_vs_wtclean$median, hapb3_vs_wtclean$lower, hapb3_vs_wtclean$upper,
              ifelse(hapb3_vs_wtclean$lower > 1.0, "[Significant]", "[Not significant]")))
}

if (!is.null(hapb3_vs_wtbiased)) {
  cat(sprintf("  HapB3 vs WT_Biased:   OR = %.4f (%.4f, %.4f)  %s\n",
              hapb3_vs_wtbiased$median, hapb3_vs_wtbiased$lower, hapb3_vs_wtbiased$upper,
              ifelse(hapb3_vs_wtbiased$lower > 1.0, "[Significant]", "[Not significant]")))
}

if (!is.null(wtclean_vs_wtbiased)) {
  cat(sprintf("  WT_Clean vs WT_Biased: OR = %.4f (%.4f, %.4f)  %s\n",
              wtclean_vs_wtbiased$median, wtclean_vs_wtbiased$lower, wtclean_vs_wtbiased$upper,
              ifelse(wtclean_vs_wtbiased$lower > 1.0 || wtclean_vs_wtbiased$upper < 1.0,
                     "[WT groups DIFFER]", "[WT groups similar]")))
}

cat("\n--- SUCRA Rankings (higher = lower toxicity) ---\n\n")
if (!is.null(sucra)) {
  sucra_sorted <- sort(sucra, decreasing = TRUE)
  for (i in seq_along(sucra_sorted)) {
    cat(sprintf("  %d. %-12s: %.1f%%\n", i, names(sucra_sorted)[i], sucra_sorted[i] * 100))
  }
}

cat("\n--- Interpretation ---\n\n")

if (!is.null(hapb3_vs_wtclean) && !is.null(hapb3_vs_wtbiased)) {
  or_diff <- abs(log(hapb3_vs_wtclean$median) - log(hapb3_vs_wtbiased$median))
  cat(sprintf("  Log-OR difference (HapB3 vs WT_Clean vs HapB3 vs WT_Biased): %.4f\n", or_diff))
  if (or_diff < 0.20) {
    cat("  --> HapB3 OR is SIMILAR regardless of WT genotyping quality.\n")
  } else {
    cat("  --> HapB3 OR DIFFERS between WT genotyping groups. Genotyping quality may matter.\n")
  }
}

if (!is.null(wtclean_vs_wtbiased)) {
  if (wtclean_vs_wtbiased$lower > 1.0 || wtclean_vs_wtbiased$upper < 1.0) {
    cat("  WT_Clean and WT_Biased show a SIGNIFICANT difference in toxicity rate.\n")
    cat("  This suggests incomplete genotyping panels may introduce misclassification bias.\n")
  } else {
    cat("  WT_Clean and WT_Biased show NO significant difference.\n")
    cat("  Genotyping panel quality does not materially affect the WT baseline estimate.\n")
  }
}

# =============================================================================
# Complete
# =============================================================================

cat("\n=============================================================\n")
cat("Genotyping Sensitivity Analysis Complete\n")
cat(sprintf("Finished: %s\n", Sys.time()))
cat("=============================================================\n")
cat("\nOutputs:\n")
cat(sprintf("  Result:       %s\n", file.path(output_dir, "genotyping_result.rds")))
cat(sprintf("  Convergence:  %s\n", file.path(output_dir, "genotyping_convergence.rds")))
cat(sprintf("  Summary:      %s\n", summary_file))
cat(sprintf("  Log:          %s\n", log_file))
