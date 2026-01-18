# =============================================================================
# R/extract_heterogeneity.R
# Extract heterogeneity parameters (sigma, R) from arm-based NMA
#
# The pcnetmeta het_cor model does not monitor Sigma parameters by default.
# This script re-runs a custom JAGS model to extract treatment-specific
# standard deviations and between-study correlations for Table S7.
# =============================================================================

cat("\n")
cat("###############################################################\n")
cat("### HETEROGENEITY PARAMETER EXTRACTION                      ###\n")
cat("###############################################################\n\n")

suppressPackageStartupMessages({
  library(rjags)
  library(coda)
})

# -----------------------------------------------------------------------------
# Load data from previous NMA run
# -----------------------------------------------------------------------------

cat("Loading NMA data and treatment names...\n")

# Use paths relative to working directory (set by run_analysis.R)
nma_data <- readRDS("output/wt_unified/wt_unified_nma_data.rds")
treatment_names <- readRDS("output/wt_unified/wt_unified_treatment_names.rds")

cat("  Data loaded:", nrow(nma_data), "arms\n")
cat("  Treatments:", paste(treatment_names, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# Custom JAGS model with sigma monitoring
# -----------------------------------------------------------------------------

cat("Building custom JAGS model with sigma/R monitoring...\n")

# Model code from pcnetmeta het_cor invwishart, modified to monitor sigma
model_string <- "
model{
  for(i in 1:len){
    p[i] <- phi(mu[t[i]] + vi[s[i], t[i]])
    r[i] ~ dbin(p[i], totaln[i])
    rhat[i] <- p[i]*totaln[i]
    dev[i] <- 2*(r[i]*(log(r[i]) - log(rhat[i])) +
      (totaln[i] - r[i])*(log(totaln[i] - r[i]) - log(totaln[i] - rhat[i])))
  }
  totresdev <- sum(dev[])
  for(j in 1:nstudy){
    vi[j, 1:ntrt] ~ dmnorm(zeros[1:ntrt], T[1:ntrt, 1:ntrt])
  }
  for(j in 1:ntrt){
    AR[j] <- phi(mu[j]/sqrt(1 + invT[j,j]))
    mu[j] ~ dnorm(0,0.001)
    sigma[j] <- sqrt(invT[j,j])
  }
  invT[1:ntrt, 1:ntrt] <- inverse(T[,])
  T[1:ntrt, 1:ntrt] ~ dwish(I[1:ntrt, 1:ntrt], ntrt + 1)

  # Correlation matrix (derived from covariance)
  for(j in 1:ntrt){
    for(k in 1:ntrt){
      R[j,k] <- invT[j,k] / (sigma[j] * sigma[k])
    }
  }

  for(j in 1:ntrt){
    for(k in 1:ntrt){
      LOR[j,k] <- log(OR[j,k])
      OR[j,k] <- AR[j]/(1 - AR[j])/AR[k]*(1 - AR[k])
    }
  }
}
"

# -----------------------------------------------------------------------------
# Prepare data for JAGS
# -----------------------------------------------------------------------------

nstudy <- length(unique(nma_data$s.id))
ntrt <- length(treatment_names)

jags_data <- list(
  s = nma_data$s.id,
  t = nma_data$t.id,
  r = nma_data$r,
  totaln = nma_data$n,
  len = nrow(nma_data),
  nstudy = nstudy,
  ntrt = ntrt,
  zeros = rep(0, ntrt),
  I = diag(ntrt)
)

cat("  JAGS data: ", jags_data$len, " arms, ",
    jags_data$nstudy, " studies, ", jags_data$ntrt, " treatments\n\n")

# Initial values function
init_func <- function() {
  list(
    mu = qnorm(rep(0.5, ntrt)),
    T = diag(ntrt)
  )
}

# Parameters to monitor
params_to_monitor <- c("sigma", "R", "AR", "LOR", "totresdev")

# -----------------------------------------------------------------------------
# Run JAGS
# -----------------------------------------------------------------------------

cat("Running JAGS model...\n")
cat("  Adaptation: 5000, Burn-in: 25000, Sampling: 50000, Thin: 5, Chains: 3\n")

set.seed(12345)

model_file <- tempfile(fileext = ".txt")
writeLines(model_string, model_file)

jags_model <- jags.model(
  file = model_file,
  data = jags_data,
  inits = init_func,
  n.chains = 3,
  n.adapt = 5000,
  quiet = TRUE
)

cat("  Burn-in...\n")
update(jags_model, n.iter = 25000, progress.bar = "none")

cat("  Sampling...\n")
samples <- coda.samples(
  model = jags_model,
  variable.names = params_to_monitor,
  n.iter = 50000,
  thin = 5,
  progress.bar = "none"
)

cat("  Sampling complete.\n\n")

# -----------------------------------------------------------------------------
# Extract and summarize sigma parameters
# -----------------------------------------------------------------------------

samples_matrix <- do.call(rbind, samples)
param_names <- colnames(samples_matrix)

# Sigma (treatment-specific SDs)
sigma_params <- param_names[grepl("^sigma\\[", param_names)]

sigma_summary <- data.frame(
  Parameter = character(),
  Treatment = character(),
  Posterior_Mean = numeric(),
  Posterior_SD = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
)

cat("Treatment-specific standard deviations (probit scale):\n")
for (i in seq_along(sigma_params)) {
  p <- sigma_params[i]
  vals <- samples_matrix[, p]

  row <- data.frame(
    Parameter = p,
    Treatment = treatment_names[i],
    Posterior_Mean = mean(vals),
    Posterior_SD = sd(vals),
    CI_Lower = quantile(vals, 0.025, names = FALSE),
    CI_Upper = quantile(vals, 0.975, names = FALSE)
  )
  sigma_summary <- rbind(sigma_summary, row)

  cat(sprintf("  %s (%s): %.4f [%.4f, %.4f]\n",
              p, treatment_names[i], mean(vals),
              quantile(vals, 0.025), quantile(vals, 0.975)))
}

# -----------------------------------------------------------------------------
# Extract correlation matrix
# -----------------------------------------------------------------------------

cat("\nBetween-study correlations:\n")

r_params <- param_names[grepl("^R\\[", param_names)]
r_summary <- data.frame(
  Parameter = character(),
  Comparison = character(),
  Posterior_Mean = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
)

R_mean <- matrix(NA, ntrt, ntrt)

for (j in 1:ntrt) {
  for (k in 1:ntrt) {
    p_name <- sprintf("R[%d,%d]", j, k)
    if (p_name %in% param_names) {
      vals <- samples_matrix[, p_name]
      R_mean[j,k] <- mean(vals)

      if (j < k) {
        row <- data.frame(
          Parameter = p_name,
          Comparison = sprintf("%s vs %s", treatment_names[j], treatment_names[k]),
          Posterior_Mean = mean(vals),
          CI_Lower = quantile(vals, 0.025, names = FALSE),
          CI_Upper = quantile(vals, 0.975, names = FALSE)
        )
        r_summary <- rbind(r_summary, row)

        cat(sprintf("  %s vs %s: %.3f [%.3f, %.3f]\n",
                    treatment_names[j], treatment_names[k],
                    mean(vals), quantile(vals, 0.025), quantile(vals, 0.975)))
      }
    }
  }
}

rownames(R_mean) <- treatment_names
colnames(R_mean) <- treatment_names

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------

cat("\nSaving heterogeneity outputs...\n")

output_dir <- "output/wt_unified"

# Combined summary CSV
full_summary <- rbind(
  data.frame(
    Parameter = sigma_summary$Parameter,
    Type = "Variance (SD)",
    Posterior_Mean = sigma_summary$Posterior_Mean,
    CI_Lower = sigma_summary$CI_Lower,
    CI_Upper = sigma_summary$CI_Upper,
    Interpretation = paste("Between-study SD for", sigma_summary$Treatment)
  ),
  data.frame(
    Parameter = r_summary$Parameter,
    Type = "Correlation",
    Posterior_Mean = r_summary$Posterior_Mean,
    CI_Lower = r_summary$CI_Lower,
    CI_Upper = r_summary$CI_Upper,
    Interpretation = paste("Correlation:", r_summary$Comparison)
  )
)

write.csv(full_summary, file.path(output_dir, "wt_unified_heterogeneity_params.csv"), row.names = FALSE)
cat("  - wt_unified_heterogeneity_params.csv\n")

# MCMC samples
saveRDS(samples, file.path(output_dir, "wt_unified_heterogeneity_mcmc.rds"))
cat("  - wt_unified_heterogeneity_mcmc.rds\n")

# Formatted text output
sink(file.path(output_dir, "wt_unified_heterogeneity_summary.txt"))
cat("=== DPYD Arm-Based NMA: Heterogeneity Parameters ===\n\n")
cat("Model: het_cor (heterogeneous variance with correlation)\n")
cat("Prior: Wishart on precision matrix T (df = ntrt + 1 = 6)\n")
cat("       Equivalent to Inverse-Wishart on covariance matrix Sigma\n\n")

cat("=== Treatment-Specific Standard Deviations (sigma) ===\n")
cat("(On probit scale; represent between-study heterogeneity for each treatment)\n\n")
for (i in 1:nrow(sigma_summary)) {
  cat(sprintf("%s (%s): %.4f [95%% CrI: %.4f, %.4f]\n",
              sigma_summary$Parameter[i],
              sigma_summary$Treatment[i],
              sigma_summary$Posterior_Mean[i],
              sigma_summary$CI_Lower[i],
              sigma_summary$CI_Upper[i]))
}

cat("\n=== Correlation Matrix ===\n")
cat("(Between-study random effect correlations across treatments)\n\n")
print(round(R_mean, 3))

cat("\n=== Off-Diagonal Correlations with 95% CrI ===\n")
for (i in 1:nrow(r_summary)) {
  cat(sprintf("%s: %.3f [%.3f, %.3f]\n",
              r_summary$Comparison[i],
              r_summary$Posterior_Mean[i],
              r_summary$CI_Lower[i],
              r_summary$CI_Upper[i]))
}

cat("\n=== Model Fit ===\n")
totresdev_vals <- samples_matrix[, "totresdev"]
cat(sprintf("Total residual deviance: %.2f [%.2f, %.2f]\n",
            mean(totresdev_vals),
            quantile(totresdev_vals, 0.025),
            quantile(totresdev_vals, 0.975)))
cat(sprintf("Number of data points: %d\n", nrow(nma_data)))
cat(sprintf("Dres/n ratio: %.2f (should be ~1 for good fit)\n",
            mean(totresdev_vals) / nrow(nma_data)))
sink()
cat("  - wt_unified_heterogeneity_summary.txt\n")

cat("\nHeterogeneity extraction complete.\n")
