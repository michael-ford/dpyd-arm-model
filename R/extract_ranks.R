# =============================================================================
# R/extract_ranks.R
# Extract treatment rank probabilities from pcnetmeta result object
# =============================================================================

cat("\n")
cat("###############################################################\n")
cat("### RANK PROBABILITY EXTRACTION                             ###\n")
cat("###############################################################\n\n")

# -----------------------------------------------------------------------------
# Load result and treatment names
# -----------------------------------------------------------------------------

cat("Loading NMA result object...\n")

result <- readRDS("output/wt_unified/wt_unified_result.rds")
treatment_names <- readRDS("output/wt_unified/wt_unified_treatment_names.rds")

cat("  Treatments:", paste(treatment_names, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# Extract rank probabilities
# -----------------------------------------------------------------------------

if (!"TrtRankProb" %in% names(result)) {
  stop("TrtRankProb not found in result object")
}

rank_probs <- result$TrtRankProb
ntrt <- length(treatment_names)

# Convert to numeric matrix
rank_numeric <- matrix(as.numeric(rank_probs), nrow = ntrt, ncol = ntrt)
rownames(rank_numeric) <- treatment_names
colnames(rank_numeric) <- paste0("Rank", 1:ntrt)

cat("Rank probability matrix:\n")
cat("(Rows = Treatments, Cols = Rank; Rank 1 = lowest toxicity/best)\n\n")
print(round(rank_numeric, 4))

# -----------------------------------------------------------------------------
# Summary statistics
# -----------------------------------------------------------------------------

cat("\n=== Most Likely Ranks ===\n")
for (i in 1:ntrt) {
  best_rank <- which.max(rank_numeric[i,])
  prob_best <- max(rank_numeric[i,])
  cat(sprintf("  %s: Rank %d (P = %.1f%%)\n",
              treatment_names[i], best_rank, prob_best * 100))
}

# -----------------------------------------------------------------------------
# Calculate SUCRA
# -----------------------------------------------------------------------------

cat("\n=== SUCRA (Surface Under Cumulative Ranking) ===\n")

sucra <- numeric(ntrt)
for (i in 1:ntrt) {
  cumulative_probs <- cumsum(rank_numeric[i, 1:(ntrt-1)])
  sucra[i] <- sum(cumulative_probs) / (ntrt - 1)
}
names(sucra) <- treatment_names

sucra_df <- data.frame(
  Treatment = treatment_names,
  SUCRA = round(sucra * 100, 2)
)
sucra_df <- sucra_df[order(-sucra_df$SUCRA), ]

for (i in 1:nrow(sucra_df)) {
  cat(sprintf("  %s: %.2f%%\n", sucra_df$Treatment[i], sucra_df$SUCRA[i]))
}

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------

cat("\nSaving rank probability outputs...\n")

output_dir <- "output/wt_unified"

# Rank probabilities CSV
write.csv(rank_numeric, file.path(output_dir, "wt_unified_rank_probabilities.csv"))
cat("  - wt_unified_rank_probabilities.csv\n")

# SUCRA CSV
sucra_out <- data.frame(
  Treatment = treatment_names,
  SUCRA = round(sucra * 100, 2)
)
write.csv(sucra_out, file.path(output_dir, "wt_unified_sucra.csv"), row.names = FALSE)
cat("  - wt_unified_sucra.csv\n")

cat("\nRank extraction complete.\n")
