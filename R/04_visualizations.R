# =============================================================================
# 04_visualizations.R
# Generate standard NMA visualizations using pcnetmeta
# =============================================================================

library(pcnetmeta)

# -----------------------------------------------------------------------------
# Load data and results
# -----------------------------------------------------------------------------

nma_data <- readRDS("output/nma_data.rds")
result <- readRDS("output/nma_result.rds")
treatment_names <- readRDS("output/treatment_names.rds")

cat("=== Generating NMA Visualizations ===\n\n")

# -----------------------------------------------------------------------------
# 1. Network Plot - Visualize network structure
# -----------------------------------------------------------------------------

cat("Creating network plot...\n")

png("output/network_plot.png", width = 800, height = 800, res = 120)
nma.networkplot(
  s.id = s.id,
  t.id = t.id,
  data = nma_data,
  trtname = treatment_names,
  title = "DPYD Network Meta-Analysis Structure"
)
dev.off()

cat("  Saved: output/network_plot.png\n")

# -----------------------------------------------------------------------------
# 2. Contrast Plot - All treatments vs WT_Clean reference
# -----------------------------------------------------------------------------

cat("Creating contrast plot...\n")

png("output/contrast_plot.png", width = 1000, height = 600, res = 120)
contrast.plot(
  result,
  reference = treatment_names[1],  # WT_Clean is the reference
  digits = 2
)
dev.off()

cat("  Saved: output/contrast_plot.png\n")

# -----------------------------------------------------------------------------
# 3. Absolute Effects Plot - Toxicity rates per treatment
# -----------------------------------------------------------------------------

cat("Creating absolute effects plot...\n")

png("output/absolute_effects_plot.png", width = 1000, height = 600, res = 120)
absolute.plot(result, digits = 2)
dev.off()

cat("  Saved: output/absolute_effects_plot.png\n")

# -----------------------------------------------------------------------------
# 4. Ranking Plot - Treatment ranking probabilities
# -----------------------------------------------------------------------------

cat("Creating ranking plot...\n")

png("output/ranking_plot.png", width = 1000, height = 600, res = 120)
rank.prob(result)
dev.off()

cat("  Saved: output/ranking_plot.png\n")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

cat("\n=== Visualizations Complete ===\n")
cat("Generated plots:\n")
cat("  - output/network_plot.png          (network structure)\n")
cat("  - output/contrast_plot.png         (ORs vs WT_Clean reference)\n")
cat("  - output/absolute_effects_plot.png (toxicity rates per treatment)\n")
cat("  - output/ranking_plot.png          (treatment ranking probabilities)\n")
