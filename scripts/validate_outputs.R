# validate_outputs.R
# Script to validate generated outputs from the DPYD arm-based NMA analysis

# Define expected files
expected_figures <- c(
  "figure2_network.png",
  "figure5_forest.png",
  "figure6_rankogram.png"
)

expected_tables <- c(
  "tableS5_absolute_risks.csv",
  "tableS6_diagnostics.csv",
  "tableS7_heterogeneity.csv",
  "table2_pairwise_or.csv"
)

# Set directory paths
fig_dir <- "output/wt_unified/figures"
table_dir <- "output/wt_unified/tables"

# Print header
cat("=== Validating Outputs ===\n\n")

# Initialize counters
fig_pass <- 0
fig_fail <- 0
table_pass <- 0
table_fail <- 0

# Check figures
cat("Figures:\n")
for (fig in expected_figures) {
  fig_path <- file.path(fig_dir, fig)
  if (file.exists(fig_path)) {
    size_kb <- round(file.info(fig_path)$size / 1024, 1)
    cat(sprintf("  ✓ %s (%s KB)\n", fig, size_kb))
    fig_pass <- fig_pass + 1
  } else {
    cat(sprintf("  ✗ %s - MISSING\n", fig))
    fig_fail <- fig_fail + 1
  }
}

cat("\n")

# Check tables
cat("Tables:\n")
for (tbl in expected_tables) {
  tbl_path <- file.path(table_dir, tbl)
  if (file.exists(tbl_path)) {
    row_count <- nrow(read.csv(tbl_path))
    cat(sprintf("  ✓ %s (%d rows)\n", tbl, row_count))
    table_pass <- table_pass + 1
  } else {
    cat(sprintf("  ✗ %s - MISSING\n", tbl))
    table_fail <- table_fail + 1
  }
}

cat("\n")

# Print summary
cat("=== Summary ===\n")
cat(sprintf("Figures: %d passed, %d failed\n", fig_pass, fig_fail))
cat(sprintf("Tables:  %d passed, %d failed\n", table_pass, table_fail))

total_pass <- fig_pass + table_pass
total_fail <- fig_fail + table_fail
cat(sprintf("Total:   %d passed, %d failed\n", total_pass, total_fail))

if (total_fail == 0) {
  cat("\nAll outputs validated successfully!\n")
} else {
  cat(sprintf("\nWARNING: %d output(s) missing!\n", total_fail))
}
