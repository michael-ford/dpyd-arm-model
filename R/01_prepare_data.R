# =============================================================================
# 01_prepare_data.R
# Prepare Excel data for pcnetmeta arm-based NMA
# =============================================================================

library(readxl)
library(dplyr)

# -----------------------------------------------------------------------------
# Pre-flight validation
# -----------------------------------------------------------------------------

data_file <- "data/Binary WT HapB3 Data for NMA (12-18-2025).xlsx"

# Check that data file exists
if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file,
       "\nPlease ensure the Excel file is in the data/ directory.")
}

# -----------------------------------------------------------------------------
# Load raw data
# -----------------------------------------------------------------------------

raw_data <- read_excel(data_file)

# Verify required columns exist
required_cols <- c("Study", "T", "R", "N")
missing_cols <- setdiff(required_cols, names(raw_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
       "\nExpected columns: ", paste(required_cols, collapse = ", "),
       "\nFound columns: ", paste(names(raw_data), collapse = ", "))
}

# Check for unexpected NA values in required columns
na_check <- sapply(required_cols, function(col) sum(is.na(raw_data[[col]])))
cols_with_na <- names(na_check[na_check > 0])
if (length(cols_with_na) > 0) {
  for (col in cols_with_na) {
    na_rows <- which(is.na(raw_data[[col]]))
    stop("Unexpected NA values in column '", col, "' at rows: ",
         paste(na_rows, collapse = ", "))
  }
}

cat("=== Raw Data Summary ===\n")
cat("Studies:", n_distinct(raw_data$Study), "\n")
cat("Total arms:", nrow(raw_data), "\n")
cat("Treatments:", paste(unique(raw_data$T), collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# Define treatment mapping
# WT_Clean = 1 (reference treatment)
# -----------------------------------------------------------------------------

treatment_map <- c(
  "WT_Clean"         = 1,
  "WT_Biased"        = 2,
  "HapB3_1129or1236" = 3,
  "2846hetho"        = 4,
  "2Ahetho"          = 5,
  "13hetho"          = 6
)

treatment_names <- c(
  "WT_Clean",
  "WT_Biased",
  "HapB3",
  "2846hetho",
  "2Ahetho",
  "13hetho"
)

# Validate all treatments in data match treatment_map (before transformation)
data_treatments <- unique(raw_data$T)
unmapped_treatments <- setdiff(data_treatments, names(treatment_map))
if (length(unmapped_treatments) > 0) {
  stop("Treatments in data not found in treatment_map: ",
       paste(unmapped_treatments, collapse = ", "),
       "\nValid treatments: ", paste(names(treatment_map), collapse = ", "),
       "\nPlease update treatment_map or check data for typos.")
}

# -----------------------------------------------------------------------------
# Transform to pcnetmeta format
# -----------------------------------------------------------------------------

nma_data <- raw_data %>%
  mutate(
    s.id = as.numeric(factor(Study, levels = unique(Study))),  # Preserve input order
    t.id = treatment_map[T],           # Numeric treatment ID
    r = R,                             # Events
    n = N                              # Total
  ) %>%
  select(s.id, t.id, r, n, Study, T) %>%
  arrange(s.id, t.id)

# -----------------------------------------------------------------------------
# Validate data
# -----------------------------------------------------------------------------

cat("=== Data Validation ===\n")

# Check for missing values
if (any(is.na(nma_data$s.id)) || any(is.na(nma_data$t.id))) {
  stop("Missing study or treatment IDs after transformation")
}

# Check r <= n (error - this should never happen)
invalid_arms <- nma_data %>% filter(r > n)
if (nrow(invalid_arms) > 0) {
  cat("ERROR: Events exceed total in these arms:\n")
  print(invalid_arms %>% select(Study, T, r, n))
  stop("Data integrity error: R > N in ", nrow(invalid_arms), " arm(s)")
}

cat("All validations passed\n\n")

# -----------------------------------------------------------------------------
# Data quality warnings
# -----------------------------------------------------------------------------

cat("=== Data Quality Warnings ===\n")
warning_count <- 0

# Warning: Single-arm studies (valid but noteworthy)
study_arm_counts <- nma_data %>%
  group_by(s.id, Study) %>%
  summarise(n_arms = n(), .groups = "drop")

single_arm_studies <- study_arm_counts %>% filter(n_arms == 1)
if (nrow(single_arm_studies) > 0) {
  warning_count <- warning_count + 1
  cat("WARNING: ", nrow(single_arm_studies), " single-arm study(ies) found:\n", sep = "")
  for (i in seq_len(nrow(single_arm_studies))) {
    cat("  - ", single_arm_studies$Study[i], "\n", sep = "")
  }
  cat("  (Single-arm studies are valid but provide less information for NMA)\n\n")
}

# Warning: Very small total N for a treatment (threshold: < 20)
small_n_threshold <- 20
treatment_totals <- nma_data %>%
  group_by(t.id, T) %>%
  summarise(total_n = sum(n), n_arms = n(), .groups = "drop")

small_treatments <- treatment_totals %>% filter(total_n < small_n_threshold)
if (nrow(small_treatments) > 0) {
  warning_count <- warning_count + 1
  cat("WARNING: Treatments with very small total N (< ", small_n_threshold, "):\n", sep = "")
  for (i in seq_len(nrow(small_treatments))) {
    cat("  - ", small_treatments$T[i], ": N = ", small_treatments$total_n[i],
        " across ", small_treatments$n_arms[i], " arm(s)\n", sep = "")
  }
  cat("  (Small samples may lead to unstable estimates)\n\n")
}

if (warning_count == 0) {
  cat("No data quality warnings.\n\n")
} else {
  cat("Total warnings: ", warning_count, "\n\n", sep = "")
}

# -----------------------------------------------------------------------------
# Report sparse data (key for arm-based model justification)
# -----------------------------------------------------------------------------

cat("=== Sparse Data Report ===\n")

zero_events <- nma_data %>% filter(r == 0)
cat("Zero-event arms:", nrow(zero_events), "\n")
if (nrow(zero_events) > 0) {
  print(zero_events %>% select(Study, T, n, r))
}

small_arms <- nma_data %>% filter(n <= 5)
cat("\nSmall sample arms (N <= 5):", nrow(small_arms), "\n")

# Check the problematic Froehlich study
froehlich <- nma_data %>% filter(grepl("Froehlich", Study))
cat("\n=== Froehlich 2015 (problematic study) ===\n")
print(froehlich %>% select(Study, T, n, r))

# -----------------------------------------------------------------------------
# Summary by treatment
# -----------------------------------------------------------------------------

cat("\n=== Treatment Summary ===\n")
treatment_summary <- nma_data %>%
  group_by(t.id, T) %>%
  summarise(
    n_arms = n(),
    total_n = sum(n),
    total_r = sum(r),
    event_rate = round(sum(r) / sum(n), 3),
    .groups = "drop"
  ) %>%
  arrange(t.id)

print(treatment_summary)

# -----------------------------------------------------------------------------
# Save prepared data
# -----------------------------------------------------------------------------

# Keep only columns needed for pcnetmeta
nma_data_clean <- nma_data %>%
  select(s.id, t.id, r, n)

saveRDS(nma_data_clean, "output/nma_data.rds")
saveRDS(treatment_names, "output/treatment_names.rds")

# Also save study mapping for reference
study_mapping <- raw_data %>%
  select(Study) %>%
  distinct() %>%
  mutate(s.id = as.numeric(factor(Study, levels = unique(Study))))

saveRDS(study_mapping, "output/study_mapping.rds")

cat("\n=== Data Preparation Complete ===\n")
cat("Saved: output/nma_data.rds\n")
cat("Saved: output/treatment_names.rds\n")
cat("Saved: output/study_mapping.rds\n")
