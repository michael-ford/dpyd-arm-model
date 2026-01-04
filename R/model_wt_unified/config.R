# =============================================================================
# R/model_wt_unified/config.R
# Configuration for WT Unified model (WT_Clean and WT_Biased merged into WT)
# =============================================================================

# Model identifier
MODEL_NAME <- "wt_unified"
MODEL_DESCRIPTION <- "WT Unified: WT_Clean and WT_Biased merged into single WT reference node"

# Treatment mapping (order determines reference: treatment 1 = reference)
# Both WT_Clean and WT_Biased map to treatment 1 (WT)
# Note: Data uses "HapB3_1129or1236" but we display as "HapB3"
TREATMENT_MAP <- c(
  "WT_Clean"         = 1,
  "WT_Biased"        = 1,  # Maps to same ID as WT_Clean (will be aggregated)
  "HapB3_1129or1236" = 2,
  "2846hetho"        = 3,
  "2Ahetho"          = 4,
  "13hetho"          = 5
)

# Treatment names vector (indexed by treatment ID)
# Only 5 treatments in unified model
TREATMENT_NAMES <- c(
  "WT",
  "HapB3",
  "2846hetho",
  "2Ahetho",
  "13hetho"
)

# Output directory for this model
OUTPUT_DIR <- "output/wt_unified"

# Data file path
DATA_FILE <- "data/Binary WT HapB3 Data for NMA (12-18-2025).xlsx"

# Key comparison of interest
KEY_COMPARISON_TREATMENT <- "13hetho"
