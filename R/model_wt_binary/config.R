# =============================================================================
# R/model_wt_binary/config.R
# Configuration for WT Binary model (WT_Clean and WT_Biased as separate nodes)
# =============================================================================

# Model identifier
MODEL_NAME <- "wt_binary"
MODEL_DESCRIPTION <- "WT Binary: WT_Clean and WT_Biased as separate treatment nodes"

# Treatment mapping (order determines reference: treatment 1 = reference)
# WT_Clean is reference for OR calculations
# Note: Data uses "HapB3_1129or1236" but we display as "HapB3"
TREATMENT_MAP <- c(
  "WT_Clean"         = 1,
  "WT_Biased"        = 2,
  "HapB3_1129or1236" = 3,
  "2846hetho"        = 4,
  "2Ahetho"          = 5,
  "13hetho"          = 6
)

# Treatment names vector (indexed by treatment ID)
# Short display names for output/visualization
TREATMENT_NAMES <- c(
  "WT_Clean",
  "WT_Biased",
  "HapB3",
  "2846hetho",
  "2Ahetho",
  "13hetho"
)

# Output directory for this model
OUTPUT_DIR <- "output/wt_binary"

# Data file path
DATA_FILE <- "data/Binary WT HapB3 Data for NMA (12-18-2025).xlsx"

# Key comparison of interest
KEY_COMPARISON_TREATMENT <- "13hetho"
