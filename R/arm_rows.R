# Shared helpers used by the final fixed-study analysis.

make_arm_rows <- function(study, treatment, n, r, duration = 1) {
  if (length(study) == 1L) study <- rep(study, length(treatment))
  if (length(study) != length(treatment)) stop("study and treatment lengths differ.", call. = FALSE)
  data.frame(
    Study = study,
    T = treatment,
    N = n,
    R = r,
    covar.duration = duration,
    stringsAsFactors = FALSE
  )
}

normalize_arm_rows <- function(rows) {
  rows <- as.data.frame(rows[, c("Study", "T", "N", "R"), drop = FALSE], stringsAsFactors = FALSE)
  rows$Study <- as.character(rows$Study)
  rows$T <- as.character(rows$T)
  rows$N <- as.numeric(rows$N)
  rows$R <- as.numeric(rows$R)
  rows <- rows[order(rows$Study, rows$T), , drop = FALSE]
  rownames(rows) <- NULL
  rows
}

bind_replacement <- function(data, replacement) {
  for (column in setdiff(names(data), names(replacement))) replacement[[column]] <- NA
  replacement <- replacement[, names(data), drop = FALSE]
  rbind(data, replacement)
}

replace_study <- function(data, study, replacement) {
  data <- data[data$Study != study, , drop = FALSE]
  bind_replacement(data, replacement)
}

verified_replacements <- function(etienne_wt_n = 232) {
  if (!etienne_wt_n %in% c(231, 232)) stop("etienne_wt_n must be 231 or 232.", call. = FALSE)
  list(
    Etienne_Grimaldi_2017 = make_arm_rows(
      "Etienne_Grimaldi_2017",
      c("13hetho", "HapB3_1129or1236", "2846hetho", "2Ahetho", "WT_Biased"),
      c(1, 4, 3, 3, etienne_wt_n), c(1, 0, 2, 2, 25), c(1, 1, 1, 1, 2)
    ),
    Jennings_2013 = make_arm_rows(
      "Jennings_2013", c("HapB3_1129or1236", "2846hetho", "2Ahetho", "WT_Clean"),
      c(10, 2, 3, 238), c(3, 2, 2, 37)
    ),
    Bozina_2021 = make_arm_rows(
      "Bozina_2021", c("HapB3_1129or1236", "13hetho", "2846hetho", "2Ahetho", "WT_Biased"),
      c(22, 1, 4, 17, 456), c(10, 1, 4, 17, 249)
    ),
    Toffoli_2015 = make_arm_rows(
      "Toffoli_2015", c("13hetho", "2846hetho", "2Ahetho", "WT_Biased"),
      c(1, 5, 12, 585), c(0, 3, 8, 84)
    ),
    Loganayagam_2013 = make_arm_rows(
      "Loganayagam_2013", c("13hetho", "HapB3_1129or1236", "2846hetho", "2Ahetho", "WT_Clean"),
      c(1, 15, 4, 3, 391), c(1, 6, 4, 3, 74)
    ),
    Salgueiro_2004 = make_arm_rows(
      "Salgueiro_2004", c("2Ahetho", "WT_Biased"), c(1, 71), c(1, 6), c(1, 2)
    ),
    Morel_2006 = make_arm_rows(
      "Morel_2006", c("13hetho", "2846hetho", "2Ahetho", "WT_Biased"),
      c(1, 10, 9, 466), c(1, 6, 5, 31), c(1, 1, 1, 2)
    ),
    Largillier_2006 = make_arm_rows(
      "Largillier_2006", c("2Ahetho", "WT_Biased"), c(1, 82), c(1, 14), c(1, 2)
    )
  )
}

validate_verified_scenario <- function(data) {
  required <- c("Study", "T", "N", "R")
  if (!all(required %in% names(data))) stop("Scenario is missing required columns.", call. = FALSE)
  if (any(is.na(data[, required]))) stop("Scenario has missing required values.", call. = FALSE)
  if (anyDuplicated(data[, c("Study", "T")])) stop("Scenario has duplicate study-treatment rows.", call. = FALSE)
  if (any(data$N < 0) || any(data$R < 0) || any(data$R > data$N)) stop("Scenario has invalid N/R values.", call. = FALSE)
  if (any(data$N %% 1 != 0) || any(data$R %% 1 != 0)) stop("Scenario N/R values must be integers.", call. = FALSE)
  invisible(TRUE)
}
