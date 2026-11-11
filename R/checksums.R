# Shared helpers used by the final fixed-study analysis.

sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) stop("File unavailable for SHA-256: ", path, call. = FALSE)
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path, call. = FALSE)
  sub(" .*", "", output[1])
}

sha256_directory <- function(path) {
  if (!dir.exists(path)) stop("Directory unavailable for SHA-256: ", path, call. = FALSE)
  command <- paste(
    "find", shQuote(path), "-type f -print | LC_ALL=C sort |",
    "while IFS= read -r f; do sha256sum \"$f\"; done | sha256sum"
  )
  output <- system2("sh", c("-c", shQuote(command)), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("Directory SHA-256 failed for ", path, call. = FALSE)
  sub(" .*", "", output[1])
}

sha256_path <- function(path, kind) {
  if (identical(kind, "directory")) sha256_directory(path) else sha256_file(path)
}

checksum_protected_artifacts <- function(artifacts) {
  hashes <- vapply(seq_len(nrow(artifacts)), function(index) {
    sha256_path(artifacts$path[index], artifacts$kind[index])
  }, character(1))
  data.frame(
    artifact = artifacts$artifact,
    path = artifacts$path,
    kind = artifacts$kind,
    expected_sha256 = artifacts$expected_sha256,
    observed_sha256 = hashes,
    matches_expected = hashes == artifacts$expected_sha256,
    stringsAsFactors = FALSE
  )
}

empty_row_delta <- function() {
  data.frame(
    scenario = character(),
    change_type = character(),
    Study = character(),
    T = character(),
    primary_N = numeric(),
    primary_R = numeric(),
    scenario_N = numeric(),
    scenario_R = numeric(),
    stringsAsFactors = FALSE
  )
}

scenario_row_delta <- function(primary, scenario, scenario_name) {
  primary_rows <- normalize_arm_rows(primary)
  scenario_rows <- normalize_arm_rows(scenario)
  merged <- merge(
    primary_rows, scenario_rows,
    by = c("Study", "T"), all = TRUE,
    suffixes = c("_primary", "_scenario"), sort = TRUE
  )
  changed <- is.na(merged$N_primary) | is.na(merged$N_scenario) |
    merged$N_primary != merged$N_scenario | merged$R_primary != merged$R_scenario
  merged <- merged[changed, , drop = FALSE]
  if (nrow(merged) == 0L) return(empty_row_delta())
  change_type <- ifelse(
    is.na(merged$N_primary), "added",
    ifelse(is.na(merged$N_scenario), "removed", "modified")
  )
  data.frame(
    scenario = scenario_name,
    change_type = change_type,
    Study = merged$Study,
    T = merged$T,
    primary_N = merged$N_primary,
    primary_R = merged$R_primary,
    scenario_N = merged$N_scenario,
    scenario_R = merged$R_scenario,
    stringsAsFactors = FALSE
  )
}

validation_row <- function(scenario, check, passed, detail) {
  data.frame(
    scenario = scenario,
    check = check,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}
source("R/publication_utils.R")
