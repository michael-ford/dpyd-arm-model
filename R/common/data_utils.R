# Shared arm transformation and validation utilities.

library(dplyr)

transform_to_nma_format <- function(raw_data, treatment_map, aggregate_wt = FALSE) {
  # Validate treatments match the map
  data_treatments <- unique(raw_data$T)
  unmapped_treatments <- setdiff(data_treatments, names(treatment_map))
  if (length(unmapped_treatments) > 0) {
    stop("Treatments in data not found in treatment_map: ",
         paste(unmapped_treatments, collapse = ", "),
         "\nValid treatments: ", paste(names(treatment_map), collapse = ", "))
  }

  nma_data <- raw_data %>%
    mutate(
      s.id = as.numeric(factor(Study, levels = unique(Study))),
      t.id = treatment_map[T],
      r = R,
      n = N
    ) %>%
    select(s.id, t.id, r, n, Study, T) %>%
    arrange(s.id, t.id)

  # If aggregating WT (for unified model), combine arms with same treatment ID
  if (aggregate_wt) {
    nma_data <- nma_data %>%
      group_by(s.id, t.id) %>%
      summarise(
        r = sum(r),
        n = sum(n),
        Study = first(Study),
        T = first(T),
        .groups = "drop"
      ) %>%
      arrange(s.id, t.id)
  }

  return(nma_data)
}

# -----------------------------------------------------------------------------
# Validate transformed data
# -----------------------------------------------------------------------------

validate_nma_data <- function(nma_data) {
  cat("=== Data Validation ===\n")

  # Check for missing values
  if (any(is.na(nma_data$s.id)) || any(is.na(nma_data$t.id))) {
    stop("Missing study or treatment IDs after transformation")
  }

  # Check r <= n
  invalid_arms <- nma_data %>% filter(r > n)
  if (nrow(invalid_arms) > 0) {
    cat("ERROR: Events exceed total in these arms:\n")
    print(invalid_arms %>% select(Study, T, r, n))
    stop("Data integrity error: R > N in ", nrow(invalid_arms), " arm(s)")
  }

  cat("All validations passed\n\n")
  return(TRUE)
}
