# Network geometry checks for the arm-level NMA input.
#
# This module deliberately uses base R only.  It is intended to be called before
# MCMC so that the evidence network, its independent designs, and the available
# inconsistency degrees of freedom are explicit and reproducible.

.network_geometry_validate <- function(data) {
  if (!is.data.frame(data)) {
    stop("Network input must be a data.frame.", call. = FALSE)
  }

  required <- c("Study", "t.id")
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Network input is missing required column(s): ",
      paste(missing_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (nrow(data) == 0L) {
    stop("Network input has no arms.", call. = FALSE)
  }

  study <- as.character(data$Study)
  treatment <- data$t.id
  if (anyNA(study) || any(!nzchar(trimws(study)))) {
    stop("Network input has a missing Study label.", call. = FALSE)
  }
  if (!is.numeric(treatment) || anyNA(treatment) || any(!is.finite(treatment))) {
    stop("Network input has a missing or non-numeric treatment node (t.id).", call. = FALSE)
  }
  if (any(treatment != floor(treatment)) || any(treatment < 1)) {
    stop("Treatment nodes (t.id) must be positive integers.", call. = FALSE)
  }
  treatment <- as.integer(treatment)

  normalized <- data.frame(
    Study = study,
    t.id = treatment,
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(normalized)) {
    stop("Network input has duplicate study-treatment arms.", call. = FALSE)
  }

  arms_per_study <- table(normalized$Study)
  single_arm <- names(arms_per_study)[arms_per_study < 2L]
  if (length(single_arm) > 0L) {
    stop(
      "Network input has single-arm study/studies: ",
      paste(sort(single_arm), collapse = ", "), ".",
      call. = FALSE
    )
  }

  normalized[order(normalized$Study, normalized$t.id), , drop = FALSE]
}

.network_geometry_study_sets <- function(data) {
  studies <- sort(unique(data$Study))
  treatment_sets <- lapply(studies, function(study) {
    sort(unique(data$t.id[data$Study == study]))
  })
  names(treatment_sets) <- studies
  treatment_sets
}

.network_geometry_rank <- function(x) {
  if (nrow(x) == 0L || ncol(x) == 0L) 0L else qr(x)$rank
}

#' Return treatment-network connectedness and its components.
#'
#' This helper is intentionally available separately from `analyze_network_geometry`
#' so callers can inspect components when the latter rejects a disconnected input.
network_connectedness <- function(data) {
  data <- .network_geometry_validate(data)
  treatment_sets <- .network_geometry_study_sets(data)
  nodes <- sort(unique(data$t.id))
  neighbours <- setNames(vector("list", length(nodes)), as.character(nodes))

  for (set in treatment_sets) {
    pairs <- utils::combn(set, 2L)
    for (i in seq_len(ncol(pairs))) {
      left <- as.character(pairs[1L, i])
      right <- as.character(pairs[2L, i])
      neighbours[[left]] <- sort(unique(c(neighbours[[left]], as.integer(right))))
      neighbours[[right]] <- sort(unique(c(neighbours[[right]], as.integer(left))))
    }
  }

  unseen <- nodes
  components <- list()
  while (length(unseen) > 0L) {
    queue <- unseen[1L]
    component <- integer()
    while (length(queue) > 0L) {
      current <- queue[1L]
      queue <- queue[-1L]
      if (current %in% component) next
      component <- c(component, current)
      queue <- c(queue, neighbours[[as.character(current)]])
      queue <- sort(unique(queue[!queue %in% component]))
    }
    component <- sort(component)
    components[[length(components) + 1L]] <- component
    unseen <- setdiff(unseen, component)
  }

  component_table <- do.call(rbind, lapply(seq_along(components), function(i) {
    data.frame(component_id = i, treatment = components[[i]])
  }))
  rownames(component_table) <- NULL

  list(
    connected = length(components) == 1L,
    n_components = length(components),
    components = components,
    component_table = component_table
  )
}

#' Identify unique treatment-set designs and their contributing studies.
network_designs <- function(data) {
  data <- .network_geometry_validate(data)
  treatment_sets <- .network_geometry_study_sets(data)
  design_keys <- vapply(treatment_sets, paste, character(1L), collapse = "-")
  unique_keys <- sort(unique(design_keys))

  rows <- lapply(seq_along(unique_keys), function(i) {
    key <- unique_keys[i]
    studies <- sort(names(treatment_sets)[design_keys == key])
    treatments <- treatment_sets[[studies[1L]]]
    baseline <- if (1L %in% treatments) 1L else min(treatments)
    data.frame(
      design_id = i,
      design = key,
      n_treatments = length(treatments),
      baseline = baseline,
      n_studies = length(studies),
      stringsAsFactors = FALSE
    )
  })
  designs <- do.call(rbind, rows)
  designs$treatments <- I(lapply(unique_keys, function(key) {
    treatment_sets[[names(treatment_sets)[match(key, design_keys)]]]
  }))
  designs$studies <- I(lapply(unique_keys, function(key) {
    sort(names(treatment_sets)[design_keys == key])
  }))
  designs
}

#' Return every direct treatment comparison and the studies contributing it.
network_direct_edges <- function(data) {
  data <- .network_geometry_validate(data)
  treatment_sets <- .network_geometry_study_sets(data)
  edge_rows <- lapply(names(treatment_sets), function(study) {
    pairs <- utils::combn(treatment_sets[[study]], 2L)
    data.frame(
      treatment1 = pairs[1L, ],
      treatment2 = pairs[2L, ],
      Study = study,
      stringsAsFactors = FALSE
    )
  })
  edge_rows <- do.call(rbind, edge_rows)
  edge_key <- paste(edge_rows$treatment1, edge_rows$treatment2, sep = "-")
  keys <- sort(unique(edge_key))

  edges <- do.call(rbind, lapply(keys, function(key) {
    rows <- edge_rows[edge_key == key, , drop = FALSE]
    data.frame(
      treatment1 = rows$treatment1[1L],
      treatment2 = rows$treatment2[1L],
      n_studies = length(unique(rows$Study)),
      stringsAsFactors = FALSE
    )
  }))
  edges$studies <- I(lapply(keys, function(key) {
    sort(unique(edge_rows$Study[edge_key == key]))
  }))
  rownames(edges) <- NULL
  edges
}

.network_geometry_canonical_cycle <- function(cycle) {
  start <- min(cycle)
  forward <- c(cycle[match(start, cycle):length(cycle)], cycle[seq_len(match(start, cycle) - 1L)])
  reversed <- rev(cycle)
  reverse_start <- match(start, reversed)
  backward <- c(reversed[reverse_start:length(reversed)], reversed[seq_len(reverse_start - 1L)])
  if (paste(forward, collapse = "-") <= paste(backward, collapse = "-")) forward else backward
}

#' Enumerate elementary undirected cycles in canonical, non-duplicated form.
network_closed_loops <- function(data) {
  data <- .network_geometry_validate(data)
  nodes <- sort(unique(data$t.id))
  edges <- network_direct_edges(data)
  neighbours <- setNames(vector("list", length(nodes)), as.character(nodes))
  for (i in seq_len(nrow(edges))) {
    left <- as.character(edges$treatment1[i])
    right <- as.character(edges$treatment2[i])
    neighbours[[left]] <- sort(unique(c(neighbours[[left]], edges$treatment2[i])))
    neighbours[[right]] <- sort(unique(c(neighbours[[right]], edges$treatment1[i])))
  }

  found <- list()
  visit <- function(start, current, path) {
    for (next_node in neighbours[[as.character(current)]]) {
      if (next_node == start && length(path) >= 3L) {
        canonical <- .network_geometry_canonical_cycle(path)
        found[[paste(canonical, collapse = "-")]] <<- canonical
      } else if (!next_node %in% path && next_node > start) {
        visit(start, next_node, c(path, next_node))
      }
    }
  }
  for (start in nodes) visit(start, start, start)

  keys <- sort(names(found))
  loops <- unname(found[keys])
  if (length(loops) == 0L) {
    return(list(
      loops = loops,
      loop_table = data.frame(loop_id = integer(), length = integer(), cycle = character())
    ))
  }
  loop_table <- data.frame(
    loop_id = seq_along(loops),
    length = lengths(loops),
    cycle = vapply(loops, function(x) paste(c(x, x[1L]), collapse = "-"), character(1L)),
    stringsAsFactors = FALSE
  )
  list(loops = loops, loop_table = loop_table)
}

.network_geometry_A <- function(designs, nodes, network_baseline) {
  treatment_columns <- setdiff(nodes, network_baseline)
  m <- sum(designs$n_treatments - 1L)
  A <- matrix(0, nrow = m, ncol = length(treatment_columns))
  colnames(A) <- paste0("treatment_", treatment_columns)
  mapping_rows <- vector("list", m)
  row <- 0L

  for (i in seq_len(nrow(designs))) {
    design_treatments <- designs$treatments[[i]]
    design_baseline <- designs$baseline[i]
    for (treatment in setdiff(design_treatments, design_baseline)) {
      row <- row + 1L
      if (treatment != network_baseline) {
        A[row, match(treatment, treatment_columns)] <- 1
      }
      if (design_baseline != network_baseline) {
        A[row, match(design_baseline, treatment_columns)] <- -1
      }
      mapping_rows[[row]] <- data.frame(
        row = row,
        design_id = designs$design_id[i],
        design = designs$design[i],
        design_baseline = design_baseline,
        treatment = treatment,
        contrast = paste0(treatment, " - ", design_baseline),
        stringsAsFactors = FALSE
      )
    }
  }
  rownames(A) <- vapply(mapping_rows, function(x) {
    paste0("design_", x$design_id, ":", x$contrast)
  }, character(1L))
  list(A = A, mapping = do.call(rbind, mapping_rows))
}

.network_geometry_arm_level_rank_check <- function(data, nodes, baseline, designs, A_mapping, q) {
  studies <- sort(unique(data$Study))
  study_matrix <- outer(data$Study, studies, `==`) * 1
  basic_treatments <- setdiff(nodes, baseline)
  treatment_matrix <- outer(data$t.id, basic_treatments, `==`) * 1
  rank_study <- .network_geometry_rank(study_matrix)
  rank_study_treatment <- .network_geometry_rank(cbind(study_matrix, treatment_matrix))
  rank_increase <- rank_study_treatment - rank_study

  # Z is the saturated, design-specific version of the independent contrast
  # rows in A.  Each column is one unique-design contrast: it is one for the
  # non-baseline arm of that design and zero for its design baseline.  Study
  # intercepts supply the corresponding local baselines.  Thus [study, X, Z]
  # adds exactly the fixed design-by-treatment directions beyond the consistent
  # treatment model [study, X].
  study_design_keys <- vapply(.network_geometry_study_sets(data), paste, character(1L), collapse = "-")
  design_by_study <- setNames(
    designs$design_id[match(study_design_keys, designs$design)],
    names(study_design_keys)
  )
  z_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(A_mapping))
  for (column in seq_len(nrow(A_mapping))) {
    matches <- design_by_study[data$Study] == A_mapping$design_id[column] &
      data$t.id == A_mapping$treatment[column]
    z_matrix[matches, column] <- 1
  }
  colnames(z_matrix) <- rownames(.network_geometry_A(designs, nodes, baseline)$A)

  rank_consistent <- rank_study_treatment
  rank_saturated <- .network_geometry_rank(cbind(study_matrix, treatment_matrix, z_matrix))
  inconsistency_rank_increase <- rank_saturated - rank_consistent

  list(
    treatment_rank_check = list(
      rank_study = rank_study,
      rank_study_treatment = rank_study_treatment,
      treatment_rank_increase = rank_increase,
      expected_rank_increase = length(nodes) - 1L,
      passes = identical(as.integer(rank_increase), as.integer(length(nodes) - 1L))
    ),
    rank_study = rank_study,
    rank_study_treatment = rank_study_treatment,
    treatment_rank_increase = rank_increase,
    expected_rank_increase = length(nodes) - 1L,
    rank_consistent = rank_consistent,
    rank_saturated = rank_saturated,
    inconsistency_rank_increase = inconsistency_rank_increase,
    expected_inconsistency_rank_increase = q,
    passes = identical(as.integer(inconsistency_rank_increase), as.integer(q))
  )
}

#' Analyse the geometry of a connected arm-level treatment network.
#'
#' `A` maps the independent contrast for each unique design (rows) to the
#' network's basic treatment effects relative to `baseline` (columns).  It is
#' the deterministic mapping required by a later fixed design-by-treatment
#' inconsistency extension.  Repeated studies with the same treatment set add
#' precision, but not rows to `A` or inconsistency degrees of freedom.
analyze_network_geometry <- function(data) {
  data <- .network_geometry_validate(data)
  connectedness <- network_connectedness(data)
  if (!connectedness$connected) {
    stop(
      "Network is disconnected; treatment components are ",
      paste(vapply(connectedness$components, paste, character(1L), collapse = ","), collapse = " | "), ".",
      call. = FALSE
    )
  }

  nodes <- sort(unique(data$t.id))
  baseline <- if (1L %in% nodes) 1L else min(nodes)
  designs <- network_designs(data)
  A_result <- .network_geometry_A(designs, nodes, baseline)
  m <- nrow(A_result$A)
  rank_A <- .network_geometry_rank(A_result$A)
  q <- m - (length(nodes) - 1L)
  arm_level_rank_check <- .network_geometry_arm_level_rank_check(
    data = data,
    nodes = nodes,
    baseline = baseline,
    designs = designs,
    A_mapping = A_result$mapping,
    q = q
  )
  loops <- network_closed_loops(data)

  list(
    treatments = nodes,
    baseline = baseline,
    connectedness = connectedness,
    designs = designs,
    direct_edges = network_direct_edges(data),
    loops = loops$loops,
    loop_table = loops$loop_table,
    m = m,
    A = A_result$A,
    A_mapping = A_result$mapping,
    rank_A = rank_A,
    q = q,
    q_from_rank = m - rank_A,
    arm_level_rank_check = arm_level_rank_check
  )
}

# Short alias for use in pre-MCMC gates.
network_geometry <- analyze_network_geometry
