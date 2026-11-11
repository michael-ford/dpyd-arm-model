#!/usr/bin/env Rscript

# Fixed-study numerical supplement tables.  This script consumes
# only the reference fixed-study primary/core reports plus the verified S9 source
# classifications; it does not run a model or alter any analysis input.

core_sha256 <- function(path) {
  value <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(value, "status")) && attr(value, "status") != 0L) {
    stop("sha256sum failed for ", path, call. = FALSE)
  }
  sub(" .*", "", value[[1L]])
}

core_stop_unless <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

core_read_csv <- function(path, label) {
  core_stop_unless(file.exists(path), paste0("Missing ", label, ": ", path))
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

core_write_csv <- function(data, path) {
  write.csv(data, path, row.names = FALSE, na = "")
}

core_md_table <- function(data, columns) {
  data <- data[, columns, drop = FALSE]
  display <- vapply(data, function(x) {
    value <- as.character(x)
    value[is.na(value)] <- ""
    gsub("\\|", "\\\\|", value)
  }, character(nrow(data)))
  if (is.null(dim(display))) display <- matrix(display, ncol = 1L)
  lines <- c(
    paste0("| ", paste(names(data), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(data)), collapse = " | "), " |")
  )
  c(lines, apply(display, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")))
}

build_core_tables <- function(project_root = getwd(), output_root = NULL) {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  core_root <- file.path(project_root, "results", "analysis", "core")
  if (is.null(output_root)) output_root <- file.path(project_root, "runs", "rebuilt", "tables", "core")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    primary_raw = file.path(project_root, "data", "core", "primary_raw.csv"),
    manifest = file.path(core_root, "scenario_input_manifest.csv"),
    run_index = file.path(core_root, "core_run_index.csv"),
    all_treatments = file.path(core_root, "all_treatment_posterior_summary.csv"),
    primary_pairwise = file.path(project_root, "results", "tables", "core", "table_s6_pairwise_or_95_cri.csv"),
    hapb3_sensitivities = file.path(core_root, "hapb3_sensitivity_summary.csv"),
    heterogeneity = file.path(core_root, "heterogeneity_fit_summary.csv"),
    diagnostics = file.path(core_root, "diagnostic_extrema.csv"),
    primary_run_metadata = file.path(core_root, "primary_run_metadata.csv"),
    evidence = file.path(project_root, "data", "reporting", "comparator_genotyping_evidence.csv")
  )
  lapply(names(paths), function(name) core_stop_unless(file.exists(paths[[name]]), paste0("Missing required ", name, " input.")))

  primary <- core_read_csv(paths$primary_raw, "fixed-study primary input")
  manifest <- core_read_csv(paths$manifest, "scenario manifest")
  run_index <- core_read_csv(paths$run_index, "core run index")
  all_treatments <- core_read_csv(paths$all_treatments, "treatment posterior summary")
  ordered_pairwise <- core_read_csv(paths$primary_pairwise, "full-precision pairwise reference")
  labels <- c("WT", "HapB3", "c.2846A>T", "*2A", "*13")
  internal <- c("WT", "HapB3", "2846hetho", "2Ahetho", "13hetho")
  keep <- match(ordered_pairwise$numerator, labels) > match(ordered_pairwise$denominator, labels)
  primary_pairwise_raw <- data.frame(
    numerator = internal[match(ordered_pairwise$numerator[keep], labels)],
    denominator = internal[match(ordered_pairwise$denominator[keep], labels)],
    or_median = ordered_pairwise$OR[keep], or_lower_95_cri = ordered_pairwise$lower_95_CrI[keep],
    or_upper_95_cri = ordered_pairwise$upper_95_CrI[keep], stringsAsFactors = FALSE)
  hapb3_sensitivities <- core_read_csv(paths$hapb3_sensitivities, "HapB3 sensitivity summary")
  heterogeneity <- core_read_csv(paths$heterogeneity, "heterogeneity summary")
  diagnostics <- core_read_csv(paths$diagnostics, "diagnostic summary")
  primary_run_metadata <- core_read_csv(paths$primary_run_metadata, "primary run metadata")
  evidence <- core_read_csv(paths$evidence, "S9 source classifications")

  core_stop_unless(nrow(manifest) == 12L && nrow(run_index) == 12L,
                            "Expected exactly 12 reference core scenarios.")
  core_stop_unless(all(run_index$status == "complete") && all(run_index$converged) &&
                              all(run_index$all_diagnostics_passed) && all(run_index$artifact_checksums_match),
                            "Core run index is not fully complete, converged, and checksum-verified.")
  primary_manifest <- manifest[manifest$scenario == "primary", , drop = FALSE]
  core_stop_unless(nrow(primary_manifest) == 1L &&
                              identical(core_sha256(paths$primary_raw), primary_manifest$raw_input_sha256),
                            "Primary input does not reconcile to the immutable manifest checksum.")

  primary_effects <- all_treatments[
    all_treatments$scenario == "primary" & all_treatments$reference == "WT",
    c("treatment_id", "treatment", "reference", "or_median", "or_lower_95_cri", "or_upper_95_cri"), drop = FALSE
  ]
  primary_effects <- primary_effects[order(primary_effects$treatment_id), , drop = FALSE]
  core_stop_unless(nrow(primary_effects) == 5L &&
                              identical(primary_effects$treatment, c("WT", "HapB3", "2846hetho", "2Ahetho", "13hetho")),
                            "Primary posterior table must contain the locked five-node network.")
  names(primary_effects) <- c("node_id", "node", "reference", "OR", "lower_95_CrI", "upper_95_CrI")

  node_labels <- c(WT = "WT", HapB3 = "HapB3", `2846hetho` = "c.2846A>T", `2Ahetho` = "*2A", `13hetho` = "*13")
  core_stop_unless(nrow(primary_pairwise_raw) == 10L &&
                              setequal(primary_pairwise_raw$numerator, names(node_labels)[-1L]) &&
                              all(primary_pairwise_raw$denominator %in% names(node_labels)),
                            "Primary pairwise result must contain all ten unique five-node comparisons.")
  pairwise_forward <- data.frame(
    numerator = unname(node_labels[primary_pairwise_raw$numerator]),
    denominator = unname(node_labels[primary_pairwise_raw$denominator]),
    OR = primary_pairwise_raw$or_median,
    lower_95_CrI = primary_pairwise_raw$or_lower_95_cri,
    upper_95_CrI = primary_pairwise_raw$or_upper_95_cri,
    stringsAsFactors = FALSE
  )
  pairwise_reverse <- data.frame(
    numerator = pairwise_forward$denominator,
    denominator = pairwise_forward$numerator,
    OR = 1 / pairwise_forward$OR,
    lower_95_CrI = 1 / pairwise_forward$upper_95_CrI,
    upper_95_CrI = 1 / pairwise_forward$lower_95_CrI,
    stringsAsFactors = FALSE
  )
  pairwise_table <- rbind(pairwise_forward, pairwise_reverse)
  pairwise_table <- pairwise_table[order(match(pairwise_table$denominator, unname(node_labels)),
                                           match(pairwise_table$numerator, unname(node_labels))), , drop = FALSE]
  rownames(pairwise_table) <- NULL
  core_stop_unless(nrow(pairwise_table) == 20L &&
                              !anyDuplicated(pairwise_table[c("numerator", "denominator")]),
                            "Ordered pairwise Table S6 must contain 20 unique non-reference comparisons.")

  node_contributions <- do.call(rbind, lapply(names(node_labels), function(node) {
    source_node <- if (node == "WT") c("WT_Clean", "WT_Biased") else if (node == "HapB3") "HapB3_1129or1236" else node
    rows <- primary[primary$T %in% source_node, , drop = FALSE]
    data.frame(
      node = unname(node_labels[[node]]),
      analytical_strata = nrow(rows),
      participants = sum(rows$N),
      events = sum(rows$R),
      stringsAsFactors = FALSE
    )
  }))
  core_stop_unless(sum(node_contributions$analytical_strata) == 79L &&
                              sum(node_contributions$participants) == 12486L &&
                              sum(node_contributions$events) == 3821L,
                            "Node contribution table does not reconcile to the primary manifest.")

  sensitivity_table <- hapb3_sensitivities[, c(
    "scenario", "definition", "studies", "arms", "total_n", "total_r",
    "or_median", "or_lower_95_cri", "or_upper_95_cri", "delta_or_median_from_primary"
  ), drop = FALSE]
  names(sensitivity_table) <- c(
    "scenario", "definition", "studies", "arms", "participants", "events",
    "HapB3_OR", "lower_95_CrI", "upper_95_CrI", "OR_difference_from_primary"
  )
  core_stop_unless(nrow(sensitivity_table) == 12L &&
                              sensitivity_table$scenario[[1L]] == "primary",
                            "Expected the primary plus 11 bounded HapB3 sensitivities.")

  heterogeneity_table <- merge(
    heterogeneity[, c("scenario", "parameter", "heterogeneity_median", "heterogeneity_lower_95_cri", "heterogeneity_upper_95_cri",
                       "residual_deviance_mean", "per_arm_mean", "DIC", "descriptive_only"), drop = FALSE],
    diagnostics[, c("scenario", "max_psrf_upper_95", "min_effective_size", "diagnostics_passed"), drop = FALSE],
    by = "scenario", all.x = TRUE, sort = FALSE
  )
  heterogeneity_table <- heterogeneity_table[order(match(heterogeneity_table$scenario, sensitivity_table$scenario), heterogeneity_table$parameter), , drop = FALSE]
  names(heterogeneity_table) <- c(
    "scenario", "parameter", "median", "lower_95_CrI", "upper_95_CrI",
    "residual_deviance_mean", "residual_deviance_per_arm", "DIC_descriptive",
    "DIC_descriptive_only", "max_PSRF_upper_95", "min_ESS", "diagnostics_passed"
  )
  core_stop_unless(nrow(heterogeneity_table) == 24L &&
                              setequal(heterogeneity_table$parameter, c("tau_arm", "tau_contrast")),
                            "Expected tau_arm and tau_contrast for every core scenario.")

  primary_heterogeneity <- heterogeneity_table[
    heterogeneity_table$scenario == "primary", , drop = FALSE
  ]
  primary_manifest <- manifest[manifest$scenario == "primary", , drop = FALSE]
  primary_metadata <- primary_run_metadata[
    primary_run_metadata$scenario == "primary" & primary_run_metadata$status == "complete", , drop = FALSE
  ]
  core_stop_unless(nrow(primary_heterogeneity) == 2L && nrow(primary_manifest) == 1L && nrow(primary_metadata) == 1L,
                            "Primary manuscript diagnostics are incomplete.")
  primary_tau_arm <- primary_heterogeneity[primary_heterogeneity$parameter == "tau_arm", , drop = FALSE]
  primary_tau_contrast <- primary_heterogeneity[primary_heterogeneity$parameter == "tau_contrast", , drop = FALSE]
  table_s4 <- data.frame(
    metric = c(
      "Analytical strata", "Analyzed arms", "Residual deviance", "Residual deviance per arm",
      "Descriptive DIC", "Maximum PSRF upper 95%", "Minimum ESS", "Chains",
      "Adaptation iterations", "Burn-in iterations", "Sampling iterations per chain",
      "Thinning interval", "Retained posterior draws"
    ),
    value = c(
      primary_manifest$studies, primary_manifest$arms,
      primary_tau_arm$residual_deviance_mean, primary_tau_arm$residual_deviance_per_arm,
      primary_tau_arm$DIC_descriptive, primary_tau_arm$max_PSRF_upper_95, primary_tau_arm$min_ESS,
      primary_metadata$n_chains, primary_metadata$n_adapt, primary_metadata$n_burnin,
      primary_metadata$n_iter, primary_metadata$n_thin, primary_metadata$posterior_draws
    ),
    note = c(
      "Fixed-study primary", "Binomial arm observations", "Posterior mean", "Residual deviance / analyzed arms",
      "Descriptive only; not a competing-model selection test", "Pass threshold <= 1.05",
      "Pass threshold >= 1,000 for reported effects and heterogeneity", rep("", 6L)
    ),
    stringsAsFactors = FALSE
  )
  table_s5 <- data.frame(
    parameter = c("Arm-level tau", "Contrast-level sqrt(2) x tau"),
    median = c(primary_tau_arm$median, primary_tau_contrast$median),
    lower_95_CrI = c(primary_tau_arm$lower_95_CrI, primary_tau_contrast$lower_95_CrI),
    upper_95_CrI = c(primary_tau_arm$upper_95_CrI, primary_tau_contrast$upper_95_CrI),
    interpretation = c(
      "Common SD of independent arm-level deviations on the log-odds scale",
      "Implied SD of a contrast between two independent arm deviations"
    ),
    stringsAsFactors = FALSE
  )

  arm_value <- function(study, treatment) {
    row <- primary[primary$Study == study & primary$T == treatment, , drop = FALSE]
    core_stop_unless(nrow(row) <= 1L, paste0("Duplicate arm: ", study, " / ", treatment))
    if (!nrow(row)) return("")
    paste0(row$N, "/", row$R)
  }
  studies <- unique(primary$Study)
  s3a <- do.call(rbind, lapply(studies, function(study) {
    comparator <- primary[primary$Study == study & primary$T %in% c("WT_Clean", "WT_Biased"), , drop = FALSE]
    data.frame(
      Study = study,
      Comparator = if (nrow(comparator)) paste0(comparator$T, " ", comparator$N, "/", comparator$R) else "None",
      HapB3 = arm_value(study, "HapB3_1129or1236"),
      `*2A` = arm_value(study, "2Ahetho"),
      `c.2846A>T` = arm_value(study, "2846hetho"),
      `*13` = arm_value(study, "13hetho"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }))
  core_stop_unless(nrow(s3a) == 25L && !any(s3a$Study == "Boisdron_Celle_2007"),
                            "S3a must contain the 25 primary strata and no Boisdron-Celle stratum.")

  wt <- primary[primary$T %in% c("WT_Clean", "WT_Biased"), c("Study", "T", "N", "R"), drop = FALSE]
  core_stop_unless(nrow(wt) == 24L && !anyDuplicated(wt$Study) &&
                              !any(wt$Study == "Boisdron_Celle_2007"),
                            "Expected 24 unique fixed-study comparator arms excluding Boisdron-Celle.")
  evidence <- evidence[evidence$Study %in% wt$Study, , drop = FALSE]
  evidence <- evidence[match(wt$Study, evidence$Study), , drop = FALSE]
  core_stop_unless(nrow(evidence) == 24L && !anyNA(evidence$Study),
                            "S9 evidence must map one-to-one to the 24 primary comparator arms.")

  reconciliation <- data.frame(
    Study = wt$Study, source_WT_T = evidence$WT_T, source_WT_N = evidence$WT_N, source_WT_R = evidence$WT_R,
    primary_WT_T = wt$T, primary_WT_N = wt$N, primary_WT_R = wt$R,
    reconciled = evidence$WT_T == wt$T,
    count_basis = ifelse(evidence$WT_N == wt$N & evidence$WT_R == wt$R,
                         "source_and_primary_agree", "primary_count_reconciled_to_source"),
    stringsAsFactors = FALSE
  )
  core_stop_unless(all(reconciliation$reconciled) &&
                              all(reconciliation$primary_WT_N == wt$N) && all(reconciliation$primary_WT_R == wt$R),
                            "Source classifications do not reconcile to primary comparator arms.")
  core_stop_unless(
    identical(as.integer(reconciliation$primary_WT_N[reconciliation$Study == "Etienne_Grimaldi_2017"]), 231L) &&
      identical(as.integer(reconciliation$primary_WT_R[reconciliation$Study == "Etienne_Grimaldi_2017"]), 25L) &&
      identical(as.integer(reconciliation$primary_WT_N[reconciliation$Study == "Largillier_2006"]), 83L) &&
      identical(as.integer(reconciliation$primary_WT_R[reconciliation$Study == "Largillier_2006"]), 14L),
    "Locked Etienne or Largillier comparator denominator is incorrect."
  )

  variants <- data.frame(
    variant = c("HapB3", "*2A", "c.2846A>T", "*13"),
    treatment = c("HapB3", "2Ahetho", "2846hetho", "13hetho"),
    category_column = c("HapB3_category", "2A_category", "2846_category", "13_category"),
    residual_column = c("HapB3_residual_fraction", "2A_residual_fraction", "2846_residual_fraction", "13_residual_fraction"),
    carrier_frequency = c(0.047, 0.016, 0.007, 0.001),
    stringsAsFactors = FALSE
  )
  lookup <- primary_effects[match(variants$treatment, primary_effects$node), c("OR", "lower_95_CrI", "upper_95_CrI"), drop = FALSE]
  core_stop_unless(!anyNA(lookup$OR), "Primary OR/CrI lookup is incomplete for QBA.")
  variants <- cbind(variants, lookup)
  allowed_categories <- c("complete_target_negative", "partial", "prescreened", "incidental_post_hoc", "no_relevant_testing")
  detail <- do.call(rbind, lapply(seq_len(nrow(evidence)), function(i) {
    do.call(rbind, lapply(seq_len(nrow(variants)), function(j) {
      category <- evidence[[variants$category_column[[j]]]][[i]]
      residual <- as.numeric(evidence[[variants$residual_column[[j]]]][[i]])
      core_stop_unless(category %in% allowed_categories && is.finite(residual) && residual >= 0 && residual <= 1,
                                paste0("Invalid S9 classification for ", evidence$Study[[i]], "."))
      data.frame(
        Study = wt$Study[[i]], WT_T = wt$T[[i]], WT_N = wt$N[[i]], WT_R = wt$R[[i]],
        variant = variants$variant[[j]], category = category, residual_fraction = residual,
        carrier_frequency = variants$carrier_frequency[[j]], OR = variants$OR[[j]],
        lower_95_CrI = variants$lower_95_CrI[[j]], upper_95_CrI = variants$upper_95_CrI[[j]],
        source_kind = evidence$source_kind[[i]], source_id_or_file = evidence$source_id_or_file[[i]],
        source_location = evidence$source_location[[i]], evidence_note = evidence$evidence_note[[i]],
        uncertainty_note = evidence$uncertainty_note[[i]], stringsAsFactors = FALSE, check.names = FALSE
      )
    }))
  }))
  core_stop_unless(nrow(detail) == 96L && !any(detail$Study == "Boisdron_Celle_2007"),
                            "S9 must contain 24 comparator arms x four targets (96 cells), without Boisdron-Celle.")

  baseline_rate <- sum(wt$R) / sum(wt$N)
  detail$implied_variant_event_rate <- with(detail, OR * baseline_rate / (1 - baseline_rate + OR * baseline_rate))
  detail$expected_misclassified <- with(detail, WT_N * carrier_frequency * residual_fraction)
  detail$expected_upward_event_rate_shift <- with(
    detail, carrier_frequency * residual_fraction * (implied_variant_event_rate - baseline_rate)
  )
  qba <- do.call(rbind, lapply(split(detail, detail$Study), function(group) {
    data.frame(
      Study = group$Study[[1L]], WT_T = group$WT_T[[1L]], WT_N = group$WT_N[[1L]], WT_R = group$WT_R[[1L]],
      expected_misclassified = sum(group$expected_misclassified), pct_of_WT_arm = 100 * sum(group$expected_misclassified) / group$WT_N[[1L]],
      expected_upward_event_rate_shift = sum(group$expected_upward_event_rate_shift),
      stringsAsFactors = FALSE
    )
  }))
  qba <- qba[match(wt$Study, qba$Study), , drop = FALSE]

  changes <- rbind(
    data.frame(artifact = "S3a", Study = s3a$Study, status = "present",
               change = "Eligible primary analytical stratum", stringsAsFactors = FALSE),
    data.frame(artifact = "S3a/S9", Study = "Boisdron_Celle_2007", status = "excluded",
               change = "Probable overlap with retained Morel; included only in the overlap sensitivity", stringsAsFactors = FALSE)
  )

  outputs <- list(
    primary = file.path(output_root, "table_primary_five_node_or_95_cri.csv"),
    network_contributions = file.path(output_root, "table_primary_network_contributions.csv"),
    pairwise = file.path(output_root, "table_s6_pairwise_or_95_cri.csv"),
    sensitivity = file.path(output_root, "table_core_hapb3_sensitivities.csv"),
    heterogeneity = file.path(output_root, "table_tau_fit_diagnostics.csv"),
    s4 = file.path(output_root, "table_s4_primary_model_diagnostics.csv"),
    s5 = file.path(output_root, "table_s5_common_heterogeneity.csv"),
    s3a = file.path(output_root, "table_s3a_arm_contributions.csv"),
    s9 = file.path(output_root, "table_s9_comparator_genotyping_coverage.csv"),
    qba = file.path(output_root, "table_s9_qba_recalculated.csv"),
    qba_lookup = file.path(output_root, "qba_primary_or_cri_lookup.csv"),
    reconciliation = file.path(output_root, "source_primary_reconciliation.csv"),
    changes = file.path(output_root, "source_eligibility.csv"),
    report = file.path(output_root, "CORE_TABLES.md"),
    checksums = file.path(output_root, "provenance_input_checksums.csv")
  )
  core_write_csv(primary_effects, outputs$primary)
  core_write_csv(node_contributions, outputs$network_contributions)
  core_write_csv(pairwise_table, outputs$pairwise)
  core_write_csv(sensitivity_table, outputs$sensitivity)
  core_write_csv(heterogeneity_table, outputs$heterogeneity)
  core_write_csv(table_s4, outputs$s4)
  core_write_csv(table_s5, outputs$s5)
  core_write_csv(s3a, outputs$s3a)
  core_write_csv(detail, outputs$s9)
  core_write_csv(qba, outputs$qba)
  core_write_csv(variants[, c("variant", "treatment", "carrier_frequency", "OR", "lower_95_CrI", "upper_95_CrI")], outputs$qba_lookup)
  core_write_csv(reconciliation, outputs$reconciliation)
  core_write_csv(changes, outputs$changes)

  report <- c(
    "# Fixed-study numerical core supplement tables", "",
    "## Primary five-node odds ratios", "",
    core_md_table(primary_effects, names(primary_effects)), "",
    "## Primary network contributions", "",
    core_md_table(node_contributions, names(node_contributions)), "",
    "The five node rows reconcile to 79 analytical arms, 12,486 arm-level participants, and 3,821 events across 25 study strata.", "",
    "## Table S6 ordered pairwise odds ratios", "",
    core_md_table(pairwise_table, names(pairwise_table)), "",
    "## HapB3 primary and 11 bounded sensitivities", "",
    core_md_table(sensitivity_table, c("scenario", "participants", "events", "HapB3_OR", "lower_95_CrI", "upper_95_CrI")), "",
    "## Heterogeneity, fit, and diagnostics", "",
    core_md_table(heterogeneity_table, c("scenario", "parameter", "median", "lower_95_CrI", "upper_95_CrI", "residual_deviance_mean", "residual_deviance_per_arm", "DIC_descriptive", "max_PSRF_upper_95", "min_ESS")), "",
    "## Manuscript-ready Table S4", "",
    core_md_table(table_s4, names(table_s4)), "",
    "## Manuscript-ready Table S5", "",
    core_md_table(table_s5, names(table_s5)), "",
    "## S3a and S9/QBA scope", "",
    paste0("S3a contains ", nrow(s3a), " primary analytical strata. S9 contains ", nrow(detail), " target-specific coverage cells across ", nrow(qba), " comparator arms; Boisdron-Celle is excluded."),
    paste0("QBA uses the primary median OR lookup with its 95% CrI recorded in `qba_primary_or_cri_lookup.csv`; the expected-carrier arithmetic remains a screening calculation."),
    "",
    "## Source eligibility", "",
    "The source eligibility table lists primary strata and the excluded overlapping cohort.",
    "",
    "## Provenance", "",
    "`provenance_input_checksums.csv` records the immutable manifest, reference core reports, primary input, source classifications."
  )
  writeLines(report, outputs$report)

  provenance_inputs <- c(unlist(paths, use.names = TRUE), builder = file.path(project_root, "R", "build_core_tables.R"))
  provenance <- data.frame(
    input = names(provenance_inputs), path = substring(normalizePath(unname(provenance_inputs), winslash = "/", mustWork = TRUE),
                     nchar(project_root) + 2L),
    bytes = file.info(unname(provenance_inputs))$size,
    sha256 = vapply(unname(provenance_inputs), core_sha256, character(1)),
    stringsAsFactors = FALSE
  )
  core_write_csv(provenance, outputs$checksums)

  invisible(list(outputs = outputs, primary = primary_effects, network_contributions = node_contributions,
                 pairwise = pairwise_table, sensitivities = sensitivity_table,
                 heterogeneity = heterogeneity_table, s4 = table_s4, s5 = table_s5,
                 s3a = s3a, s9 = detail, qba = qba,
                 reconciliation = reconciliation, changes = changes, provenance = provenance))
}

if (sys.nframe() == 0L) {
  result <- build_core_tables()
  cat(sprintf("Built fixed-study core tables: 5-node primary, %d sensitivities, %d S3a strata, and %d S9 cells.\n",
              nrow(result$sensitivities), nrow(result$s3a), nrow(result$s9)))
}
