source("R/build_core_tables.R")

project_root <- normalizePath(getwd(), mustWork = TRUE)
temporary_root <- tempfile("fixed-study-core-tables-")
dir.create(temporary_root, recursive = TRUE)
on.exit(unlink(temporary_root, recursive = TRUE), add = TRUE)

result <- build_core_tables(project_root, temporary_root)

required <- c(
  "table_primary_five_node_or_95_cri.csv", "table_core_hapb3_sensitivities.csv",
  "table_primary_network_contributions.csv", "table_s6_pairwise_or_95_cri.csv",
  "table_tau_fit_diagnostics.csv", "table_s4_primary_model_diagnostics.csv",
  "table_s5_common_heterogeneity.csv", "table_s3a_arm_contributions.csv",
  "table_s9_comparator_genotyping_coverage.csv", "table_s9_qba_recalculated.csv",
  "qba_primary_or_cri_lookup.csv", "source_primary_reconciliation.csv",
  "source_eligibility.csv", "CORE_TABLES.md", "provenance_input_checksums.csv"
)
stopifnot(all(file.exists(file.path(temporary_root, required))))

primary <- read.csv(file.path(temporary_root, "table_primary_five_node_or_95_cri.csv"), stringsAsFactors = FALSE)
network_contributions <- read.csv(file.path(temporary_root, "table_primary_network_contributions.csv"), stringsAsFactors = FALSE)
pairwise <- read.csv(file.path(temporary_root, "table_s6_pairwise_or_95_cri.csv"), stringsAsFactors = FALSE)
sensitivities <- read.csv(file.path(temporary_root, "table_core_hapb3_sensitivities.csv"), stringsAsFactors = FALSE)
heterogeneity <- read.csv(file.path(temporary_root, "table_tau_fit_diagnostics.csv"), stringsAsFactors = FALSE)
s4 <- read.csv(file.path(temporary_root, "table_s4_primary_model_diagnostics.csv"), stringsAsFactors = FALSE)
s5 <- read.csv(file.path(temporary_root, "table_s5_common_heterogeneity.csv"), stringsAsFactors = FALSE)
s3a <- read.csv(file.path(temporary_root, "table_s3a_arm_contributions.csv"), stringsAsFactors = FALSE, check.names = FALSE)
s9 <- read.csv(file.path(temporary_root, "table_s9_comparator_genotyping_coverage.csv"), stringsAsFactors = FALSE)
qba <- read.csv(file.path(temporary_root, "table_s9_qba_recalculated.csv"), stringsAsFactors = FALSE)
reconciliation <- read.csv(file.path(temporary_root, "source_primary_reconciliation.csv"), stringsAsFactors = FALSE)
changes <- read.csv(file.path(temporary_root, "source_eligibility.csv"), stringsAsFactors = FALSE)

stopifnot(
  nrow(primary) == 5L,
  identical(primary$node, c("WT", "HapB3", "2846hetho", "2Ahetho", "13hetho")),
  nrow(network_contributions) == 5L,
  sum(network_contributions$analytical_strata) == 79L,
  sum(network_contributions$participants) == 12486L,
  sum(network_contributions$events) == 3821L,
  nrow(pairwise) == 20L,
  !anyDuplicated(pairwise[c("numerator", "denominator")]),
  pairwise$OR[pairwise$numerator == "HapB3" & pairwise$denominator == "WT"] > 1.72,
  pairwise$OR[pairwise$numerator == "HapB3" & pairwise$denominator == "WT"] < 1.73,
  nrow(sensitivities) == 12L,
  nrow(heterogeneity) == 24L,
  setequal(heterogeneity$parameter, c("tau_arm", "tau_contrast")),
  all(heterogeneity$DIC_descriptive_only),
  nrow(s4) == 13L,
  s4$value[s4$metric == "Analyzed arms"] == 79,
  s4$value[s4$metric == "Chains"] == 4,
  nrow(s5) == 2L,
  identical(s5$parameter, c("Arm-level tau", "Contrast-level sqrt(2) x tau")),
  nrow(s3a) == 25L,
  !any(s3a$Study == "Boisdron_Celle_2007"),
  nrow(qba) == 24L,
  nrow(s9) == 96L,
  length(table(s9$Study)) == 24L,
  all(table(s9$Study) == 4L),
  !any(s9$Study == "Boisdron_Celle_2007"),
  all(reconciliation$reconciled),
  reconciliation$primary_WT_N[reconciliation$Study == "Etienne_Grimaldi_2017"] == 231L,
  reconciliation$primary_WT_R[reconciliation$Study == "Etienne_Grimaldi_2017"] == 25L,
  reconciliation$primary_WT_N[reconciliation$Study == "Largillier_2006"] == 83L,
  reconciliation$primary_WT_R[reconciliation$Study == "Largillier_2006"] == 14L,
  any(changes$Study == "Boisdron_Celle_2007" & changes$status == "excluded")
)

output_names <- list.files(temporary_root, full.names = FALSE)
stopifnot(!any(grepl("absolute|sucra|correlation", output_names, ignore.case = TRUE)))

cat("Fixed-study core supplement table tests passed.\n")
