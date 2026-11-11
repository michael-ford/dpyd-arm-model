# DPYD fixed-study network meta-analysis

This repository reproduces the final arm-level analysis of paper-defined severe fluoropyrimidine toxicity: 25 analytical strata, 79 mutually exclusive arms, 12,486 participants, and 3,821 events. The Bayesian binomial-logit model uses fixed study effects, independent deviations on every arm, common arm-level heterogeneity, and wild type (WT) as the global reference. Genotype effects use concurrent within-stratum information (White et al., 2019; Piepho et al., 2024).

The primary odds ratios versus WT, from `results/tables/core/table_primary_five_node_or_95_cri.csv`, are:

| Variant | OR | 95% credible interval |
|---|---:|---:|
| HapB3 | 1.73 | 1.00–2.91 |
| c.2846A>T | 7.30 | 3.55–16.14 |
| *2A | 8.20 | 4.42–16.01 |
| *13 | 2.81 | 0.80–10.08 |

The HapB3 lower bound before rounding is 1.000412; several sensitivity intervals include 1. Only relative odds ratios are reported: there are no clinical absolute-risk estimates or treatment rankings. The comparator-classification QBA is a screening calculation, not corrected patient data. The global inconsistency diagnostic converged but failed its information gate; inconsistency could not be assessed reliably.

## Layout

- `R/fixed_study_model.R`, `R/model_fixed_study/`: model and JAGS definitions.
- `R/design_by_treatment.R`, `R/network_geometry.R`: inconsistency diagnostic and network geometry.
- `R/core_*`, `R/dependent_*`, `R/inconsistency_runner.R`: scenario preparation, checkpoints, and execution.
- `R/arm_rows.R`, `R/scenario_helpers.R`, `R/checksums.R`, `R/common/data_utils.R`: shared row, integrity, and transformation helpers.
- `R/build_core_tables.R`, `R/build_sensitivity_tables.R`, `R/render_core_figures.R`, `R/render_sensitivity_figures.R`: publication assets.
- `R/run_scenario.R`, `R/compare_to_reference.R`, `R/rebuild_outputs.R`: public commands.
- `R/generate_core_inputs.R`, `R/generate_dependent_inputs.R`: optional checkpoint generation for inspection under `runs/checkpoints/`.
- `R/summarize_dependent.R`: internal checkpoint helpers; it expects checkpoint-local attempts and is not a public command for aggregating repository-level reruns.
- `tests/`: focused validation and model tests.
- `data/core/`, `data/dependent/`: exact scenario inputs; `data/reporting/`: published comparator classifications.
- `data/input_checksums.csv`: SHA-256 manifest of public inputs, code, and reference assets.
- `results/`: reference tables, figures, reports, diagnostics, scenario specifications, and both committed seed schedules under `results/analysis/`.
- `runs/`: ignored reruns and rebuilt assets; existing attempts are preserved.

## Reproduction

Run from the repository root with Docker installed:

```sh
docker build --build-arg GIT_COMMIT="$(git rev-parse HEAD)" -t dpyd-fixed-study .

# Validate every runnable scenario without MCMC.
docker run --rm -v "$PWD:/analysis" -w /analysis dpyd-fixed-study Rscript R/run_scenario.R --all --preflight

# Primary analysis.
docker run --rm -v "$PWD:/analysis" -w /analysis dpyd-fixed-study Rscript R/run_scenario.R primary

# All 12 core fits, 29 dependent fits, and two inconsistency diagnostics.
# Use a separate output directory/clean checkout if the primary was already run;
# or record a reason for another preserved attempt:
docker run --rm -e RERUN_REASON="Full reproduction" -v "$PWD:/analysis" -w /analysis dpyd-fixed-study Rscript R/run_scenario.R --all

# Rebuild tables and figures from the committed reference summaries.
docker run --rm -v "$PWD:/analysis" -w /analysis dpyd-fixed-study Rscript R/rebuild_outputs.R

# Compare the latest primary rerun with the committed reference.
docker run --rm -v "$PWD:/analysis" -w /analysis dpyd-fixed-study Rscript R/compare_to_reference.R primary

# Focused tests, each in a fresh R process.
docker run --rm -v "$PWD:/analysis" -w /analysis dpyd-fixed-study sh -c 'for f in tests/test_*.R; do Rscript "$f" || exit 1; done'
```

`--list` lists runnable scenario names. Three grouped exclusions alias core fits; the N≤5 sparse scenario is explicitly infeasible because it removes *13. Outputs go to `runs/<scenario>/attempt_NNN/`; rebuilt assets go to `runs/rebuilt/`. A nonconverged fit stops without extension or fallback. An inconsistency information-gate failure is retained and reported as non-estimable.

Every production fit uses four chains, 50,000 adaptation iterations, 125,000 burn-in iterations, 250,000 sampling iterations per chain, and thinning 10 (100,000 retained draws). The master seed is 12345; scenario-specific chain seeds come from the committed schedules. Every reported effect, pairwise contrast, and heterogeneity parameter must have PSRF upper 95% ≤1.05 and ESS ≥1,000. Run metadata records input hashes, Git commit, R, JAGS, rjags, and coda versions.

The recorded dependent fits took about 3 minutes to 7.4 hours each; WT-split fits were the longest. Hardware and concurrent execution affect these times. No primary-specific timing was included in the reference package. Allow many hours for all scenarios.

| Runtime component | Version or build policy |
|---|---|
| R | 4.4.1 (`rocker/r-ver:4.4.1`) |
| JAGS | 4.3.0; installed through APT and asserted at build time |
| rjags | 4-17 (displayed by R as 4.17) |
| coda | 0.19-4.1 (displayed by R as 0.19.4.1) |
| dplyr | CRAN at image build; arm transformations |
| ggplot2 | CRAN at image build; figures |

The build installs normal transitive dependencies. `remotes` is an installation helper and is removed afterward. Supporting-package and operating-system revisions were not recorded in the reference execution log and are not pinned here. Retain the verified Docker image for exact runtime reproduction; a future rebuild may resolve newer supporting packages and change plot appearance.

The comparison prints maximum absolute differences in OR/CrI, tau, and diagnostic extrema. Core references include OR/CrI and tau; dependent references contain selected HapB3 contrasts and diagnostic extrema, but no tau summaries. Per-parameter reference PSRF/ESS are unavailable, so comparisons use maximum PSRF upper 95% and minimum ESS. Each rerun still records and gates every monitored parameter. Inconsistency references contain gate results and diagnostic extrema, not full OR/tau summaries.

Reproduction status of this release: the primary analysis was rerun in the release image and reproduced the reference results exactly. The other 40 fits were not rerun for this release; their reference results are carried over unchanged, and their inputs, seeds and model settings were verified identical to those of the runs that produced them.

Table and figure rebuilding uses committed lightweight reference summaries; it does not reconstruct those summaries from posterior draws. Rebuilt assets include repository-relative input provenance ledgers and source eligibility. Optional checkpoint generation is for comparison and inspection; it does not replace the committed inputs or release metadata.

This replaces the earlier pcnetmeta analysis, which remains in Git history. The final model is described in the manuscript and `results/METHODS.md`, with numerical reports under `results/analysis/`.

White et al. (2019): https://doi.org/10.1002/sim.8360

Piepho et al. (2024): https://doi.org/10.1002/jrsm.1721
