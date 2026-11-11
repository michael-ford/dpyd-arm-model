# Fixed-study core results verification

## Verification decision

**PASS — all 12 manifested core scenarios have exactly one `RUN_COMPLETE` attempt. No core attempt is failed or nonconverged.**

Every core attempt used the locked four-chain configuration with 100,000 retained draws (50,000 adaptation; 125,000 burn-in; 250,000 iterations; thinning 10). All recorded diagnostics passed.
The verifier rehashed 216 runtime-recorded attempt artifacts; each matched its recorded checksum. All manifested input checksums also matched.

## Prespecified primary network totals

Primary manifest totals: 25 studies, 79 arms, N=12486, and 3821 events.
Network totals: 8 designs, 10 direct comparisons, 37 enumerated elementary closed loops, and inconsistency df=14.

## Numerical stability

Across the 12 prespecified scenarios, HapB3-versus-WT posterior median ORs range from 1.6498 to 2.0337. Scenario definitions and primary-relative median deltas are provided without additional interpretation in `hapb3_sensitivity_summary.csv`.
Diagnostic extrema: maximum PSRF upper 95% = 1.004941; minimum effective sample size = 3268.0.

## Files

- `core_run_index.csv` — immutable-attempt and contract verification index.
- `all_treatment_posterior_summary.csv` — posterior summaries for every treatment in every core scenario.
- `hapb3_sensitivity_summary.csv` — HapB3-versus-WT summaries, scenario definitions, and deltas from primary.
- `heterogeneity_fit_summary.csv` and `diagnostic_extrema.csv` — numerical fit and convergence summaries.

These are reference summaries. Public reruns record their own checksums, diagnostics, and software versions under `runs/`. Posterior chains are not included in this repository.
