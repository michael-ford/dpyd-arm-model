# Fixed-study inconsistency report

Both reference design-by-treatment fits completed the locked 100,000-draw MCMC contract and passed combined MCMC convergence diagnostics. They are nevertheless non-reportable because each failed the pre-specified information gate. The empirical omnibus output was calculable in each run, but is non-reportable; it must not be interpreted as evidence of no inconsistency.

Inconsistency could not be reliably estimated.

No node splitting was triggered because no reportable global conflict was identified.

## Verified execution contract

- Exactly one reference `attempt_001` was recorded for each prior scale; both carry `RUN_NONCONVERGED` solely because the information gate failed, not because MCMC failed.
- The manifested primary raw and transformed input checksums, deterministic four-chain seed schedules, production configuration (50,000 adaptation; 125,000 burn-in; 250,000 iterations; thinning 10), saved before/after code hashes, saved before/after protected-artifact hashes, and every indexed attempt artifact were revalidated.

## Information and diagnostics

| Prior SD | Max posterior/prior variance ratio | Threshold | Combined MCMC | Max PSRF upper 95% | Min ESS | Empirical omnibus | Residual deviance/arm | Descriptive DIC |
|---:|---:|---:|---|---:|---:|---|---:|---:|
| 2 | 0.640433 | 0.64 | passed | 1.002082 | 2942.6 | calculated but non-reportable (no_evidence_of_inconsistency) | 1.079797 | 369.607723 |
| 1 | 0.848845 | 0.64 | passed | 1.002020 | 3897.0 | calculated but non-reportable (no_evidence_of_inconsistency) | 1.081702 | 367.397762 |

The SD=2 ratio (0.640433) exceeds the 0.64 threshold by 0.000433; the SD=1 ratio is 0.848845. Thus the stable MCMC samples do not establish that the design-interaction directions are sufficiently data-informed for inferential reporting.

Machine-readable details: `inconsistency_run_gate_summary.csv`, `inconsistency_artifact_checksum_path_manifest.csv`, and `inconsistency_protected_artifact_snapshot_manifest.csv`.
