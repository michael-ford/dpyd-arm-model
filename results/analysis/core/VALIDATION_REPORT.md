# Fixed-study pre-MCMC validation report

Status: **PASS**

Scenarios: 12; checks: 86; passed: 86.

## Primary network geometry

- Treatments: `1,2,3,4,5`.
- Connected: `TRUE`.
- Unique designs: `8`.
- Direct comparisons: `10`.
- Enumerated elementary closed loops: `37`.
- Design contrasts `m`: `18`; rank(A): `4`; inconsistency df `q`: `14`.
- Independent arm-level rank increase: `14` (expected `14`).

## Scenario inputs

| Scenario | Source strata | Source arms | Source N | Source events | Raw SHA-256 |
|---|---:|---:|---:|---:|---|
| `primary` | 25 | 79 | 12486 | 3821 | `bae7903764c071945e3910c56f25f36a85e98d7e6e46e694d1aa05448fa2680a` |
| `sensitivity_exclude_cevik` | 24 | 77 | 12383 | 3795 | `f684daafb71bb91a606fa2a631e513c5518dcc78c788368fe1400df0351263da` |
| `sensitivity_include_iachetta_kleibl` | 27 | 85 | 12982 | 4047 | `30c9be0b8000930a94e5d3bd1858ec43c2aecd9e9b45fafcac3323b2957e78b8` |
| `sensitivity_include_boisdron_overlap` | 26 | 82 | 12658 | 3836 | `07cf155881cf41d39c253a44407dcaac6e6989005b1bd64dc082cb8f3db17ba7` |
| `sensitivity_exclude_onesti` | 24 | 77 | 12360 | 3793 | `e7adeac9d2db171dd14a4f21d9bb54fdd8cf7c688076effc1152453f157de918` |
| `sensitivity_exclude_genotype_informed_carriers` | 25 | 71 | 12395 | 3778 | `b02d8dbac427e4d0c24a660bcd46d1fc19ec78aeaa113a5391affce6482fec74` |
| `sensitivity_exclude_pk_informed_morel` | 24 | 75 | 12000 | 3778 | `438f58fb5f3426b03cff0a86bb8fc2f0c6ce47fc362315bd37bf1713b6eede10` |
| `sensitivity_exclude_all_treatment_moderated` | 24 | 67 | 11909 | 3735 | `4dedea050d29585d5b4dfe95ad64e19dac07484fa754b51d00331bb146d32379` |
| `sensitivity_wigle_early` | 25 | 79 | 12486 | 3682 | `0135729877e18a4adaa8125fafc2acfd3203ebd4ca26f8de662f28b7f201bb84` |
| `sensitivity_cremolini_first_four` | 25 | 79 | 12486 | 3771 | `67858e6b9c33be093fa0573d9bad853709fa6e1cd397c61a6df1d99b30cc13e8` |
| `sensitivity_compound_lower_function` | 25 | 79 | 12489 | 3824 | `7090e58a5a19db5b0b1be160fc892ab676fd931653931e1502d1b3a44f5f908d` |
| `sensitivity_compound_reciprocal` | 25 | 79 | 12489 | 3824 | `3437e81751fc14a15813c93456e5d1c1c4400b3a10ad9b654f4becbd9a8ac724` |

## Validation coverage

The checkpoint verifies all twelve approved source transformations, exact scenario totals and row deltas, source-specific primary rows, treatment mappings, unique study/treatment arms, integer binomial validity, fixed-study identifiability, presence of all five nodes, network connectivity, designs, direct comparisons, closed loops, two independent inconsistency-df calculations, immutable input checksums, deterministic seed scheduling, and input stability.

The genotype-informed and combined treatment-moderation source scenarios retain the prospective Wigle and Bozina standard-care comparators for provenance. Once their carrier arms are removed, those comparator rows are single-arm strata and therefore cannot identify a relative effect under fixed study effects. The transformed model inputs mechanically omit those two orphan strata; the manifest records both source-level and fitted-input totals.

MCMC is not authorized unless this report and every machine-readable check pass.
