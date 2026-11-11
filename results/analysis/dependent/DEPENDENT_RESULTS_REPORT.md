# Fixed-study dependent-analysis verification

**PASS — every planned dependent fit has exactly one selected complete production attempt; aliases resolve to verified core runs and the infeasible sparse configuration remains unfitted.**

New fits: 29; core aliases: 3; infeasible records: 1.
Maximum PSRF upper 95%: 1.006142; minimum ESS: 1081.7.
Rehashed runtime artifacts: 551; all matched their recorded checksums.

## Non-WT-split HapB3 results by family

| family | scenarios | minimum_OR | maximum_OR | minimum_lower_95_CrI | maximum_upper_95_CrI |
| --- | --- | --- | --- | --- | --- |
| ld | 6 | 1.68237575092161 | 1.74123691921239 | 0.964570154820795 | 2.96405806681801 |
| loo | 23 | 1.50373878453741 | 1.9175959526877 | 0.845846633752714 | 3.28083042831482 |
| primary | 1 | 1.72650158394933 | 1.72650158394933 | 1.00041234292026 | 2.9104162655257 |
| sparse | 1 | 1.70981391317515 | 1.70981391317515 | 1.02722978331838 | 2.82318710417696 |

## WT-split key comparisons

| scenario | comparison | or_median | or_lower_95_cri | or_upper_95_cri |
| --- | --- | --- | --- | --- |
| genotyping_wt_split_loose | HapB3 vs WT_Clean | 1.32130061143995 | 0.62360665454948 | 2.67186393339199 |
| genotyping_wt_split_loose | HapB3 vs WT_Biased | 2.09629825600375 | 1.08464071972078 | 4.08083026815852 |
| genotyping_wt_split_strict | HapB3 vs WT_Clean | 1.16427545232862 | 0.51018290999277 | 2.61600993495622 |
| genotyping_wt_split_strict | HapB3 vs WT_Biased | 2.06184202196731 | 1.11967384706465 | 3.77262482697924 |


These are reference summaries. Public reruns record their own checksums, diagnostics, and software versions under `runs/`. Posterior chains are not included in this repository.
