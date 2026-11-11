# Fixed-study final sensitivity tables

Core rows: 12; WT-split comparisons: 4; LD bounds: 6; supported sparse-data rows: 1; grouped source exclusions: 23.
The disconnected N<=5 sparse scenario remains an infeasible record and is not represented as a fitted result.

## WT-split

| scenario | comparison | OR | lower_95_CrI | upper_95_CrI |
| --- | --- | --- | --- | --- |
| genotyping_wt_split_loose | HapB3 vs WT_Clean | 1.32130061143995 | 0.62360665454948 | 2.67186393339199 |
| genotyping_wt_split_loose | HapB3 vs WT_Biased | 2.09629825600375 | 1.08464071972078 | 4.08083026815852 |
| genotyping_wt_split_strict | HapB3 vs WT_Clean | 1.16427545232862 | 0.51018290999277 | 2.61600993495622 |
| genotyping_wt_split_strict | HapB3 vs WT_Biased | 2.06184202196731 | 1.11967384706465 | 3.77262482697924 |

## LD bounds

| scenario | OR | lower_95_CrI | upper_95_CrI |
| --- | --- | --- | --- |
| ld_best_bozina | 1.73589229690812 | 1.0114684657691 | 2.92052397171756 |
| ld_worst_bozina | 1.70306261969217 | 0.977836425927392 | 2.90473680647433 |
| ld_best_wigle_prospective | 1.73802974475205 | 1.01222885237093 | 2.93980569994584 |
| ld_worst_wigle_prospective | 1.68237575092161 | 0.964570154820795 | 2.88410754496204 |
| ld_best_rosmarin_arm_b | 1.74123691921239 | 1.00874394931601 | 2.96405806681801 |
| ld_worst_rosmarin_arm_b | 1.70751832476899 | 0.994587193084426 | 2.89798025127308 |

## Supported sparse-data result

| scenario | OR | lower_95_CrI | upper_95_CrI |
| --- | --- | --- | --- |
| sparse_exclude_n2 | 1.70981391317515 | 1.02722978331838 | 2.82318710417696 |

## Grouped source-influence envelope

Median OR range: 1.5037–1.9176; minimum lower 95% CrI: 0.8458.

All rows are linked to a verified new attempt or an explicitly named completed-core alias.
