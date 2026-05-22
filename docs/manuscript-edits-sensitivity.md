# Manuscript Edit Suggestions: Sensitivity Analysis

**Generated from:** sensitivity_summary.csv, sensitivity_forest.png
**Date:** 2026-03-02
**Applies to:** docs/FINAL DRAFT HapB3.pdf

## Instructions for corresponding author

These are suggested additions to address PRISMA 2020 Items 13f (sensitivity
analysis methods) and 20d (sensitivity analysis results). Each edit below
specifies the insertion point and the exact text to add. All numerical values
are drawn from the completed sensitivity analysis.

---

## Edit 1: Methods -- Add "Sensitivity Analyses" subsection

**Location:** Section 2.5 "Network Meta-Analysis", after the paragraph ending
with "...thinning interval of 10" (the MCMC configuration paragraph, page 6).

**Insert new subsection:**

> **2.6 Sensitivity Analyses**
>
> To assess robustness, we conducted 37 sensitivity analyses varying the model specification (homogeneous equicorrelation), link function (logit vs probit), prior specification (inverse-gamma hyperparameters with diffuse and weakly informative settings), and sparse data handling (excluding arms with N ≤ 2 and N ≤ 5). We also performed a complete leave-one-out analysis across all 31 studies to identify influential observations. For each scenario, we compared the HapB3 odds ratio (median, 95% CrI), SUCRA rankings, and model fit statistics against the primary analysis. Reduced MCMC sampling (100,000 iterations, 50,000 burn-in) was used for sensitivity runs, with convergence assessed via the Gelman-Rubin potential scale reduction factor (PSRF ≤ 1.10); scenarios exceeding this threshold were re-run with the primary model's full sampling configuration.

---

## Edit 2: Results -- Add "Sensitivity Analyses" subsection

**Location:** Section 3 "Results", after subsection 3.6 "Study Quality
Assessment" and before Section 4 "Discussion" (page 10).

**Insert new subsection:**

> **3.7 Sensitivity Analyses**
>
> Of 37 sensitivity scenarios, two were infeasible — the heterogeneous correlation model with inverse-gamma prior (unsupported by pcnetmeta) and exclusion of arms with N ≤ 5 (network disconnection, isolating the *13 node and indicating that all *13 evidence derives from very small study arms). Of the 35 feasible scenarios, 34 converged (PSRF ≤ 1.10) and all 34 were consistent with the primary result (Supplementary Table S1, Supplementary Figure S1). Among converged scenarios, HapB3 OR ranged from 1.84 to 2.12 (primary: 2.01), with SUCRA consistently near 74.5% (range: 73.8%–74.8%). In all converged scenarios, HapB3 maintained lowest predicted toxicity among variant genotypes, and no converged scenario produced a 95% CrI that included 1.0. One scenario (exclusion of arms with N ≤ 2) did not converge after escalation to the primary model's full sampling configuration (PSRF 1.17), likely reflecting structural data sparsity in the reduced network rather than insufficient sampling, and should be interpreted with caution.

---

## Edit 3: Discussion -- Add limitations acknowledgment

**Location:** Section 4.3 "Methodological Considerations", in the paragraph
beginning "In terms of limitations regarding study design..." (page 12).
Append to the end of the existing limitations discussion.

**Insert:**

> While sensitivity analyses addressed model specification, link function, prior choice, and data sparsity, we did not perform formal publication bias assessment, as standard methods (e.g., funnel plots, Egger's test) are designed for pairwise meta-analysis of direct comparisons and are not directly applicable to arm-based network meta-analysis, where treatment effects are estimated simultaneously across the network rather than as independent pairwise contrasts. Similarly, a formal certainty-of-evidence assessment using the CINeMA framework was not conducted; future work should address this gap. The inverse-Wishart prior is fixed for the heterogeneous correlation model in pcnetmeta; prior sensitivity testing was therefore limited to the homogeneous equicorrelation model, confounding model structure with prior choice in those comparisons.

---

## Edit 4: Supplementary materials

**Add to supplementary materials:**

- **Table S1:** Sensitivity analysis summary
  - File: output/sensitivity/sensitivity_summary.csv
  - Caption: "Summary of 37 sensitivity analyses testing robustness of HapB3 odds ratio. Each row represents a model perturbation with OR (median, 95% CrI), SUCRA, and convergence diagnostics. Of 35 feasible scenarios, 34 converged (PSRF ≤ 1.10)."

- **Figure S1:** Sensitivity analysis forest plot
  - File: output/sensitivity/sensitivity_forest.png
  - Caption: "Forest plot of HapB3 odds ratio vs wild-type across sensitivity scenarios. Filled circles represent converged main perturbation scenarios; open circles represent leave-one-out iterations; orange triangles represent non-converged scenarios. Horizontal lines represent 95% credible intervals. Dashed vertical line: OR = 1.0 (null). Dotted vertical line: primary analysis HapB3 OR."

---

## Edit 5: Abstract update (optional)

**Location:** Abstract, after the sentence reporting primary results.

**Consider adding:** "Sensitivity analyses (37 scenarios including complete leave-one-out analysis) confirmed robustness of findings; all 34 converged scenarios were consistent with the primary result."

(Check journal word limit before adding.)

---

## Edit 6: Data availability statement

**Location:** End of manuscript, data availability section.

**Append:** "Analytic code, including sensitivity analysis scripts, is
available at ___(insert repository URL)."
