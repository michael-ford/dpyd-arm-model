# Fixed-study model and reporting scope

The analysis fits arm-level severe-toxicity counts with a binomial likelihood and logit link. Each analytical stratum has its own fixed intercept, so genotype effects are estimated from concurrent within-stratum comparisons without a shared population distribution for baseline study risks. White et al. (2019) and Piepho et al. (2024) support this fixed-study principle. Piepho's exact arm-based/contrast-based equivalence result concerns linear models; exact equivalence is not claimed for this binomial-logit model.

WT is the global reference. An independent normal deviation is placed on every arm, with one common arm-level heterogeneity parameter tau. Piepho et al. provide the direct network-meta-analysis support for symmetric arm deviations. Jackson et al. (2018) provide supporting pairwise-meta-analysis evidence against assigning residual heterogeneity only to a designated non-reference arm.

Non-WT log-odds effects have independent N(0, 2²) priors. Each study intercept has an independent diffuse N(0, 100²) prior. Arm-level tau has a half-normal prior with scale 1; contrast-level heterogeneity is sqrt(2) × tau. Friede et al. (2017) support the half-normal prior family, not this exact arm-scale choice.

Eligible zero-event and all-event arms enter the binomial likelihood without continuity corrections (Stijnen et al., 2010). Outcome-sampled studies are excluded from the primary because their arm proportions do not estimate cohort incidence. Connected multi-arm strata without WT contribute relative genotype information; isolated single-arm strata cannot identify a genotype effect separately from their fixed intercept. Pure-variant primary arms are mutually exclusive; the committed sensitivities bound compound-carrier allocation, overlap, endpoint windows, clinical management, comparator classification, sparse data, and grouped source influence.

Four chains use 50,000 adaptation, 125,000 burn-in, and 250,000 sampling iterations per chain, thinning 10, and master seed 12345. The committed schedules specify deterministic chain seeds. Production settings permit no fallback or extension. All non-reference effects, pairwise log contrasts, and heterogeneity parameters require PSRF upper 95% ≤1.05 and ESS ≥1,000.

Only relative odds ratios and 95% credible intervals are reported. The fitted model does not independently anchor a transportable clinical baseline risk (Dias and Ades, 2016). Rankings and a treatment-correlation matrix are not estimated. The comparator-classification QBA is a prevalence-based screening calculation, not an input correction or a clinical absolute-risk estimate.

The fixed design-by-treatment extension follows the global inconsistency framework of Higgins et al. (2012). Both prior-scale fits converged, but failed the information criteria. Inconsistency could not be reliably assessed; this is not evidence of consistency. The numerical gates are retained in `analysis/inconsistency/INCONSISTENCY_REPORT.md`.

## References

- White IR, Turner RM, Karahalios A, Salanti G. A comparison of arm-based and contrast-based models for network meta-analysis. Statistics in Medicine. 2019;38:5197–5213. https://doi.org/10.1002/sim.8360
- Piepho HP, Madden LV, Williams ER. The use of fixed study main effects in arm-based network meta-analysis. Research Synthesis Methods. 2024;15:747–750. https://doi.org/10.1002/jrsm.1721
- Jackson D, Law M, Stijnen T, Viechtbauer W, White IR. A comparison of 7 random-effects models for meta-analyses that estimate the summary odds ratio. Statistics in Medicine. 2018;37:1059–1085. https://doi.org/10.1002/sim.7588
- Stijnen T, Hamza TH, Özdemir P. Random effects meta-analysis of event outcome in the framework of the generalized linear mixed model with applications in sparse data. Statistics in Medicine. 2010;29:3046–3067. https://doi.org/10.1002/sim.4040
- Friede T, Röver C, Wandel S, Neuenschwander B. Meta-analysis of few small studies in orphan diseases. Research Synthesis Methods. 2017;8:79–91. https://doi.org/10.1002/jrsm.1217
- Dias S, Ades AE. Absolute or relative effects? Arm-based synthesis of trial data. Research Synthesis Methods. 2016;7:23–28. https://doi.org/10.1002/jrsm.1184
- Higgins JPT, Jackson D, Barrett JK, Lu G, Ades AE, White IR. Consistency and inconsistency in network meta-analysis: concepts and models for multi-arm studies. Research Synthesis Methods. 2012;3:98–110. https://doi.org/10.1002/jrsm.1044
- R Core Team. R: A Language and Environment for Statistical Computing. 2024. https://www.R-project.org/
- Plummer M. JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling. Proceedings of the 3rd International Workshop on Distributed Statistical Computing. 2003. https://www.r-project.org/conferences/DSC-2003/Proceedings/Plummer.pdf
- Plummer M. rjags: Bayesian Graphical Models using MCMC. R package version 4-17. 2025. https://doi.org/10.32614/CRAN.package.rjags
- Plummer M, Best N, Cowles K, Vines K. CODA: convergence diagnosis and output analysis for MCMC. R News. 2006;6:7–11. https://journal.r-project.org/articles/RN-2006-002/
