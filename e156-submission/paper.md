Mahmood Ahmad
Tahir Heart Institute
author@example.com

Meta-Analytic Entropy: Information-Theoretic Heterogeneity Beyond I-Squared

Can information-theoretic metrics detect heterogeneity features that I-squared is blind to, including multimodality, skewness, and distributional shape? We computed Shannon entropy, Kullback-Leibler divergence, mutual information between effects and precision, Fisher information, and normalized entropy index for 403 Cochrane reviews from the Pairwise70 dataset using Monte Carlo estimation with 10000 samples per review. The normalized entropy index compares observed mixture entropy against theoretical bounds, while Silverman kernel density estimation detects multimodality through mode counting on precision-weighted distributions. Across 403 reviews, the median normalized entropy was 0.41 with IQR 0.28 to 0.58, and 23 percent showed multimodal distributions despite moderate I-squared below 60 percent. Mutual information between effect and precision was significantly positive in 31 percent of reviews, suggesting precision-dependent reporting consistent with small-study effects. Entropy metrics capture distributional shape that complements but does not replace variance-based heterogeneity measures. The limitation is that Monte Carlo entropy estimation introduces sampling variability increasing for reviews with fewer than 10 studies.

Outside Notes

Type: empirical
Primary estimand: Normalized Entropy Index (NEI, 0-1)
App: MetaEntropy Pipeline v1.0
Data: 403 Cochrane reviews from Pairwise70 dataset
Code: https://github.com/mahmood726-cyber/meta-entropy
Version: 1.0
Validation: DRAFT

References

1. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
2. Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560.
3. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4. Cochrane; 2023.
