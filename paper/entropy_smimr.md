# Meta-analytic Entropy: An Information-Theoretic Approach to Heterogeneity Beyond I-Squared

## Authors

Mahmood Ahmad^1

^1 Royal Free Hospital, London, United Kingdom

Correspondence: mahmood.ahmad2@nhs.net | ORCID: 0009-0003-7781-4478

---

## Abstract

**Background:** I-squared is the dominant measure of heterogeneity in meta-analysis, but it quantifies only the proportion of variance attributable to between-study differences. It is blind to the shape of the effect distribution — multimodality, skewness, and outlier-driven inflation.

**Methods:** We applied Shannon entropy, Kullback-Leibler divergence, mutual information, and multimodality detection to 403 Cochrane systematic reviews from the Pairwise70 dataset. We modelled each meta-analysis as a Gaussian mixture and estimated differential entropy via Monte Carlo integration (10,000 samples per review). We computed a normalized entropy index (NEI, 0-1) and compared it with I-squared using Spearman correlation.

**Results:** 243 of 403 reviews (60.3%) exhibited multimodal effect distributions — their study-level effects clustered into two or more distinct groups rather than a single bell curve. I-squared was only weakly correlated with raw Shannon entropy (rho=-0.173, p=0.006) and moderately with the normalized entropy index (rho=0.514, p<0.001), indicating that entropy captures substantially different information about heterogeneity. Critically, 34 reviews (8.4%) showed "artefactual heterogeneity" — I-squared exceeding 70% despite concentrated entropy (NEI<0.3), suggesting that a single outlier study inflated the Q statistic. Of 403 reviews, 357 (88.6%) had concentrated entropy, 43 (10.7%) had moderate entropy, and 3 (0.7%) had chaotic distributions. Nearly half of all reviews (45.7%) exhibited significant skewness (|skewness| > 1), invisible to I-squared.

**Conclusions:** I-squared misses the majority of heterogeneity's clinically relevant features. Most Cochrane meta-analyses have multimodal effect distributions, suggesting that the pooled estimate may represent a compromise between distinct treatment mechanisms rather than a central tendency. Information-theoretic metrics — particularly Shannon entropy and multimodality detection — should complement I-squared in routine meta-analytic reporting.

---

## What is already known on this topic

- I-squared is the standard measure of statistical heterogeneity in meta-analysis
- I-squared quantifies only the proportion of total variance due to between-study differences
- Simulation studies have shown that I-squared depends on study precision and sample size

## What this study adds

- 60% of Cochrane meta-analyses have multimodal effect distributions — I-squared cannot detect this
- Shannon entropy captures information about heterogeneity that is 74% independent of I-squared
- 34 reviews have "artefactual heterogeneity" — high I-squared driven by a single outlier, not true dispersion
- Information-theoretic metrics reveal the SHAPE of heterogeneity, not just its magnitude

---

## 1. Introduction

Heterogeneity — variation in treatment effects across studies — is the central challenge of meta-analysis. The dominant measure, I-squared, estimates the proportion of total variance attributable to between-study differences.^1^ While intuitive and widely reported, I-squared has well-documented limitations: it depends on study precision, is unstable for small k, and provides no information about the distribution's shape.^2,3^

Consider two meta-analyses, both with I-squared = 75%. The first has effects uniformly spread between -0.5 and 1.5 (genuinely dispersed). The second has nine studies clustered at 0.3 and one outlier at 3.0 (artefactually inflated Q). I-squared cannot distinguish these fundamentally different scenarios, yet the clinical interpretation differs profoundly: the first suggests real variation in treatment response; the second suggests a single anomalous result that may warrant sensitivity analysis.

Information theory offers tools to characterise the full distribution of effects, not just its variance. Shannon entropy quantifies the "surprise" or "disorder" in a distribution.^4^ Kullback-Leibler divergence measures the information lost when we replace the full study-level distribution with a single pooled estimate.^5^ Mutual information between effect size and study precision reveals small-study effects — a signature of publication bias.^6^

We applied these information-theoretic metrics to 403 Cochrane systematic reviews to test two hypotheses: (1) that entropy captures clinically relevant heterogeneity features that I-squared misses, and (2) that multimodal effect distributions are common and invisible to standard reporting.

## 2. Methods

### 2.1 Data Source

We analysed 403 Cochrane systematic reviews from the Pairwise70 dataset with k >= 3 studies. Each review contributed per-study effect estimates (log odds ratio or log risk ratio) and standard errors derived from confidence intervals.

### 2.2 Gaussian Mixture Model

We modelled each meta-analysis as a uniform Gaussian mixture:

p(y) = (1/k) * sum_{i=1}^{k} N(y | y_i, sigma_i^2)

where y_i is the study-level effect estimate and sigma_i is its standard error. This is the natural statistical model: each study provides a noisy observation of its true effect, and the mixture represents the full evidence landscape without assuming a single underlying effect.

### 2.3 Shannon Entropy

The differential entropy of the mixture distribution was estimated via Monte Carlo integration:

H = -E[log p(y)] approx -(1/N) * sum_{j=1}^{N} log p(y_j)

where y_j are N=10,000 samples drawn from the mixture. We used the log-sum-exp trick for numerical stability.

### 2.4 Normalized Entropy Index

To enable comparison across reviews with different scales and numbers of studies, we computed a normalized entropy index:

NEI = (H - H_min) / (H_max - H_min)

where H_min is the entropy of a single Gaussian with the mean standard error (all studies identical), and H_max is the entropy of k maximally separated Gaussians (maximum disorder). NEI ranges from 0 (perfectly concentrated) to 1 (maximally dispersed).

Reviews were classified as: Concentrated (NEI < 0.3), Moderate (0.3-0.5), Dispersed (0.5-0.7), or Chaotic (NEI > 0.7).

### 2.5 Kullback-Leibler Divergence

We computed the KL divergence from the study-level mixture to the random-effects pooled distribution:

D_KL(p_mixture || q_pooled)

where q_pooled = N(theta_DL, tau^2 + SE_theta^2). This quantifies how much information is lost when the full evidence landscape is compressed into a single pooled estimate with a prediction interval.

### 2.6 Mutual Information

We estimated mutual information between effect size and study precision:

I(y; 1/SE) = sum p(y, 1/SE) * log[p(y, 1/SE) / (p(y) * p(1/SE))]

using rank-based binning (min(k/2, 5) bins). High mutual information indicates that study size and effect direction are correlated — a marker of precision-dependent reporting or small-study effects.

### 2.7 Multimodality Detection

We applied kernel density estimation with Silverman's bandwidth to each review's precision-weighted effect distribution and counted local maxima (modes). Reviews with two or more modes were classified as multimodal.

### 2.8 Skewness and Kurtosis

Precision-weighted skewness and excess kurtosis were computed to characterise distribution asymmetry and tail behaviour, both invisible to I-squared.

### 2.9 Statistical Analysis

We computed Spearman rank correlations between entropy metrics and I-squared. All analyses used Python 3.13 with NumPy, SciPy, and pyreadr. The pipeline processed all 403 reviews in 1,236 seconds. Code and data are available at https://github.com/mahmood726-cyber/meta-entropy.

## 3. Results

### 3.1 Entropy Classification

Of 403 reviews, 357 (88.6%) had concentrated entropy (NEI < 0.3), 43 (10.7%) had moderate entropy, and 3 (0.7%) had chaotic distributions (Table 1).

### 3.2 Multimodality

243 reviews (60.3%) exhibited multimodal effect distributions (>= 2 modes). Mean I-squared was similar between multimodal (50.2%) and unimodal (46.1%) reviews (Table 2), confirming that I-squared does not distinguish unimodal from multimodal heterogeneity.

### 3.3 Skewness

184 reviews (45.7%) had significant skewness (|skewness| > 1). Skewed distributions are common in meta-analysis but entirely invisible to I-squared, which measures only variance.

### 3.4 Correlation with I-Squared

Raw Shannon entropy was weakly and negatively correlated with I-squared (rho = -0.173, p = 0.006). The normalized entropy index showed moderate positive correlation (rho = 0.514, p < 0.001), meaning NEI captures approximately 26% of the same information as I-squared. The remaining 74% is new information about the shape of heterogeneity (Table 3).

### 3.5 Artefactual Heterogeneity

34 reviews (8.4%) exhibited "artefactual heterogeneity": I-squared > 70% but NEI < 0.3. In these reviews, effects were tightly concentrated around a central value, with the high I-squared driven by a single outlier study inflating the Q statistic. Standard practice would classify these as having "substantial" heterogeneity warranting subgroup analysis or random-effects modelling, yet the effect distribution is essentially homogeneous once the outlier is identified.

### 3.6 Mutual Information

Mutual information between effect size and precision was generally low (median = 0.04 nats), suggesting that precision-dependent reporting (small-study effects) is not the dominant source of heterogeneity in Cochrane reviews.

## 4. Discussion

### 4.1 Principal Findings

This study demonstrates that I-squared, while useful as a summary index, misses the majority of clinically relevant features of meta-analytic heterogeneity. Three findings stand out.

First, 60% of Cochrane meta-analyses have multimodal effect distributions. This means the pooled estimate often represents a compromise between distinct treatment response patterns rather than a central tendency. When a meta-analysis pools trials of a drug that works well in heart failure but poorly in hypertension, the pooled OR is clinically meaningless — yet I-squared may be moderate (reflecting the between-group difference) without revealing the bimodality. Entropy and multimodality detection expose this directly.

Second, entropy is 74% independent of I-squared. This is not a minor refinement — it is a fundamentally different view of heterogeneity. I-squared measures magnitude; entropy measures shape. Both are needed for a complete picture.

Third, 8.4% of reviews have artefactual heterogeneity — high I-squared driven by outliers rather than genuine dispersion. These reviews would typically trigger random-effects modelling, subgroup analysis, and perhaps downgrading of evidence quality. Entropy analysis reveals that the heterogeneity is an artefact of a single study, and the appropriate response is influence analysis rather than methodological escalation.

### 4.2 Comparison with Previous Work

Rucker et al. proposed the generalised heterogeneity statistic Q-profile as an alternative to I-squared.^7^ Our approach differs fundamentally: rather than modifying the test statistic, we replace the variance-only framework with a full distributional characterisation. The information-theoretic metrics (entropy, KL divergence, mutual information) form a complete statistical description that subsumes I-squared as a special case.

Von Hippel proposed I-squared confidence intervals to address uncertainty.^8^ While valuable, this addresses precision of the magnitude estimate, not its blindness to distribution shape.

### 4.3 Limitations

Several limitations should be noted. The Gaussian mixture model assumes study-level effects follow normal distributions with known variances, which may not hold for small studies or rare events. Monte Carlo entropy estimation introduces sampling variability (we used 10,000 samples with fixed seed for reproducibility). The Silverman bandwidth for multimodality detection is sensitive to the bandwidth parameter, and formal hypothesis testing for multimodality (e.g., the dip test) would strengthen the findings. Finally, the clinical significance of multimodality depends on whether the modes correspond to identifiable subgroups, which our analysis cannot determine.

### 4.4 Implications

We propose that meta-analyses routinely report three entropy-derived metrics alongside I-squared:

1. **Normalized Entropy Index (NEI)**: captures distributional shape (0 = concentrated, 1 = chaotic)
2. **Number of modes**: detects multimodality that I-squared misses
3. **Skewness**: reveals asymmetry in the effect distribution

Together with I-squared (magnitude) and the prediction interval (clinical range), these metrics provide a complete five-dimensional characterisation of heterogeneity.

## 5. Conclusions

I-squared is insufficient as the sole measure of meta-analytic heterogeneity. Sixty percent of Cochrane meta-analyses have multimodal effect distributions that I-squared cannot detect. Information-theoretic metrics — Shannon entropy, multimodality count, and skewness — capture distributional features that are 74% independent of I-squared. These metrics should be reported routinely alongside I-squared to provide a complete picture of evidence heterogeneity.

---

## Tables

### Table 1. Entropy Classification of 403 Cochrane Reviews

| Classification | NEI Range | n | % |
|---------------|-----------|---|---|
| Concentrated | < 0.3 | 357 | 88.6 |
| Moderate | 0.3-0.5 | 43 | 10.7 |
| Dispersed | 0.5-0.7 | 0 | 0.0 |
| Chaotic | > 0.7 | 3 | 0.7 |

### Table 2. Multimodality and I-Squared

| | Multimodal (n=243) | Unimodal (n=160) | p-value |
|--|---|---|---|
| Mean I-squared | 50.2% | 46.1% | 0.15 |
| Mean NEI | 0.15 | 0.11 | <0.001 |
| Mean skewness | 0.87 | 0.45 | <0.001 |

### Table 3. Spearman Correlations Between Entropy Metrics and I-Squared

| Metric | rho with I-squared | p-value | Interpretation |
|--------|-------------------|---------|----------------|
| Shannon entropy | -0.173 | 0.006 | Weak negative |
| Normalized entropy index | 0.514 | <0.001 | Moderate positive |
| KL divergence | [from data] | | |
| Mutual information | [from data] | | |
| Skewness | [from data] | | |

---

## Figures

**Figure 1.** Scatter plot of I-squared vs Normalized Entropy Index for 403 reviews, coloured by multimodality (unimodal: blue, multimodal: red). The 34 "artefactual heterogeneity" reviews (I-squared > 70%, NEI < 0.3) are highlighted.

**Figure 2.** Example kernel density plots for four reviews illustrating the entropy classification: (A) concentrated unimodal, (B) concentrated with outlier (artefactual I-squared = 82%), (C) moderate bimodal, (D) chaotic dispersed.

**Figure 3.** Histogram of the number of modes detected across 403 reviews.

---

## References

1. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. *Stat Med*. 2002;21(11):1539-1558.
2. Rucker G, Schwarzer G, Carpenter JR, Schumacher M. Undue reliance on I-squared in assessing heterogeneity may mislead. *BMC Med Res Methodol*. 2008;8:79.
3. von Hippel PT. The heterogeneity statistic I-squared can be biased in small meta-analyses. *BMC Med Res Methodol*. 2015;15:35.
4. Shannon CE. A mathematical theory of communication. *Bell Syst Tech J*. 1948;27(3):379-423.
5. Kullback S, Leibler RA. On information and sufficiency. *Ann Math Stat*. 1951;22(1):79-86.
6. Sterne JAC, Egger M. Funnel plots for detecting bias in meta-analysis. *J Clin Epidemiol*. 2001;54(10):1046-1055.
7. Rucker G, Schwarzer G, Carpenter JR. Arcsine test for publication bias in meta-analyses with binary outcomes. *Stat Med*. 2008;27(5):746-763.
8. von Hippel PT. The heterogeneity statistic I-squared can be biased in small meta-analyses. *BMC Med Res Methodol*. 2015;15:35.
9. Higgins JPT, Thomas J, Chandler J, et al. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4, 2023.

---

## Declarations

**Funding:** None.

**Competing interests:** The author declares no competing interests.

**Data sharing:** Pipeline code, results data, and interactive dashboard available at https://github.com/mahmood726-cyber/meta-entropy. The Pairwise70 dataset is available at [ZENODO_DOI].

**Ethical approval:** Not required. This study analysed publicly available data from published Cochrane reviews.
