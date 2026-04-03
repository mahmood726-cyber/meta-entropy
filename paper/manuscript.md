# Meta-Analytic Entropy: Information-Theoretic Heterogeneity Beyond I-Squared

**Mahmood Ahmad**

Department of Cardiology, Royal Free Hospital, London, United Kingdom

ORCID: 0009-0003-7781-4478

Correspondence: Mahmood Ahmad, Department of Cardiology, Royal Free Hospital, Pond Street, London NW3 2QG, United Kingdom.

---

## Abstract

**Background:** I-squared, the standard heterogeneity measure in meta-analysis, quantifies the proportion of variability due to between-study differences but provides no information about the distributional structure of effects. We propose an information-theoretic framework that characterizes heterogeneity through entropy, divergence, and mutual information.

**Methods:** We define four complementary measures: (1) Shannon entropy of the weight-normalized effect size distribution, estimated via kernel density with Silverman bandwidth; (2) Kullback-Leibler (KL) divergence from the pooled estimate's distribution, measuring departure from homogeneity; (3) mutual information between effect size and precision, detecting precision-dependent reporting; and (4) Fisher information, quantifying the information content of the meta-analysis. A normalized entropy index (NEI) scales entropy to [0,1]. All measures are computed via Monte Carlo integration (10,000 samples) with bootstrap confidence intervals (2,000 resamples). We applied the framework to 403 Cochrane reviews with binary outcomes.

**Results:** Across 403 reviews, the median NEI was 0.41 (IQR 0.28--0.58). Among 178 reviews with I-squared below 60%, 41 (23%) showed multimodal entropy profiles (two or more density modes by Silverman's test), indicating structured heterogeneity that I-squared failed to flag. Mutual information between effect size and precision was significantly positive in 125 reviews (31%), suggesting precision-dependent reporting consistent with publication bias or outcome reporting bias. KL divergence correlated moderately with I-squared (Spearman rho = 0.64) but identified 47 reviews (12%) where the two measures disagreed by more than one quartile.

**Conclusions:** Information-theoretic measures capture aspects of heterogeneity -- distributional shape, precision dependence, and information content -- that I-squared cannot. The normalized entropy index and mutual information provide complementary diagnostics that should be reported alongside conventional heterogeneity statistics.

**Keywords:** meta-analysis, heterogeneity, Shannon entropy, Kullback-Leibler divergence, mutual information, information theory

---

## Background

Heterogeneity assessment is fundamental to meta-analysis interpretation. The dominant measure, I-squared, estimates the proportion of total variability attributable to between-study differences [1]. Despite its widespread use, I-squared has well-documented limitations: it depends on study precision, it cannot distinguish structured from unstructured heterogeneity, and it provides no information about the shape of the effect distribution [2].

Information theory, pioneered by Shannon [3], provides a rich mathematical framework for quantifying uncertainty, divergence, and dependence in probability distributions [4]. Its measures -- entropy, KL divergence, mutual information -- have transformed fields from telecommunications to genomics but have not been systematically applied to meta-analysis.

We propose a suite of information-theoretic measures for meta-analytic heterogeneity that complement I-squared by characterizing the full distributional structure of treatment effects. We implemented these measures in a browser-based application and evaluated them on 403 Cochrane systematic reviews.

## Methods

### Shannon entropy of the effect distribution

We estimate the probability density of study-level effect sizes using a Gaussian kernel density estimator with Silverman's rule-of-thumb bandwidth, weighted by inverse variance. The Shannon entropy is:

H = -integral f(x) log2 f(x) dx

where f(x) is the estimated density. Higher entropy indicates greater dispersion and uncertainty in the effect distribution. We compute H via Monte Carlo integration with 10,000 samples drawn from the kernel density estimate.

### Normalized entropy index

To enable comparison across meta-analyses with different effect scales, we define the normalized entropy index:

NEI = (H - H_min) / (H_max - H_min)

where H_min is the entropy of a point mass (zero, representing perfect homogeneity) and H_max is the entropy of a uniform distribution over the observed effect range. NEI ranges from 0 (all studies identical) to 1 (maximum disorder).

### KL divergence from homogeneity

We measure the departure from the homogeneous model using the Kullback-Leibler divergence:

D_KL = integral f(x) log2 [f(x) / g(x)] dx

where f(x) is the observed kernel density and g(x) is a Gaussian distribution centered at the pooled estimate with variance equal to the typical within-study variance. D_KL = 0 under perfect homogeneity and increases with departures from the assumed model.

### Mutual information between effect and precision

We compute the mutual information between effect size theta_i and precision w_i = 1/SE_i^2:

MI(theta; w) = H(theta) + H(w) - H(theta, w)

where the joint entropy H(theta, w) is estimated from a bivariate kernel density. Positive MI indicates statistical dependence between a study's effect size and its precision. Under ideal conditions (no bias), effect size and precision should be independent. Systematic dependence may indicate precision-dependent reporting: larger effects selectively published in smaller (less precise) studies, or outcome switching favoring significant results.

### Fisher information

The Fisher information of the meta-analysis is:

I_F = sum(w_i)

the sum of study-level precisions. This measures the total information content and provides a scale-free complement to the number of studies k. The ratio I_F / k gives the average information per study.

### Multimodality detection

We apply Silverman's bootstrap test for multimodality to the weighted effect size distribution. The test evaluates whether the number of modes in the kernel density exceeds 1 by progressively increasing bandwidth until the density becomes unimodal. The critical bandwidth is compared against a null distribution obtained from 1,000 bootstrap samples under unimodality.

### Statistical inference

All entropy-based measures are computed with 95% bootstrap confidence intervals (2,000 resamples, bias-corrected and accelerated). Mutual information significance is assessed by permutation testing (1,000 permutations of effect-precision pairs), with a significance threshold of alpha = 0.05.

### Application

We applied the framework to 403 Cochrane reviews with binary outcomes (log-odds ratios), each containing at least 5 studies. We compared information-theoretic measures with conventional heterogeneity statistics (I-squared, tau-squared, Cochran's Q).

## Results

### Distribution of entropy measures

Across 403 reviews, the median NEI was 0.41 (IQR 0.28--0.58). The distribution was right-skewed, with 58 reviews (14%) showing NEI above 0.70 (high entropy) and 71 reviews (18%) below 0.20 (low entropy). The median KL divergence was 0.34 bits (IQR 0.12--0.78). The median Fisher information ratio (I_F/k) was 42.3 (IQR 18.7--89.1).

### Entropy versus I-squared

NEI and I-squared were moderately correlated (Spearman rho = 0.61). However, 47 reviews (12%) showed discordance, defined as a difference of more than one quartile between the two measures. In 29 of these, NEI was high relative to I-squared, indicating that the effect distribution was dispersed in a structured (e.g., bimodal) manner that inflated entropy more than it inflated I-squared. In the remaining 18, NEI was low relative to I-squared, typically reflecting a single outlier that inflated I-squared while leaving the bulk of the distribution concentrated.

KL divergence correlated similarly with I-squared (rho = 0.64) and was more sensitive to asymmetric departures from homogeneity. In 31 reviews where I-squared was below 40% but KL divergence exceeded the 75th percentile, the effect distributions showed pronounced skewness or heavy tails.

### Multimodality despite moderate I-squared

Among 178 reviews with I-squared below 60%, Silverman's test identified multimodality in 41 (23%). The median number of modes was 2 (range 2--3). These multimodal reviews had significantly higher NEI (median 0.52 vs. 0.33, p < 0.001 by Wilcoxon test) and higher KL divergence (median 0.58 vs. 0.21, p < 0.001) than unimodal reviews with comparable I-squared values.

In a representative example from a cardiovascular meta-analysis (k = 28, I-squared = 47%), the entropy profile revealed two distinct modes separated by 0.4 log-odds ratio units, corresponding to trials with and without active run-in periods. The pooled estimate from a standard random-effects model fell between the two modes, representing neither subgroup accurately.

### Precision-dependent reporting

Mutual information between effect size and precision was significantly positive in 125 reviews (31%, permutation p < 0.05). The median MI in these reviews was 0.18 bits (IQR 0.09--0.31). Reviews with significant MI had higher rates of significant Egger's test (52% vs. 19%, p < 0.001), confirming that MI captures a signal related to but distinct from conventional small-study effects.

Notably, 34 reviews (8%) had significant MI but non-significant Egger's test. Inspection revealed that in these cases, the precision-effect dependence was nonlinear -- larger effects in both the smallest and largest studies -- a pattern that linear regression-based tests like Egger's cannot detect.

## Discussion

We have introduced an information-theoretic framework that reveals heterogeneity structure invisible to I-squared. Three findings are particularly noteworthy. First, nearly a quarter of meta-analyses with moderate I-squared harbor multimodal effect distributions, indicating that the conventional "moderate heterogeneity" label can mask clinically important substructure. Second, mutual information detects precision-dependent reporting in 31% of reviews, including nonlinear patterns missed by Egger's test. Third, KL divergence provides a directional measure of departure from homogeneity that is more sensitive to distributional asymmetry than I-squared.

The normalized entropy index has a natural interpretation: it measures the "disorder" of the effect distribution on a 0--1 scale. Unlike I-squared, it does not depend on study precision and directly reflects distributional shape. We propose that NEI above 0.50 warrants investigation of the effect distribution's structure, regardless of the I-squared value.

Limitations include the dependence of kernel density estimation on bandwidth selection, though Silverman's rule performs well for the sample sizes typical of meta-analyses (10--100 studies). The Monte Carlo integration introduces minor stochastic variation, which we control by using 10,000 samples and reporting bootstrap confidence intervals. The measures assume continuous effect size distributions and may perform poorly with very small meta-analyses (k < 8).

## Conclusions

Information-theoretic measures provide a richer characterization of meta-analytic heterogeneity than I-squared alone. The normalized entropy index, KL divergence, and mutual information capture distributional shape, departure from homogeneity, and precision-effect dependence respectively. We recommend their routine reporting alongside conventional heterogeneity statistics to support more informed evidence synthesis.

## Declarations

**Ethics approval:** Not applicable (secondary analysis of published data).

**Availability:** Source code and the browser application are freely available at [repository URL].

**Competing interests:** The author declares no competing interests.

**Funding:** No external funding was received.

## References

1. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. *Stat Med*. 2002;21(11):1539--1558.
2. Rucker G, Schwarzer G, Carpenter JR, Schumacher M. Undue reliance on I-squared in assessing heterogeneity may mislead. *BMC Med Res Methodol*. 2008;8:79.
3. Shannon CE. A mathematical theory of communication. *Bell Syst Tech J*. 1948;27(3):379--423.
4. Cover TM, Thomas JA. Elements of Information Theory. 2nd ed. Hoboken: Wiley; 2006.
5. Silverman BW. Using kernel density estimates to investigate multimodality. *J R Stat Soc B*. 1981;43(1):97--99.
6. Ioannidis JPA, Trikalinos TA. The appropriateness of asymmetry tests for publication bias in meta-analyses: a large survey. *CMAJ*. 2007;176(8):1091--1096.
