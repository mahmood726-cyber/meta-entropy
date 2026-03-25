"""Meta-analytic Entropy: Information-theoretic heterogeneity for 403 Cochrane reviews.

For each review, computes:
  1. Shannon entropy of the effect size distribution (continuous: differential entropy)
  2. Kullback-Leibler divergence from pooled estimate (information loss from pooling)
  3. Mutual information between study effects and study precision
  4. Fisher information (total precision of the evidence)
  5. Normalized entropy index (0-1, comparable across reviews)

Key thesis: I-squared captures the proportion of variance due to heterogeneity,
but entropy captures the SHAPE of the heterogeneity — multimodality, skewness,
heavy tails — things I-squared is blind to.

Usage: python pipeline.py
"""

import csv
import json
import math
import time
import numpy as np
import pyreadr
from pathlib import Path
from scipy import stats as sp_stats

PAIRWISE_DIR = Path(r'C:\Models\Pairwise70\data')
OUTPUT_DIR = Path(r'C:\Models\MetaEntropy\data\output')


# ═══════════════════════════════════════
# INFORMATION-THEORETIC METRICS
# ═══════════════════════════════════════

def shannon_entropy(yi, sei):
    """Differential (continuous) Shannon entropy of the effect distribution.

    Models the meta-analytic distribution as a mixture of Gaussians:
      p(y) = sum_i (1/k) * N(y | yi, sei^2)

    Entropy estimated via Monte Carlo sampling from this mixture.
    """
    k = len(yi)
    if k < 3:
        return None

    # Sample from the mixture
    n_samples = 10000
    rng = np.random.RandomState(42)

    # Choose which component to sample from (uniform mixture)
    components = rng.randint(0, k, n_samples)
    samples = np.array([rng.normal(yi[c], sei[c]) for c in components])

    # Evaluate mixture density at each sample
    log_densities = np.zeros(n_samples)
    for i in range(n_samples):
        # log p(x) = log( (1/k) * sum_j N(x | yj, sej^2) )
        log_components = np.array([
            -0.5 * ((samples[i] - yi[j]) / sei[j])**2 - np.log(sei[j]) - 0.5 * np.log(2 * np.pi)
            for j in range(k)
        ])
        log_densities[i] = np.log(1.0 / k) + log_sum_exp(log_components)

    # H = -E[log p(x)]
    entropy = -np.mean(log_densities)
    return float(entropy)


def log_sum_exp(x):
    """Numerically stable log-sum-exp."""
    x = np.asarray(x)
    mx = np.max(x)
    return mx + np.log(np.sum(np.exp(x - mx)))


def kl_divergence_from_pooled(yi, sei, theta, se_theta, tau2):
    """KL divergence: D(mixture || N(theta, tau2+se^2)).

    Measures information lost when we replace the full study-level
    distribution with the single pooled estimate.
    """
    k = len(yi)
    if k < 3:
        return None

    # Sample from the mixture
    n_samples = 10000
    rng = np.random.RandomState(42)
    components = rng.randint(0, k, n_samples)
    samples = np.array([rng.normal(yi[c], sei[c]) for c in components])

    # Pooled distribution: N(theta, tau2 + se_theta^2)
    sigma_pooled = math.sqrt(tau2 + se_theta**2) if (tau2 + se_theta**2) > 0 else 0.01

    kl = 0.0
    for i in range(n_samples):
        # log p_mixture(x)
        log_mix = np.log(1.0 / k) + log_sum_exp([
            -0.5 * ((samples[i] - yi[j]) / sei[j])**2 - np.log(sei[j]) - 0.5 * np.log(2 * np.pi)
            for j in range(k)
        ])

        # log q_pooled(x)
        log_pooled = -0.5 * ((samples[i] - theta) / sigma_pooled)**2 - math.log(sigma_pooled) - 0.5 * math.log(2 * math.pi)

        kl += (log_mix - log_pooled)

    kl /= n_samples
    return max(0.0, float(kl))  # KL >= 0 by definition


def mutual_information_effect_precision(yi, sei):
    """Mutual information I(effect; precision).

    Measures whether larger/smaller studies systematically report
    different effects. High MI suggests precision-dependent reporting
    (small study effects / publication bias signature).
    """
    k = len(yi)
    if k < 5:
        return None

    # Discretize into bins
    n_bins = min(k // 2, 5)
    if n_bins < 2:
        return None

    # Rank-based binning to handle outliers
    yi_ranks = np.argsort(np.argsort(yi))
    sei_ranks = np.argsort(np.argsort(sei))

    yi_bins = np.minimum(yi_ranks * n_bins // k, n_bins - 1)
    sei_bins = np.minimum(sei_ranks * n_bins // k, n_bins - 1)

    # Joint and marginal counts
    joint = np.zeros((n_bins, n_bins))
    for i in range(k):
        joint[yi_bins[i], sei_bins[i]] += 1

    joint /= k  # normalize to probabilities
    py = joint.sum(axis=1)
    ps = joint.sum(axis=0)

    # MI = sum p(y,s) * log(p(y,s) / (p(y)*p(s)))
    mi = 0.0
    for i in range(n_bins):
        for j in range(n_bins):
            if joint[i, j] > 0 and py[i] > 0 and ps[j] > 0:
                mi += joint[i, j] * math.log(joint[i, j] / (py[i] * ps[j]))

    return float(mi)


def fisher_information(sei):
    """Total Fisher information: sum(1/sei^2).

    Measures the total precision of the evidence base.
    """
    return float(np.sum(1.0 / sei**2))


def normalized_entropy_index(entropy, k, sei):
    """Normalize entropy to [0, 1] by comparing to maximum entropy
    distribution (uniform over the observed range) and minimum entropy
    (all studies identical).

    NEI = (H - H_min) / (H_max - H_min)
    """
    if entropy is None or k < 3:
        return None

    # H_min: all studies at the same value → entropy of a single Gaussian with mean SE
    mean_se = float(np.mean(sei))
    h_min = 0.5 * math.log(2 * math.pi * math.e * mean_se**2)

    # H_max: uniform distribution over observed range
    # For uniform on [a, b]: H = log(b - a)
    # But we're comparing mixture entropy, so use the entropy of k widely-spaced Gaussians
    # Approximate: H_max ≈ log(k) + 0.5*log(2*pi*e) + log(mean_se)
    h_max = math.log(max(k, 2)) + 0.5 * math.log(2 * math.pi * math.e * mean_se**2)

    if h_max <= h_min:
        return 0.5

    nei = (entropy - h_min) / (h_max - h_min)
    return float(max(0.0, min(1.0, nei)))


def detect_multimodality(yi, sei):
    """Silverman's test for multimodality.

    Uses kernel density estimation with increasing bandwidth to find
    the number of modes. Returns estimated number of modes.
    """
    k = len(yi)
    if k < 5:
        return 1

    # Weighted KDE (weight by precision)
    weights = 1.0 / sei**2
    weights /= weights.sum()

    # Silverman's rule of thumb bandwidth
    std = np.sqrt(np.average((yi - np.average(yi, weights=weights))**2, weights=weights))
    h = 0.9 * std * k**(-0.2)

    if h < 1e-10:
        return 1

    # Evaluate density on a grid
    grid = np.linspace(np.min(yi) - 3*h, np.max(yi) + 3*h, 200)
    density = np.zeros(len(grid))
    for i in range(k):
        density += weights[i] * np.exp(-0.5 * ((grid - yi[i]) / h)**2) / (h * np.sqrt(2 * np.pi))

    # Count modes (local maxima)
    modes = 0
    for i in range(1, len(density) - 1):
        if density[i] > density[i-1] and density[i] > density[i+1]:
            modes += 1

    return max(1, modes)


def skewness_kurtosis(yi, sei):
    """Precision-weighted skewness and excess kurtosis of effect distribution."""
    k = len(yi)
    if k < 4:
        return 0, 0

    wi = 1.0 / sei**2
    wi_norm = wi / wi.sum()

    mu = float(np.sum(wi_norm * yi))
    m2 = float(np.sum(wi_norm * (yi - mu)**2))
    m3 = float(np.sum(wi_norm * (yi - mu)**3))
    m4 = float(np.sum(wi_norm * (yi - mu)**4))

    if m2 < 1e-15:
        return 0, 0

    skew = m3 / m2**1.5
    kurt = m4 / m2**2 - 3  # excess kurtosis

    return float(skew), float(kurt)


# ═══════════════════════════════════════
# DATA LOADING (reuses Fragility Atlas loader pattern)
# ═══════════════════════════════════════

def load_review(rda_path):
    """Load one RDA, select primary analysis, return yi/sei/metadata."""
    result = pyreadr.read_r(str(rda_path))
    df = list(result.values())[0].copy()
    df.columns = df.columns.str.replace(' ', '.', regex=False)

    review_id = rda_path.stem.split('_')[0]

    import pandas as pd
    groups = []
    for (grp, num), sub in df.groupby(['Analysis.group', 'Analysis.number']):
        has_binary = (sub['Experimental.cases'].notna() & (sub['Experimental.cases'] > 0)).any()
        groups.append({'grp': grp, 'num': num, 'k': len(sub), 'binary': has_binary})

    if not groups:
        return None

    gdf = pd.DataFrame(groups)
    binary = gdf[gdf['binary']]
    best = binary.loc[binary['k'].idxmax()] if len(binary) > 0 else gdf.loc[gdf['k'].idxmax()]
    primary = df[(df['Analysis.group'] == best['grp']) & (df['Analysis.number'] == best['num'])]

    has_binary = (primary['Experimental.cases'].notna() & (primary['Experimental.cases'] > 0)).any()
    if has_binary:
        scale = 'ratio'
    else:
        means = primary['Mean'].dropna()
        scale = 'ratio' if len(means) > 0 and (means > 0).all() else 'difference'

    if scale == 'ratio':
        valid = (primary['Mean'].notna() & (primary['Mean'] > 0) &
                 primary['CI.start'].notna() & (primary['CI.start'] > 0) &
                 primary['CI.end'].notna() & (primary['CI.end'] > 0))
        sub = primary[valid]
        if len(sub) < 3:
            return None
        yi = np.log(sub['Mean'].values.astype(float))
        sei = (np.log(sub['CI.end'].values.astype(float)) - np.log(sub['CI.start'].values.astype(float))) / (2 * 1.96)
    else:
        valid = primary['Mean'].notna() & primary['CI.start'].notna() & primary['CI.end'].notna()
        sub = primary[valid]
        if len(sub) < 3:
            return None
        yi = sub['Mean'].values.astype(float)
        sei = (sub['CI.end'].values.astype(float) - sub['CI.start'].values.astype(float)) / (2 * 1.96)

    ok = (sei > 0) & np.isfinite(yi) & np.isfinite(sei)
    yi, sei = yi[ok], sei[ok]

    if len(yi) < 3:
        return None

    # DL pooling
    wi = 1.0 / sei**2
    sw = np.sum(wi)
    theta_fe = np.sum(wi * yi) / sw
    Q = float(np.sum(wi * (yi - theta_fe)**2))
    C = float(sw - np.sum(wi**2) / sw)
    tau2 = max(0, (Q - (len(yi) - 1)) / C) if C > 0 else 0
    ws = 1.0 / (sei**2 + tau2)
    sws = np.sum(ws)
    theta = float(np.sum(ws * yi) / sws)
    se_theta = float(1.0 / math.sqrt(sws))
    I2 = max(0, (Q - (len(yi) - 1)) / Q * 100) if Q > 0 else 0

    return {
        'review_id': review_id,
        'yi': yi, 'sei': sei,
        'k': len(yi), 'scale': scale,
        'theta': theta, 'se_theta': se_theta,
        'tau2': tau2, 'I2': I2, 'Q': Q,
    }


# ═══════════════════════════════════════
# MAIN
# ═══════════════════════════════════════

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Meta-analytic Entropy Pipeline")
    print("=" * 40)

    t0 = time.time()
    rda_files = sorted(PAIRWISE_DIR.glob('*.rda'))
    print(f"  Found {len(rda_files)} RDA files")

    results = []
    n_loaded = 0

    for rda in rda_files:
        review = load_review(rda)
        if review is None:
            continue

        yi, sei = review['yi'], review['sei']
        k = review['k']
        n_loaded += 1

        # Compute all information-theoretic metrics
        H = shannon_entropy(yi, sei)
        KL = kl_divergence_from_pooled(yi, sei, review['theta'], review['se_theta'], review['tau2'])
        MI = mutual_information_effect_precision(yi, sei)
        FI = fisher_information(sei)
        NEI = normalized_entropy_index(H, k, sei)
        n_modes = detect_multimodality(yi, sei)
        skew, kurt = skewness_kurtosis(yi, sei)

        row = {
            'review_id': review['review_id'],
            'k': k,
            'scale': review['scale'],
            'theta': round(review['theta'], 4),
            'I2': round(review['I2'], 1),
            'tau2': round(review['tau2'], 4),
            'entropy': round(H, 4) if H is not None else '',
            'kl_divergence': round(KL, 4) if KL is not None else '',
            'mutual_info': round(MI, 4) if MI is not None else '',
            'fisher_info': round(FI, 2),
            'norm_entropy_index': round(NEI, 4) if NEI is not None else '',
            'n_modes': n_modes,
            'skewness': round(skew, 3),
            'kurtosis': round(kurt, 3),
            # Entropy-based classification
            'entropy_class': '',
        }

        # Classify
        if NEI is not None:
            if NEI < 0.3:
                row['entropy_class'] = 'Concentrated'
            elif NEI < 0.5:
                row['entropy_class'] = 'Moderate'
            elif NEI < 0.7:
                row['entropy_class'] = 'Dispersed'
            else:
                row['entropy_class'] = 'Chaotic'

        results.append(row)

    elapsed = time.time() - t0
    print(f"  Loaded: {n_loaded} reviews in {elapsed:.1f}s")

    # ═══════════════════════════════════════
    # KEY ANALYSIS: Does entropy capture what I² misses?
    # ═══════════════════════════════════════

    # Spearman correlation: entropy vs I²
    valid = [r for r in results if r['entropy'] != '' and r['I2'] > 0]
    if len(valid) > 10:
        from scipy.stats import spearmanr
        ent_vals = [r['entropy'] for r in valid]
        i2_vals = [r['I2'] for r in valid]
        rho_ent_i2, p_ent_i2 = spearmanr(ent_vals, i2_vals)

        nei_vals = [r['norm_entropy_index'] for r in valid]
        rho_nei_i2, p_nei_i2 = spearmanr(nei_vals, i2_vals)

        # Multimodality vs I²
        multi = [r for r in valid if r['n_modes'] >= 2]
        uni = [r for r in valid if r['n_modes'] == 1]

    # ═══════════════════════════════════════
    # HEADLINE FINDINGS
    # ═══════════════════════════════════════

    from collections import Counter
    ent_classes = Counter(r['entropy_class'] for r in results if r['entropy_class'])
    multi_count = sum(1 for r in results if r['n_modes'] >= 2)
    skewed = sum(1 for r in results if abs(r['skewness']) > 1)

    print(f"\n{'='*50}")
    print("HEADLINE FINDINGS")
    print(f"{'='*50}")
    print(f"  Reviews analyzed: {len(results)}")
    print(f"  Entropy classes: {dict(ent_classes)}")
    print(f"  Multimodal (>=2 modes): {multi_count} ({100*multi_count/len(results):.1f}%)")
    print(f"  Skewed (|skew|>1): {skewed} ({100*skewed/len(results):.1f}%)")

    if len(valid) > 10:
        print(f"\n  Spearman(Entropy, I2): rho={rho_ent_i2:.3f}, p={p_ent_i2:.4f}")
        print(f"  Spearman(NEI, I2):     rho={rho_nei_i2:.3f}, p={p_nei_i2:.4f}")

        # Key insight: reviews where I² is low but entropy is high (hidden heterogeneity)
        hidden_het = [r for r in valid if r['I2'] < 30 and r['norm_entropy_index'] > 0.5]
        print(f"\n  HIDDEN HETEROGENEITY (I2<30% but NEI>0.5): {len(hidden_het)} reviews")
        print(f"  These reviews appear homogeneous by I2 but have dispersed/chaotic effect distributions")

        # And the reverse: I² high but entropy low (artefactual heterogeneity)
        artefact = [r for r in valid if r['I2'] > 70 and r['norm_entropy_index'] < 0.3]
        print(f"  ARTEFACTUAL HETEROGENEITY (I2>70% but NEI<0.3): {len(artefact)} reviews")
        print(f"  These reviews look heterogeneous by I2 but effects are actually concentrated (one outlier inflates Q)")

        if multi_count > 0:
            # Multimodal reviews: I² vs entropy
            multi_i2 = [r['I2'] for r in results if r['n_modes'] >= 2 and r['I2'] > 0]
            uni_i2 = [r['I2'] for r in results if r['n_modes'] == 1 and r['I2'] > 0]
            print(f"\n  Multimodal reviews: mean I2={np.mean(multi_i2):.1f}%")
            print(f"  Unimodal reviews:   mean I2={np.mean(uni_i2):.1f}%")

    # ═══════════════════════════════════════
    # EXPORT
    # ═══════════════════════════════════════

    fields = list(results[0].keys())
    with open(OUTPUT_DIR / 'entropy_results.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)

    summary = {
        'n_reviews': len(results),
        'entropy_classes': dict(ent_classes),
        'n_multimodal': multi_count,
        'pct_multimodal': round(100 * multi_count / len(results), 1),
        'n_skewed': skewed,
        'pct_skewed': round(100 * skewed / len(results), 1),
        'correlation_entropy_i2': round(rho_ent_i2, 3) if len(valid) > 10 else None,
        'correlation_nei_i2': round(rho_nei_i2, 3) if len(valid) > 10 else None,
        'n_hidden_heterogeneity': len(hidden_het) if len(valid) > 10 else 0,
        'n_artefactual_heterogeneity': len(artefact) if len(valid) > 10 else 0,
        'elapsed_seconds': round(elapsed, 1),
    }
    with open(OUTPUT_DIR / 'entropy_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Saved to {OUTPUT_DIR}/")


if __name__ == '__main__':
    main()
