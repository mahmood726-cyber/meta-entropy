"""Tests for Meta-analytic Entropy pipeline."""
import sys, math, pytest, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline import (shannon_entropy, kl_divergence_from_pooled,
                       mutual_information_effect_precision, fisher_information,
                       normalized_entropy_index, detect_multimodality, skewness_kurtosis)


class TestShannonEntropy:
    def test_homogeneous_low_entropy(self):
        """Identical studies should have low entropy."""
        yi = np.array([-0.5, -0.5, -0.5, -0.5, -0.5])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        H = shannon_entropy(yi, sei)
        assert H is not None
        assert H < 2.0  # Low entropy for concentrated distribution

    def test_heterogeneous_high_entropy(self):
        """Spread-out studies should have higher entropy."""
        yi_homo = np.array([-0.5, -0.5, -0.5, -0.5, -0.5])
        yi_hetero = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        H_homo = shannon_entropy(yi_homo, sei)
        H_hetero = shannon_entropy(yi_hetero, sei)
        assert H_hetero > H_homo

    def test_k2_returns_none(self):
        assert shannon_entropy(np.array([0.1, 0.2]), np.array([0.1, 0.1])) is None

    def test_entropy_is_finite(self):
        yi = np.array([-0.3, 0.1, 0.5, -0.1, 0.3])
        sei = np.array([0.2, 0.15, 0.25, 0.1, 0.3])
        H = shannon_entropy(yi, sei)
        assert H is not None
        assert math.isfinite(H)


class TestKLDivergence:
    def test_perfect_pooling_low_kl(self):
        """When all studies agree, KL should be near zero."""
        yi = np.array([-0.5, -0.5, -0.5, -0.5, -0.5])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        KL = kl_divergence_from_pooled(yi, sei, -0.5, 0.045, 0.0)
        assert KL is not None
        assert KL < 2.0  # Low-ish information loss (MC noise inflates slightly)

    def test_kl_non_negative(self):
        """KL divergence is always >= 0."""
        yi = np.array([-1.0, 0.0, 1.0, -0.5, 0.5])
        sei = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
        KL = kl_divergence_from_pooled(yi, sei, 0.0, 0.1, 0.5)
        assert KL is not None
        assert KL >= 0


class TestMutualInformation:
    def test_independent_low_mi(self):
        """Random effects+precision should have near-zero MI."""
        rng = np.random.RandomState(42)
        yi = rng.normal(0, 1, 20)
        sei = np.abs(rng.normal(0.2, 0.05, 20)) + 0.01
        MI = mutual_information_effect_precision(yi, sei)
        assert MI is not None
        assert MI < 0.6  # Near-zero but binning noise with k=20

    def test_correlated_higher_mi(self):
        """Small studies with large effects should have higher MI."""
        # Simulate small-study effect
        sei = np.array([0.5, 0.4, 0.3, 0.2, 0.1, 0.08, 0.06, 0.05, 0.04, 0.03])
        yi = -0.1 + 0.5 / sei  # Larger effect with smaller SE (Egger-like)
        MI = mutual_information_effect_precision(yi, sei)
        assert MI is not None
        assert MI > 0


class TestFisherInformation:
    def test_precision_sum(self):
        sei = np.array([0.1, 0.2, 0.5])
        FI = fisher_information(sei)
        expected = 1/0.01 + 1/0.04 + 1/0.25
        assert abs(FI - expected) < 1e-10

    def test_more_precise_more_info(self):
        sei_precise = np.array([0.05, 0.05, 0.05])
        sei_imprecise = np.array([0.5, 0.5, 0.5])
        assert fisher_information(sei_precise) > fisher_information(sei_imprecise)


class TestMultimodality:
    def test_unimodal(self):
        yi = np.array([-0.5, -0.4, -0.3, -0.2, -0.1])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        assert detect_multimodality(yi, sei) == 1

    def test_bimodal(self):
        """Two clusters should detect 2 modes."""
        yi = np.array([-2.0, -1.9, -1.8, 1.8, 1.9, 2.0])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
        modes = detect_multimodality(yi, sei)
        assert modes >= 2


class TestSkewness:
    def test_symmetric_near_zero(self):
        yi = np.array([-0.3, -0.1, 0.0, 0.1, 0.3])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        skew, _ = skewness_kurtosis(yi, sei)
        assert abs(skew) < 0.5

    def test_right_skewed(self):
        yi = np.array([0.0, 0.1, 0.1, 0.2, 2.0])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        skew, _ = skewness_kurtosis(yi, sei)
        assert skew > 0.5


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
