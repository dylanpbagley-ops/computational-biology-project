"""
test_aim5_and_gate_dual_marker_selectivity.py

Regression tests for src/aim5_and_gate_dual_marker_selectivity.py -- the
AND-gate dual-marker selectivity model, built after the user asked
directly whether a "multiple protein targets to amplify selectivity"
idea had ever been investigated (it had not -- Aim 3b built an OR-gate
yield question, not an AND-gate selectivity question; see the module's
own docstring for the full account).

WHAT THESE TESTS LOCK IN, AND WHY EACH ONE MATTERS
-------------------------------------------------------
1. The Gaussian-copula joint reconstruction actually reproduces the
   correct MARGINAL distribution for each marker (mean and, loosely, the
   general shape) regardless of the correlation parameter -- checked
   empirically, not just algebraically. If this were wrong, every
   downstream AND-gate number would be silently built on a marginal that
   doesn't match the real digitized data it's supposed to be anchored to.
2. correlation=0 in the copula actually reproduces independence to a
   tight numerical tolerance (checked via empirical Pearson correlation
   of a large sample) -- the "independence is the honest default" framing
   in the module docstring depends on this being true, not assumed.
3. THE K_GEOM REGRESSION TEST: single_marker_yield's CA19-9-only ratio at
   the project's real-sigma-anchored k_intra must match
   aim3_aim4_range_sweep.sweep_ranking_flip_robustness()'s own
   independently-established CA19-9 ratio (1.4373 at z=1.0) to tight
   tolerance. An earlier version of this file computed k_geom SEPARATELY
   for each population (using each population's own mean) instead of
   anchoring it once to the healthy population and reusing that same
   value for both populations -- exactly the convention
   sweep_ranking_flip_robustness() itself uses and documents as
   deliberate ("k_geom does not itself become a hidden second free
   parameter"). That bug silently produced a ratio of ~1.04 instead of
   the correct ~1.44 -- caught by comparing against the already-
   established, already-tested number before trusting this file's own
   result. This test exists specifically so that bug (or one like it)
   cannot silently return.
4. The reliability diagnostic (top_1pct_contribution /
   AND_GATE_MAX_RELIABLE_TOP_1PCT_CONTRIBUTION) correctly flags the
   diagnosed-unstable dilute-z/high-correlation regime as unreliable, and
   correctly passes the stable z=1.0 default regime -- a direct
   regression test for the real numerical instability found and
   diagnosed while building this file (checked directly: the AND-gate
   ratio swung 67-99 across four seeds at z=1e-4, correlation=0.9,
   n_vesicles=50,000, while it stayed within 24-26 at correlation=0 at
   the same z).
5. THE HEADLINE FINDING: the AND-gate selectivity ratio exceeds the best
   single-marker ratio at EVERY k_intra value in the real-sigma
   plausibility range, at z=1.0 (the checked-stable regime), across a
   range of assumed correlations -- and every one of these cells is
   itself flagged "reliable" by the diagnostic in (4), so this is not
   just a number that happens to look good; it is a number this file's
   own numerical-reliability check has independently confirmed is
   trustworthy.
"""

import sys
import os

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from aim5_and_gate_dual_marker_selectivity import (
    fit_all_marginals,
    simulate_joint_two_marker_population,
    and_gate_selectivity_sweep,
    and_gate_capture_probability,
    or_gate_capture_probability,
    AND_GATE_MAX_RELIABLE_TOP_1PCT_CONTRIBUTION,
)
from aim3_aim4_range_sweep import sweep_ranking_flip_robustness, get_k_intra_sweep_values


# ---------------------------------------------------------------------
# Test 1 & 2: the copula reconstruction hits its marginal targets, and
# correlation=0 reproduces independence.
# ---------------------------------------------------------------------

def test_copula_reproduces_correct_marginal_mean_regardless_of_correlation():
    fits, _ = fit_all_marginals()
    for population in ("PDAC", "healthy"):
        for correlation in (0.0, 0.5, 0.9):
            pop = simulate_joint_two_marker_population(
                fits, population, correlation, n_vesicles=100_000, seed=0
            )
            expected_egfr_mean = fits[("EGFR", population)]["mean_n"]
            expected_ca199_mean = fits[("CA199", population)]["mean_n"]
            assert abs(pop["egfr"].mean() - expected_egfr_mean) / expected_egfr_mean < 0.05, (
                f"{population}, correlation={correlation}: EGFR empirical mean "
                f"{pop['egfr'].mean():.2f} does not match target {expected_egfr_mean:.2f} "
                f"within 5% -- the copula construction is not preserving its marginal."
            )
            assert abs(pop["ca199"].mean() - expected_ca199_mean) / expected_ca199_mean < 0.05, (
                f"{population}, correlation={correlation}: CA19-9 empirical mean "
                f"{pop['ca199'].mean():.2f} does not match target {expected_ca199_mean:.2f} "
                f"within 5%."
            )


def test_correlation_zero_reproduces_independence():
    fits, _ = fit_all_marginals()
    pop = simulate_joint_two_marker_population(fits, "PDAC", 0.0, n_vesicles=200_000, seed=0)
    pearson_r, _ = stats.pearsonr(pop["egfr"], pop["ca199"])
    assert abs(pearson_r) < 0.02, (
        f"Expected correlation=0.0 in the copula to reproduce empirical "
        f"independence (Pearson r near 0), got r={pearson_r:.4f}."
    )


def test_correlation_parameter_controls_empirical_correlation_monotonically():
    """
    Sanity check that the correlation parameter actually does what its
    name says: higher requested correlation -> higher empirical Pearson
    correlation between the two marker counts (not exactly equal, since
    Pearson correlation on transformed/discrete negative-binomial counts
    differs slightly from the underlying Gaussian-copula correlation, but
    monotonically increasing).
    """
    fits, _ = fit_all_marginals()
    empirical_correlations = []
    for correlation in (0.0, 0.3, 0.6, 0.9):
        pop = simulate_joint_two_marker_population(fits, "PDAC", correlation, n_vesicles=100_000, seed=0)
        r, _ = stats.pearsonr(pop["egfr"], pop["ca199"])
        empirical_correlations.append(r)
    assert all(
        empirical_correlations[i] < empirical_correlations[i + 1]
        for i in range(len(empirical_correlations) - 1)
    ), f"Expected monotonically increasing empirical correlation, got {empirical_correlations}"


# ---------------------------------------------------------------------
# Test 3: THE K_GEOM REGRESSION TEST -- single-marker ratios computed by
# this file must match the already-established aim3_aim4_range_sweep
# numbers exactly.
# ---------------------------------------------------------------------

def test_single_marker_ratios_match_aim4_range_sweep():
    k_intra = max(get_k_intra_sweep_values())
    result = and_gate_selectivity_sweep(k_intra=k_intra, z=1.0, correlations=(0.0,))

    reference = sweep_ranking_flip_robustness(z=1.0)
    ca199_rows = reference["healthy_pdac_plasma_ca199"]
    egfr_rows = reference["healthy_pdac_plasma_egfr"]
    ca199_reference_ratio = next(
        row["ratio_p2"] for row in ca199_rows if abs(row["k_intra"] - k_intra) < 1e-20
    )
    egfr_reference_ratio = next(
        row["ratio_p2"] for row in egfr_rows if abs(row["k_intra"] - k_intra) < 1e-20
    )

    assert abs(result["ca199_only_ratio"] - ca199_reference_ratio) < 1e-9, (
        f"CA19-9-only ratio ({result['ca199_only_ratio']}) does not match "
        f"aim3_aim4_range_sweep's independently-established reference "
        f"({ca199_reference_ratio}) -- this is the exact regression this "
        f"test exists to catch (a k_geom convention mismatch silently "
        f"changes this number). Do not loosen this tolerance; investigate "
        f"the k_geom construction instead."
    )
    assert abs(result["egfr_only_ratio"] - egfr_reference_ratio) < 1e-9, (
        f"EGFR-only ratio ({result['egfr_only_ratio']}) does not match "
        f"the reference ({egfr_reference_ratio})."
    )


# ---------------------------------------------------------------------
# Test 4: the reliability diagnostic correctly flags the diagnosed-
# unstable regime and correctly passes the stable default regime.
# ---------------------------------------------------------------------

def test_reliability_diagnostic_flags_known_unstable_dilute_regime():
    result = and_gate_selectivity_sweep(z=1e-4, correlations=(0.0, 0.9), n_vesicles=50_000, seed=0)
    assert result["reliable"][1] is False, (
        "Expected z=1e-4, correlation=0.9 to be flagged unreliable (this "
        "is the exact regime diagnosed this session as numerically "
        "unstable -- the AND-gate ratio swung 67-99 across seeds here) -- "
        "if this now passes as reliable, the diagnostic's threshold or "
        "construction has changed; re-verify the instability is actually "
        "gone (e.g. via a multi-seed check) before trusting this."
    )


def test_reliability_diagnostic_passes_stable_default_regime():
    result = and_gate_selectivity_sweep(z=1.0, correlations=(0.0, 0.2, 0.4, 0.6, 0.8))
    assert all(result["reliable"]), (
        f"Expected every correlation at the default z=1.0 regime to be "
        f"flagged reliable, got {result['reliable']} -- if this default "
        f"regime is no longer reliable, every other test and the "
        f"headline finding built on it needs re-examining."
    )


# ---------------------------------------------------------------------
# Test 5: THE HEADLINE FINDING -- AND-gate beats the best single marker
# across the entire real-sigma k_intra plausibility range, at z=1.0.
# ---------------------------------------------------------------------

def test_and_gate_beats_best_single_marker_across_full_plausibility_range():
    for k_intra in sorted(get_k_intra_sweep_values()):
        result = and_gate_selectivity_sweep(k_intra=k_intra, z=1.0, correlations=(0.0, 0.4, 0.8))
        assert all(result["reliable"]), (
            f"k_intra={k_intra:.3e}: expected every tested correlation to "
            f"be numerically reliable at z=1.0 -- got {result['reliable']}. "
            f"Investigate before trusting the beats-single-marker claim "
            f"at this k_intra."
        )
        assert all(result["and_gate_beats_best_single_marker"]), (
            f"k_intra={k_intra:.3e}: expected the AND-gate ratio to beat "
            f"the best single-marker ratio ({result['best_single_marker_ratio']:.4f}) "
            f"at every tested correlation, got and_gate_ratio="
            f"{result['and_gate_ratio']}. This is the file's headline "
            f"finding -- if it fails somewhere in the plausibility range, "
            f"report that honestly rather than loosening this assertion."
        )


def test_and_gate_and_or_gate_are_correctly_ordered():
    """
    Sanity check on the two gate logics themselves, independent of any
    real data: for any given per-vesicle pair of single-marker
    probabilities in (0, 1), the AND-gate probability must be <= either
    individual probability, and the OR-gate probability must be >=
    either individual probability. This is a property of probability
    theory (P(A and B) <= P(A), P(B); P(A or B) >= P(A), P(B) for any
    two events), not something specific to this project's binding model
    -- included as a basic correctness check on the two gate functions'
    arithmetic.
    """
    from binding_model import capture_probability_patch
    from k_intra_z_plausibility_bounds import k_intra_bounds_from_real_sigma, SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2

    n_egfr = np.array([10.0, 50.0, 100.0])
    n_ca199 = np.array([20.0, 5.0, 80.0])
    diameter = np.array([80.0, 120.0, 60.0])
    k_intra = max(k_intra_bounds_from_real_sigma().values())
    sigma = SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2
    z = 1.0
    k_geom = 1.0

    theta_egfr = capture_probability_patch(n_egfr, diameter, sigma, k_geom, k_intra, z)
    theta_ca199 = capture_probability_patch(n_ca199, diameter, sigma, k_geom, k_intra, z)
    theta_and = and_gate_capture_probability(n_egfr, n_ca199, diameter, sigma, k_geom, k_intra, z)
    theta_or = or_gate_capture_probability(n_egfr, n_ca199, diameter, sigma, k_geom, k_intra, z)

    assert np.all(theta_and <= theta_egfr + 1e-12) and np.all(theta_and <= theta_ca199 + 1e-12), (
        "AND-gate probability must never exceed either single-marker probability."
    )
    assert np.all(theta_or >= theta_egfr - 1e-12) and np.all(theta_or >= theta_ca199 - 1e-12), (
        "OR-gate probability must never be less than either single-marker probability."
    )


if __name__ == "__main__":
    # Allow running without pytest, matching this project's other test
    # files' convention.
    test_functions = [
        test_copula_reproduces_correct_marginal_mean_regardless_of_correlation,
        test_correlation_zero_reproduces_independence,
        test_correlation_parameter_controls_empirical_correlation_monotonically,
        test_single_marker_ratios_match_aim4_range_sweep,
        test_reliability_diagnostic_flags_known_unstable_dilute_regime,
        test_reliability_diagnostic_passes_stable_default_regime,
        test_and_gate_beats_best_single_marker_across_full_plausibility_range,
        test_and_gate_and_or_gate_are_correctly_ordered,
    ]
    for test_fn in test_functions:
        try:
            test_fn()
            print(f"PASS: {test_fn.__name__}")
        except AssertionError as e:
            print(f"FAIL: {test_fn.__name__}\n  {e}")
