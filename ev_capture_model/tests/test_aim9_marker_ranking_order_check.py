"""
test_aim9_marker_ranking_order_check.py

Regression tests for src/aim9_marker_ranking_order_check.py -- the
project's first "new claim" check (does the nonlinear model, fed real
marginal data, predict a real, externally-published AUC ranking?). See
that file's own module docstring for the full account of what data
actually exists, what Check 1 and Check 2 each honestly do and don't
claim, and why Check 2 reports an invalid ratio-transfer as a finding
rather than forcing a number through it.

WHAT THESE TESTS LOCK IN, AND WHY EACH ONE MATTERS
-------------------------------------------------------
1. The real ground-truth data (Guo et al.'s AUCs, the four real CRC
   positivity fractions) is correctly exposed and internally consistent
   with aim7/aim8's own data (no silent divergence between what aim9
   imports and what aim7/aim8 actually define).
2. auc_from_two_samples() is mathematically correct -- checked against
   hand-constructed samples with a KNOWN true AUC (perfect separation,
   perfect anti-separation, and no separation), not just checked for being
   "in range."
3. linear_ranking_from_real_positivity() matches an independent hand
   computation from crc_empirical_marginal_positivity()'s own real data.
4. nonlinear_yield_ranking() returns valid probabilities and a ranking
   over exactly the same four markers.
5. sweep_reordering_frequency() is internally consistent (n_reordered /
   n_grid_points arithmetic is correct; every grid point's ranking is a
   permutation of the same four markers; the itgam_rank_positions list
   has one entry per grid point, all in [1, 4]).
6. itgam_cluster_ratio_transfer_check() correctly reproduces the REAL
   caught finding -- the naive transfer IS invalid (produces a fraction
   above 1.0) for the real numbers currently in this project. If this
   test ever finds is_valid=True, that means the underlying real numbers
   changed and this file's Check 2 narrative needs re-writing, not that
   the test should just be loosened.
7. print_summary() runs end-to-end without raising.
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from aim9_marker_ranking_order_check import (
    MARKERS_WITH_REAL_CRC_POSITIVITY,
    MARKER_WITH_BOTH_REAL_AUC_AND_POSITIVITY,
    GENERIC_MEAN_COPIES_IF_POSITIVE,
    _reconstruct_crc_population,
    auc_from_two_samples,
    linear_ranking_from_real_positivity,
    nonlinear_yield_ranking,
    sweep_reordering_frequency,
    itgam_cluster_ratio_transfer_check,
    print_summary,
)
from aim8_raw_vesicle_correlation_check import crc_empirical_marginal_positivity, CRC_MARKERS
from aim7_crc_plasma_exosome_cluster_replication import (
    GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC,
    GUO_ET_AL_CLUSTER11_HEALTHY_TO_CANCER_DECLINE,
)
from k_intra_z_plausibility_bounds import SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2


# ---------------------------------------------------------------------
# Test 1: real ground-truth data is exposed correctly and consistently.
# ---------------------------------------------------------------------

def test_markers_with_real_crc_positivity_matches_aim8():
    assert MARKERS_WITH_REAL_CRC_POSITIVITY == CRC_MARKERS, (
        "aim9's MARKERS_WITH_REAL_CRC_POSITIVITY should be exactly aim8's "
        "CRC_MARKERS (ITGB3, ITGAM, ITGAL, ITGB2) -- if these have "
        "diverged, aim9 may be using a stale or incorrect marker list."
    )


def test_guo_auc_data_has_four_markers_including_itgam():
    assert set(GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC.keys()) == {
        "ITGAM", "ITGB3", "EPS8", "CD9"
    }
    assert MARKER_WITH_BOTH_REAL_AUC_AND_POSITIVITY == "ITGAM"
    assert MARKER_WITH_BOTH_REAL_AUC_AND_POSITIVITY in GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC
    assert MARKER_WITH_BOTH_REAL_AUC_AND_POSITIVITY in MARKERS_WITH_REAL_CRC_POSITIVITY
    # EPS8 and CD9 must NOT have real positivity data anywhere in this
    # project -- confirmed directly against the raw CRC file's own columns.
    assert "EPS8" not in MARKERS_WITH_REAL_CRC_POSITIVITY
    assert "CD9" not in MARKERS_WITH_REAL_CRC_POSITIVITY
    # ITGAM should have the HIGHEST real AUC among the four Guo markers --
    # this is a real fact about the ground-truth data, checked here so any
    # future edit to GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC that silently
    # breaks this is caught.
    best_marker = max(GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC,
                       key=GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC.get)
    assert best_marker == "ITGAM"


# ---------------------------------------------------------------------
# Test 2: auc_from_two_samples() is mathematically correct.
# ---------------------------------------------------------------------

def test_auc_perfect_separation():
    sample_a = np.array([10.0, 11.0, 12.0])
    sample_b = np.array([1.0, 2.0, 3.0])
    result = auc_from_two_samples(sample_a, sample_b)
    assert abs(result["auc"] - 1.0) < 1e-9, (
        f"Every value in sample_a exceeds every value in sample_b -- "
        f"expected AUC=1.0 exactly, got {result['auc']}."
    )


def test_auc_perfect_anti_separation():
    sample_a = np.array([1.0, 2.0, 3.0])
    sample_b = np.array([10.0, 11.0, 12.0])
    result = auc_from_two_samples(sample_a, sample_b)
    assert abs(result["auc"] - 0.0) < 1e-9, (
        f"Every value in sample_a is LESS than every value in sample_b -- "
        f"expected AUC=0.0 exactly, got {result['auc']}."
    )


def test_auc_no_separation_with_identical_distributions():
    rng = np.random.default_rng(42)
    sample_a = rng.normal(loc=5.0, scale=1.0, size=2000)
    sample_b = rng.normal(loc=5.0, scale=1.0, size=2000)
    result = auc_from_two_samples(sample_a, sample_b)
    assert abs(result["auc"] - 0.5) < 0.03, (
        f"Two samples from the SAME distribution should give AUC close to "
        f"0.5 (random/no separation), got {result['auc']:.4f}."
    )


def test_auc_bounded_and_matches_manual_u_statistic():
    rng = np.random.default_rng(7)
    sample_a = rng.uniform(0, 1, size=50)
    sample_b = rng.uniform(0, 1, size=40)
    result = auc_from_two_samples(sample_a, sample_b)
    assert 0.0 <= result["auc"] <= 1.0
    assert 0.0 <= result["p_value"] <= 1.0
    # Manual, brute-force AUC computation (independent of scipy's U
    # statistic path) -- counts, for every pair, whether a > b.
    count = 0
    for a in sample_a:
        for b in sample_b:
            if a > b:
                count += 1
            elif a == b:
                count += 0.5
    manual_auc = count / (len(sample_a) * len(sample_b))
    assert abs(result["auc"] - manual_auc) < 1e-9, (
        f"auc_from_two_samples() gives {result['auc']}, independent "
        f"brute-force pairwise computation gives {manual_auc}."
    )


# ---------------------------------------------------------------------
# Test 3: linear ranking matches real data directly.
# ---------------------------------------------------------------------

def test_linear_ranking_matches_real_marginal_positivity():
    marg = crc_empirical_marginal_positivity()
    result = linear_ranking_from_real_positivity()
    for marker in MARKERS_WITH_REAL_CRC_POSITIVITY:
        assert result["positivity_fractions"][marker] == marg[marker], (
            f"{marker}: aim9's linear_ranking_from_real_positivity() "
            f"reports {result['positivity_fractions'][marker]}, but "
            f"aim8's crc_empirical_marginal_positivity() (the real data "
            f"source) gives {marg[marker]}."
        )
    # Ranking must be a valid permutation of the four markers, sorted
    # descending by the real fractions.
    assert set(result["ranking_best_first"]) == set(MARKERS_WITH_REAL_CRC_POSITIVITY)
    fractions_in_rank_order = [
        result["positivity_fractions"][m] for m in result["ranking_best_first"]
    ]
    assert fractions_in_rank_order == sorted(fractions_in_rank_order, reverse=True), (
        "ranking_best_first should list markers in descending order of "
        "real positivity fraction."
    )


# ---------------------------------------------------------------------
# Test 4: nonlinear ranking is valid.
# ---------------------------------------------------------------------

def test_nonlinear_yield_ranking_valid_probabilities_and_markers():
    result = nonlinear_yield_ranking(
        sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, k_intra=1e-13, z=1e-4
    )
    assert set(result["yields"].keys()) == set(MARKERS_WITH_REAL_CRC_POSITIVITY)
    for marker, yield_val in result["yields"].items():
        assert 0.0 <= yield_val <= 1.0, (
            f"{marker}: capture probability yield {yield_val} out of "
            f"bounds [0, 1]."
        )
    assert set(result["ranking_best_first"]) == set(MARKERS_WITH_REAL_CRC_POSITIVITY)
    yields_in_rank_order = [result["yields"][m] for m in result["ranking_best_first"]]
    assert yields_in_rank_order == sorted(yields_in_rank_order, reverse=True)


def test_reconstructed_population_mean_roughly_matches_target_positivity():
    """
    Sanity check on _reconstruct_crc_population(): the fraction of
    synthetic vesicles with n > 0 should be close to the real target
    positivity fraction it was built from (not exact, due to sampling
    noise, but close at n_samples=5000).
    """
    marg = crc_empirical_marginal_positivity()
    for marker in MARKERS_WITH_REAL_CRC_POSITIVITY:
        target = marg[marker]
        n_values = _reconstruct_crc_population(marker, target, n_samples=5000, seed=0)
        empirical_positivity = float((n_values > 0).mean())
        assert abs(empirical_positivity - target) < 0.03, (
            f"{marker}: target positivity {target:.4f}, but reconstructed "
            f"population's empirical positivity is {empirical_positivity:.4f} "
            f"-- fit_negative_binomial_to_positivity may not be matching "
            f"its target closely enough."
        )


# ---------------------------------------------------------------------
# Test 5: sweep_reordering_frequency() internal consistency.
# ---------------------------------------------------------------------

def test_sweep_reordering_frequency_internal_consistency():
    result = sweep_reordering_frequency(n_samples=1000)  # smaller n_samples
    # for test speed -- the reordering-frequency FINDING itself is
    # reported by print_summary()/the module's own __main__ run at the
    # full n_samples=5000; this test only checks internal arithmetic
    # consistency, which does not depend on sample size.

    assert result["n_grid_points"] == len(result["grid_results"])
    assert result["n_grid_points"] > 0

    manual_n_reordered = sum(
        1 for g in result["grid_results"] if not g["matches_linear"]
    )
    assert manual_n_reordered == result["n_reordered"], (
        f"Reported n_reordered={result['n_reordered']} does not match a "
        f"direct recount from grid_results ({manual_n_reordered})."
    )

    expected_fraction = result["n_reordered"] / result["n_grid_points"]
    assert abs(result["reordering_fraction"] - expected_fraction) < 1e-9

    assert 0.0 <= result["reordering_fraction"] <= 1.0

    # Every grid point's nonlinear_ranking must be a permutation of the
    # four real-data markers.
    for g in result["grid_results"]:
        assert set(g["nonlinear_ranking"]) == set(MARKERS_WITH_REAL_CRC_POSITIVITY)
        assert g["matches_linear"] == (
            g["nonlinear_ranking"] == result["linear_ranking"]["ranking_best_first"]
        )

    # itgam_rank_positions: one entry per grid point, every value in [1, 4].
    assert len(result["itgam_rank_positions"]) == result["n_grid_points"]
    for pos in result["itgam_rank_positions"]:
        assert 1 <= pos <= len(MARKERS_WITH_REAL_CRC_POSITIVITY)


# ---------------------------------------------------------------------
# Test 6: itgam_cluster_ratio_transfer_check() correctly reproduces the
# REAL caught invalid-ratio-transfer finding.
# ---------------------------------------------------------------------

def test_itgam_ratio_transfer_is_correctly_found_invalid():
    """
    THE HEADLINE CHECK 2 FINDING: for the real numbers currently in this
    project (ITGAM's real CRC positivity ~14.62%, the cluster's real
    13.47%/1.13% healthy/CRC decline), naively transferring the cluster's
    healthy:CRC ratio to ITGAM's own marginal produces a fraction ABOVE
    1.0 -- confirmed here directly (not just trusting the module's own
    is_valid flag) via an independent hand recomputation. If this test
    ever finds the transferred fraction is LEGITIMATELY <= 1.0, that means
    the real underlying numbers changed and this file's Check 2 narrative
    (and its is_valid=False framing throughout print_summary()) needs to
    be revisited -- not that this test should just be loosened to pass.
    """
    marg = crc_empirical_marginal_positivity()
    decline = GUO_ET_AL_CLUSTER11_HEALTHY_TO_CANCER_DECLINE

    hand_ratio = decline["healthy_fraction"] / decline["primary_and_metastatic_crc_fraction"]
    hand_transferred = marg["ITGAM"] * hand_ratio

    result = itgam_cluster_ratio_transfer_check()

    assert abs(result["naive_decline_ratio"] - hand_ratio) < 1e-9
    assert abs(result["naive_transferred_fraction"] - hand_transferred) < 1e-9
    assert abs(result["itgam_real_crc_fraction"] - marg["ITGAM"]) < 1e-9
    assert result["cluster_healthy_fraction"] == decline["healthy_fraction"]
    assert result["cluster_crc_fraction"] == decline["primary_and_metastatic_crc_fraction"]

    assert result["naive_transferred_fraction"] > 1.0, (
        f"Expected the naive ratio-transfer to produce an IMPOSSIBLE "
        f"fraction above 1.0 (this is the real, caught finding this "
        f"function exists to report) -- got "
        f"{result['naive_transferred_fraction']:.4f}, which would be a "
        f"valid probability. This means the real underlying numbers "
        f"changed; the Check 2 narrative needs to be revisited, not this "
        f"test loosened."
    )
    assert result["is_valid"] is False


def test_itgam_ratio_transfer_check_does_not_fabricate_a_healthy_population():
    """
    Confirms itgam_cluster_ratio_transfer_check() returns only the
    arithmetic/diagnostic numbers, and does NOT return any key suggesting
    a synthetic healthy population or AUC was constructed from this
    invalid approach (e.g. no "synthetic_healthy_fraction" or "auc_result"
    key) -- guards against a future edit silently reintroducing the
    abandoned, invalid approach.
    """
    result = itgam_cluster_ratio_transfer_check()
    forbidden_keys = {"synthetic_healthy_fraction", "auc_result", "auc"}
    assert not (forbidden_keys & set(result.keys())), (
        f"itgam_cluster_ratio_transfer_check() should only report the "
        f"diagnostic arithmetic proving the transfer is invalid, not "
        f"construct a synthetic healthy population or AUC from it -- "
        f"found unexpected key(s): {forbidden_keys & set(result.keys())}."
    )


# ---------------------------------------------------------------------
# Test 7: print_summary() runs end-to-end without raising.
# ---------------------------------------------------------------------

def test_print_summary_runs_without_error():
    print_summary()  # Just confirms no exception -- output correctness is
    # covered by the more targeted tests above.


if __name__ == "__main__":
    test_functions = [
        test_markers_with_real_crc_positivity_matches_aim8,
        test_guo_auc_data_has_four_markers_including_itgam,
        test_auc_perfect_separation,
        test_auc_perfect_anti_separation,
        test_auc_no_separation_with_identical_distributions,
        test_auc_bounded_and_matches_manual_u_statistic,
        test_linear_ranking_matches_real_marginal_positivity,
        test_nonlinear_yield_ranking_valid_probabilities_and_markers,
        test_reconstructed_population_mean_roughly_matches_target_positivity,
        test_sweep_reordering_frequency_internal_consistency,
        test_itgam_ratio_transfer_is_correctly_found_invalid,
        test_itgam_ratio_transfer_check_does_not_fabricate_a_healthy_population,
        test_print_summary_runs_without_error,
    ]
    for test_fn in test_functions:
        try:
            test_fn()
            print(f"PASS: {test_fn.__name__}")
        except AssertionError as e:
            print(f"FAIL: {test_fn.__name__}\n  {e}")
