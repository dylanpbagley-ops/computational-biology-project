"""
test_stem_cell_marker_ranking_extension.py

Regression tests for src/stem_cell_marker_ranking_extension.py -- the
v11->v12 handoff transition's first extension of the range-sweep treatment
(aim3_aim4_range_sweep.py) to a literature-table system OTHER than the two
real-digitized Lennon comparisons, per v11's "Immediate next steps" item 5.

WHAT THESE TESTS LOCK IN, AND WHY EACH ONE MATTERS
-------------------------------------------------------
1. The negative-binomial fitting procedure (fit_negative_binomial_to_
   positivity) actually hits both of its two targets (the zero-fraction
   AND the positive-conditional mean) to tight tolerance, for every
   marker/tissue combination actually used in this file -- not just
   "runs without crashing." This is the foundation everything else is
   built on; if the fit is wrong, every downstream yield and ranking is
   silently wrong too.

2. stable_joint_population_yield() NEVER returns a non-finite value,
   across the full real-sigma K_INTRA plausibility range, for every
   marker/tissue population used in this file. This is the DIRECT
   regression test for the overflow bug this session found and fixed
   while building this file (see that function's own docstring for the
   full account: the original population_model.joint_population_yield(),
   applied to this file's heavy-tailed reconstructed populations at the
   largest plausibility-range K_INTRA, silently produced nan, which then
   silently produced a MEANINGLESS ranking via Python's sorted() --
   nan comparisons are always False, so sorted() does not raise an error,
   it just returns an arbitrary order). This test is written specifically
   so that if the overflow-safe rewrite is ever "simplified" back to a
   direct call to the original function, or the overflow_cutoff constant
   is loosened incorrectly, this test fails immediately rather than
   silently producing a wrong ranking again.

3. The marker-level ranking crossover (CD81 overtakes CD63 at the largest
   K_INTRA in the real-sigma plausibility range, despite CD63 having the
   higher naive fold-change score at every one of the three tissues) is
   locked in as the corrected, independently-reproduced finding -- not
   just printed once and trusted. This finding was independently
   re-verified across four different random seeds before being reported
   (see the handoff for the numbers) specifically because a single-seed
   result on a reconstructed (not real digitized) population is more
   vulnerable to being a sampling fluke than a result on real data; this
   test checks the SAME crossover reproduces with the file's own default
   (seed=0) reconstruction, deterministically.

4. CD73 stays last in the marker-level ranking at every tested K_INTRA --
   the "always dominated" half of the finding, worth locking in
   separately from the CD81/CD63 crossover since it's a different kind of
   claim (a ranking that does NOT change, rather than one that does).

HOW TO RUN THIS FILE
-----------------------
From the tests/ folder (with src/ next to it):
    pytest test_stem_cell_marker_ranking_extension.py -v
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from stem_cell_marker_ranking_extension import (
    fit_negative_binomial_to_positivity,
    stable_joint_population_yield,
    reconstruct_marker_population,
    naive_fold_change_ranking,
    marker_level_ranking,
    sweep_marker_ranking_across_plausibility_range,
    POSITIVITY_FRACTION_RANGE,
    COPIES_PER_POSITIVE_PARTICLE,
)
from k_intra_z_plausibility_bounds import k_intra_bounds_from_real_sigma, SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2


# ---------------------------------------------------------------------
# Test 1: the negative-binomial fit actually hits both targets, for
# every marker/tissue combination this file actually uses.
# ---------------------------------------------------------------------

def test_negative_binomial_fit_hits_both_targets_for_every_marker_tissue():
    for marker, tissue_copies in COPIES_PER_POSITIVE_PARTICLE.items():
        lo, hi = POSITIVITY_FRACTION_RANGE[marker]
        fraction_positive = 0.5 * (lo + hi)
        target_p0 = 1.0 - fraction_positive
        for tissue, mean_copies_if_positive in tissue_copies.items():
            fit = fit_negative_binomial_to_positivity(fraction_positive, mean_copies_if_positive)

            achieved_p0 = fit["probabilities"][0]
            assert abs(achieved_p0 - target_p0) < 1e-4, (
                f"{marker} ({tissue}): fitted P(n=0)={achieved_p0:.6f} does "
                f"not match target {target_p0:.6f}"
            )

            n_values = fit["n_values"]
            probabilities = fit["probabilities"]
            positive_mask = n_values > 0
            achieved_conditional_mean = (
                np.sum(n_values[positive_mask] * probabilities[positive_mask])
                / probabilities[positive_mask].sum()
            )
            assert abs(achieved_conditional_mean - mean_copies_if_positive) < 0.05, (
                f"{marker} ({tissue}): fitted E[n|n>0]={achieved_conditional_mean:.4f} "
                f"does not match target {mean_copies_if_positive}"
            )


# ---------------------------------------------------------------------
# Test 2: stable_joint_population_yield() never returns a non-finite
# value, across the whole real-sigma plausibility range, for every
# marker/tissue population -- the direct regression test for the
# overflow bug this session found and fixed.
# ---------------------------------------------------------------------

def test_stable_yield_never_returns_nan_across_full_plausibility_range():
    k_intra_values = sorted(set(k_intra_bounds_from_real_sigma().values()))
    z = 1e-4

    pair_index = 0
    for marker, tissue_copies in COPIES_PER_POSITIVE_PARTICLE.items():
        for tissue in tissue_copies:
            n_vals, d_vals = reconstruct_marker_population(
                marker, tissue, n_samples=5000, seed=pair_index
            )
            pair_index += 1
            for k_intra in k_intra_values:
                y = stable_joint_population_yield(
                    n_vals, d_vals, SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, 1.0, k_intra, z,
                    patch_exponent=2.0,
                )
                assert np.isfinite(y), (
                    f"{marker} ({tissue}) at k_intra={k_intra:.3e}: yield is "
                    f"non-finite ({y}) -- the overflow-safe rewrite has "
                    f"regressed; see this test's module docstring for why "
                    f"this must never happen silently."
                )
                assert 0.0 <= y <= 1.0, (
                    f"{marker} ({tissue}) at k_intra={k_intra:.3e}: yield "
                    f"{y} is outside the physically valid [0, 1] range."
                )


def test_stable_yield_matches_original_where_original_does_not_overflow():
    """
    Cross-check: at K_INTRA values SMALL enough that the original
    population_model.joint_population_yield() does not overflow,
    stable_joint_population_yield() must give the SAME answer (to
    floating-point precision) -- confirming the overflow-safe rewrite is
    a genuine drop-in replacement, not a different model, everywhere the
    original is trustworthy.
    """
    from population_model import joint_population_yield

    n_vals, d_vals = reconstruct_marker_population("CD73", "AD", n_samples=2000, seed=0)
    # CD73's reconstruction has the smallest mean_n of the three markers,
    # so its smallest-and-mid-range K_INTRA values are the safest place to
    # check exact agreement without touching the overflow regime at all.
    small_k_intra_values = sorted(set(k_intra_bounds_from_real_sigma().values()))[:5]
    z = 1e-4

    for k_intra in small_k_intra_values:
        y_original = joint_population_yield(n_vals, d_vals, SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
                                              1.0, k_intra, z, patch_exponent=2.0)
        y_stable = stable_joint_population_yield(n_vals, d_vals, SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
                                                    1.0, k_intra, z, patch_exponent=2.0)
        assert np.isfinite(y_original), (
            "Test setup assumption violated: the original function was "
            "expected NOT to overflow at this k_intra for CD73(AD) -- if "
            "this fails, the k_intra grid or sample has changed and this "
            "test's small_k_intra_values slice needs revisiting."
        )
        assert abs(y_original - y_stable) < 1e-9, (
            f"k_intra={k_intra:.3e}: stable_joint_population_yield "
            f"({y_stable}) disagrees with the original, non-overflowing "
            f"joint_population_yield ({y_original}) -- the rewrite should "
            f"be numerically identical wherever the original is valid."
        )


# ---------------------------------------------------------------------
# Test 3: the marker-level ranking crossover (CD81 overtakes CD63 at the
# largest plausibility-range K_INTRA) reproduces deterministically.
# ---------------------------------------------------------------------

def test_cd81_overtakes_cd63_at_largest_k_intra_despite_lower_naive_score():
    naive = naive_fold_change_ranking()
    naive_marker_ranking = marker_level_ranking(naive["ranking_best_first"])
    assert naive_marker_ranking == ["CD63", "CD81", "CD73"], (
        f"Naive marker ranking changed from the expected ['CD63', 'CD81', "
        f"'CD73'] to {naive_marker_ranking} -- check "
        f"POSITIVITY_FRACTION_RANGE / COPIES_PER_POSITIVE_PARTICLE for an "
        f"unintended edit before trusting anything else in this test file."
    )

    result = sweep_marker_ranking_across_plausibility_range()
    marker_rankings = [marker_level_ranking(r) for r in result["nonlinear_rankings"]]

    largest_k_intra_ranking = marker_rankings[int(np.argmax(result["k_intra_values"]))]
    assert largest_k_intra_ranking[0] == "CD81", (
        f"Expected CD81 to lead the marker-level ranking at the largest "
        f"plausibility-range K_INTRA (the corrected, four-seed-verified "
        f"finding -- see handoff), but got ranking "
        f"{largest_k_intra_ranking} instead. If this is a genuine, "
        f"re-verified change (e.g. from a deliberate change to the "
        f"reconstructed copy numbers or diameter assumption), update this "
        f"test AND the handoff together -- do not just loosen this "
        f"assertion."
    )

    # And confirm the naive ranking is NOT what the nonlinear model gives
    # at this same K_INTRA -- i.e., this is a genuine crossover, not a
    # tautology.
    assert largest_k_intra_ranking != naive_marker_ranking, (
        "Expected a genuine ranking crossover at the largest K_INTRA "
        "(nonlinear ranking should differ from naive here) but the "
        "nonlinear ranking matched naive -- the crossover finding has "
        "not reproduced."
    )


def test_marker_level_crossover_is_robust_across_random_seeds():
    """
    The crossover found in test_cd81_overtakes_cd63_at_largest_k_intra_
    despite_lower_naive_score used the file's default seed=0
    reconstruction. Because this is a RECONSTRUCTED (not real digitized)
    population, a single-seed result is more vulnerable to being a
    sampling fluke than a result on real data -- this was checked
    directly across four different seed offsets before being reported as
    a finding (see the handoff). This test locks in that same
    seed-robustness check so it does not need to be re-verified by hand
    in every future session.
    """
    from k_intra_z_plausibility_bounds import k_intra_bounds_from_real_sigma

    k_intra_max = max(k_intra_bounds_from_real_sigma().values())
    z = 1e-4

    for seed_offset in (0, 100, 500, 1000):
        yields = {}
        pair_index = 0
        for marker, tissue_copies in COPIES_PER_POSITIVE_PARTICLE.items():
            for tissue in tissue_copies:
                n_vals, d_vals = reconstruct_marker_population(
                    marker, tissue, n_samples=5000, seed=seed_offset + pair_index
                )
                pair_index += 1
                yields[(marker, tissue)] = stable_joint_population_yield(
                    n_vals, d_vals, SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, 1.0, k_intra_max, z,
                    patch_exponent=2.0,
                )
        ranking = sorted(yields, key=yields.get, reverse=True)
        marker_ranking = marker_level_ranking(ranking)
        assert marker_ranking[0] == "CD81", (
            f"seed_offset={seed_offset}: expected CD81 to lead at the "
            f"largest K_INTRA (robust across seeds per the handoff), got "
            f"{marker_ranking} instead."
        )


# ---------------------------------------------------------------------
# Test 4: CD73 stays last in the marker-level ranking at EVERY tested
# K_INTRA -- the "always dominated" half of the finding.
# ---------------------------------------------------------------------

def test_cd73_stays_last_at_every_k_intra():
    result = sweep_marker_ranking_across_plausibility_range()
    for k_intra, ranking in zip(result["k_intra_values"], result["nonlinear_rankings"]):
        marker_ranking = marker_level_ranking(ranking)
        assert marker_ranking[-1] == "CD73", (
            f"k_intra={k_intra:.3e}: expected CD73 to remain last in the "
            f"marker-level ranking (its naive score is ~4-5x lower than "
            f"either CD81 or CD63 at every tissue), got {marker_ranking} "
            f"instead -- if this is a genuine finding, describe the "
            f"nonlinearity closing that large a gap explicitly in the "
            f"handoff rather than silently updating this test."
        )


if __name__ == "__main__":
    # Allow running without pytest, matching this project's other test
    # files' convention.
    test_functions = [
        test_negative_binomial_fit_hits_both_targets_for_every_marker_tissue,
        test_stable_yield_never_returns_nan_across_full_plausibility_range,
        test_stable_yield_matches_original_where_original_does_not_overflow,
        test_cd81_overtakes_cd63_at_largest_k_intra_despite_lower_naive_score,
        test_marker_level_crossover_is_robust_across_random_seeds,
        test_cd73_stays_last_at_every_k_intra,
    ]
    for test_fn in test_functions:
        try:
            test_fn()
            print(f"PASS: {test_fn.__name__}")
        except AssertionError as e:
            print(f"FAIL: {test_fn.__name__}\n  {e}")
