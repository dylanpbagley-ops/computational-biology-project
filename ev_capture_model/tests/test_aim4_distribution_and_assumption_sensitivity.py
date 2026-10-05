"""
test_aim4_distribution_and_assumption_sensitivity.py

Regression tests for src/aim4_distribution_and_assumption_sensitivity.py --
the CURRENT.md handoff's "Immediate next steps" item 1: distribution-family
sensitivity and a fuller plausibility-bounds assumption sensitivity, beyond
the closed target_exponents check.

WHAT THESE TESTS LOCK IN, AND WHY EACH ONE MATTERS
-------------------------------------------------------
1. fit_zero_inflated_poisson_to_positivity() actually hits both of its two
   targets (empirically, via a large sample, not just algebraically) --
   the same discipline test_stem_cell_marker_ranking_extension.py already
   applies to the negative-binomial fit, applied here to its ZIP
   counterpart so a future refactor of the ZIP construction cannot
   silently drift away from matching the same two summary statistics.
2. The stem-cell marker-ranking crossover (CD81 overtakes CD63 at the
   largest plausibility-range K_INTRA) survives being reconstructed with
   a zero-inflated Poisson instead of a negative binomial -- the actual
   new finding this file exists to check, locked in as a real assertion
   so a future change to either reconstruction cannot silently reverse it
   without a test failing.
3. The superselectivity-reachability finding (any_superselective=True)
   is stable across three different n_values grid choices for the
   plausibility bounds, and across the full documented sigma range
   (40% coverage operating point through 100% close-packed ceiling).
4. The same crossover survives four different diameter-model variants
   (two alternate CVs and a genuinely different family, gamma).

Each of these was actually run and its printed output inspected before
being written as an assertion here -- not assumed to hold from the
module's own docstring claims.
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from aim4_distribution_and_assumption_sensitivity import (
    fit_zero_inflated_poisson_to_positivity,
    sample_zero_inflated_poisson,
    reconstruct_marker_population_zip,
    check_crossover_survives_zip_family,
    sweep_superselectivity_across_n_values_grids,
    sweep_superselectivity_across_sigma_choices,
    check_crossover_survives_diameter_assumption,
)
from stem_cell_marker_ranking_extension import (
    POSITIVITY_FRACTION_RANGE,
    COPIES_PER_POSITIVE_PARTICLE,
    fit_negative_binomial_to_positivity,
)


# ---------------------------------------------------------------------
# Test 1: the ZIP fit hits both targets empirically, for every
# marker/tissue combination actually used in this project.
# ---------------------------------------------------------------------

def test_zip_fit_hits_both_targets_for_every_marker_tissue():
    for marker, tissue_copies in COPIES_PER_POSITIVE_PARTICLE.items():
        lo, hi = POSITIVITY_FRACTION_RANGE[marker]
        fraction_positive = 0.5 * (lo + hi)
        for tissue, mean_copies_if_positive in tissue_copies.items():
            samp = sample_zero_inflated_poisson(
                fraction_positive, mean_copies_if_positive, n_samples=200_000, seed=1
            )
            emp_frac_positive = float((samp > 0).mean())
            emp_mean_if_positive = float(samp[samp > 0].mean())

            assert abs(emp_frac_positive - fraction_positive) < 0.01, (
                f"{marker}/{tissue}: ZIP empirical fraction_positive "
                f"{emp_frac_positive:.4f} does not match target "
                f"{fraction_positive:.4f} within tolerance."
            )
            assert abs(emp_mean_if_positive - mean_copies_if_positive) < 0.5, (
                f"{marker}/{tissue}: ZIP empirical mean_copies_if_positive "
                f"{emp_mean_if_positive:.2f} does not match target "
                f"{mean_copies_if_positive} within tolerance."
            )


def test_zip_and_negative_binomial_have_matching_overall_mean():
    """
    Both reconstructions are matched to the SAME two summary statistics
    (fraction_positive, mean_copies_if_positive), so their resulting
    OVERALL mean_n should be very close -- this is the property that
    makes the family comparison fair (both families represent
    approximately the same "naive" mean, only the SHAPE around that mean
    differs).

    CHECKED DIRECTLY, NOT ASSUMED EXACT: the ZIP construction is exact by
    algebraic construction (mean_n = fraction_positive * mean_copies_if_
    positive, verified in test_zip_fit_hits_both_targets_for_every_
    marker_tissue). The negative-binomial fit, by contrast, finds
    dispersion_r via a 1D root-find (fit_negative_binomial_to_positivity's
    own docstring) that exactly hits P(n=0)=1-fraction_positive (confirmed
    directly: scipy.stats.nbinom.pmf(0, r, p) matches to ~1e-9), but its
    resulting mean_n is only approximately fraction_positive *
    mean_copies_if_positive -- checked directly for CD63/AD: nb mean_n =
    26.5021 vs. the exact target product 26.5, a real ~0.002 (0.008%)
    discrepancy inherent to how fit_negative_binomial_to_positivity solves
    its nested root-find (it matches the zero-fraction exactly, then
    separately solves for the mean_n that hits the CONDITIONAL mean given
    that dispersion, which does not algebraically guarantee the overall
    mean lands exactly on the naive product). This is a real property of
    the existing negative-binomial fit, not a bug introduced by this
    file -- the tolerance below (0.05, comfortably above the ~0.008%
    relative discrepancy actually observed) reflects that, rather than
    hiding a larger, more concerning mismatch.
    """
    for marker, tissue_copies in COPIES_PER_POSITIVE_PARTICLE.items():
        lo, hi = POSITIVITY_FRACTION_RANGE[marker]
        fraction_positive = 0.5 * (lo + hi)
        for tissue, mean_copies_if_positive in tissue_copies.items():
            zip_fit = fit_zero_inflated_poisson_to_positivity(fraction_positive, mean_copies_if_positive)
            nb_fit = fit_negative_binomial_to_positivity(fraction_positive, mean_copies_if_positive)
            assert abs(zip_fit["mean_n"] - nb_fit["mean_n"]) < 0.05, (
                f"{marker}/{tissue}: ZIP mean_n ({zip_fit['mean_n']}) and "
                f"negative-binomial mean_n ({nb_fit['mean_n']}) differ by "
                f"more than the expected small root-find discrepancy "
                f"(~0.002 typical) -- investigate before trusting the "
                f"family comparison as apples-to-apples."
            )


# ---------------------------------------------------------------------
# Test 2: the stem-cell crossover survives the zero-inflated-Poisson
# distribution-family alternative.
# ---------------------------------------------------------------------

def test_crossover_survives_zip_family():
    result = check_crossover_survives_zip_family()
    assert result["naive_marker_ranking"] == ["CD63", "CD81", "CD73"], (
        f"Naive marker ranking changed from the expected baseline -- "
        f"got {result['naive_marker_ranking']}, check "
        f"POSITIVITY_FRACTION_RANGE / COPIES_PER_POSITIVE_PARTICLE for an "
        f"unintended edit before trusting this test."
    )
    assert all(result["cd81_leads_under_nb"]), (
        f"CD81 did not lead under negative binomial at every seed "
        f"({result['nb_marker_rankings']}) -- the baseline finding itself "
        f"did not reproduce, so the family-robustness check cannot be "
        f"meaningfully evaluated."
    )
    assert all(result["cd81_leads_under_zip"]), (
        f"CD81 did NOT lead under the zero-inflated-Poisson "
        f"reconstruction at every seed ({result['zip_marker_rankings']}) "
        f"-- the crossover finding does not survive this distribution-"
        f"family alternative. This is a genuine result to report "
        f"honestly, not a bug to silently fix by loosening this "
        f"assertion."
    )
    assert result["crossover_survives_zip_family"] is True


# ---------------------------------------------------------------------
# Test 3: superselectivity-reachability is stable across n_values grid
# choices and across the documented sigma range.
# ---------------------------------------------------------------------

def test_superselectivity_stable_across_n_values_grids():
    result = sweep_superselectivity_across_n_values_grids()
    for grid, r in result.items():
        if grid == "qualitative_finding_stable":
            continue
        assert r["any_superselective"] is True, (
            f"n_values grid {grid}: superselectivity (alpha>1) was NOT "
            f"reached anywhere in the plausibility grid under this "
            f"alternative receptor-count assumption -- this would be a "
            f"genuine, newly-discovered fragility in the headline "
            f"superselectivity finding, worth reporting directly rather "
            f"than adjusting the grid choice to hide it."
        )
    assert result["qualitative_finding_stable"] is True


def test_superselectivity_stable_across_sigma_choices():
    result = sweep_superselectivity_across_sigma_choices()
    for label, r in result.items():
        if label == "qualitative_finding_stable":
            continue
        assert r["any_superselective"] is True, (
            f"sigma choice {label!r}: superselectivity (alpha>1) was NOT "
            f"reached anywhere in the plausibility grid at this sigma -- "
            f"a genuine finding worth reporting directly."
        )
    assert result["qualitative_finding_stable"] is True

    # Sanity check on the specific numbers, not just the boolean: the
    # 40%-coverage and 100%-max-packed sigma choices should give
    # IDENTICAL max_alpha, because k_intra_bounds_from_real_sigma()
    # rescales K_INTRA inversely with sigma, and q_bound only depends on
    # the product sigma*K_INTRA (the identifiability property documented
    # in Limitations_and_Calibration_Gap.md) -- so changing sigma alone,
    # with K_INTRA re-derived consistently against it, must reproduce
    # the exact same swept product grid.
    forty_pct = result["lopez_baltazar_40pct"]["max_alpha_found"]
    hundred_pct = result["lopez_baltazar_100pct_max"]["max_alpha_found"]
    assert abs(forty_pct - hundred_pct) < 1e-6, (
        f"Expected max_alpha to be IDENTICAL across sigma choices (since "
        f"K_INTRA is re-derived consistently against each sigma, so the "
        f"product sigma*K_INTRA sweep is unchanged -- see the "
        f"identifiability proof in Limitations_and_Calibration_Gap.md), "
        f"but got {forty_pct} vs {hundred_pct}."
    )


# ---------------------------------------------------------------------
# Test 4: the stem-cell crossover survives alternative diameter-model
# assumptions (wider/narrower lognormal CV, and a different family).
# ---------------------------------------------------------------------

def test_crossover_survives_diameter_assumption_variants():
    result = check_crossover_survives_diameter_assumption()
    for variant, rankings in result.items():
        if variant == "crossover_survives_all_diameter_variants":
            continue
        assert all(r[0] == "CD81" for r in rankings), (
            f"Diameter variant {variant}: CD81 did not lead at every "
            f"seed ({rankings}) -- the crossover finding does not "
            f"survive this diameter assumption. Report this honestly; "
            f"do not narrow the tested variants to avoid it."
        )
    assert result["crossover_survives_all_diameter_variants"] is True


if __name__ == "__main__":
    # Allow running without pytest, matching this project's other test
    # files' convention.
    test_functions = [
        test_zip_fit_hits_both_targets_for_every_marker_tissue,
        test_zip_and_negative_binomial_have_matching_overall_mean,
        test_crossover_survives_zip_family,
        test_superselectivity_stable_across_n_values_grids,
        test_superselectivity_stable_across_sigma_choices,
        test_crossover_survives_diameter_assumption_variants,
    ]
    for test_fn in test_functions:
        try:
            test_fn()
            print(f"PASS: {test_fn.__name__}")
        except AssertionError as e:
            print(f"FAIL: {test_fn.__name__}\n  {e}")
