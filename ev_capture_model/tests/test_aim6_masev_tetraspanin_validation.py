"""
test_aim6_masev_tetraspanin_validation.py

Regression tests for src/aim6_masev_tetraspanin_validation.py -- the
project's first real external-data validation (see that file's own module
docstring for the full account of why this file exists and what it does
and does not claim).

WHAT THESE TESTS LOCK IN, AND WHY EACH ONE MATTERS
-------------------------------------------------------
1. The real data itself is internally consistent: each cell line's 8
   categories sum to 1.0 (a basic sanity check on the transcription from
   the paper's figure -- if this fails, a number was mistyped).
2. marginals_from_real_data() correctly derives each marker's marginal
   positivity fraction by summing the right categories -- checked against
   an independently hand-computed value, not just re-deriving the same
   arithmetic the function itself does.
3. The marginal negative-binomial fit actually reproduces the real
   positivity fraction it was fit to (empirically, by simulating from it),
   for every marker/cell-line combination actually used in this file --
   the same kind of "does the fit actually hit its target" check this
   project's other reconstruction fits (aim5, aim3b,
   stem_cell_marker_ranking_extension) all include for themselves.
4. The 3-marker copula preserves each marginal regardless of the shared
   correlation parameter (the direct 3-marker analog of aim5's own
   test_copula_reproduces_correct_marginal_mean_regardless_of_correlation).
5. correlation=0 reproduces pairwise independence to a tight empirical
   tolerance (the 3-marker analog of aim5's
   test_correlation_zero_reproduces_independence).
6. categorize_joint_population()'s 8 category fractions always sum to
   1.0 (every vesicle falls into exactly one category, by construction --
   a basic correctness check on the categorization logic itself).
7. total_variation_distance() is 0 for identical distributions and
   correctly bounded in [0, 1] -- a basic correctness check on the metric,
   independent of any real data.
8. THE HEADLINE FINDING: validate_cell_line() finds at least one swept
   correlation, for BOTH real cell lines, whose predicted joint
   distribution is within a generously-set TVD tolerance of the real
   data -- i.e., this project's existing joint-reconstruction machinery,
   fed the real marginals, is actually CAPABLE of predicting something
   close to the real published joint distribution shape, not just
   producing some number that happens to be reported. The tolerance
   (0.05) is set deliberately loose relative to the actual best-fit TVDs
   found (0.011 for PANC-1, 0.019 for CAPAN-2, both well under this bound)
   so that this test is a real regression guard against the fit
   meaningfully breaking, not a tight bound that would be brittle to
   incidental numerical changes.
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from aim6_masev_tetraspanin_validation import (
    REAL_MASEV_DATA,
    CATEGORY_ORDER,
    CATEGORY_MARKERS_POSITIVE,
    ASPC1_A549_REAL_MARGINALS,
    marginals_from_real_data,
    fit_marginal_to_positivity_fraction,
    fit_all_marginals_for_cell_line,
    simulate_joint_tetraspanin_population,
    categorize_joint_population,
    total_variation_distance,
    validate_cell_line,
    independence_baseline_tvd,
    fit_committed_correlation,
    held_out_prediction_for_cell_line,
)


# ---------------------------------------------------------------------
# Test 1: real data sums to 1.0 per cell line.
# ---------------------------------------------------------------------

def test_real_data_categories_sum_to_one():
    """
    Tolerance is 0.2 percentage points, not exact -- checked directly
    against the transcribed numbers (both cell lines re-verified against
    the buildlog this test suite was written from): CAPAN-2 sums to
    exactly 100.0%, but PANC-1 sums to 99.9% (24.8+37.6+4.6+4.4+11.2+1.6+
    10.6+5.1 = 99.9), a genuine rounding artifact in the published
    figure's own one-decimal-place percentages, not a transcription typo
    here. A tight 1e-6 tolerance would incorrectly flag a real feature of
    the source data as a bug.
    """
    for cell_line, data in REAL_MASEV_DATA.items():
        total = sum(data.values())
        assert abs(total - 1.0) < 2e-3, (
            f"{cell_line}: real Figure 4 categories sum to {total:.4f}, "
            f"expected ~1.0 (within rounding) -- check the transcribed "
            f"numbers for a typo if this is off by more than a rounding "
            f"artifact."
        )


# ---------------------------------------------------------------------
# Test 2: marginals_from_real_data() matches an independently
# hand-computed value.
# ---------------------------------------------------------------------

def test_marginals_from_real_data_match_hand_computation():
    # PANC-1, hand-computed directly from REAL_MASEV_DATA's own numbers:
    # P(CD9+) = CD9_only + CD81_CD9 + CD63_CD9 + all_three
    #         = 0.376 + 0.112 + 0.106 + 0.051 = 0.645
    # P(CD63+) = CD63_only + CD63_CD81 + CD63_CD9 + all_three
    #          = 0.044 + 0.016 + 0.106 + 0.051 = 0.217
    # P(CD81+) = CD81_only + CD81_CD9 + CD63_CD81 + all_three
    #          = 0.046 + 0.112 + 0.016 + 0.051 = 0.225
    marg = marginals_from_real_data("PANC-1")
    assert abs(marg["CD9"] - 0.645) < 1e-9, marg
    assert abs(marg["CD63"] - 0.217) < 1e-9, marg
    assert abs(marg["CD81"] - 0.225) < 1e-9, marg


# ---------------------------------------------------------------------
# Test 3: the marginal negative-binomial fit reproduces its target
# positivity fraction, empirically.
# ---------------------------------------------------------------------

def test_marginal_fit_reproduces_target_positivity_fraction():
    from scipy import stats as scipy_stats
    for cell_line in REAL_MASEV_DATA:
        marg = marginals_from_real_data(cell_line)
        for marker, target_frac in marg.items():
            fit = fit_marginal_to_positivity_fraction(target_frac)
            p_param = fit["dispersion_r"] / (fit["dispersion_r"] + fit["mean_n"])
            p0 = scipy_stats.nbinom.pmf(0, fit["dispersion_r"], p_param)
            achieved_frac = 1.0 - p0
            assert abs(achieved_frac - target_frac) < 1e-4, (
                f"{cell_line}/{marker}: fit achieves positivity fraction "
                f"{achieved_frac:.6f}, target was {target_frac:.6f}."
            )


# ---------------------------------------------------------------------
# Test 4 & 5: the 3-marker copula preserves marginals and reproduces
# independence at correlation=0.
# ---------------------------------------------------------------------

def test_copula_preserves_marginal_positivity_fraction_regardless_of_correlation():
    marginal_fits = fit_all_marginals_for_cell_line("PANC-1")
    real_marg = marginals_from_real_data("PANC-1")
    for correlation in (0.0, 0.4, 0.8):
        counts = simulate_joint_tetraspanin_population(
            marginal_fits, correlation, n_vesicles=200_000, seed=0
        )
        for marker in ("CD9", "CD63", "CD81"):
            empirical_frac = float((counts[marker] > 0).mean())
            target_frac = real_marg[marker]
            assert abs(empirical_frac - target_frac) < 0.02, (
                f"correlation={correlation}, marker={marker}: empirical "
                f"positivity fraction {empirical_frac:.4f} does not match "
                f"target {target_frac:.4f} within 0.02 -- the copula is not "
                f"preserving its marginal."
            )


def test_correlation_zero_reproduces_pairwise_independence():
    from scipy import stats as scipy_stats
    marginal_fits = fit_all_marginals_for_cell_line("PANC-1")
    counts = simulate_joint_tetraspanin_population(marginal_fits, 0.0, n_vesicles=200_000, seed=0)
    pairs = [("CD9", "CD63"), ("CD9", "CD81"), ("CD63", "CD81")]
    for a, b in pairs:
        r, _ = scipy_stats.pearsonr(counts[a], counts[b])
        assert abs(r) < 0.02, (
            f"Expected correlation=0.0 to reproduce pairwise independence "
            f"between {a} and {b}, got Pearson r={r:.4f}."
        )


# ---------------------------------------------------------------------
# Test 6: categorize_joint_population() sums to 1.0.
# ---------------------------------------------------------------------

def test_categorize_joint_population_sums_to_one():
    marginal_fits = fit_all_marginals_for_cell_line("CAPAN-2")
    counts = simulate_joint_tetraspanin_population(marginal_fits, 0.3, n_vesicles=50_000, seed=1)
    predicted = categorize_joint_population(counts)
    total = sum(predicted.values())
    assert abs(total - 1.0) < 1e-9, (
        f"Categorized fractions sum to {total}, expected 1.0 -- every "
        f"vesicle should fall into exactly one of the 8 categories."
    )


# ---------------------------------------------------------------------
# Test 7: total_variation_distance() correctness.
# ---------------------------------------------------------------------

def test_tvd_zero_for_identical_distributions():
    d = {cat: 1.0 / len(CATEGORY_ORDER) for cat in CATEGORY_ORDER}
    assert total_variation_distance(d, d) == 0.0


def test_tvd_bounded_in_zero_one():
    for cell_line in REAL_MASEV_DATA:
        result = validate_cell_line(cell_line, correlations=(0.0, 0.8), n_vesicles=20_000, seed=0)
        for tvd in result["tvd"]:
            assert 0.0 <= tvd <= 1.0, f"TVD {tvd} out of bounds [0, 1]."


# ---------------------------------------------------------------------
# Test 8: THE HEADLINE FINDING -- the model's best-fit correlation gets
# genuinely close to the real data for both cell lines.
# ---------------------------------------------------------------------

def test_best_fit_tvd_is_small_for_both_real_cell_lines():
    for cell_line in REAL_MASEV_DATA:
        result = validate_cell_line(cell_line, n_vesicles=100_000, seed=0)
        assert result["best_tvd"] < 0.05, (
            f"{cell_line}: best-fit TVD is {result['best_tvd']:.4f}, "
            f"expected < 0.05 -- this is the headline validation finding; "
            f"if this now fails, the joint-reconstruction approach no "
            f"longer predicts a close match to the real MASEV data, which "
            f"is a real and important result to report honestly, not to "
            f"paper over by loosening this tolerance."
        )


def test_independence_baseline_is_worse_than_or_similar_to_best_fit():
    """
    The independence-only (correlation=0) prediction should never be
    BETTER than the best swept correlation, by construction (correlation=0
    is itself one of the swept values in CORRELATION_SWEEP, so the best
    TVD is always <= the independence TVD). This is a basic internal-
    consistency check, not a claim about which correlation is "true."
    """
    for cell_line in REAL_MASEV_DATA:
        best = validate_cell_line(cell_line, n_vesicles=100_000, seed=0)["best_tvd"]
        indep = independence_baseline_tvd(cell_line, n_vesicles=100_000, seed=0)
        assert best <= indep + 1e-9, (
            f"{cell_line}: best_tvd ({best:.4f}) should never exceed the "
            f"independence baseline ({indep:.4f}), since correlation=0.0 "
            f"is itself included in the swept range."
        )


# ---------------------------------------------------------------------
# Tests 9-13: the ASPC1/A549 held-out prediction -- fit_committed_
# correlation() and held_out_prediction_for_cell_line(). See that
# module's own docstring ("HELD-OUT PREDICTION DATA" section) for the
# full account of what this validates and what it explicitly does not.
# ---------------------------------------------------------------------

def test_committed_correlation_is_from_the_swept_grid():
    """
    fit_committed_correlation() must return one of the values actually
    tested in CORRELATION_SWEEP -- it selects among candidates, it does
    not invent a new one via continuous optimization. A value outside the
    tested grid would mean the selection logic is broken.
    """
    from aim6_masev_tetraspanin_validation import CORRELATION_SWEEP
    committed = fit_committed_correlation()
    assert committed["committed_correlation"] in CORRELATION_SWEEP, (
        f"Committed correlation {committed['committed_correlation']} is not "
        f"one of the swept candidates {CORRELATION_SWEEP}."
    )


def test_committed_correlation_minimizes_average_training_tvd():
    """
    Directly re-derive the average-TVD-minimizing correlation by brute
    force across the same grid, and check fit_committed_correlation()
    agrees -- an independent re-computation, not just re-reading the
    function's own internal logic back.
    """
    from aim6_masev_tetraspanin_validation import CORRELATION_SWEEP

    training = ("PANC-1", "CAPAN-2")
    tvd_by_corr = {corr: [] for corr in CORRELATION_SWEEP}
    for cell_line in training:
        result = validate_cell_line(cell_line, correlations=CORRELATION_SWEEP, n_vesicles=100_000, seed=0)
        for corr, tvd in zip(result["correlations"], result["tvd"]):
            tvd_by_corr[corr].append(tvd)
    average_tvd = {corr: np.mean(tvds) for corr, tvds in tvd_by_corr.items()}
    expected_best = min(average_tvd, key=average_tvd.get)

    committed = fit_committed_correlation(training_cell_lines=training, n_vesicles=100_000, seed=0)
    assert committed["committed_correlation"] == expected_best, (
        f"fit_committed_correlation() chose {committed['committed_correlation']}, "
        f"but an independent brute-force recomputation across the same grid "
        f"found {expected_best} minimizes average training TVD."
    )


def test_held_out_prediction_rejects_unknown_cell_line():
    committed = fit_committed_correlation()
    try:
        held_out_prediction_for_cell_line("PANC-1", committed["committed_correlation"])
        assert False, (
            "Expected held_out_prediction_for_cell_line() to reject 'PANC-1' "
            "(a TRAINING cell line with real joint data -- use "
            "validate_cell_line() for it instead), but it did not raise."
        )
    except (ValueError, KeyError):
        pass


def test_held_out_prediction_categories_sum_to_one():
    committed = fit_committed_correlation()
    for cell_line in ASPC1_A549_REAL_MARGINALS:
        result = held_out_prediction_for_cell_line(
            cell_line, committed["committed_correlation"], n_vesicles=50_000, seed=0
        )
        total = sum(result["predicted_categories"].values())
        assert abs(total - 1.0) < 1e-9, (
            f"{cell_line}: predicted categories sum to {total}, expected 1.0."
        )


def test_held_out_predicted_marginals_recovered_closely():
    """
    THE ONE GENUINELY CHECKABLE PIECE of the held-out prediction (see
    module docstring): the marginals recomputed from the predicted joint
    distribution should closely match the real input marginals. This is
    not circular -- the categorization step (thresholding simulated counts
    at >0 per marker) is a separate computation from the marginal fit
    itself, so a construction bug (e.g. a mismatched marker order when
    building the copula's covariance matrix, or an off-by-one in which
    column corresponds to which marker) COULD make these drift apart even
    though the marginal fit was individually correct. A tight tolerance
    (2 percentage points) is used because this is checking arithmetic
    correctness, not a scientific hypothesis -- a real bug here should
    fail loudly, not be tuned away with a loose tolerance.
    """
    committed = fit_committed_correlation()
    for cell_line in ASPC1_A549_REAL_MARGINALS:
        result = held_out_prediction_for_cell_line(
            cell_line, committed["committed_correlation"], n_vesicles=200_000, seed=0
        )
        for marker, error in result["marginal_recovery_error"].items():
            assert error < 0.02, (
                f"{cell_line}/{marker}: marginal recovery error {error:.4f} "
                f"(as a fraction) exceeds 0.02 -- the predicted joint "
                f"distribution's own implied marginal has drifted from the "
                f"real input marginal it was fit to. This points to a "
                f"construction bug in the copula or categorization step, "
                f"not a scientific finding to report as-is."
            )


def test_held_out_predicted_categories_match_marker_definitions():
    """
    Sanity check that the predicted categories for a held-out cell line
    are internally consistent with CATEGORY_MARKERS_POSITIVE's own
    definitions -- e.g., "CD9_only" should never have a HIGHER predicted
    fraction than the overall CD9+ marginal (since CD9_only is a SUBSET
    of CD9+ vesicles, specifically those also negative for CD63 and
    CD81). This is a basic set-theoretic property that must hold for any
    valid categorization, independent of whether the correlation or
    marginals happen to be scientifically correct.
    """
    committed = fit_committed_correlation()
    for cell_line in ASPC1_A549_REAL_MARGINALS:
        result = held_out_prediction_for_cell_line(
            cell_line, committed["committed_correlation"], n_vesicles=100_000, seed=0
        )
        predicted = result["predicted_categories"]
        recomputed_marginals = result["predicted_marginals_recomputed"]
        for marker in ("CD9", "CD63", "CD81"):
            only_category = f"{marker}_only"
            assert predicted[only_category] <= recomputed_marginals[marker] + 1e-9, (
                f"{cell_line}: predicted {only_category} fraction "
                f"({predicted[only_category]:.4f}) exceeds the overall "
                f"{marker}+ marginal ({recomputed_marginals[marker]:.4f}) -- "
                f"{only_category} must be a subset of {marker}+ vesicles."
            )


if __name__ == "__main__":
    # Allow running without pytest, matching this project's other test
    # files' convention.
    test_functions = [
        test_real_data_categories_sum_to_one,
        test_marginals_from_real_data_match_hand_computation,
        test_marginal_fit_reproduces_target_positivity_fraction,
        test_copula_preserves_marginal_positivity_fraction_regardless_of_correlation,
        test_correlation_zero_reproduces_pairwise_independence,
        test_categorize_joint_population_sums_to_one,
        test_tvd_zero_for_identical_distributions,
        test_tvd_bounded_in_zero_one,
        test_best_fit_tvd_is_small_for_both_real_cell_lines,
        test_independence_baseline_is_worse_than_or_similar_to_best_fit,
        test_committed_correlation_is_from_the_swept_grid,
        test_committed_correlation_minimizes_average_training_tvd,
        test_held_out_prediction_rejects_unknown_cell_line,
        test_held_out_prediction_categories_sum_to_one,
        test_held_out_predicted_marginals_recovered_closely,
        test_held_out_predicted_categories_match_marker_definitions,
    ]
    for test_fn in test_functions:
        try:
            test_fn()
            print(f"PASS: {test_fn.__name__}")
        except AssertionError as e:
            print(f"FAIL: {test_fn.__name__}\n  {e}")
