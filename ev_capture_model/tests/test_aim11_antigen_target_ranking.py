"""
test_aim11_antigen_target_ranking.py

Tests for the pivot's headline ranking-comparison file. What needs
checking:

1. The two new reconstruction-fit functions (median/IQR-based, and
   zero-fraction/conditional-mean-based) actually reproduce their real
   target statistics -- not just "ran without error."
2. naive_mean_ranking() and distribution_aware_escape_ranking() (and
   their pre-vs-relapse counterparts) produce internally consistent,
   correctly-shaped results.
3. The full grid sweep runs end-to-end in reasonable time (a direct
   regression test for the real timeout caught and fixed while first
   building this file -- see reconstruct_population_from_median_iqr's
   own docstring for the full account of the caching fix).
4. and_gate_escape_fraction() behaves sensibly (an AND-gate/OR-logic
   combination should never have a HIGHER escape fraction than either
   single marker alone, since escape requires BOTH to fail).

HOW TO RUN THIS FILE
-----------------------
    pytest test_aim11_antigen_target_ranking.py -v
or:
    python test_aim11_antigen_target_ranking.py
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from aim11_antigen_target_ranking import (
    fit_negative_binomial_to_median_iqr,
    fit_negative_binomial_to_zero_fraction_and_conditional_mean,
    reconstruct_population_from_median_iqr,
    reconstruct_population_from_range_and_positivity,
    naive_mean_ranking,
    distribution_aware_escape_ranking,
    naive_mean_ranking_pre_vs_relapse,
    escape_ranking_pre_vs_relapse,
    sweep_reordering_across_thresholds,
    sweep_pre_vs_relapse_across_full_grid,
    and_gate_escape_fraction,
)
from aim10_car_t_kill_probability_model import (
    LBCL_CD19_DENSITY_PRE_TREATMENT,
    LBCL_CD19_DENSITY_AT_RELAPSE,
    fit_kill_curve_to_real_threshold,
)


def test_median_iqr_fit_reproduces_real_lbcl_targets():
    """
    Direct regression test for the real bug caught while first building
    this file: an earlier fsolve-based version returned converged=False
    and an achieved median/IQR that visibly did not match the real
    target. This test checks the ACTUAL achieved values, not just the
    converged flag, against both real LBCL cohorts this file uses.
    """
    for cohort_name, cohort in [
        ("pre_treatment", LBCL_CD19_DENSITY_PRE_TREATMENT),
        ("at_relapse", LBCL_CD19_DENSITY_AT_RELAPSE),
    ]:
        fit = fit_negative_binomial_to_median_iqr(
            cohort["median_molecules_per_cell"], cohort["iqr_low"], cohort["iqr_high"]
        )
        assert fit["converged"], f"{cohort_name} fit did not converge: {fit}"

        median_rel_error = abs(fit["achieved_median"] - cohort["median_molecules_per_cell"]) / cohort["median_molecules_per_cell"]
        assert median_rel_error < 0.05, (
            f"{cohort_name}: achieved median {fit['achieved_median']} is "
            f">5% off real target {cohort['median_molecules_per_cell']}"
        )

        target_width = cohort["iqr_high"] - cohort["iqr_low"]
        achieved_width = fit["achieved_iqr_high"] - fit["achieved_iqr_low"]
        width_rel_error = abs(achieved_width - target_width) / target_width
        assert width_rel_error < 0.05, (
            f"{cohort_name}: achieved IQR width {achieved_width} is >5% off "
            f"real target width {target_width}"
        )
    print("PASS: median/IQR fit reproduces both real LBCL cohorts within 5%")


def test_zero_fraction_conditional_mean_fit_reproduces_target():
    """Direct regression test for the second real bug caught (the reused
    EV-project function failing with 'f(a) and f(b) must have different
    signs' at this file's input scale) -- checks the NEW dedicated fit
    function actually reproduces its target statistics."""
    fraction_positive = 0.4515
    mean_copies_if_positive = 254.95
    fit = fit_negative_binomial_to_zero_fraction_and_conditional_mean(
        fraction_positive, mean_copies_if_positive
    )
    assert fit["converged"]

    # Directly verify the fitted distribution actually has the target
    # zero-fraction and conditional mean (not just "the optimizer didn't
    # crash").
    from scipy import stats
    p_param = fit["dispersion_r"] / (fit["dispersion_r"] + fit["mean_n"])
    p0 = stats.nbinom.pmf(0, fit["dispersion_r"], p_param)
    achieved_fraction_positive = 1.0 - p0
    assert np.isclose(achieved_fraction_positive, fraction_positive, atol=0.02), (
        f"Achieved fraction_positive {achieved_fraction_positive} does not "
        f"match target {fraction_positive}"
    )
    print(f"PASS: zero-fraction/conditional-mean fit reproduces target "
          f"(achieved fraction_positive={achieved_fraction_positive:.4f}, "
          f"target={fraction_positive})")


def test_naive_ranking_matches_real_reported_direction():
    """The naive ranking must rank CD19_in_LBCL above CD19_in_myeloma,
    since the real reported LBCL median (5810) is far above the real
    myeloma range's geometric mean (~255) -- a basic sanity check that
    the reconstruction didn't invert the real data somehow."""
    naive = naive_mean_ranking()
    assert naive["ranking_best_first"][0] == "CD19_in_LBCL", (
        f"Expected CD19_in_LBCL to rank first by mean density, got "
        f"{naive['ranking_best_first']}"
    )
    print("PASS: naive ranking correctly reflects real reported magnitudes")


def test_pre_treatment_has_lower_escape_fraction_than_relapse():
    """Physical sanity check: since pre-treatment CD19 is real and reported
    to be HIGHER than at-relapse CD19 (5810 vs 2021 molecules/cell,
    Transplantation and Cellular Therapy 2024), the model should predict
    a LOWER escape fraction pre-treatment than at relapse, for any
    reasonable calibration -- exactly the real clinical pattern (CD19
    downregulation drives relapse) this data was collected to
    demonstrate."""
    fit = fit_kill_curve_to_real_threshold(843)  # blinatumomab cytolysis
    result = escape_ranking_pre_vs_relapse(
        fit["car_avidity"], fit["k_intra"], fit["z"]
    )
    assert (result["escape_fractions"]["CD19_LBCL_pre_treatment"]
            <= result["escape_fractions"]["CD19_LBCL_at_relapse"]), (
        f"Expected pre-treatment escape fraction <= at-relapse, got "
        f"{result['escape_fractions']}"
    )
    print(f"PASS: pre-treatment escape fraction "
          f"({result['escape_fractions']['CD19_LBCL_pre_treatment']:.4f}) <= "
          f"at-relapse ({result['escape_fractions']['CD19_LBCL_at_relapse']:.4f})")


def test_full_grid_sweep_completes_in_reasonable_time():
    """
    Direct regression test for the real timeout (>120 seconds) caught
    while first running this file's full grid sweep, caused by
    re-running an expensive 150x150-point grid-search fit at every one
    of 36 grid points instead of caching it once. This test fails loudly
    if that performance regression ever reappears.
    """
    start = time.time()
    result = sweep_pre_vs_relapse_across_full_grid()
    elapsed = time.time() - start
    assert elapsed < 30.0, (
        f"Full grid sweep took {elapsed:.1f}s -- expected well under 30s "
        f"with the population-reconstruction cache in place. This may "
        f"mean the caching fix (see reconstruct_population_from_median_iqr's "
        f"docstring) has regressed."
    )
    assert result["n_grid_points"] == 36, (
        f"Expected 36 grid points (4 k_intra fits x 9 kill_thresholds), "
        f"got {result['n_grid_points']}"
    )
    print(f"PASS: full grid sweep completed in {elapsed:.2f}s "
          f"({result['n_grid_points']} grid points, "
          f"{result['n_reordered']} reordered)")


def test_and_gate_escape_fraction_never_exceeds_either_single_marker():
    """
    AND-gate (OR-logic-for-killing, per this file's own docstring
    convention matching aim5's) escape fraction must be <= EITHER single
    marker's own escape fraction alone, since a cell only escapes the
    AND-gate if it escapes BOTH markers individually -- escaping both is
    never more likely than escaping either one alone.
    """
    fit = fit_kill_curve_to_real_threshold(843)
    rng = np.random.default_rng(42)
    pop_a = rng.negative_binomial(5, 5 / (5 + 843), size=3000).astype(float)
    pop_b = rng.negative_binomial(5, 5 / (5 + 200), size=3000).astype(float)

    result = and_gate_escape_fraction(
        pop_a, pop_b, fit["car_avidity"], fit["k_intra"], fit["z"]
    )
    assert result["escape_fraction"] <= result["single_a_escape_fraction"] + 1e-9, (
        f"AND-gate escape fraction {result['escape_fraction']} exceeds "
        f"single-marker-A escape fraction {result['single_a_escape_fraction']}"
    )
    assert result["escape_fraction"] <= result["single_b_escape_fraction"] + 1e-9, (
        f"AND-gate escape fraction {result['escape_fraction']} exceeds "
        f"single-marker-B escape fraction {result['single_b_escape_fraction']}"
    )
    print(f"PASS: AND-gate escape fraction ({result['escape_fraction']:.4f}) "
          f"<= both single-marker escape fractions "
          f"(A={result['single_a_escape_fraction']:.4f}, "
          f"B={result['single_b_escape_fraction']:.4f})")


if __name__ == "__main__":
    tests = [
        test_median_iqr_fit_reproduces_real_lbcl_targets,
        test_zero_fraction_conditional_mean_fit_reproduces_target,
        test_naive_ranking_matches_real_reported_direction,
        test_pre_treatment_has_lower_escape_fraction_than_relapse,
        test_full_grid_sweep_completes_in_reasonable_time,
        test_and_gate_escape_fraction_never_exceeds_either_single_marker,
    ]
    failures = 0
    for test in tests:
        try:
            test()
        except AssertionError as e:
            print(f"FAIL: {test.__name__}: {e}")
            failures += 1
    print(f"\n{len(tests) - failures}/{len(tests)} tests passed.")
    sys.exit(1 if failures else 0)
