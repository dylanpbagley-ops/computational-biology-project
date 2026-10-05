"""
test_aim10_car_t_kill_probability_model.py

Tests for the pivot's core model. Three things need checking:

1. kill_probability() is genuinely identical to
   binding_model.capture_probability() (it should be -- it's a thin
   wrapper, not a reimplementation; a test that isn't identical would mean
   the wrapper introduced a real bug).
2. fit_kill_curve_to_real_threshold() actually reproduces the real
   threshold it was fit to (i.e. plugging the fitted k_intra back into
   kill_probability() at the real threshold density gives back
   threshold_kill_probability, not some other number) -- and does so
   WITHOUT raising the NaN/overflow error that was caught and fixed three
   times while first building this file (see that function's own
   docstring for the full account) for ALL three real thresholds, not
   just one.
3. escape_fraction() behaves physically sensibly: a population entirely
   above a calibrated threshold should show ~0% escape; a population
   entirely below should show ~100% escape; a population straddling the
   threshold should show something in between.

HOW TO RUN THIS FILE
-----------------------
From the command line, from inside the tests/ folder (with src/ next to
it):
    pytest test_aim10_car_t_kill_probability_model.py -v
or, if pytest isn't installed:
    python test_aim10_car_t_kill_probability_model.py
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from binding_model import capture_probability
from aim10_car_t_kill_probability_model import (
    kill_probability,
    escape_fraction,
    fit_kill_curve_to_real_threshold,
    fit_all_real_thresholds,
    REAL_TCE_ACTIVATION_THRESHOLDS,
)


def test_kill_probability_is_identical_to_capture_probability():
    """kill_probability() must be an exact, unmodified wrapper."""
    n = np.array([1.0, 10.0, 100.0, 1000.0])
    sigma, k_intra, z = 1.0, 0.01, 1.0
    expected = capture_probability(n, sigma, k_intra, z)
    actual = kill_probability(n, sigma, k_intra, z)
    assert np.allclose(actual, expected), (
        f"kill_probability diverged from capture_probability: "
        f"{actual} vs {expected}"
    )
    print("PASS: kill_probability is identical to capture_probability")


def test_fit_reproduces_all_three_real_thresholds_without_error():
    """
    Every one of the three real threshold values this project actually
    uses (81, 843, 8690 molecules/cell) must fit WITHOUT raising, and the
    fitted k_intra must reproduce theta=0.5 at that exact density. This is
    the direct regression test for the NaN/overflow bug caught three times
    while building this file -- see fit_kill_curve_to_real_threshold's own
    docstring for the full account of why a single spot-check (just n=81)
    was not enough to catch it the first two times.
    """
    thresholds_to_test = [81, 843, 8690]
    for threshold in thresholds_to_test:
        fit = fit_kill_curve_to_real_threshold(threshold)
        assert fit["converged"], f"Fit did not converge for threshold={threshold}"

        theta_at_threshold = kill_probability(
            threshold, fit["car_avidity"], fit["k_intra"], fit["z"]
        )
        assert not np.isnan(theta_at_threshold), (
            f"theta_at_threshold is NaN for threshold={threshold} -- the "
            f"exact bug caught and fixed three times while building this "
            f"file has reappeared."
        )
        assert np.isclose(theta_at_threshold, 0.5, atol=1e-6), (
            f"Fitted k_intra does not reproduce theta=0.5 at threshold="
            f"{threshold}: got {theta_at_threshold}"
        )
    print(f"PASS: fit reproduces theta=0.5 at all {len(thresholds_to_test)} "
          f"real thresholds, no NaN/overflow")


def test_fit_all_real_thresholds_runs_end_to_end():
    """fit_all_real_thresholds() (used directly by print_summary()) must
    run without error and produce a converged fit for every real
    threshold in REAL_TCE_ACTIVATION_THRESHOLDS."""
    results = fit_all_real_thresholds()
    expected_count = sum(
        1 for thresholds in REAL_TCE_ACTIVATION_THRESHOLDS.values()
        for key in thresholds if key.endswith("_molecules_per_cell")
    )
    assert len(results) == expected_count, (
        f"Expected {expected_count} fits, got {len(results)}"
    )
    for name, fit in results.items():
        assert fit["converged"], f"Fit '{name}' did not converge"
        assert not np.isnan(fit["k_intra"]), f"Fit '{name}' produced NaN k_intra"
    print(f"PASS: fit_all_real_thresholds produced {len(results)} converged, "
          f"non-NaN fits")


def test_escape_fraction_high_density_population_rarely_escapes():
    """A population with antigen density far ABOVE a calibrated threshold
    should show escape_fraction close to 0."""
    fit = fit_kill_curve_to_real_threshold(843)
    rng = np.random.default_rng(0)
    high_density_population = rng.negative_binomial(
        5, 5 / (5 + 843 * 20), size=5000
    ).astype(float)
    result = escape_fraction(
        high_density_population, fit["car_avidity"], fit["k_intra"], fit["z"]
    )
    assert result["escape_fraction"] < 0.05, (
        f"Expected near-zero escape fraction for a high-density population, "
        f"got {result['escape_fraction']}"
    )
    print(f"PASS: high-density population escape_fraction = "
          f"{result['escape_fraction']:.4f} (expected near 0)")


def test_escape_fraction_low_density_population_mostly_escapes():
    """A population with antigen density far BELOW a calibrated threshold
    should show escape_fraction close to 1."""
    fit = fit_kill_curve_to_real_threshold(843)
    rng = np.random.default_rng(0)
    low_density_population = rng.negative_binomial(
        5, 5 / (5 + 843 / 20), size=5000
    ).astype(float)
    result = escape_fraction(
        low_density_population, fit["car_avidity"], fit["k_intra"], fit["z"]
    )
    assert result["escape_fraction"] > 0.95, (
        f"Expected near-total escape fraction for a low-density population, "
        f"got {result['escape_fraction']}"
    )
    print(f"PASS: low-density population escape_fraction = "
          f"{result['escape_fraction']:.4f} (expected near 1)")


def test_escape_fraction_rejects_empty_population():
    """escape_fraction() must raise, not silently return a meaningless
    result, for an empty population."""
    fit = fit_kill_curve_to_real_threshold(843)
    try:
        escape_fraction(np.array([]), fit["car_avidity"], fit["k_intra"], fit["z"])
        assert False, "Expected ValueError for empty population, got no error"
    except ValueError:
        print("PASS: escape_fraction correctly rejects an empty population")


def test_escape_fraction_threshold_is_configurable():
    """A stricter (higher) kill_threshold should never produce a LOWER
    escape fraction than a looser (lower) one, for the same population --
    this is a basic monotonicity sanity check on the threshold parameter."""
    fit = fit_kill_curve_to_real_threshold(843)
    rng = np.random.default_rng(1)
    population = rng.negative_binomial(5, 5 / (5 + 843), size=5000).astype(float)

    loose = escape_fraction(
        population, fit["car_avidity"], fit["k_intra"], fit["z"], kill_threshold=0.1
    )
    strict = escape_fraction(
        population, fit["car_avidity"], fit["k_intra"], fit["z"], kill_threshold=0.9
    )
    assert strict["escape_fraction"] >= loose["escape_fraction"], (
        f"A stricter kill_threshold produced a LOWER escape fraction "
        f"({strict['escape_fraction']}) than a looser one "
        f"({loose['escape_fraction']}) -- this should never happen."
    )
    print(f"PASS: escape_fraction is monotonic in kill_threshold "
          f"(loose={loose['escape_fraction']:.4f}, strict={strict['escape_fraction']:.4f})")


if __name__ == "__main__":
    tests = [
        test_kill_probability_is_identical_to_capture_probability,
        test_fit_reproduces_all_three_real_thresholds_without_error,
        test_fit_all_real_thresholds_runs_end_to_end,
        test_escape_fraction_high_density_population_rarely_escapes,
        test_escape_fraction_low_density_population_mostly_escapes,
        test_escape_fraction_rejects_empty_population,
        test_escape_fraction_threshold_is_configurable,
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
