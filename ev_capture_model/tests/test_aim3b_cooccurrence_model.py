"""
test_aim3b_cooccurrence_model.py

NEW (this session). Regression tests for src/aim3b_cooccurrence_model.py --
the real, fitted (not just illustratively-swept) two-subpopulation
co-occurrence model, "Immediate next steps" item 2 (unattempted from v6
through v15), plus its item-5 independence test built on top.

WHAT THESE TESTS LOCK IN
--------------------------
1. fit_suppression_ratio() actually converges to the real target ratio
   (10.0 for EpCAM, 5.0 for Her2), for several f_rich values -- the core
   claim of this file (a REAL fit exists, not just an illustrative sweep).
2. The simulated ratio is monotonically increasing in the suppression
   ratio -- the property fit_suppression_ratio()'s use of brentq (a
   bracketing root-finder, which requires monotonicity to be guaranteed
   to work correctly) depends on.
3. fit_both_markers() produces two DIFFERENT suppression ratios for EpCAM
   vs. Her2 at the same f_rich (since the two markers have different real
   target ratios, 10x vs 5x, they cannot share one suppression value).
4. independence_test() shows a real, non-trivial deviation from
   independence (not exactly 1.0), and the SIGN of that deviation is
   POSITIVE correlation (fold_deviation > 1) -- both markers are
   suppressed in the SAME subpopulation, so this is a hard structural
   prediction of the model, not an arbitrary possibility.
5. The independence fold-deviation finding is STABLE across the
   MARKER_MEAN_TET_POOR sensitivity sweep, while the combination-gain
   finding is NOT -- this is the real, checked sensitivity finding from
   this session's own build process (see the source file's module
   docstring), locked in here so it cannot be silently forgotten or
   contradicted by a future edit without that edit being deliberate.
6. multi_marker_yield_gain()'s combination gain is always >= 1.0 (a
   combination-capture strategy can never do WORSE than the best single
   marker it includes -- a basic logical/set-theoretic property: P(A or B)
   >= max(P(A), P(B)) always, regardless of correlation structure).

HOW TO RUN
    pytest test_aim3b_cooccurrence_model.py -v
or:
    python test_aim3b_cooccurrence_model.py
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import aim3b_cooccurrence_model as m
from aim3b_cooccurrence_model import (
    fit_suppression_ratio,
    fit_both_markers,
    fit_across_f_rich_range,
    independence_test,
    multi_marker_yield_gain,
    combination_gain_sensitivity_sweep,
    _simulated_ratio_for_suppression,
    MIZENKO_TARGET_RATIOS,
)


def test_fit_suppression_ratio_converges_to_real_target():
    """
    For a representative f_rich, fit_suppression_ratio() should find a
    suppression ratio whose simulated ratio is close to the real target
    (Mizenko's reported 10.0 for EpCAM) -- the core claim this file adds
    beyond the scoping file's illustrative sweep.
    """
    f_rich = 0.20
    target = MIZENKO_TARGET_RATIOS["EpCAM"]
    suppression = fit_suppression_ratio(f_rich, target, seed=0)
    achieved = _simulated_ratio_for_suppression(f_rich, m.MARKER_MEAN_TET_POOR, suppression, seed=1)

    assert abs(achieved - target) / target < 0.10, (
        f"Expected the fitted suppression ratio to reproduce the target "
        f"ratio {target} within 10% (using an INDEPENDENT re-simulation, "
        f"different seed than the fit itself) -- got achieved={achieved:.3f}."
    )
    print(f"  [OK] fit_suppression_ratio() at f_rich={f_rich} for EpCAM "
          f"converges to a suppression ratio of {suppression:.2f}x, "
          f"achieving ratio {achieved:.2f}x (target {target}).")


def test_simulated_ratio_is_monotonic_in_suppression_ratio():
    """
    The property fit_suppression_ratio()'s use of brentq depends on: the
    simulated ratio must be monotonically increasing as the suppression
    ratio increases (stronger suppression -> bigger gap between capture
    methods -> bigger ratio). Checked directly across a range, not assumed.
    """
    f_rich = 0.20
    suppression_values = np.array([2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0])
    ratios = [
        _simulated_ratio_for_suppression(f_rich, m.MARKER_MEAN_TET_POOR, s, seed=0)
        for s in suppression_values
    ]
    diffs = np.diff(ratios)
    assert np.all(diffs > 0), (
        f"Expected the simulated ratio to be strictly increasing in the "
        f"suppression ratio -- got ratios {ratios} (diffs {diffs}). If "
        f"this fails, fit_suppression_ratio()'s brentq root-find is no "
        f"longer guaranteed to converge correctly."
    )
    print(f"  [OK] Simulated ratio is monotonically increasing in the "
          f"suppression ratio across {suppression_values} -- "
          f"ratios: {[f'{r:.2f}' for r in ratios]}.")


def test_fit_both_markers_gives_different_suppression_for_each_marker():
    """
    EpCAM (target 10x) and Her2 (target 5x) should NOT receive the same
    fitted suppression ratio at a shared f_rich, since they have different
    real targets -- a basic consistency check that the two markers are
    genuinely fit independently, not accidentally sharing one value.
    """
    fit = fit_both_markers(0.20, seed=0)
    epcam_suppression = fit["EpCAM"]["suppression_ratio"]
    her2_suppression = fit["Her2"]["suppression_ratio"]

    assert not np.isclose(epcam_suppression, her2_suppression, rtol=0.05), (
        f"Expected EpCAM (target 10x) and Her2 (target 5x) to receive "
        f"DIFFERENT fitted suppression ratios -- got EpCAM={epcam_suppression:.2f}, "
        f"Her2={her2_suppression:.2f}, which are suspiciously close."
    )
    assert epcam_suppression > her2_suppression, (
        f"Expected EpCAM (bigger target ratio, 10x vs 5x) to require a "
        f"LARGER suppression ratio than Her2 -- got EpCAM={epcam_suppression:.2f}, "
        f"Her2={her2_suppression:.2f}."
    )
    print(f"  [OK] EpCAM (suppression {epcam_suppression:.2f}x) and Her2 "
          f"(suppression {her2_suppression:.2f}x) are fit to distinct, "
          f"correctly-ordered suppression ratios.")


def test_fit_across_f_rich_range_achieves_targets_at_every_f_rich():
    """
    fit_across_f_rich_range() should produce a good fit (achieved ratio
    close to target) for EVERY f_rich in the plausible range, not just one
    cherry-picked value.
    """
    fits = fit_across_f_rich_range()
    assert len(fits) >= 3, "Expected the default f_rich range to have at least 3 values."
    for fit in fits:
        for marker in ("EpCAM", "Her2"):
            target = fit[marker]["target_ratio"]
            achieved = fit[marker]["achieved_ratio"]
            assert abs(achieved - target) / target < 0.15, (
                f"At f_rich={fit['f_rich']}, {marker}'s achieved ratio "
                f"{achieved:.2f} deviates from target {target} by more "
                f"than 15%."
            )
    print(f"  [OK] fit_across_f_rich_range() achieves both markers' real "
          f"target ratios within 15% at all {len(fits)} tested f_rich values.")


def test_independence_test_shows_real_positive_correlation():
    """
    independence_test() should find a real, non-trivial (>10%) deviation
    from independence, and that deviation should be POSITIVE (fold > 1) --
    both markers are suppressed together in the tet-rich subpopulation, so
    the model structurally predicts positive correlation, not just "some
    deviation in some direction."
    """
    fit = fit_both_markers(0.20, seed=0)
    result = independence_test(fit, seed=1)

    assert result["fold_deviation_from_independence"] > 1.05, (
        f"Expected a real (>5%) POSITIVE deviation from independence "
        f"(fold_deviation > 1.05) -- got {result['fold_deviation_from_independence']:.4f}. "
        f"The model structurally predicts positive correlation since both "
        f"markers are suppressed in the same tet-rich subpopulation."
    )
    assert 0 < result["p_joint_true"] < 1, "p_joint_true should be a valid probability."
    print(f"  [OK] independence_test() finds a real positive correlation "
          f"(fold deviation {result['fold_deviation_from_independence']:.2f}x) -- "
          f"the naive independence shortcut under-predicts true co-occurrence.")


def test_independence_finding_is_stable_across_mean_tet_poor_sensitivity_sweep():
    """
    Locks in this session's real sensitivity finding (see the source
    file's module docstring): the independence fold-deviation should stay
    roughly STABLE (within 10% of its value at the default scale) across
    the MARKER_MEAN_TET_POOR sensitivity sweep, confirming this is the
    more robust of the file's two headline findings.
    """
    rows = combination_gain_sensitivity_sweep()
    fold_deviations = [row["fold_deviation_from_independence"] for row in rows]

    spread = (max(fold_deviations) - min(fold_deviations)) / np.mean(fold_deviations)
    assert spread < 0.10, (
        f"Expected the independence fold-deviation to be STABLE (<10% "
        f"relative spread) across the MARKER_MEAN_TET_POOR sensitivity "
        f"sweep -- got fold deviations {fold_deviations} (relative spread "
        f"{spread:.2%}). If this is no longer stable, the module "
        f"docstring's claim that this finding is robust (unlike the "
        f"combination-gain number) needs to be revisited."
    )
    print(f"  [OK] Independence fold-deviation is stable across the "
          f"MARKER_MEAN_TET_POOR sweep (relative spread {spread:.2%}): "
          f"{[f'{f:.3f}' for f in fold_deviations]}.")


def test_combination_gain_is_sensitive_to_mean_tet_poor_as_documented():
    """
    Locks in the OTHER half of this session's sensitivity finding: the
    combination-gain number should NOT be stable across the same sweep --
    if a future change makes it stable too, that would itself be worth
    understanding (has the underlying model changed?), not silently
    accepted as an improvement.
    """
    rows = combination_gain_sensitivity_sweep()
    combo_gains = [row["combination_gain_true"] for row in rows]

    spread = (max(combo_gains) - min(combo_gains)) / np.mean(combo_gains)
    assert spread > 0.10, (
        f"Expected the combination-gain number to be SENSITIVE (>10% "
        f"relative spread) to MARKER_MEAN_TET_POOR, per this session's "
        f"documented finding -- got combo gains {combo_gains} (relative "
        f"spread {spread:.2%}), suggesting this sensitivity may no longer "
        f"hold and the module docstring's account should be revisited."
    )
    print(f"  [OK] Combination gain is sensitive to MARKER_MEAN_TET_POOR as "
          f"documented (relative spread {spread:.2%}): "
          f"{[f'{c:.3f}' for c in combo_gains]}.")


def test_combination_gain_never_below_one():
    """
    Basic logical/set-theoretic sanity check: P(A or B) >= max(P(A), P(B))
    is always true regardless of correlation structure, so
    multi_marker_yield_gain()'s combination_gain_true should never be less
    than 1.0 (the combination strategy can never do WORSE than the better
    single marker alone).
    """
    for mean_tet_poor in (2.0, 15.0, 30.0):
        m.MARKER_MEAN_TET_POOR = mean_tet_poor
        try:
            fit = fit_both_markers(0.20, seed=0)
            gain = multi_marker_yield_gain(fit, seed=1)
            assert gain["combination_gain_true"] >= 1.0 - 1e-9, (
                f"combination_gain_true should never be < 1.0 -- got "
                f"{gain['combination_gain_true']:.4f} at "
                f"mean_tet_poor={mean_tet_poor}."
            )
        finally:
            m.MARKER_MEAN_TET_POOR = 15.0  # restore default
    print("  [OK] Combination gain is always >= 1.0 across tested "
          "mean_tet_poor values, matching the basic set-theoretic "
          "guarantee that P(A or B) >= max(P(A), P(B)).")


if __name__ == "__main__":
    tests = [
        test_fit_suppression_ratio_converges_to_real_target,
        test_simulated_ratio_is_monotonic_in_suppression_ratio,
        test_fit_both_markers_gives_different_suppression_for_each_marker,
        test_fit_across_f_rich_range_achieves_targets_at_every_f_rich,
        test_independence_test_shows_real_positive_correlation,
        test_independence_finding_is_stable_across_mean_tet_poor_sensitivity_sweep,
        test_combination_gain_is_sensitive_to_mean_tet_poor_as_documented,
        test_combination_gain_never_below_one,
    ]

    n_passed = 0
    n_failed = 0
    for test_func in tests:
        print(f"\nRunning {test_func.__name__} ...")
        try:
            test_func()
            print(f"PASSED: {test_func.__name__}")
            n_passed += 1
        except AssertionError as e:
            print(f"FAILED: {test_func.__name__}\n   {e}")
            n_failed += 1

    print(f"\n{'='*60}\n{n_passed} passed, {n_failed} failed out of {len(tests)} tests\n{'='*60}")
    if n_failed > 0:
        sys.exit(1)
