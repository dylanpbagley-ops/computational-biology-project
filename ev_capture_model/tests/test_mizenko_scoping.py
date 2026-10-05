"""
test_mizenko_scoping.py

Regression test for src/mizenko_cooccurrence_scoping.py's headline finding
(see that file's module docstring for the full context): a two-subpopulation
mixture model needs BOTH a fairly small tetraspanin-rich fraction AND a
strong (~50-100x), specific suppression of EpCAM within that subpopulation
to reproduce Mizenko et al.'s reported ~10x biotin/tetraspanin EpCAM-
positive-fraction ratio -- mild/moderate anti-correlation alone falls far
short (~1.0-1.6x). This is exploratory SCOPING work for Aim 3b, not a
validation of binding_model.py, and this test exists only so the specific
numbers reported in the docstring/handoff stay checked rather than becoming
a stale, unverified claim if the script changes later.

HOW TO RUN
    pytest test_mizenko_scoping.py -v
or:
    python test_mizenko_scoping.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from mizenko_cooccurrence_scoping import epcam_positive_fraction_ratio


def test_mild_anticorrelation_falls_far_short_of_mizenko_magnitude():
    """
    At a large tet-rich fraction (0.5) and only mild EpCAM suppression
    (15 vs 5 mean copies, a 3x difference), the resulting EpCAM-positive-
    fraction ratio should stay well below 2x -- nowhere near Mizenko's
    reported ~10x. This locks in the "mild anti-correlation is not enough"
    half of the scoping finding.
    """
    result = epcam_positive_fraction_ratio(
        frac_tet_rich=0.5, epcam_mean_tet_poor=15, epcam_mean_tet_rich=5, seed=0
    )
    assert result["ratio"] < 2.0, (
        f"Expected mild anti-correlation to produce a small ratio (<2x), "
        f"got {result['ratio']:.2f}x -- if this now exceeds 2x, the scoping "
        f"model's behavior has changed and the docstring claim needs "
        f"revisiting."
    )
    print(f"  [OK] mild anti-correlation (tet_rich=0.5, EpCAM 15/5) gives "
          f"ratio={result['ratio']:.2f}x, well below Mizenko's ~10x -- "
          f"confirms mild correlation alone cannot explain the magnitude.")


def test_strong_structural_anticorrelation_can_reach_mizenko_magnitude():
    """
    At a small tet-rich fraction (0.3) and strong EpCAM suppression (12 vs
    0.1 mean copies, a 120x difference), the ratio should approach or
    exceed Mizenko's reported ~10x (allowing this specific configuration to
    land at ~7-8x, since this locks in "the right ballpark," not an exact
    digit match). This confirms the mixture-model FORM is at least capable
    of reaching the right order of magnitude when the anti-correlation is
    strong and structural, not just capable of small effects.
    """
    result = epcam_positive_fraction_ratio(
        frac_tet_rich=0.3, epcam_mean_tet_poor=12, epcam_mean_tet_rich=0.1, seed=0
    )
    assert result["ratio"] > 5.0, (
        f"Expected strong, structural anti-correlation to reach at least "
        f"5x (approaching Mizenko's reported ~10x), got "
        f"{result['ratio']:.2f}x -- if this now falls below 5x, the "
        f"scoping model's behavior has changed and the docstring claim "
        f"needs revisiting."
    )
    print(f"  [OK] strong structural anti-correlation (tet_rich=0.3, EpCAM "
          f"12/0.1) gives ratio={result['ratio']:.2f}x, in the right "
          f"ballpark of Mizenko's ~10x -- confirms the mixture-model form "
          f"CAN reach the right order of magnitude, but only under strong, "
          f"specific assumptions, not mild ones.")


def test_zero_anticorrelation_gives_ratio_near_one():
    """
    Sanity/null check: if EpCAM mean copy number is IDENTICAL in both
    subpopulations (no anti-correlation at all), the ratio must be very
    close to 1.0 -- capture selection on tetraspanin alone should not bias
    EpCAM detection if the two markers are truly independent. This is the
    null-model check that makes the other two tests meaningful.
    """
    result = epcam_positive_fraction_ratio(
        frac_tet_rich=0.2, epcam_mean_tet_poor=10, epcam_mean_tet_rich=10, seed=0
    )
    assert abs(result["ratio"] - 1.0) < 0.05, (
        f"Expected ratio very close to 1.0 when EpCAM is independent of "
        f"tetraspanin status, got {result['ratio']:.3f} -- a real "
        f"deviation here would suggest a bug in the capture-selection "
        f"logic, not a physical finding."
    )
    print(f"  [OK] identical EpCAM distributions in both subpopulations "
          f"(true independence) give ratio={result['ratio']:.3f}, "
          f"essentially 1.0 as expected -- confirms the model doesn't "
          f"spuriously manufacture a ratio when there's no real "
          f"anti-correlation.")


if __name__ == "__main__":
    tests = [
        test_zero_anticorrelation_gives_ratio_near_one,
        test_mild_anticorrelation_falls_far_short_of_mizenko_magnitude,
        test_strong_structural_anticorrelation_can_reach_mizenko_magnitude,
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
