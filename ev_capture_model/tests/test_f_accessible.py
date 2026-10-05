"""
test_f_accessible.py

NEW (this session). Regression tests for binding_model.py's new
accessibility/protein-corona multiplier: accessible_receptors(),
capture_probability_accessible(), and
selectivity_ratio_sensitivity_to_accessibility() -- "Immediate next steps"
item 3 (unattempted from v4 through v15).

WHAT THESE TESTS LOCK IN
--------------------------
1. accessible_receptors() at the default f_accessible=1.0 returns n
   UNCHANGED -- the core backward-compatibility guarantee this addition
   depends on (no existing call site or test needed to change).
2. accessible_receptors() at f_accessible<1.0 correctly scales n down,
   and rejects out-of-range f_accessible values explicitly.
3. capture_probability_accessible() at f_accessible=1.0 matches
   capture_probability() EXACTLY (bit-for-bit), confirming the wrapper
   really is a no-op at the default, not just approximately so.
4. capture_probability_accessible() at f_accessible<1.0 gives a LOWER
   capture probability than at f_accessible=1.0 (fewer accessible
   receptors should never increase capture probability -- a basic
   monotonicity requirement of the underlying model, since q_bound is
   monotonically increasing in n).
5. selectivity_ratio_sensitivity_to_accessibility() with EQUAL tumor/
   normal f_accessible gives distortion_factor == 1.0 exactly (no
   distortion when accessibility doesn't differ between populations --
   the "errors cancel" case the handoff describes).
6. selectivity_ratio_sensitivity_to_accessibility() with DIFFERENT tumor/
   normal f_accessible gives a real, non-trivial distortion_factor,
   confirming the handoff's flagged risk ("if f_accessible differs... it
   contaminates the selectivity ratio") is a real, checkable property of
   this model, not just a hypothetical concern in prose.

HOW TO RUN
    pytest test_f_accessible.py -v
or:
    python test_f_accessible.py
"""

import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from binding_model import (
    accessible_receptors,
    capture_probability,
    capture_probability_accessible,
    selectivity_ratio_sensitivity_to_accessibility,
)


def test_default_f_accessible_leaves_n_unchanged():
    """
    accessible_receptors(n, f_accessible=1.0) (the default) must return n
    exactly unchanged -- the core backward-compatibility guarantee.
    """
    n = np.array([0, 1, 10, 50, 200], dtype=float)
    n_accessible = accessible_receptors(n)
    assert np.array_equal(n_accessible, n), (
        "Expected accessible_receptors() at the default f_accessible=1.0 "
        "to return n bit-for-bit unchanged."
    )
    print("  [OK] Default f_accessible=1.0 leaves n exactly unchanged.")


def test_f_accessible_below_one_scales_n_down_correctly():
    """
    accessible_receptors(n, f_accessible=0.3) should return exactly
    0.3 * n.
    """
    n = np.array([10.0, 100.0, 200.0])
    n_accessible = accessible_receptors(n, f_accessible=0.3)
    assert np.allclose(n_accessible, 0.3 * n), (
        f"Expected n_accessible = 0.3 * n = {0.3 * n}, got {n_accessible}."
    )
    print("  [OK] f_accessible=0.3 correctly scales n down to 0.3*n.")


def test_f_accessible_out_of_range_raises():
    """
    f_accessible values outside (0, 1] should raise ValueError explicitly,
    per accessible_receptors()'s own documented contract -- this is very
    likely a units/convention mistake at the call site, not a real
    parameter choice, and should be caught loudly rather than silently
    producing a physically nonsensical n_accessible.
    """
    with pytest.raises(ValueError):
        accessible_receptors(100.0, f_accessible=1.5)
    with pytest.raises(ValueError):
        accessible_receptors(100.0, f_accessible=0.0)
    with pytest.raises(ValueError):
        accessible_receptors(100.0, f_accessible=-0.2)
    print("  [OK] f_accessible values outside (0, 1] correctly raise ValueError.")


def test_capture_probability_accessible_matches_original_at_default():
    """
    capture_probability_accessible(n, sigma, k_intra, z) at the default
    f_accessible=1.0 must match capture_probability(n, sigma, k_intra, z)
    EXACTLY (bit-for-bit) -- confirming the wrapper is a true no-op at the
    default, not merely a close numerical approximation.
    """
    n, sigma, k_intra, z = 50.0, 1e11, 1e-12, 1.0
    original = capture_probability(n, sigma, k_intra, z)
    wrapped = capture_probability_accessible(n, sigma, k_intra, z)
    assert original == wrapped, (
        f"Expected capture_probability_accessible() at f_accessible=1.0 "
        f"to match capture_probability() EXACTLY -- got {wrapped} vs "
        f"{original}."
    )
    print("  [OK] capture_probability_accessible() at the default matches "
          "capture_probability() bit-for-bit.")


def test_lower_accessibility_never_increases_capture_probability():
    """
    Basic monotonicity requirement: since q_bound(n, ...) is monotonically
    increasing in n (more receptors = more binding opportunity, never
    fewer), reducing the accessible fraction of receptors should never
    INCREASE capture probability -- checked across several f_accessible
    values.
    """
    n, sigma, k_intra, z = 100.0, 1e11, 1e-12, 1.0
    theta_full = capture_probability_accessible(n, sigma, k_intra, z, f_accessible=1.0)

    for f_acc in (0.9, 0.5, 0.1, 0.01):
        theta_reduced = capture_probability_accessible(n, sigma, k_intra, z, f_accessible=f_acc)
        assert theta_reduced <= theta_full, (
            f"Expected capture probability at f_accessible={f_acc} "
            f"({theta_reduced}) to be <= capture probability at "
            f"f_accessible=1.0 ({theta_full})."
        )
    print("  [OK] Reducing f_accessible never increases capture probability, "
          "across several tested values.")


def test_equal_accessibility_gives_no_distortion():
    """
    When tumor and normal populations have the SAME f_accessible, the
    selectivity ratio should be UNCHANGED by accessibility at all --
    distortion_factor should be exactly 1.0 (the "errors cancel" case).
    """
    result = selectivity_ratio_sensitivity_to_accessibility(
        n_tumor=[80, 90, 100], n_normal=[10, 15, 20],
        sigma=1e11, k_intra=1e-12, z=1.0,
        f_accessible_tumor=0.6, f_accessible_normal=0.6,
    )
    assert np.isclose(result["distortion_factor"], 1.0, atol=1e-9), (
        f"Expected distortion_factor == 1.0 exactly when tumor and normal "
        f"f_accessible are equal -- got {result['distortion_factor']}."
    )
    print("  [OK] Equal tumor/normal f_accessible gives distortion_factor "
          "== 1.0 exactly, confirming no distortion when accessibility "
          "doesn't differ between populations.")


def test_unequal_accessibility_produces_real_distortion():
    """
    Confirms the handoff's flagged risk directly: when tumor and normal
    f_accessible DIFFER substantially, the selectivity ratio should be
    meaningfully distorted (distortion_factor far from 1.0) -- this is
    the actual, numeric answer to "if f_accessible differs between tumor
    and normal EVs... it contaminates the selectivity ratio."
    """
    result = selectivity_ratio_sensitivity_to_accessibility(
        n_tumor=[80, 90, 100], n_normal=[10, 15, 20],
        sigma=1e11, k_intra=1e-12, z=1.0,
        f_accessible_tumor=0.9, f_accessible_normal=0.3,
    )
    assert abs(result["distortion_factor"] - 1.0) > 0.1, (
        f"Expected a real (>10%) distortion when tumor f_accessible=0.9 "
        f"and normal f_accessible=0.3 differ substantially -- got "
        f"distortion_factor={result['distortion_factor']:.4f}."
    )
    print(f"  [OK] Unequal tumor/normal f_accessible (0.9 vs 0.3) produces "
          f"a real distortion_factor of {result['distortion_factor']:.3f} -- "
          f"confirming the handoff's flagged risk is a checkable, real "
          f"property of this model.")


if __name__ == "__main__":
    tests = [
        test_default_f_accessible_leaves_n_unchanged,
        test_f_accessible_below_one_scales_n_down_correctly,
        test_f_accessible_out_of_range_raises,
        test_capture_probability_accessible_matches_original_at_default,
        test_lower_accessibility_never_increases_capture_probability,
        test_equal_accessibility_gives_no_distortion,
        test_unequal_accessibility_produces_real_distortion,
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
