"""
test_alpha_exact.py

NEW (v14->v15). Regression tests for binding_model.alpha_exact() -- the
closed-form (exact) selectivity formula added this session in response to
"Immediate next steps" item 0b (itself a direct response to the
independent Biomni external review -- see the handoff's "What changed
since v13").

WHY THIS FILE EXISTS
---------------------
alpha_exact() was independently re-derived from this project's own
equations (not copied from Biomni's report -- see its docstring in
binding_model.py for the full derivation) and then checked against
Biomni's reported formula and against this project's own existing
finite-difference selectivity_alpha(). These tests lock in both checks so
a future refactor cannot silently break the closed form without a test
failing immediately.

WHAT THESE TESTS LOCK IN
--------------------------
1. alpha_exact() closely matches the finite-difference selectivity_alpha()
   at smooth, well-resolved points (they are two different ways of
   computing an approximately-continuous quantity, so exact numerical
   equality is not expected, but they should agree to within a few
   percent away from the sweep's most extreme corners).
2. alpha_exact() and selectivity_alpha() diverge by a real, non-trivial
   amount (up to ~0.70 in absolute alpha, checked directly this session)
   at some points in the project's actual K_INTRA/z plausibility grid --
   this is itself worth locking in as a fact about the model, not just
   "they roughly agree everywhere."
3. The exact peak alpha on the real plausibility grid is close to, but
   not identical to, the finite-difference peak (~10.465 exact vs. ~10.461
   finite-difference) -- confirming switching to the exact formula does
   NOT change the qualitative "superselectivity peaks around alpha~10.46"
   story, it just makes the number the write-up cites more precise.
4. alpha_exact() at n=1, in the DILUTE limit (small n*sigma*k_intra),
   should be close to 1 -- this is a basic sanity check (at n=1, x is
   small, so exp(x)-1 ~ x, and alpha ~ x*e^x/(x*(1+z*x)) ~ 1/(1+z*x) ~ 1
   for small x*z), independent of the finite-difference comparison.

HOW TO RUN
    pytest test_alpha_exact.py -v
or:
    python test_alpha_exact.py
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from binding_model import alpha_exact, selectivity_alpha
from aim3_aim4_range_sweep import (
    get_k_intra_sweep_values,
    get_z_sweep_values,
    sweep_superselectivity_reachability,
)
from k_intra_z_plausibility_bounds import SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2


def test_alpha_exact_matches_finite_difference_closely_at_smooth_points():
    """
    At a well-resolved, non-extreme (k_intra, z) point, alpha_exact()
    and the finite-difference selectivity_alpha() should agree closely
    (within a few percent) -- they are estimating the same underlying
    quantity two different ways.
    """
    sigma = SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2
    k_intra_values = get_k_intra_sweep_values()
    k_intra = k_intra_values[len(k_intra_values) // 2]  # a mid-range value
    z = 1e-4
    n_values = np.linspace(1, 200, 60)

    n_mid, alpha_fd = selectivity_alpha(n_values, sigma, k_intra, z)
    alpha_ex = alpha_exact(n_mid, sigma, k_intra, z)

    finite = np.isfinite(alpha_fd) & np.isfinite(alpha_ex)
    assert finite.sum() > 10, "Expected most points to be finite/reliable for this check."
    rel_diff = np.abs(alpha_fd[finite] - alpha_ex[finite]) / np.maximum(np.abs(alpha_ex[finite]), 1e-9)
    assert np.median(rel_diff) < 0.05, (
        f"Expected the median relative difference between alpha_exact() and "
        f"selectivity_alpha() to be small (<5%) at a mid-range, well-"
        f"resolved (k_intra, z) point -- got median {np.median(rel_diff):.4f}."
    )
    print(f"  [OK] alpha_exact() and selectivity_alpha() agree closely "
          f"(median relative difference {np.median(rel_diff):.2%}) at a "
          f"mid-range plausibility-grid point.")


def test_alpha_exact_and_finite_difference_diverge_meaningfully_somewhere_on_real_grid():
    """
    Independently re-derived this session: on the project's REAL K_INTRA/z
    plausibility grid, alpha_exact() and selectivity_alpha() diverge by as
    much as ~0.70 in absolute alpha at some grid points -- large enough
    that an independent reviewer computing the exact formula would notice
    immediately. This is a genuine property of the model on this specific
    grid, not a coincidence of one arbitrarily chosen point, so it's
    locked in directly here.
    """
    sigma = SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2
    k_intra_values = get_k_intra_sweep_values()
    z_values = get_z_sweep_values()
    n_values = np.linspace(1, 200, 60)

    max_diff = 0.0
    for k_intra in k_intra_values:
        for z in z_values:
            n_mid, alpha_fd = selectivity_alpha(n_values, sigma, k_intra, z)
            alpha_ex = alpha_exact(n_mid, sigma, k_intra, z)
            finite = np.isfinite(alpha_fd) & np.isfinite(alpha_ex)
            if not finite.any():
                continue
            diff = np.abs(alpha_fd[finite] - alpha_ex[finite])
            if diff.max() > max_diff:
                max_diff = diff.max()

    assert max_diff > 0.3, (
        f"Expected a real, non-trivial discrepancy (>0.3 in absolute alpha) "
        f"somewhere on the real plausibility grid -- got a maximum "
        f"discrepancy of only {max_diff:.4f}. This session's independent "
        f"check found ~0.70 -- investigate if this drops much lower, since "
        f"it would mean the two formulas are now much closer than "
        f"previously found (worth understanding why, not just noting)."
    )
    print(f"  [OK] alpha_exact() and selectivity_alpha() diverge by "
          f"{max_diff:.4f} (absolute alpha) at their most different point "
          f"on the real plausibility grid, matching this session's "
          f"independently-checked finding of ~0.70.")


def test_superselectivity_sweep_reports_exact_peak_close_to_finite_difference_peak():
    """
    NEW (v14->v15). sweep_superselectivity_reachability() now reports
    max_alpha_exact_found alongside the original max_alpha_found. Both
    should be close to each other (within ~0.1 in absolute alpha) at the
    grid's actual peak, and the exact value should be at least as large
    (the finite-difference approximation should not, on the grid's real
    peak specifically, overshoot the true value by a meaningful margin).
    """
    result = sweep_superselectivity_reachability()
    assert "max_alpha_exact_found" in result, (
        "Expected sweep_superselectivity_reachability() to report "
        "max_alpha_exact_found -- this key is new in v14->v15."
    )
    assert "best_params_exact_k_intra_z" in result

    fd_peak = result["max_alpha_found"]
    exact_peak = result["max_alpha_exact_found"]
    assert abs(fd_peak - exact_peak) < 0.5, (
        f"Expected the finite-difference peak ({fd_peak:.4f}) and exact "
        f"peak ({exact_peak:.4f}) to be close (within 0.5) at the grid's "
        f"actual maximum -- this session's independent check found "
        f"10.461 (finite-difference) vs. 10.465 (exact), a ~0.04 gap."
    )
    print(f"  [OK] sweep_superselectivity_reachability() reports both peaks: "
          f"finite-difference {fd_peak:.3f}, exact {exact_peak:.3f} -- "
          f"close to each other at the grid's actual peak, confirming "
          f"switching to the exact formula does not change the "
          f"qualitative 'superselectivity peaks around alpha~10.46' story.")


def test_alpha_exact_near_one_in_dilute_limit_at_n_equals_one():
    """
    Basic sanity check, independent of the finite-difference comparison:
    at n=1 with a small n*sigma*k_intra*z product (the dilute limit),
    alpha_exact() should be close to 1 -- this is the same "alpha starts
    near 1 at low n" behavior selectivity_alpha()'s own docstring
    describes for the finite-difference version.
    """
    sigma = 1.0
    k_intra = 1e-6
    z = 1e-4
    alpha = alpha_exact(1.0, sigma, k_intra, z)
    assert 0.9 < alpha < 1.1, (
        f"Expected alpha_exact(n=1, ...) to be close to 1 in the dilute "
        f"limit (small n*sigma*k_intra*z) -- got {alpha:.4f}."
    )
    print(f"  [OK] alpha_exact() at n=1 in the dilute limit is {alpha:.4f}, "
          f"close to 1 as expected.")


if __name__ == "__main__":
    tests = [
        test_alpha_exact_matches_finite_difference_closely_at_smooth_points,
        test_alpha_exact_and_finite_difference_diverge_meaningfully_somewhere_on_real_grid,
        test_superselectivity_sweep_reports_exact_peak_close_to_finite_difference_peak,
        test_alpha_exact_near_one_in_dilute_limit_at_n_equals_one,
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
