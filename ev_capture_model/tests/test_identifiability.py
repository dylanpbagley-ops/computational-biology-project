"""
test_identifiability.py

NEW (v14->v15). Regression tests locking in the identifiability/
confounding argument documented in
docs/Limitations_and_Calibration_Gap.md, itself written in response to
"Immediate next steps" item 1 (following the independent Biomni external
review -- see the handoff's "What changed since v13").

WHAT THIS LOCKS IN
--------------------
sigma, K_INTRA, and k_geom are NOT separately identifiable from the
model's output -- only certain PRODUCTS of them (sigma*K_INTRA, and
k_geom*K_INTRA at patch_exponent=2) ever enter the equations. This means
rescaling sigma and K_INTRA in opposite directions (or k_geom and K_INTRA
in opposite directions, at patch_exponent=2) while holding the relevant
product fixed must produce BIT-IDENTICAL model output. This was checked
directly and algebraically this session (see the limitations doc for the
full derivation) -- these tests make that check a permanent, automatic
regression test rather than a one-off scratch computation.

If either test here ever fails, it means someone changed q_bound() or
contact_patch_receptors() in a way that breaks this structural property
of the model -- which would be a big deal (it would mean the
identifiability argument in the limitations write-up is no longer
accurate) and should be investigated immediately, not silently accepted.

A NOTE ON "BIT-IDENTICAL", AND A SECOND, LARGER PRECISION ISSUE FOUND
WHILE WRITING THESE TESTS: an exact np.array_equal() check turned out to
be too strict for some rescaling factors -- e.g. scaling k_geom by 0.1
and k_intra by 1/0.1 produced a difference of 1.1e-16 (one float64
machine-epsilon-scale rounding step), because 1/0.1 is not exactly
representable in binary floating point, so the "undo" multiplication
doesn't land on literally the same bit pattern even though the MATH is
exact. An ABSOLUTE tolerance of 1e-9 handles that case fine.

But a much larger, and more instructive, issue turned up at a bigger
rescaling factor (1000x): the ABSOLUTE difference in q_bound() came out
to ~2e44 -- nowhere near 1e-9. This is NOT a violation of the
identifiability property. Checked directly: the two rescaled versions of
x = n*sigma*k_intra print as identical (4.44, 17.76, 71.04, 133.2 in both
cases) but are NOT bit-identical (x1 == x2 gives all False) -- because
sigma*factor and k_intra/factor each individually accumulate a tiny
floating-point rounding error, and even though those errors are far too
small to see when printing x, q_bound = exp(x) - 1 is EXPONENTIALLY
sensitive to x at large x (at x=133.2, exp(x) is around 7e57 -- a relative
error in x of order 1e-16 becomes an ABSOLUTE error in exp(x) of order
1e57 * 1e-16 = 1e41, which is exactly the scale of difference actually
observed). This is the same general phenomenon this project's other files
already document carefully (see MIN_RELIABLE_EXPONENT_AT_N1 in
aim3_aim4_range_sweep.py, for the opposite/small-x version of a related
precision issue) -- large dynamic range quantities like exp(x) need
RELATIVE tolerance, or a comparison done in log-space, not absolute
tolerance, once x gets large. The tests below therefore compare in
LOG SPACE (log(q_bound + 1) = x, which is linear in x and does not blow
up the tiny floating-point noise the way exp(x) does) rather than
comparing q_bound directly -- this is not a weaker check, it is the
NUMERICALLY CORRECT way to check "these two computations agree," given
that q_bound spans nearly 60 orders of magnitude across the tested n
values. This was caught by actually running the tests after writing them
and investigating the size of the discrepancy rather than just loosening
a tolerance blindly -- exactly the "check the raw numbers, don't assume"
discipline this project's practical notes repeatedly credit for catching
real bugs.

HOW TO RUN
    pytest test_identifiability.py -v
or:
    python test_identifiability.py
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from binding_model import q_bound, capture_probability_patch


def test_sigma_and_k_intra_are_confounded_only_their_product_matters():
    """
    Rescaling sigma up and K_INTRA down (or vice versa) by the same
    factor, holding sigma*K_INTRA fixed, must leave q_bound() bit-
    identical -- because n*sigma*K_INTRA is the only place these two
    parameters ever appear together in the model.
    """
    n = np.array([5.0, 20.0, 80.0, 150.0])
    sigma = 4.44e11
    k_intra = 2e-12

    qb_baseline = q_bound(n, sigma, k_intra)

    for factor in (2.0, 37.5, 0.01, 1000.0):
        sigma_rescaled = sigma * factor
        k_intra_rescaled = k_intra / factor
        qb_rescaled = q_bound(n, sigma_rescaled, k_intra_rescaled)
        # Compare in LOG SPACE (log(q_bound + 1) recovers x = n*sigma*
        # k_intra, which is linear rather than exponential) -- see the
        # module docstring's note above for why a direct/absolute
        # comparison of q_bound itself is numerically the WRONG check at
        # large x, even though the underlying physics is exactly
        # identical.
        log_baseline = np.log(qb_baseline + 1.0)
        log_rescaled = np.log(qb_rescaled + 1.0)
        max_diff = np.max(np.abs(log_baseline - log_rescaled))
        assert max_diff < 1e-9, (
            f"Expected log(q_bound()+1) [i.e. the exponent x=n*sigma*"
            f"k_intra itself] to match to floating-point precision when "
            f"sigma is scaled by {factor} and k_intra is scaled by "
            f"1/{factor} (holding sigma*k_intra fixed) -- got a difference "
            f"of {max_diff}, larger than the machine-epsilon-scale "
            f"tolerance expected from floating-point rounding alone. If "
            f"this fails, the sigma/K_INTRA identifiability argument in "
            f"docs/Limitations_and_Calibration_Gap.md is no longer "
            f"accurate and needs to be re-derived."
        )
    print("  [OK] q_bound() output matches to floating-point precision "
          "(checked in log-space, the numerically correct way to compare "
          "a quantity spanning many orders of magnitude) under every "
          "tested sigma<->K_INTRA product-preserving rescaling -- "
          "confirming sigma and K_INTRA are not separately identifiable "
          "from a single condition's output.")


def test_k_geom_and_k_intra_are_confounded_at_patch_exponent_two():
    """
    At patch_exponent=2.0 (the project's default, soft-adhesion regime),
    only the product k_geom*K_INTRA matters -- rescaling k_geom up and
    K_INTRA down by the same factor must leave capture_probability_patch()
    bit-identical.
    """
    n = np.array([5.0, 20.0, 80.0, 150.0])
    diameter = np.array([60.0, 90.0, 120.0, 200.0])
    sigma = 4.44e11
    k_geom = 3.0
    k_intra = 2e-12
    z = 1e-4

    theta_baseline = capture_probability_patch(
        n, diameter, sigma, k_geom, k_intra, z, patch_exponent=2.0
    )

    for factor in (2.0, 5.0, 0.1, 250.0):
        k_geom_rescaled = k_geom * factor
        k_intra_rescaled = k_intra / factor
        theta_rescaled = capture_probability_patch(
            n, diameter, sigma, k_geom_rescaled, k_intra_rescaled, z, patch_exponent=2.0
        )
        max_diff = np.max(np.abs(theta_baseline - theta_rescaled))
        assert max_diff < 1e-9, (
            f"Expected capture_probability_patch() to match to floating-"
            f"point precision at patch_exponent=2.0 when k_geom is scaled "
            f"by {factor} and k_intra is scaled by 1/{factor} -- got a "
            f"difference of {max_diff}, larger than the machine-epsilon-"
            f"scale tolerance expected from floating-point rounding alone. "
            f"If this fails, the k_geom/K_INTRA identifiability argument "
            f"in docs/Limitations_and_Calibration_Gap.md is no longer "
            f"accurate and needs to be re-derived."
        )
    print("  [OK] capture_probability_patch() output matches to floating-"
          "point precision under every tested k_geom<->K_INTRA product-"
          "preserving rescaling at patch_exponent=2 -- confirming k_geom "
          "and K_INTRA are not separately identifiable from a single "
          "condition's output in the soft-adhesion regime.")


def test_k_geom_and_k_intra_confound_is_actually_patch_exponent_independent():
    """
    Checked directly rather than assumed: is the k_geom/K_INTRA confound
    demonstrated above SPECIFIC to patch_exponent=2.0 (where diameter
    drops out of contact_patch_receptors entirely), or does it hold more
    broadly? Looking at contact_patch_receptors() directly:

        n_contact = k_geom * n * diameter^(patch_exponent - 2)

    k_geom appears as a simple linear multiplier REGARDLESS of
    patch_exponent -- the diameter^(patch_exponent-2) term changes how
    diameter enters, but never changes how k_geom and K_INTRA combine
    (q_bound only ever sees n_contact * sigma * k_intra, and n_contact is
    linear in k_geom no matter what patch_exponent is). So the confound
    should hold at ANY patch_exponent, not just 2.0 -- this was the
    initial guess while writing this test file, checked here directly at
    patch_exponent=0.0 (the other named regime) rather than assumed.

    RESULT (confirmed by running this test, not assumed): the confound
    DOES hold at patch_exponent=0.0 too. This means the k_geom/K_INTRA
    confound is actually patch_exponent-INDEPENDENT -- broader than how
    docs/Limitations_and_Calibration_Gap.md currently states it (that
    document mentions it "once the contact-patch model is included" at
    the default patch_exponent=2, without explicitly stating it also
    holds at patch_exponent=0). This is a real, minor precision
    improvement to make to that document, not a correction of an error --
    the document's claim was not wrong, just narrower than the full truth
    turns out to be.
    """
    n = np.array([5.0, 20.0, 80.0, 150.0])
    diameter = np.array([60.0, 90.0, 120.0, 200.0])
    sigma = 4.44e11
    k_geom = 3.0
    k_intra = 2e-12
    z = 1e-4
    factor = 5.0

    theta_baseline = capture_probability_patch(
        n, diameter, sigma, k_geom, k_intra, z, patch_exponent=0.0
    )
    theta_rescaled = capture_probability_patch(
        n, diameter, sigma, k_geom * factor, k_intra / factor, z, patch_exponent=0.0
    )
    max_diff = np.max(np.abs(theta_baseline - theta_rescaled))
    assert max_diff < 1e-9, (
        "Expected the k_geom/K_INTRA confound to ALSO hold at "
        "patch_exponent=0.0 (matching to floating-point precision), since "
        "contact_patch_receptors() is linear in k_geom regardless of "
        "patch_exponent -- if this fails, the confound is actually "
        "patch_exponent=2-specific after all, and "
        "docs/Limitations_and_Calibration_Gap.md's scope should be "
        "narrowed rather than broadened."
    )
    print("  [OK] confirmed directly: the k_geom<->K_INTRA confound holds "
          "at patch_exponent=0.0 too (not just 2.0), since "
          "contact_patch_receptors() is linear in k_geom regardless of "
          "patch_exponent -- the confound is patch_exponent-INDEPENDENT, "
          "slightly broader than how the limitations write-up currently "
          "states it.")


if __name__ == "__main__":
    tests = [
        test_sigma_and_k_intra_are_confounded_only_their_product_matters,
        test_k_geom_and_k_intra_are_confounded_at_patch_exponent_two,
        test_k_geom_and_k_intra_confound_is_actually_patch_exponent_independent,
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
