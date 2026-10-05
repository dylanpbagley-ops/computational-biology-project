"""
test_k_intra_z_plausibility_bounds.py

Regression test for src/k_intra_z_plausibility_bounds.py (see that file's
module docstring for the full context: this is the documented FALLBACK
used after searching for, and not finding, a clean single-marker absolute
calibration source for K_INTRA/z/k_geom this session -- not a fit, not a
calibration, just a checked plausibility range and a record of the real
numbers found). This test exists so the specific numbers reported stay
checked rather than becoming stale if the module changes later, and so
the internal consistency of the dimensionless-exponent table is verified
(target = n * sigma_K, by construction).

HOW TO RUN
    pytest test_k_intra_z_plausibility_bounds.py -v
or:
    python test_k_intra_z_plausibility_bounds.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from k_intra_z_plausibility_bounds import (
    dimensionless_exponent_bounds,
    concentration_axis_cross_check,
    k_intra_bounds_from_real_sigma,
    LOPEZ_BALTAZAR_KD_M,
    GENERAL_ANTIBODY_KD_RANGE_M,
    MIZENKO_DILUTION_RANGE_PARTICLES_PER_ML,
    KAKAVOULIA_EC50_ENGINEERED_PARTICLES_PER_ML,
    KAKAVOULIA_EC50_NATURAL_PARTICLES_PER_ML,
    SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2,
    SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
    LOPEZ_BALTAZAR_COVERAGE_FRACTION_AT_KD,
    LOPEZ_BALTAZAR_CONTROL_CONCENTRATIONS_PARTICLES_PER_ML,
)


def test_dimensionless_exponent_table_is_internally_consistent():
    """
    For every (n, target_exponent) pair, the returned sigma*K_INTRA value
    must satisfy n * sigma_K == target_exponent exactly (this is how the
    table is constructed -- solving n*sigma*K_INTRA = target for
    sigma*K_INTRA). This is a pure consistency check, not a physical claim:
    it catches a broken formula, not a wrong physical assumption.
    """
    bounds = dimensionless_exponent_bounds()
    for (n, target), sigma_k in bounds.items():
        recovered = n * sigma_k
        assert abs(recovered - target) < 1e-9, (
            f"n={n}, target={target}: n*sigma_K = {recovered}, expected "
            f"{target} -- the table's own defining relationship is broken."
        )
    print(f"  [OK] all {len(bounds)} (n, target_exponent) entries in the "
          f"plausibility table satisfy n*sigma_K == target exactly.")


def test_dimensionless_exponent_bounds_span_expected_range():
    """
    Sanity check that the honest width of the bound is what the module
    docstring claims: roughly four orders of magnitude between the
    smallest sigma*K_INTRA (large n, small target) and the largest
    (small n, large target). This test exists so a future edit that
    quietly narrows the claimed range (making it look falsely more
    precise) gets caught.
    """
    bounds = dimensionless_exponent_bounds()
    values = list(bounds.values())
    ratio = max(values) / min(values)
    assert ratio > 1000, (
        f"Expected the plausibility bound to span at least 3 orders of "
        f"magnitude (matching the module's honest 'this is wide, not a "
        f"fit' framing), got a max/min ratio of {ratio:.1f} -- if this "
        f"has narrowed a lot, check whether the module docstring's "
        f"'four orders of magnitude' claim still matches."
    )
    print(f"  [OK] plausibility bound spans a max/min ratio of {ratio:.0f}x, "
          f"confirming this is a wide honest range, not a tight fit "
          f"dressed up as one.")


def test_lopez_baltazar_kd_is_within_general_antibody_range():
    """
    The one real EV-specific Kd found this session (Lopez Baltazar,
    anti-CD81, SPR) should fall within the broader general antibody-antigen
    Kd range also cited -- if it didn't, that would be worth flagging as
    an inconsistency between the specific and general anchors used in the
    module docstring's reasoning.
    """
    lo, hi = GENERAL_ANTIBODY_KD_RANGE_M
    assert lo <= LOPEZ_BALTAZAR_KD_M <= hi, (
        f"Lopez Baltazar's Kd ({LOPEZ_BALTAZAR_KD_M:.2e} M) falls outside "
        f"the general antibody-antigen range ({lo:.2e}-{hi:.2e} M) cited "
        f"alongside it -- this specific anchor and the general range "
        f"should agree, or the docstring's framing needs revisiting."
    )
    print(f"  [OK] Lopez Baltazar's Kd ({LOPEZ_BALTAZAR_KD_M*1e9:.2f} nM) "
          f"sits within the general antibody-antigen range "
          f"({lo*1e9:.1f}-{hi*1e9:.1f} nM), confirming the specific "
          f"EV anchor is consistent with the broader literature range.")


def test_concentration_cross_check_ranges_overlap():
    """
    The module's headline concentration-axis finding is that Mizenko's own
    dilution range and Kakavoulia's EC50 values converge on the same order
    of magnitude (~1e8-1e9 particles/mL). This test checks that claim
    directly: Kakavoulia's EC50 values (both engineered and natural) must
    fall within (or very close to) Mizenko's tested dilution range, since
    that overlap is the actual evidence for "convergence" being asserted
    in the docstring -- if a future edit changes either number such that
    they no longer overlap, this claim would become false and should be
    caught here rather than silently going stale.
    """
    mizenko_lo, mizenko_hi = MIZENKO_DILUTION_RANGE_PARTICLES_PER_ML
    kakavoulia_lo, kakavoulia_hi = KAKAVOULIA_EC50_NATURAL_PARTICLES_PER_ML

    assert mizenko_lo <= KAKAVOULIA_EC50_ENGINEERED_PARTICLES_PER_ML <= mizenko_hi, (
        f"Kakavoulia's engineered-EV EC50 "
        f"({KAKAVOULIA_EC50_ENGINEERED_PARTICLES_PER_ML:.2e} particles/mL) "
        f"falls outside Mizenko's tested dilution range "
        f"({mizenko_lo:.2e}-{mizenko_hi:.2e} particles/mL) -- the "
        f"'convergent order of magnitude' claim in the docstring would "
        f"no longer be supported."
    )
    assert mizenko_lo <= kakavoulia_hi and kakavoulia_lo <= mizenko_hi, (
        f"Kakavoulia's natural-EV EC50 range "
        f"({kakavoulia_lo:.2e}-{kakavoulia_hi:.2e} particles/mL) does not "
        f"overlap Mizenko's tested dilution range "
        f"({mizenko_lo:.2e}-{mizenko_hi:.2e} particles/mL) at all -- the "
        f"'convergent order of magnitude' claim in the docstring would "
        f"no longer be supported."
    )
    print(f"  [OK] Kakavoulia's EC50 values (engineered="
          f"{KAKAVOULIA_EC50_ENGINEERED_PARTICLES_PER_ML:.1e}, natural="
          f"{kakavoulia_lo:.1e}-{kakavoulia_hi:.1e} particles/mL) overlap "
          f"Mizenko's own tested dilution range "
          f"({mizenko_lo:.1e}-{mizenko_hi:.1e} particles/mL), confirming "
          f"the cross-platform convergence claim.")


def test_concentration_axis_cross_check_returns_expected_keys():
    """
    Structural check that concentration_axis_cross_check() returns the
    fields the module docstring promises, so a future refactor can't
    silently drop one of the cited numbers without a test noticing.
    """
    result = concentration_axis_cross_check()
    expected_keys = {
        "mizenko_dilution_range_particles_per_mL",
        "mizenko_practical_ceiling_particles_per_mL",
        "kakavoulia_ec50_engineered_particles_per_mL",
        "kakavoulia_ec50_natural_range_particles_per_mL",
        "convergent_order_of_magnitude",
        "caveat",
    }
    assert expected_keys.issubset(result.keys()), (
        f"Missing expected keys: {expected_keys - result.keys()}"
    )
    print(f"  [OK] concentration_axis_cross_check() returns all "
          f"{len(expected_keys)} expected fields.")


def test_sigma_at_kd_equals_coverage_fraction_times_sigma_max():
    """
    NEW (v7->v8): SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2 must equal
    LOPEZ_BALTAZAR_COVERAGE_FRACTION_AT_KD * SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2
    exactly -- this is how it's defined (40% of the real max close-packed
    density found in the dissertation's Appendix A.1), and a future edit
    that breaks this relationship would silently make the "real sigma"
    claim wrong.
    """
    expected = LOPEZ_BALTAZAR_COVERAGE_FRACTION_AT_KD * SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2
    assert abs(SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2 - expected) < 1e6, (
        f"SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2 ({SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2:.3e}) "
        f"does not equal coverage_fraction * sigma_max ({expected:.3e})."
    )
    print(f"  [OK] sigma at Kd ({SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2:.2e} Abs/cm^2) = "
          f"{LOPEZ_BALTAZAR_COVERAGE_FRACTION_AT_KD:.0%} of sigma_max "
          f"({SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2:.2e} Abs/cm^2), as reported "
          f"in the dissertation.")


def test_k_intra_bounds_from_real_sigma_matches_dimensionless_bound_divided_by_sigma():
    """
    NEW (v7->v8): k_intra_bounds_from_real_sigma() should be exactly the
    dimensionless_exponent_bounds() table divided by the real sigma value
    -- this is a pure unit-conversion function, not a new physical claim,
    so it should reproduce that arithmetic exactly for every (n, target)
    key.
    """
    dimensionless = dimensionless_exponent_bounds()
    k_intra_bounds = k_intra_bounds_from_real_sigma()
    assert set(dimensionless.keys()) == set(k_intra_bounds.keys()), (
        "k_intra_bounds_from_real_sigma() should cover the same (n, target) "
        "keys as dimensionless_exponent_bounds() by default."
    )
    for key, sigma_k in dimensionless.items():
        expected_k_intra = sigma_k / SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2
        actual_k_intra = k_intra_bounds[key]
        rel_err = abs(actual_k_intra - expected_k_intra) / expected_k_intra
        assert rel_err < 1e-9, (
            f"key={key}: expected K_INTRA={expected_k_intra:.3e}, got "
            f"{actual_k_intra:.3e} (relative error {rel_err:.2e})."
        )
    print(f"  [OK] all {len(k_intra_bounds)} K_INTRA bound entries equal "
          f"the dimensionless bound divided by the real sigma "
          f"({SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2:.2e} Abs/cm^2), confirming "
          f"this is a pure unit conversion, not a new fit.")


def test_k_intra_bounds_still_span_wide_range():
    """
    NEW (v7->v8): converting to real physical units (cm^2/Ab) via a real
    sigma must NOT narrow the honest width of the bound -- dividing every
    entry by the same constant preserves the max/min ratio exactly. This
    test guards against a future edit accidentally presenting the
    real-sigma version as more precise than the dimensionless version it
    comes from.
    """
    k_intra_bounds = k_intra_bounds_from_real_sigma()
    values = list(k_intra_bounds.values())
    ratio = max(values) / min(values)
    assert ratio > 1000, (
        f"Expected the K_INTRA bound (like the dimensionless bound it's "
        f"derived from) to span at least 3 orders of magnitude, got a "
        f"max/min ratio of {ratio:.1f} -- if this has narrowed, check "
        f"whether k_intra_bounds_from_real_sigma() is still a pure "
        f"division by a constant sigma."
    )
    print(f"  [OK] K_INTRA bound (real sigma, physical units) still spans "
          f"a max/min ratio of {ratio:.0f}x, confirming the real-sigma "
          f"conversion did not narrow the honest width of the bound.")


def test_lopez_baltazar_control_concentrations_are_above_linear_range():
    """
    NEW (v7->v8): the control-experiment concentrations retrieved from
    Appendix A of the dissertation (used for nonspecific-binding controls,
    not the headline titration) are all reported as ABOVE the linear
    dynamic range's upper bound -- this matches the dissertation's own
    framing (these are higher-concentration checks, not points on the
    working titration curve). If a future edit introduces a value that
    contradicts this, it's worth checking against the source again.
    """
    from k_intra_z_plausibility_bounds import LOPEZ_BALTAZAR_LINEAR_RANGE_PARTICLES_PER_ML
    _, linear_hi = LOPEZ_BALTAZAR_LINEAR_RANGE_PARTICLES_PER_ML
    for conc in LOPEZ_BALTAZAR_CONTROL_CONCENTRATIONS_PARTICLES_PER_ML:
        assert conc >= linear_hi, (
            f"Control concentration {conc:.2e} particles/mL is below the "
            f"linear range's upper bound ({linear_hi:.2e}) -- expected "
            f"control experiments to probe concentrations at or above the "
            f"working linear range, per the dissertation's Appendix A."
        )
    print(f"  [OK] all {len(LOPEZ_BALTAZAR_CONTROL_CONCENTRATIONS_PARTICLES_PER_ML)} "
          f"control concentrations are at or above the linear range's "
          f"upper bound ({linear_hi:.2e} particles/mL), consistent with "
          f"the dissertation's framing of these as higher-concentration "
          f"control checks.")


if __name__ == "__main__":
    tests = [
        test_dimensionless_exponent_table_is_internally_consistent,
        test_dimensionless_exponent_bounds_span_expected_range,
        test_lopez_baltazar_kd_is_within_general_antibody_range,
        test_concentration_cross_check_ranges_overlap,
        test_concentration_axis_cross_check_returns_expected_keys,
        test_sigma_at_kd_equals_coverage_fraction_times_sigma_max,
        test_k_intra_bounds_from_real_sigma_matches_dimensionless_bound_divided_by_sigma,
        test_k_intra_bounds_still_span_wide_range,
        test_lopez_baltazar_control_concentrations_are_above_linear_range,
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
