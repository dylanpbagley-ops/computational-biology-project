"""
test_aim1.py

Aim 2's unit tests: the two "limiting cases" the handoff requires before
this model can be trusted for anything downstream.

WHAT IS A "LIMITING CASE" AND WHY DOES IT MATTER?
----------------------------------------------------
A limiting case is a special, simplified version of a general model where
we ALREADY KNOW (from independent, well-established theory) what the
correct answer should be. If our general model doesn't reduce to the known
correct answer in that special case, something is wrong with the general
model -- a bug, a sign error, a misremembered formula, etc. This is one of
the most reliable ways to catch mistakes in a physics/math model, because
it doesn't depend on having real experimental data (which we don't have
much of yet) -- it only depends on internal mathematical consistency.

THE TWO LIMITING CASES REQUIRED (from the handoff)
-----------------------------------------------------
1. Poisson-distributed receptor counts -> mean-field theory must be exact.
   As explained in population_model.mean_field_qb(), this is exact
   specifically at the level of the partition function q_b, not at the
   final capture-probability level (see that function's docstring for the
   full explanation of why). This test checks the part that actually is
   provably exact.

2. Valency = 1 (n=1 receptor) -> must recover the ordinary single-site
   Langmuir binding isotherm. We test this by comparing
   binding_model.capture_probability(n=1, ...) against
   binding_model.langmuir_probability(...), which was written completely
   independently (from the textbook Langmuir formula) so this is a
   genuine cross-check, not comparing a formula to itself.

   IMPORTANT: while writing this test, it turned out the two do NOT match
   exactly at every ligand density -- only in the dilute / weak-binding
   limit (small sigma * k_intra). See the long comment in
   binding_model.py (in the module docstring, under "WHY WE ALSO NEED THE
   n=1 (LANGMUIR) CASE") for the full explanation: this is because our
   q_bound() formula is exp(n*sigma*k_intra) - 1, which only approximates
   the linear Langmuir form (k_intra*sigma) for small sigma*k_intra (first
   term of a Taylor series). This is a known, understood property of the
   chosen physical regime, not a bug -- so the test below checks for
   convergence to Langmuir as density decreases, rather than exact
   equality everywhere.

HOW TO RUN THIS FILE
-----------------------
From the command line, from inside the tests/ folder (with src/ next to
it), run:
    pytest test_aim1.py -v
or, if pytest isn't installed:
    python test_aim1.py
(the __main__ block at the bottom runs all tests manually and prints
PASS/FAIL, so this file works even without pytest installed.)
"""

import sys
import os
import numpy as np

# Make src/ importable from this tests/ folder without needing to install
# anything as a package -- this just tells Python "also look in this other
# folder for modules to import."
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from binding_model import capture_probability, langmuir_probability, q_bound
from population_model import mean_field_qb, poisson_distribution, population_yield


# ---------------------------------------------------------------------
# Test 1: Poisson mean-field exactness (at the q_b level -- see the
# module docstring above and population_model.mean_field_qb for why this
# specific level, not the full theta level, is the one that is exactly
# true).
# ---------------------------------------------------------------------

def test_poisson_mean_field_exact_at_qb_level():
    """
    E[q_b(n)] computed by directly summing over a Poisson(mean_n)
    distribution must match the closed-form mean_field_qb() formula to
    essentially machine precision (i.e., differences only from
    floating-point rounding, not from any real approximation).
    """
    test_cases = [
        # (mean_n, sigma, k_intra)
        (5.0, 1.0, 0.02),
        (20.0, 1.0, 0.05),
        (50.0, 2.0, 0.01),
        (100.0, 0.5, 0.03),
    ]

    for mean_n, sigma, k_intra in test_cases:
        lam = sigma * k_intra
        # n_max needs to comfortably cover the Poisson distribution's
        # support; mean + 10 standard deviations is generous for a Poisson
        # distribution (std dev = sqrt(mean)).
        n_max = int(mean_n + 10 * np.sqrt(mean_n) + 20)

        n_values, probabilities = poisson_distribution(mean_n, n_max)
        qb_values = q_bound(n_values, sigma, k_intra)
        E_qb_direct_sum = np.sum(probabilities * qb_values)

        E_qb_closed_form = mean_field_qb(mean_n, sigma, k_intra)

        diff = abs(E_qb_direct_sum - E_qb_closed_form)
        assert diff < 1e-6, (
            f"Poisson mean-field theorem failed for mean_n={mean_n}, "
            f"sigma={sigma}, k_intra={k_intra}: direct sum gives "
            f"{E_qb_direct_sum}, closed form gives {E_qb_closed_form}, "
            f"difference={diff}"
        )
        print(
            f"  [OK] mean_n={mean_n:>6}, sigma={sigma}, k_intra={k_intra}: "
            f"direct sum={E_qb_direct_sum:.10f}, closed form={E_qb_closed_form:.10f}, "
            f"diff={diff:.2e}"
        )


def test_poisson_mean_field_theta_level_is_only_approximate():
    """
    This is a companion test to the one above -- it does NOT assert
    exactness, on purpose. It documents (and checks) that the FULL
    capture-probability version (after the Langmuir wrapper) is only
    approximately equal between "average theta over the population" and
    "theta evaluated at the mean-field q_b" -- and that the size of that
    approximation error stays within an expected, bounded range for
    reasonable parameters. If this test starts failing because the
    difference got MUCH bigger than before, that's worth a fresh look;
    if it fails because the difference became exactly zero, that would
    actually be suspicious too (it would mean something about the model
    changed in a way that accidentally made the two levels coincide).
    """
    z = 1.0
    test_cases = [
        (20.0, 1.0, 0.05),
        (20.0, 1.0, 0.1),
        (50.0, 0.5, 0.05),
    ]

    for mean_n, sigma, k_intra in test_cases:
        n_max = int(mean_n + 10 * np.sqrt(mean_n) + 20)
        n_values, probabilities = poisson_distribution(mean_n, n_max)

        # True population average of theta, computed directly.
        E_theta = population_yield(n_values, probabilities, sigma, k_intra, z)

        # Mean-field approximation: apply the Langmuir wrapper to the
        # exact mean-field q_b, rather than averaging theta directly.
        qb_mf = mean_field_qb(mean_n, sigma, k_intra)
        theta_mean_field = z * qb_mf / (1 + z * qb_mf)

        relative_diff = abs(E_theta - theta_mean_field) / E_theta

        # We expect this approximation to be in the same ballpark we found
        # during development (roughly 1-5%), not exact (0%) and not wildly
        # off (e.g. >20%, which would suggest a real bug rather than the
        # known Langmuir-step nonlinearity).
        assert 0 < relative_diff < 0.20, (
            f"Expected a small but non-zero approximation gap for "
            f"mean_n={mean_n}, sigma={sigma}, k_intra={k_intra}, got "
            f"relative_diff={relative_diff:.4f} (E_theta={E_theta:.6f}, "
            f"theta_mean_field={theta_mean_field:.6f})"
        )
        print(
            f"  [OK] mean_n={mean_n}, sigma={sigma}, k_intra={k_intra}: "
            f"E[theta]={E_theta:.6f}, theta(mean-field q_b)={theta_mean_field:.6f}, "
            f"relative gap={relative_diff:.2%} (expected, not a bug)"
        )


# ---------------------------------------------------------------------
# Test 2: valency = 1 must recover ordinary Langmuir binding.
# ---------------------------------------------------------------------

def test_valency_one_recovers_langmuir_in_dilute_limit():
    """
    With n=1 receptor, capture_probability() must converge to
    langmuir_probability() (an independently-written textbook Langmuir
    formula) as ligand density gets small (the "dilute limit").

    As explained above and in binding_model.py's module docstring, the two
    formulas are NOT identical at every density -- q_bound() uses
    exp(sigma*k_intra) - 1 while Langmuir uses the linear sigma*k_intra,
    and these only agree closely when sigma*k_intra is small. So this test
    checks two things instead of one blanket "must match" assertion:
        (a) at small sigma*k_intra, the relative gap between the two must
            be small (a strict, tight tolerance);
        (b) the relative gap must shrink as sigma*k_intra shrinks further
            -- confirming this is a controlled, well-behaved approximation
            rather than an arbitrary discrepancy.
    """
    k_intra_values = [0.01, 0.05, 0.2]
    z_values = [0.1, 1.0, 10.0]

    # A sequence of shrinking sigma values -- each 10x smaller than the
    # last -- to check that the gap shrinks correspondingly.
    shrinking_sigma = np.array([0.1, 0.01, 0.001, 0.0001])

    for k_intra in k_intra_values:
        for z in z_values:
            theta_general = capture_probability(1, shrinking_sigma, k_intra, z)
            theta_langmuir = langmuir_probability(shrinking_sigma, k_intra, z)
            relative_gap = np.abs(theta_general - theta_langmuir) / theta_langmuir

            # (a) smallest sigma tested should already be a very tight match
            assert relative_gap[-1] < 1e-4, (
                f"Expected near-exact Langmuir recovery at very low density "
                f"for k_intra={k_intra}, z={z}, but relative gap was "
                f"{relative_gap[-1]:.2e}"
            )
            # (b) gap should shrink monotonically as sigma shrinks (allowing
            # for tiny floating point wobble at the smallest scales)
            assert np.all(np.diff(relative_gap) <= 1e-12), (
                f"Expected the Langmuir approximation gap to shrink "
                f"monotonically as sigma decreases for k_intra={k_intra}, "
                f"z={z}, but it did not: relative gaps were {relative_gap}"
            )

    print("  [OK] n=1 converges to independent Langmuir formula in the "
          "dilute limit, with the gap shrinking as expected (not exact at "
          "all densities -- see binding_model.py docstring for why).")


# ---------------------------------------------------------------------
# Test 3 (bonus / sanity check): capture probability must be a
# monotonically increasing function of n. This isn't one of the two
# required limiting cases from the handoff, but it's a basic physical
# sanity check that's cheap to test and would catch a wide class of bugs
# (e.g. an accidental sign flip).
# ---------------------------------------------------------------------

def test_capture_probability_increases_with_receptor_count():
    """More receptors should never decrease capture probability."""
    sigma, k_intra, z = 1.0, 0.05, 1.0
    n_values = np.arange(0, 100)
    theta_values = capture_probability(n_values, sigma, k_intra, z)

    differences = np.diff(theta_values)
    assert np.all(differences >= 0), (
        "capture_probability is not monotonically increasing in n -- "
        "found a decrease somewhere, which should be physically impossible."
    )
    # Also check the n=0 boundary case explicitly: zero receptors should
    # mean zero capture probability.
    assert theta_values[0] == 0.0, (
        f"Expected theta=0 for n=0 receptors, got {theta_values[0]}"
    )
    print("  [OK] capture probability is monotonically increasing in n, "
          "and theta(n=0)=0 exactly.")




# ---------------------------------------------------------------------
# Test 4: contact-patch-aware model (n vs. areal-density resolution).
# Added to resolve "Immediate next steps" item 1 from the v3 handoff.
# ---------------------------------------------------------------------

def test_patch_exponent_2_reproduces_original_model_exactly():
    """
    q_bound_patch()/capture_probability_patch() at patch_exponent=2.0 (the
    soft-adhesion default) must be EXACTLY identical to the original,
    size-blind q_bound()/capture_probability(), with k_geom*k_intra playing
    the role of the original model's single k_intra. This is the backward-
    compatibility guarantee: adding the diameter-aware machinery must not
    change any previously-validated result (analytical vs. Monte Carlo,
    Poisson mean-field, Langmuir limiting case) unless patch_exponent=0 is
    explicitly requested.
    """
    from binding_model import q_bound_patch, capture_probability_patch

    n = np.array([1.0, 5.0, 20.0, 50.0, 100.0])
    # Diameters are arbitrary and DELIBERATELY varied -- if patch_exponent=2
    # is working correctly, the result must not depend on these at all.
    diameters = np.array([30.0, 55.0, 90.0, 150.0, 300.0])
    sigma, k_intra, k_geom, z = 1.0, 0.05, 0.03, 1.0

    qb_original = q_bound(n, sigma, k_geom * k_intra)
    qb_patch = q_bound_patch(n, diameters, sigma, k_geom, k_intra, patch_exponent=2.0)
    assert np.allclose(qb_original, qb_patch, atol=1e-12), (
        f"q_bound_patch at patch_exponent=2 should exactly match q_bound "
        f"with k_geom*k_intra, but got original={qb_original}, patch={qb_patch}"
    )

    theta_original = capture_probability(n, sigma, k_geom * k_intra, z)
    theta_patch = capture_probability_patch(n, diameters, sigma, k_geom, k_intra, z, patch_exponent=2.0)
    assert np.allclose(theta_original, theta_patch, atol=1e-12), (
        "capture_probability_patch at patch_exponent=2 should exactly "
        "match capture_probability with k_geom*k_intra."
    )

    print("  [OK] patch_exponent=2.0 reproduces the original size-blind "
          "model exactly (diameter-independent, as required for backward "
          "compatibility).")


def test_patch_exponent_0_depends_only_on_areal_density():
    """
    At patch_exponent=0.0 (fixed-patch / density-driven regime), two
    vesicles with the SAME areal density (n / (pi * D^2)) but DIFFERENT
    individual n and D must give the SAME q_b -- confirming this regime
    genuinely reduces to a pure density dependence, not some intermediate
    mix of n and D.
    """
    from binding_model import contact_patch_receptors, q_bound_patch

    # Same density (n/D^2 ratio held fixed), different absolute n and D:
    n_pair = np.array([20.0, 80.0])    # 4x more receptors...
    d_pair = np.array([55.0, 110.0])   # ...but 2x diameter -> 4x more area
    sigma, k_intra, k_geom = 1.0, 0.05, 0.03

    density = n_pair / (np.pi * d_pair ** 2)
    assert np.isclose(density[0], density[1], rtol=1e-10), (
        "Test setup error: these two (n, D) pairs should have equal "
        "areal density by construction."
    )

    qb = q_bound_patch(n_pair, d_pair, sigma, k_geom, k_intra, patch_exponent=0.0)
    assert np.isclose(qb[0], qb[1], rtol=1e-10), (
        f"At patch_exponent=0, equal-density vesicles should give equal "
        f"q_b regardless of absolute size, but got {qb}"
    )
    print("  [OK] patch_exponent=0.0 depends only on areal density, "
          "confirmed with two same-density, different-size vesicles "
          "giving identical q_b.")


def test_joint_population_yield_matches_manual_mean():
    """
    joint_population_yield() should be exactly the plain mean of
    capture_probability_patch() over the sample -- no hidden weighting or
    normalization. Checked against a manual calculation.
    """
    from binding_model import capture_probability_patch
    from population_model import joint_population_yield

    rng = np.random.default_rng(seed=7)
    n_sample = rng.uniform(5, 100, size=50)
    d_sample = rng.uniform(30, 200, size=50)
    sigma, k_geom, k_intra, z = 1.0, 0.03, 0.05, 1.0

    for p in (2.0, 0.0):
        expected = np.mean(capture_probability_patch(n_sample, d_sample, sigma, k_geom, k_intra, z, patch_exponent=p))
        got = joint_population_yield(n_sample, d_sample, sigma, k_geom, k_intra, z, patch_exponent=p)
        assert np.isclose(expected, got, rtol=1e-12), (
            f"joint_population_yield (patch_exponent={p}) should equal a "
            f"plain mean of per-vesicle capture probability; expected "
            f"{expected}, got {got}"
        )
    print("  [OK] joint_population_yield matches a manual mean calculation "
          "at both patch_exponent=2 and patch_exponent=0.")


def test_lennon_ranking_flips_between_patch_regimes():
    """
    THE HEADLINE CHECK for 'Immediate next steps' item 1: running the
    diameter-aware model on the actual digitized Lennon Figure 3c
    healthy-vs-PDAC populations, the PDAC/healthy yield ratio should be
    ABOVE 1 (PDAC captures more) under patch_exponent=2 (size-blind,
    matches the original model) and BELOW 1 (PDAC captures less) under
    patch_exponent=0 (density-driven) -- i.e., the qualitative ranking
    genuinely flips depending on this modeling choice. This is not a
    "should be true" assertion in the normal sense -- it's confirming
    that the sensitivity finding reported in the v4 handoff is real and
    reproducible from the data files, not a one-off script output that
    could silently stop being true if the data or code changes.

    Uses a separately-calibrated k_geom per regime (see
    patch_regime_sensitivity.py for why a single shared k_geom is not a
    fair comparison -- it has different physical units in each regime).
    """
    import csv
    import os

    data_dir = os.path.join(os.path.dirname(__file__), "..", "..", "ev_capture_data")
    if not os.path.isdir(data_dir):
        print("  [SKIPPED] ev_capture_data/ not found relative to this test "
              "file -- this check requires the digitized Lennon data, which "
              "lives in the workspace folder, not the code repo. Not a "
              "failure; just not runnable in isolation.")
        return

    def load_pairs(fn):
        diam, copies = [], []
        with open(os.path.join(data_dir, fn)) as f:
            r = csv.DictReader(f)
            for row in r:
                d = float(row["diameter_nm"])
                c = float(row["copies_per_ev"])
                if d > 0:
                    diam.append(d)
                    copies.append(c)
        return np.array(copies), np.array(diam)

    from population_model import joint_population_yield

    n_h, d_h = load_pairs("fig3c_EGFR_healthy_digitized.csv")
    n_p, d_p = load_pairs("fig3c_EGFR_PDAC_digitized.csv")

    sigma, k_intra, z = 1.0, 0.05, 1.0
    ratios = {}
    for p in (2.0, 0.0):
        k_geom = 1.0 / (n_h.mean() * d_h.mean() ** (p - 2.0))
        y_h = joint_population_yield(n_h, d_h, sigma, k_geom, k_intra, z, patch_exponent=p)
        y_p = joint_population_yield(n_p, d_p, sigma, k_geom, k_intra, z, patch_exponent=p)
        ratios[p] = y_p / y_h

    assert ratios[2.0] > 1.0, (
        f"Expected PDAC/healthy ratio > 1 under patch_exponent=2 (size-"
        f"blind), got {ratios[2.0]:.3f}"
    )
    assert ratios[0.0] < 1.0, (
        f"Expected PDAC/healthy ratio < 1 under patch_exponent=0 "
        f"(density-driven), got {ratios[0.0]:.3f}"
    )
    print(f"  [OK] Ranking genuinely flips on real Lennon EGFR data: "
          f"patch_exponent=2 gives PDAC/healthy={ratios[2.0]:.2f}x (more), "
          f"patch_exponent=0 gives PDAC/healthy={ratios[0.0]:.2f}x (less).")
if __name__ == "__main__":
    tests = [
        test_poisson_mean_field_exact_at_qb_level,
        test_poisson_mean_field_theta_level_is_only_approximate,
        test_valency_one_recovers_langmuir_in_dilute_limit,
        test_capture_probability_increases_with_receptor_count,
        test_patch_exponent_2_reproduces_original_model_exactly,
        test_patch_exponent_0_depends_only_on_areal_density,
        test_joint_population_yield_matches_manual_mean,
        test_lennon_ranking_flips_between_patch_regimes,
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
