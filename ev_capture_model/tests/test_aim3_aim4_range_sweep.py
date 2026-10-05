"""
test_aim3_aim4_range_sweep.py

Regression tests for src/aim3_aim4_range_sweep.py -- the v7->v8 handoff
transition's range-sweep work (v7's "Immediate next steps" items 5 and 7:
sweep Aim 3/4 findings across the K_INTRA plausibility range rather than
using arbitrary placeholder values), UPDATED in the v8->v9 transition
after a confirmed bug in sweep_superselectivity_reachability() was found
and fixed (see that function's docstring in aim3_aim4_range_sweep.py, and
handoff v9's "What changed since v8").

WHAT CHANGED IN THIS FILE, AND WHY (v8->v9)
---------------------------------------------
v8's version of this file asserted, as a locked-in finding,
"any_superselective is False" and "max_alpha_found stays below the
superselectivity threshold." That finding does NOT reproduce: it was the
output of sweep_superselectivity_reachability() silently defaulting its
`sigma` parameter to the OLD placeholder value 1.0 instead of
SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, the real antibody surface density every
swept K_INTRA value was actually derived against. With the bug fixed
(sigma now correctly defaults to the real value), the corrected,
independently-reproduced finding is the OPPOSITE: alpha>1 IS reached,
across roughly 60% of the numerically-reliable (K_INTRA, z) grid, peaking
around alpha~10.46. The two tests that encoded the old (wrong) finding as
an assertion have been REPLACED below with tests that lock in the
CORRECTED behavior, plus a new explicit regression test that checks the
function's default argument directly -- so a future refactor cannot
silently reintroduce the sigma=1.0 default bug without a test failing
immediately, even before anyone inspects the numerical output.

This is deliberately the same discipline the original (v7->v8) version of
this file used for its own first two bugs: locking in the corrected
answer with a test, not just fixing the code and moving on.

These tests now lock in:

1. The healthy/PDAC-plasma (size-mismatched) ranking flip is robust
   across the ENTIRE plausibility range (ratio_p2>1 and ratio_p0<1 at
   every tested K_INTRA value). Unchanged from v8 -- independently
   re-confirmed in the v8->v9 transition, not just assumed still valid.
2. The PANC-1/HPDEC (size-matched) comparison does NOT show this large,
   clean flip -- it stays close to ratio~1 throughout, dipping only
   slightly below 1 at the single largest (most-saturating) K_INTRA
   value tested. Unchanged from v8 -- independently re-confirmed.
3. sweep_ranking_flip_robustness() defaults to the real sigma. Unchanged
   from v8 -- independently re-confirmed, not assumed correct by analogy
   just because a bug turned up in its sibling function.
4. NEW -- sweep_superselectivity_reachability() ALSO defaults to the
   real sigma (SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2), not the old placeholder
   1.0. This is the direct regression test for the v8->v9 bug fix.
5. NEW, CORRECTED -- alpha>1 (genuine superselectivity) IS reached in a
   substantial fraction (checked: roughly 50-70%, see tolerance below) of
   the numerically-reliable portion of the (K_INTRA, z) plausibility
   grid, with a peak alpha well above the 1.01 threshold used to
   distinguish a real signal from float64 noise. This REPLACES v8's
   "not reached anywhere" test, which was itself the bug's output.
6. NEW -- with the real sigma, few or zero grid cells need to be excluded
   as numerically unreliable (unlike the old buggy sigma=1.0 run, which
   excluded 91 of 104 cells). This is also new coverage: it distinguishes
   "the corrected run is trustworthy across nearly the whole grid" from
   "the corrected run just moved the numerical-artifact problem
   somewhere else," which was checked directly rather than assumed.

HOW TO RUN
    pytest test_aim3_aim4_range_sweep.py -v
or:
    python test_aim3_aim4_range_sweep.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from aim3_aim4_range_sweep import (
    sweep_ranking_flip_robustness,
    sweep_superselectivity_reachability,
    strict_resolved_egfr_ca199_comparison,
    ca199_selectivity_sensitivity_to_accessibility,
    get_k_intra_sweep_values,
    get_z_sweep_values,
    ALPHA_SUPERSELECTIVE_THRESHOLD,
)
from k_intra_z_plausibility_bounds import SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2


def test_healthy_pdac_plasma_flip_is_robust_across_full_k_intra_range():
    """
    The size-mismatched comparison must show ratio_p2>1 AND ratio_p0<1 at
    EVERY tested K_INTRA value -- this is the headline "robust flip"
    finding, independently re-confirmed in the v8->v9 transition.
    """
    results = sweep_ranking_flip_robustness()
    rows = results["healthy_pdac_plasma"]
    assert len(rows) >= 4, "Expected the sweep to cover multiple K_INTRA values."
    for r in rows:
        assert r["ratio_p2"] > 1.0, (
            f"k_intra={r['k_intra']:.3e}: expected ratio_p2>1 (size-blind "
            f"patch regime), got {r['ratio_p2']:.4f}."
        )
        assert r["ratio_p0"] < 1.0, (
            f"k_intra={r['k_intra']:.3e}: expected ratio_p0<1 (density-"
            f"driven patch regime), got {r['ratio_p0']:.4f}."
        )
    print(f"  [OK] healthy/PDAC-plasma ranking flip (ratio_p2>1, ratio_p0<1) "
          f"holds at all {len(rows)} tested K_INTRA values, confirming "
          f"this finding is robust across the plausibility range, not an "
          f"artifact of one arbitrary parameter choice.")


def test_panc1_hpdec_stays_near_one_not_a_clean_flip():
    """
    The size-matched comparison should NOT show the same kind of large,
    clean flip -- ratios should stay in a narrow band close to 1
    (checked: within roughly 0.99-1.5), unlike the mismatched comparison's
    ratio_p0 which drops as low as ~0.17-0.27. This distinguishes a real
    "roughly matched populations behave similarly" finding from a
    coincidence.
    """
    results = sweep_ranking_flip_robustness()
    rows = results["panc1_hpdec"]
    ratios_p0 = [r["ratio_p0"] for r in rows]
    ratios_p2 = [r["ratio_p2"] for r in rows]

    assert min(ratios_p0) > 0.9, (
        f"Expected the size-matched comparison's ratio_p0 to stay well "
        f"above the ~0.17-0.27 seen for the mismatched comparison (i.e. "
        f"no large flip) -- got a minimum of {min(ratios_p0):.4f}, which "
        f"would suggest a much larger flip than this session's checked "
        f"finding."
    )
    assert max(ratios_p2) < 2.0, (
        f"ratio_p2 for the size-matched comparison should stay in a "
        f"modest range (checked this session: up to ~1.48), not blow up -- "
        f"got a maximum of {max(ratios_p2):.4f}."
    )
    print(f"  [OK] PANC-1/HPDEC ratios stay in the checked narrow band "
          f"(ratio_p0 in [{min(ratios_p0):.3f}, {max(ratios_p0):.3f}], "
          f"ratio_p2 in [{min(ratios_p2):.3f}, {max(ratios_p2):.3f}]) -- "
          f"consistent with 'no large clean flip for this roughly "
          f"size-matched comparison,' not a coincidence of one parameter "
          f"choice.")


def test_ranking_flip_sweep_uses_real_sigma_not_placeholder():
    """
    Regression test for the unit-mismatch bug caught in the v7->v8
    session: the default sigma in sweep_ranking_flip_robustness() must be
    the REAL Lopez Baltazar sigma, not the old placeholder sigma=1.0 --
    because the K_INTRA values it's paired with were derived assuming
    that real sigma. This test checks the function signature's actual
    default value directly (not just its behavior), so a future refactor
    that reintroduces sigma=1.0 as the default is caught immediately.

    Independently re-confirmed in the v8->v9 transition (not assumed
    correct just because this function's docstring said so) after a bug
    of the identical class turned up in this function's sibling,
    sweep_superselectivity_reachability() -- see the next test.
    """
    import inspect
    sig = inspect.signature(sweep_ranking_flip_robustness)
    default_sigma = sig.parameters["sigma"].default
    assert default_sigma == SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, (
        f"sweep_ranking_flip_robustness()'s default sigma is "
        f"{default_sigma}, expected it to be "
        f"SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2 "
        f"({SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2:.3e}) -- using a different "
        f"sigma than the one the K_INTRA plausibility values were derived "
        f"against silently breaks the n*sigma*K_INTRA product, as caught "
        f"and fixed in the v7->v8 session."
    )
    print(f"  [OK] sweep_ranking_flip_robustness() defaults to the real "
          f"sigma ({SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2:.3e} Abs/cm^2), "
          f"matching the K_INTRA values it's paired with.")


def test_superselectivity_sweep_uses_real_sigma_not_placeholder():
    """
    NEW in v8->v9 -- direct regression test for the confirmed bug found
    in this session: sweep_superselectivity_reachability()'s default
    sigma previously was 1.0 (the old placeholder), NOT
    SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, despite being called with no
    arguments at its one call site in print_summary(). This produced
    v8's entire "superselectivity not reached anywhere" finding, which
    was consequently a numerical artifact, not real physics -- see
    handoff v9's "What changed since v8," item 1, for the full numeric
    confirmation (calling the function both ways with everything else
    identical: sigma=1.0 gives max_alpha=1.000013 / any_superselective=
    False / 13 of 104 cells reliable; the real sigma gives max_alpha=
    10.46 / any_superselective=True / all 104 cells reliable).

    This test checks the function signature's actual default value
    directly, exactly mirroring the equivalent test for
    sweep_ranking_flip_robustness() above -- so if this exact bug is ever
    reintroduced (e.g. by a careless refactor of the two sibling
    functions), it is caught immediately by a failing test, not by a
    user pressure-testing a suspiciously clean result in a later session.
    """
    import inspect
    sig = inspect.signature(sweep_superselectivity_reachability)
    default_sigma = sig.parameters["sigma"].default
    assert default_sigma == SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, (
        f"sweep_superselectivity_reachability()'s default sigma is "
        f"{default_sigma}, expected it to be "
        f"SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2 "
        f"({SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2:.3e}). Using sigma=1.0 (the "
        f"old placeholder) here silently invalidates the entire "
        f"superselectivity-reachability finding -- this is the exact bug "
        f"found and fixed in the v8->v9 transition. Do not reintroduce it."
    )
    print(f"  [OK] sweep_superselectivity_reachability() defaults to the "
          f"real sigma ({SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2:.3e} Abs/cm^2), "
          f"matching the K_INTRA values it's paired with -- the v8->v9 "
          f"bug fix is in place.")


def test_genuine_superselectivity_is_reached_in_current_plausibility_range():
    """
    CORRECTED finding, replacing v8's (bug-caused, retracted) opposite
    claim. With the sigma bug fixed, alpha>1 IS reached across a
    substantial fraction of the numerically-reliable (K_INTRA, z)
    plausibility grid -- independently re-derived in the v8->v9
    transition (ad hoc in a review chat: max_alpha~10.46, ~60% of grid)
    and re-confirmed here via the actual fixed function, not just copied
    from the review chat's interactive numbers.

    Tolerances are set loosely around the checked values (peak alpha
    comfortably above 2x the superselectivity threshold; superselective
    fraction between 30% and 90%) rather than pinned to the exact
    max_alpha=10.461... and fraction=59.6% found this session, so that
    small, expected floating-point differences across machines/numpy
    versions don't spuriously fail this test, while still catching a
    real regression (e.g. the bug being reintroduced, which would drop
    the fraction to 0% and max_alpha to ~1.0).
    """
    result = sweep_superselectivity_reachability()
    assert result["grid_size"] > 0, "Expected at least one numerically reliable grid cell."
    assert result["any_superselective"] is True, (
        "Expected any_superselective=True with the real sigma correctly "
        "applied -- if this is False, check whether the sigma-default "
        "bug (fixed in v8->v9) has been reintroduced before concluding "
        "superselectivity is genuinely absent."
    )
    assert result["max_alpha_found"] > 2 * ALPHA_SUPERSELECTIVE_THRESHOLD, (
        f"Expected max_alpha_found well above the superselectivity "
        f"threshold ({ALPHA_SUPERSELECTIVE_THRESHOLD}) -- this session's "
        f"corrected run found ~10.46. Got {result['max_alpha_found']:.4f}, "
        f"which is suspiciously close to the threshold and could indicate "
        f"a partial regression (e.g. the sigma bug reintroduced for only "
        f"part of the grid)."
    )
    assert 0.3 < result["fraction_of_grid_superselective"] < 0.9, (
        f"Expected roughly 30-90% of the reliable grid to be "
        f"superselective (this session's corrected run found ~59.6%) -- "
        f"got {result['fraction_of_grid_superselective']:.1%}, outside "
        f"the loose tolerance band. Investigate before trusting either "
        f"this run or the original 59.6% finding."
    )
    print(f"  [OK] genuine superselectivity (alpha>1) IS reached: max "
          f"alpha {result['max_alpha_found']:.3f}, in "
          f"{result['fraction_of_grid_superselective']:.1%} of the "
          f"reliable (K_INTRA, z) grid -- the corrected, bug-fixed "
          f"finding, replacing v8's retracted opposite claim.")


def test_superselectivity_sweep_excludes_few_or_no_cells_with_real_sigma():
    """
    UPDATED from v8's version of this test (which asserted the OPPOSITE:
    that SOME cells must be excluded). That was correct under the buggy
    sigma=1.0 default (which pushed most cells to float64 machine
    epsilon, requiring 91 of 104 to be excluded), but is no longer the
    expected behavior under the corrected, real sigma -- since the real
    sigma is ~4.4e11x larger, n*sigma*K_INTRA clears machine epsilon for
    nearly the entire swept range, and 0 cells were excluded when this
    was directly re-checked in the v8->v9 transition and again in this
    fix.

    This test checks that the exclusion COUNT dropped substantially
    relative to the full grid size, and specifically that it is much
    smaller than what the old buggy sigma produced (91 of 104) --
    without hard-pinning it to exactly 0, since a future change to the
    swept K_INTRA range could legitimately reintroduce a few edge-case
    exclusions without that being a bug.
    """
    result = sweep_superselectivity_reachability()
    k_intra_values = get_k_intra_sweep_values()
    full_grid_size = len(k_intra_values) * len(get_z_sweep_values())

    assert result["n_cells_excluded_numerically_unreliable"] < 0.2 * full_grid_size, (
        f"Expected fewer than 20% of the full grid "
        f"({0.2 * full_grid_size:.0f} of {full_grid_size}) to be excluded "
        f"as numerically unreliable with the real sigma -- got "
        f"{result['n_cells_excluded_numerically_unreliable']}, close to "
        f"or matching the ~87% (91 of 104) excluded under the OLD buggy "
        f"sigma=1.0 default. This could indicate the sigma bug has been "
        f"reintroduced."
    )
    assert result["grid_size"] > 0.8 * full_grid_size, (
        f"Expected the reliable grid ({result['grid_size']}) to cover "
        f"more than 80% of the full grid ({full_grid_size}) with the "
        f"real sigma correctly applied."
    )
    print(f"  [OK] with the real sigma, only "
          f"{result['n_cells_excluded_numerically_unreliable']} of "
          f"{full_grid_size} grid cells are excluded as numerically "
          f"unreliable (reliable grid size {result['grid_size']}) -- "
          f"confirming the real sigma clears float64 machine epsilon "
          f"across nearly the whole plausibility range, unlike the old "
          f"buggy sigma=1.0 default which excluded 91 of 104 cells.")


def test_egfr_flip_does_not_survive_strict_resolved_only_filter():
    """
    NEW (v14->v15). This is the direct regression test for the single
    most consequential finding from this project's independent diagnosis
    of the external Biomni report (see the handoff's "What changed since
    v13" and "Immediate next steps," item 0a). Prior to this session, the
    EGFR healthy/PDAC-plasma ranking flip was cited as a second robust
    demonstration alongside CA19-9. Independently re-derived this session:
    when only pixel-exact (non-density-estimated) points are used, EGFR's
    patch_exponent=2 ratio stays BELOW 1 at every tested K_INTRA value --
    the flip does NOT hold on strict-resolved-only data.

    Tolerance is set loosely around the checked range (0.885-0.935) so
    small floating-point differences don't spuriously fail this test,
    while still catching a real regression (e.g. if a future change
    accidentally fixed or masked the underlying data-provenance issue,
    which would be worth knowing about immediately, not silently).
    """
    result = strict_resolved_egfr_ca199_comparison()
    egfr = result["egfr"]
    assert egfr["flip_survives_strict_filter"] is False, (
        f"Expected the EGFR ranking flip to NOT survive strict-resolved-"
        f"only filtering (ratio_p2 should stay below 1 at every tested "
        f"K_INTRA value) -- got ratio_p2 range "
        f"[{egfr['ratio_p2_min']:.4f}, {egfr['ratio_p2_max']:.4f}], which "
        f"implies flip_survives_strict_filter=True. This would contradict "
        f"this session's (and Biomni's) independently-checked finding -- "
        f"investigate before citing EGFR as a robust flip again."
    )
    assert 0.7 < egfr["ratio_p2_min"] < 1.0, (
        f"Expected EGFR's minimum strict-resolved ratio_p2 in a loose "
        f"band around the checked value (~0.885) -- got "
        f"{egfr['ratio_p2_min']:.4f}."
    )
    assert 0.7 < egfr["ratio_p2_max"] < 1.0, (
        f"Expected EGFR's maximum strict-resolved ratio_p2 in a loose "
        f"band around the checked value (~0.935), still below 1.0 -- got "
        f"{egfr['ratio_p2_max']:.4f}."
    )
    print(f"  [OK] EGFR ranking flip does NOT survive strict-resolved-only "
          f"filtering: ratio_p2 stays in "
          f"[{egfr['ratio_p2_min']:.4f}, {egfr['ratio_p2_max']:.4f}], "
          f"below 1.0 throughout -- confirming EGFR should be retired as "
          f"a demonstrated ranking-flip case.")


def test_ca199_flip_survives_strict_resolved_only_filter():
    """
    NEW (v14->v15). The counterpart to the EGFR test above: CA19-9's
    ranking flip DOES survive the same strict-resolved-only filter,
    independently re-derived this session (ratio_p2 range approximately
    1.227-1.449, all clearly above 1). This is what makes CA19-9 the
    project's primary, headline Lennon-based ranking-flip demonstration
    per "Immediate next steps," item 0a -- it is the one that survives the
    strictest, most skeptical version of the check.
    """
    result = strict_resolved_egfr_ca199_comparison()
    ca199 = result["ca199"]
    assert ca199["flip_survives_strict_filter"] is True, (
        f"Expected the CA19-9 ranking flip to SURVIVE strict-resolved-only "
        f"filtering (ratio_p2 should stay above 1 at every tested K_INTRA "
        f"value) -- got ratio_p2 range "
        f"[{ca199['ratio_p2_min']:.4f}, {ca199['ratio_p2_max']:.4f}], "
        f"which implies flip_survives_strict_filter=False. If this is no "
        f"longer true, CA19-9 can no longer be described as the project's "
        f"robust headline finding -- investigate immediately, this is the "
        f"project's primary demonstrated result."
    )
    assert 1.0 < ca199["ratio_p2_min"] < 1.6, (
        f"Expected CA19-9's minimum strict-resolved ratio_p2 in a loose "
        f"band around the checked value (~1.227) -- got "
        f"{ca199['ratio_p2_min']:.4f}."
    )
    assert 1.0 < ca199["ratio_p2_max"] < 1.6, (
        f"Expected CA19-9's maximum strict-resolved ratio_p2 in a loose "
        f"band around the checked value (~1.449) -- got "
        f"{ca199['ratio_p2_max']:.4f}."
    )
    print(f"  [OK] CA19-9 ranking flip DOES survive strict-resolved-only "
          f"filtering: ratio_p2 stays in "
          f"[{ca199['ratio_p2_min']:.4f}, {ca199['ratio_p2_max']:.4f}], "
          f"above 1.0 throughout -- confirming CA19-9 as the project's "
          f"robust, primary Lennon-based ranking-flip demonstration.")


def test_panc1_hpdec_unaffected_by_strict_resolved_only_flag():
    """
    NEW (v14->v15). The PANC-1/HPDEC comparison's CSVs (panel_a/d_*) have
    no density_estimated rows to filter (unlike the fig3c_* healthy/PDAC-
    plasma CSVs) -- the digitization pipeline never needed a density-
    estimation step for those smaller, less-overlapping point clouds. This
    test confirms strict_resolved_only=True leaves "panc1_hpdec" results
    numerically identical to the full-reconstruction run, so the new flag
    cannot silently and unexpectedly change this unrelated comparison.
    """
    full = sweep_ranking_flip_robustness(strict_resolved_only=False)
    strict = sweep_ranking_flip_robustness(strict_resolved_only=True)

    full_rows = full["panc1_hpdec"]
    strict_rows = strict["panc1_hpdec"]
    assert len(full_rows) == len(strict_rows)
    for r_full, r_strict in zip(full_rows, strict_rows):
        assert r_full["k_intra"] == r_strict["k_intra"]
        assert abs(r_full["ratio_p2"] - r_strict["ratio_p2"]) < 1e-9, (
            "panc1_hpdec's ratio_p2 should be unaffected by "
            "strict_resolved_only, since panel_a/d_* CSVs have no "
            "density_estimated rows."
        )
        assert abs(r_full["ratio_p0"] - r_strict["ratio_p0"]) < 1e-9, (
            "panc1_hpdec's ratio_p0 should be unaffected by "
            "strict_resolved_only, since panel_a/d_* CSVs have no "
            "density_estimated rows."
        )
    print("  [OK] panc1_hpdec results are numerically identical with "
          "strict_resolved_only=True vs False, as expected (no "
          "density_estimated rows exist in those CSVs).")


def test_healthy_pdac_plasma_backward_compatible_key_unaffected_by_flag():
    """
    NEW (v14->v15). "healthy_pdac_plasma" (the original key, predating
    this session's EGFR/CA19-9 split) is documented to always return
    full-reconstruction EGFR data regardless of strict_resolved_only, so
    that any existing caller reading that exact key is not silently
    changed by adding the new flag. This locks in that backward-
    compatibility guarantee directly.
    """
    full = sweep_ranking_flip_robustness(strict_resolved_only=False)
    strict = sweep_ranking_flip_robustness(strict_resolved_only=True)

    for r_full, r_strict in zip(full["healthy_pdac_plasma"], strict["healthy_pdac_plasma"]):
        assert abs(r_full["ratio_p2"] - r_strict["ratio_p2"]) < 1e-9, (
            "'healthy_pdac_plasma' must stay full-reconstruction "
            "regardless of strict_resolved_only, for backward "
            "compatibility with existing callers/tests."
        )
    # And it should match the explicit "_egfr" full-reconstruction alias.
    for r_a, r_b in zip(full["healthy_pdac_plasma"], full["healthy_pdac_plasma_egfr"]):
        assert abs(r_a["ratio_p2"] - r_b["ratio_p2"]) < 1e-12
    print("  [OK] 'healthy_pdac_plasma' key remains full-reconstruction "
          "EGFR data regardless of strict_resolved_only, and matches the "
          "explicit 'healthy_pdac_plasma_egfr' alias exactly.")


def test_ca199_accessibility_equal_case_matches_no_effect_case():
    """
    NEW (this session). binding_model.py's own proof (see
    selectivity_ratio_sensitivity_to_accessibility()'s docstring) is that
    EQUAL accessibility between tumor and normal populations cancels out
    of the tumor/normal ratio exactly, regardless of its absolute value.
    This locks in that same guarantee at the level of the actual CA19-9
    real-data comparison: "no_effect" (f_accessible=1.0 for both) and
    "illustrative_equal_direction" (f_accessible=0.7 for both) must give
    IDENTICAL ratios at every swept K_INTRA value.
    """
    results = ca199_selectivity_sensitivity_to_accessibility()
    no_effect = results["no_effect"]
    equal_dir = results["illustrative_equal_direction"]

    for r_no, r_eq in zip(no_effect, equal_dir):
        assert abs(r_no["k_intra"] - r_eq["k_intra"]) < 1e-30
        assert abs(r_no["ratio"] - r_eq["ratio"]) < 1e-9, (
            f"Equal accessibility (0.7 vs 0.7) should give an IDENTICAL "
            f"ratio to no accessibility effect (1.0 vs 1.0) at "
            f"k_intra={r_no['k_intra']:.3e} -- got {r_no['ratio']:.6f} vs "
            f"{r_eq['ratio']:.6f}."
        )
    print("  [OK] Equal tumor/normal accessibility leaves the real CA19-9 "
          "ratio completely unchanged from the no-effect baseline, "
          "confirming binding_model.py's cancellation proof holds against "
          "real digitized data too.")


def test_ca199_accessibility_unequal_case_weakens_but_does_not_erase_flip():
    """
    NEW (this session). The illustrative unequal-accessibility scenario
    (tumor 30% less accessible than normal, on top of a shared 0.7
    baseline -- see the source function's docstring for why these are
    explicitly labeled illustrative, not measured) should substantially
    REDUCE the CA19-9 ratio compared to the no-effect baseline, since
    reduced tumor accessibility directly reduces tumor yield relative to
    normal. This is checked directly rather than assumed, and the ratio
    is confirmed to stay above 1.0 in this specific tested scenario (the
    flip survives, though weakened) -- not asserted to always survive
    ANY accessibility difference, since a large enough one could reverse
    it entirely (a genuinely open, real-data-dependent question this
    project cannot currently answer without real tumor-vs-normal
    accessibility measurements -- see the source function's docstring).
    """
    results = ca199_selectivity_sensitivity_to_accessibility()
    no_effect_ratios = [r["ratio"] for r in results["no_effect"]]
    unequal_ratios = [r["ratio"] for r in results["illustrative_unequal_direction"]]

    for r_no, r_un in zip(no_effect_ratios, unequal_ratios):
        assert r_un < r_no, (
            "Reduced tumor accessibility (relative to normal) should "
            "REDUCE the tumor/normal ratio compared to the no-effect "
            "baseline at every k_intra value."
        )
    assert all(r > 1.0 for r in unequal_ratios), (
        "In this specific tested illustrative scenario, the CA19-9 flip "
        "should still survive (ratio>1) though weakened -- if this now "
        "fails, that's a real, worth-reporting change in the sensitivity "
        "finding, not something to silently loosen this test for."
    )
    print("  [OK] Illustrative unequal accessibility (tumor 30% less "
          "accessible than normal) substantially weakens but does not "
          "erase the CA19-9 ratio in this specific tested scenario -- "
          "ratios dropped from ~1.23-1.45 (no effect) to ~1.01-1.03.")


def test_ca199_accessibility_note_documents_the_real_data_gap():
    """
    NEW (this session). The function's returned 'note' must actually state
    the honest finding (real corona-masking evidence exists, but no real
    tumor-vs-normal comparison was found) rather than silently becoming
    generic boilerplate in a future edit -- checked for the key factual
    content, not exact wording.
    """
    results = ca199_selectivity_sensitivity_to_accessibility()
    note = results["note"]
    assert "10.1002/jev2.70124" in note, "Note should cite the real DOI this session retrieved."
    assert "not" in note.lower() and ("compare" in note.lower() or "comparison" in note.lower()), (
        "Note should state plainly that the cited paper does not compare "
        "tumor vs. normal EVs."
    )
    print("  [OK] The accessibility-sensitivity function's note honestly "
          "documents what real literature does and does not support here.")


if __name__ == "__main__":
    tests = [
        test_healthy_pdac_plasma_flip_is_robust_across_full_k_intra_range,
        test_panc1_hpdec_stays_near_one_not_a_clean_flip,
        test_ranking_flip_sweep_uses_real_sigma_not_placeholder,
        test_superselectivity_sweep_uses_real_sigma_not_placeholder,
        test_genuine_superselectivity_is_reached_in_current_plausibility_range,
        test_superselectivity_sweep_excludes_few_or_no_cells_with_real_sigma,
        test_egfr_flip_does_not_survive_strict_resolved_only_filter,
        test_ca199_flip_survives_strict_resolved_only_filter,
        test_panc1_hpdec_unaffected_by_strict_resolved_only_flag,
        test_healthy_pdac_plasma_backward_compatible_key_unaffected_by_flag,
        test_ca199_accessibility_equal_case_matches_no_effect_case,
        test_ca199_accessibility_unequal_case_weakens_but_does_not_erase_flip,
        test_ca199_accessibility_note_documents_the_real_data_gap,
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
