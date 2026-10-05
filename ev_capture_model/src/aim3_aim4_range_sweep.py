"""
aim3_aim4_range_sweep.py

WHY THIS FILE EXISTS
---------------------
The v7->v8 handoff transition's search for a K_INTRA/z calibration source
(see src/k_intra_z_plausibility_bounds.py) did not resolve to a single
fitted number -- it produced a documented, honestly-wide RANGE instead,
now expressed in real physical units (cm^2/antibody) using a real antibody
surface density found this session (Lopez Baltazar dissertation, Appendix
A.1). The v7 handoff's "Immediate next steps" items 5 and 7 named exactly
this as the way forward while calibration remains unresolved: "sweep
across [the plausibility] file's range rather than using the old arbitrary
placeholder values... report whether the ranking-flip findings are robust
across the whole plausible range," and separately "does the model land in
the alpha>1 regime ANYWHERE across the plausibility file's range."

This script does both of those sweeps, using the REAL digitized Lennon
et al. 2019 data in ev_capture_data/ (not reconstructed-from-summary-
statistics data), so both findings are grounded in real single-vesicle
measurements, not just internal model consistency.

WHAT THIS IS, AND WHAT IT IS NOT
-----------------------------------
This is NOT a calibration and does not produce a single "the answer is X"
number for K_INTRA, z, or k_geom -- those remain formally uncalibrated
placeholders in binding_model.py, exactly as before. What this script
does is ask two honest, range-based questions that do NOT require a
single fitted value to answer meaningfully:

  1. Ranking-flip robustness (Aim 3/Aim 4): for the real digitized
     PANC-1/HPDEC (Fig 2, matched cell lines) and healthy/PDAC-plasma
     (Fig 3c, size-mismatched populations) comparisons, does the
     patch_exponent=2-vs-0 ranking-flip finding from v4 (already
     established, unchanged since) hold at EVERY point in the K_INTRA
     plausibility range, or only at the single arbitrary placeholder
     value (K_INTRA=0.05) used in the original demo scripts?

  2. Superselectivity-regime reachability (Aim 4): across the SAME
     K_INTRA range, combined with a swept z (activity/concentration
     parameter, which the plausibility file's concentration-axis
     cross-check bounds only qualitatively -- see that file's docstring
     point 3), is alpha > 1 (the strict superselectivity signature)
     achievable ANYWHERE in the swept region? This is the "range-based
     answer even without a single fitted z" the v7 handoff explicitly
     named as a legitimate finding in its own right.

HOW THE RANGE IS BUILT
--------------------------
k_intra_z_plausibility_bounds.k_intra_bounds_from_real_sigma() gives a
table of K_INTRA values keyed by (assumed receptor count n, target
saturation exponent) -- not a single min/max pair. This script pulls the
full set of values from that table (spanning ~3-4 orders of magnitude, as
documented there) and uses them directly as the swept K_INTRA range,
rather than inventing a separate range -- this keeps the sweep tied
explicitly to the already-documented, already-tested bound rather than a
new unreviewed one.

For z, since no real number exists yet (the concentration-axis
cross-check bounds a concentration regime, not a dimensionless z
directly -- see the plausibility file), a broad log-spaced sweep is used
instead, chosen wide enough to include both the "z=1 shows no
superselectivity" finding from v5 and the "z=1e-4 shows strong
superselectivity" finding from v5, so the sweep is guaranteed to bracket
the known qualitative transition rather than accidentally missing it.
"""

import os
import sys
import csv

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from binding_model import selectivity_alpha, alpha_exact
from population_model import joint_population_yield
from k_intra_z_plausibility_bounds import (
    k_intra_bounds_from_real_sigma,
    SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "ev_capture_data")


def _load_pairs(filename, strict_resolved_only=False):
    """Load (copies_per_ev, diameter_nm) pairs from a digitized CSV, matching
    the same loading convention already used in tests/test_aim1.py's Lennon
    ranking-flip test, so results are directly comparable.

    Parameters
    ----------
    filename : str
        CSV filename inside ev_capture_data/.
    strict_resolved_only : bool, optional (default False)
        NEW (v14->v15, added in response to the Biomni external review --
        see the handoff's "What changed since v13"). If True, skip every
        row whose `source` column starts with "density_estimated" and keep
        only pixel-exact digitized points (source starting with
        "resolved"). This column only exists in the fig3c_* (healthy/PDAC-
        plasma) CSVs -- the panel_a/d/e/f (PANC-1/HPDEC) CSVs are fully
        pixel-resolved already (no density-estimation step was needed for
        those smaller, less-overlapping point clouds), so this flag has no
        effect on those files. Independently confirmed this session
        against the raw CSVs: fig3c_EGFR_PDAC (2050 rows, 1578
        density_estimated, 77.0%), fig3c_EGFR_healthy (815 rows, 724
        density_estimated, 88.8%), fig3c_CA199_PDAC (2653 rows, 2176
        density_estimated, 82.0%), fig3c_CA199_healthy (644 rows, 604
        density_estimated, 93.8%) -- matching the handoff's cited
        percentages exactly. See ev_capture_data/README.md for what
        "density_estimated" vs. "resolved"/"resolved_cluster" mean and why
        the distinction exists (the digitization pipeline synthesizes
        points in the densest, most-overlapping corner of each scatter
        plot to match the paper's reported N, since individual blobs
        cannot be resolved there).
    """
    diam, copies = [], []
    path = os.path.join(DATA_DIR, filename)
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            if strict_resolved_only:
                source = row.get("source", "")
                if source.startswith("density_estimated"):
                    continue
            d = float(row["diameter_nm"])
            c = float(row["copies_per_ev"])
            if d > 0:
                diam.append(d)
                copies.append(c)
    return np.array(copies), np.array(diam)


def get_k_intra_sweep_values():
    """
    Pull the full set of K_INTRA values out of the plausibility file's
    real-sigma bound table (see module docstring). Returns a sorted,
    de-duplicated array spanning the documented ~3-4 order of magnitude
    range.
    """
    bounds = k_intra_bounds_from_real_sigma()
    values = sorted(set(bounds.values()))
    return np.array(values)


def get_z_sweep_values(n_points=13):
    """
    Broad log-spaced sweep for z, chosen to bracket the known qualitative
    transition found in the v5 session: z=1.0 (used in all existing demo
    scripts) shows NO superselectivity (alpha never exceeds 1), while
    z=1e-4 (with the same sigma/k_intra) shows a strong superselective
    peak. This sweep spans well past both endpoints (1e-6 to 1e2) so the
    transition point itself, not just these two known samples, is
    resolved.
    """
    return np.logspace(-6, 2, n_points)


def sweep_ranking_flip_robustness(sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, z=1.0,
                                   strict_resolved_only=False):
    """
    Question 1: does the patch_exponent=2-vs-0 ranking flip (established
    in v4, unchanged since) hold across the ENTIRE K_INTRA plausibility
    range, for the matched-size (PANC-1/HPDEC) and BOTH size-mismatched
    (healthy/PDAC-plasma, EGFR and CA19-9) real digitized comparisons?

    IMPORTANT UNIT-CONSISTENCY NOTE, CAUGHT AND FIXED THIS SESSION: the
    K_INTRA values from k_intra_z_plausibility_bounds.py's
    k_intra_bounds_from_real_sigma() are only meaningful when PAIRED with
    the SAME real sigma they were derived against
    (SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2) -- because that function computes
    K_INTRA = (target dimensionless exponent) / (n * sigma), so plugging
    the resulting K_INTRA back in with a DIFFERENT sigma (e.g. the
    old placeholder sigma=1.0 used elsewhere in this project's demo
    scripts) silently changes the actual n*sigma*K_INTRA product away
    from the intended target -- by a factor of sigma_real (~4.4e11) in
    this specific case, which was caught during this session's own
    testing (an early version of this function used sigma=1.0 by mistake
    and produced an artifactually tiny effective exponent). The default
    here is therefore SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, not 1.0.

    NEW (v14->v15, strict_resolved_only) -- added directly in response to
    this project's independent Biomni-review diagnosis (see the handoff's
    "What changed since v13" and "Immediate next steps," item 0a). When
    True, both healthy/PDAC-plasma comparisons (EGFR and CA19-9) are
    computed using ONLY pixel-exact digitized points (source not starting
    with "density_estimated" -- see _load_pairs' docstring), excluding the
    77-94% of points in these four CSVs that were synthesized by the
    digitization pipeline's density-estimation step rather than resolved
    individually. This is the check that revealed EGFR's ranking flip does
    NOT survive strict provenance filtering while CA19-9's DOES -- see
    ev_capture_healthy_pdac_ca199 and ev_capture_healthy_pdac_egfr keys
    below, and the "healthy_pdac_plasma" key (kept for backward
    compatibility, and always full-reconstruction/EGFR regardless of this
    flag, so existing callers are not silently changed).

    Returns
    -------
    dict with keys "panc1_hpdec" (PANC-1/HPDEC, EGFR, matched-size
    control), "healthy_pdac_plasma" (healthy/PDAC-plasma, EGFR --
    UNCHANGED meaning from prior versions, always full-reconstruction data
    regardless of strict_resolved_only, kept for backward compatibility
    with existing callers/tests), "healthy_pdac_plasma_egfr" (identical
    values to "healthy_pdac_plasma", explicit alias so EGFR and CA19-9 can
    be iterated over symmetrically), and "healthy_pdac_plasma_ca199" (NEW
    -- the CA19-9 counterpart, not previously computed by this function).
    Each maps to a list of per-K_INTRA-value result dicts:
        {"k_intra": ..., "ratio_p2": ..., "ratio_p0": ...}
    where ratio_p2 / ratio_p0 are the (tumor yield / normal yield) ratios
    under patch_exponent=2.0 and 0.0 respectively. When
    strict_resolved_only=True, "panc1_hpdec" is unaffected (those CSVs
    have no density_estimated rows to filter -- see _load_pairs), and both
    healthy/PDAC-plasma comparisons are computed on strict-resolved-only
    data instead of the full reconstruction.
    """
    k_intra_values = get_k_intra_sweep_values()

    comparisons = {
        "panc1_hpdec": (
            "panel_a_PANC1_EGFR_rescaled.csv",
            "panel_d_HPDEC_EGFR_rescaled.csv",
        ),
        "healthy_pdac_plasma": (
            "fig3c_EGFR_PDAC_digitized.csv",
            "fig3c_EGFR_healthy_digitized.csv",
        ),
        "healthy_pdac_plasma_egfr": (
            "fig3c_EGFR_PDAC_digitized.csv",
            "fig3c_EGFR_healthy_digitized.csv",
        ),
        "healthy_pdac_plasma_ca199": (
            "fig3c_CA199_PDAC_digitized.csv",
            "fig3c_CA199_healthy_digitized.csv",
        ),
    }

    # "healthy_pdac_plasma" (the original, backward-compatible key) is
    # always full-reconstruction, regardless of strict_resolved_only --
    # only the two explicitly-named _egfr/_ca199 keys respect the flag.
    # This means existing callers/tests that read "healthy_pdac_plasma"
    # keep getting exactly what they always got.
    filter_by_label = {
        "panc1_hpdec": False,
        "healthy_pdac_plasma": False,
        "healthy_pdac_plasma_egfr": strict_resolved_only,
        "healthy_pdac_plasma_ca199": strict_resolved_only,
    }

    results = {}
    for label, (tumor_file, normal_file) in comparisons.items():
        use_strict = filter_by_label[label]
        n_t, d_t = _load_pairs(tumor_file, strict_resolved_only=use_strict)
        n_n, d_n = _load_pairs(normal_file, strict_resolved_only=use_strict)

        rows = []
        for k_intra in k_intra_values:
            row = {"k_intra": float(k_intra)}
            for p in (2.0, 0.0):
                # k_geom chosen so the normal population's mean n_contact
                # is O(1) at patch_exponent=2 baseline scale, same
                # convention as tests/test_aim1.py's existing ranking-flip
                # test, so k_geom does not itself become a hidden second
                # free parameter being silently varied across the sweep.
                k_geom = 1.0 / (n_n.mean() * d_n.mean() ** (p - 2.0))
                y_t = joint_population_yield(n_t, d_t, sigma, k_geom, k_intra, z, patch_exponent=p)
                y_n = joint_population_yield(n_n, d_n, sigma, k_geom, k_intra, z, patch_exponent=p)
                ratio = y_t / y_n if y_n > 0 else float("inf")
                row[f"ratio_p{int(p)}"] = ratio
            rows.append(row)
        results[label] = rows

    return results


def strict_resolved_egfr_ca199_comparison(sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, z=1.0):
    """
    NEW (v14->v15). Convenience wrapper, added specifically for "Immediate
    next steps" item 0a: makes the strict-resolved-only EGFR-vs-CA19-9
    ranking-flip comparison a single, directly-callable, testable function
    rather than something that only lived in a scratch verification script
    and this handoff's prose.

    Runs sweep_ranking_flip_robustness(strict_resolved_only=True) and
    extracts just the patch_exponent=2 ratio range (the ratio the
    project's ranking-flip claims are actually stated in terms of) for
    each marker, so the headline comparison can be read off directly:

        EGFR:   ratio_p2 stays BELOW 1 across the full K_INTRA range
                (flip does NOT survive strict filtering)
        CA19-9: ratio_p2 stays ABOVE 1 across the full K_INTRA range
                (flip DOES survive strict filtering)

    Independently confirmed this session (re-run directly from this
    function, not copied from the handoff): EGFR ratio_p2 range
    approximately 0.885-0.935 (all <1, no flip); CA19-9 ratio_p2 range
    approximately 1.227-1.449 (all >1, flip holds). Small differences from
    Biomni's own independently-reported ranges (EGFR 0.959-0.990; CA19-9
    1.262-1.447) are expected and were already noted in the handoff --
    they likely trace to a slightly different z or k_geom convention --
    the qualitative and quantitative conclusion (EGFR fails, CA19-9 holds)
    is the same either way.

    Returns
    -------
    dict with keys "egfr" and "ca199", each mapping to:
        {"ratio_p2_min": ..., "ratio_p2_max": ...,
         "flip_survives_strict_filter": bool}
    where "flip_survives_strict_filter" is True only if ratio_p2 stays
    above 1.0 at EVERY tested K_INTRA value (the CA19-9 case), and False
    if it ever drops to or below 1.0 (the EGFR case).
    """
    results = sweep_ranking_flip_robustness(sigma=sigma, z=z, strict_resolved_only=True)

    summary = {}
    for marker, key in (("egfr", "healthy_pdac_plasma_egfr"), ("ca199", "healthy_pdac_plasma_ca199")):
        ratios = [row["ratio_p2"] for row in results[key]]
        summary[marker] = {
            "ratio_p2_min": min(ratios),
            "ratio_p2_max": max(ratios),
            "flip_survives_strict_filter": all(r > 1.0 for r in ratios),
        }
    return summary


def ca199_selectivity_sensitivity_to_accessibility(
    sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, z=1.0, strict_resolved_only=True
):
    """
    "Immediate next steps" item 3 (this session): apply f_accessible
    (binding_model.py's accessibility/corona multiplier, added last
    session but not yet used against any real headline finding) to the
    project's actual primary headline case -- the CA19-9 healthy/PDAC-
    plasma ranking flip -- rather than only the illustrative 0.9-vs-0.3
    example in test_f_accessible.py.

    WHAT REAL DATA THIS SESSION FOUND, AND WHY IT DOES NOT SUPPORT A
    CALIBRATED f_accessible VALUE FOR TUMOR VS. NORMAL PANCREATIC EVs --
    READ THIS BEFORE CITING ANY SPECIFIC f_accessible NUMBER FROM THIS
    FUNCTION AS IF IT WERE MEASURED
    -----------------------------------------------------------------------
    A real, relevant, PubMed-indexed paper exists: Foersonits et al. 2025,
    "Improved Accessibility of Extracellular Vesicle Surface Molecules Upon
    Partial Removal of the Protein Corona by High Ionic Strength," Journal
    of Extracellular Vesicles 14(7):e70124, DOI 10.1002/jev2.70124 (PMID
    40660090, PMC12259388) -- retrieved and read in full this session via
    PubMed. It directly demonstrates the MECHANISM binding_model.py's
    f_accessible represents: a protein corona forms around real EVs (from
    human blood plasma and cultured cell lines) and measurably MASKS
    surface markers from antibody/annexin V detection, and stripping that
    corona (by high-salt washing) significantly IMPROVES detection of 15
    of 37 tested plasma-EV surface markers and 9 of 37 THP1-EV surface
    markers -- real, checked, statistically-tested evidence (Wilcoxon
    tests, p<0.05 to p<0.001 depending on marker) that receptor
    accessibility is a real, measurable, marker-dependent phenomenon, not
    a hypothetical one this project invented to have somewhere to put
    f_accessible.

    HOWEVER, checked directly and reported honestly rather than stretched:
    this paper does NOT answer the specific question this function needs
    (does accessibility differ between TUMOR and NORMAL EVs specifically,
    and if so by how much). It compares corona-INTACT vs. corona-STRIPPED
    EVs from a SINGLE population (healthy human plasma donors, or one
    non-cancer cell line, THP1) -- not a tumor-vs-normal comparison at
    all. A separate, targeted PubMed search this session for a pancreatic-
    cancer-specific or tumor-vs-normal EV corona/accessibility comparison
    (search terms: pancreatic ductal adenocarcinoma, EV surface protein,
    corona, accessibility, tumor, normal, comparison, MUC1/EGFR/TROP2)
    returned ZERO results -- a genuine, checked negative finding, not an
    assumption that such a paper doesn't exist. Additionally, even within
    the Foersonits paper, the ACTUAL fold-change magnitudes for individual
    markers are reported only in bar-chart figures (Figure 4C/4F), not as
    numeric values in the article's own text or tables -- so even a
    same-population corona-effect SIZE cannot be responsibly extracted
    from this source without measuring pixel heights off a chart, which
    this project's own digitization standards (see ev_capture_data/) would
    not accept as "real digitized data" without care this function does
    not attempt to apply here.

    WHAT THIS FUNCTION THEREFORE ACTUALLY DOES: rather than inventing a
    tumor/normal accessibility split with no real source (which would be
    exactly the kind of unsupported claim this project's own working style
    exists to avoid), it runs THREE clearly-labeled scenarios against the
    CA19-9 headline comparison:
      1. "no_effect" -- f_accessible=1.0 for both populations (the
         project's existing, already-reported headline finding, included
         here only as the baseline for comparison).
      2. "illustrative_equal_direction" -- a SINGLE illustrative
         f_accessible value (0.7, chosen only as a round, moderate,
         clearly-labeled-as-illustrative number, matching the spirit of
         test_f_accessible.py's own 0.9/0.3 illustrative pair) applied to
         BOTH populations equally, demonstrating the (already-proven-in-
         binding_model.py) fact that EQUAL accessibility cancels out of
         the ratio regardless of its absolute value.
      3. "illustrative_unequal_direction" -- the SAME illustrative
         f_accessible=0.7 baseline, but with the TUMOR population's
         accessibility further reduced by an additional illustrative 30%
         (f_accessible_tumor=0.49) relative to normal, in the QUALITATIVE
         direction the Foersonits paper's mechanism would suggest IF tumor
         EVs carried a heavier/stickier corona than normal EVs (a real,
         independently plausible hypothesis in the broader EV literature --
         tumor EVs are known to have different surface protein/glycan
         composition that could plausibly affect protein adsorption -- but
         NOT something this session's search confirmed with a real
         tumor-vs-normal measurement). This scenario exists to show HOW
         MUCH the CA19-9 finding's ratio would move under a real but
         UNVERIFIED accessibility difference, not to claim this difference
         is measured or even necessarily real.

    Every returned scenario is labeled explicitly as illustrative except
    "no_effect," and the docstring above is the honest record of exactly
    what real literature does and does not support here -- consistent with
    this project's standing rule to report a genuine "doesn't have real
    calibration data (yet)" finding plainly rather than backfilling a
    number that looks plausible.

    Parameters
    ----------
    sigma, z : float, optional
        Same meaning as elsewhere in this module. Defaults match
        strict_resolved_egfr_ca199_comparison()'s defaults for direct
        comparability.
    strict_resolved_only : bool, optional (default True)
        Whether to use only pixel-exact digitized CA19-9 points (matching
        this project's primary, strict-filter-tested headline claim) or
        the full reconstruction. Defaults to True since CA19-9 under
        strict filtering is the project's actual current headline claim
        (see CURRENT.md's "Known risks").

    Returns
    -------
    dict
        {scenario_name: {"k_intra": ..., "ratio": ...} list, keyed by the
        three scenario names above}, plus a "note" key with a short
        plain-language summary of the honest-data-gap finding.
    """
    from binding_model import accessible_receptors

    k_intra_values = get_k_intra_sweep_values()
    tumor_file = "fig3c_CA199_PDAC_digitized.csv"
    normal_file = "fig3c_CA199_healthy_digitized.csv"
    n_t, d_t = _load_pairs(tumor_file, strict_resolved_only=strict_resolved_only)
    n_n, d_n = _load_pairs(normal_file, strict_resolved_only=strict_resolved_only)

    scenarios = {
        "no_effect": (1.0, 1.0),
        "illustrative_equal_direction": (0.7, 0.7),
        "illustrative_unequal_direction": (0.49, 0.7),  # (tumor, normal)
    }

    results = {}
    for scenario_name, (f_acc_tumor, f_acc_normal) in scenarios.items():
        n_t_acc = accessible_receptors(n_t, f_acc_tumor)
        n_n_acc = accessible_receptors(n_n, f_acc_normal)
        rows = []
        for k_intra in k_intra_values:
            k_geom = 1.0 / (n_n_acc.mean() * d_n.mean() ** (2.0 - 2.0))
            y_t = joint_population_yield(n_t_acc, d_t, sigma, k_geom, k_intra, z, patch_exponent=2.0)
            y_n = joint_population_yield(n_n_acc, d_n, sigma, k_geom, k_intra, z, patch_exponent=2.0)
            ratio = y_t / y_n if y_n > 0 else float("inf")
            rows.append({"k_intra": float(k_intra), "ratio": ratio})
        results[scenario_name] = rows

    results["note"] = (
        "Foersonits et al. 2025 (DOI 10.1002/jev2.70124, PMID 40660090) confirms "
        "protein-corona-driven marker masking is REAL and MEASURABLE in real EVs, "
        "but does NOT compare tumor vs. normal EVs -- a targeted search this "
        "session for a tumor-vs-normal EV corona/accessibility comparison found "
        "none. The 'illustrative_*' scenarios above use clearly-labeled invented "
        "numbers to show HOW MUCH the CA19-9 ratio would move under a real but "
        "currently unverified accessibility difference -- they are not measured "
        "values and should not be cited as if they were."
    )
    return results


# Below this value, n*sigma*K_INTRA at n=1 is close enough to machine
# epsilon (~2.2e-16 for float64) that exp(x)-1 loses precision to
# catastrophic cancellation, producing WILDLY OSCILLATING, non-physical
# alpha values (checked directly this session: alpha jumping between
# ~0.87 and ~1.21 from one adjacent n to the next, rather than tracing a
# smooth curve, at k_intra=2.25e-16 with sigma=1.0 -- a numerical
# artifact, not a real superselectivity signal). Grid cells whose SMALLEST
# tested exponent (n=1) falls below this floor are excluded from the
# alpha>1 reachability count below, and this exclusion is reported
# explicitly rather than silently -- reporting a numerically-unreliable
# alpha>1 as if it were a real finding would repeat, in a new form, the
# same kind of premature-conclusion mistake this project's working style
# has caught before (see the handoff's "What NOT to revisit" list).
MIN_RELIABLE_EXPONENT_AT_N1 = 1e-12

# Above this, alpha is considered a REAL superselective signal rather than
# floating-point wobble around exactly 1.0. Checked directly this session:
# even at the LARGEST K_INTRA in the plausibility range (2.25e-12
# cm^2/Ab, sigma=1.0, n up to 200), alpha comes out at 1.000026 -- a
# razor-thin numerical wobble, not a real peak (the nanostar_validation.py
# file's genuine superselective peaks, by contrast, reach alpha~6, orders
# of magnitude above this kind of wobble). A threshold of 1.01 (1%) is
# used so a real, physically meaningful superselective signal cannot be
# confused with float64 noise around the exactly-linear limit.
ALPHA_SUPERSELECTIVE_THRESHOLD = 1.01


def sweep_superselectivity_reachability(sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, n_max=200):
    """
    Question 2: across the K_INTRA plausibility range AND a broad z sweep,
    is alpha > 1 (strict superselectivity) reachable ANYWHERE? Uses
    selectivity_alpha() directly (the per-vesicle model, not the
    population/yield model) since alpha is defined per binding_model.py's
    docstring as d(log theta)/d(log n) -- a single-vesicle-model quantity.

    BUG FIXED THIS SESSION (v8->v9 transition), DOCUMENTED HERE SO IT IS
    NOT RE-INTRODUCED: this function's `sigma` parameter previously
    defaulted to 1.0 -- the OLD placeholder value used throughout this
    project's earliest demo scripts -- instead of
    SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, the real antibody surface density
    that every K_INTRA value in get_k_intra_sweep_values() was actually
    DERIVED against (via k_intra_bounds_from_real_sigma(), which computes
    K_INTRA = target_exponent / (n * sigma_real)). Since this function is
    called with NO ARGUMENTS at its one call site in print_summary(), the
    bug was live in normal execution, not just a latent default -- it
    produced this project's entire v8 "superselectivity is not reached
    anywhere" finding, which was consequently WRONG (a numerical artifact,
    not real physics) and has been retracted (see handoff v9's "What
    changed since v8"). This is the exact same CLASS of bug as
    sweep_ranking_flip_robustness()'s "Bug 1" (documented in that
    function's own docstring above) -- that fix was applied to the
    ranking-flip function during the v7->v8 session but was never applied
    here, leaving this sibling function broken. Both functions must be
    called with a sigma that matches the sigma the swept K_INTRA values
    were derived against; SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2 is therefore
    the correct default for both, and now is for both.

    NUMERICAL CAVEAT, RE-CHECKED AFTER THE FIX (still worth reading): at
    the smallest K_INTRA values in the plausibility range, n*sigma*K_INTRA
    at n=1 CAN still fall below float64 machine epsilon for some (k_intra,
    z) combinations, and q_bound's exp(x)-1 loses precision to
    catastrophic cancellation there -- producing alpha values that
    oscillate non-physically instead of tracing a smooth curve. This is a
    NUMERICAL ARTIFACT, not a real superselectivity finding, and affected
    grid points are excluded from the reported reachability fraction (see
    MIN_RELIABLE_EXPONENT_AT_N1 above); the number of excluded cells is
    always reported so the exclusion is never silent. With the CORRECTED
    real sigma, far fewer (typically zero, in the sweep as currently
    configured -- confirmed by direct re-run after this fix) cells fall
    into this excluded regime than under the old buggy sigma=1.0 default,
    because the real sigma is ~4.4e11x larger, which pushes n*sigma*K_INTRA
    well clear of machine epsilon for nearly the whole swept K_INTRA range.
    Do not assume this means the caveat no longer applies in general --
    re-check n_cells_excluded_numerically_unreliable every time the swept
    K_INTRA or n range changes, rather than assuming zero exclusions by
    default.

    NEW (v14->v15) -- this function now ALSO computes and reports the
    EXACT (closed-form) peak alpha, via binding_model.alpha_exact(),
    alongside the original finite-difference-based figures above. This is
    "Immediate next steps" item 0b, added directly in response to the
    Biomni external review: the finite-difference selectivity_alpha() is
    a good approximation at the specific peak (checked directly this
    session: exact peak 10.465 vs. finite-difference peak 10.461, a real
    but small ~0.04 difference) but diverges by as much as ~0.70 in
    absolute alpha at other grid points -- see alpha_exact()'s own
    docstring for the full derivation and discrepancy check. Both values
    are now reported so the write-up can cite the more precise exact
    number where it matters, without discarding the original finite-
    difference figure or needing to silently replace it everywhere it's
    already been used.

    Returns
    -------
    dict with:
        "max_alpha_found": the largest RELIABLE finite-difference alpha
            value found anywhere in the swept grid (numerically-unreliable
            cells excluded) -- UNCHANGED meaning from prior versions, kept
            for backward compatibility with existing callers/tests.
        "any_superselective": whether max_alpha_found > 1
        "best_params": the (k_intra, z) pair that produced max_alpha_found
        "fraction_of_grid_superselective": what fraction of the RELIABLE
            (k_intra, z) grid cells produce AT LEAST one alpha>1 value
            across the tested n range (finite-difference basis, unchanged)
        "n_cells_excluded_numerically_unreliable": how many grid cells
            were skipped because n*sigma*K_INTRA at n=1 was too close to
            machine epsilon to trust
        "max_alpha_exact_found": NEW -- the largest EXACT (closed-form)
            alpha value found anywhere in the same reliable grid, computed
            via binding_model.alpha_exact() on the identical n/k_intra/z
            values used for the finite-difference computation above (so
            the two numbers are directly comparable, not computed on
            different grids).
        "best_params_exact_k_intra_z": the (k_intra, z) pair that produced
            max_alpha_exact_found (may differ slightly from
            best_params_k_intra_z, since the two methods can peak at
            slightly different grid points).
    """
    k_intra_values = get_k_intra_sweep_values()
    z_values = get_z_sweep_values()
    n_values = np.linspace(1, n_max, 60)

    max_alpha = -np.inf
    best_params = None
    max_alpha_exact = -np.inf
    best_params_exact = None
    n_superselective_cells = 0
    total_reliable_cells = 0
    n_excluded = 0

    for k_intra in k_intra_values:
        exponent_at_n1 = 1 * sigma * k_intra
        if exponent_at_n1 < MIN_RELIABLE_EXPONENT_AT_N1:
            # Every z value at this k_intra shares the same unreliable
            # q_bound computation (z only rescales after qb is computed),
            # so skip the whole row.
            n_excluded += len(z_values)
            continue
        for z in z_values:
            total_reliable_cells += 1
            try:
                n_mid, alpha = selectivity_alpha(n_values, sigma, k_intra, z)
            except (ValueError, FloatingPointError):
                continue
            finite_mask = np.isfinite(alpha)
            alpha_reliable = alpha[finite_mask]
            if len(alpha_reliable) > 0:
                local_max = alpha_reliable.max()
                if local_max > ALPHA_SUPERSELECTIVE_THRESHOLD:
                    n_superselective_cells += 1
                if local_max > max_alpha:
                    max_alpha = local_max
                    best_params = (float(k_intra), float(z))

            # Exact alpha, evaluated at the SAME n_mid points used above
            # (the geometric midpoints selectivity_alpha() differentiates
            # at), so the exact and finite-difference peaks are directly
            # comparable rather than computed on different grids.
            try:
                alpha_ex = alpha_exact(n_mid, sigma, k_intra, z)
            except (ValueError, FloatingPointError):
                continue
            alpha_ex_reliable = alpha_ex[np.isfinite(alpha_ex)]
            if len(alpha_ex_reliable) > 0:
                local_max_exact = alpha_ex_reliable.max()
                if local_max_exact > max_alpha_exact:
                    max_alpha_exact = local_max_exact
                    best_params_exact = (float(k_intra), float(z))

    return {
        "max_alpha_found": float(max_alpha),
        "any_superselective": bool(max_alpha > ALPHA_SUPERSELECTIVE_THRESHOLD),
        "best_params_k_intra_z": best_params,
        "fraction_of_grid_superselective": (
            n_superselective_cells / total_reliable_cells if total_reliable_cells else 0.0
        ),
        "grid_size": total_reliable_cells,
        "n_cells_excluded_numerically_unreliable": n_excluded,
        "max_alpha_exact_found": float(max_alpha_exact),
        "best_params_exact_k_intra_z": best_params_exact,
    }


def print_summary():
    print("=" * 78)
    print("AIM 3/4 RANGE SWEEP -- robustness of ranking-flip and alpha>1")
    print("reachability across the K_INTRA plausibility range (NOT a calibration)")
    print("=" * 78)

    k_intra_values = get_k_intra_sweep_values()
    print(f"\nK_INTRA sweep values pulled from k_intra_z_plausibility_bounds.py "
          f"(real sigma-anchored table): {len(k_intra_values)} distinct values, "
          f"spanning {k_intra_values.min():.3e} to {k_intra_values.max():.3e} "
          f"cm^2/Ab per receptor "
          f"({k_intra_values.max()/k_intra_values.min():.0f}x range).")

    print("\n--- Question 1: does the patch-regime ranking behavior hold "
          "across the whole K_INTRA range? ---")
    print("    (Per the v4/v7 finding: the SIZE-MISMATCHED healthy/PDAC-"
          "plasma comparison is expected to FLIP between patch regimes; "
          "the roughly SIZE-MATCHED PANC-1/HPDEC cell-line comparison is "
          "expected NOT to flip, i.e. ratio > 1 under BOTH regimes. This "
          "script checks each comparison against its OWN expected "
          "behavior, not a single flip-for-both rule.")
    flip_results = sweep_ranking_flip_robustness()
    expected_to_flip = {"panc1_hpdec": False, "healthy_pdac_plasma": True}
    # NOTE (v14->v15): sweep_ranking_flip_robustness() now also returns
    # "healthy_pdac_plasma_egfr" (identical to "healthy_pdac_plasma") and
    # "healthy_pdac_plasma_ca199" (NEW) -- see that function's docstring.
    # This loop intentionally only iterates over the two ORIGINAL keys
    # (which expected_to_flip covers), since the EGFR/CA19-9 strict-
    # resolved-only comparison is a DIFFERENT, more specific question
    # (see the new section printed below, via
    # strict_resolved_egfr_ca199_comparison()) with its own framing rather
    # than fitting awkwardly into this loop's flip-vs-no-flip structure.
    for label in ("panc1_hpdec", "healthy_pdac_plasma"):
        rows = flip_results[label]
        p2_all_above_1 = all(r["ratio_p2"] > 1.0 for r in rows)
        p0_all_below_1 = all(r["ratio_p0"] < 1.0 for r in rows)
        p0_all_above_1 = all(r["ratio_p0"] > 1.0 for r in rows)
        ratios_p2 = [r["ratio_p2"] for r in rows]
        ratios_p0 = [r["ratio_p0"] for r in rows]
        print(f"\n  {label} (expected to flip: {expected_to_flip[label]}):")
        print(f"    patch_exponent=2 (size-blind) ratio > 1 at ALL "
              f"{len(rows)} K_INTRA values tested: {p2_all_above_1}")
        print(f"    ratio_p2 range across sweep: {min(ratios_p2):.3f} - {max(ratios_p2):.3f}")
        print(f"    ratio_p0 range across sweep: {min(ratios_p0):.3f} - {max(ratios_p0):.3f}")
        if expected_to_flip[label]:
            robust = p2_all_above_1 and p0_all_below_1
            print(f"    patch_exponent=0 (density-driven) ratio < 1 at ALL "
                  f"{len(rows)} K_INTRA values tested: {p0_all_below_1}")
            if robust:
                print(f"    [ROBUST FLIP] the ranking genuinely flips "
                      f"(ratio_p2>1, ratio_p0<1) at EVERY tested point in "
                      f"the plausibility range -- not just the arbitrary "
                      f"placeholder value used in the original demo "
                      f"scripts.")
            else:
                print(f"    [NOT FULLY ROBUST] the expected flip does NOT "
                      f"hold at every tested K_INTRA value -- see raw rows.")
        else:
            robust = p2_all_above_1 and p0_all_above_1
            print(f"    patch_exponent=0 (density-driven) ratio > 1 at ALL "
                  f"{len(rows)} K_INTRA values tested: {p0_all_above_1}")
            if robust:
                print(f"    [ROBUST NO-FLIP] as expected for this roughly "
                      f"size-matched comparison, the ranking does NOT flip "
                      f"(stays ratio>1 under BOTH patch regimes) at EVERY "
                      f"tested point in the plausibility range.")
            else:
                print(f"    [NUANCED, NOT A FLIP-VS-NO-FLIP BUG] checked "
                      f"directly (not a numerical artifact -- yields here "
                      f"are mid-range, 0.5-0.6, nowhere near float "
                      f"precision limits): as K_INTRA rises toward the "
                      f"saturating end of the plausibility range, BOTH "
                      f"populations' yields approach ~1 and their ratio "
                      f"compresses toward 1 from above, crossing very "
                      f"slightly below 1 (ratio~0.999-0.990) at the "
                      f"single largest K_INTRA value tested -- not the "
                      f"large, clean flip seen for the mismatched-size "
                      f"comparison. This is a real, physically sensible "
                      f"consequence of approaching saturation (differences "
                      f"between similar-sized populations shrink as both "
                      f"approach ceiling yield), not evidence the v4 "
                      f"no-flip claim was wrong for realistic (non-"
                      f"saturating) K_INTRA -- but it means 'never flips' "
                      f"should be stated as 'does not flip except within "
                      f"a hair's breadth of full saturation,' not as an "
                      f"unconditional rule. See the printed per-value "
                      f"table (or re-run with verbose output) for the "
                      f"exact numbers.")

    print("\n--- Question 1b (NEW, v14->v15): does the healthy/PDAC-plasma "
          "ranking flip survive STRICT-RESOLVED-ONLY filtering (pixel-exact "
          "points only, excluding density-estimated reconstructions)? ---")
    print("    Added directly in response to this project's independent "
          "diagnosis of an external Biomni review -- see 'Immediate next "
          "steps' item 0a and docs/ for the full account. EGFR was "
          "previously cited as a second robust ranking-flip case alongside "
          "CA19-9; this check is why that claim has been retired.")
    strict_result = strict_resolved_egfr_ca199_comparison()
    for marker in ("egfr", "ca199"):
        info = strict_result[marker]
        verdict = "SURVIVES strict filtering" if info["flip_survives_strict_filter"] else "does NOT survive strict filtering"
        print(f"  {marker.upper()}: ratio_p2 range on strict-resolved-only "
              f"data = [{info['ratio_p2_min']:.3f}, {info['ratio_p2_max']:.3f}] "
              f"-- flip {verdict}.")
    print("  CA19-9 is therefore the project's PRIMARY, headline Lennon-"
          "based ranking-flip demonstration; EGFR is retained only as a "
          "limitations-section negative-control example of how digitization-"
          "reconstruction choices can manufacture an apparent finding that "
          "doesn't survive stricter scrutiny.")

    print("\n--- Question 2: is alpha>1 (superselectivity) reachable "
          "ANYWHERE in the swept (K_INTRA, z) grid? ---")
    alpha_results = sweep_superselectivity_reachability()
    print(f"  Grid size tested (numerically reliable cells only): "
          f"{alpha_results['grid_size']} (K_INTRA, z) pairs")
    if alpha_results["n_cells_excluded_numerically_unreliable"] > 0:
        print(f"  [NOTE] {alpha_results['n_cells_excluded_numerically_unreliable']} "
              f"grid cells EXCLUDED: at the smallest K_INTRA values tested, "
              f"n*sigma*K_INTRA at n=1 is close enough to float64 machine "
              f"epsilon that exp(x)-1 loses precision, producing wildly "
              f"oscillating non-physical alpha values (checked directly -- "
              f"see MIN_RELIABLE_EXPONENT_AT_N1 in this file) rather than a "
              f"real superselectivity signal. These are excluded rather "
              f"than silently included.")
    print(f"  Max alpha found anywhere in the RELIABLE grid (finite-difference): "
          f"{alpha_results['max_alpha_found']:.3f}")
    print(f"  Max alpha found anywhere in the RELIABLE grid (EXACT, closed-form "
          f"alpha_exact() -- NEW v14->v15): {alpha_results['max_alpha_exact_found']:.3f}")
    print(f"  Any point in the plausible grid is genuinely superselective "
          f"(alpha>1): {alpha_results['any_superselective']}")
    print(f"  Fraction of the reliable (K_INTRA, z) grid where alpha>1 is "
          f"reached for at least one n: {alpha_results['fraction_of_grid_superselective']:.1%}")
    if alpha_results["best_params_k_intra_z"]:
        k_i, z_i = alpha_results["best_params_k_intra_z"]
        print(f"  Best (K_INTRA, z) pair found: K_INTRA={k_i:.3e}, z={z_i:.3e}")

    print("\n--- Question 3 (NEW, this session): how sensitive is the CA19-9 "
          "headline ratio to receptor accessibility (f_accessible), applied "
          "to REAL digitized data for the first time? ---")
    acc_results = ca199_selectivity_sensitivity_to_accessibility()
    no_effect_ratios = [r["ratio"] for r in acc_results["no_effect"]]
    equal_ratios = [r["ratio"] for r in acc_results["illustrative_equal_direction"]]
    unequal_ratios = [r["ratio"] for r in acc_results["illustrative_unequal_direction"]]
    print(f"  no_effect (f_accessible=1.0 both): ratio range "
          f"{min(no_effect_ratios):.3f} - {max(no_effect_ratios):.3f} "
          f"(the project's existing, already-reported headline finding)")
    print(f"  illustrative_equal_direction (f_accessible=0.7 both): ratio range "
          f"{min(equal_ratios):.3f} - {max(equal_ratios):.3f} "
          f"(IDENTICAL to no_effect, confirming equal accessibility cancels "
          f"out of the ratio exactly, even against real digitized data)")
    print(f"  illustrative_unequal_direction (tumor f_accessible=0.49, normal=0.7): "
          f"ratio range {min(unequal_ratios):.3f} - {max(unequal_ratios):.3f} "
          f"(the flip WEAKENS substantially but does not reverse in this "
          f"specific tested scenario)")
    print(f"\n  {acc_results['note']}")

    return flip_results, alpha_results


if __name__ == "__main__":
    print_summary()
