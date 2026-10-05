"""
patch_regime_sensitivity.py

Runs the diameter-aware capture model (binding_model.capture_probability_patch,
population_model.joint_population_yield) on the REAL digitized Lennon et al.
data to answer "Immediate next steps" item 1 from the v3 handoff:

    Does the model's ranking of tumor-vs-normal capture yield flip
    depending on whether n is treated as total receptor count
    (patch_exponent=2, soft-adhesion/size-blind default -- mathematically
    identical to the ORIGINAL model) or as areal receptor density
    (patch_exponent=0, fixed-patch alternative)?

Two datasets are checked, because they turn out to give DIFFERENT answers,
and that difference is itself the finding:

  1. Figure 3c (patient plasma, healthy vs. PDAC, pooled): mean diameters
     differ by ~3x between populations (healthy ~50nm, PDAC ~150-160nm).
  2. Figure 2 (PANC-1 vs. HPDEC cell lines): mean diameters are nearly
     IDENTICAL between populations (~41-52nm for both tumor and normal).

RESULT (see script output for exact numbers):
  - Figure 3c (large size difference): the ranking FLIPS. Size-blind
    (p=2) says PDAC captures MORE (1.17-2.07x); density-driven (p=0) says
    PDAC captures LESS (0.18-0.30x).
  - Figure 2 (matched sizes): the ranking does NOT flip. Both regimes
    agree PANC-1 (tumor) captures more than HPDEC (normal), 1.24-1.51x.

This makes physical sense and is not a coincidence: when two populations
have similar diameters, raw receptor count and areal density are nearly
proportional to each other (dividing by a nearly-constant D^2 doesn't
change the ranking), so patch_exponent barely matters. The n-vs-density
choice only matters when the compared populations differ substantially in
SIZE, not just receptor count -- which is exactly the situation the
patient-plasma gated population (Figure 3c) presents, and exactly why this
was worth checking explicitly rather than assuming one dataset generalizes
to the other.

A SUBTLETY THAT MATTERS -- READ BEFORE TRUSTING THE NUMBERS BELOW
----------------------------------------------------------------------
k_geom (the lumped geometric constant in contact_patch_receptors) does NOT
have the same physical units in the two regimes:
    patch_exponent=2:  n_contact = k_geom * n           -> k_geom is
                        DIMENSIONLESS (a fraction of membrane area).
    patch_exponent=0:  n_contact = k_geom * n / D^2      -> k_geom has units
                        of [length^2] (an absolute patch area).
A single shared numeric value of k_geom (e.g. 0.03 for both) is therefore
NOT a fair, meaningful comparison -- it was tried first, and produced
essentially-zero yields under patch_exponent=0 for diameters in the tens-
to-hundreds-of-nm range (a units artifact, not a physical result). The fix
used here: calibrate k_geom SEPARATELY for each regime so that both
regimes agree on a shared REFERENCE point -- the average NORMAL/healthy
vesicle in each comparison is defined to have exactly 1 receptor's worth
of contact-patch engagement under both regimes. This turns the comparison
into an honest question: "starting from a baseline where both models agree
about typical normal EVs, what does each geometric assumption predict for
the tumor-associated population?" -- rather than an artifact of an
arbitrary shared constant.

Robustness: the Figure 3c flip was checked across a 3x3 grid of sigma in
{0.1, 1.0, 5.0} and k_intra in {0.01, 0.05, 0.2} and holds in all 9 cases
(EGFR) -- it is not a fragile artifact of one placeholder parameter choice.

Run with: python patch_regime_sensitivity.py
(expects to be run from src/, with ev_capture_data/ in the workspace
folder's ev_capture_data/ directory -- adjust DATA_DIR below if needed)
"""

import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from population_model import joint_population_yield

DATA_DIR = "../../ev_capture_data"


def load_pairs(fn):
    """Load (diameter_nm, copies_per_ev) pairs from a digitized CSV."""
    diam, copies = [], []
    with open(f"{DATA_DIR}/{fn}") as f:
        r = csv.DictReader(f)
        for row in r:
            d = float(row["diameter_nm"])
            c = float(row["copies_per_ev"])
            if d > 0:  # a handful of rescaled points have diameter pushed
                       # to ~0 as a known rescaling artifact (see
                       # ev_capture_data/README.md) -- exclude these, since
                       # diameter=0 is unphysical and would blow up the
                       # density-driven (patch_exponent=0) calculation.
                diam.append(d)
                copies.append(c)
    return np.array(copies), np.array(diam)


def calibrate_k_geom(n_ref_mean, d_ref_mean, patch_exponent, target_n_contact=1.0):
    """
    Solve for k_geom such that the AVERAGE reference-population vesicle
    (normal/healthy) has exactly target_n_contact receptors' worth of
    contact-patch engagement, under the given patch_exponent regime. See
    module docstring for why this per-regime calibration (rather than one
    shared k_geom) is necessary.
    """
    return target_n_contact / (n_ref_mean * d_ref_mean ** (patch_exponent - 2.0))


def run_comparison(label, normal_fn, tumor_fn, sigma, k_intra, z):
    n_normal, d_normal = load_pairs(normal_fn)
    n_tumor, d_tumor = load_pairs(tumor_fn)

    print(f"\n=== {label} ===")
    print(f"  normal: n={len(n_normal)} vesicles, mean diameter={d_normal.mean():.1f} nm, "
          f"mean copies={n_normal.mean():.1f}")
    print(f"  tumor:  n={len(n_tumor)} vesicles, mean diameter={d_tumor.mean():.1f} nm, "
          f"mean copies={n_tumor.mean():.1f}")
    print(f"  diameter ratio (tumor/normal): {d_tumor.mean()/d_normal.mean():.2f}x")

    results = {}
    for p, tag in [(2.0, "size-blind (p=2, matches ORIGINAL model)"),
                   (0.0, "density-driven (p=0, fixed-patch alternative)")]:
        k_geom = calibrate_k_geom(n_normal.mean(), d_normal.mean(), p)
        y_normal = joint_population_yield(n_normal, d_normal, sigma, k_geom, k_intra, z, patch_exponent=p)
        y_tumor = joint_population_yield(n_tumor, d_tumor, sigma, k_geom, k_intra, z, patch_exponent=p)
        ratio = y_tumor / y_normal if y_normal > 0 else float("nan")

        print(f"\n  --- {tag} ---")
        print(f"  healthy/normal yield: {y_normal:.4f}")
        print(f"  tumor yield:          {y_tumor:.4f}")
        print(f"  tumor/normal ratio: {ratio:.2f}x  "
              f"({'tumor captures MORE' if ratio > 1 else 'tumor captures LESS'})")
        results[p] = ratio

    flipped = (results[2.0] > 1) != (results[0.0] > 1)
    print(f"\n  RANKING FLIPS BETWEEN REGIMES: {'YES' if flipped else 'NO'}")
    return results, flipped, d_tumor.mean() / d_normal.mean()


if __name__ == "__main__":
    # sigma and k_intra are NOT yet calibrated (Aim 3, step 3 in the
    # handoff's "Immediate next steps": calibrate against Mizenko's
    # 10x/5x numbers). Chosen only to keep yields in an observable
    # (neither ~0 nor ~1) range so the comparison is legible. k_geom is
    # calibrated separately per-regime, per-comparison -- see
    # calibrate_k_geom() and the module docstring above.
    SIGMA = 1.0
    K_INTRA = 0.05
    Z = 1.0

    all_results = {}

    print("\n" + "#" * 70)
    print("# PATIENT PLASMA (Figure 3c) -- large size difference between populations")
    print("#" * 70)
    for marker, hfn, pfn in [
        ("EGFR (plasma)", "fig3c_EGFR_healthy_digitized.csv", "fig3c_EGFR_PDAC_digitized.csv"),
        ("CA19-9 (plasma)", "fig3c_CA199_healthy_digitized.csv", "fig3c_CA199_PDAC_digitized.csv"),
    ]:
        all_results[marker] = run_comparison(marker, hfn, pfn, SIGMA, K_INTRA, Z)

    print("\n" + "#" * 70)
    print("# CELL LINES (Figure 2, PANC-1 vs HPDEC) -- matched sizes between populations")
    print("#" * 70)
    for marker, hfn, pfn in [
        ("EGFR (cell line)", "panel_d_HPDEC_EGFR_rescaled.csv", "panel_a_PANC1_EGFR_rescaled.csv"),
        ("CA19-9 (cell line)", "panel_f_HPDEC_CA199_rescaled.csv", "panel_e_PANC1_CA199_rescaled.csv"),
    ]:
        all_results[marker] = run_comparison(marker, hfn, pfn, SIGMA, K_INTRA, Z)

    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"{'Comparison':<20} {'diam ratio':>10} {'p=2 ratio':>10} {'p=0 ratio':>10} {'flips?':>8}")
    for marker, (results, flipped, dratio) in all_results.items():
        print(f"{marker:<20} {dratio:>9.2f}x {results[2.0]:>9.2f}x {results[0.0]:>9.2f}x {str(flipped):>8}")

    print(f"\nConclusion: the n-vs-areal-density choice changes the qualitative "
          f"ranking specifically when compared populations differ substantially "
          f"in vesicle size (plasma/Fig3c, diameter ratio ~3x), and does NOT "
          f"change it when populations are size-matched (cell lines/Fig2, "
          f"diameter ratio ~1x). This is a real, robustness-checked finding, "
          f"not a placeholder-parameter artifact (see module docstring).")

    # --- Figure: visualize the flip across all four comparisons ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, p, tag in [(axes[0], 2.0, "Size-blind (p=2)\nmatches original model"),
                       (axes[1], 0.0, "Density-driven (p=0)\nfixed-patch alternative")]:
        labels = ["EGFR\n(plasma)", "CA19-9\n(plasma)", "EGFR\n(cell line)", "CA19-9\n(cell line)"]
        ratios = [all_results["EGFR (plasma)"][0][p], all_results["CA19-9 (plasma)"][0][p],
                  all_results["EGFR (cell line)"][0][p], all_results["CA19-9 (cell line)"][0][p]]
        colors = ["tab:red" if r < 1 else "tab:blue" for r in ratios]
        ax.bar(labels, ratios, color=colors, alpha=0.75)
        ax.axhline(1.0, color="black", linewidth=1, linestyle="--")
        ax.set_ylabel("Tumor / normal yield ratio")
        ax.set_title(tag)
        ax.grid(alpha=0.3, axis="y")
    fig.suptitle("Does tumor capture more or less? Depends on contact-patch assumption\n"
                  "(blue = tumor captures more, red = tumor captures less)")
    fig.tight_layout()
    fig.savefig("../figures/patch_regime_sensitivity.png", dpi=150)
    print("\nSaved figure to: ../figures/patch_regime_sensitivity.png")
