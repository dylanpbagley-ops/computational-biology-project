"""
mizenko_cooccurrence_scoping.py

WHY THIS FILE EXISTS
---------------------
This is NOT a validation script, NOT a calibration, and NOT part of Aim 2's
"reproduce external signatures" leg. It is exploratory SCOPING work for
Aim 3b (multi-marker joint distributions / the independence question),
done to answer one specific question before deciding how to build Aim 3b's
real model: "how strong would tetraspanin/tumor-marker anti-correlation
need to be, physically, to produce Mizenko et al.'s observed ~10x (EpCAM)
and ~5x (Her2) fold-differences between single-tetraspanin capture and
non-specific (biotin) capture?"

WHY THIS QUESTION MATTERS -- THE LOGIC ERROR THIS SESSION CAUGHT
-------------------------------------------------------------------
The v5 handoff's "Immediate next steps" item 1 said to calibrate
K_INTRA, z, and k_geom (binding_model.py's single-ligand-pair binding
parameters) against Mizenko's 10x/5x numbers. Before doing that, this
session re-read Mizenko et al. 2021's actual methodology directly (RSC/
BioMed Central full text, not just the handoff's summary) and found this
doesn't work, for a specific, checkable reason:

Mizenko's 10x/5x numbers do NOT measure how efficiently one ligand
(tetraspanin) captures EVs bearing that ligand. They measure something
different: EVs were captured by TWO DIFFERENT METHODS from the SAME
starting population -- (a) single-tetraspanin capture (anti-CD9, -CD63, or
-CD81 coated chip) and (b) genuinely non-specific capture (EVs chemically
biotinylated at any surface amine, captured on an anti-biotin chip,
independent of any particular protein) -- and THEN both captured
populations were stained for EpCAM and Her2. The paper's own words: "EpCAM+
EVs could be detected on all three tetraspanin capture spots at similar
rates, [but] biotin capture identified more EVs by approximately a
tenfold difference" -- and their own interpretation: "there may be a large
percentage of non-tetraspanin expressing EVs that contain the bulk of
these markers."

This is a CO-OCCURRENCE / CORRELATION measurement between TWO DIFFERENT
surface markers (tetraspanin vs. EpCAM/Her2) across the same EV
population -- not a measurement of how strongly ONE marker binds its own
ligand. binding_model.py's K_INTRA, z, and k_geom describe single
ligand-receptor pair thermodynamics (how likely IS a vesicle with n copies
of marker X to bind a surface coated with ligand for marker X). There is
no principled mapping from a cross-marker co-occurrence ratio onto those
single-marker parameters -- forcing a fit would produce numbers that don't
mean what a calibration claims to mean. This was surfaced directly and the
user chose (see handoff) to reassign Mizenko's numbers to Aim 3b (the
multi-marker/independence question, where the handoff had ALREADY
correctly flagged co-occurrence and independence as central) rather than
force a single-marker calibration.

WHAT THIS SCRIPT DOES
------------------------
A quick, honest scoping exercise (NOT a fitted model, NOT using real
per-vesicle EpCAM/tetraspanin joint data -- none exists for this
comparison) to answer: "is a 5-10x co-occurrence deficit even PHYSICALLY
PLAUSIBLE from a simple two-subpopulation mixture model, and if so, how
strong does the anti-correlation between tetraspanin and EpCAM/Her2 need
to be?" This matters because if even extreme anti-correlation assumptions
couldn't reach 5-10x, that would suggest either a different mechanism is
needed, or Mizenko's numbers reflect something beyond simple population
mixture structure (e.g., a detection/normalization artifact worth
re-reading the paper for). Both possibilities are worth knowing BEFORE
building a real Aim 3b model.

THE FINDING
--------------
Modeled as a two-subpopulation mixture: a "tetraspanin-rich" fraction
(captured by single-tetraspanin methods) and a "tetraspanin-poor" fraction
(effectively not captured by single-tetraspanin methods, only by
non-specific/biotin capture), with EpCAM copy number drawn from different
zero-inflated negative binomial distributions in each subpopulation.

Key result, swept across several parameter combinations (see
run_scoping_sweep() below): getting anywhere near a 5-10x ratio in
EpCAM-positive fraction between non-specific and tetraspanin capture
requires BOTH (a) a fairly small tetraspanin-rich fraction (roughly
10-30% of the population) AND (b) a strong, specific suppression of EpCAM
in the tetraspanin-rich subpopulation (roughly 50-100x lower mean EpCAM
copy number in tet-rich vesicles than in tet-poor vesicles) -- mild or
moderate anti-correlation (a Gaussian-copula-style correlation coefficient
of even 0.5, or a 5-10x EpCAM suppression) only produces ratios in the
1.0-1.6x range, nowhere near Mizenko's reported 5-10x. This is a real,
useful finding for scoping Aim 3b: the Mizenko result implies a STRUCTURAL,
near-exclusive relationship between tetraspanin-positive and EpCAM/Her2-
positive status (something closer to two largely-distinct subpopulations
than a graded correlation) -- consistent with, and a quantitative
sharpening of, the paper's own "different EVs carry different markers"
interpretation and the Lennon non-uniform-EGFR-distribution finding
already noted elsewhere in this project as evidence for clustering.

WHAT THIS IS NOT
-------------------
This is not calibrated against any real joint (tetraspanin, EpCAM)
per-vesicle dataset -- none was available this session. It's a plausibility
sweep over an assumed mixture-model FORM, meant to bound what Aim 3b's
real model needs to explain, not a fitted or validated result. Treat the
specific numbers (10-30% tet-rich fraction, 50-100x suppression) as
"this is the ballpark of what would be needed," not as claims about the
real biology.
"""

import numpy as np
from scipy import stats


def epcam_positive_fraction_ratio(
    frac_tet_rich,
    epcam_mean_tet_poor,
    epcam_mean_tet_rich,
    tet_mean_rich=18.0,
    epcam_dispersion=2.0,
    n_vesicles=300_000,
    seed=0,
):
    """
    Simulate a two-subpopulation mixture (tetraspanin-rich vs. -poor) with
    different EpCAM copy-number distributions in each, and compute the
    ratio of EpCAM-positive fraction between a "non-specific capture"
    (uniform sample of the whole population) and a "tetraspanin capture"
    (only the tet-rich subpopulation, an idealized hard-threshold capture)
    -- analogous to the biotin/tetraspanin comparison in Mizenko et al.

    Parameters
    ----------
    frac_tet_rich : float
        Fraction of the population that is "tetraspanin-rich" (i.e., would
        be captured by single-tetraspanin methods). The rest are treated
        as having zero tetraspanin (not captured by tetraspanin methods at
        all -- an idealized hard threshold, not claiming real capture is
        this sharp).
    epcam_mean_tet_poor, epcam_mean_tet_rich : float
        Mean EpCAM copy number in the tet-poor and tet-rich subpopulations
        respectively. Set these unequal to create anti-correlation between
        tetraspanin and EpCAM status.
    tet_mean_rich : float
        Mean tetraspanin copy number within the tet-rich subpopulation
        (not otherwise used in this simplified capture model beyond
        defining who's "tet-rich").
    epcam_dispersion : float
        Negative-binomial dispersion parameter for EpCAM counts (same
        convention as population_model.negative_binomial_distribution).
    n_vesicles : int
        Simulated population size.
    seed : int
        RNG seed for reproducibility.

    Returns
    -------
    dict with keys:
        "frac_tet_rich", "P_EpCAM_pos_given_captured" (tetraspanin capture),
        "P_EpCAM_pos_overall" (non-specific/biotin capture), "ratio"
        (overall / captured -- analogous to Mizenko's biotin/tetraspanin
        fold-difference).
    """
    rng = np.random.default_rng(seed)
    is_tet_rich = rng.random(n_vesicles) < frac_tet_rich

    r = epcam_dispersion
    epcam_poor = stats.nbinom.rvs(
        r, r / (r + epcam_mean_tet_poor), size=n_vesicles, random_state=rng
    )
    epcam_rich = stats.nbinom.rvs(
        r, r / (r + max(epcam_mean_tet_rich, 1e-6)), size=n_vesicles, random_state=rng
    )
    epcam = np.where(is_tet_rich, epcam_rich, epcam_poor)
    epcam_positive = epcam > 0

    captured_by_tetraspanin = is_tet_rich
    p_epcam_given_captured = epcam_positive[captured_by_tetraspanin].mean()
    p_epcam_overall = epcam_positive.mean()
    ratio = p_epcam_overall / p_epcam_given_captured

    return {
        "frac_tet_rich": frac_tet_rich,
        "P_EpCAM_pos_given_tet_capture": p_epcam_given_captured,
        "P_EpCAM_pos_overall": p_epcam_overall,
        "ratio": ratio,
    }


def run_scoping_sweep():
    """
    Reproduce the sweep reported in the module docstring: how does the
    biotin/tetraspanin EpCAM-positive-fraction ratio depend on (a) how
    small the tetraspanin-rich fraction is, and (b) how strongly EpCAM is
    suppressed in the tetraspanin-rich subpopulation?
    """
    print(f"{'frac_tet_rich':>14} | {'EpCAM mean (poor/rich)':>24} | "
          f"{'P(EpCAM+|tet-cap)':>18} | {'P(EpCAM+|overall)':>18} | {'ratio':>7}")
    print("-" * 95)

    configs = [
        (0.5, 15, 5),
        (0.2, 15, 5),
        (0.2, 15, 1),
        (0.1, 15, 0.5),
        (0.1, 20, 0.2),
        (0.3, 12, 0.1),
    ]
    for frac_tet_rich, epcam_poor, epcam_rich in configs:
        r = epcam_positive_fraction_ratio(frac_tet_rich, epcam_poor, epcam_rich)
        print(f"{frac_tet_rich:>14.2f} | {f'{epcam_poor}/{epcam_rich}':>24} | "
              f"{r['P_EpCAM_pos_given_tet_capture']:>18.3f} | "
              f"{r['P_EpCAM_pos_overall']:>18.3f} | {r['ratio']:>7.2f}x")

    print("\nTAKEAWAY: reaching Mizenko's reported ~10x (EpCAM) requires BOTH a "
          "fairly small tetraspanin-rich fraction (~10-30%) AND strong "
          "(~50-100x) suppression of EpCAM specifically within the "
          "tetraspanin-rich subpopulation -- mild/moderate anti-correlation "
          "alone (the first two rows above) falls far short, only reaching "
          "~1.0-1.2x. This is a real constraint on what Aim 3b's real "
          "co-occurrence model needs to reproduce, not a claim about the "
          "true underlying biology (no real joint tetraspanin/EpCAM "
          "per-vesicle dataset was used here -- see module docstring).")


if __name__ == "__main__":
    print(__doc__.strip().split("\n\n")[0])
    print()
    run_scoping_sweep()
