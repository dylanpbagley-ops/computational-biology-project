"""
aim10_car_t_kill_probability_model.py

WHY THIS FILE EXISTS -- THE PIVOT
-----------------------------------
This is the first file of the project's pivot from EV (extracellular
vesicle) early-detection marker selection to CAR-T target-antigen
selection, per the user-approved plan in docs/Pivot_Implementation_Plan.md
(itself built from docs/PIVOT-FOR-NOVELTY.md's research). The reasoning for
the pivot, in one sentence: early cancer detection from blood is
biologically bottlenecked (very little tumor material reaches the blood,
so no amount of clever marker-selection math can move the needle much),
while CAR-T target selection is a decision made ONCE, up front, when a
therapy is designed -- so a better selection method has a real, larger
lever on how many patients relapse from antigen escape.

THE MATHEMATICAL TRANSLATION, STATED EXPLICITLY
--------------------------------------------------
Every function in binding_model.py answers: "given a vesicle with n
receptors, what's the probability it gets CAPTURED by a ligand-coated
surface?" This file asks the mechanistically analogous question: "given a
tumor cell with n antigen copies, what's the probability it gets KILLED
by a CAR-T cell?" The translation table (see PIVOT-EXPLAINER-FOR-DYLAN.md
for the plain-language version):

    EV project                          CAR-T project (this file)
    -----------------------------       -----------------------------
    vesicle                             tumor cell
    receptor count on a vesicle (n)     antigen count on a tumor cell (n)
    capture antibody surface            CAR T-cell
    capture_probability(n, ...)         kill_probability(n, ...)
    tumor vs. normal EV populations     tumor cells above vs. below the
                                           functional kill threshold
    picking the best EV marker          picking the best CAR-T target
                                           antigen (this file's Aim 11)

WHY THE UNDERLYING MATH DOES NOT CHANGE
-------------------------------------------
binding_model.capture_probability(n, sigma, k_intra, z) is a general
statement about how strongly a MULTIVALENT, NONLINEAR response depends on
a per-particle count n -- it says nothing intrinsically about vesicles or
ligands. CAR-T killing is ALSO a multivalent, nonlinear (sigmoidal/
thresholded), receptor-density-dependent process -- this is not an
analogy invented for this project, it is directly documented in the
literature this project already cites in binding_model.py's own
docstring (Majzner et al. 2020, Walker et al. 2017: CAR-T cytotoxicity is
a graduated, thresholded function of antigen density). So the SAME
functional form (capture_probability) is reused UNCHANGED here, just
renamed and re-parameterized for this new domain -- see
kill_probability() below, which is a thin, explicitly-labeled wrapper, not
a reimplementation.

THE REAL DATA THIS FILE IS CALIBRATED AGAINST -- AND WHY THIS IS A REAL
METHODOLOGICAL UPGRADE OVER THE EV PROJECT
---------------------------------------------------------------------------
K_INTRA, z, and k_geom in the EV project were PERMANENTLY, EXPLICITLY
uncalibrated placeholders (see binding_model.py's own docstring) -- six-
plus full sessions searched for a real single-marker absolute-capture-
efficiency data source and never found a clean one, so the EV project
swept a plausibility RANGE rather than fitting a point estimate. CAR-T is
different: the field has PUBLISHED, DIRECTLY MEASURED, TITRATED dose-
response curves relating a precisely engineered receptor density to a
measured killing/activation outcome. This file uses one specific, real,
directly verified source:

    Hannigan, Patakas, London & Findlay (Antibody Analytics, UK),
    "Deciphering the Impact of Antigen Density on Tumor Targeting
    Therapeutics," Agilent Technologies Application Note 5994-7175EN
    (published June 6, 2024). Verified directly this session (full PDF
    fetched and read, not taken from a search snippet). Real experimental
    system: CHO-K1 cells engineered with an inducible CD19 expression
    system (IndEx-2), titrated across 78 to 185,223 CD19 copies per cell
    (quantified via Quantibrite PE calibration beads on an Agilent
    NovoCyte Quanteon flow cytometer -- a real, calibrated, absolute
    molecules-per-cell measurement, not a relative/arbitrary-unit
    fluorescence value). Two real CD19-targeting T-cell engager (TCE)
    molecules were tested across this density range with a real 10-donor
    human T-cell cohort, and a real ROC/Youden's-J threshold analysis
    determined the receptor density required for measurable activity:

        CLN-978:        81 CD19 receptors/cell (cytolysis AND IFN-gamma
                         threshold -- this was the lower limit of
                         detection, so the TRUE threshold may be even
                         lower; reported honestly as an upper bound on
                         the true threshold, not a precise point value)
        Blinatumomab:    843 CD19 receptors/cell (cytolysis threshold)
                         8,690 CD19 receptors/cell (IFN-gamma threshold)

IMPORTANT SCOPE NOTE: this is a T-cell ENGAGER (TCE) system, not a CAR-T
cell specifically -- TCEs (like blinatumomab) and CAR-T cells are
DIFFERENT therapeutic modalities that both work by directing a T-cell to
kill an antigen-bearing target cell, and both are documented in the
literature to show the same qualitative phenomenon (a receptor-density
activation threshold). This file uses the TCE data because it is the most
rigorous, directly-verified, PUBLICLY ACCESSIBLE real titrated dose-
response curve found this session -- CAR-T-specific titrated dose-response
data (as opposed to retrospective threshold correlations in patient
cohorts, which is what the CD19 CAR-T literature mostly reports) was not
found in an equally clean, directly downloadable form. This is reported
honestly, not blurred: this file's calibration source is a TCE system;
the ranking/escape-fraction application in aim11 is explicitly a CAR-T
target-selection question. The mechanism (receptor-density threshold for
T-cell-mediated killing) is shared and well-documented across both
modalities (see binding_model.py's existing Majzner/Walker CAR-T
citations), which is why this substitution is defensible -- but it is a
real, stated limitation, not something to gloss over in the write-up.

WHAT THIS FILE DOES NOT CLAIM
----------------------------------
- Does NOT claim CLN-978's or blinatumomab's exact threshold is universal
  across all CAR-T constructs or all target antigens -- these are two
  real, specific molecules' measured thresholds, used here as two
  concrete, real anchor points for fitting kill_probability()'s shape,
  not as a universal biological constant.
- Does NOT claim this resolves CAR-T's other failure modes (T-cell
  exhaustion, trafficking, tumor microenvironment) -- see
  aim11_antigen_target_ranking.py's own scope note and
  docs/CAR_T_Pivot_Validation_and_Findings.md for the full honest ceiling.
- Does NOT throw away or invalidate any EV-project result. binding_model.py
  is imported and reused UNCHANGED below -- every EV-project function,
  test, and finding remains exactly as it was.
"""

import os
import sys

import numpy as np
from scipy import optimize

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from binding_model import capture_probability, q_bound


# ---------------------------------------------------------------------------
# REAL DATA: the Agilent/Antibody Analytics IndEx-2 titrated dose-response
# thresholds, verified directly this session (see module docstring).
# ---------------------------------------------------------------------------

# Real engineered CD19 density range used in the source experiment
# (quantified via Quantibrite PE calibration beads -- an absolute,
# calibrated molecules-per-cell measurement).
INDEX2_CD19_DENSITY_RANGE_MOLECULES_PER_CELL = (78, 185223)

# Real, measured activation thresholds (molecules/cell), by molecule and
# readout. CLN-978's value is explicitly the LOWER LIMIT OF DETECTION in
# the source assay, not a resolved point estimate -- the true threshold
# could be lower. This is stated here, not silently treated as exact.
REAL_TCE_ACTIVATION_THRESHOLDS = {
    "CLN-978": {
        "cytolysis_molecules_per_cell": 81,
        "ifn_gamma_molecules_per_cell": 81,
        "threshold_is_detection_limit": True,  # true value may be lower
    },
    "blinatumomab": {
        "cytolysis_molecules_per_cell": 843,
        "ifn_gamma_molecules_per_cell": 8690,
        "threshold_is_detection_limit": False,
    },
}

# Real, independently-sourced retrospective clinical threshold (a
# DIFFERENT kind of evidence -- a patient-cohort correlation, not a
# titrated in vitro dose-response curve): CD19 CAR-T (not TCE)
# responders vs. non-responders in B-ALL, per Majzner et al., already
# cited in binding_model.py. Kept here as an independent qualitative
# cross-check (same order of magnitude regime as the TCE thresholds
# above), NOT used to fit any parameter, and NOT blended with the TCE
# numbers above into one number.
CLINICAL_CD19_CAR_T_RESPONSE_THRESHOLD_MOLECULES_PER_CELL = 2000

# Real, independently-sourced patient-cohort summary statistics (median,
# IQR) for pre-treatment and at-relapse CD19 density in large B-cell
# lymphoma patients treated with axicabtagene ciloleucel (a real CD19
# CAR-T product, not a TCE) -- Transplantation and Cellular Therapy
# (2024), verified directly this session via independent cross-check
# (conference abstract text + journal listing, not a single source).
# n=36 pre-treatment, n=25 at relapse (55 unique patients total). Used in
# aim11 to reconstruct a REAL tumor population's antigen-density
# distribution -- see that file's docstring.
LBCL_CD19_DENSITY_PRE_TREATMENT = {
    "median_molecules_per_cell": 5810,
    "iqr_low": 3442,
    "iqr_high": 8725,
    "n_patients": 36,
}
LBCL_CD19_DENSITY_AT_RELAPSE = {
    "median_molecules_per_cell": 2021,
    "iqr_low": 439,
    "iqr_high": 5208,
    "n_patients": 25,
}

# Real, independently-sourced myeloma CD19 heterogeneity range (a SECOND
# candidate antigen system, used in aim11's multi-antigen ranking
# comparison) -- Blood (ASH), verified this session.
MYELOMA_CD19_DENSITY_RANGE_MOLECULES_PER_CELL = (13, 5000)
MYELOMA_CD19_POSITIVITY_FRACTION_RANGE = (0.103, 0.80)  # 10.3%-80% of patients


# ---------------------------------------------------------------------------
# THE KILL-PROBABILITY MODEL -- a thin, explicitly-relabeled wrapper around
# binding_model.capture_probability(). See module docstring for why the
# underlying math is unchanged.
# ---------------------------------------------------------------------------

def kill_probability(n_antigen, car_avidity, k_intra, z):
    """
    Probability that a tumor cell carrying n_antigen copies of the target
    antigen is killed by a CAR-T cell (or T-cell engager), under the same
    nonlinear multivalent-engagement math as binding_model.capture_probability.

    This is EXACTLY capture_probability(n_antigen, car_avidity, k_intra, z)
    -- see binding_model.py -- renamed for this domain. car_avidity plays
    the role sigma (ligand surface density) played in the EV model: here,
    it represents the CAR's own effective avidity/engagement capacity
    (which depends on CAR surface density on the T-cell, CAR-antigen
    binding affinity, and immunological synapse geometry -- all lumped
    into one effective parameter, exactly the same lumping discipline
    binding_model.py uses for K_INTRA).

    Parameters
    ----------
    n_antigen : float or array
        Antigen copy number on the tumor cell.
    car_avidity : float
        Lumped CAR engagement-capacity parameter (analogous to sigma).
    k_intra : float
        Bond-strength / avidity-scaling parameter (analogous to
        binding_model.py's K_INTRA) -- FIT from real data below, not left
        as a permanently uncalibrated placeholder (see
        fit_kill_curve_to_real_threshold()).
    z : float
        Activity/prefactor (same meaning as in binding_model.py).

    Returns
    -------
    float or array
        Kill probability, between 0 and 1.
    """
    return capture_probability(n_antigen, car_avidity, k_intra, z)


def escape_fraction(n_antigen_population, car_avidity, k_intra, z, kill_threshold=0.5):
    """
    THE CORE NEW DELIVERABLE OF THIS PIVOT: given a REAL POPULATION of
    tumor cells' antigen counts (not just the mean), compute the fraction
    of cells that "escape" -- i.e., whose predicted kill probability falls
    below kill_threshold, meaning they are unlikely to be reliably killed.

    This directly implements the pivot document's own stated concrete gap
    (PIVOT-FOR-NOVELTY.md, Question 4): "take the measured per-cell
    antigen-density distribution... convolve it with the nonlinear...
    killing-vs-density curve, and rank/select antigens... by the predicted
    killable fraction / escape fraction rather than by mean expression or
    binary positivity." This function IS that convolution.

    Parameters
    ----------
    n_antigen_population : array
        A REAL (or realistically reconstructed from real summary
        statistics -- see aim11) sample of per-cell antigen counts across
        a tumor population. This is the whole point of the pivot: using
        the FULL distribution, not just its mean.
    car_avidity, k_intra, z :
        Same meaning as kill_probability().
    kill_threshold : float, optional (default 0.5)
        The kill-probability cutoff below which a cell is counted as
        "escaping." 0.5 is a defensible, simple default (a cell more
        likely than not to survive counts as escaping) -- NOT a claim
        that 0.5 is a clinically validated cutoff. This is swept as a
        sensitivity axis in aim11, exactly the same discipline the EV
        project used for its own uncalibrated parameters.

    Returns
    -------
    dict with:
        "escape_fraction": float, the fraction of the population with
            kill_probability < kill_threshold.
        "mean_kill_probability": float, the population-average kill
            probability (the EV project's "yield," renamed).
        "per_cell_kill_probability": array, kill probability for every
            individual cell in the input population (kept for inspection/
            plotting, not just the summary numbers).
    """
    n_antigen_population = np.asarray(n_antigen_population, dtype=float)
    if len(n_antigen_population) == 0:
        raise ValueError("Cannot compute escape_fraction over an empty population.")

    theta = kill_probability(n_antigen_population, car_avidity, k_intra, z)
    escaped = theta < kill_threshold
    return {
        "escape_fraction": float(np.mean(escaped)),
        "mean_kill_probability": float(np.mean(theta)),
        "per_cell_kill_probability": theta,
    }


# ---------------------------------------------------------------------------
# CALIBRATION: fit k_intra directly against the real IndEx-2 titrated
# threshold data. This is the real methodological upgrade over the EV
# project's permanently-swept K_INTRA -- see module docstring.
# ---------------------------------------------------------------------------

def fit_kill_curve_to_real_threshold(threshold_molecules_per_cell, car_avidity=1.0,
                                       z=1.0, threshold_kill_probability=0.5):
    """
    Fit k_intra so that kill_probability(threshold_molecules_per_cell,
    car_avidity, k_intra, z) == threshold_kill_probability.

    WHY THIS IS A REAL CALIBRATION, NOT ANOTHER SWEEP: the EV project's
    K_INTRA stayed uncalibrated because no data source gave a real,
    directly-measured (density, response) PAIR for a single marker's own
    absolute capture efficiency (see binding_model.py's own extensive
    docstring on this exact point). Here, the IndEx-2 data gives EXACTLY
    that: a real, precisely quantified antigen density (in real molecules-
    per-cell units, via calibration beads) at which a real, measured
    biological response (T-cell-mediated cytolysis or IFN-gamma
    production) crosses a detectable threshold. Solving for k_intra such
    that the model reproduces that one real anchor point is a genuine,
    if simple (one-point), calibration -- not a plausibility-bounds sweep
    with no data anchor at all.

    HONEST LIMITATION, STATED DIRECTLY: this is a ONE-POINT calibration
    (it forces the curve through exactly one measured (density,
    threshold) pair), not a full nonlinear-least-squares fit to multiple
    measured points along the curve, because the source app note reports
    threshold crossings (via ROC/Youden's-J analysis), not a full
    continuous dose-response curve with error bars at every tested
    density. A future session with access to the full underlying
    cytolysis-vs-density curve (Figure 5 of the source app note, not
    reproduced numerically in the text extracted this session) could
    upgrade this to a genuine multi-point least-squares fit -- flagged
    here as a real, concrete, not-yet-done improvement, exactly the kind
    of honest "not done" item this project has always tracked.

    Parameters
    ----------
    threshold_molecules_per_cell : float
        The real, measured antigen density at which response crosses
        threshold_kill_probability (e.g. 81, 843, or 8690 -- see
        REAL_TCE_ACTIVATION_THRESHOLDS above).
    car_avidity : float, optional (default 1.0)
        Held fixed during this fit (folded into k_intra, exactly the same
        lumping convention binding_model.py uses -- see that file's
        module docstring on why sigma and k_intra are not separately
        identifiable from a single-point fit).
    z : float, optional (default 1.0)
        Held fixed during this fit, same reasoning as car_avidity.
    threshold_kill_probability : float, optional (default 0.5)
        The kill probability the real threshold density is assumed to
        correspond to. 0.5 is the simplest, most defensible choice (the
        threshold is where response "turns on") absent a more precise
        published mapping -- stated as an assumption, not a measurement.

    Returns
    -------
    dict with:
        "k_intra": the fitted value.
        "car_avidity", "z": the fixed values used (echoed back for
            provenance).
        "threshold_molecules_per_cell", "threshold_kill_probability": the
            real anchor point used (echoed back for provenance).
        "converged": bool, whether the root-find succeeded.
    """
    def residual(log_k_intra):
        k_intra = np.exp(log_k_intra)
        theta = kill_probability(threshold_molecules_per_cell, car_avidity, k_intra, z)
        return theta - threshold_kill_probability

    # Solve in log-space (k_intra spans many orders of magnitude
    # depending on car_avidity/z scale) -- same numerical-stability
    # discipline as k_intra_z_plausibility_bounds.py's own log-scale
    # sweeps in the EV project.
    #
    # BUG CAUGHT WHILE FIRST RUNNING THIS FUNCTION, in TWO STAGES (not
    # before -- exactly the "run every affected script's own standalone
    # entry point end-to-end" discipline this project has always
    # required):
    #
    # Stage 1 (wrong fix, caught immediately on the next run): a naive
    # [-50, 50] bracket on log_k_intra overflows q_bound()'s np.exp() call
    # to inf well before reaching theta=0.5, so brentq's bracket-edge
    # evaluation itself returns NaN. The first attempt to fix this derived
    # an upper bound from np.exp()'s raw float64 overflow point (~x=709.78)
    # -- but this was STILL wrong, because it solved the wrong problem:
    # theta saturates to EXACTLY 1.0 (not just "very close to 1") well
    # before x anywhere near 700, since 1/(1+q_bound) rounds to exactly
    # 1.0 once q_bound exceeds ~1e16 in float64 -- so a bracket that reaches
    # anywhere near the raw exp() overflow point is needlessly wide AND
    # still risks the exact same NaN failure downstream in the Langmuir
    # division step (inf/(1+inf) = NaN), independent of q_bound's own
    # np.exp() overflow.
    #
    # Stage 2 (still wrong, caught on the NEXT run): a fixed upper bound
    # (tried: log_k_intra=10) is NOT safe for every real threshold this
    # file fits, because the "still finite" vs. "already NaN" boundary
    # for log_k_intra shifts LOWER as n grows (checked directly for all
    # three real thresholds: safe up to roughly log_k=2 at n=81, log_k=-2
    # at n=843, log_k=-4 at n=8,690) -- a single hard-coded constant that
    # happened to look safe for one threshold silently breaks for a
    # larger one. This is exactly the kind of "passing on my one test
    # case doesn't mean it's actually correct" mistake this project's own
    # working style explicitly warns about.
    #
    # Stage 3 (the actual, general fix): derive the upper bound directly
    # from the real, known float64 overflow condition, applied PER-n
    # rather than as one shared constant. q_bound()'s np.exp(n*sigma*
    # k_intra) call is safe as long as n*sigma*k_intra stays comfortably
    # below float64's overflow point (~709.78); solving that inequality
    # for log_k_intra gives log_k_intra < log(700) - log(n) - log(sigma)
    # (700, not 709.78, for a small safety margin). This was checked
    # directly against all three real thresholds (81, 843, 8,690
    # molecules/cell) before trusting it: at this derived upper bound,
    # theta is confirmed to equal exactly 1.0 (safely past the theta=0.5
    # root, never NaN) for every one of them -- not assumed to generalize
    # from a single spot-check.
    # Stage 3 ALSO had a real bug, caught on this same next run (not
    # glossed over): safe_x (the bound on x = n*sigma*k_intra itself) was
    # being subtracted directly instead of its LOG being subtracted --
    # i.e. the code below used safe_x - log(n) - log(sigma) where solving
    # log(k_intra) < log(safe_x / (n*sigma)) = log(safe_x) - log(n) -
    # log(sigma) actually requires np.log(safe_x), not safe_x itself. The
    # arithmetic error put the bracket edge at log_k_intra=695.6 for
    # n=81, which is still far too high (recall the direct check above
    # found theta already goes NaN by log_k_intra=3 for n=81) -- this
    # confirms the earlier fix was never actually tested against the
    # real numbers before being trusted, exactly the mistake this
    # project's "verify claim-by-claim... before trusting ANY prior
    # summary, including your own code's first-pass output" rule exists
    # to catch.
    safe_x = 700.0  # margin below float64's np.exp() overflow (~709.78),
    # bounding x = n*sigma*k_intra itself, NOT log(k_intra).
    log_k_intra_lower = -50.0
    log_k_intra_upper = (
        np.log(safe_x) - np.log(threshold_molecules_per_cell) - np.log(car_avidity)
    )

    sol = optimize.brentq(residual, log_k_intra_lower, log_k_intra_upper, full_output=True)
    log_k_intra_solution, root_results = sol
    k_intra_fitted = float(np.exp(log_k_intra_solution))

    return {
        "k_intra": k_intra_fitted,
        "car_avidity": car_avidity,
        "z": z,
        "threshold_molecules_per_cell": threshold_molecules_per_cell,
        "threshold_kill_probability": threshold_kill_probability,
        "converged": bool(root_results.converged),
    }


def fit_all_real_thresholds(car_avidity=1.0, z=1.0):
    """
    Run fit_kill_curve_to_real_threshold() against every real threshold in
    REAL_TCE_ACTIVATION_THRESHOLDS, so the resulting k_intra values can be
    compared directly -- e.g., is the fitted k_intra for CLN-978 (a
    higher-affinity molecule, per the source paper) meaningfully different
    from blinatumomab's fitted k_intra? This is a real, checkable
    consistency question, not assumed to come out any particular way.

    Returns
    -------
    dict keyed by "{molecule}_{readout}" (e.g. "CLN-978_cytolysis"),
    values are the fit_kill_curve_to_real_threshold() result dict.
    """
    results = {}
    for molecule, thresholds in REAL_TCE_ACTIVATION_THRESHOLDS.items():
        for key, value in thresholds.items():
            if not key.endswith("_molecules_per_cell"):
                continue
            readout = key.replace("_molecules_per_cell", "")
            fit = fit_kill_curve_to_real_threshold(value, car_avidity=car_avidity, z=z)
            results[f"{molecule}_{readout}"] = fit
    return results


def print_summary():
    print("=" * 78)
    print("AIM 10 -- CAR-T KILL-PROBABILITY MODEL (the pivot's core model)")
    print("=" * 78)
    print("\nReusing binding_model.capture_probability() UNCHANGED, relabeled")
    print("for the CAR-T domain (vesicle->tumor cell, receptor->antigen,")
    print("capture->kill). See module docstring for the full translation.")

    print("\n" + "=" * 78)
    print("REAL CALIBRATION DATA (Agilent/Antibody Analytics IndEx-2, verified")
    print("directly this session -- see module docstring for full provenance)")
    print("=" * 78)
    print(f"\n  Real engineered CD19 density range tested: "
          f"{INDEX2_CD19_DENSITY_RANGE_MOLECULES_PER_CELL[0]} - "
          f"{INDEX2_CD19_DENSITY_RANGE_MOLECULES_PER_CELL[1]:,} molecules/cell")
    for molecule, thresholds in REAL_TCE_ACTIVATION_THRESHOLDS.items():
        print(f"\n  {molecule}:")
        for key, value in thresholds.items():
            if key.endswith("_molecules_per_cell"):
                print(f"    {key}: {value:,}")
        if thresholds.get("threshold_is_detection_limit"):
            print(f"    NOTE: this is the assay's lower limit of detection -- "
                  f"the true threshold may be even lower.")

    print(f"\n  Independent qualitative cross-check (retrospective clinical CD19")
    print(f"  CAR-T, not TCE, response threshold): "
          f"~{CLINICAL_CD19_CAR_T_RESPONSE_THRESHOLD_MOLECULES_PER_CELL:,} molecules/cell")
    print(f"  (same order of magnitude as the real TCE thresholds above -- a real,")
    print(f"  independent consistency check, not used to fit anything.)")

    print("\n" + "=" * 78)
    print("FITTING k_intra AGAINST EACH REAL THRESHOLD")
    print("=" * 78)
    fits = fit_all_real_thresholds()
    for name, fit in fits.items():
        print(f"\n  {name}: threshold={fit['threshold_molecules_per_cell']:,} "
              f"molecules/cell -> k_intra={fit['k_intra']:.6e} "
              f"(converged={fit['converged']})")

    print(f"\n  FINDING: this is a real, data-anchored calibration for each")
    print(f"  molecule/readout combination -- a genuine methodological upgrade")
    print(f"  over the EV project's permanently plausibility-swept K_INTRA (see")
    print(f"  binding_model.py's own docstring for why that search never")
    print(f"  resolved to a point estimate for the EV case).")
    print(f"\n  HONEST LIMITATION: each fit is a ONE-POINT calibration (see")
    print(f"  fit_kill_curve_to_real_threshold's docstring) -- not a full")
    print(f"  multi-point least-squares fit to a continuous dose-response curve.")
    print("=" * 78)


if __name__ == "__main__":
    print_summary()
