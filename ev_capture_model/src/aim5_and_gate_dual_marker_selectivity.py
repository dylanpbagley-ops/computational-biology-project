"""
aim5_and_gate_dual_marker_selectivity.py

WHY THIS FILE EXISTS
---------------------
This project's existing multi-marker work (aim3b_cooccurrence_model.py)
answers an OR-GATE question: if a capture surface catches a vesicle that
is positive for EITHER of two markers, how much more total YIELD do you
get versus your best single marker? That work found a real but modest
answer (roughly 1.02x-1.26x combination gain, well below the project's
original informal 2-4x expectation), fit directly to Mizenko et al.'s
tetraspanin/EpCAM/Her2 data.

This file asks a DIFFERENT, mechanistically distinct question that the
project had discussed (as far back as the v5 handoff's "multi-marker
capture is a core sub-aim") but never actually built: an AND-GATE
question. If a capture surface requires a vesicle to be positive for
BOTH of two markers simultaneously (rather than either one), does that
sharpen the TUMOR-VS-NORMAL SELECTIVITY RATIO -- not just change total
yield -- beyond what either single marker achieves alone?

This is a real, independently motivated strategy in the EV field (dual-
aptamer AND-logic-gate detection schemes and dual-surface-protein
orthogonal-barcoding approaches both exist in the published literature
for exactly this reason -- see the project's handoff for citations found
during the novelty-scoping conversation that led to this file), and it is
mechanistically different from the OR-gate question: an AND-gate can, in
principle, suppress off-target (normal) capture much more than it
suppresses on-target (tumor) capture, IF the tumor and normal populations
differ more on the JOINT two-marker distribution than on either marker's
marginal distribution alone.

THE HONEST DATA PROBLEM, STATED UP FRONT
-------------------------------------------
Lennon et al.'s real digitized data (ev_capture_data/) has real per-
vesicle EGFR counts (from one scatter plot) and real per-vesicle CA19-9
counts (from a DIFFERENT scatter plot) for the same two populations
(healthy plasma, PDAC plasma) -- but NOT both markers measured on the
SAME individual vesicle. This is the exact same structural gap Aim 3b
hit with Mizenko's data (see that file's own module docstring), and it is
handled the same honest way here: rather than either (a) naively
assuming independence and presenting one number as if it were measured,
or (b) refusing to build anything because the ideal dataset doesn't
exist, this file builds a JOINT reconstruction anchored to the two REAL
marginal distributions (real mean, real dispersion, fit directly from the
real digitized per-vesicle data for each marker/population), with an
EXPLICIT, SWEPT rank-correlation parameter between the two markers. The
independence case (correlation=0) is reported as the conservative
default; a plausible positive-correlation range is swept as a first-class
sensitivity axis, not hidden. This mirrors exactly the discipline
Limitations_and_Calibration_Gap.md and the Aim 3b write-up already use
for other unconstrained parameters in this project: report the number
that's supportable, and report the range for the number that isn't.

THE PHYSICAL MODEL, STATED EXPLICITLY
-----------------------------------------
An AND-gate capture surface is modeled as co-immobilizing BOTH antibody
types (anti-EGFR and anti-CA19-9) on the same surface. This file treats
binding through each marker type as an independent event (a standard,
defensible simplification -- not a claim that this is the only possible
physical model): a vesicle is captured only if EITHER marker's own
independent capture_probability_patch() draw succeeds for AND-logic, or
the corresponding OR-logic combination for comparison. This is
DELIBERATELY the simpler of two possible AND-gate physical models: a
true single multivalent bond-counting process across BOTH marker types
simultaneously (bonds of either type contributing to one shared
q_bound-style partition function) would be a larger physics change to
the core model and is named explicitly in "Immediate next steps" below
as a real, unbuilt extension -- not silently assumed away.

WHAT THIS FILE DOES NOT DO
------------------------------
- Does NOT claim a measured, real dual-marker selectivity number. Every
  reported ratio is either the independence case (an honest, clearly-
  labeled default) or a swept function of an assumed correlation.
- Does NOT touch or recalibrate K_INTRA, z, or k_geom -- same permanently
  uncalibrated placeholders as everywhere else in this project.
- Does NOT build the "true single shared multivalent bond-counting
  across two marker types" alternative physical model (bonds of EITHER
  type contributing to one shared q_bound-style partition function,
  rather than two independent capture_probability_patch() draws
  multiplied together) -- a real, unbuilt extension, left for a future
  session if this line of work is pursued further.
- Is numerically UNRELIABLE at very dilute z combined with high assumed
  correlation (found and diagnosed this session -- see
  AND_GATE_MAX_RELIABLE_TOP_1PCT_CONTRIBUTION's docstring for the full
  account). The z=1.0 default this file reports by default is the
  checked-stable regime; any result computed at a different z should be
  checked against its own "reliable" flag before being trusted or cited.
"""

import os
import sys
import csv

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from binding_model import capture_probability_patch
from k_intra_z_plausibility_bounds import (
    k_intra_bounds_from_real_sigma,
    SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
)
from aim3_aim4_range_sweep import get_k_intra_sweep_values

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "ev_capture_data")


# ---------------------------------------------------------------------------
# PART 1: load the two REAL marginal distributions (EGFR, CA19-9) for the
# healthy and PDAC-plasma populations, and fit each to a negative binomial
# (same distribution family used everywhere else in this project for
# reconstructed marginals -- see population_model.py).
# ---------------------------------------------------------------------------

def _load_pairs(filename):
    """Load (copies_per_ev, diameter_nm) pairs from a digitized CSV -- same
    loading convention as aim3_aim4_range_sweep.py's _load_pairs(), kept
    as a separate local copy here (rather than importing the other
    module's private helper) so this file's data dependencies are fully
    self-contained and explicit.
    """
    diam, copies = [], []
    path = os.path.join(DATA_DIR, filename)
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            d = float(row["diameter_nm"])
            c = float(row["copies_per_ev"])
            if d > 0:
                diam.append(d)
                copies.append(c)
    return np.array(copies), np.array(diam)


# Real digitized files for the two markers, both populations (Lennon et
# al. 2019, Fig 3c, healthy/PDAC plasma -- the same pair used as this
# project's primary CA19-9 headline finding and its retired EGFR case).
MARKER_FILES = {
    "EGFR": {"PDAC": "fig3c_EGFR_PDAC_digitized.csv", "healthy": "fig3c_EGFR_healthy_digitized.csv"},
    "CA199": {"PDAC": "fig3c_CA199_PDAC_digitized.csv", "healthy": "fig3c_CA199_healthy_digitized.csv"},
}


def fit_marginal(copies, dispersion_r=2.0):
    """
    Fit a negative binomial to a real digitized copies_per_ev array by
    matching its MEAN exactly (the real, measured mean of the actual
    digitized data) at a fixed dispersion_r. Unlike stem_cell_marker_
    ranking_extension.py's fit_negative_binomial_to_positivity() (which
    solves for dispersion_r to hit a target zero-fraction, because that
    file only has POSITIVITY FRACTION + mean as its two real summary
    statistics), this function's job is simpler: the real data here is a
    full per-vesicle sample, so the real empirical mean can be used
    directly as the negative binomial's target mean, with dispersion_r
    fixed at a standard value (2.0, matching the shared default already
    used in aim3b_cooccurrence_model.simulate_joint_population() and
    population_model.negative_binomial_distribution()'s own module
    discussion of "moderate overdispersion") rather than fit from the
    data. This keeps this file's marginals directly traceable to the
    real digitized means, while making the (necessary, since a copula
    construction needs an explicit parametric marginal) distribution-
    family choice explicit and consistent with the rest of the project.

    Parameters
    ----------
    copies : array
        Real digitized copies_per_ev values for one marker/population.
    dispersion_r : float
        Fixed negative-binomial dispersion parameter.

    Returns
    -------
    dict with "mean_n" (the real empirical mean, used directly),
    "dispersion_r".
    """
    mean_n = float(np.mean(copies))
    return {"mean_n": mean_n, "dispersion_r": dispersion_r}


def fit_all_marginals(dispersion_r=2.0):
    """
    Fit every (marker, population) marginal in MARKER_FILES, returning a
    dict keyed by (marker, population) -> fit dict (see fit_marginal).
    Also returns the real diameter array for each population (used
    directly, per-vesicle, matching the convention everywhere else in
    this project of using real digitized diameters rather than an
    assumed diameter distribution wherever real diameter data exists).
    """
    fits = {}
    diameters = {}
    for marker, pop_files in MARKER_FILES.items():
        for pop, filename in pop_files.items():
            copies, diam = _load_pairs(filename)
            fits[(marker, pop)] = fit_marginal(copies, dispersion_r=dispersion_r)
            # Diameter is per-population, not per-marker (a given
            # population of vesicles has one size distribution regardless
            # of which marker's plot digitized it) -- but Lennon's two
            # marker plots for the SAME population were digitized from
            # DIFFERENT scatter panels, so their diameter arrays are not
            # literally the same array. Keyed by population only, using
            # whichever marker's diameters were loaded LAST for that
            # population is arbitrary, so instead each (marker, population)
            # diameter array is kept separately and the correlation-
            # construction step below draws its OWN diameter directly
            # from the same real array it draws marker counts from,
            # never mixing a diameter from one marker's plot with a
            # count from the other's.
            diameters[(marker, pop)] = diam
    return fits, diameters


# ---------------------------------------------------------------------------
# PART 2: build a joint (EGFR, CA19-9) reconstruction via a Gaussian
# copula, at an explicit, swept rank correlation, anchored to the two
# real marginals above.
# ---------------------------------------------------------------------------

def simulate_joint_two_marker_population(fits, population, correlation,
                                           n_vesicles=50_000, seed=0):
    """
    Simulate a population of vesicles with jointly-drawn (EGFR, CA19-9)
    counts, at a given Gaussian-copula rank correlation, anchored to the
    REAL marginal fits for the given population (healthy or PDAC).

    HOW THE COPULA WORKS, IN PLAIN TERMS: draw two correlated standard
    normal variables (correlation = the requested `correlation`
    parameter, exactly), then convert each one to a negative-binomial
    count by matching its percentile under the standard normal CDF to
    the same percentile under the target marker's negative-binomial CDF
    (the standard "Gaussian copula" construction). This produces two
    count variables whose INDIVIDUAL marginal distributions are exactly
    the two real fitted negative binomials (by construction, checked
    directly in tests), while their JOINT dependence structure is
    controlled by the single `correlation` parameter -- correlation=0
    reproduces independence exactly; correlation>0 makes vesicles that
    are high in one marker more likely to also be high in the other,
    matching the qualitative direction Aim 3b's real, robust finding
    (EpCAM/Her2 are positively correlated, not independent) already
    established for a different marker pair on real data.

    Parameters
    ----------
    fits : dict
        Output of fit_all_marginals() -- keyed by (marker, population).
    population : str
        "PDAC" or "healthy".
    correlation : float
        Target Gaussian-copula correlation, in [-1, 1]. 0 = independence.
    n_vesicles : int
        Simulated population size.
    seed : int
        RNG seed.

    Returns
    -------
    dict with "egfr" (int array), "ca199" (int array), both length
    n_vesicles, jointly drawn as described above.
    """
    if not (-1.0 <= correlation <= 1.0):
        raise ValueError(f"correlation must be in [-1, 1], got {correlation}")

    rng = np.random.default_rng(seed)
    mean_cov = [0.0, 0.0]
    cov = [[1.0, correlation], [correlation, 1.0]]
    z1, z2 = rng.multivariate_normal(mean_cov, cov, size=n_vesicles).T

    u1 = stats.norm.cdf(z1)
    u2 = stats.norm.cdf(z2)

    egfr_fit = fits[("EGFR", population)]
    ca199_fit = fits[("CA199", population)]

    egfr_p = egfr_fit["dispersion_r"] / (egfr_fit["dispersion_r"] + egfr_fit["mean_n"])
    ca199_p = ca199_fit["dispersion_r"] / (ca199_fit["dispersion_r"] + ca199_fit["mean_n"])

    egfr = stats.nbinom.ppf(u1, egfr_fit["dispersion_r"], egfr_p)
    ca199 = stats.nbinom.ppf(u2, ca199_fit["dispersion_r"], ca199_p)

    return {"egfr": egfr.astype(float), "ca199": ca199.astype(float)}


# ---------------------------------------------------------------------------
# PART 3: the AND-gate (and OR-gate, for comparison) capture model, and
# the tumor-vs-normal selectivity-ratio sweep across correlation.
# ---------------------------------------------------------------------------

def and_gate_capture_probability(n_egfr, n_ca199, diameter, sigma, k_geom, k_intra, z,
                                   patch_exponent=2.0):
    """
    Per-vesicle AND-gate capture probability: the vesicle is captured
    only if it would ALSO be captured through EGFR-only binding AND
    ALSO through CA19-9-only binding, treating the two marker-antibody
    interactions as INDEPENDENT binding events on a surface that
    co-immobilizes both antibody types (see module docstring's "THE
    PHYSICAL MODEL" section for the explicit statement of this
    simplification and its limits).

    theta_AND = theta_EGFR(n_egfr) * theta_CA199(n_ca199)

    Both single-marker probabilities use the SAME diameter, sigma,
    k_geom, k_intra, z, and patch_exponent -- i.e., this does not model
    differential antibody affinity or surface density between the two
    marker types, another explicit simplification (a real capture
    surface's two antibody spots could well have different k_intra/sigma
    values; that is a further, unbuilt sensitivity axis -- see
    "Immediate next steps" below).

    Returns
    -------
    float or array, same shape as inputs: the joint AND-gate capture
    probability, in [0, 1] by construction (a product of two [0,1]
    probabilities).
    """
    theta_egfr = capture_probability_patch(n_egfr, diameter, sigma, k_geom, k_intra, z, patch_exponent)
    theta_ca199 = capture_probability_patch(n_ca199, diameter, sigma, k_geom, k_intra, z, patch_exponent)
    return theta_egfr * theta_ca199


def or_gate_capture_probability(n_egfr, n_ca199, diameter, sigma, k_geom, k_intra, z,
                                  patch_exponent=2.0):
    """
    Per-vesicle OR-gate capture probability, for direct comparison against
    the AND-gate above: theta_OR = 1 - (1-theta_EGFR)*(1-theta_CA199),
    i.e. captured if EITHER marker's independent binding event succeeds.
    This is the multi-marker logic conceptually closest to Aim 3b's
    existing OR-gate yield-gain question, but computed here through the
    SAME per-vesicle multivalent binding model (capture_probability_patch)
    used throughout Aims 1-2/4, rather than Aim 3b's simpler
    "positive/not positive" threshold logic -- included here as a
    reference point, not as this file's main question.
    """
    theta_egfr = capture_probability_patch(n_egfr, diameter, sigma, k_geom, k_intra, z, patch_exponent)
    theta_ca199 = capture_probability_patch(n_ca199, diameter, sigma, k_geom, k_intra, z, patch_exponent)
    return 1.0 - (1.0 - theta_egfr) * (1.0 - theta_ca199)


def population_and_gate_yield(pop, diameter_egfr, diameter_ca199, sigma, k_geom, k_intra, z,
                                patch_exponent=2.0):
    """
    Population-average AND-gate yield for a simulated joint population
    (output of simulate_joint_two_marker_population()). Diameter is drawn
    per-vesicle from the population's REAL digitized diameter array (see
    _draw_matched_diameters() below for how a diameter is matched to each
    synthetic vesicle) -- separately for the "as if measured via the EGFR
    plot" and "as if measured via the CA19-9 plot" diameter arrays, since
    (per the module docstring) these are not literally the same array in
    Lennon's real data. The mean of the two per-vesicle patch-adjusted
    contact terms is used, keeping the model anchored to real diameters
    from both source panels rather than silently picking one arbitrarily.

    Returns
    -------
    dict with "yield" (the population-mean AND-gate capture probability)
    and "effective_contributing_fraction" (see
    and_gate_reliability_diagnostic() below for what this measures and
    why it matters -- checked directly this session after a real
    numerical instability was found and diagnosed at very dilute z).
    """
    theta = and_gate_capture_probability(
        pop["egfr"], pop["ca199"], (diameter_egfr + diameter_ca199) / 2.0,
        sigma, k_geom, k_intra, z, patch_exponent
    )
    mean_theta = float(np.mean(theta))
    # See and_gate_reliability_diagnostic() docstring for why this
    # specific diagnostic (fraction of the mean contributed by the
    # top 1% of vesicles) is the right numerical-reliability check here.
    if mean_theta > 0:
        sorted_theta = np.sort(theta)[::-1]
        top_1pct_n = max(1, len(theta) // 100)
        top_1pct_contribution = sorted_theta[:top_1pct_n].sum() / (mean_theta * len(theta))
    else:
        top_1pct_contribution = 0.0
    return {"yield": mean_theta, "top_1pct_contribution": float(top_1pct_contribution)}


def _draw_matched_diameters(diameters, marker, population, n_vesicles, seed):
    """
    Draw n_vesicles diameters by resampling WITH REPLACEMENT from the
    real digitized diameter array for the given (marker, population) --
    used to give each synthetic AND-gate vesicle a real, empirically-
    grounded diameter rather than an assumed distribution, while still
    allowing the synthetic sample size (n_vesicles) to differ from the
    real digitized N (typically much larger, for smoother sweep
    statistics, matching the convention already used in
    stem_cell_marker_ranking_extension.py's reconstruct_marker_population,
    which also draws a larger synthetic n_samples than any single real
    digitized panel).
    """
    rng = np.random.default_rng(seed + 500_000)  # offset so diameter
                                                    # resampling doesn't
                                                    # share a stream with
                                                    # the copula draws
    real_diam = diameters[(marker, population)]
    return rng.choice(real_diam, size=n_vesicles, replace=True)


# A cell (correlation, k_intra, z combination) is flagged numerically
# unreliable if either population's AND-gate yield has more than this
# fraction of its total contribution coming from the top 1% of vesicles
# in the simulated sample -- i.e., the population mean is being driven
# by a small handful of extreme-tail synthetic vesicles rather than a
# broad, stable majority of the sample. FOUND AND DIAGNOSED THIS
# SESSION: at very dilute z (e.g. z=1e-4, the same dilute regime this
# project's OWN superselectivity-reachability finding lives in) combined
# with high assumed EGFR/CA19-9 correlation, the AND-gate ratio was
# checked directly across four different random seeds and swung by as
# much as 50% (67 to 99 at correlation=0.9, z=1e-4, n_vesicles=50,000) --
# a real, reproducible numerical instability, not a transient fluke. This
# is the SAME general failure mode this project's MIN_RELIABLE_EXPONENT_
# AT_N1 guard (in aim3_aim4_range_sweep.py) protects against for a
# different reason (there, floating-point precision at the low end of an
# exponential; here, Monte Carlo sampling noise in a genuinely
# rare-event tail) -- both are "don't trust a mean that's secretly being
# computed from a handful of extreme points" checks, and this project's
# standing practice is to make that check explicit and reported, not
# silently absorbed into a single trusted number. A threshold of 0.5 (the
# top 1% of vesicles contributing more than half the population mean) was
# chosen as a conservative, round trigger point -- checked directly
# against the diagnosed instability above (which showed top-1%
# contributions well past this threshold) and against the STABLE z=1.0
# default case (which showed top-1% contributions far below it).
AND_GATE_MAX_RELIABLE_TOP_1PCT_CONTRIBUTION = 0.5


def and_gate_selectivity_sweep(correlations=(0.0, 0.2, 0.4, 0.6, 0.8),
                                 k_intra=None, z=1.0,
                                 sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
                                 n_vesicles=50_000, seed=0):
    """
    THE MAIN QUESTION THIS FILE ANSWERS: across a range of assumed
    EGFR/CA19-9 correlation strengths, how does the AND-gate tumor-vs-
    normal selectivity ratio (PDAC AND-gate yield / healthy AND-gate
    yield) compare to the best SINGLE-marker selectivity ratio (the
    better of the EGFR-only or CA19-9-only ratio) at the same k_intra/z?

    Uses the SAME real-sigma-anchored k_intra default as every other
    range-sweep function in this project (the largest value in
    k_intra_bounds_from_real_sigma()'s range) and the SAME z=1.0 default
    aim3_aim4_range_sweep.sweep_ranking_flip_robustness() uses, unless
    overridden. IMPORTANT, FOUND THIS SESSION: unlike that function, this
    one's result becomes numerically UNRELIABLE at very dilute z (e.g.
    z=1e-4) combined with high correlation -- see
    AND_GATE_MAX_RELIABLE_TOP_1PCT_CONTRIBUTION above for the full
    diagnosis. Every returned row is checked against this reliability
    criterion and flagged, not silently trusted -- see "reliable" in the
    Returns section below. The z=1.0 default used here is the STABLE
    regime (checked directly: top-1% contribution stays under 5% at
    every tested correlation at z=1.0), so the headline numbers this
    file reports by default are numerically trustworthy; the dilute-
    regime instability is documented as a known risk for anyone who
    reruns this sweep at a different z, not swept under the rug.

    Parameters
    ----------
    correlations : tuple of float
        Gaussian-copula correlations to sweep, applied IDENTICALLY to
        both the PDAC and healthy populations (a simplifying assumption,
        stated explicitly: this does not test the possibility that
        correlation strength itself differs between tumor and normal
        vesicles, which is a further, unbuilt sensitivity axis -- see
        "Immediate next steps" below).
    k_intra, z, sigma : float
        Binding-model parameters, same meaning as everywhere else in this
        project. k_intra defaults to the largest plausibility-range value
        if None.
    n_vesicles, seed : int
        Passed through to the joint-population simulation.

    Returns
    -------
    dict with:
        "correlations": the swept correlation values (echoed back).
        "and_gate_ratio": list of (PDAC AND-gate yield / healthy AND-gate
            yield), one per correlation.
        "egfr_only_ratio", "ca199_only_ratio": the corresponding single-
            marker selectivity ratios (constant across correlation, since
            a single marker's own capture probability does not depend on
            the OTHER marker's correlation with it -- computed once,
            included at every row purely for side-by-side comparison).
        "best_single_marker_ratio": max(egfr_only_ratio, ca199_only_ratio).
        "and_gate_beats_best_single_marker": list of bool, one per
            correlation -- whether the AND-gate ratio exceeds the best
            achievable single-marker ratio at that correlation.
        "reliable": list of bool, one per correlation -- NEW, added after
            this session's diagnosed numerical instability. False means
            EITHER population's AND-gate yield at that correlation had
            more than AND_GATE_MAX_RELIABLE_TOP_1PCT_CONTRIBUTION of its
            mean driven by the top 1% of simulated vesicles -- treat that
            row's and_gate_ratio as numerically untrustworthy, not as a
            real finding, regardless of how large or small it looks.
    """
    if k_intra is None:
        k_intra = max(k_intra_bounds_from_real_sigma().values())

    fits, diameters = fit_all_marginals()

    # k_geom convention, matching aim3_aim4_range_sweep.sweep_ranking_
    # flip_robustness() EXACTLY (this is a correctness-critical detail,
    # not a style choice -- see this file's own regression test,
    # test_single_marker_ratios_match_aim4_range_sweep, which was added
    # specifically because an earlier version of this function computed
    # k_geom SEPARATELY for each population from its own mean, which
    # silently introduces k_geom as a hidden second free parameter that
    # partially cancels the real diameter/copy-number difference between
    # populations -- caught by comparing against aim3_aim4_range_sweep's
    # already-established, already-tested CA19-9 ratio before trusting
    # this file's own number). The correct convention: k_geom is anchored
    # ONCE, to the NORMAL (healthy) population's mean, and the SAME
    # k_geom value is then used for BOTH populations' yield.
    def single_marker_yield(marker, population, k_geom):
        copies, diam = _load_pairs(MARKER_FILES[marker][population])
        theta = capture_probability_patch(copies, diam, sigma, k_geom, k_intra, z, patch_exponent=2.0)
        return float(np.mean(theta))

    egfr_healthy_copies, _ = _load_pairs(MARKER_FILES["EGFR"]["healthy"])
    ca199_healthy_copies, _ = _load_pairs(MARKER_FILES["CA199"]["healthy"])
    egfr_k_geom = 1.0 / egfr_healthy_copies.mean()
    ca199_k_geom = 1.0 / ca199_healthy_copies.mean()

    egfr_pdac_yield = single_marker_yield("EGFR", "PDAC", egfr_k_geom)
    egfr_healthy_yield = single_marker_yield("EGFR", "healthy", egfr_k_geom)
    ca199_pdac_yield = single_marker_yield("CA199", "PDAC", ca199_k_geom)
    ca199_healthy_yield = single_marker_yield("CA199", "healthy", ca199_k_geom)

    egfr_only_ratio = egfr_pdac_yield / egfr_healthy_yield if egfr_healthy_yield > 0 else float("inf")
    ca199_only_ratio = ca199_pdac_yield / ca199_healthy_yield if ca199_healthy_yield > 0 else float("inf")
    best_single = max(egfr_only_ratio, ca199_only_ratio)

    # AND-gate k_geom: a single shared value is needed since the AND-gate
    # probability multiplies two capture_probability_patch() calls that
    # must use consistent geometry, and (per the fix above) each single-
    # marker k_geom is itself anchored to the healthy population. The
    # mean of the two (correctly healthy-anchored) single-marker k_geom
    # values is used -- a new modeling choice this file introduces (no
    # single-marker sweep in this project has needed to reconcile two
    # different k_geom values before), flagged explicitly rather than
    # silently picked.
    and_gate_k_geom = 0.5 * (egfr_k_geom + ca199_k_geom)

    and_gate_ratios = []
    reliable_flags = []
    for correlation in correlations:
        pdac_pop = simulate_joint_two_marker_population(
            fits, "PDAC", correlation, n_vesicles=n_vesicles, seed=seed
        )
        healthy_pop = simulate_joint_two_marker_population(
            fits, "healthy", correlation, n_vesicles=n_vesicles, seed=seed + 1
        )
        pdac_diam_egfr = _draw_matched_diameters(diameters, "EGFR", "PDAC", n_vesicles, seed)
        pdac_diam_ca199 = _draw_matched_diameters(diameters, "CA199", "PDAC", n_vesicles, seed)
        healthy_diam_egfr = _draw_matched_diameters(diameters, "EGFR", "healthy", n_vesicles, seed + 1)
        healthy_diam_ca199 = _draw_matched_diameters(diameters, "CA199", "healthy", n_vesicles, seed + 1)

        pdac_result = population_and_gate_yield(
            pdac_pop, pdac_diam_egfr, pdac_diam_ca199, sigma, and_gate_k_geom, k_intra, z
        )
        healthy_result = population_and_gate_yield(
            healthy_pop, healthy_diam_egfr, healthy_diam_ca199, sigma, and_gate_k_geom, k_intra, z
        )
        pdac_yield = pdac_result["yield"]
        healthy_yield = healthy_result["yield"]
        ratio = pdac_yield / healthy_yield if healthy_yield > 0 else float("inf")
        and_gate_ratios.append(ratio)

        is_reliable = (
            pdac_result["top_1pct_contribution"] < AND_GATE_MAX_RELIABLE_TOP_1PCT_CONTRIBUTION
            and healthy_result["top_1pct_contribution"] < AND_GATE_MAX_RELIABLE_TOP_1PCT_CONTRIBUTION
        )
        reliable_flags.append(is_reliable)

    return {
        "correlations": list(correlations),
        "and_gate_ratio": and_gate_ratios,
        "egfr_only_ratio": egfr_only_ratio,
        "ca199_only_ratio": ca199_only_ratio,
        "best_single_marker_ratio": best_single,
        "and_gate_beats_best_single_marker": [r > best_single for r in and_gate_ratios],
        "reliable": reliable_flags,
    }


def print_summary():
    print("=" * 78)
    print("AIM 5 -- AND-GATE DUAL-MARKER SELECTIVITY (EGFR + CA19-9)")
    print("=" * 78)
    result = and_gate_selectivity_sweep()
    print(f"\nBest single-marker selectivity ratio (PDAC/healthy): {result['best_single_marker_ratio']:.4f}")
    print(f"  EGFR-only ratio: {result['egfr_only_ratio']:.4f}")
    print(f"  CA19-9-only ratio: {result['ca199_only_ratio']:.4f}")
    print("\nAND-gate ratio across swept EGFR/CA19-9 correlation (z=1.0, the numerically stable regime):")
    for corr, ratio, beats, reliable in zip(result["correlations"], result["and_gate_ratio"],
                                              result["and_gate_beats_best_single_marker"],
                                              result["reliable"]):
        marker = "BEATS best single marker" if beats else "does not beat best single marker"
        flag = "" if reliable else "  [UNRELIABLE -- see AND_GATE_MAX_RELIABLE_TOP_1PCT_CONTRIBUTION]"
        print(f"  correlation={corr:+.1f}: AND-gate ratio={ratio:.4f}  ({marker}){flag}")

    print("\n--- Sensitivity check: does the AND-gate advantage hold across the full k_intra range? ---")
    for k_intra in sorted(get_k_intra_sweep_values()):
        r = and_gate_selectivity_sweep(k_intra=k_intra, correlations=(0.0, 0.8))
        print(f"  k_intra={k_intra:.3e}: best_single={r['best_single_marker_ratio']:.4f}, "
              f"AND@corr0={r['and_gate_ratio'][0]:.4f}, AND@corr0.8={r['and_gate_ratio'][1]:.4f}, "
              f"reliable={r['reliable']}")

    print("\n" + "=" * 78)


if __name__ == "__main__":
    print_summary()
