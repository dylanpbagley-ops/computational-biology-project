"""
aim3b_cooccurrence_model.py

WHY THIS FILE EXISTS
---------------------
"Immediate next steps" item 2 (formerly item 3, unchanged and unattempted
from v6 through v15): "Build Aim 3b's real co-occurrence model, using
src/mizenko_cooccurrence_scoping.py's scoping result as a starting target."

mizenko_cooccurrence_scoping.py answered a narrower, prior question
("is a 5-10x co-occurrence deficit even PHYSICALLY PLAUSIBLE from a simple
two-subpopulation mixture model, and roughly how strong does the
anti-correlation need to be") with a hand-picked table of six illustrative
(frac_tet_rich, epcam_mean_poor, epcam_mean_rich) triples, none of which
was actually FIT to Mizenko's real reported numbers -- the scoping file's
own docstring says plainly "this is a plausibility sweep... not a fitted
or validated result."

THIS file is the real model the scoping work was always pointed toward:
it actually FITS the same two-subpopulation mixture structure to
Mizenko et al. 2021's real reported ratios (~10x for EpCAM, ~5x for Her2),
using both markers TOGETHER (not one at a time, since both are measured on
the SAME tetraspanin-defined subpopulation structure in the paper), and
then uses the fitted model to answer the actual Aim 3b questions the
project's "Aim 3" section poses: how much does combining multiple markers
improve yield over any single marker, and does the naive INDEPENDENCE
assumption (multiply single-marker positive-fractions together) hold, or
does the real co-occurrence structure break it?

WHAT THIS FILE DOES NOT DO -- READ THIS BEFORE USING ANY NUMBER HERE FOR
A K_INTRA/z/k_geom CLAIM
--------------------------------------------------------------------------
This file NEVER touches binding_model.py's K_INTRA, z, or k_geom. Per the
handoff's "What NOT to revisit": Mizenko's 10x/5x numbers are a CROSS-
MARKER CO-OCCURRENCE measurement (how does EpCAM/Her2 positivity differ
between a tetraspanin-selected subpopulation and the whole population),
not a single-marker capture-efficiency measurement, and there is no
principled way to map a co-occurrence ratio onto single-bond binding-
strength parameters. This file's mixture-model parameters (the tet-rich
fraction, and each marker's suppression factor within that subpopulation)
describe POPULATION STRUCTURE -- how markers cluster together across
vesicles -- which is a genuinely different, and genuinely necessary,
model component from binding_model.py's per-vesicle biophysics. Both are
needed for a complete Aim 3 multi-marker analysis: this file supplies the
population structure P(tetraspanin, EpCAM, Her2) jointly; binding_model.py
(unchanged, reused as-is) supplies how any GIVEN vesicle's marker counts
translate into a capture probability once a joint distribution exists.

THE MODEL, BUILT UP EXPLICITLY
-----------------------------------
A population of vesicles is modeled as a two-component mixture:
  - a "tetraspanin-rich" subpopulation, fraction f_rich, which is what
    single-tetraspanin capture (anti-CD9/CD63/CD81) actually captures
    (idealized as: ALL of this subpopulation is captured, NONE of the
    tet-poor subpopulation is -- an idealized hard threshold, the same
    simplification mizenko_cooccurrence_scoping.py already used and
    flagged explicitly as not claiming real capture is this sharp).
  - a "tetraspanin-poor" subpopulation, fraction (1 - f_rich), which
    single-tetraspanin capture misses but non-specific (biotin) capture
    catches along with everything else.

Each subpopulation has its OWN mean copy number for EpCAM and for Her2
(independently settable per marker, since Mizenko's ratios differ between
the two markers -- 10x vs 5x -- so a single shared suppression factor
across both markers would not fit both numbers at once). Within each
subpopulation, marker counts are drawn from a negative binomial
distribution (matching population_model.py's own "zero-inflated,
overdispersed" convention elsewhere in this project, for consistency, and
because a real single copy-number-per-subpopulation is itself an
idealization already -- adding more distributional structure than that
would not be justified by what this ratio-only data can actually
constrain).

WHAT "FITTING" MEANS HERE, PRECISELY
-----------------------------------------
Given a candidate f_rich, and a candidate suppression ratio for one
marker (mean copies in tet-poor / mean copies in tet-rich -- ALWAYS >1 by
construction, matching Mizenko's finding that non-specific capture finds
MORE EpCAM+/Her2+ vesicles than single-tetraspanin capture), simulate the
resulting P(marker positive | tetraspanin capture) and P(marker positive |
non-specific capture), and compute their ratio -- exactly
epcam_positive_fraction_ratio() from the scoping file, reused directly
(not reimplemented) here. "Fitting" means finding, for f_rich fixed, the
suppression ratio that makes this simulated ratio match Mizenko's real
reported ratio (10.0 for EpCAM, 5.0 for Her2) as closely as possible --
a one-dimensional root-find (scipy.optimize.brentq), not a multi-parameter
optimization, because f_rich is shared across both markers (it describes
the SAME tetraspanin-defined subpopulation structure) while each marker's
own suppression ratio is fit independently GIVEN that f_rich. This
directly answers "Aim 3b's co-occurrence model is scoped, not built" from
the handoff's own "Current state: done vs. not done" section.

A REAL SENSITIVITY FINDING, CAUGHT WHILE FIRST RUNNING THIS FILE END TO
END -- READ THIS BEFORE CITING ANY "COMBINATION GAIN" NUMBER
--------------------------------------------------------------------------
MARKER_MEAN_TET_POOR (the absolute copy-number scale, held fixed since
Mizenko's RATIO-only data cannot constrain it -- see "WHY f_rich ITSELF IS
NOT ALSO FIT" below for the same reasoning applied to f_rich) turns out to
matter A LOT for one of this file's two headline outputs, but NOT AT ALL
for the other -- checked directly by sweeping MARKER_MEAN_TET_POOR over
2.0-30.0, not assumed:

  - independence_test()'s fold-deviation-from-independence result is
    STABLE across this whole range (~1.18x-1.18x, i.e. essentially
    unchanged) -- this is a robust, reportable structural finding.
  - multi_marker_yield_gain()'s combination-gain result is HIGHLY
    SENSITIVE to this same choice (1.26x at mean_tet_poor=2.0, versus only
    1.02x at mean_tet_poor=30.0) -- because a larger mean_tet_poor pushes
    the tet-poor subpopulation's own single-marker positivity fraction
    toward saturation (at mean=15, dispersion=2, P(positive)=98.6% --
    checked directly), leaving little room for ANY combination strategy to
    add anything, regardless of the real co-occurrence structure.

This is reported honestly rather than picking whichever MARKER_MEAN_TET_POOR
makes the combination-gain number look most like the handoff's own
existing "expect 2-4x, not 20x" expectation. See
combination_gain_sensitivity_sweep() below, which makes this sensitivity a
first-class, directly-callable result rather than something only found by
manually varying a module constant -- and see the "Immediate next steps"
entry in the handoff for what real data (an independent measurement of
the true absolute EpCAM/Her2 copy-number scale, not just Mizenko's
capture-method ratio) would be needed to pin this down further.

WHY f_rich ITSELF IS NOT ALSO FIT, AND WHAT THIS MEANS
------------------------------------------------------------
There is no independent, real measurement of f_rich (the tetraspanin-rich
fraction of the true EV population) in Mizenko et al. 2021 -- the paper
reports the AFTER-CAPTURE ratios, not the population's own true tetraspanin
positivity fraction as a fraction of an unbiased reference sample. Fitting
f_rich itself would require a THIRD piece of information (e.g. an
independent measurement of what fraction of a REAL EV population is
tetraspanin-positive by some method other than the two being compared)
that this paper does not supply. This file therefore treats f_rich as a
SWEPT, not fit, parameter -- exactly analogous to how the rest of this
project treats K_INTRA/z/k_geom as swept plausibility ranges rather than
fit point values when no calibration source exists (see
k_intra_z_plausibility_bounds.py) -- and reports, for EACH f_rich in a
plausible range (10-30%, per the scoping file's own finding), what
suppression ratio for EACH marker is needed to hit its real target, plus
whether that required suppression ratio is itself physically plausible
(e.g., not absurdly large).

INDEPENDENCE TEST -- "Immediate next steps" item 5, built directly on
TOP of this fitted model
--------------------------------------------------------------------------
Once EpCAM and Her2 are BOTH fit within the SAME two-subpopulation
structure (same f_rich, each with its own fitted suppression ratio), the
model makes a further, testable prediction: because both markers are
suppressed in the SAME tet-rich subpopulation and enhanced in the SAME
tet-poor subpopulation, EpCAM-positive and Her2-positive status are
predicted to be POSITIVELY CORRELATED across vesicles (both markers are
more likely to co-occur in tet-poor vesicles) -- even though nothing in
the model explicitly links EpCAM and Her2 to each other directly (they are
conditionally independent GIVEN tetraspanin-rich/poor status, but not
unconditionally independent). This is exactly the kind of "does the
independence assumption hold" question the handoff's "Multi-marker capture
is a core sub-aim" section flags: "the 9x math assumes the proteins appear
independently across vesicles. If they cluster, the multiplication fails."
independence_test() below computes exactly this: the TRUE joint
EpCAM+/Her2+ double-positive fraction (from the fitted mixture) versus the
NAIVE independence-assumption prediction (P(EpCAM+) x P(Her2+) multiplied
together, the standard field practice) -- and reports how far apart they
are, i.e., exactly how wrong the "multiply the single-marker rates
together" shortcut is for THIS system.
"""

import numpy as np
from scipy import stats
from scipy.optimize import brentq

# Reuse the scoping file's simulation function directly rather than
# reimplementing it -- this file builds ON TOP of that scoping work, not
# alongside a duplicate of it.
from mizenko_cooccurrence_scoping import epcam_positive_fraction_ratio


# ---------------------------------------------------------------------------
# Mizenko et al. 2021's real reported target ratios (see this project's own
# module docstring in mizenko_cooccurrence_scoping.py and the handoff's
# "Key literature" section for the full methodology and direct quotes).
# ---------------------------------------------------------------------------
MIZENKO_TARGET_RATIOS = {
    "EpCAM": 10.0,
    "Her2": 5.0,
}

# Plausible range for the tetraspanin-rich fraction, per the scoping file's
# own finding ("roughly 10-30% of the population"). Swept, not fit -- see
# module docstring for why.
F_RICH_PLAUSIBLE_RANGE = (0.10, 0.15, 0.20, 0.25, 0.30)

# A generic, held-fixed tet-poor mean copy number for each marker -- the
# ABSOLUTE scale is not constrained by a ratio-only measurement (only the
# RATIO of tet-rich to tet-poor means is), so this is fixed at a single
# representative value (matching the scoping file's own convention, e.g.
# its epcam_mean_tet_poor=15 configs) and only the SUPPRESSED (tet-rich)
# mean is solved for, given a target ratio -- see fit_suppression_ratio().
MARKER_MEAN_TET_POOR = 15.0


def _simulated_ratio_for_suppression(f_rich, mean_tet_poor, suppression_ratio,
                                       n_vesicles=300_000, seed=0):
    """
    Thin wrapper around mizenko_cooccurrence_scoping.epcam_positive_fraction_ratio(),
    parameterized by a SUPPRESSION RATIO (mean_tet_poor / mean_tet_rich)
    instead of the two means directly -- more convenient for root-finding,
    since the quantity being searched over (how many-fold suppressed is
    the marker in the tet-rich subpopulation) is the ratio itself.

    Parameters
    ----------
    f_rich : float
        Tetraspanin-rich fraction.
    mean_tet_poor : float
        Mean marker copy number in the tet-poor subpopulation (held fixed
        -- see module docstring, MARKER_MEAN_TET_POOR).
    suppression_ratio : float
        mean_tet_poor / mean_tet_rich. Must be > 1 for MEANINGFUL
        suppression (mean_tet_rich < mean_tet_poor, matching Mizenko's
        finding that single-tetraspanin capture UNDER-counts these
        markers relative to non-specific capture).
    n_vesicles, seed :
        Passed through to epcam_positive_fraction_ratio().

    Returns
    -------
    float
        The simulated P(marker+ | non-specific) / P(marker+ | tetraspanin
        capture) ratio -- directly comparable to Mizenko's reported 10x/5x.
    """
    mean_tet_rich = mean_tet_poor / suppression_ratio
    result = epcam_positive_fraction_ratio(
        frac_tet_rich=f_rich,
        epcam_mean_tet_poor=mean_tet_poor,
        epcam_mean_tet_rich=mean_tet_rich,
        n_vesicles=n_vesicles,
        seed=seed,
    )
    return result["ratio"]


def fit_suppression_ratio(f_rich, target_ratio, mean_tet_poor=MARKER_MEAN_TET_POOR,
                            suppression_search_range=(1.01, 500.0), n_vesicles=300_000, seed=0):
    """
    Given a fixed tetraspanin-rich fraction f_rich and a target co-
    occurrence ratio (Mizenko's real reported 10.0 for EpCAM or 5.0 for
    Her2), find the suppression ratio (mean_tet_poor / mean_tet_rich) that
    makes the simulated ratio match the target, via a 1-D root-find
    (scipy.optimize.brentq on the monotonic function suppression_ratio ->
    simulated_ratio - target_ratio).

    WHY brentq (BISECTION-STYLE ROOT-FINDING), NOT A GENERIC OPTIMIZER
    ------------------------------------------------------------------------
    The simulated ratio is a MONOTONICALLY INCREASING function of the
    suppression ratio (stronger suppression in the tet-rich subpopulation
    -> a bigger gap between tetraspanin-capture and non-specific-capture
    positivity -> a bigger ratio) -- checked directly before relying on
    this (see test suite). A monotonic 1-D function with a known target is
    exactly what a bracketing root-finder like brentq is built for: it is
    guaranteed to converge (unlike a generic gradient-based optimizer,
    which could in principle get stuck or overshoot) as long as the target
    is bracketed by the search range's endpoints, which is checked
    explicitly below rather than assumed.

    Parameters
    ----------
    f_rich : float
        Tetraspanin-rich fraction (held fixed for this fit).
    target_ratio : float
        The real ratio to fit to (10.0 for EpCAM, 5.0 for Her2).
    mean_tet_poor : float
        Held-fixed tet-poor mean copy number (see module docstring).
    suppression_search_range : tuple of float
        (low, high) bounds to search the suppression ratio within.
    n_vesicles, seed :
        Passed through to the simulation.

    Returns
    -------
    float
        The fitted suppression ratio.

    Raises
    ------
    ValueError
        If the target ratio is not bracketed by the search range (i.e.,
        even the maximum tested suppression ratio can't reach the target
        at this f_rich) -- reported explicitly rather than silently
        returning a boundary value that doesn't actually fit.
    """
    def objective(suppression_ratio):
        return _simulated_ratio_for_suppression(
            f_rich, mean_tet_poor, suppression_ratio, n_vesicles=n_vesicles, seed=seed
        ) - target_ratio

    low, high = suppression_search_range
    f_low, f_high = objective(low), objective(high)
    if f_low > 0 or f_high < 0:
        raise ValueError(
            f"Target ratio {target_ratio} is not bracketed by the search "
            f"range {suppression_search_range} at f_rich={f_rich} "
            f"(objective at range endpoints: {f_low + target_ratio:.3f}, "
            f"{f_high + target_ratio:.3f}). This f_rich may be too large "
            f"or too small to ever reach this target ratio -- widen the "
            f"search range or check whether this f_rich is physically "
            f"reasonable before assuming a fit exists."
        )

    return brentq(objective, low, high, xtol=1e-3)


def fit_both_markers(f_rich, n_vesicles=300_000, seed=0):
    """
    Fit BOTH EpCAM and Her2 suppression ratios at a shared f_rich, using
    Mizenko's real reported targets (MIZENKO_TARGET_RATIOS). This is the
    core "real co-occurrence model" fit this file exists to provide --
    the scoping file never fit anything to the real numbers; this does.

    Parameters
    ----------
    f_rich : float
        Shared tetraspanin-rich fraction.
    n_vesicles, seed :
        Passed through to the fits.

    Returns
    -------
    dict
        {"f_rich": ..., "EpCAM": {"suppression_ratio": ..., "mean_tet_poor": ...,
         "mean_tet_rich": ..., "target_ratio": ..., "achieved_ratio": ...},
         "Her2": {...same keys...}}
    """
    result = {"f_rich": f_rich}
    for marker, target in MIZENKO_TARGET_RATIOS.items():
        suppression = fit_suppression_ratio(
            f_rich, target, mean_tet_poor=MARKER_MEAN_TET_POOR,
            n_vesicles=n_vesicles, seed=seed,
        )
        achieved = _simulated_ratio_for_suppression(
            f_rich, MARKER_MEAN_TET_POOR, suppression, n_vesicles=n_vesicles, seed=seed + 1
        )  # re-simulated with a DIFFERENT seed than the fit itself, so
           # "achieved" is an independent check, not just re-reporting the
           # exact simulation the root-finder converged on.
        result[marker] = {
            "suppression_ratio": suppression,
            "mean_tet_poor": MARKER_MEAN_TET_POOR,
            "mean_tet_rich": MARKER_MEAN_TET_POOR / suppression,
            "target_ratio": target,
            "achieved_ratio": achieved,
        }
    return result


def fit_across_f_rich_range(f_rich_values=F_RICH_PLAUSIBLE_RANGE, n_vesicles=300_000, seed=0):
    """
    Run fit_both_markers() across the plausible f_rich range (swept, not
    fit -- see module docstring), so the write-up can report how the
    required suppression ratio for each marker depends on the assumed
    tetraspanin-rich fraction, rather than committing to one arbitrary
    f_rich value.

    Returns
    -------
    list of dict
        One fit_both_markers() result per f_rich value.
    """
    return [fit_both_markers(f_rich, n_vesicles=n_vesicles, seed=seed) for f_rich in f_rich_values]


def simulate_joint_population(f_rich, epcam_mean_tet_poor, epcam_mean_tet_rich,
                                her2_mean_tet_poor, her2_mean_tet_rich,
                                dispersion=2.0, n_vesicles=300_000, seed=0):
    """
    Simulate a population of vesicles with BOTH EpCAM and Her2 counts
    drawn jointly, under the two-subpopulation mixture structure -- i.e.,
    for each vesicle, first decide tet-rich or tet-poor, then draw ITS
    EpCAM count AND its Her2 count from the corresponding subpopulation's
    distributions. This is what makes EpCAM and Her2 co-occurrence emerge
    naturally from the shared subpopulation assignment (conditionally
    independent given tet-status, unconditionally correlated) rather than
    being imposed directly as a marker-marker correlation parameter.

    Parameters
    ----------
    f_rich : float
        Tetraspanin-rich fraction.
    epcam_mean_tet_poor, epcam_mean_tet_rich : float
        EpCAM mean copy number in each subpopulation.
    her2_mean_tet_poor, her2_mean_tet_rich : float
        Her2 mean copy number in each subpopulation.
    dispersion : float
        Shared negative-binomial dispersion parameter for both markers
        (same convention as population_model.negative_binomial_distribution
        -- smaller means more overdispersed/zero-inflated).
    n_vesicles : int
        Simulated population size.
    seed : int
        RNG seed.

    Returns
    -------
    dict with keys:
        "is_tet_rich" : boolean array, which subpopulation each vesicle is in
        "epcam" : int array, EpCAM copy number per vesicle
        "her2" : int array, Her2 copy number per vesicle
    """
    rng = np.random.default_rng(seed)
    is_tet_rich = rng.random(n_vesicles) < f_rich

    def draw(mean_poor, mean_rich):
        r = dispersion
        poor = stats.nbinom.rvs(r, r / (r + mean_poor), size=n_vesicles, random_state=rng)
        rich = stats.nbinom.rvs(r, r / (r + max(mean_rich, 1e-6)), size=n_vesicles, random_state=rng)
        return np.where(is_tet_rich, rich, poor)

    epcam = draw(epcam_mean_tet_poor, epcam_mean_tet_rich)
    her2 = draw(her2_mean_tet_poor, her2_mean_tet_rich)

    return {"is_tet_rich": is_tet_rich, "epcam": epcam, "her2": her2}


def independence_test(fit_result, n_vesicles=300_000, seed=2):
    """
    "Immediate next steps" item 5 (Aim 3b multi-marker joint distributions
    and the independence test), built directly on top of this file's
    fitted two-marker mixture model.

    Simulates the FULL joint (EpCAM, Her2) population implied by a
    fit_both_markers() result, then compares:
      - TRUE joint double-positive fraction: P(EpCAM+ AND Her2+), measured
        directly from the simulated joint population.
      - NAIVE INDEPENDENCE prediction: P(EpCAM+) x P(Her2+), the standard
        field-practice shortcut the handoff explicitly flags as untested
        ("the 9x math assumes the proteins appear independently across
        vesicles").

    If the two are close, the independence assumption is a reasonable
    approximation FOR THIS SYSTEM. If they differ substantially, that is
    itself the finding: real marker co-occurrence structure breaks the
    naive multiplication rule, and by how much.

    Parameters
    ----------
    fit_result : dict
        Output of fit_both_markers() -- must have "f_rich", "EpCAM", "Her2"
        keys with "mean_tet_poor"/"mean_tet_rich" sub-keys.
    n_vesicles, seed :
        Passed through to simulate_joint_population().

    Returns
    -------
    dict with keys:
        "p_epcam_positive", "p_her2_positive" : marginal positive fractions
        "p_joint_true" : true P(EpCAM+ AND Her2+) from the simulated joint population
        "p_joint_independence_prediction" : P(EpCAM+) * P(Her2+)
        "fold_deviation_from_independence" : p_joint_true / p_joint_independence_prediction
            (>1 means POSITIVE correlation -- co-occurrence is MORE common
            than independence predicts; <1 would mean anti-correlation)
    """
    epcam_info = fit_result["EpCAM"]
    her2_info = fit_result["Her2"]

    pop = simulate_joint_population(
        f_rich=fit_result["f_rich"],
        epcam_mean_tet_poor=epcam_info["mean_tet_poor"],
        epcam_mean_tet_rich=epcam_info["mean_tet_rich"],
        her2_mean_tet_poor=her2_info["mean_tet_poor"],
        her2_mean_tet_rich=her2_info["mean_tet_rich"],
        n_vesicles=n_vesicles,
        seed=seed,
    )

    epcam_pos = pop["epcam"] > 0
    her2_pos = pop["her2"] > 0

    p_epcam = epcam_pos.mean()
    p_her2 = her2_pos.mean()
    p_joint_true = (epcam_pos & her2_pos).mean()
    p_joint_independence = p_epcam * p_her2

    fold_deviation = p_joint_true / p_joint_independence if p_joint_independence > 0 else float("inf")

    return {
        "p_epcam_positive": float(p_epcam),
        "p_her2_positive": float(p_her2),
        "p_joint_true": float(p_joint_true),
        "p_joint_independence_prediction": float(p_joint_independence),
        "fold_deviation_from_independence": float(fold_deviation),
    }


def multi_marker_yield_gain(fit_result, n_vesicles=300_000, seed=3):
    """
    A second, directly practical Aim 3b question this fitted model can now
    answer: how much does capturing on EITHER marker (EpCAM OR Her2 --
    the simplest possible "combination" capture strategy, an OR-gate
    across markers) improve yield over the better of the two SINGLE-marker
    strategies alone? This is the "insurance if single-marker ceilings are
    low" argument from the handoff's "Multi-marker capture is a core
    sub-aim" section, made concrete and numeric for the first time using
    an actual fitted population rather than an assumed independent-OR
    calculation.

    Parameters
    ----------
    fit_result : dict
        Output of fit_both_markers().
    n_vesicles, seed :
        Passed through to simulate_joint_population().

    Returns
    -------
    dict with keys:
        "p_epcam_positive", "p_her2_positive" : single-marker yields
        "p_either_positive_true" : true P(EpCAM+ OR Her2+) from the joint population
        "p_either_positive_independence_prediction" : the naive
            independence-based OR formula, 1-(1-p_epcam)*(1-p_her2)
        "best_single_marker_yield" : max(p_epcam, p_her2)
        "combination_gain_true" : p_either_positive_true / best_single_marker_yield
            (how much the REAL combination strategy improves over the
            best single marker)
    """
    epcam_info = fit_result["EpCAM"]
    her2_info = fit_result["Her2"]

    pop = simulate_joint_population(
        f_rich=fit_result["f_rich"],
        epcam_mean_tet_poor=epcam_info["mean_tet_poor"],
        epcam_mean_tet_rich=epcam_info["mean_tet_rich"],
        her2_mean_tet_poor=her2_info["mean_tet_poor"],
        her2_mean_tet_rich=her2_info["mean_tet_rich"],
        n_vesicles=n_vesicles,
        seed=seed,
    )
    epcam_pos = pop["epcam"] > 0
    her2_pos = pop["her2"] > 0

    p_epcam = epcam_pos.mean()
    p_her2 = her2_pos.mean()
    p_either_true = (epcam_pos | her2_pos).mean()
    p_either_independence = 1.0 - (1.0 - p_epcam) * (1.0 - p_her2)
    best_single = max(p_epcam, p_her2)

    return {
        "p_epcam_positive": float(p_epcam),
        "p_her2_positive": float(p_her2),
        "p_either_positive_true": float(p_either_true),
        "p_either_positive_independence_prediction": float(p_either_independence),
        "best_single_marker_yield": float(best_single),
        "combination_gain_true": float(p_either_true / best_single) if best_single > 0 else float("inf"),
    }


def combination_gain_sensitivity_sweep(mean_tet_poor_values=(2.0, 5.0, 10.0, 15.0, 20.0, 30.0),
                                         f_rich=0.20, seed=3):
    """
    Sweep MARKER_MEAN_TET_POOR (the absolute copy-number scale, which
    Mizenko's ratio-only data cannot itself constrain -- see module
    docstring) and report both independence_test() and
    multi_marker_yield_gain() at each value, so the write-up can show
    directly which of this file's two headline results is robust to this
    unconstrained choice (the independence fold-deviation) and which is
    NOT (the combination-gain multiplier) -- rather than reporting only
    one arbitrarily-chosen-scale number for the latter as if it were as
    solid as the former.

    Parameters
    ----------
    mean_tet_poor_values : tuple of float
        Absolute tet-poor mean copy-number values to test.
    f_rich : float
        Held fixed at a representative mid-range value for this sweep
        (the f_rich-dependence is already characterized separately by
        fit_across_f_rich_range()).
    seed : int
        Passed through to the fits and simulations.

    Returns
    -------
    list of dict
        One entry per mean_tet_poor value:
        {"mean_tet_poor": ..., "best_single_marker_yield": ...,
         "combination_gain_true": ..., "fold_deviation_from_independence": ...}
    """
    global MARKER_MEAN_TET_POOR
    original = MARKER_MEAN_TET_POOR
    rows = []
    try:
        for mean_poor in mean_tet_poor_values:
            MARKER_MEAN_TET_POOR = mean_poor
            fit = fit_both_markers(f_rich, seed=seed)
            gain = multi_marker_yield_gain(fit, seed=seed)
            indep = independence_test(fit, seed=seed)
            rows.append({
                "mean_tet_poor": mean_poor,
                "best_single_marker_yield": gain["best_single_marker_yield"],
                "combination_gain_true": gain["combination_gain_true"],
                "fold_deviation_from_independence": indep["fold_deviation_from_independence"],
            })
    finally:
        MARKER_MEAN_TET_POOR = original  # always restore the module-level
                                          # default, even if an exception
                                          # occurs mid-sweep, so this
                                          # function has no lasting side
                                          # effect on the module's state.
    return rows


if __name__ == "__main__":
    print(__doc__.strip().split("\n\n")[0])
    print()
    print("Fitting the two-subpopulation mixture model DIRECTLY to Mizenko "
          "et al. 2021's real reported ratios (EpCAM ~10x, Her2 ~5x) -- "
          "not just sweeping illustrative configurations as the scoping "
          "file (mizenko_cooccurrence_scoping.py) did.\n")

    print(f"{'f_rich':>8} | {'EpCAM suppression':>18} | {'EpCAM achieved':>15} | "
          f"{'Her2 suppression':>18} | {'Her2 achieved':>14}")
    print("-" * 84)
    fits = fit_across_f_rich_range()
    for fit in fits:
        e = fit["EpCAM"]
        h = fit["Her2"]
        print(f"{fit['f_rich']:>8.2f} | {e['suppression_ratio']:>18.2f}x | "
              f"{e['achieved_ratio']:>14.2f}x | {h['suppression_ratio']:>18.2f}x | "
              f"{h['achieved_ratio']:>14.2f}x")

    print("\nTAKEAWAY: for each tested tetraspanin-rich fraction, a specific, "
          "fittable suppression ratio for EACH marker reproduces Mizenko's "
          "real reported ratio closely -- confirming the scoping file's "
          "qualitative finding with an actual fit, not just an illustrative "
          "sweep. Smaller f_rich requires LESS extreme per-marker "
          "suppression to reach the same target ratio (fewer 'hidden' "
          "vesicles need to carry a larger fraction of the signal); larger "
          "f_rich requires MORE extreme suppression.")

    print(f"\n{'='*74}")
    print("INDEPENDENCE TEST ('Immediate next steps' item 5) -- does the "
          "naive P(EpCAM+)*P(Her2+) shortcut hold?")
    print(f"{'='*74}")
    mid_fit = fit_both_markers(0.20)  # a representative mid-range f_rich
    indep = independence_test(mid_fit)
    print(f"  P(EpCAM+) = {indep['p_epcam_positive']:.4f}")
    print(f"  P(Her2+)  = {indep['p_her2_positive']:.4f}")
    print(f"  TRUE P(EpCAM+ AND Her2+)              = {indep['p_joint_true']:.4f}")
    print(f"  NAIVE INDEPENDENCE prediction (product) = {indep['p_joint_independence_prediction']:.4f}")
    print(f"  Fold deviation from independence: {indep['fold_deviation_from_independence']:.2f}x")
    if indep["fold_deviation_from_independence"] > 1.1:
        print("  FINDING: EpCAM and Her2 positivity are POSITIVELY "
              "CORRELATED across vesicles (more than independence would "
              "predict) -- both markers are suppressed in the SAME "
              "tetraspanin-rich subpopulation, so 'clustering,' exactly as "
              "the handoff's 'independence assumption is itself a finding' "
              "section anticipated, is real in this fitted model. The "
              "naive multiplication shortcut UNDER-predicts true double-"
              "positive yield.")
    elif indep["fold_deviation_from_independence"] < 0.9:
        print("  FINDING: EpCAM and Her2 positivity are NEGATIVELY "
              "correlated -- the naive multiplication shortcut OVER-"
              "predicts true double-positive yield.")
    else:
        print("  FINDING: independence is a reasonable approximation for "
              "this fitted model at this f_rich.")

    print(f"\n{'='*74}")
    print("MULTI-MARKER YIELD GAIN -- does combining EpCAM+Her2 (OR-gate "
          "capture) meaningfully beat the better single marker? AND: how "
          "sensitive is this to the arbitrary absolute copy-number scale?")
    print(f"{'='*74}")
    gain = multi_marker_yield_gain(mid_fit)
    print(f"  At the default MARKER_MEAN_TET_POOR={MARKER_MEAN_TET_POOR}:")
    print(f"    Best single-marker yield: {gain['best_single_marker_yield']:.4f}")
    print(f"    TRUE combination (EpCAM OR Her2) yield: {gain['p_either_positive_true']:.4f}")
    print(f"    Naive independence-OR prediction: {gain['p_either_positive_independence_prediction']:.4f}")
    print(f"    Real combination gain over best single marker: "
          f"{gain['combination_gain_true']:.2f}x")

    print(f"\n  IMPORTANT, checked directly (not assumed) -- see module "
          f"docstring's sensitivity-finding section: this combination-gain "
          f"number is HIGHLY SENSITIVE to MARKER_MEAN_TET_POOR, an absolute "
          f"scale that Mizenko's ratio-only data cannot itself constrain, "
          f"while the independence fold-deviation is NOT:\n")
    print(f"  {'mean_tet_poor':>14} | {'best single':>12} | {'combo gain':>11} | {'indep. fold-dev':>16}")
    print("  " + "-" * 62)
    for row in combination_gain_sensitivity_sweep():
        print(f"  {row['mean_tet_poor']:>14.1f} | {row['best_single_marker_yield']:>12.4f} | "
              f"{row['combination_gain_true']:>10.2f}x | {row['fold_deviation_from_independence']:>15.2f}x")

    print(f"\n  HONEST TAKEAWAY: the combination-gain multiplier ranges from "
          f"~1.02x to ~1.26x across this sweep -- BELOW the handoff's own "
          f"prior expectation of '2-4x, not 20x,' and highly dependent on "
          f"an absolute copy-number scale this fit cannot pin down. The "
          f"independence fold-deviation (~1.18x, positive correlation), by "
          f"contrast, is STABLE across the same sweep -- it is the more "
          f"trustworthy of this file's two headline numbers. Reporting "
          f"'combination gain is somewhere between 1x and 1.3x, and depends "
          f"on the (currently unmeasured) absolute marker density' is the "
          f"honest version of this finding -- not a single point estimate.")
