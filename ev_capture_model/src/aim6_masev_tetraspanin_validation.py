"""
aim6_masev_tetraspanin_validation.py

WHY THIS FILE EXISTS
---------------------
This is the project's pivot module (decided in the session documented in
docs/Pivot_Decision_and_MASEV_Validation_Scoping_Buildlog.md, and pointed to
directly by EV-Capture-CURRENT.md's "Immediate next steps," item 0). Every
finding this project has produced up to this point has been checked
against ITSELF -- analytical vs. Monte Carlo agreement, sensitivity sweeps,
robustness checks -- but never against a real measurement made by someone
else, on real vesicles, that this project did not produce. This file is
the first place that happens.

THE REAL DATA THIS FILE VALIDATES AGAINST
--------------------------------------------
Spitzberg, Ferguson, Yang, Peterson, Carlson & Weissleder, "Multiplexed
analysis of EV reveals specific biomarker composition with diagnostic
impact" (the "MASEV" method), Nat Commun 14:1239 (2023), DOI
10.1038/s41467-023-36932-z. This paper tracks the SAME individual vesicle
across 5 cycles of stain/image/cleave-and-restain (image-registration-
based identity tracking, >99% retention validated) -- genuine per-vesicle
JOINT multi-marker data, not population summary statistics reconstructed
from two separate, unlinked measurements (unlike this project's own
EGFR/CA19-9 joint reconstruction in aim5_and_gate_dual_marker_selectivity.py,
which is honest about being built from two DIFFERENT scatter plots that
never measured both markers on the same vesicle).

The paper's Figure 4 reports, for each of several cell lines, the exact
percentage of vesicles falling into each of the 8 possible combinations of
CD9/CD63/CD81 presence-absence (0, 1, 2, or all 3 tetraspanins positive,
broken out by WHICH ones). Real numbers for two cell lines (the only two
with a full 8-category breakdown available this session -- see "Known
risks" in EV-Capture-CURRENT.md for why ASPC1/A549 are not included) are
hard-coded in REAL_MASEV_DATA below, transcribed directly from the paper's
own figure.

WHAT THIS FILE ASKS, PRECISELY
----------------------------------
Can this project's existing joint-reconstruction machinery -- the SAME
Gaussian-copula approach already built and tested in
aim5_and_gate_dual_marker_selectivity.py for EGFR/CA19-9 -- predict the
correct SHAPE of a real joint 3-marker distribution, when its three
marginal inputs (P(CD9+), P(CD63+), P(CD81+)) are set to match the real
paper's own reported marginals for that cell line?

This is NOT a free-parameter fit in the sense of "tune parameters until
the numbers match." The three marginals are read directly off the real
data (see marginals_from_real_data() below -- they are not assumed or
guessed). The ONE thing this file sweeps, rather than fits, is the
pairwise correlation between the three tetraspanins' presence across
vesicles -- exactly the same "swept, not fit" discipline this project
already applies to K_INTRA/z/k_geom (k_intra_z_plausibility_bounds.py)
and to Aim 3b's f_rich (aim3b_cooccurrence_model.py), for the same
reason: there is no independent measurement in this paper (or anywhere
else checked so far) of what the TRUE pairwise tetraspanin correlation is,
only the resulting joint category percentages that a given correlation
would produce. Reporting a swept range, with an honest goodness-of-fit
number at each swept value, is the right level of confidence for what
this data can actually support -- not a single cherry-picked "best-fit"
correlation presented as if it were measured.

THE SCOPE THIS VALIDATES -- AND, EXPLICITLY, WHAT IT DOES NOT
------------------------------------------------------------------
This validates: "this project's joint multi-marker reconstruction
approach, when given the correct real marginals, predicts a joint
distribution SHAPE close to a real one, for markers with no tumor-vs-
normal biology story." CD9/CD63/CD81 are generic EV-membrane structural
proteins (tetraspanins), not tumor markers -- there is no "tumor vesicles
should have more CD9" story here the way there is for EGFR/CA19-9.

This does NOT validate: "this project's tumor-vs-normal selectivity
claims are now supported by real data." Those claims (the CA19-9
healthy/PDAC-plasma flip, the AND-gate dual-marker selectivity finding,
etc.) remain exactly as unvalidated/range-based as they were before this
file existed. Do not blur these two claims together in any future
write-up -- see EV-Capture-CURRENT.md's "Known risks" section, which
states this same warning for the pivot decision as a whole.

THE GOODNESS-OF-FIT METRIC USED, AND WHY
---------------------------------------------
Total variation distance (TVD): half the sum of |predicted_pct -
real_pct| across the 8 categories, expressed as a fraction (0 = perfect
match, 1 = the two distributions share no probability mass at all). TVD
is used (rather than, say, chi-square) because: (a) it has an intuitive
meaning directly in the same units as the percentages being compared --
"the predicted and real distributions differ by X percentage points of
vesicles, summed appropriately" -- which is easy to report honestly to a
non-expert audience (see this project's own "plain language" standing
requirement); (b) it does not require an assumed sample-size-dependent
null distribution the way a chi-square GOODNESS-of-fit p-value would,
which would invite an unjustified claim of statistical significance from
what are, after all, published percentages without the underlying raw
counts (see EV-Capture-CURRENT.md's "Known risks": no raw per-vesicle
data is available from this paper, only Figure 4's summary percentages).
"""

import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))


# ---------------------------------------------------------------------------
# PART 1: the real data itself, transcribed directly from Spitzberg et al.
# 2023, Figure 4. Percentages are fractions of the total vesicle population
# for that cell line (all 8 categories for a given cell line sum to 100%,
# checked directly in tests).
# ---------------------------------------------------------------------------

# Category order, fixed throughout this file: (CD9, CD63, CD81) presence
# tuples, in the exact order the real data below is given.
CATEGORY_ORDER = (
    "none",       # no tetraspanin positive
    "CD9_only",
    "CD81_only",
    "CD63_only",
    "CD81_CD9",   # CD9 and CD81 positive, CD63 negative
    "CD63_CD81",  # CD63 and CD81 positive, CD9 negative
    "CD63_CD9",   # CD63 and CD9 positive, CD81 negative
    "all_three",
)

# Which markers are positive in each category -- used to compute each
# marker's own marginal positivity fraction directly from the real data
# (see marginals_from_real_data() below), and to build the corresponding
# category from a simulated joint population on equal footing.
CATEGORY_MARKERS_POSITIVE = {
    "none":       {"CD9": False, "CD63": False, "CD81": False},
    "CD9_only":   {"CD9": True,  "CD63": False, "CD81": False},
    "CD81_only":  {"CD9": False, "CD63": False, "CD81": True},
    "CD63_only":  {"CD9": False, "CD63": True,  "CD81": False},
    "CD81_CD9":   {"CD9": True,  "CD63": False, "CD81": True},
    "CD63_CD81":  {"CD9": False, "CD63": True,  "CD81": True},
    "CD63_CD9":   {"CD9": True,  "CD63": True,  "CD81": False},
    "all_three":  {"CD9": True,  "CD63": True,  "CD81": True},
}

# Real percentages, Spitzberg et al. 2023, Figure 4 -- read directly at
# full resolution from the bioRxiv preprint's figure (see
# Pivot_Decision_and_MASEV_Validation_Scoping_Buildlog.md for the full
# retrieval account). Expressed as fractions (0-1), not percentages, for
# direct use in probability arithmetic below.
REAL_MASEV_DATA = {
    "PANC-1": {
        "none": 0.248, "CD9_only": 0.376, "CD81_only": 0.046, "CD63_only": 0.044,
        "CD81_CD9": 0.112, "CD63_CD81": 0.016, "CD63_CD9": 0.106, "all_three": 0.051,
    },
    "CAPAN-2": {
        "none": 0.526, "CD9_only": 0.237, "CD81_only": 0.075, "CD63_only": 0.048,
        "CD81_CD9": 0.044, "CD63_CD81": 0.014, "CD63_CD9": 0.035, "all_three": 0.020,
    },
}

# Number of vesicles analyzed per cell line in the real MASEV experiment
# (Spitzberg et al. 2023) -- not used in any calculation below, kept here
# purely for reporting/context (a reader may reasonably want to know the
# real sample size behind these percentages).
REAL_MASEV_N_VESICLES = {
    "PANC-1": 8500,
    "CAPAN-2": 4600,
    "ASPC1": 14_000,
    "A549": 8_100,
}

# ---------------------------------------------------------------------------
# NEW: single-marker vs. cocktail capture-loss numbers -- real, already-
# published data from THIS SAME PAPER (Spitzberg et al. 2023), never
# previously extracted into this codebase. Added following the second
# diagnosis-followthrough investigation (docs/Diagnosis_Followthrough_and_
# Recommendation.md, Part 6.2 item 2, and Category 19's finding that mining
# the project's own already-downloaded source papers was the single highest-
# value, lowest-risk category of the entire re-investigation).
#
# The paper states directly (in the same text discussing Figure 4/6A) that
# purifying EVs using a SINGLE tetraspanin antibody misses a large fraction
# of KRAS-mutant EVs -- up to 80% lost -- while a COCKTAIL of multiple
# tetraspanin antibodies still misses a real, substantial fraction (36-47%,
# depending on cell line/condition). This is real, quantified evidence that
# multi-marker capture strategies outperform single-marker capture but do
# NOT eliminate the problem -- directly relevant context for this project's
# own Aim 5 (the EGFR/CA19-9 AND-gate dual-marker selectivity model), which
# makes the qualitative case that combining markers changes what gets
# captured, without previously citing any real quantified example of how
# much loss single-marker capture actually causes in practice.
#
# SCOPE NOTE: these are KRAS-mutant-EV recovery-loss percentages (a
# different marker/mechanism than the tumor-vs-normal tetraspanin
# comparisons Aim 6 validates), so this is cited as independent
# corroborating context for the general "single-marker capture loses real
# signal, multi-marker capture recovers some but not all of it" story --
# not blended into Aim 6's own TVD-based validation numbers, and not used
# in any calculation in this file.
SPITZBERG_SINGLE_MARKER_MAX_LOSS_FRACTION = 0.80
SPITZBERG_COCKTAIL_MISS_FRACTION_RANGE = (0.36, 0.47)

# ---------------------------------------------------------------------------
# HELD-OUT PREDICTION DATA -- ASPC1 and A549 (NEW, added after the user
# manually downloaded the paper's Supplementary Information PDF, which
# TWO prior automated fetch attempts could not retrieve).
# ---------------------------------------------------------------------------
"""
WHY THIS SECTION EXISTS, AND WHAT IT DOES AND DOES NOT LET US DO
--------------------------------------------------------------------
Reading the retrieved Supplementary Information PDF directly (Supplementary
Fig. 12, a bar chart of single-marker EV positivity rates across 6 cell
lines, and Supplementary Fig. 15, a per-vesicle heatmap across 4 cell
lines) resolved PART of the original gap, not all of it:

- Supplementary Fig. 12 DOES give real, directly-readable single-marker
  MARGINAL positivity percentages for CD9, CD63, and CD81 in ASPC1 and
  A549 (read directly off the bar chart's y-axis, at the resolution a
  bar chart honestly supports -- to the nearest percentage point, not
  sub-percent precision). These are hard-coded below as
  ASPC1_A549_REAL_MARGINALS.

- Supplementary Fig. 15 does NOT give a numeric joint 8-category
  breakdown (the ASPC1/A549 equivalent of Figure 4's PANC-1/CAPAN-2
  numbers). It is a raw per-vesicle heatmap image (12,000 individual
  vesicles as rows, markers as columns, presence/absence shown by pixel
  shading) with no axis giving category percentages -- genuinely not
  transcribable into the same 8-category format REAL_MASEV_DATA uses,
  confirmed by inspecting the image directly rather than assumed. No
  other table or figure in the 19-page Supplementary Information gives
  this breakdown for ASPC1/A549 either (checked directly, page by page).

WHAT THIS MEANS FOR THE "HELD-OUT PREDICTION" STEP
-------------------------------------------------------
Because only the marginals (not the joint 8-category numbers) are
available for ASPC1/A549, this file cannot do for them what it did for
PANC-1/CAPAN-2 -- there is no real joint distribution to compute a TVD
against. What CAN honestly be done, and what
held_out_prediction_for_cell_line() below does:

1. Fit the shared copula correlation using ONLY the TRAINING data
   (PANC-1 and CAPAN-2's real joint 8-category numbers) -- i.e., commit
   to ONE correlation value derived entirely from the two cell lines
   with full joint data, rather than re-sweeping correlation separately
   for ASPC1/A549 (which would silently let the "prediction" cheat by
   picking whatever correlation looks best for the held-out cell line,
   defeating the point of a held-out test).
2. Feed that ONE committed correlation, together with ASPC1/A549's OWN
   real marginals (which ARE available), into the same joint
   reconstruction machinery, producing a full predicted 8-category
   distribution for each.
3. Report that prediction plainly as NOT independently checkable against
   a real joint number, because no such real number exists in this
   paper for these two cell lines.
4. As the one genuinely checkable piece: recompute the MARGINALS implied
   by the predicted joint distribution (by summing the right predicted
   categories, exactly as marginals_from_real_data() does for real data)
   and compare them to the real observed marginals from Supplementary
   Fig. 12. This is a real, if narrower, held-out consistency check --
   it verifies the predicted joint distribution's own internal marginals
   land close to reality, even though the full joint SHAPE prediction
   itself has no real number to check against directly.

This is a strictly weaker claim than the PANC-1/CAPAN-2 validation, and
should always be described that way: "the model's committed correlation,
applied out-of-sample to two new cell lines' real marginals, produces
marginally-consistent predictions" is NOT the same claim as "the model's
full joint-distribution prediction was checked against real joint data" --
only the latter is true for PANC-1/CAPAN-2, only the former is available
for ASPC1/A549. Do not blur these two strengths of evidence together.
"""

# Real single-marker positivity percentages, read directly off
# Supplementary Fig. 12's bar chart (y-axis: "% of EV pos for a
# biomarker"). Read to the nearest percentage point, the resolution a
# bar chart on a 0-70% axis honestly supports -- not claimed as more
# precise than that. Expressed as fractions (0-1) for direct use here.
ASPC1_A549_REAL_MARGINALS = {
    "ASPC1": {"CD9": 0.195, "CD63": 0.07, "CD81": 0.07},
    "A549":  {"CD9": 0.36,  "CD63": 0.14, "CD81": 0.03},
}

# Shared negative-binomial dispersion parameter, matching this project's
# standing convention elsewhere (aim3b_cooccurrence_model.py's
# simulate_joint_population, aim5_and_gate_dual_marker_selectivity.py's
# fit_marginal) for "moderate overdispersion" when no data-specific
# dispersion estimate is available -- which is the case here, since the
# real data gives only positivity fractions, not full copy-number
# distributions, so dispersion itself cannot be fit from this data any
# more than it could be in either of those two prior files.
DEFAULT_DISPERSION_R = 2.0

# Correlation values swept across all three marker pairs simultaneously
# (a single shared pairwise correlation, not three independent ones -- see
# "WHY A SINGLE SHARED CORRELATION, NOT THREE INDEPENDENT ONES" in
# simulate_joint_tetraspanin_population()'s docstring for why this
# simplification is used and what it does not claim).
CORRELATION_SWEEP = (0.0, 0.2, 0.4, 0.6, 0.8)


def marginals_from_real_data(cell_line):
    """
    Compute each tetraspanin's own marginal positivity fraction --
    P(CD9+), P(CD63+), P(CD81+) -- directly from the real 8-category
    Figure 4 breakdown, by summing the categories in which that marker is
    positive. This is arithmetic on the real published numbers, not a fit
    or an assumption -- e.g. P(CD9+) = CD9_only + CD81_CD9 + CD63_CD9 +
    all_three.

    Parameters
    ----------
    cell_line : str
        "PANC-1" or "CAPAN-2" (keys of REAL_MASEV_DATA).

    Returns
    -------
    dict with keys "CD9", "CD63", "CD81", each the real marginal
    positivity fraction (0-1) for that marker in that cell line.
    """
    data = REAL_MASEV_DATA[cell_line]
    marginals = {}
    for marker in ("CD9", "CD63", "CD81"):
        marginals[marker] = sum(
            frac for category, frac in data.items()
            if CATEGORY_MARKERS_POSITIVE[category][marker]
        )
    return marginals


# ---------------------------------------------------------------------------
# PART 2: fit a negative-binomial marginal to each marker's real positivity
# fraction (reusing the same one-free-parameter zero-fraction-matching
# approach as stem_cell_marker_ranking_extension.py's
# fit_negative_binomial_to_positivity(), specialized here since this data
# gives ONLY a positivity fraction per marker -- no separate "mean copies
# if positive" number exists in this paper the way it does for the stem-
# cell sources, so dispersion is held at DEFAULT_DISPERSION_R and only the
# mean is solved for, rather than solving for both mean and dispersion
# jointly as that file's two-constraint version does).
# ---------------------------------------------------------------------------

def fit_marginal_to_positivity_fraction(fraction_positive, dispersion_r=DEFAULT_DISPERSION_R):
    """
    Find the negative-binomial mean (at a FIXED dispersion_r) whose
    resulting P(n=0) matches (1 - fraction_positive) -- i.e., whose
    resulting P(n>0) matches the real marginal positivity fraction for
    one tetraspanin marker.

    WHY DISPERSION IS HELD FIXED HERE, UNLIKE stem_cell_marker_ranking_
    extension.py's TWO-CONSTRAINT VERSION
    ------------------------------------------------------------------------
    That file has TWO real summary statistics per marker (a positivity
    fraction AND a separate mean-copies-if-positive number), so it can
    solve for both the mean and the dispersion jointly. This paper's
    Figure 4 gives only a categorical presence/absence breakdown -- ONE
    real number per marker (the positivity fraction) -- so only one
    negative-binomial parameter can be solved for without inventing a
    second piece of information this data doesn't contain. Dispersion is
    therefore fixed at this project's standing default (2.0, "moderate
    overdispersion" -- see DEFAULT_DISPERSION_R above) and only the mean
    is fit, exactly the same "one real constraint, one free parameter"
    discipline aim3b_cooccurrence_model.py already applies when Mizenko's
    data likewise supplies only a ratio, not a full distribution.

    Parameters
    ----------
    fraction_positive : float
        Real P(marker+) for this marker/cell line, in (0, 1).
    dispersion_r : float
        Fixed negative-binomial dispersion (default: this project's
        standard 2.0).

    Returns
    -------
    dict with "mean_n", "dispersion_r".
    """
    if not (0.0 < fraction_positive < 1.0):
        raise ValueError(
            f"fraction_positive must be strictly between 0 and 1, got "
            f"{fraction_positive} -- a boundary value would make the "
            f"zero-fraction root-find degenerate."
        )

    target_p0 = 1.0 - fraction_positive

    def p0_error(mean_n):
        p_param = dispersion_r / (dispersion_r + mean_n)
        return stats.nbinom.pmf(0, dispersion_r, p_param) - target_p0

    # Bracket: mean_n=1e-6 (P(n=0) -> ~1) to mean_n=1000 (P(n=0) -> ~0).
    # Comfortably brackets any target_p0 in (0, 1) for dispersion_r in the
    # range this project uses -- checked directly in this file's test
    # suite across all real marginals actually used here.
    from scipy.optimize import brentq
    mean_n = brentq(p0_error, 1e-6, 1000.0, xtol=1e-8)

    return {"mean_n": mean_n, "dispersion_r": dispersion_r}


def fit_all_marginals_for_cell_line(cell_line, dispersion_r=DEFAULT_DISPERSION_R):
    """
    Fit all three tetraspanin marginals (CD9, CD63, CD81) for one cell
    line, using the real positivity fractions from marginals_from_real_data().

    Returns
    -------
    dict keyed by marker name ("CD9", "CD63", "CD81") -> fit dict (see
    fit_marginal_to_positivity_fraction).
    """
    marginals = marginals_from_real_data(cell_line)
    return {
        marker: fit_marginal_to_positivity_fraction(frac, dispersion_r=dispersion_r)
        for marker, frac in marginals.items()
    }


# ---------------------------------------------------------------------------
# PART 3: simulate a joint (CD9, CD63, CD81) population via a 3-variable
# Gaussian copula -- the direct 3-marker extension of aim5's 2-marker
# copula (simulate_joint_two_marker_population()), reusing the exact same
# construction technique.
# ---------------------------------------------------------------------------

def simulate_joint_tetraspanin_population(marginal_fits, correlation,
                                             n_vesicles=100_000, seed=0):
    """
    Simulate a population of vesicles with jointly-drawn (CD9, CD63, CD81)
    counts, via a 3-variable Gaussian copula at a single SHARED pairwise
    correlation, anchored to the three real marginal fits.

    WHY A SINGLE SHARED CORRELATION, NOT THREE INDEPENDENT PAIRWISE ONES
    --------------------------------------------------------------------------
    A full 3x3 correlation matrix would need three independent pairwise
    correlations (CD9-CD63, CD9-CD81, CD63-CD81), each itself an
    unmeasured free parameter -- this paper's Figure 4 does not report
    pairwise correlations directly, only the resulting joint category
    percentages, so there is no more basis for fitting three independent
    correlations than for fitting one. Using a SINGLE shared correlation
    for all three pairs is the simplest version of "these markers are all
    biogenesis-linked structural tetraspanins, plausibly co-varying
    together for a shared reason (e.g. overall EV membrane-protein
    density, or exosome-vs-microvesicle subtype identity), not via three
    independently-different pairwise relationships" -- an explicit
    simplifying assumption, stated here rather than silently built in,
    exactly the same honesty standard this project applies to aim5's
    "identical correlation assumed for tumor and normal populations"
    simplification (see that file's own docstring, "Immediate next
    steps"). If this validation's fit turns out to be poor, unequal
    pairwise correlations are a natural, explicitly-named next thing to
    try -- not silently assumed away.

    correlation=0 reproduces independence exactly (checked directly in
    this file's test suite, the same way aim5's copula's independence
    case is checked).

    Parameters
    ----------
    marginal_fits : dict
        Output of fit_all_marginals_for_cell_line() -- keyed by "CD9",
        "CD63", "CD81", each a fit dict with "mean_n", "dispersion_r".
    correlation : float
        Shared pairwise Gaussian-copula correlation, in [-1, 1]. Applied
        identically to all three marker pairs (see above).
    n_vesicles : int
        Simulated population size.
    seed : int
        RNG seed.

    Returns
    -------
    dict with "CD9", "CD63", "CD81" (int arrays, length n_vesicles) --
    the simulated per-vesicle copy numbers for each marker, jointly drawn.
    """
    if not (-1.0 <= correlation <= 1.0):
        raise ValueError(f"correlation must be in [-1, 1], got {correlation}")

    rng = np.random.default_rng(seed)
    mean_vec = [0.0, 0.0, 0.0]
    cov = [
        [1.0, correlation, correlation],
        [correlation, 1.0, correlation],
        [correlation, correlation, 1.0],
    ]
    z = rng.multivariate_normal(mean_vec, cov, size=n_vesicles)
    u = stats.norm.cdf(z)  # shape (n_vesicles, 3), columns = CD9, CD63, CD81

    counts = {}
    for i, marker in enumerate(("CD9", "CD63", "CD81")):
        fit = marginal_fits[marker]
        p_param = fit["dispersion_r"] / (fit["dispersion_r"] + fit["mean_n"])
        counts[marker] = stats.nbinom.ppf(u[:, i], fit["dispersion_r"], p_param).astype(float)

    return counts


def categorize_joint_population(counts):
    """
    Convert a simulated joint population (output of
    simulate_joint_tetraspanin_population(), or any dict with "CD9",
    "CD63", "CD81" per-vesicle count arrays) into the same 8-category
    breakdown REAL_MASEV_DATA uses, as fractions of the total population.

    A marker is treated as "positive" if its simulated count is > 0 --
    matching exactly how the real paper's own presence/absence categories
    are defined (a binary detected/not-detected call per marker per
    vesicle, not a copy-number threshold), and matching how
    marginals_from_real_data() computes real marginals from the same real
    category breakdown.

    Parameters
    ----------
    counts : dict
        "CD9", "CD63", "CD81" -> per-vesicle count arrays (same shape).

    Returns
    -------
    dict keyed by category name (CATEGORY_ORDER) -> fraction of the
    simulated population in that category. Sums to 1.0 by construction
    (every vesicle falls into exactly one of the 8 mutually exclusive,
    collectively exhaustive categories).
    """
    cd9_pos = counts["CD9"] > 0
    cd63_pos = counts["CD63"] > 0
    cd81_pos = counts["CD81"] > 0

    n_total = len(cd9_pos)
    fractions = {}
    for category, markers in CATEGORY_MARKERS_POSITIVE.items():
        mask = np.ones(n_total, dtype=bool)
        mask &= cd9_pos if markers["CD9"] else ~cd9_pos
        mask &= cd63_pos if markers["CD63"] else ~cd63_pos
        mask &= cd81_pos if markers["CD81"] else ~cd81_pos
        fractions[category] = float(mask.sum()) / n_total
    return fractions


# ---------------------------------------------------------------------------
# PART 4: goodness-of-fit against the real data, and the correlation sweep.
# ---------------------------------------------------------------------------

def total_variation_distance(predicted, real):
    """
    Total variation distance between two category-fraction dicts (same
    keys): half the sum of absolute differences. 0 = identical
    distributions; 1 = no shared probability mass at all. See module
    docstring "THE GOODNESS-OF-FIT METRIC USED, AND WHY" for why this
    metric was chosen over a chi-square-style test statistic.

    Parameters
    ----------
    predicted, real : dict
        Category fraction dicts, same keys (CATEGORY_ORDER).

    Returns
    -------
    float
        Total variation distance, in [0, 1].
    """
    return 0.5 * sum(abs(predicted[cat] - real[cat]) for cat in CATEGORY_ORDER)


def validate_cell_line(cell_line, correlations=CORRELATION_SWEEP,
                          n_vesicles=100_000, seed=0):
    """
    THE MAIN FUNCTION THIS FILE EXISTS TO PROVIDE: for one real cell line,
    fit the three tetraspanin marginals to the real data, then -- across a
    swept range of assumed shared correlations -- simulate the joint
    population, categorize it into the same 8 real-data categories, and
    report the total variation distance against the real numbers.

    Parameters
    ----------
    cell_line : str
        "PANC-1" or "CAPAN-2".
    correlations : tuple of float
        Shared pairwise correlations to sweep (default: CORRELATION_SWEEP).
    n_vesicles, seed : int
        Passed through to the simulation.

    Returns
    -------
    dict with:
        "cell_line": echoed back.
        "real_marginals": output of marginals_from_real_data().
        "correlations": the swept correlation values.
        "predicted_categories": list of dict, one 8-category breakdown per
            swept correlation.
        "tvd": list of float, total variation distance at each correlation.
        "best_correlation": the correlation with the LOWEST tvd (reported
            as the swept value that fits best -- NOT claimed as a
            measured true correlation, see module docstring).
        "best_tvd": the tvd at best_correlation.
    """
    marginal_fits = fit_all_marginals_for_cell_line(cell_line)
    real_data = REAL_MASEV_DATA[cell_line]

    predicted_list = []
    tvd_list = []
    for correlation in correlations:
        counts = simulate_joint_tetraspanin_population(
            marginal_fits, correlation, n_vesicles=n_vesicles, seed=seed
        )
        predicted = categorize_joint_population(counts)
        predicted_list.append(predicted)
        tvd_list.append(total_variation_distance(predicted, real_data))

    best_idx = int(np.argmin(tvd_list))

    return {
        "cell_line": cell_line,
        "real_marginals": marginals_from_real_data(cell_line),
        "correlations": list(correlations),
        "predicted_categories": predicted_list,
        "tvd": tvd_list,
        "best_correlation": correlations[best_idx],
        "best_tvd": tvd_list[best_idx],
    }


def independence_baseline_tvd(cell_line, n_vesicles=100_000, seed=0):
    """
    Convenience function returning JUST the correlation=0 (independence)
    TVD for a cell line -- the conservative default this project reports
    everywhere else when no correlation measurement exists (see aim5's
    module docstring, "the independence case (correlation=0) is reported
    as the conservative default"). Kept as its own named function so the
    write-up can cite "the independence-only prediction misses real
    tetraspanin clustering by X" as a standalone, simple number, separate
    from the full correlation sweep.
    """
    result = validate_cell_line(cell_line, correlations=(0.0,), n_vesicles=n_vesicles, seed=seed)
    return result["tvd"][0]


# ---------------------------------------------------------------------------
# PART 5: the held-out prediction for ASPC1/A549 -- committing to ONE
# correlation from the training cell lines (PANC-1, CAPAN-2), then
# applying it out-of-sample.
# ---------------------------------------------------------------------------

def fit_committed_correlation(training_cell_lines=("PANC-1", "CAPAN-2"),
                                 correlations=CORRELATION_SWEEP,
                                 n_vesicles=100_000, seed=0):
    """
    Commit to ONE shared copula correlation value, derived ONLY from the
    training cell lines (the ones with real, full joint 8-category data
    available) -- the correlation that minimizes the AVERAGE total
    variation distance across both training cell lines simultaneously
    (not each cell line's own individually-best correlation, which would
    already be two different numbers -- 0.2 happens to be individually
    best for both PANC-1 and CAPAN-2 in the CORRELATION_SWEEP grid tested
    here, but this function computes the shared-average-minimizing choice
    explicitly rather than relying on that coincidence, so it remains
    correct if the swept grid or training set ever changes).

    THIS IS THE KEY DISCIPLINE THAT MAKES ASPC1/A549 A REAL HELD-OUT TEST
    --------------------------------------------------------------------------
    If ASPC1/A549 were each allowed their own re-swept, independently-
    best-fit correlation, the "prediction" for them would not be a
    prediction at all -- it would be a second in-sample fit, exactly like
    PANC-1/CAPAN-2's own validate_cell_line() sweep. Locking the
    correlation to a SINGLE value chosen using ONLY the training data,
    then applying that exact value to a cell line whose joint data was
    never used to choose it, is what makes the resulting ASPC1/A549
    prediction a genuine out-of-sample test rather than a second fit
    dressed up as one.

    Parameters
    ----------
    training_cell_lines : tuple of str
        Cell lines with real joint data to fit the shared correlation on
        (default: the two with full REAL_MASEV_DATA breakdowns).
    correlations : tuple of float
        Candidate correlation values to choose among (default:
        CORRELATION_SWEEP, the same grid validate_cell_line() sweeps).
    n_vesicles, seed : int
        Passed through to the simulation.

    Returns
    -------
    dict with:
        "committed_correlation": the single correlation value chosen.
        "average_tvd_on_training": the average TVD across
            training_cell_lines at this correlation.
        "per_cell_line_tvd_at_committed_correlation": dict, cell line ->
            its own TVD at the committed correlation (for transparency --
            shows whether the shared choice is a good fit for BOTH
            training cell lines individually, not just on average).
    """
    tvd_by_correlation = {corr: [] for corr in correlations}
    for cell_line in training_cell_lines:
        result = validate_cell_line(cell_line, correlations=correlations,
                                       n_vesicles=n_vesicles, seed=seed)
        for corr, tvd in zip(result["correlations"], result["tvd"]):
            tvd_by_correlation[corr].append(tvd)

    average_tvd = {corr: float(np.mean(tvds)) for corr, tvds in tvd_by_correlation.items()}
    committed_correlation = min(average_tvd, key=average_tvd.get)

    per_cell_line_tvd = {}
    for i, cell_line in enumerate(training_cell_lines):
        per_cell_line_tvd[cell_line] = tvd_by_correlation[committed_correlation][i]

    return {
        "committed_correlation": committed_correlation,
        "average_tvd_on_training": average_tvd[committed_correlation],
        "per_cell_line_tvd_at_committed_correlation": per_cell_line_tvd,
    }


def held_out_prediction_for_cell_line(cell_line, committed_correlation,
                                         n_vesicles=100_000, seed=0):
    """
    Predict the full joint 8-category tetraspanin distribution for a
    HELD-OUT cell line (ASPC1 or A549 -- one whose real joint data does
    not exist in this paper, only its real marginals), using the
    committed correlation (fit_committed_correlation()'s output, derived
    ONLY from the training cell lines) and this cell line's OWN real
    marginal positivity fractions (ASPC1_A549_REAL_MARGINALS).

    Then performs the ONE genuinely checkable comparison available: the
    predicted joint distribution's own implied marginals (recomputed by
    summing categories, exactly as marginals_from_real_data() does) are
    compared against the real observed marginals -- NOT because this is
    circular (the real marginals were an INPUT to the marginal fits, so
    of course the fitted marginal distributions' own means match by
    construction) but because the copula's categorization step
    (thresholding at count > 0 per marker) is a genuinely separate
    computation from the marginal fit itself, and a bug in the joint
    construction (e.g. a sign error in the correlation, a mismatched
    marker ordering) COULD make the recomputed marginals drift from the
    fitted targets even though the marginal fit itself was correct --
    see this file's own test suite for a check that exercises exactly
    this failure mode.

    Parameters
    ----------
    cell_line : str
        "ASPC1" or "A549" (keys of ASPC1_A549_REAL_MARGINALS).
    committed_correlation : float
        The single correlation value to apply (from
        fit_committed_correlation(), NOT re-fit for this cell line).
    n_vesicles, seed : int
        Passed through to the simulation.

    Returns
    -------
    dict with:
        "cell_line": echoed back.
        "committed_correlation": echoed back.
        "real_marginals": the real input marginals (from
            ASPC1_A549_REAL_MARGINALS).
        "predicted_categories": the predicted full 8-category joint
            distribution (NOT independently checkable -- no real joint
            number exists for this cell line, see module discussion).
        "predicted_marginals_recomputed": the marginals implied by the
            predicted joint distribution, recomputed by summing
            categories -- the one genuinely checkable output.
        "marginal_recovery_error": dict, marker -> absolute difference
            between real_marginals and predicted_marginals_recomputed
            (should be small; a real check on the joint-construction
            arithmetic, not a claim about the correlation's correctness).
    """
    if cell_line not in ASPC1_A549_REAL_MARGINALS:
        raise ValueError(
            f"cell_line must be one of {list(ASPC1_A549_REAL_MARGINALS.keys())} "
            f"(held-out cell lines with real marginals but no real joint "
            f"8-category data) -- got {cell_line!r}. Use validate_cell_line() "
            f"instead for PANC-1/CAPAN-2, which DO have real joint data."
        )

    real_marginals = ASPC1_A549_REAL_MARGINALS[cell_line]
    marginal_fits = {
        marker: fit_marginal_to_positivity_fraction(frac)
        for marker, frac in real_marginals.items()
    }

    counts = simulate_joint_tetraspanin_population(
        marginal_fits, committed_correlation, n_vesicles=n_vesicles, seed=seed
    )
    predicted_categories = categorize_joint_population(counts)

    predicted_marginals_recomputed = {}
    for marker in ("CD9", "CD63", "CD81"):
        predicted_marginals_recomputed[marker] = sum(
            frac for category, frac in predicted_categories.items()
            if CATEGORY_MARKERS_POSITIVE[category][marker]
        )

    marginal_recovery_error = {
        marker: abs(real_marginals[marker] - predicted_marginals_recomputed[marker])
        for marker in ("CD9", "CD63", "CD81")
    }

    return {
        "cell_line": cell_line,
        "committed_correlation": committed_correlation,
        "real_marginals": real_marginals,
        "predicted_categories": predicted_categories,
        "predicted_marginals_recomputed": predicted_marginals_recomputed,
        "marginal_recovery_error": marginal_recovery_error,
    }


def print_summary():
    print("=" * 78)
    print("AIM 6 -- MASEV TETRASPANIN JOINT-DISTRIBUTION VALIDATION")
    print("(first real external-data validation this project has ever had)")
    print("=" * 78)

    for cell_line in REAL_MASEV_DATA:
        print(f"\n--- {cell_line} (real N={REAL_MASEV_N_VESICLES[cell_line]} vesicles, "
              f"Spitzberg et al. 2023, Fig. 4) ---")
        real_marg = marginals_from_real_data(cell_line)
        print(f"  Real marginals: P(CD9+)={real_marg['CD9']:.3f}, "
              f"P(CD63+)={real_marg['CD63']:.3f}, P(CD81+)={real_marg['CD81']:.3f}")

        result = validate_cell_line(cell_line)
        print(f"  {'correlation':>12} | {'TVD':>8}")
        print("  " + "-" * 24)
        for corr, tvd in zip(result["correlations"], result["tvd"]):
            marker = "  <- best" if corr == result["best_correlation"] else ""
            print(f"  {corr:>12.1f} | {tvd:>8.4f}{marker}")

        print(f"\n  Independence-only (correlation=0) TVD: {result['tvd'][0]:.4f}")
        print(f"  Best swept correlation: {result['best_correlation']:.1f} "
              f"(TVD={result['best_tvd']:.4f})")
        print(f"\n  Category-by-category comparison at best-fit correlation "
              f"({result['best_correlation']:.1f}):")
        best_idx = result["correlations"].index(result["best_correlation"])
        predicted = result["predicted_categories"][best_idx]
        real = REAL_MASEV_DATA[cell_line]
        print(f"  {'category':>12} | {'real %':>8} | {'predicted %':>12} | {'diff (pp)':>10}")
        print("  " + "-" * 50)
        for cat in CATEGORY_ORDER:
            diff = (predicted[cat] - real[cat]) * 100
            print(f"  {cat:>12} | {real[cat]*100:>7.1f}% | {predicted[cat]*100:>11.1f}% | {diff:>+9.1f}")

    print(f"\n{'='*78}")
    print("HELD-OUT PREDICTION -- ASPC1 and A549")
    print("(real marginals available; real joint 8-category data is NOT --")
    print(" see module docstring, 'HELD-OUT PREDICTION DATA' section)")
    print("=" * 78)

    committed = fit_committed_correlation()
    print(f"\nCommitted correlation (fit using ONLY PANC-1/CAPAN-2's real joint "
          f"data): {committed['committed_correlation']:.1f}")
    print(f"  Average TVD on training cell lines at this correlation: "
          f"{committed['average_tvd_on_training']:.4f}")
    for cell_line, tvd in committed["per_cell_line_tvd_at_committed_correlation"].items():
        print(f"    {cell_line}: TVD={tvd:.4f}")

    print(f"\n{'='*78}")
    print("NEW: single-marker vs. cocktail capture-loss (same paper, real data,")
    print("never previously extracted -- see Diagnosis_Followthrough doc)")
    print("=" * 78)
    print(f"  Single tetraspanin antibody: up to "
          f"{SPITZBERG_SINGLE_MARKER_MAX_LOSS_FRACTION*100:.0f}% of KRAS-mutant EVs missed.")
    lo, hi = SPITZBERG_COCKTAIL_MISS_FRACTION_RANGE
    print(f"  Tetraspanin COCKTAIL (multi-marker): still misses {lo*100:.0f}-{hi*100:.0f}%.")
    print("  Real, quantified evidence that multi-marker capture helps but does")
    print("  not solve the problem -- relevant context for Aim 5's AND-gate model.")

    for cell_line in ASPC1_A549_REAL_MARGINALS:
        print(f"\n--- {cell_line} (real N={REAL_MASEV_N_VESICLES[cell_line]} vesicles) ---")
        held_out = held_out_prediction_for_cell_line(
            cell_line, committed["committed_correlation"]
        )
        print(f"  Real marginals (read from Supplementary Fig. 12): "
              f"CD9={held_out['real_marginals']['CD9']*100:.1f}%, "
              f"CD63={held_out['real_marginals']['CD63']*100:.1f}%, "
              f"CD81={held_out['real_marginals']['CD81']*100:.1f}%")
        print(f"  Predicted joint 8-category distribution (using the committed "
              f"correlation={held_out['committed_correlation']:.1f} -- NOT "
              f"independently checkable, no real joint number exists for "
              f"this cell line):")
        for cat in CATEGORY_ORDER:
            print(f"    {cat:>12}: {held_out['predicted_categories'][cat]*100:>5.1f}%")
        print(f"  Held-out consistency check -- predicted marginals recomputed "
              f"from the predicted joint distribution vs. the real input "
              f"marginals (this IS checkable):")
        for marker in ("CD9", "CD63", "CD81"):
            real_m = held_out["real_marginals"][marker] * 100
            pred_m = held_out["predicted_marginals_recomputed"][marker] * 100
            err = held_out["marginal_recovery_error"][marker] * 100
            print(f"    {marker}: real={real_m:.1f}%, recomputed-from-joint={pred_m:.1f}%, "
                  f"error={err:.2f} pp")

    print(f"\n{'='*78}")
    print("SCOPE REMINDER (read before citing this result anywhere):")
    print("This validates that the project's joint multi-marker reconstruction")
    print("approach can predict a real published joint-distribution SHAPE for")
    print("generic EV structural proteins (tetraspanins). It does NOT validate")
    print("this project's tumor-vs-normal selectivity claims (EGFR/CA19-9),")
    print("which remain exactly as range-based/unvalidated as before this file.")
    print()
    print("The ASPC1/A549 section above is a WEAKER form of evidence than the")
    print("PANC-1/CAPAN-2 validation: it checks that a correlation committed to")
    print("from OTHER cell lines produces internally-consistent marginals on")
    print("new cell lines, NOT that the full predicted joint shape matches a")
    print("real joint number (no such number exists in this paper for these")
    print("two cell lines). Keep these two claims separate when citing this.")
    print("=" * 78)


if __name__ == "__main__":
    print_summary()
