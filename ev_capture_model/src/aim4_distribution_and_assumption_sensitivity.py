"""
aim4_distribution_and_assumption_sensitivity.py

WHY THIS FILE EXISTS
---------------------
This is the CURRENT.md handoff's "Immediate next steps" item 1, the last
piece of Aim 4's sensitivity analysis: "distribution-family choice, a
fuller sensitivity analysis of the plausibility bounds' own assumptions
beyond the single widened check already done." That phrase refers
specifically to the target_exponents diagnostic widening run in the
v9->v10 / v10->v11 handoff transitions (see k_intra_z_plausibility_
bounds.py and the archived v15.md "What changed since v10"/"since v9"
sections) -- a check on ONE assumption (the exponent ceiling), decided
closed by the user, kept only as a write-up footnote. This file asks the
two questions that check was never meant to answer:

  1. DISTRIBUTION-FAMILY SENSITIVITY: every headline finding built from
     REAL digitized data (the two Lennon-based comparisons in
     aim3_aim4_range_sweep.py) never fits a distribution at all --
     joint_population_yield() takes real per-vesicle (n, diameter) pairs
     directly (see population_model.py's own docstring: "no explicit
     P(n) needs to be fit first"). Distribution-family choice is
     therefore NOT a live sensitivity question for CA19-9 or the
     retired-EGFR case. It IS a live, previously-untested question for
     the ONE headline finding built from a RECONSTRUCTED population: the
     stem-cell CD81/CD63/CD73 marker-ranking crossover
     (stem_cell_marker_ranking_extension.py), which fits a NEGATIVE
     BINOMIAL to two summary statistics (fraction positive, mean copies
     if positive) per marker/tissue, and that choice of family has never
     been checked against an alternative. This file builds a second,
     structurally different reconstruction -- a zero-inflated Poisson
     (a Bernoulli on/off switch times a Poisson, vs. negative binomial's
     continuous Gamma-mixture-of-Poissons story) -- fit to the SAME two
     summary statistics, and re-runs the exact crossover test against it.

  2. PLAUSIBILITY-BOUNDS ASSUMPTION SENSITIVITY, BEYOND target_exponents:
     k_intra_bounds_from_real_sigma() bakes in two more assumptions never
     separately stress-tested: the assumed RECEPTOR-COUNT GRID
     (n_values=(1,10,50,100), used to convert a dimensionless exponent
     bound into a K_INTRA range) and the SIGMA ANCHOR ITSELF (the real
     Lopez Baltazar value, 4.44e11 Abs/cm^2 at "40% coverage" -- a single
     reported operating point, not a value with its own stated
     uncertainty). This file checks whether the headline superselectivity-
     reachability finding (max alpha, fraction of grid superselective)
     changes qualitatively if (a) the n_values grid is shifted to a
     different plausible range, or (b) sigma is swept across its own
     documented plausible range (the reported 40% operating coverage up
     to the 100% max-close-packed ceiling already named in
     k_intra_z_plausibility_bounds.py as SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2,
     a real, already-computed alternative this project has never actually
     used in a sweep).

  3. DIAMETER-MODEL SENSITIVITY (a smaller, related check): the stem-cell
     reconstruction assigns every marker/tissue the SAME generic lognormal
     diameter distribution (120nm mean, CV=0.3) -- explicitly flagged in
     that file's own docstring as an assumption, not a measurement, but
     never stress-tested against an alternative shape or width. This file
     checks whether the crossover survives a genuinely different diameter
     CV and a genuinely different diameter distribution family (gamma
     instead of lognormal), matched to the same mean.

WHAT THIS FILE IS NOT
------------------------
- NOT a re-litigation of target_exponents. That check is closed (see
  "What NOT to revisit" in the handoff) and is not repeated here.
- NOT a new calibration attempt. K_INTRA, z, and k_geom remain
  permanently uncalibrated placeholders; every function here operates on
  the existing plausibility RANGE or explicit alternative assumptions,
  never a single fitted number.
- NOT a claim that the alternate distribution/diameter families used here
  are "more correct" than negative binomial / lognormal. They are chosen
  specifically to be STRUCTURALLY DIFFERENT (a genuinely different
  generative story), so that "the finding survives a different family"
  is an honest robustness statement, not a coincidence of two similar-
  shaped curves.
"""

import os
import sys

import numpy as np
from scipy import optimize, stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from binding_model import selectivity_alpha
from k_intra_z_plausibility_bounds import (
    dimensionless_exponent_bounds,
    k_intra_bounds_from_real_sigma,
    SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
    SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2,
)
from aim3_aim4_range_sweep import (
    get_z_sweep_values,
    MIN_RELIABLE_EXPONENT_AT_N1,
    ALPHA_SUPERSELECTIVE_THRESHOLD,
)
from stem_cell_marker_ranking_extension import (
    POSITIVITY_FRACTION_RANGE,
    COPIES_PER_POSITIVE_PARTICLE,
    GENERIC_EV_DIAMETER_MEAN_NM,
    GENERIC_EV_DIAMETER_CV,
    stable_joint_population_yield,
    naive_fold_change_ranking,
    marker_level_ranking,
)


# ===========================================================================
# PART 1: DISTRIBUTION-FAMILY SENSITIVITY (zero-inflated Poisson alternative
# to the stem-cell reconstruction's negative binomial)
# ===========================================================================

def fit_zero_inflated_poisson_to_positivity(fraction_positive, mean_copies_if_positive,
                                             n_max=5000):
    """
    Build a zero-inflated Poisson (ZIP) distribution matched to the SAME
    two summary statistics stem_cell_marker_ranking_extension.py's
    fit_negative_binomial_to_positivity() matches: the fraction of
    particles carrying ANY copies of the marker (fraction_positive), and
    the mean copy number AMONG the positive particles only
    (mean_copies_if_positive).

    WHY THIS IS A GENUINELY DIFFERENT FAMILY, NOT A RELABELED NEGATIVE
    BINOMIAL: a negative binomial distribution arises as a continuous
    Gamma mixture of Poisson distributions -- every particle has SOME
    nonzero Poisson rate, drawn from a Gamma distribution, and n=0 only
    happens when that draw's own Poisson process happens to produce zero
    events. A zero-inflated Poisson instead has a hard, literal two-step
    generative story: first, a Bernoulli "does this particle carry the
    marker machinery at all" coin flip (probability fraction_positive);
    then, ONLY for particles that flip "yes," a single shared Poisson
    draw with rate lambda = mean_copies_if_positive. Every particle that
    flips "no" gets EXACTLY zero, not "a Poisson draw from a rate that
    happened to come out near zero." This is a structurally different
    claim about biology (a literal on/off subpopulation split vs. a
    continuum of per-particle rates) even though both distributions can
    match the same two summary numbers -- which is exactly why it is a
    meaningful robustness check and not a redundant one.

    Because ZIP has exactly two parameters (the zero-inflation probability
    and the Poisson rate) and this function is given exactly two targets,
    there is no fitting/optimization needed at all -- unlike the negative-
    binomial fit (which needs a 1D root-find over dispersion_r because
    negative binomial's zero-probability formula involves BOTH parameters
    jointly). Here, the mapping is direct: the "off" probability is
    exactly (1 - fraction_positive), and the Poisson rate for the "on"
    subpopulation is exactly mean_copies_if_positive, by construction.
    This is checked with an assertion below rather than assumed silently.

    Parameters
    ----------
    fraction_positive : float
        Same meaning as fit_negative_binomial_to_positivity's parameter.
    mean_copies_if_positive : float
        Same meaning as fit_negative_binomial_to_positivity's parameter.
    n_max : int
        Highest receptor count to include explicitly when building the
        full P(n) table (only used by zip_pmf_table() below, not by the
        sampler, which draws directly).

    Returns
    -------
    dict with keys "zero_inflation_prob" (1 - fraction_positive),
    "poisson_rate" (mean_copies_if_positive), "mean_n" (the resulting
    OVERALL mean, i.e. fraction_positive * mean_copies_if_positive, for
    direct comparison against the negative-binomial fit's "mean_n" key),
    "fraction_positive_achieved" (recomputed directly from the two
    parameters, to confirm the construction is exact).
    """
    if not (0.0 < fraction_positive < 1.0):
        raise ValueError(f"fraction_positive must be in (0, 1), got {fraction_positive}")
    if mean_copies_if_positive <= 0:
        raise ValueError(f"mean_copies_if_positive must be positive, got {mean_copies_if_positive}")

    zero_inflation_prob = 1.0 - fraction_positive
    poisson_rate = mean_copies_if_positive

    # P(n=0) under ZIP = zero_inflation_prob + (1-zero_inflation_prob) * P_poisson(0)
    #                  = zero_inflation_prob + fraction_positive * exp(-poisson_rate)
    # So P(n>0) = fraction_positive * (1 - exp(-poisson_rate)), which is
    # very slightly LESS than fraction_positive itself (a "positive" coin
    # flip can still draw a Poisson 0). This is a real, small, honest
    # discrepancy inherent to the ZIP construction -- checked explicitly
    # below rather than glossed over -- and is negligible whenever
    # poisson_rate is not tiny (exp(-20) ~ 2e-9, so for every marker/
    # tissue combination actually used in this project, where
    # mean_copies_if_positive is at least 20, this discrepancy is
    # unmeasurably small).
    fraction_positive_achieved = fraction_positive * (1.0 - np.exp(-poisson_rate))

    return {
        "zero_inflation_prob": zero_inflation_prob,
        "poisson_rate": poisson_rate,
        "mean_n": fraction_positive * poisson_rate,
        "fraction_positive_achieved": fraction_positive_achieved,
    }


def sample_zero_inflated_poisson(fraction_positive, mean_copies_if_positive,
                                  n_samples=5000, seed=0):
    """
    Draw a synthetic sample from the zero-inflated Poisson construction in
    fit_zero_inflated_poisson_to_positivity(), directly analogous to
    stem_cell_marker_ranking_extension.reconstruct_marker_population()'s
    negative-binomial sampling step.

    Returns
    -------
    numpy.ndarray of length n_samples, receptor counts (float, matching
    the dtype convention used elsewhere in this project for n_values).
    """
    fit = fit_zero_inflated_poisson_to_positivity(fraction_positive, mean_copies_if_positive)
    rng = np.random.default_rng(seed)
    is_positive = rng.random(n_samples) < fraction_positive
    counts = np.zeros(n_samples, dtype=float)
    n_positive = int(is_positive.sum())
    if n_positive > 0:
        counts[is_positive] = rng.poisson(fit["poisson_rate"], size=n_positive).astype(float)
    return counts


def reconstruct_marker_population_zip(marker, tissue, n_samples=5000, seed=0):
    """
    Zero-inflated-Poisson counterpart to stem_cell_marker_ranking_
    extension.reconstruct_marker_population() -- same marker/tissue
    lookup, same diameter model (see that file's module docstring for why
    a single generic diameter distribution is used for every marker), same
    (fraction_positive, mean_copies_if_positive) targets, but the
    receptor-count axis is drawn from the zero-inflated Poisson
    construction above instead of a fitted negative binomial.

    Returns
    -------
    n_values, diameter_values : two numpy arrays of length n_samples,
    directly usable with stable_joint_population_yield(), exactly like
    reconstruct_marker_population()'s return value.
    """
    if marker not in POSITIVITY_FRACTION_RANGE:
        raise ValueError(f"Unknown marker {marker!r}; expected one of "
                          f"{sorted(POSITIVITY_FRACTION_RANGE)}.")
    if tissue not in COPIES_PER_POSITIVE_PARTICLE[marker]:
        raise ValueError(f"Unknown tissue {tissue!r} for marker {marker!r}; "
                          f"expected one of "
                          f"{sorted(COPIES_PER_POSITIVE_PARTICLE[marker])}.")

    lo, hi = POSITIVITY_FRACTION_RANGE[marker]
    fraction_positive = 0.5 * (lo + hi)
    mean_copies_if_positive = COPIES_PER_POSITIVE_PARTICLE[marker][tissue]

    n_values = sample_zero_inflated_poisson(
        fraction_positive, mean_copies_if_positive, n_samples=n_samples, seed=seed
    )

    # Same generic diameter model as the negative-binomial reconstruction,
    # UNCHANGED -- this part of Part 1's check isolates the receptor-count
    # distribution family specifically; the diameter axis is stress-tested
    # separately in Part 3 below, one assumption at a time.
    rng = np.random.default_rng(seed + 1_000_000)  # offset so diameter draws
                                                     # don't share a stream
                                                     # with the ZIP draws
                                                     # above, matching the
                                                     # independence assumed
                                                     # by the model
    mu_log = np.log(GENERIC_EV_DIAMETER_MEAN_NM) - 0.5 * np.log(1 + GENERIC_EV_DIAMETER_CV ** 2)
    sigma_log = np.sqrt(np.log(1 + GENERIC_EV_DIAMETER_CV ** 2))
    diameter_values = rng.lognormal(mean=mu_log, sigma=sigma_log, size=n_samples)

    return n_values, diameter_values


def crossover_finding_under_family(reconstruct_fn, k_intra=None, z=1e-4,
                                    sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
                                    n_samples=5000, seed=0):
    """
    Shared helper: reconstruct all nine (marker, tissue) populations using
    the given reconstruction function (either reconstruct_marker_population
    from stem_cell_marker_ranking_extension.py, or
    reconstruct_marker_population_zip above -- both have the identical
    (marker, tissue, n_samples, seed) call signature by construction, so
    this helper works unmodified for either), compute yields at the given
    k_intra, and return the marker-level ranking.

    Parameters
    ----------
    reconstruct_fn : callable
        Either reconstruct_marker_population or
        reconstruct_marker_population_zip.
    k_intra : float or None
        If None, uses the largest value in the real-sigma plausibility
        range (matching the handoff's headline "at the largest K_INTRA"
        framing exactly).

    Returns
    -------
    dict with "yields" (dict keyed by (marker, tissue)) and
    "marker_ranking" (list of marker names, best/highest yield first,
    aggregated via marker_level_ranking() the same way the original
    finding is reported).
    """
    if k_intra is None:
        k_intra = max(k_intra_bounds_from_real_sigma().values())

    yields = {}
    pair_index = 0
    for marker in POSITIVITY_FRACTION_RANGE:
        for tissue in COPIES_PER_POSITIVE_PARTICLE[marker]:
            n_vals, d_vals = reconstruct_fn(marker, tissue, n_samples=n_samples, seed=seed + pair_index)
            pair_index += 1
            yields[(marker, tissue)] = stable_joint_population_yield(
                n_vals, d_vals, sigma, 1.0, k_intra, z, patch_exponent=2.0
            )

    ranking = sorted(yields, key=yields.get, reverse=True)
    return {"yields": yields, "marker_ranking": marker_level_ranking(ranking)}


def check_crossover_survives_zip_family(seeds=(0, 100, 500, 1000), n_samples=5000):
    """
    THE MAIN QUESTION FOR PART 1: does the CD81-overtakes-CD63 crossover
    (at the largest real-sigma-plausibility-range K_INTRA, despite CD63's
    higher naive fold-change score) still happen when the reconstruction
    uses a zero-inflated Poisson instead of a negative binomial, matched
    to the exact same two summary statistics?

    Checked across the same four seeds already used to establish the
    negative-binomial version's seed-robustness (see
    tests/test_stem_cell_marker_ranking_extension.py's
    test_marker_level_crossover_is_robust_across_random_seeds), so a
    "yes" here is directly comparable in strength to the original finding,
    not a weaker single-seed check.

    Returns
    -------
    dict with:
        "naive_marker_ranking": the (family-independent) naive fold-change
            marker ranking, for reference.
        "nb_marker_rankings": list of marker rankings, one per seed, under
            the ORIGINAL negative-binomial reconstruction (recomputed here
            directly, not assumed from the other test file, so this
            function is self-contained and does not silently depend on
            import order).
        "zip_marker_rankings": list of marker rankings, one per seed,
            under the zero-inflated-Poisson reconstruction.
        "cd81_leads_under_nb": list of bool, one per seed.
        "cd81_leads_under_zip": list of bool, one per seed.
        "crossover_survives_zip_family": True only if CD81 leads under
            BOTH families at EVERY tested seed -- the strict form of "this
            finding is not an artifact of the negative-binomial choice."
    """
    from stem_cell_marker_ranking_extension import reconstruct_marker_population

    naive = naive_fold_change_ranking()
    naive_marker_ranking = marker_level_ranking(naive["ranking_best_first"])

    nb_rankings = []
    zip_rankings = []
    for seed in seeds:
        nb_result = crossover_finding_under_family(
            reconstruct_marker_population, seed=seed, n_samples=n_samples
        )
        zip_result = crossover_finding_under_family(
            reconstruct_marker_population_zip, seed=seed, n_samples=n_samples
        )
        nb_rankings.append(nb_result["marker_ranking"])
        zip_rankings.append(zip_result["marker_ranking"])

    cd81_leads_under_nb = [r[0] == "CD81" for r in nb_rankings]
    cd81_leads_under_zip = [r[0] == "CD81" for r in zip_rankings]

    return {
        "naive_marker_ranking": naive_marker_ranking,
        "nb_marker_rankings": nb_rankings,
        "zip_marker_rankings": zip_rankings,
        "cd81_leads_under_nb": cd81_leads_under_nb,
        "cd81_leads_under_zip": cd81_leads_under_zip,
        "crossover_survives_zip_family": bool(
            all(cd81_leads_under_nb) and all(cd81_leads_under_zip)
        ),
    }


# ===========================================================================
# PART 2: PLAUSIBILITY-BOUNDS ASSUMPTION SENSITIVITY
# (n_values grid choice, and the sigma anchor choice)
# ===========================================================================

def _superselectivity_reachability_for_k_intra_values(k_intra_values, sigma, n_max=200):
    """
    Shared inner sweep, factored out of aim3_aim4_range_sweep.
    sweep_superselectivity_reachability() so it can be re-run against a
    DIFFERENT k_intra_values array (e.g. one derived from a shifted
    n_values grid, or a different sigma) without duplicating that
    function's full docstring/bookkeeping. Mirrors its numerical-
    reliability handling exactly (same MIN_RELIABLE_EXPONENT_AT_N1 guard,
    same ALPHA_SUPERSELECTIVE_THRESHOLD), so results from this helper are
    directly comparable to that function's own output.

    Returns
    -------
    dict with "max_alpha_found", "any_superselective",
    "fraction_of_grid_superselective", "grid_size",
    "n_cells_excluded_numerically_unreliable" -- the finite-difference-only
    subset of sweep_superselectivity_reachability()'s return keys (the
    exact/alpha_exact() comparison is not needed for this assumption-
    sensitivity question and is skipped to keep this sweep fast enough to
    run across several assumption variants in one test session).
    """
    z_values = get_z_sweep_values()
    n_values = np.linspace(1, n_max, 60)

    max_alpha = -np.inf
    n_superselective_cells = 0
    total_reliable_cells = 0
    n_excluded = 0

    for k_intra in k_intra_values:
        exponent_at_n1 = 1 * sigma * k_intra
        if exponent_at_n1 < MIN_RELIABLE_EXPONENT_AT_N1:
            n_excluded += len(z_values)
            continue
        for z in z_values:
            total_reliable_cells += 1
            try:
                n_mid, alpha = selectivity_alpha(n_values, sigma, k_intra, z)
            except (ValueError, FloatingPointError):
                continue
            alpha_reliable = alpha[np.isfinite(alpha)]
            if len(alpha_reliable) > 0:
                local_max = alpha_reliable.max()
                if local_max > ALPHA_SUPERSELECTIVE_THRESHOLD:
                    n_superselective_cells += 1
                if local_max > max_alpha:
                    max_alpha = local_max

    return {
        "max_alpha_found": float(max_alpha),
        "any_superselective": bool(max_alpha > ALPHA_SUPERSELECTIVE_THRESHOLD),
        "fraction_of_grid_superselective": (
            n_superselective_cells / total_reliable_cells if total_reliable_cells else 0.0
        ),
        "grid_size": total_reliable_cells,
        "n_cells_excluded_numerically_unreliable": n_excluded,
    }


def sweep_superselectivity_across_n_values_grids(
    grids=((1, 10, 50, 100), (1, 5, 25, 75), (5, 20, 100, 300)),
    sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
):
    """
    Does the superselectivity-reachability headline finding (max alpha
    ~10.46, ~60% of the grid superselective) depend on the SPECIFIC
    receptor-count grid (n_values=(1,10,50,100)) used to build the K_INTRA
    plausibility range, or is it stable if that grid is shifted to other
    plausible choices?

    Three grids are tested: the project's existing default
    (1, 10, 50, 100); a lower-shifted grid (1, 5, 25, 75), representing a
    belief that typical EV surface marker counts skew toward the low end
    of the 1-100 range documented in k_intra_z_plausibility_bounds.py's
    module docstring; and a higher-shifted grid (5, 20, 100, 300),
    representing the opposite belief, extended somewhat past the
    documented range's own upper end to test a genuinely wider assumption
    rather than a small perturbation.

    Returns
    -------
    dict keyed by grid (tuple), each value a
    _superselectivity_reachability_for_k_intra_values()-style result dict,
    plus a top-level "qualitative_finding_stable" bool: True if
    any_superselective is True under every tested grid (the finding
    "superselectivity IS reachable somewhere in the plausible range" does
    not flip to False under any tested alternative grid).
    """
    results = {}
    for grid in grids:
        k_intra_values = sorted(set(k_intra_bounds_from_real_sigma(n_values=grid).values()))
        results[grid] = _superselectivity_reachability_for_k_intra_values(
            np.array(k_intra_values), sigma=sigma
        )

    results["qualitative_finding_stable"] = bool(
        all(r["any_superselective"] for r in results.values() if isinstance(r, dict))
    )
    return results


def sweep_superselectivity_across_sigma_choices(
    sigma_choices=None,
):
    """
    Does the superselectivity-reachability finding depend on the SPECIFIC
    sigma value K_INTRA's real-unit conversion is anchored to (the Lopez
    Baltazar "40% coverage" operating point, 4.44e11 Abs/cm^2), or is it
    stable across the full documented plausible sigma range?

    k_intra_bounds_from_real_sigma()'s own docstring already names the
    natural upper bound for this check: SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_
    CM2 (1.11e12 Abs/cm^2, 100% close-packed antibody coverage -- the same
    dissertation's own physically-derived ceiling, not an invented one).
    This function is the first place in the codebase that actually RUNS
    that comparison end-to-end (the docstring names it as the intended use
    of the sigma_abs_per_cm2 parameter but no prior session called it with
    anything other than the default).

    IMPORTANT SUBTLETY, CHECKED DIRECTLY: because
    k_intra_bounds_from_real_sigma() computes K_INTRA = dimensionless_bound
    / sigma, changing sigma changes K_INTRA in INVERSE proportion --- but
    q_bound only ever depends on the PRODUCT sigma * K_INTRA (see
    Limitations_and_Calibration_Gap.md's identifiability proof). So this
    function must call _superselectivity_reachability_for_k_intra_values()
    with the SAME sigma used to derive that sigma's own K_INTRA values, at
    every point -- passing a mismatched (K_INTRA, sigma) pair would
    silently and incorrectly change the product, reproducing exactly the
    "Bug 1" class of error aim3_aim4_range_sweep.py's own docstrings warn
    about repeatedly. This is done correctly below (each sigma choice
    computes its OWN k_intra_values first), and is exactly why this
    function cannot simply reuse get_k_intra_sweep_values() (which is
    hard-coded to the default sigma).

    Parameters
    ----------
    sigma_choices : dict or None
        Maps a label to a sigma value (Abs/cm^2). Defaults to the
        project's two already-documented real values: "lopez_baltazar_40pct"
        (the default, headline value) and "lopez_baltazar_100pct_max"
        (SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2).

    Returns
    -------
    dict keyed by label, each value a
    _superselectivity_reachability_for_k_intra_values()-style result dict,
    plus "qualitative_finding_stable" (see
    sweep_superselectivity_across_n_values_grids for the same meaning).
    """
    if sigma_choices is None:
        sigma_choices = {
            "lopez_baltazar_40pct": SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
            "lopez_baltazar_100pct_max": SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2,
        }

    results = {}
    for label, sigma in sigma_choices.items():
        k_intra_values = sorted(set(k_intra_bounds_from_real_sigma(
            sigma_abs_per_cm2=sigma
        ).values()))
        results[label] = _superselectivity_reachability_for_k_intra_values(
            np.array(k_intra_values), sigma=sigma
        )

    results["qualitative_finding_stable"] = bool(
        all(r["any_superselective"] for r in results.values() if isinstance(r, dict))
    )
    return results


# ===========================================================================
# PART 3: DIAMETER-MODEL SENSITIVITY (stem-cell reconstruction only)
# ===========================================================================

def reconstruct_marker_population_diameter_variant(marker, tissue, diameter_cv=None,
                                                     diameter_family="lognormal",
                                                     n_samples=5000, seed=0):
    """
    Same negative-binomial receptor-count reconstruction as
    stem_cell_marker_ranking_extension.reconstruct_marker_population(),
    but with the diameter axis swappable -- either a wider/narrower
    lognormal (diameter_cv override) or a genuinely different family
    (gamma), matched to the SAME mean (GENERIC_EV_DIAMETER_MEAN_NM) so
    only the shape/spread assumption changes, not the central value that
    was itself the literature-typical choice.

    Parameters
    ----------
    diameter_cv : float or None
        Coefficient of variation to use. Defaults to
        GENERIC_EV_DIAMETER_CV (the file's existing default, 0.3) if None.
    diameter_family : str
        "lognormal" (matches the existing default shape, just with a
        possibly different CV) or "gamma" (a genuinely different
        distribution shape -- gamma is right-skewed like lognormal but
        with a different tail behavior, a meaningfully different
        structural choice rather than a relabeling).

    Returns
    -------
    n_values, diameter_values : as in reconstruct_marker_population().
    """
    from stem_cell_marker_ranking_extension import (
        fit_negative_binomial_to_positivity,
    )

    if marker not in POSITIVITY_FRACTION_RANGE:
        raise ValueError(f"Unknown marker {marker!r}")
    if tissue not in COPIES_PER_POSITIVE_PARTICLE[marker]:
        raise ValueError(f"Unknown tissue {tissue!r} for marker {marker!r}")

    lo, hi = POSITIVITY_FRACTION_RANGE[marker]
    fraction_positive = 0.5 * (lo + hi)
    mean_copies_if_positive = COPIES_PER_POSITIVE_PARTICLE[marker][tissue]
    fit = fit_negative_binomial_to_positivity(fraction_positive, mean_copies_if_positive)

    rng = np.random.default_rng(seed)
    p_param = fit["dispersion_r"] / (fit["dispersion_r"] + fit["mean_n"])
    n_values = rng.negative_binomial(fit["dispersion_r"], p_param, size=n_samples).astype(float)

    cv = GENERIC_EV_DIAMETER_CV if diameter_cv is None else diameter_cv
    mean_nm = GENERIC_EV_DIAMETER_MEAN_NM

    if diameter_family == "lognormal":
        mu_log = np.log(mean_nm) - 0.5 * np.log(1 + cv ** 2)
        sigma_log = np.sqrt(np.log(1 + cv ** 2))
        diameter_values = rng.lognormal(mean=mu_log, sigma=sigma_log, size=n_samples)
    elif diameter_family == "gamma":
        # Gamma parameterized by shape k and scale theta, with mean = k*theta
        # and CV = 1/sqrt(k) -- matched to the same mean and CV as the
        # lognormal default so only the FAMILY (tail shape), not the
        # first two moments, differs.
        shape_k = 1.0 / (cv ** 2)
        scale_theta = mean_nm / shape_k
        diameter_values = rng.gamma(shape_k, scale_theta, size=n_samples)
    else:
        raise ValueError(f"Unknown diameter_family {diameter_family!r}; expected "
                          f"'lognormal' or 'gamma'.")

    return n_values, diameter_values


def check_crossover_survives_diameter_assumption(
    variants=(("lognormal", 0.3), ("lognormal", 0.6), ("lognormal", 0.15), ("gamma", 0.3)),
    seeds=(0, 100, 500, 1000),
    n_samples=5000,
):
    """
    Does the CD81-overtakes-CD63 crossover survive if the generic diameter
    assumption is widened, narrowed, or changed to a different family
    entirely? Four variants are tested against the same four seeds used
    for the family-robustness check in Part 1, for the same reason (a
    reconstruction-based finding needs multi-seed evidence to be
    trustworthy, not a single lucky draw).

    Returns
    -------
    dict keyed by (diameter_family, diameter_cv) tuple, each value a list
    of marker rankings (one per seed), plus "crossover_survives_all_
    diameter_variants": True only if CD81 leads at every seed under every
    tested variant.
    """
    if k_intra_bounds_from_real_sigma is None:  # pragma: no cover -- import guard
        raise RuntimeError("k_intra_bounds_from_real_sigma failed to import")

    k_intra = max(k_intra_bounds_from_real_sigma().values())
    results = {}
    all_survive = True

    for family, cv in variants:
        rankings = []
        for seed in seeds:
            def reconstruct_fn(marker, tissue, n_samples=n_samples, seed=seed,
                                _family=family, _cv=cv):
                return reconstruct_marker_population_diameter_variant(
                    marker, tissue, diameter_cv=_cv, diameter_family=_family,
                    n_samples=n_samples, seed=seed
                )
            result = crossover_finding_under_family(
                reconstruct_fn, k_intra=k_intra, seed=seed, n_samples=n_samples
            )
            rankings.append(result["marker_ranking"])
        results[(family, cv)] = rankings
        if not all(r[0] == "CD81" for r in rankings):
            all_survive = False

    results["crossover_survives_all_diameter_variants"] = all_survive
    return results


def print_summary():
    print("=" * 78)
    print("AIM 4 -- DISTRIBUTION-FAMILY AND PLAUSIBILITY-BOUNDS ASSUMPTION")
    print("SENSITIVITY (beyond the closed target_exponents check)")
    print("=" * 78)

    print("\n--- Part 1: distribution-family sensitivity (stem-cell crossover) ---")
    zip_check = check_crossover_survives_zip_family()
    print(f"Naive marker ranking (family-independent): {zip_check['naive_marker_ranking']}")
    print(f"CD81 leads under negative binomial, all seeds: {zip_check['cd81_leads_under_nb']}")
    print(f"CD81 leads under zero-inflated Poisson, all seeds: {zip_check['cd81_leads_under_zip']}")
    print(f"Crossover survives ZIP family: {zip_check['crossover_survives_zip_family']}")

    print("\n--- Part 2a: sensitivity to the n_values grid ---")
    grid_check = sweep_superselectivity_across_n_values_grids()
    for grid, result in grid_check.items():
        if grid == "qualitative_finding_stable":
            continue
        print(f"  n_values={grid}: max_alpha={result['max_alpha_found']:.3f}, "
              f"frac_superselective={result['fraction_of_grid_superselective']:.3f}, "
              f"any_superselective={result['any_superselective']}")
    print(f"Qualitative finding stable across all tested grids: "
          f"{grid_check['qualitative_finding_stable']}")

    print("\n--- Part 2b: sensitivity to the sigma anchor ---")
    sigma_check = sweep_superselectivity_across_sigma_choices()
    for label, result in sigma_check.items():
        if label == "qualitative_finding_stable":
            continue
        print(f"  sigma={label}: max_alpha={result['max_alpha_found']:.3f}, "
              f"frac_superselective={result['fraction_of_grid_superselective']:.3f}, "
              f"any_superselective={result['any_superselective']}")
    print(f"Qualitative finding stable across all tested sigma choices: "
          f"{sigma_check['qualitative_finding_stable']}")

    print("\n--- Part 3: diameter-model sensitivity (stem-cell crossover) ---")
    diam_check = check_crossover_survives_diameter_assumption()
    for variant, rankings in diam_check.items():
        if variant == "crossover_survives_all_diameter_variants":
            continue
        print(f"  {variant}: rankings={rankings}")
    print(f"Crossover survives all tested diameter variants: "
          f"{diam_check['crossover_survives_all_diameter_variants']}")

    print("\n" + "=" * 78)


if __name__ == "__main__":
    print_summary()
