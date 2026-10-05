"""
population_model.py

Aim 1's core novelty: instead of asking "what fraction of vesicles WITH
EXACTLY n RECEPTORS get captured" (that's binding_model.capture_probability,
a per-vesicle question), we ask "what fraction of a REAL POPULATION of
vesicles -- which all have DIFFERENT receptor counts -- gets captured."

WHY THIS IS THE IMPORTANT STEP
---------------------------------
A real blood sample doesn't contain vesicles that all have exactly n=20
receptors. It contains a whole population: some vesicles with 0 receptors,
some with 5, some with 40, etc. The published literature almost always
reports just the MEAN receptor count across this population (a single
number). The standard practice in the field ("fold-change ranking") assumes
that plugging the mean into a capture model is a good enough substitute for
knowing the whole distribution.

Because our capture_probability() formula is NONLINEAR in n (n sits inside
an exponential -- see binding_model.py), that substitution is NOT generally
valid. This is a basic but important fact from probability: for a nonlinear
function g(n),

    average of g(n)   is NOT THE SAME AS   g(average of n)

in general. (This inequality even has a name for a wide class of functions:
Jensen's inequality.) The only case where the two ARE guaranteed to match
exactly is a special one we'll test for explicitly below (the Poisson case).

WHAT THIS FILE COMPUTES
--------------------------
Given:
    - a distribution over receptor counts, P(n) -- represented here simply
      as two paired arrays: n_values (the possible receptor counts) and
      probabilities (how likely each count is, summing to 1)
    - a ligand density sigma, k_intra, and z (same meaning as before)

we compute the POPULATION-AVERAGE capture probability (also called "yield"
in the handoff):

    Yield(sigma) = sum over n of [ P(n) * theta(n, sigma) ]

This is just a weighted average: for each possible receptor count n, look
up how likely a random vesicle is to have that many receptors (P(n)), and
how likely it is to be captured if it does (theta(n, sigma)), multiply
those together, and add up over all possible n. That's it -- there's no
new physics here beyond binding_model.py, just correctly averaging over
the real distribution instead of using a single "typical" n.

REPRESENTING P(n)
-------------------
We build P(n) using scipy.stats distributions where possible (e.g. Poisson,
negative binomial for "zero-inflated and overdispersed" -- see handoff),
which is more robust and less error-prone than manually computing
distribution formulas by hand.
"""

import numpy as np
from scipy import stats

from binding_model import capture_probability, q_bound


def population_yield(n_values, probabilities, sigma, k_intra, z):
    """
    Compute the population-average capture probability (yield), given an
    explicit discrete distribution over receptor counts.

    Parameters
    ----------
    n_values : array of int
        The possible receptor counts (e.g. [0, 1, 2, ..., 200]).
    probabilities : array of float
        P(n) for each value in n_values. Must be the same length as
        n_values and should sum to (approximately) 1. (We do NOT
        silently renormalize -- if these don't sum close to 1, something
        upstream is probably wrong, and we want that to be visible rather
        than hidden.)
    sigma : float
        Ligand surface density.
    k_intra : float
        Binding-strength parameter.
    z : float
        Activity/prefactor.

    Returns
    -------
    float
        The population-average capture probability (yield), between 0
        and 1.
    """
    n_values = np.asarray(n_values, dtype=float)
    probabilities = np.asarray(probabilities, dtype=float)

    total_prob = probabilities.sum()
    if not np.isclose(total_prob, 1.0, atol=1e-3):
        raise ValueError(
            f"probabilities sum to {total_prob:.6f}, not 1.0 -- check the "
            "distribution before using it (this is not auto-normalized "
            "on purpose, so mistakes don't get silently hidden)."
        )

    theta_per_n = capture_probability(n_values, sigma, k_intra, z)

    # Vesicles with n=0 receptors can never be captured (there is nothing to
    # form a bond with). capture_probability(0, ...) should already
    # naturally come out to 0 -- q_bound(0, ...) = exp(0) - 1 = 0, so
    # theta = 0. We don't need a special case, but it's worth stating
    # explicitly since it's an important physical checkpoint: the model
    # should never predict that receptor-less vesicles get captured.

    yield_value = np.sum(probabilities * theta_per_n)
    return yield_value


def mean_field_qb(mean_n, sigma, k_intra):
    """
    Closed-form expression for the POPULATION-AVERAGE bond-count partition
    function, E[q_b(n)], when receptor counts n are Poisson-distributed
    with mean mean_n.

    IMPORTANT SCOPE NOTE -- read this before using this function for a
    "mean-field is exact" claim elsewhere in the project.

    While working out the Poisson mean-field theorem for this project I
    found that "mean-field theory is exact for Poisson receptors" is TRUE,
    but only at THIS specific level -- the average of q_b(n) over a Poisson
    population of receptor counts exactly equals this closed-form
    expression, with no approximation:

        E[q_b(n)] = exp( mean_n * (exp(sigma * k_intra) - 1) ) - 1

    This can be checked directly: for a Poisson random variable n with
    mean mean_n, its moment generating function is
        E[exp(n * t)] = exp( mean_n * (exp(t) - 1) )
    Setting t = sigma * k_intra and subtracting 1 from both sides gives
    exactly the formula above, since q_b(n) = exp(n * sigma * k_intra) - 1
    by definition (see binding_model.q_bound). This is genuinely exact --
    not a numerical coincidence or an approximation that happens to be
    close -- and is confirmed numerically in test_population_model.py to
    machine precision (differences around 1e-10, i.e. floating-point
    rounding only).

    WHAT IS *NOT* EXACT: converting q_b into an actual capture probability
    requires the additional Langmuir step theta = z*q_b / (1 + z*q_b) (see
    binding_model.capture_probability). That division is ALSO nonlinear,
    and averaging theta(n) over the Poisson population does NOT exactly
    equal evaluating theta at this mean-field q_b -- the two differ by
    roughly 1-3% for typical parameter ranges tested here, more at
    extreme z. This is a second, separate source of approximation error
    stacked on top of the first (exact) step. In other words:

        E[q_b(n)]                                  <- EXACT closed form (this function)
        theta_mean_field := z*E[q_b(n)] / (1+z*E[q_b(n)])   <- one more step, APPROXIMATE
        E[theta(n)]  (true population yield)        <- what population_yield() computes directly

    We keep these as two separate, clearly-named functions
    (mean_field_qb here, and population_yield elsewhere) specifically so
    the code cannot accidentally conflate "exact" with "approximate."

    Parameters
    ----------
    mean_n : float
        Mean receptor count of the Poisson-distributed population.
    sigma : float
        Ligand surface density.
    k_intra : float
        Binding-strength parameter.

    Returns
    -------
    float
        The exact population-average bond-count partition function,
        E[q_b(n)].
    """
    return np.exp(mean_n * (np.exp(sigma * k_intra) - 1.0)) - 1.0


def poisson_distribution(mean_n, n_max):
    """
    Build a Poisson distribution over receptor counts 0, 1, ..., n_max,
    with the given mean.

    A Poisson distribution describes counts that arise from many
    independent, low-probability opportunities for an event to happen
    (here: many independent chances for a receptor protein to end up on
    the vesicle membrane during its formation) -- it's the "textbook
    default" model for count data when there's no reason to expect
    clustering or extra variability. We use it here specifically because
    it is the SPECIAL CASE where population averaging and using the mean
    directly are proven to give identical answers -- see
    poisson_mean_field_theorem_holds() below. Real EV data is NOT well
    described by a Poisson distribution (see negative_binomial_distribution
    below), but Poisson is exactly the right choice for a unit test.

    Parameters
    ----------
    mean_n : float
        The desired mean receptor count.
    n_max : int
        Highest receptor count to include explicitly. Should be chosen
        large enough that the Poisson distribution's probability mass
        above n_max is negligible (a good rule of thumb: n_max should be
        at least mean_n + 10 * sqrt(mean_n)).

    Returns
    -------
    n_values : numpy.ndarray
        Array [0, 1, 2, ..., n_max].
    probabilities : numpy.ndarray
        P(n) for each value, from the Poisson probability mass function.
    """
    n_values = np.arange(0, n_max + 1)
    probabilities = stats.poisson.pmf(n_values, mu=mean_n)

    # The PMF values for n=0..n_max won't sum to EXACTLY 1 (a true Poisson
    # distribution has infinite support), just very close to 1 as long as
    # n_max is large enough. We check this rather than silently
    # renormalizing, so that a poorly chosen n_max produces a visible error
    # instead of quietly biased results.
    total = probabilities.sum()
    if total < 0.999:
        raise ValueError(
            f"n_max={n_max} is too small for mean_n={mean_n}: only "
            f"{total:.4f} of the Poisson probability mass is captured. "
            "Increase n_max."
        )

    return n_values, probabilities


def negative_binomial_distribution(mean_n, dispersion_r, n_max):
    """
    Build a negative binomial distribution over receptor counts.

    WHY NEGATIVE BINOMIAL, NOT POISSON, FOR REAL DATA:
    The handoff notes real EV receptor distributions are "zero-inflated and
    overdispersed" -- meaning (a) an unusually large fraction of vesicles
    carry exactly zero of the marker (more than a Poisson distribution with
    the same mean would predict), and (b) the spread/variance of receptor
    counts across vesicles is bigger than a Poisson distribution allows
    (for a Poisson distribution, variance always exactly equals the mean --
    that's actually its defining mathematical property). The negative
    binomial distribution is a standard, well-established way to model
    "Poisson-like count data with extra variance and often extra zeros" --
    it's widely used for exactly this kind of biological count data (e.g.
    RNA-seq read counts have the same issue and the field standardized on
    negative binomial models for that reason).

    The negative binomial distribution has two parameters instead of
    Poisson's one: the mean (mean_n), and a "dispersion" parameter
    (dispersion_r) that controls how much EXTRA variance there is beyond
    Poisson. As dispersion_r gets very large, the negative binomial
    distribution actually converges to a Poisson distribution (This is
    a useful sanity check we'll also use as a secondary consistency test.)
    Smaller dispersion_r means more overdispersion/more zero-inflation.

    Parameters
    ----------
    mean_n : float
        Desired mean receptor count.
    dispersion_r : float
        Dispersion parameter (sometimes called "size" or "r" in different
        conventions). Must be > 0. Smaller = more overdispersed/more
        zero-inflated.
    n_max : int
        Highest receptor count to include explicitly (same purpose as in
        poisson_distribution).

    Returns
    -------
    n_values : numpy.ndarray
        Array [0, 1, 2, ..., n_max].
    probabilities : numpy.ndarray
        P(n) for each value, from the negative binomial probability mass
        function.
    """
    # scipy.stats.nbinom is parameterized by (n_param, p_param), NOT
    # directly by mean and dispersion, so we need to convert.
    # Using the mean/dispersion convention: mean = r * (1-p) / p
    # Solving for p given r=dispersion_r and mean=mean_n:
    #   p = r / (r + mean)
    p_param = dispersion_r / (dispersion_r + mean_n)

    n_values = np.arange(0, n_max + 1)
    probabilities = stats.nbinom.pmf(n_values, dispersion_r, p_param)

    total = probabilities.sum()
    if total < 0.999:
        raise ValueError(
            f"n_max={n_max} is too small for mean_n={mean_n}, "
            f"dispersion_r={dispersion_r}: only {total:.4f} of the "
            "probability mass is captured. Increase n_max."
        )

    return n_values, probabilities


# ---------------------------------------------------------------------------
# JOINT (DIAMETER, COPY-NUMBER) POPULATION MODEL
# ---------------------------------------------------------------------------
"""
WHY THIS SECTION EXISTS
---------------------------
Everything above assumes a population is described by a single distribution
P(n) over receptor counts, with no reference to vesicle size. That was fine
before the Lennon Figure 2/3 digitization, because no dataset had per-
vesicle DIAMETER paired with per-vesicle COPY NUMBER. Now that
ev_capture_data/ has exactly that (each row is one real or reconstructed
vesicle, with both a diameter_nm and a copies_per_ev value), the population
model can be extended to use the JOINT distribution directly, and to run
the diameter-aware capture_probability_patch() from binding_model.py across
an actual empirical (diameter, copy-number) point cloud rather than an
assumed P(n) shape.

This directly answers "Immediate next steps" item 1 from the v3 handoff:
whether the contact-patch-adjusted model (which can run in either the
n-only regime, patch_exponent=2, or the density-driven regime,
patch_exponent=0) changes the yield/selectivity comparison between two
real, matched vesicle populations (e.g. the gated PDAC-enriched plasma
population vs. the healthy plasma population).
"""

from binding_model import capture_probability_patch


def joint_population_yield(n_values, diameter_values, sigma, k_geom, k_intra, z,
                            patch_exponent=2.0):
    """
    Compute the population-average capture probability (yield) for an
    EMPIRICAL population of vesicles, each with its own (diameter,
    receptor-count) pair -- e.g. one row per vesicle from a digitized
    ev_capture_data/ CSV -- rather than an assumed analytic P(n)
    distribution.

    Unlike population_yield() above (which takes a probability-weighted
    discrete distribution), this takes a plain SAMPLE of individual
    vesicles (as you'd get directly from a digitized scatter plot or real
    single-vesicle measurement) and averages capture_probability_patch()
    over that sample directly -- each vesicle counts equally, so this is
    just a plain mean over the data, no explicit P(n) needs to be fit
    first.

    Parameters
    ----------
    n_values : array
        Receptor (copy) count for each vesicle in the sample. Same length
        as diameter_values -- n_values[i] and diameter_values[i] describe
        the SAME vesicle (i.e., this expects paired data, not two
        independent marginal distributions).
    diameter_values : array
        Diameter for each vesicle in the sample, same units used
        elsewhere (ev_capture_data/ uses nanometers).
    sigma : float
        Ligand surface density.
    k_geom : float
        Contact-patch geometric constant (see binding_model.py).
    k_intra : float
        Bond-strength parameter.
    z : float
        Activity/prefactor.
    patch_exponent : float, optional
        Contact-patch scaling exponent (default 2.0 = soft-adhesion limit,
        diameter-independent, matches the original size-blind model
        exactly; use 0.0 for the fixed-patch/density-driven alternative).
        See binding_model.py's module discussion for the full physical
        reasoning behind these two regimes.

    Returns
    -------
    float
        Population-average capture probability (yield), between 0 and 1 --
        the plain mean of capture_probability_patch() over every vesicle
        in the sample.
    """
    n_values = np.asarray(n_values, dtype=float)
    diameter_values = np.asarray(diameter_values, dtype=float)

    if len(n_values) != len(diameter_values):
        raise ValueError(
            f"n_values and diameter_values must be paired (same vesicle, "
            f"same index) and therefore the same length: got "
            f"{len(n_values)} vs {len(diameter_values)}."
        )
    if len(n_values) == 0:
        raise ValueError("Cannot compute a yield over an empty population.")

    theta_per_vesicle = capture_probability_patch(
        n_values, diameter_values, sigma, k_geom, k_intra, z, patch_exponent
    )
    return float(np.mean(theta_per_vesicle))


def compare_patch_regimes(n_values, diameter_values, sigma, k_geom, k_intra, z):
    """
    Convenience function: compute joint_population_yield() for the SAME
    population under both named patch_exponent regimes (2.0 = soft-
    adhesion/size-blind, 0.0 = fixed-patch/density-driven), and return
    both, so the two can be compared directly. This is the concrete tool
    for checking whether a conclusion (e.g. "population A captures better
    than population B") depends on which contact-patch assumption is used
    -- exactly the check the v3 handoff's "Immediate next steps" item 1
    calls for.

    Parameters
    ----------
    n_values, diameter_values, sigma, k_geom, k_intra, z :
        Same meaning as in joint_population_yield().

    Returns
    -------
    dict with keys:
        "yield_size_blind" : yield at patch_exponent=2.0
        "yield_density_driven" : yield at patch_exponent=0.0
    """
    return {
        "yield_size_blind": joint_population_yield(
            n_values, diameter_values, sigma, k_geom, k_intra, z, patch_exponent=2.0
        ),
        "yield_density_driven": joint_population_yield(
            n_values, diameter_values, sigma, k_geom, k_intra, z, patch_exponent=0.0
        ),
    }
