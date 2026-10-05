"""
monte_carlo.py

A brute-force stochastic simulation of the same physical situation
binding_model.py describes analytically, used purely as a cross-check
(Aim 2, validation #1: "analytical vs. Monte Carlo must agree").

WHAT "BRUTE FORCE" MEANS HERE
-------------------------------
Instead of solving the math in closed form, we directly simulate the random
process that the math is describing, many times, and see what fraction of
the time a vesicle ends up "captured." If our closed-form formula in
binding_model.py is correct, the two methods should agree (within the
statistical noise you'd expect from doing something random a finite number
of times).

WHAT RANDOM PROCESS ARE WE SIMULATING?
-----------------------------------------
Recall from binding_model.py that we're using the "large colloid / mobile
receptors AND mobile ligands" regime from Curk, Dobnikar & Frenkel (their
Eq 1.15). In that regime, the bound-state partition function has the form

    q_b(n, sigma) = exp(n * sigma * K_INTRA) - 1

A partition function of the form exp(mean) - 1 is exactly what you get from
a POISSON PROCESS: if the number of successful bond-formation "attempts"
that actually turn into bonds is a Poisson random variable with mean
lambda = n * sigma * K_INTRA, then:

    P(zero bonds form)     = exp(-lambda)
    P(at least one bond)   = 1 - exp(-lambda)

This is the physical picture: each of the vesicle's n receptors has many
chances to encounter a ligand as it explores the contact patch (because
receptors and ligands are both mobile -- see the handoff's mechanism (a),
lateral diffusion). The NUMBER of bond-forming encounters that succeed is
random, and its average scales with both how many receptors there are (n)
and how much ligand is available (sigma). That's exactly a Poisson process
with mean n * sigma * K_INTRA.

So, in this Monte Carlo, for a single simulated vesicle we:
    1. Compute lambda = n * sigma * K_INTRA (the average number of bonds
       we'd expect to form).
    2. Draw a random number of actual bonds formed from a Poisson
       distribution with that mean.
    3. Say the vesicle is "attached" (has at least one bond) if that number
       is >= 1.

Then, separately, the Langmuir wrapper theta = z*q_b / (1 + z*q_b) accounts
for whether the vesicle even gets the opportunity to start forming bonds in
the first place (transport/concentration/incubation-time effects bundled
into z). We simulate that as a separate coin flip: with probability
z*q_b / (1 + z*q_b + ...) -- see the docstring of simulate_capture() below
for exactly how the two pieces combine, because getting this right matters:
plugging the two steps together naively would NOT reproduce the Langmuir
formula. We work out the correct way to combine them below.

HOW TO COMBINE THE TWO STEPS CORRECTLY
-----------------------------------------
We want P(capture) = z*q_b / (1 + z*q_b). A clean way to simulate exactly
this probability with two sequential random draws is:

    Draw U ~ Uniform(0, 1).
    Capture happens if U < z*q_b / (1 + z*q_b).

This is the most direct simulation: draw one random number, compare to the
formula's own probability. This does NOT yet test anything, though --
it would "agree" by construction no matter what formula we typed in,
because we are just re-implementing the same equation with a random number
generator glued on.

To make this an honest, independent check, our Monte Carlo needs to
build P(capture) up from more primitive random ingredients -- ones that do
NOT already assume the Langmuir formula -- and then check that the
resulting frequency matches theta. Here is the decomposition we use,
which follows directly from the algebra of z*q_b/(1+z*q_b):

    Let q_b = exp(lambda) - 1, where lambda = n * sigma * k_intra.
    z*q_b / (1 + z*q_b)
        = z*(exp(lambda)-1) / (1 + z*(exp(lambda)-1))

There isn't a simpler "two independent coin flips" decomposition of this
exact expression in general (the Langmuir wrapper couples z and q_b
nonlinearly). So instead, our Monte Carlo tests the part of the model that
actually contains the physics we care about -- the exp(lambda)-1 bond-count
term -- by simulating the Poisson bond-counting process explicitly (steps
1-3 above) across MANY vesicles and MANY random ligand "encounter" trials,
and comparing the resulting fraction with >=1 bond against the analytical
q_b-derived quantity q_b/(1+q_b) (i.e., we set z=1 for this specific check,
which isolates the combinatorial/multivalency part of the model from the
transport/concentration part). This is a meaningful, non-circular test:
the Monte Carlo never "knows" the formula exp(lambda)-1 -- it only knows
"draw a Poisson-distributed number of bonds with mean n*sigma*k_intra, and
check if it's >= 1." If our analytical formula for q_b is right, the two
must agree.

We separately test the z-dependent Langmuir wrapper on its own in the unit
tests (validating against the textbook Langmuir isotherm), since that part
of the model is standard and not itself in question.
"""

import numpy as np


def simulate_bond_formation(n, sigma, k_intra, n_trials, rng=None):
    """
    Simulate whether a vesicle with n receptors forms at least one bond
    with a surface of ligand density sigma, using a Poisson bond-count
    model, repeated n_trials times (i.e., n_trials independent vesicles,
    all with the same n and sigma).

    This function does NOT use the q_bound() formula from binding_model.py
    at all -- it only uses the physical assumption behind it (bonds forming
    via a Poisson process with mean n*sigma*k_intra). That's what makes
    comparing its output to q_bound() a genuine check rather than a
    circular one.

    Parameters
    ----------
    n : int
        Number of receptors on the vesicle (must be a whole number here --
        this is a real simulated vesicle, not the continuous math version).
    sigma : float
        Ligand surface density.
    k_intra : float
        Binding-strength parameter (same meaning as in binding_model.py).
    n_trials : int
        How many independent simulated vesicles to run. More trials =
        less statistical noise, at the cost of more compute time.
    rng : numpy.random.Generator, optional
        A random number generator. If not given, a fresh default one is
        created. Passing your own rng makes results reproducible (same
        rng seed -> same simulated outcome every time), which is useful
        for debugging.

    Returns
    -------
    fraction_bound : float
        The fraction of the n_trials simulated vesicles that formed at
        least one bond. This is our Monte Carlo estimate of q_b/(1+q_b)
        (see module docstring for why that's the right comparison).
    bond_counts : numpy.ndarray of shape (n_trials,)
        The actual number of bonds simulated for each vesicle, in case
        you want to inspect the full distribution (not just the
        bound/unbound outcome).
    """
    if rng is None:
        rng = np.random.default_rng()

    lam = n * sigma * k_intra  # mean number of bonds expected, per vesicle
    if lam < 0:
        raise ValueError(
            f"Computed a negative mean bond count (lambda={lam}). "
            "Check that n, sigma, and k_intra are all non-negative."
        )

    bond_counts = rng.poisson(lam=lam, size=n_trials)
    is_bound = bond_counts >= 1
    fraction_bound = np.mean(is_bound)

    return fraction_bound, bond_counts


def simulate_capture_probability(n, sigma, k_intra, n_trials, rng=None):
    """
    Convenience wrapper: run simulate_bond_formation() and return just the
    fraction_bound (the Monte Carlo estimate of theta at z=1, i.e. of
    q_b / (1 + q_b)).

    Parameters are identical to simulate_bond_formation().

    Returns
    -------
    float
        Monte Carlo estimate of q_b / (1 + q_b) for the given n, sigma,
        k_intra.
    """
    fraction_bound, _ = simulate_bond_formation(n, sigma, k_intra, n_trials, rng=rng)
    return fraction_bound


def monte_carlo_standard_error(fraction_bound, n_trials):
    """
    Estimate the statistical uncertainty (standard error) on a Monte Carlo
    fraction_bound estimate, using the standard formula for the standard
    error of a proportion:

        SE = sqrt( p * (1 - p) / n_trials )

    This tells us how much wiggle room to allow when comparing the Monte
    Carlo result to the analytical formula -- if the two disagree by much
    more than a few standard errors, that's a sign of a real bug, not just
    random noise.

    Parameters
    ----------
    fraction_bound : float
        The estimated fraction (between 0 and 1) from a Monte Carlo run.
    n_trials : int
        How many trials were used to produce that estimate.

    Returns
    -------
    float
        The estimated standard error.
    """
    p = fraction_bound
    return np.sqrt(p * (1 - p) / n_trials)
