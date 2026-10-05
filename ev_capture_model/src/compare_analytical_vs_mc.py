"""
compare_analytical_vs_mc.py

Aim 2, validation #1: run the analytical formula (binding_model.q_bound)
and the Monte Carlo simulation (monte_carlo.simulate_capture_probability)
across a range of receptor counts (n) and ligand densities (sigma), and
check that they agree.

As explained in monte_carlo.py's docstring, we compare the two methods at
z=1, where the analytical prediction simplifies to:

    theta_analytical = q_b / (1 + q_b),   q_b = exp(n * sigma * k_intra) - 1

This script:
    1. Picks a grid of (n, sigma) combinations to test.
    2. Computes the analytical prediction for each.
    3. Runs the Monte Carlo simulation for each (with enough trials to keep
       statistical noise small).
    4. Prints a table comparing the two, with the Monte Carlo standard
       error, so we can see by eye whether disagreements are just noise or
       a real problem.
    5. Makes a plot: capture probability vs. n, both curves overlaid, so
       the agreement (or lack of it) is visually obvious.

Run this from the command line with:  python compare_analytical_vs_mc.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")  # write plots to file rather than needing a display
import matplotlib.pyplot as plt

from binding_model import q_bound
from monte_carlo import simulate_capture_probability, monte_carlo_standard_error


def run_comparison(n_values, sigma, k_intra, n_trials, rng=None):
    """
    Compute analytical and Monte Carlo capture probabilities (at z=1) for
    each receptor count in n_values, at a fixed sigma and k_intra.

    Returns a list of dicts, one per n value, each containing:
        n, theta_analytical, theta_mc, mc_standard_error,
        difference, difference_in_standard_errors
    """
    if rng is None:
        rng = np.random.default_rng(seed=42)  # fixed seed = reproducible run

    results = []
    for n in n_values:
        qb = q_bound(n, sigma, k_intra)
        theta_analytical = qb / (1.0 + qb)

        theta_mc = simulate_capture_probability(n, sigma, k_intra, n_trials, rng=rng)
        se = monte_carlo_standard_error(theta_mc, n_trials)

        difference = theta_mc - theta_analytical
        # Avoid division by zero if se happens to be exactly 0 (can happen
        # if theta_mc is exactly 0 or 1, e.g. for extreme parameter values)
        diff_in_se = difference / se if se > 0 else float("nan")

        results.append({
            "n": n,
            "theta_analytical": theta_analytical,
            "theta_mc": theta_mc,
            "mc_standard_error": se,
            "difference": difference,
            "difference_in_standard_errors": diff_in_se,
        })

    return results


def print_comparison_table(results):
    """Pretty-print the comparison results as a table."""
    header = f"{'n':>8} | {'analytical':>12} | {'monte carlo':>12} | {'MC std err':>10} | {'diff (in SE)':>12}"
    print(header)
    print("-" * len(header))
    for r in results:
        print(
            f"{r['n']:>8.1f} | {r['theta_analytical']:>12.5f} | "
            f"{r['theta_mc']:>12.5f} | {r['mc_standard_error']:>10.5f} | "
            f"{r['difference_in_standard_errors']:>12.2f}"
        )


def plot_comparison(results, sigma, k_intra, save_path):
    """Overlay analytical and Monte Carlo curves and save the figure."""
    n_values = [r["n"] for r in results]
    theta_analytical = [r["theta_analytical"] for r in results]
    theta_mc = [r["theta_mc"] for r in results]
    mc_err = [r["mc_standard_error"] for r in results]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(n_values, theta_analytical, "-", color="tab:blue", linewidth=2,
             label="Analytical formula")
    ax.errorbar(n_values, theta_mc, yerr=mc_err, fmt="o", color="tab:orange",
                markersize=5, capsize=3, label="Monte Carlo simulation")

    ax.set_xlabel("Number of receptors on vesicle (n)")
    ax.set_ylabel("Capture probability (theta, at z=1)")
    ax.set_title(f"Analytical vs. Monte Carlo\n(sigma={sigma}, k_intra={k_intra})")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    print(f"\nSaved comparison plot to: {save_path}")


if __name__ == "__main__":
    # --- Parameters for this test run ---
    # These are arbitrary "reasonable" numbers chosen just to exercise the
    # model across a range where capture probability goes from near 0 to
    # near 1 -- not yet calibrated to any real biological system (that's
    # Aim 3's job). sigma and k_intra are chosen together so that
    # n * sigma * k_intra ranges from small to several, over the n values
    # tested below.
    SIGMA = 1.0
    K_INTRA = 0.05
    N_VALUES = np.array([1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 96, 128], dtype=float)
    N_TRIALS = 200_000  # Monte Carlo trials per n value; large enough to
                         # keep standard error small (see printed table)

    results = run_comparison(N_VALUES, SIGMA, K_INTRA, N_TRIALS)
    print_comparison_table(results)

    max_diff_se = max(
        abs(r["difference_in_standard_errors"]) for r in results
        if not np.isnan(r["difference_in_standard_errors"])
    )
    print(f"\nLargest disagreement across all n values: {max_diff_se:.2f} standard errors.")
    if max_diff_se < 4:
        print("PASS: analytical and Monte Carlo agree within statistical noise.")
    else:
        print("WARNING: disagreement larger than expected from noise alone -- investigate.")

    plot_comparison(results, SIGMA, K_INTRA, save_path="../figures/analytical_vs_mc.png")
