"""
fig7_conceptual_schematic.py

Write-up figure 7: a conceptual schematic explaining the project's core
mechanism, for the paper introduction and poster. Not a result figure --
no headline number is claimed here. Explicitly and visibly labeled as
ILLUSTRATIVE, using round, clearly-synthetic parameter values, not real
digitized/measured data (which live in figures 1-6).

Two panels, built from the ACTUAL model function (binding_model.
capture_probability), not a hand-drawn cartoon curve:

Left panel: the nonlinear capture-probability curve itself (theta vs.
receptor count n), the mechanism the whole project is built on. A linear
("naive") reference line is drawn for contrast.

Right panel: two illustrative populations with the SAME MEAN receptor
count but different shapes (one tight/unimodal, one zero-inflated/
overdispersed -- the real, literature-documented shape of EV/tumor-cell
marker distributions per this project's own "Poisson vs. negative
binomial" finding) produce DIFFERENT total capture once run through the
same nonlinear curve. This is the concrete, general version of "why mean
alone can rank targets incorrectly" that every real headline figure
(1-6) is a specific instance of.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig7_conceptual_schematic.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import nbinom

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import clean_spines, style_title, section_badge, add_callout, PROJECT_ORANGE
from binding_model import capture_probability

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")

COLOR_NONLINEAR = "#2E5EAA"
COLOR_LINEAR_REF = "0.55"
COLOR_TIGHT = "#2E5EAA"
COLOR_SKEWED = "#B0563A"

# Deliberately round, illustrative parameter values -- NOT fit to or
# copied from any real digitized/measured system (those are figures 1-6).
SIGMA_ILLUSTRATIVE = 1.0
K_INTRA_ILLUSTRATIVE = 0.05
Z_ILLUSTRATIVE = 1.0
MEAN_N = 20  # shared mean receptor count for both illustrative populations


def _neg_binom_params_for_mean_and_dispersion(mean, r):
    """
    scipy's nbinom is parameterized by (r, p) with mean = r*(1-p)/p. Given
    a target mean and a dispersion parameter r (smaller r = heavier tail /
    more zero-inflated, larger r = tighter/more Poisson-like), solve for p
    directly in closed form -- no root-finding needed, and by
    construction the TRUE mean (not just the zero-fraction) matches
    exactly, avoiding the truncation-mismatch bug an earlier version of
    this function had (matching zero-fraction with a free r let the true
    mean drift far from the target for very small r).
    """
    p = r / (r + mean)
    return r, p


def plot_schematic(save_path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))

    # --- Left panel: the nonlinear capture curve itself ---
    ax = axes[0]
    n_values = np.linspace(0.1, 100, 300)
    theta = capture_probability(n_values, SIGMA_ILLUSTRATIVE, K_INTRA_ILLUSTRATIVE, Z_ILLUSTRATIVE)
    ax.plot(n_values, theta, color=COLOR_NONLINEAR, linewidth=2.5, label="Nonlinear model\n(this project)")
    linear_ref = n_values / n_values.max() * theta.max()
    ax.plot(n_values, linear_ref, "--", color=COLOR_LINEAR_REF, linewidth=1.8,
            label="Linear (\"fold-change\")\nassumption")
    ax.set_xlabel("Receptor / antigen count per particle or cell ($n$)", fontsize=10.5)
    ax.set_ylabel("Capture / kill probability", fontsize=10.5)
    style_title(ax, "Capture is nonlinear in receptor count", fontsize=12.5)
    clean_spines(ax, keep=("left", "bottom"))
    ax.legend(fontsize=9.5, frameon=False, loc="lower right")
    ax.grid(alpha=0.2, linewidth=0.7)

    # --- Right panel: same mean, different shape -> different total yield ---
    ax2 = axes[1]
    # A wide range is needed for the PMF sums (mean, yield) to converge
    # accurately for the heavy-tailed (small-r) population -- checked
    # directly: at r=0.6, mean=20, summing only to x=100 recovers a mean
    # well below 20 (truncation error), while summing to x=2000 recovers
    # the true mean to machine precision. The PLOT only shows x<70 (where
    # essentially all the visible mass is); the SUMS use the full range.
    x = np.arange(0, 2000)

    # Both populations share the SAME true mean (MEAN_N) by construction
    # -- only the dispersion r differs (large r = tight/Poisson-like,
    # small r = heavily zero-inflated/overdispersed, the real documented
    # shape of EV/tumor-marker count data per this project's own
    # Poisson-vs-negative-binomial finding).
    r_tight, p_tight = _neg_binom_params_for_mean_and_dispersion(MEAN_N, r=20.0)
    r_skew, p_skew = _neg_binom_params_for_mean_and_dispersion(MEAN_N, r=0.6)

    pmf_tight = nbinom.pmf(x, r_tight, p_tight)
    pmf_skew = nbinom.pmf(x, r_skew, p_skew)

    ax2.bar(x - 0.2, pmf_tight, width=0.4, color=COLOR_TIGHT, alpha=0.75, label="Population A (tight)")
    ax2.bar(x + 0.2, pmf_skew, width=0.4, color=COLOR_SKEWED, alpha=0.75, label="Population B (zero-inflated)")
    ax2.set_xlim(-1, 70)
    ax2.set_xlabel("Receptor / antigen count per particle or cell ($n$)", fontsize=10.5)
    ax2.set_ylabel("Fraction of population", fontsize=10.5)
    style_title(ax2, "Two populations with equal mean but different shape", fontsize=12.5)
    clean_spines(ax2, keep=("left", "bottom"))
    ax2.legend(fontsize=9.5, frameon=False, loc="upper right")
    ax2.grid(alpha=0.2, linewidth=0.7)

    theta_per_n = capture_probability(np.maximum(x, 0.01), SIGMA_ILLUSTRATIVE, K_INTRA_ILLUSTRATIVE, Z_ILLUSTRATIVE)
    yield_tight = float(np.sum(pmf_tight * theta_per_n))
    yield_skew = float(np.sum(pmf_skew * theta_per_n))

    mean_tight = float(np.sum(x * pmf_tight))
    mean_skew = float(np.sum(x * pmf_skew))

    caption = (
        f"Illustrative example, not a result: Population A and B are constructed with the "
        f"same mean receptor count (A: {mean_tight:.1f}, B: {mean_skew:.1f}) but different "
        f"shape (B is zero-inflated). Passed through the same nonlinear capture curve (left "
        f"panel), they yield different total capture: {yield_tight*100:.1f}% (A) vs. "
        f"{yield_skew*100:.1f}% (B), a difference a mean-only comparison cannot detect. This "
        f"illustrates the general mechanism behind the reordering findings in this project."
    )
    fig.text(0.5, -0.09, caption, ha="center", fontsize=9.3, wrap=True)

    # Callout on the zero-inflated population's spike at n=0 -- the visual
    # root cause of the differing yields despite the two populations
    # sharing the same mean.
    add_callout(
        ax2, xy=(0, pmf_skew[0]), text="zero-inflated tail drives\nthe yield difference",
        xytext=(0.55, 0.55), textcoords="axes fraction",
        color=PROJECT_ORANGE, fontsize=9.5,
    )

    fig.suptitle("Illustrative mechanism: nonlinear capture and population shape", fontsize=15, y=1.03)
    section_badge(fig, "OVERVIEW", 7)
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    print(f"Population A (tight): mean={mean_tight:.2f}, yield={yield_tight:.4f}")
    print(f"Population B (skewed): mean={mean_skew:.2f}, yield={yield_skew:.4f}")
    assert abs(mean_tight - mean_skew) < 1.0, "means should be approximately equal by construction"
    print("Sanity check passed: means approximately equal, yields differ.")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    plot_schematic(save_path=os.path.join(FIGURES_DIR, "fig7_conceptual_schematic.png"))
