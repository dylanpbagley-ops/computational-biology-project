"""
fig3_stemcell_marker_crossover.py

Write-up figure 3 (Claim 1, EV foundation, third independent dataset):
the CD81/CD63/CD73 marker-level yield crossover on reconstructed
stem-cell (MSC-EV) data. This is the figure version of
sweep_marker_ranking_across_plausibility_range()'s headline finding: the
naive fold-change ranking (CD63 > CD81 > CD73) holds at 7 of 8 tested
K_INTRA values, but INVERTS at the top (CD81 > CD63) at the single
largest, most-saturating K_INTRA -- independently verified across 4
random seeds in the project's own prior work.

Chart choice: a line plot of marker-level yield vs. K_INTRA (log x-axis),
one line per marker (3 series -- comfortably within the "color alone is
fine" 1-3 series band), so the crossing point is directly visible rather
than only implied by a boolean ranking-match table. "Marker-level" means,
at each K_INTRA, each marker's yield is its BEST-PLACED tissue's yield --
the same convention marker_level_ranking() uses, because CD73's three
tissues are statistically indistinguishable (identical reconstruction
parameters) and their relative jitter among themselves is sampling noise,
not signal (see that function's own docstring). Plotting all nine
(marker, tissue) populations would bury the real 3-marker crossover this
finding is actually about inside noise from a comparison this project's
own documentation explicitly warns is not real.

Uses reconstruct_marker_population() and stable_joint_population_yield()
directly, at the SAME k_intra grid and z value
sweep_marker_ranking_across_plausibility_range() uses -- no new
computation.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig3_stemcell_marker_crossover.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import clean_spines, style_title, section_badge
from stem_cell_marker_ranking_extension import (
    reconstruct_marker_population,
    stable_joint_population_yield,
    naive_fold_change_ranking,
    sweep_marker_ranking_across_plausibility_range,
    marker_level_ranking,
    POSITIVITY_FRACTION_RANGE,
    COPIES_PER_POSITIVE_PARTICLE,
    SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
)
from k_intra_z_plausibility_bounds import k_intra_bounds_from_real_sigma

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")
Z_VALUE = 1e-4  # same representative dilute-regime value used throughout this file

# Fixed categorical color assignment, in the field's naive-ranking order
# (CD63 first/highest naive score, then CD81, then CD73) -- NOT reassigned
# by the nonlinear model's ranking, so a reader tracking "the orange line"
# can follow it whether it's winning or losing at a given K_INTRA.
MARKER_COLORS = {
    "CD63": "#B0563A",  # slot 1 -- naive-ranking leader
    "CD81": "#2E5EAA",  # slot 2 -- the one that overtakes at high K_INTRA
    "CD73": "#5B8C5A",  # slot 3 -- stays last throughout
}


def compute_marker_level_yields(n_samples=5000):
    k_intra_values = sorted(set(k_intra_bounds_from_real_sigma().values()))

    # Reconstruct every (marker, tissue) population once, same seeding
    # convention as sweep_marker_ranking_across_plausibility_range().
    populations = {}
    pair_index = 0
    for marker in POSITIVITY_FRACTION_RANGE:
        for tissue in COPIES_PER_POSITIVE_PARTICLE[marker]:
            populations[(marker, tissue)] = reconstruct_marker_population(
                marker, tissue, n_samples=n_samples, seed=pair_index
            )
            pair_index += 1

    # marker -> list of yields (one per k_intra), each the BEST tissue's
    # yield at that k_intra (marker_level_ranking's convention).
    marker_yields = {m: [] for m in POSITIVITY_FRACTION_RANGE}

    for k_intra in k_intra_values:
        per_pair_yield = {}
        for pair, (n_vals, d_vals) in populations.items():
            y = stable_joint_population_yield(
                n_vals, d_vals, SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, 1.0, k_intra, Z_VALUE,
                patch_exponent=2.0,
            )
            per_pair_yield[pair] = y
        for marker in POSITIVITY_FRACTION_RANGE:
            best = max(per_pair_yield[(marker, t)] for t in COPIES_PER_POSITIVE_PARTICLE[marker])
            marker_yields[marker].append(best)

    return k_intra_values, marker_yields


def plot_crossover(k_intra_values, marker_yields, save_path):
    naive = naive_fold_change_ranking()
    naive_marker_order = []
    for marker, tissue in naive["ranking_best_first"]:
        if marker not in naive_marker_order:
            naive_marker_order.append(marker)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for marker in ["CD63", "CD81", "CD73"]:
        ax.plot(
            k_intra_values, marker_yields[marker],
            "-o", color=MARKER_COLORS[marker], linewidth=2, markersize=6,
            label=marker, zorder=3,
        )

    # Mark the crossover point directly: the largest k_intra where CD81
    # overtakes CD63 (per this project's own documented finding). Placed
    # well inside the axes (below-left of the point, not above the plot)
    # so it can never collide with the title -- checked visually after an
    # earlier version placed it above the axes and it overlapped the
    # second title line.
    crossover_idx = len(k_intra_values) - 1  # documented: only at the single largest K_INTRA
    if marker_yields["CD81"][crossover_idx] > marker_yields["CD63"][crossover_idx]:
        ax.annotate(
            "CD81 overtakes CD63\n(only at the largest $K_{INTRA}$)",
            xy=(k_intra_values[crossover_idx], marker_yields["CD81"][crossover_idx]),
            xytext=(-190, -55), textcoords="offset points",
            fontsize=9.5, color="0.25",
            arrowprops=dict(arrowstyle="->", color="0.4", linewidth=1.2,
                             connectionstyle="arc3,rad=0.15"),
        )

    ax.set_xscale("log")
    ax.set_xlabel(r"$K_{INTRA}$  (per-bond binding strength, cm$^2$/antibody)", fontsize=11)
    ax.set_ylabel("Marker-level population yield\n(best tissue source)", fontsize=11)
    style_title(
        ax,
        "Marker-level population yield vs. $K_{INTRA}$ for CD63, CD81, CD73",
        fontsize=12.5,
    )
    clean_spines(ax, keep=("left", "bottom"))
    ax.legend(loc="upper left", fontsize=10, frameon=False, title="Marker (best tissue)")
    ax.grid(axis="y", alpha=0.25, linewidth=0.8)
    ymax = max(max(v) for v in marker_yields.values())
    ax.set_ylim(top=ymax * 1.18)

    result = sweep_marker_ranking_across_plausibility_range()
    naive_marker_level = marker_level_ranking(result["naive_ranking"])
    marker_level_matches = [
        marker_level_ranking(r) == naive_marker_level for r in result["nonlinear_rankings"]
    ]
    n_match = sum(marker_level_matches)
    n_total = len(marker_level_matches)
    caption = (
        f"Naive fold-change ranking (positivity fraction $\\times$ mean copies/positive "
        f"particle): {' > '.join(naive_marker_order)}. Nonlinear model matches this "
        f"full ranking at {n_match}/{n_total} tested $K_{{INTRA}}$ values; inverts only "
        f"at the single largest (most-saturating) value. Reconstructed from Cytotherapy "
        f"2025 summary statistics (AD/UC/BM tissue sources), 4-seed-verified."
    )
    fig.text(0.5, -0.06, caption, ha="center", fontsize=9.2, wrap=True)

    section_badge(fig, "EV", 3)
    # No add_callout() here: the crossover point already has a hand-placed
    # arrowed annotation above (positioned specifically to avoid the title),
    # which does the same "box/arrow at the key data feature" job -- a
    # second callout would be redundant.

    fig.tight_layout()
    fig.subplots_adjust(top=0.84)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    print(f"Naive marker order: {naive_marker_order}")
    print(f"Nonlinear ranking matches naive at {n_match}/{n_total} K_INTRA values")
    for marker in ["CD63", "CD81", "CD73"]:
        print(f"{marker} yields: {[f'{y:.4f}' for y in marker_yields[marker]]}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    k_intra_values, marker_yields = compute_marker_level_yields()
    plot_crossover(k_intra_values, marker_yields,
                    save_path=os.path.join(FIGURES_DIR, "fig3_stemcell_marker_crossover.png"))
