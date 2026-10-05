"""
fig1_superselectivity_heatmap_poster.py

POSTER version of Figure 1. Same underlying data/computation as
fig1_superselectivity_heatmap.py (imports build_alpha_grid/_log_edges from
it directly -- no recomputation, no duplicated logic), but re-rendered for
the physical board: bigger fonts throughout, and the small italic
technical caption replaced by ONE big bold plain-language takeaway banner
(the technical caption stays in the paper version).

Run from ev_capture_model/src/figures_writeup/:
    python3 fig1_superselectivity_heatmap_poster.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import (
    clean_spines, section_badge, add_callout, add_takeaway_banner,
    PROJECT_BLUE,
)
from fig1_superselectivity_heatmap import build_alpha_grid, _log_edges
from aim3_aim4_range_sweep import ALPHA_SUPERSELECTIVE_THRESHOLD

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")


def build(save_path):
    k_intra_values, z_values, grid, excluded = build_alpha_grid()

    fig, ax = plt.subplots(figsize=(11, 8.5))
    z_edges = _log_edges(z_values)
    k_edges = _log_edges(k_intra_values)
    masked_grid = np.ma.masked_invalid(grid)

    im = ax.pcolormesh(z_edges, k_edges, masked_grid, cmap="viridis",
                        shading="flat", edgecolors="white", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")

    Zg, Kg = np.meshgrid(z_values, k_intra_values)
    ax.scatter(Zg, Kg, s=8, color="white", alpha=0.55, linewidths=0, zorder=3)

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(r"Maximum selectivity $\alpha_{exact}$ reached", fontsize=15)
    cbar.ax.tick_params(labelsize=13)

    ax.set_xlabel(r"$z$  (vesicle concentration / transport activity)", fontsize=16)
    ax.set_ylabel(r"$K_{INTRA}$  (per-bond binding strength)", fontsize=16)
    ax.tick_params(labelsize=13)
    ax.set_title("Selectivity across the physically-plausible parameter grid",
                 fontsize=22, fontweight="bold", pad=14)
    clean_spines(ax, keep=("left", "bottom"))
    section_badge(fig, "EV", 1)

    n_reliable = int(np.sum(~np.isnan(grid)))
    n_superselective = int(np.sum(grid[~np.isnan(grid)] > ALPHA_SUPERSELECTIVE_THRESHOLD))
    frac = n_superselective / n_reliable if n_reliable else 0.0
    peak = float(np.nanmax(grid))
    peak_idx = np.unravel_index(np.nanargmax(grid), grid.shape)
    peak_k, peak_z = k_intra_values[peak_idx[0]], z_values[peak_idx[1]]

    add_callout(
        ax, xy=(peak_z, peak_k), text=f"peak α={peak:.1f}",
        xytext=(0.32, 0.74), textcoords="axes fraction",
        color=PROJECT_BLUE, fontsize=15,
    )

    fig.tight_layout()
    fig.subplots_adjust(top=0.86)  # leave room for the section badge above the title
    add_takeaway_banner(
        fig,
        f"Peak selectivity reaches {peak:.1f}x the baseline, exceeding the superselectivity "
        f"threshold in {frac*100:.0f}% of physically-plausible (K_INTRA, z) conditions.",
        color=PROJECT_BLUE, y=-0.10,
    )
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    build(os.path.join(FIGURES_DIR, "fig1_superselectivity_heatmap_poster.png"))
