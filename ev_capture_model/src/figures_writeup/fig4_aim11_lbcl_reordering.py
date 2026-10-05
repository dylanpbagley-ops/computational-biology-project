"""
fig4_aim11_lbcl_reordering.py

Write-up figure 4 (Claim 2, CAR-T pivot): the Aim 11 pre-treatment vs.
at-relapse LBCL CD19 reordering check. This is the figure version of
sweep_pre_vs_relapse_across_full_grid()'s headline finding: across a
4 (real k_intra fit) x 9 (kill_threshold, 0.1-0.9) = 36-point grid,
naive mean-density tracking is reordered by the distribution-aware
escape-fraction model at 2/36 (5.6%) of grid points.

Chart choice: a STATUS GRID (2 discrete outcomes -- "matches naive" vs.
"reordered" -- not a continuous heatmap). Using a continuous colormap
here would misstate the data: this is a binary per-cell outcome, not a
magnitude field, so a sequential/viridis fill would imply a false
continuity between cells. Status colors (a muted "matches" gray-blue vs.
a flagged "reordered" warm accent) are used instead, reserved for exactly
this good/notable-state distinction and not reused as a general
categorical palette elsewhere in this figure set.

Uses sweep_pre_vs_relapse_across_full_grid() directly -- no new
computation, just re-shaped into a grid for plotting.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig4_aim11_lbcl_reordering.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import clean_spines, style_title, section_badge, add_callout, PROJECT_ORANGE
from aim11_antigen_target_ranking import sweep_pre_vs_relapse_across_full_grid

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")

COLOR_MATCHES = "#C7CDD6"   # neutral/status "no change" -- matches naive ranking
COLOR_REORDERED = "#C0392B"  # status "flagged" -- distribution-aware model disagrees

FIT_LABELS = {
    "CLN-978_cytolysis": "CLN-978\n(cytolysis)",
    "CLN-978_ifn_gamma": "CLN-978\n(IFN-γ)",
    "blinatumomab_cytolysis": "blinatumomab\n(cytolysis)",
    "blinatumomab_ifn_gamma": "blinatumomab\n(IFN-γ)",
}


def build_status_grid():
    result = sweep_pre_vs_relapse_across_full_grid()
    grid_results = result["grid_results"]

    fit_names = sorted(set(r["fit_name"] for r in grid_results), key=lambda x: x)
    thresholds = sorted(set(r["kill_threshold"] for r in grid_results))

    status = np.zeros((len(fit_names), len(thresholds)), dtype=int)  # 0=matches, 1=reordered
    for r in grid_results:
        i = fit_names.index(r["fit_name"])
        j = thresholds.index(r["kill_threshold"])
        status[i, j] = 0 if r["matches_naive"] else 1

    return fit_names, thresholds, status, result


def plot_status_grid(fit_names, thresholds, status, result, save_path):
    fig, ax = plt.subplots(figsize=(8.5, 4.5))

    cmap = ListedColormap([COLOR_MATCHES, COLOR_REORDERED])
    ax.pcolormesh(
        np.arange(len(thresholds) + 1), np.arange(len(fit_names) + 1),
        status, cmap=cmap, vmin=0, vmax=1,
        edgecolors="white", linewidth=2,
    )

    ax.set_xticks(np.arange(len(thresholds)) + 0.5)
    ax.set_xticklabels([f"{t:.1f}" for t in thresholds], fontsize=9)
    ax.set_yticks(np.arange(len(fit_names)) + 0.5)
    ax.set_yticklabels([FIT_LABELS.get(f, f) for f in fit_names], fontsize=9)

    ax.set_xlabel("Kill-probability threshold (fraction of cells that must be\n"
                   "reliably killed to count as \"not escaped\")", fontsize=10.5)
    ax.set_ylabel("Calibration fit", fontsize=10.5)
    style_title(
        ax,
        "Reordering status by kill-probability threshold and calibration fit\n"
        "(pre-treatment vs. at-relapse CD19, LBCL cohort)",
        fontsize=13,
    )
    clean_spines(ax, keep=())

    # Callout on one of the reordered cells -- the single most important
    # feature of a status grid whose whole point is "which cells flip."
    reordered_ij = np.argwhere(status == 1)
    if len(reordered_ij) > 0:
        i0, j0 = reordered_ij[0]
        add_callout(
            ax, xy=(j0 + 0.5, i0 + 0.5), text="reordered",
            xytext=(0.82, 0.90), textcoords="axes fraction",
            color=PROJECT_ORANGE, fontsize=10.5,
        )

    legend_elements = [
        Patch(facecolor=COLOR_MATCHES, edgecolor="white", label="Matches naive ranking"),
        Patch(facecolor=COLOR_REORDERED, edgecolor="white", label="Reordered by the model"),
    ]
    ax.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, -0.32),
              ncol=2, fontsize=10, frameon=False)

    n_reordered = result["n_reordered"]
    n_total = result["n_grid_points"]
    frac = result["reordering_fraction"]
    caption = (
        f"{n_reordered}/{n_total} grid points ({frac*100:.1f}%) show the distribution-aware "
        f"model producing a different ranking than tracking mean CD19 density alone "
        f"(5,810 to 2,021 molecules/cell). Based on summary-statistic reconstructions, "
        f"not a validated clinical outcome prediction."
    )
    fig.text(0.5, -0.16, caption, ha="center", fontsize=9.2, wrap=True)

    section_badge(fig, "CAR-T", 4)

    fig.tight_layout()
    fig.subplots_adjust(top=0.80)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    print(f"Reordered: {n_reordered}/{n_total} ({frac*100:.2f}%)")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    fit_names, thresholds, status, result = build_status_grid()
    plot_status_grid(fit_names, thresholds, status, result,
                      save_path=os.path.join(FIGURES_DIR, "fig4_aim11_lbcl_reordering.png"))
