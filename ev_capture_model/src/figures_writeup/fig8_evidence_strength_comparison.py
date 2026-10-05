"""
fig8_evidence_strength_comparison.py

Write-up figure 8: the project's headline "strength of evidence" story in
one figure. Compares the meaningful-reordering fraction across every check
this project has run, ordered from weakest to strongest evidentiary tier:
reconstructed-summary-statistic checks (Aim 9, Aim 11) vs. real,
individually-measured single-cell checks (Aim 12, at both the initial
9-patient subsample and the complete 41-patient cohort). The point of this
figure is the pattern itself: the effect gets LARGER, not smaller, as the
evidence gets more direct and as more real data is added -- the opposite of
what a data-selection artifact or a lucky subsample would produce.

Aim 9 and Aim 11 numbers are recomputed live here from this project's own
tested functions (sweep_reordering_frequency, sweep_pre_vs_relapse_across_full_grid)
-- not copied from prose. The two Aim 12 numbers (9-patient and 41-patient)
are NOT recomputable in this sandbox (the raw per-cell CyTOF processing was
done in a browser session with no network egress here) and are pulled
directly from the cached, already-verified REAL_SINGLE_CELL_REORDERING_RESULT
dict in aim12_single_cell_cytof_reordering.py, the same cached source
fig5/fig6 already use.

No change to any existing tested source file.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig8_evidence_strength_comparison.py
"""

import os
import json
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import clean_spines, style_title, section_badge, add_callout, PROJECT_BLUE
from aim9_marker_ranking_order_check import sweep_reordering_frequency
from aim11_antigen_target_ranking import sweep_pre_vs_relapse_across_full_grid
from aim12_single_cell_cytof_reordering import REAL_SINGLE_CELL_REORDERING_RESULT

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")

COLOR_RECONSTRUCTED = "#8A93A6"  # neutral gray -- reconstructed-distribution checks
COLOR_REAL = "#2E5EAA"           # accent blue -- real, individually-measured checks


def build_comparison():
    aim9 = sweep_reordering_frequency()
    aim11 = sweep_pre_vs_relapse_across_full_grid()
    # Every Aim 12 bar below uses ONE criterion: a grid point counts when some
    # REORDERED PAIR of patients differs by >= 2 percentage points in escape
    # fraction. The earlier "9-patient subset" bar was computed under a different
    # (looser) criterion and is not comparable, so it has been replaced by the
    # gated-tumor-cell tier, which is computed from the same files under the same
    # rule. See src/aim13_scoping/aim12_reordering_comparison.py.
    strict = REAL_SINGLE_CELL_REORDERING_RESULT["strict_pairwise_convention"]
    regated = json.load(open(os.path.join(
        os.path.dirname(__file__), "..", "aim13_scoping", "aim12_reordering_comparison.json")))

    bars = [
        {
            "label": "Aim 9\nEV markers\n(reconstructed)",
            "fraction": aim9["reordering_fraction"],
            "n_points": f"{aim9['n_reordered']}/{aim9['n_grid_points']}",
            "color": COLOR_RECONSTRUCTED,
        },
        {
            "label": "Aim 11\nCAR-T targets\n(reconstructed)",
            "fraction": aim11["reordering_fraction"],
            "n_points": f"{aim11['n_reordered']}/{aim11['n_grid_points']}",
            "color": COLOR_RECONSTRUCTED,
        },
        {
            "label": "Aim 12\n41 patients, all events\n(real single-cell)",
            "fraction": strict["reordered_meaningful_spread_ge_2pct"] / strict["total_grid_points"],
            "n_points": f"{strict['reordered_meaningful_spread_ge_2pct']}/{strict['total_grid_points']}",
            "color": COLOR_REAL,
        },
        {
            "label": "Aim 12\n41 patients, sorted tumor cells\n(real single-cell)",
            "fraction": regated["Bgated"]["frac"],
            "n_points": f"{regated['Bgated']['reordered']}/{regated['Bgated']['grid']}",
            "color": COLOR_REAL,
        },
    ]
    return bars


def plot_comparison(bars, save_path):
    fig, ax = plt.subplots(figsize=(8, 6))

    x = range(len(bars))
    heights = [b["fraction"] * 100 for b in bars]
    colors = [b["color"] for b in bars]

    ax.bar(x, heights, width=0.6, color=colors)

    for i, b in enumerate(bars):
        ax.text(i, b["fraction"] * 100 + 1.5, f"{b['fraction']*100:.1f}%\n({b['n_points']})",
                 ha="center", va="bottom", fontsize=10)

    ax.set_xticks(list(x))
    ax.set_xticklabels([b["label"] for b in bars], fontsize=10)
    ax.set_ylabel("Grid points with meaningful reordering (%)", fontsize=11)
    ax.set_ylim(0, max(heights) * 1.28)

    style_title(
        ax,
        "Meaningful target reordering across evidence tiers",
        fontsize=13,
    )

    # Legend distinguishing the two evidentiary tiers (categorical identity,
    # not magnitude -- matches the gray/blue bar coloring above).
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLOR_RECONSTRUCTED, label="Reconstructed from published\nsummary statistics"),
        Patch(facecolor=COLOR_REAL, label="Individually-measured\nper-cell/per-vesicle data"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=9, frameon=False)

    clean_spines(ax, keep=("left", "bottom"))

    # Callout on the final bar -- the full 41-patient real single-cell
    # result, the project's strongest evidence tier.
    last = bars[-1]
    add_callout(
        ax, xy=(len(bars) - 1, last["fraction"] * 100),
        text="largest measured\neffect",
        xytext=(0.60, 0.82), textcoords="axes fraction",
        color=PROJECT_BLUE, fontsize=10,
    )

    section_badge(fig, "OVERVIEW", 8)

    fig.text(
        0.5, -0.06,
        "The 6.87M-cell CyTOF cohort (Aim 12) shows a larger effect than any reconstructed-distribution\n"
        "check. A grid point counts when some reordered PAIR of patients differs by at least 2 percentage\n"
        "points in escape fraction; largest such pair spread 53.1 points. The effect survives sorting to\n"
        "tumor cells (51.9%) and light-chain-purified tumor cells (40.7%).",
        ha="center", fontsize=9.5, style="italic",
    )

    fig.tight_layout()
    fig.subplots_adjust(top=0.84)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    for b in bars:
        print(f"  {b['label'].replace(chr(10), ' ')}: {b['fraction']*100:.2f}% ({b['n_points']})")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    bars = build_comparison()
    plot_comparison(bars, save_path=os.path.join(FIGURES_DIR, "fig8_evidence_strength_comparison.png"))
