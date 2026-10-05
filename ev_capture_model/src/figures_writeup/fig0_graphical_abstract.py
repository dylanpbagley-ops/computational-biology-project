"""
fig0_graphical_abstract.py

A drawn overview/graphical-abstract figure (icons, arrows, flow -- NOT a
data plot) meant to orient a viewer before they hit the 8 technical
figures. Built after Dylan flagged that the figure set "didn't look like
a winning ISEF project's" and design research confirmed a real, specific
gap: no overview figure existed, and graphical abstracts (icons +
arrows + simple compare/contrast, letting a viewer "understand findings
before they actually understand what's going on") are a standard,
well-documented part of how strong scientific communication -- and
strong posters -- work. See BioRender's own description of graphical
abstracts, and EV-Capture-CURRENT.md's September 3, 2026 entries for the
full sourcing.

Where this figure uses real project content (not just decoration):
  - The little "model" curve in stage 2 IS binding_model.capture_probability(),
    the actual tested function every other figure in this project is built
    on -- not a fake illustrative curve.
  - The per-cell-density icon in stage 1 is sampled from the real cached
    REAL_PATIENT_CD19_SUMMARY statistics (a real patient's mean/zero-frac),
    not arbitrary numbers.
  - The headline number in the takeaway banner (48.1%, 13/27, the strict
    pairwise criterion) is pulled
    directly from REAL_SINGLE_CELL_REORDERING_RESULT, the same cached,
    already-verified source fig5/fig6/fig8 use.

No change to any existing tested source file.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig0_graphical_abstract.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import PROJECT_BLUE, PROJECT_ORANGE, NEUTRAL_GRAY, INK, section_badge
from binding_model import capture_probability
from aim12_single_cell_cytof_reordering import (
    REAL_PATIENT_CD19_SUMMARY,
    REAL_SINGLE_CELL_REORDERING_RESULT,
)

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")


def draw_stage_box(ax, x0, x1, y0, y1, color, title):
    box = FancyBboxPatch(
        (x0, y0), x1 - x0, y1 - y0,
        boxstyle="round,pad=0.04,rounding_size=0.12",
        linewidth=2.2, edgecolor=color, facecolor="white", zorder=1,
    )
    ax.add_patch(box)
    ax.text((x0 + x1) / 2, y1 - 0.35, title, ha="center", va="top",
             fontsize=13.5, fontweight="bold", color=color, zorder=3)


def draw_arrow(ax, x0, x1, y, color):
    arrow = FancyArrowPatch(
        (x0, y), (x1, y), arrowstyle="-|>", mutation_scale=28,
        linewidth=3.5, color=color, zorder=2,
    )
    ax.add_patch(arrow)


def stage1_cell_icon(ax, cx, cy, w, h):
    """Cluster of circles whose sizes/shades are sampled from REAL patient
    CD19 summary stats (log-scaled means across the real 41-patient
    cohort) -- illustrating real per-patient variation, not arbitrary
    decoration."""
    means = np.array([v["mean"] for v in REAL_PATIENT_CD19_SUMMARY.values()])
    rng = np.random.default_rng(0)
    sample = rng.choice(means, size=26, replace=False)
    log_sample = np.log10(sample)
    norm = (log_sample - log_sample.min()) / (log_sample.max() - log_sample.min())

    n_cols, n_rows = 6, 5
    xs = np.linspace(cx - w * 0.38, cx + w * 0.38, n_cols)
    ys = np.linspace(cy + h * 0.05, cy - h * 0.32, n_rows)
    k = 0
    for yi in ys:
        for xi in xs:
            if k >= len(norm):
                break
            r = 0.035 + 0.075 * norm[k]
            shade = plt.cm.viridis(0.15 + 0.7 * norm[k])
            ax.add_patch(Circle((xi + rng.uniform(-0.03, 0.03),
                                  yi + rng.uniform(-0.03, 0.03)),
                                 r, facecolor=shade, edgecolor="white",
                                 linewidth=0.5, zorder=3))
            k += 1


def stage2_model_icon(ax, cx, cy, w, h):
    """The real capture_probability() curve -- the actual function every
    other figure in this project is built on -- as a small inset sparkline,
    with a linear reference for contrast (same visual grammar as fig7)."""
    inset = ax.inset_axes([
        (cx - w * 0.36 - (ax.get_xlim()[0])) / (ax.get_xlim()[1] - ax.get_xlim()[0]),
        (cy - h * 0.22 - (ax.get_ylim()[0])) / (ax.get_ylim()[1] - ax.get_ylim()[0]),
        (w * 0.72) / (ax.get_xlim()[1] - ax.get_xlim()[0]),
        (h * 0.42) / (ax.get_ylim()[1] - ax.get_ylim()[0]),
    ])
    n = np.linspace(0, 100, 200)
    theta = capture_probability(n, sigma=3.5e10, k_intra=1e-11, z=1e-6)  # chosen to visibly saturate within n=[0,100], same qualitative shape as fig7
    inset.plot(n, theta, color=PROJECT_ORANGE, linewidth=2.8)
    inset.plot([0, 100], [0, theta.max()], "--", color=NEUTRAL_GRAY, linewidth=1.3)
    inset.set_xticks([]); inset.set_yticks([])
    for spine in inset.spines.values():
        spine.set_visible(False)


def stage3_ranking_icon(ax, cx, cy, w, h):
    """Two tiny 5-row rank columns -- naive vs. model -- with ONE
    highlighted swap, the same visual grammar as fig5's dumbbell chart,
    shrunk to icon scale."""
    left_x = cx - w * 0.30
    right_x = cx + w * 0.30
    top_y = cy - h * 0.02
    row_h = h * 0.11
    bar_w = w * 0.20
    bar_h = row_h * 0.55

    naive_order = [0, 1, 2, 3, 4]
    model_order = [0, 2, 3, 1, 5]  # index 4 (bottom-naive) swapped out for a new entrant (5)
    colors = [NEUTRAL_GRAY] * 4 + [PROJECT_ORANGE]

    for i, idx in enumerate(naive_order):
        y = top_y - i * row_h
        ax.add_patch(plt.Rectangle((left_x - bar_w / 2, y - bar_h / 2), bar_w, bar_h,
                                    facecolor=NEUTRAL_GRAY if i < 4 else PROJECT_ORANGE,
                                    edgecolor="none", zorder=3))
    for i, idx in enumerate(model_order):
        y = top_y - i * row_h
        color = PROJECT_BLUE if idx == 5 else NEUTRAL_GRAY
        ax.add_patch(plt.Rectangle((right_x - bar_w / 2, y - bar_h / 2), bar_w, bar_h,
                                    facecolor=color, edgecolor="none", zorder=3))

    # Connector for the displaced->promoted swap.
    ax.add_patch(FancyArrowPatch(
        (left_x + bar_w / 2, top_y - 4 * row_h),
        (right_x - bar_w / 2, top_y - 4 * row_h),
        arrowstyle="-|>", mutation_scale=10, linewidth=1.4,
        color=PROJECT_ORANGE, linestyle=(0, (3, 2)), zorder=2, alpha=0.7,
    ))
    ax.text(left_x, top_y - 5 * row_h - 0.12, "naive\n(mean)", ha="center",
            va="top", fontsize=9, color=INK)
    ax.text(right_x, top_y - 5 * row_h - 0.12, "model\n(distribution)", ha="center",
            va="top", fontsize=9, color=INK)


def build(save_path):
    fig, ax = plt.subplots(figsize=(15, 6.2))
    ax.set_xlim(0, 15.5)
    ax.set_ylim(0, 6.2)
    ax.axis("off")

    stage_y0, stage_y1 = 0.9, 5.35
    boxes = [
        (0.4, 4.6, "SINGLE-CELL DATA", PROJECT_BLUE),
        (5.6, 9.8, "NONLINEAR BINDING MODEL", PROJECT_ORANGE),
        (10.8, 15.1, "REORDERED RANKING", PROJECT_BLUE),
    ]
    for x0, x1, title, color in boxes:
        draw_stage_box(ax, x0, x1, stage_y0, stage_y1, color, title)

    cy = (stage_y0 + stage_y1) / 2 - 0.15
    stage1_cell_icon(ax, (0.4 + 4.6) / 2, cy, 4.2, stage_y1 - stage_y0)
    stage2_model_icon(ax, (5.6 + 9.8) / 2, cy, 4.2, stage_y1 - stage_y0)
    stage3_ranking_icon(ax, (10.8 + 15.1) / 2, cy, 4.3, stage_y1 - stage_y0)

    draw_arrow(ax, 4.7, 5.5, (stage_y0 + stage_y1) / 2, INK)
    draw_arrow(ax, 9.9, 10.7, (stage_y0 + stage_y1) / 2, INK)

    ax.text((0.4 + 4.6) / 2, stage_y0 - 0.25, "41 patients · 6.87M cells (CyTOF)",
            ha="center", va="top", fontsize=10.5, color=INK)
    ax.text((5.6 + 9.8) / 2, stage_y0 - 0.25, "distribution-aware escape fraction",
            ha="center", va="top", fontsize=10.5, color=INK)
    ax.text((10.8 + 15.1) / 2, stage_y0 - 0.25, "not just the average density",
            ha="center", va="top", fontsize=10.5, color=INK)

    ax.text(7.75, 5.85, "Overview: from single-cell data to a reordered CAR-T target ranking",
            ha="center", va="bottom", fontsize=17, fontweight="bold", color=INK)

    strict = REAL_SINGLE_CELL_REORDERING_RESULT["strict_pairwise_convention"]
    n_reordered = strict["reordered_meaningful_spread_ge_2pct"]
    n_total = strict["total_grid_points"]
    frac = n_reordered / n_total
    takeaway = (
        f"Across 41 DLBCL patients, distribution-aware scoring changes the CAR-T target ranking "
        f"in {frac*100:.1f}% ({n_reordered}/{n_total}) of tested conditions (Figure 8)."
    )
    fig.text(0.5, -0.06, takeaway, ha="center", fontsize=11.5, fontweight="bold",
              color="white",
              bbox=dict(boxstyle="round,pad=0.6", facecolor=PROJECT_BLUE, edgecolor="none"))

    section_badge(fig, "OVERVIEW", 0)

    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    print(f"Headline fraction used in banner: {frac*100:.1f}% ({n_reordered}/{n_total})")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    build(os.path.join(FIGURES_DIR, "fig0_graphical_abstract.png"))
