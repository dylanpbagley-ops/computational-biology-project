"""
fig5_aim12_cytof_reordering_example_poster.py

POSTER version of Figure 5 -- same cached REAL_SINGLE_CELL_REORDERING_RESULT
data as the paper version, re-rendered bigger/bolder with a plain-language
takeaway banner replacing the technical caption.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig5_aim12_cytof_reordering_example_poster.py
"""

import os
import sys
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import clean_spines, section_badge, add_takeaway_banner, PROJECT_BLUE, PROJECT_ORANGE, NEUTRAL_GRAY
import json
GATED_JSON = os.path.join(os.path.dirname(__file__), '..', 'aim13_scoping', 'gated_reordering_example.json')

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")


def _parse(entry):
    m = re.match(r"([A-Z0-9]+)\s*\(.*?=([\d.]+)%?\)", entry)
    return m.group(1), float(m.group(2))


def build(save_path):
    ex = json.load(open(GATED_JSON))
    naive = [(pid, float(v)) for pid, v in ex["naive_top5"]]
    model = [(pid, 100*float(v)) for pid, v in ex["model_top5"]]
    naive_ids = [p[0] for p in naive]
    model_ids = [p[0] for p in model]
    displaced = set(naive_ids) - set(model_ids)
    promoted = set(model_ids) - set(naive_ids)

    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-0.2, 6.6)

    for i, (pid, val) in enumerate(naive):
        y = 5.5 - i
        color = PROJECT_ORANGE if pid in displaced else NEUTRAL_GRAY
        ax.scatter(0, y, s=420, color=color, zorder=3)
        ax.text(-0.08, y, f"{pid}", ha="right", va="center", fontsize=17, fontweight="bold")
    for i, (pid, val) in enumerate(model):
        y = 5.5 - i
        color = PROJECT_BLUE if pid in promoted else NEUTRAL_GRAY
        ax.scatter(1, y, s=420, color=color, zorder=3)
        ax.text(1.08, y, f"{pid}", ha="left", va="center", fontsize=17, fontweight="bold")

    naive_pos = {pid: 5.5 - i for i, (pid, _) in enumerate(naive)}
    model_pos = {pid: 5.5 - i for i, (pid, _) in enumerate(model)}
    for pid in naive_ids:
        if pid in model_pos:
            ax.plot([0, 1], [naive_pos[pid], model_pos[pid]], color=NEUTRAL_GRAY, linewidth=2.2, zorder=1)
    for pid in displaced:
        ax.plot([0, 0.35], [naive_pos[pid], naive_pos[pid] - 0.9], "--", color=PROJECT_ORANGE, linewidth=3, zorder=2)
    for pid in promoted:
        ax.plot([0.65, 1], [model_pos[pid] + 0.9, model_pos[pid]], "--", color=PROJECT_BLUE, linewidth=3, zorder=2)

    ax.text(0, 6.3, "Ranked by\naverage density", ha="center", va="bottom", fontsize=16, fontweight="bold")
    ax.text(1, 6.3, "Ranked by\nescape fraction", ha="center", va="bottom", fontsize=16, fontweight="bold")

    ax.set_title("Naive vs. model-based top-5 CAR-T target ranking",
                 fontsize=22, fontweight="bold", pad=20)
    ax.axis("off")
    clean_spines(ax, keep=())
    section_badge(fig, "CAR-T", 5)

    fig.tight_layout()
    fig.subplots_adjust(top=0.84)  # leave room for the section badge above the title
    add_takeaway_banner(
        fig,
        "Two patients enter the model's top 5 while two with HIGHER average antigen density "
        "drop out. Using the full per-cell distribution instead of the average changes the "
        "ranking. Computed on sorted tumor cells, not whole biopsies.",
        color=PROJECT_ORANGE, y=-0.06,
    )
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    build(os.path.join(FIGURES_DIR, "fig5_aim12_cytof_reordering_example_poster.png"))
