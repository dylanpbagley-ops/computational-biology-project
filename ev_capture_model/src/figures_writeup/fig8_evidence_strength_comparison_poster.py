"""
fig8_evidence_strength_comparison_poster.py

POSTER version of Figure 8 -- reuses build_comparison() from the paper
version directly (same live-recomputed Aim 9/11 numbers, same cached Aim
12 numbers), re-rendered bigger/bolder with a plain-language takeaway
banner replacing the technical caption. This is the poster's HERO figure.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig8_evidence_strength_comparison_poster.py
"""

import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import clean_spines, section_badge, add_callout, add_takeaway_banner, PROJECT_BLUE
from fig8_evidence_strength_comparison import build_comparison

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")


def build(save_path):
    bars = build_comparison()

    fig, ax = plt.subplots(figsize=(12.5, 9))
    x = range(len(bars))
    heights = [b["fraction"] * 100 for b in bars]
    colors = [b["color"] for b in bars]
    ax.bar(x, heights, width=0.62, color=colors)

    for i, b in enumerate(bars):
        ax.text(i, b["fraction"] * 100 + 1.8, f"{b['fraction']*100:.1f}%",
                 ha="center", va="bottom", fontsize=19, fontweight="bold")

    ax.set_xticks(list(x))
    ax.set_xticklabels([b["label"] for b in bars], fontsize=15)
    ax.set_ylabel("Targets meaningfully reordered (%)", fontsize=17)
    ax.tick_params(axis="y", labelsize=14)
    ax.set_ylim(0, max(heights) * 1.3)
    ax.set_title("Meaningful target reordering across evidence tiers",
                 fontsize=24, fontweight="bold", pad=16)

    legend_elements = [
        Patch(facecolor=bars[0]["color"], label="Estimated from published averages"),
        Patch(facecolor=bars[-1]["color"], label="Measured directly, patient by patient"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=13.5, frameon=False)
    clean_spines(ax, keep=("left", "bottom"))
    section_badge(fig, "OVERVIEW", 8)

    add_callout(
        ax, xy=(3, bars[3]["fraction"] * 100), text="largest measured\neffect",
        xytext=(0.62, 0.86), textcoords="axes fraction",
        color=PROJECT_BLUE, fontsize=15,
    )

    fig.tight_layout()
    fig.subplots_adjust(top=0.84)  # leave room for the section badge above the title
    add_takeaway_banner(
        fig,
        "Single-cell data (41 patients, 6.87M cells) shows several times more target reordering "
        "than any estimate built from summary statistics, and the effect persists after sorting "
        "to tumor cells.",
        color=PROJECT_BLUE, y=-0.08,
    )
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    for b in bars:
        print(f"  {b['label'].replace(chr(10), ' ')}: {b['fraction']*100:.2f}%")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    build(os.path.join(FIGURES_DIR, "fig8_evidence_strength_comparison_poster.png"))
