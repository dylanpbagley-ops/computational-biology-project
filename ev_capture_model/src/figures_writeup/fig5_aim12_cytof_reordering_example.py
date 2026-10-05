"""
fig5_aim12_cytof_reordering_example.py

Write-up figure 5 (Claim 2, CAR-T pivot, the project's strongest result):
a concrete, patient-level illustration of the Aim 12 reordering effect --
the real 41-patient CyTOF DLBCL cohort's naive top-5 (by mean CD19
density) vs. the distribution-aware model's real top-5 (by escape
fraction), at the CLN-978_cytolysis calibration, kill_threshold=0.5.

This figure exists specifically because a grid-wide reordering percentage
is an abstract, hard-to-visualize summary number -- this shows concrete
instances of what a reordering actually looks like.

REBUILT ON GATED TUMOR CELLS (September 7, 2026). The earlier version used
per-patient means computed over EVERY event in each file -- whole biopsies,
which are 5%-87% B cells. Its worked example (P7676 mean 109.1 entering the
top 5 and displacing Q7915 mean 111.8) inverts once cells are gated to B
cells: Q7915 is then 230.3 and P7676 is 139.0. The figure now uses
QC-gated, CD3-negative, CD19/CD20-positive B cells, computed live from the
local FCS files rather than from a cached constant. The gated example is
cleaner: TWO patients enter the top 5 (O8155 mean 176.9, Q5971 mean 178.8)
displacing TWO with HIGHER means (Q7915 mean 230.3, R5097 mean 205.2), so
the point -- lower mean, better predicted outcome -- holds unambiguously
for every promoted/displaced pair rather than for one pair only.

Chart choice: paired dumbbell/slope chart (2 columns: naive rank vs.
model rank), the standard form for "before -> after per item".
Data source: src/aim13_scoping/gated_reordering_example.json, produced by
gated_reordering_example.py from the real local FCS files.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig5_aim12_cytof_reordering_example.py
"""

import os
import sys
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import json
from figure_style import clean_spines, section_badge, TITLE_WEIGHT

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")

COLOR_STABLE = "#5B6472"      # patients that appear in both top-5 lists
COLOR_DISPLACED = "#B0563A"   # drops out of the top 5 despite a higher mean
COLOR_PROMOTED = "#2E5EAA"    # enters the top 5 despite a lower mean

GATED_JSON = os.path.join(os.path.dirname(__file__), "..", "aim13_scoping",
                          "gated_reordering_example.json")


def _parse_entry(entry):
    """'O6357 (mean=348.5)' or 'O6357 (escape=19.1%)' -> (id, value_str)."""
    m = re.match(r"(\S+)\s+\((?:mean|escape)=([^)]+)\)", entry)
    return m.group(1), m.group(2)


def plot_reordering_example(save_path):
    with open(GATED_JSON) as fh:
        ex = json.load(fh)
    naive = [(pid, f"{v:.1f}") for pid, v in ex["naive_top5"]]
    model = [(pid, f"{100*v:.1f}%") for pid, v in ex["model_top5"]]

    naive_ids = [pid for pid, _ in naive]
    model_ids = [pid for pid, _ in model]

    fig, ax = plt.subplots(figsize=(7.5, 6))

    y_naive = {pid: 5.5 - i for i, pid in enumerate(naive_ids)}
    y_model = {pid: 5.5 - i for i, pid in enumerate(model_ids)}

    all_ids = sorted(set(naive_ids) | set(model_ids))
    for pid in all_ids:
        y0 = y_naive.get(pid)
        y1 = y_model.get(pid)
        if y0 is not None and y1 is not None:
            color = COLOR_STABLE
            lw = 1.5
        elif y0 is not None and y1 is None:
            color = COLOR_DISPLACED
            lw = 2.5
        else:
            color = COLOR_PROMOTED
            lw = 2.5
        # Stable patients (appear in both top-5 lists) get a solid connecting
        # line. Displaced/promoted patients do NOT get a line at all -- an
        # earlier version used a short dashed "stub" trailing off each one,
        # but that was genuinely ambiguous (it read as if the two stubs were
        # connected to each other). A small colored badge next to each dot
        # (added in the naive/model plotting loops below) replaces it.
        if y0 is not None and y1 is not None:
            ax.plot([0, 1], [y0, y1], "-", color=color, linewidth=lw, zorder=2)

    badge_style = dict(fontsize=8.5, color="white", fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.3", edgecolor="none"))

    for pid, val in naive:
        y = y_naive[pid]
        is_displaced = pid in naive_ids and pid not in model_ids
        color = COLOR_DISPLACED if is_displaced else COLOR_STABLE
        ax.scatter(0, y, s=90, color=color, zorder=3)
        ax.text(-0.06, y, f"{pid}  (mean={val})", ha="right", va="center", fontsize=9.5)
        if is_displaced:
            # Small badge instead of a trailing dashed line (which read as
            # ambiguously connected to the "promoted" patient's own stub).
            ax.text(0.10, y - 0.34, "drops out of top 5", ha="left", va="center",
                     bbox=dict(boxstyle="round,pad=0.3", facecolor=COLOR_DISPLACED, edgecolor="none"),
                     **{k: v for k, v in badge_style.items() if k != "bbox"})

    for pid, val in model:
        y = y_model[pid]
        is_promoted = pid in model_ids and pid not in naive_ids
        color = COLOR_PROMOTED if is_promoted else COLOR_STABLE
        ax.scatter(1, y, s=90, color=color, zorder=3)
        ax.text(1.06, y, f"{pid}  (escape={val})", ha="left", va="center", fontsize=9.5)
        if is_promoted:
            ax.text(0.90, y + 0.34, "new entry to top 5", ha="right", va="center",
                     bbox=dict(boxstyle="round,pad=0.3", facecolor=COLOR_PROMOTED, edgecolor="none"),
                     **{k: v for k, v in badge_style.items() if k != "bbox"})

    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-0.2, 6.6)

    # Column headers as text annotations ABOVE the topmost dot (not as
    # x-axis tick labels), so they can never collide with the caption
    # below -- an earlier version used set_xticklabels() at y=0 and the
    # labels landed on top of the caption text.
    ax.text(0, 6.3, "Naive ranking\n(mean CD19 density)", ha="center", va="bottom", fontsize=11)
    ax.text(1, 6.3, "Model ranking\n(escape fraction)", ha="center", va="bottom", fontsize=11)

    ax.set_xticks([])
    ax.set_yticks([])
    clean_spines(ax, keep=())

    # Not using style_title() here: this figure needs a much larger pad
    # (55, vs. the shared default of 12) so the title clears the "Naive
    # ranking" / "Model ranking" column headers placed just above the
    # topmost dot -- so the bold weight is applied directly via the same
    # TITLE_WEIGHT constant instead, keeping the same look with the pad
    # this layout actually needs.
    ax.set_title(
        "Naive vs. model-based top-5 CAR-T target ranking\n"
        "(mean CD19 density vs. escape fraction)",
        fontsize=13, fontweight=TITLE_WEIGHT, pad=55,
    )

    caption = (
        "41-patient DLBCL CyTOF cohort, gated tumor B cells, CLN-978_cytolysis calibration, "
        "kill_threshold=0.5 (antigen-density threshold n* = 81 molecules/cell). O8155 (mean 176.9) "
        "and Q5971 (mean 178.8) enter the model's top 5 by escape fraction, displacing Q7915 "
        "(mean 230.3) and R5097 (mean 205.2), both of which have higher averages. R5097 falls from "
        "4th to 14th. One instance of the grid-wide reordering finding (48.1% of grid points)."
    )
    fig.text(0.5, -0.06, caption, ha="center", fontsize=9.2, wrap=True)

    section_badge(fig, "CAR-T", 5)
    # No add_callout() here: the displaced ("drops out of top 5", orange)
    # and promoted ("new entry to top 5", blue) patients already carry their
    # own small colored badges, plus color-coded dots distinct from the
    # gray "stable" ones -- the whole chart IS the callout on this figure's
    # key feature, so a boxed arrow on top of it would be redundant.

    fig.tight_layout()
    fig.subplots_adjust(top=0.72)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    print(f"Naive top 5: {naive_ids}")
    print(f"Model top 5: {model_ids}")
    displaced = set(naive_ids) - set(model_ids)
    promoted = set(model_ids) - set(naive_ids)
    print(f"Displaced: {displaced}, Promoted: {promoted}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    plot_reordering_example(save_path=os.path.join(FIGURES_DIR, "fig5_aim12_cytof_reordering_example.png"))
