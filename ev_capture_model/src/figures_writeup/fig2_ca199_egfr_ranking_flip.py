"""
fig2_ca199_egfr_ranking_flip.py

Write-up figure 2 (Claim 1, EV foundation): the CA19-9 vs. EGFR ranking-
flip robustness comparison, on strict-resolved-only (pixel-exact-only)
digitized Lennon healthy/PDAC-plasma data. This is the figure version of
strict_resolved_egfr_ca199_comparison()'s headline result: EGFR's
tumor/normal yield ratio stays BELOW 1 (no flip) across the whole K_INTRA
plausibility range, while CA19-9's stays ABOVE 1 (flip holds) across the
same range -- the project's primary EV-side headline finding, with EGFR
retained only as a negative-control comparison.

Chart choice: a line plot (K_INTRA, log x-axis) vs. ratio_p2 (y-axis),
one line per marker, with a horizontal reference line at ratio=1 (the
ranking-flip threshold: above 1 means tumor yield exceeds normal yield).
This is the natural form for "does a ratio cross a threshold across an
ordered parameter sweep" -- not a bar chart (K_INTRA is continuous/ordered,
not categorical) and not a heatmap (only 2 series, not a 2D field).

Uses sweep_ranking_flip_robustness(strict_resolved_only=True) directly --
no new computation, just re-shaped for plotting.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig2_ca199_egfr_ranking_flip.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import clean_spines, style_title, section_badge, add_callout, PROJECT_BLUE
from aim3_aim4_range_sweep import sweep_ranking_flip_robustness

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")

# Fixed categorical color assignment (not cycled): CA19-9 is the primary
# headline finding, gets slot 1 (a confident blue); EGFR is the retired
# negative-control comparison, gets slot 2 (a muted orange/gray-warm) so
# it visually reads as "the one that didn't hold up," not as an equally-
# weighted second result.
COLOR_CA199 = "#2E5EAA"   # slot 1 -- primary finding
COLOR_EGFR = "#B0563A"    # slot 2 -- retired/negative-control


def plot_flip_comparison(save_path):
    results = sweep_ranking_flip_robustness(strict_resolved_only=True)

    k_ca199 = [row["k_intra"] for row in results["healthy_pdac_plasma_ca199"]]
    r_ca199 = [row["ratio_p2"] for row in results["healthy_pdac_plasma_ca199"]]
    k_egfr = [row["k_intra"] for row in results["healthy_pdac_plasma_egfr"]]
    r_egfr = [row["ratio_p2"] for row in results["healthy_pdac_plasma_egfr"]]

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.axhline(1.0, color="0.5", linewidth=1, zorder=1)
    ax.text(k_ca199[0], 1.0, "  ratio = 1  (no ranking flip)", va="bottom",
            ha="left", fontsize=8.5, color="0.35", style="italic")

    ax.plot(k_ca199, r_ca199, "-o", color=COLOR_CA199, linewidth=2,
            markersize=6, label="CA19-9 (flip holds)", zorder=3)
    ax.plot(k_egfr, r_egfr, "-o", color=COLOR_EGFR, linewidth=2,
            markersize=6, label="EGFR (flip does not survive)", zorder=3)

    ax.set_xscale("log")
    ax.set_xlabel(r"$K_{INTRA}$  (per-bond binding strength, cm$^2$/antibody)", fontsize=11)
    ax.set_ylabel("Tumor : normal yield ratio\n(strict pixel-exact points only)", fontsize=11)
    style_title(ax, "Tumor:normal yield ratio vs. $K_{INTRA}$ for CA19-9 and EGFR")
    clean_spines(ax, keep=("left", "bottom"))
    ax.legend(loc="center left", fontsize=10, frameon=False)
    ax.grid(axis="y", alpha=0.25, linewidth=0.8)

    ca199_min, ca199_max = min(r_ca199), max(r_ca199)
    egfr_min, egfr_max = min(r_egfr), max(r_egfr)

    # Callout on the K_INTRA value where CA19-9 and EGFR are furthest
    # apart -- the clearest single point demonstrating the ranking flip
    # that CA19-9 survives and EGFR does not.
    gaps = [rc - re for rc, re in zip(r_ca199, r_egfr)]
    peak_gap_idx = int(np.argmax(gaps))
    add_callout(
        ax, xy=(k_ca199[peak_gap_idx], r_ca199[peak_gap_idx]),
        text="largest CA19-9/EGFR\nseparation",
        xytext=(0.68, 0.88), textcoords="axes fraction",
        color=PROJECT_BLUE, fontsize=10,
    )

    section_badge(fig, "EV", 2)

    caption = (
        f"CA19-9 ratio stays in [{ca199_min:.2f}, {ca199_max:.2f}] (above 1 at every "
        f"tested $K_{{INTRA}}$); EGFR ratio stays in [{egfr_min:.2f}, {egfr_max:.2f}] "
        f"(below 1 at every tested $K_{{INTRA}}$). Healthy/PDAC-plasma comparison, "
        f"pixel-exact digitized points only (Lennon et al. 2019)."
    )
    fig.text(0.5, -0.05, caption, ha="center", fontsize=9.5, wrap=True)

    fig.tight_layout()
    fig.subplots_adjust(top=0.84)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    print(f"CA19-9 ratio_p2 range: [{ca199_min:.4f}, {ca199_max:.4f}]")
    print(f"EGFR ratio_p2 range: [{egfr_min:.4f}, {egfr_max:.4f}]")
    assert all(r > 1.0 for r in r_ca199), "CA19-9 should hold above 1 at every K_INTRA"
    assert all(r < 1.0 for r in r_egfr), "EGFR should stay below 1 at every K_INTRA"
    print("Sanity check passed: CA19-9 always >1, EGFR always <1.")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    plot_flip_comparison(save_path=os.path.join(FIGURES_DIR, "fig2_ca199_egfr_ranking_flip.png"))
