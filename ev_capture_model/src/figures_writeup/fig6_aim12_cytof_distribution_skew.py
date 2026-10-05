"""
fig6_aim12_cytof_distribution_skew.py

Write-up figure 6 (Claim 2, CAR-T pivot): WHY the Aim 12 reordering effect
is so much larger than every reconstructed-distribution check this
project has run. The project's own documented mechanistic explanation:
real per-patient mean/median ratios range 1.08x-29.78x (a 27.5x spread in
skewness) and near-zero-cell fractions range 0.28%-21.68% (a 77x spread)
across the real 41-patient cohort -- variation a 2-parameter
(mean, variance) reconstruction from summary statistics cannot represent,
but which a distribution-SHAPE-sensitive model is specifically able to
respond to.

Chart choice: a scatter plot of mean vs. median CD19 density, one point
per real patient (41 total, comfortably a "many items, show the
relationship" job -- not a categorical/identity chart, so a single
sequential-toned scatter is correct, not 41 distinct hues). A patient
sitting far from the diagonal (mean >> median) is a strongly right-skewed
population -- exactly the shape information a naive mean-only ranking
throws away. The diagonal reference line marks "mean = median" (a
symmetric distribution); distance from it is the story.

Uses ONLY REAL_PATIENT_CD19_SUMMARY, the cached real per-patient
statistics already reported in aim12_single_cell_cytof_reordering.py --
no new computation.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig6_aim12_cytof_distribution_skew.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import clean_spines, style_title, section_badge
from aim12_single_cell_cytof_reordering import REAL_PATIENT_CD19_SUMMARY

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")


def plot_skew_scatter(save_path):
    patients = list(REAL_PATIENT_CD19_SUMMARY.keys())
    means = np.array([REAL_PATIENT_CD19_SUMMARY[p]["mean"] for p in patients])
    medians = np.array([REAL_PATIENT_CD19_SUMMARY[p]["median"] for p in patients])
    zero_frac = np.array([REAL_PATIENT_CD19_SUMMARY[p]["zero_frac"] for p in patients])

    # Avoid division by (near-)zero median for the ratio metric; only used
    # for identifying the extremes to label, not for the plotted axes.
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(medians > 0, means / medians, np.nan)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))

    lo = min(means.min(), medians.min()) * 0.7
    hi = max(means.max(), medians.max()) * 1.15
    ax.plot([lo, hi], [lo, hi], color="0.6", linewidth=1.2, zorder=1)
    ax.text(hi * 0.022, hi * 0.60, "mean = median line\n(symmetric distribution)",
            fontsize=8.5, color="0.4", style="italic", ha="left")

    sc = ax.scatter(medians, means, s=55, c=zero_frac, cmap="viridis",
                     edgecolors="white", linewidths=0.6, zorder=3)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Fraction of near-zero-CD19 cells", fontsize=10.5)

    # Label the two most extreme skew ratios, since 41 direct labels would
    # be chartjunk -- selective labeling per the marks-and-anatomy rule.
    valid = ~np.isnan(ratio)
    top2_idx = np.argsort(np.where(valid, ratio, -np.inf))[-2:]
    for idx in top2_idx:
        pid = patients[idx]
        ax.annotate(
            f"{pid}\n({ratio[idx]:.1f}$\\times$ mean/median)",
            xy=(medians[idx], means[idx]),
            xytext=(12, 8), textcoords="offset points", fontsize=8.5,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Median per-cell CD19 signal (real, per patient)", fontsize=11)
    ax.set_ylabel("Mean per-cell CD19 signal (real, per patient)", fontsize=11)
    style_title(
        ax,
        "Mean vs. median per-cell CD19 signal across 41 patients",
        fontsize=12.5,
    )
    clean_spines(ax, keep=("left", "bottom"))

    finite_ratio = ratio[valid]
    caption = (
        f"41 patients, complete DLBCL CyTOF cohort. Mean/median ratio spans "
        f"{finite_ratio.min():.2f}$\\times$ to {finite_ratio.max():.2f}$\\times$ "
        f"({finite_ratio.max()/finite_ratio.min():.1f}$\\times$ spread in skewness); "
        f"near-zero-cell fraction spans {zero_frac.min()*100:.2f}% to {zero_frac.max()*100:.2f}% "
        f"({zero_frac.max()/zero_frac.min():.0f}$\\times$ spread). This distributional "
        f"heterogeneity is discarded by a mean-only ranking, consistent with directly-measured "
        f"single-cell data showing a larger reordering effect (48.1% of grid points) than "
        f"reconstructed-distribution checks elsewhere in this project (5.6%)."
    )
    fig.text(0.5, -0.1, caption, ha="center", fontsize=9, wrap=True)

    section_badge(fig, "CAR-T", 6)
    # No add_callout() here: the two most extreme-skew patients are already
    # individually labeled with their own annotations above -- an
    # add_callout() box would just duplicate that.

    fig.tight_layout()
    fig.subplots_adjust(top=0.84)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    print(f"Mean/median ratio range: {finite_ratio.min():.4f} - {finite_ratio.max():.4f}")
    print(f"Zero-fraction range: {zero_frac.min():.4f} - {zero_frac.max():.4f}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    plot_skew_scatter(save_path=os.path.join(FIGURES_DIR, "fig6_aim12_cytof_distribution_skew.png"))
