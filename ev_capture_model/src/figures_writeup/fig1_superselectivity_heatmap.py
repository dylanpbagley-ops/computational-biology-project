"""
fig1_superselectivity_heatmap.py

Write-up figure 1 (Claim 1, EV foundation): heatmap of the maximum
selectivity alpha reached across the full (K_INTRA, z) plausibility grid,
for the exact closed-form alpha_exact(). This is the figure version of
sweep_superselectivity_reachability()'s headline number ("alpha>1 reached
in ~59.6% of the numerically-reliable grid, peaking at alpha~10.46") --
that function only returns the single global max; this script recomputes
the SAME per-cell values it computes internally, but keeps the full 2D
grid instead of collapsing it to one number, purely for plotting. No
change to any existing tested source file.

REVISED after a design review (see docs -- the first version used
imshow(), which draws every grid cell at equal size regardless of its
actual K_INTRA value. The 8 K_INTRA values are NOT evenly log-spaced (they
come from a small n x target_exponent cross-product, not a clean
np.logspace sweep -- confirmed directly: 2.25e-16, 4.50e-16, 2.25e-15,
4.50e-15, 2.25e-14, 4.50e-14, 2.25e-13, 2.25e-12), so equal-sized cells
visually implied an even sweep that isn't real -- a real distortion, not
just a style choice. Fixed by using pcolormesh with true log-spaced cell
EDGES (computed as geometric midpoints between consecutive values), so
each cell's size on the page reflects its real spacing on the K_INTRA
axis. The z axis (13 points) genuinely is an even np.logspace, so it
didn't have this problem, but is now drawn the same correct way for
consistency.

Run from ev_capture_model/src/figures_writeup/:
    python3 fig1_superselectivity_heatmap.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from figure_style import clean_spines, style_title, section_badge, add_callout, PROJECT_BLUE
from binding_model import selectivity_alpha, alpha_exact
from aim3_aim4_range_sweep import (
    get_k_intra_sweep_values,
    get_z_sweep_values,
    SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
    MIN_RELIABLE_EXPONENT_AT_N1,
    ALPHA_SUPERSELECTIVE_THRESHOLD,
)

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "Poster, Paper, and Write Up", "Figures")


def build_alpha_grid(sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2, n_max=200):
    """
    Reproduce sweep_superselectivity_reachability()'s internal per-cell
    computation, but return the full 2D grid of max-alpha-per-cell values
    (exact, closed-form) instead of collapsing to one global max. Mirrors
    that function's logic exactly (same exclusion rule, same n grid) so
    the figure is a faithful visualization of the same tested result, not
    a new computation.
    """
    k_intra_values = get_k_intra_sweep_values()
    z_values = get_z_sweep_values()
    n_values = np.linspace(1, n_max, 60)

    grid = np.full((len(k_intra_values), len(z_values)), np.nan)
    excluded = np.zeros((len(k_intra_values), len(z_values)), dtype=bool)

    for i, k_intra in enumerate(k_intra_values):
        exponent_at_n1 = 1 * sigma * k_intra
        if exponent_at_n1 < MIN_RELIABLE_EXPONENT_AT_N1:
            excluded[i, :] = True
            continue
        for j, z in enumerate(z_values):
            try:
                n_mid, _ = selectivity_alpha(n_values, sigma, k_intra, z)
                alpha_ex = alpha_exact(n_mid, sigma, k_intra, z)
            except (ValueError, FloatingPointError):
                continue
            reliable = alpha_ex[np.isfinite(alpha_ex)]
            if len(reliable) > 0:
                grid[i, j] = reliable.max()

    return k_intra_values, z_values, grid, excluded


def _log_edges(values):
    """
    Given N sorted, positive values (possibly unevenly log-spaced), return
    N+1 cell EDGES in log space: each interior edge is the geometric mean
    of its two neighbors, and the two outer edges extrapolate the same
    half-step used at the nearest interior edge. This is what makes
    pcolormesh draw each cell at a width that reflects its REAL spacing
    on a log axis, instead of imshow's equal-size-per-cell distortion.
    """
    log_v = np.log10(values)
    interior = (log_v[:-1] + log_v[1:]) / 2.0
    first_half_step = interior[0] - log_v[0]
    last_half_step = log_v[-1] - interior[-1]
    edges = np.concatenate((
        [log_v[0] - first_half_step],
        interior,
        [log_v[-1] + last_half_step],
    ))
    return 10 ** edges


def plot_heatmap(k_intra_values, z_values, grid, excluded, save_path):
    fig, ax = plt.subplots(figsize=(8.5, 6))

    z_edges = _log_edges(z_values)
    k_edges = _log_edges(k_intra_values)

    masked_grid = np.ma.masked_invalid(grid)

    im = ax.pcolormesh(
        z_edges, k_edges, masked_grid,
        cmap="viridis", shading="flat",
        edgecolors="white", linewidth=0.6,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Mark the actual sampled (z, K_INTRA) points so the reader can see
    # these are 104 discrete tested combinations, not a continuous field.
    Zg, Kg = np.meshgrid(z_values, k_intra_values)
    ax.scatter(Zg, Kg, s=4, color="white", alpha=0.5, linewidths=0, zorder=3)

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(r"Maximum selectivity $\alpha_{exact}$ reached", fontsize=11)

    ax.set_xlabel(r"$z$  (vesicle concentration / transport activity)", fontsize=11)
    ax.set_ylabel(r"$K_{INTRA}$  (per-bond binding strength, cm$^2$/antibody)", fontsize=11)
    style_title(ax, "Selectivity $\\alpha_{exact}$ across the physically-plausible\n"
                "$(K_{INTRA}, z)$ parameter grid")
    clean_spines(ax, keep=("left", "bottom"))
    section_badge(fig, "EV", 1)

    n_reliable = int(np.sum(~np.isnan(grid)))
    n_superselective = int(np.sum(grid[~np.isnan(grid)] > ALPHA_SUPERSELECTIVE_THRESHOLD))
    frac = n_superselective / n_reliable if n_reliable else 0.0
    peak = float(np.nanmax(grid))

    # Callout on the peak cell -- the ISEF poster guide's own "box or arrow
    # highlighting important data" recommendation; the original version
    # relied entirely on the caption text below to convey this.
    peak_idx = np.unravel_index(np.nanargmax(grid), grid.shape)
    peak_k, peak_z = k_intra_values[peak_idx[0]], z_values[peak_idx[1]]
    add_callout(
        ax, xy=(peak_z, peak_k), text=f"peak \u03b1={peak:.2f}",
        xytext=(0.30, 0.72), textcoords="axes fraction",
        color=PROJECT_BLUE, fontsize=10,
    )

    caption = (
        f"{frac*100:.1f}% of the {n_reliable} numerically-reliable grid points "
        f"reach $\\alpha$ > 1 (superselectivity); peak $\\alpha_{{exact}}$ = {peak:.2f}. "
        f"White dots mark the {n_reliable} tested (K_INTRA, z) combinations."
    )
    fig.text(0.5, -0.04, caption, ha="center", fontsize=9.5, wrap=True)

    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    print(f"Reliable grid cells: {n_reliable} / {grid.size}")
    print(f"Superselective fraction: {frac:.4f} ({n_superselective}/{n_reliable})")
    print(f"Peak alpha_exact: {peak:.4f}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    k_intra_values, z_values, grid, excluded = build_alpha_grid()
    plot_heatmap(
        k_intra_values, z_values, grid, excluded,
        save_path=os.path.join(FIGURES_DIR, "fig1_superselectivity_heatmap.png"),
    )
