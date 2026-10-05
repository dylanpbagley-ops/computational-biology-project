"""
figure_style.py

Shared visual-identity module for every write-up figure (paper + poster).
Built after a direct design review (see EV-Capture-CURRENT.md, September 3
2026 entry) found that reusing plain journal-style matplotlib figures
one-for-one as poster figures was a real gap, not just a feeling -- sourced
against real poster-design guidance:
  - "Recreate figures specifically for your poster" (not paper figures
    reused as-is) -- https://conceptviz.app/blog/how-to-design-award-winning-scientific-poster
  - Poster figures need a unified color/font/line-weight system across the
    whole set -- same source.
  - Add boxes/arrows highlighting the important part of a chart; unmarked
    charts are explicitly flagged as a judge-facing problem -- Idaho STEM
    ISEF Poster Do's/Don'ts.
  - Poster captions should run ~16-18pt, not paper-scale ~9-10pt -- same
    source.
  - One figure should read as the poster's "hero" -- largest, most
    prominent, carrying the main finding.

This module gives every figure script ONE shared way to apply that, so the
whole set reads as one designed system instead of 9 independent one-off
matplotlib defaults.

Two identity colors run through the whole project, used as a small
top-left section badge (not as a forced recolor of each chart's own
correct-by-job internal data encoding -- a heatmap stays sequential, a
3-way categorical chart keeps 3 distinguishable hues, etc.):
    PROJECT_BLUE   -- Aims 1-9,  the EV foundation
    PROJECT_ORANGE -- Aims 10-12, the CAR-T pivot (the project's strongest
                       evidence tier -- see fig8)
"""

import matplotlib.patches as mpatches

PROJECT_BLUE = "#2E5EAA"
PROJECT_ORANGE = "#B0563A"
NEUTRAL_GRAY = "#8A93A6"
INK = "#1C2530"

TITLE_FONTSIZE = 15
TITLE_WEIGHT = "bold"
CAPTION_FONTSIZE_PAPER = 10
CAPTION_FONTSIZE_POSTER = 22  # matches/exceeds ISEF's own 16-18pt poster-caption guidance
TAKEAWAY_FONTSIZE_POSTER = 26

SECTION_LABELS = {
    "EV": ("EV FOUNDATION", PROJECT_BLUE),
    "CAR-T": ("CAR-T PIVOT", PROJECT_ORANGE),
    "OVERVIEW": ("PROJECT OVERVIEW", INK),
}


def clean_spines(ax, keep=("left", "bottom")):
    """Strip chart borders down to the ISEF-guide-recommended minimal
    frame (bottom + left only, by default) instead of a full boxy frame."""
    for spine in ("top", "right", "left", "bottom"):
        ax.spines[spine].set_visible(spine in keep)


def section_badge(fig, section, number):
    """Small rounded badge, top-left of the figure: FIG <n> - <SECTION>.
    Gives every figure an instant, consistent visual identity tag so the
    whole 9-figure set reads as one system rather than 9 independent
    charts. `section` is one of SECTION_LABELS's keys."""
    label_text, color = SECTION_LABELS[section]
    fig.text(
        0.01, 0.99, f"FIG {number}   ·   {label_text}",
        fontsize=9.5, fontweight="bold", color="white",
        va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.45", facecolor=color, edgecolor="none"),
    )


def style_title(ax, text, fontsize=TITLE_FONTSIZE):
    ax.set_title(text, fontsize=fontsize, fontweight=TITLE_WEIGHT, pad=12)


def add_takeaway_banner(fig, text, color=PROJECT_BLUE, y=-0.05):
    """Big, bold, plain-language one/two-line takeaway banner for the
    POSTER version of a figure -- stands in for (or sits above) the small
    italic technical caption used in the paper version, at a font size
    that clears ISEF's own 16-18pt poster-legibility guidance."""
    fig.text(
        0.5, y, text, ha="center", va="top",
        fontsize=TAKEAWAY_FONTSIZE_POSTER, fontweight="bold", color="white",
        bbox=dict(boxstyle="round,pad=0.7", facecolor=color, edgecolor="none"),
        wrap=True,
    )


def add_callout(ax, xy, text, xytext, color=PROJECT_ORANGE, fontsize=16,
                 textcoords="data"):
    """A boxed, arrowed callout pointing at the specific data feature that
    matters -- the "boxes or arrows highlighting important data" ISEF's
    own poster guide calls for, which none of the original 8 figures had.
    Pass textcoords="axes fraction" (with xytext a 0-1 (x, y) pair) to
    place the box at a fixed spot regardless of data scale -- safer than
    data coordinates when the box needs to avoid axis labels/edges."""
    ax.annotate(
        text, xy=xy, xytext=xytext, fontsize=fontsize, fontweight="bold",
        color="white", ha="center", xycoords="data", textcoords=textcoords,
        bbox=dict(boxstyle="round,pad=0.5", facecolor=color, edgecolor="none"),
        arrowprops=dict(arrowstyle="-|>", color=color, linewidth=2.5,
                         shrinkA=0, shrinkB=8),
    )

# ---------------------------------------------------------------------------
# POSTER-MODE helpers. A poster figure is a genuinely different artifact
# from a paper figure, not a bigger export of the same PNG -- see the
# module-level rationale above. Poster variants exist for the 4 figures
# actually selected for the physical board (fig0, fig1, fig5, fig8);
# everything else stays paper-only. Poster mode means: bigger fonts
# throughout, ONE big bold plain-language takeaway banner replacing the
# small italic technical caption (the technical detail stays in the
# paper), and the same badge/callout/clean-spine treatment as the paper
# version.
# ---------------------------------------------------------------------------

POSTER_SCALE = 1.55  # applied to every paper-version fontsize for poster figures


def poster_scale(fontsize):
    return fontsize * POSTER_SCALE
