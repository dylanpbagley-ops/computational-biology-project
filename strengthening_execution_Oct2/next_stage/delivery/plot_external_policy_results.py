"""Render the locked GSE305971 10% policy comparison from verified TSV outputs."""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
NEXT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "strengthening_execution_Oct2" / "dependencies"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, PercentFormatter


RESULTS = NEXT / "analysis" / "results"
OUT = NEXT / "delivery" / "figures"
PLOT_DATA = NEXT / "delivery" / "figure_policy_plot_values.tsv"
POLICIES = ("pooled", "equal_code", "donor_envelope")
POLICY_LABELS = {
    "pooled": "Pooled",
    "equal_code": "Equal code",
    "donor_envelope": "Donor envelope",
}
COLORS = {"pooled": "#38598c", "equal_code": "#299d8f", "donor_envelope": "#d97835"}
MARKERS = ("TNFRSF17", "GPRC5D", "SLAMF7")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    summary_rows = read_tsv(RESULTS / "external_full_newcode_thresholds.tsv")
    bound_rows = read_tsv(RESULTS / "external_full_newcode_missingness_bounds.tsv")
    primary = [
        r for r in summary_rows
        if r["analysis_set"] == "fixed_pilot_reference_new_codes" and float(r["budget"]) == 0.1
    ]
    primary_bounds = [
        r for r in bound_rows
        if r["analysis_set"] == "fixed_pilot_reference_new_codes" and float(r["budget"]) == 0.1
    ]
    if len(primary) != 27 or len(primary_bounds) != 3 * 3 * (9 + 90 + 7):
        raise ValueError(f"Unexpected fixed 10% row counts: summary={len(primary)}, bounds={len(primary_bounds)}")

    summary_by_key = {(r["marker"], r["policy"], r["heldout_group"]): r for r in primary}
    bounds_by_key: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in primary_bounds:
        bounds_by_key[(row["marker"], row["policy"], row["heldout_group"])].append(row)

    values: list[dict[str, str | int | float]] = []
    for marker in MARKERS:
        for policy in POLICIES:
            healthy_key = (marker, policy, "heldout_new_healthy_normal_pc")
            mm_key = (marker, policy, "heldout_new_mm_malignant_label")
            healthy_rows = bounds_by_key[healthy_key]
            mm_rows = bounds_by_key[mm_key]
            healthy_summary = summary_by_key[healthy_key]
            healthy_obs = int(healthy_summary["healthy_codes_above_normal_budget"])
            healthy_max = sum(float(r["fraction_upper_all_missing_positive"]) > 0.1 for r in healthy_rows)
            mm_summary = summary_by_key[mm_key]
            mm_n = sum(int(r["matched_cells"]) for r in mm_rows)
            mm_annotation_n = sum(int(r["annotated_cells"]) for r in mm_rows)
            mm_positive = sum(int(r["positive_matched_cells"]) for r in mm_rows)
            mm_missing = sum(int(r["unmatched_annotation_cells"]) for r in mm_rows)
            mm_matched_fraction = mm_positive / mm_n
            mm_lower = mm_positive / mm_annotation_n
            mm_upper = (mm_positive + mm_missing) / mm_annotation_n
            row = {
                "marker": marker,
                "policy": policy,
                "healthy_codes_observed_over_budget": int(healthy_obs),
                "healthy_codes_max_possible_over_budget": int(healthy_max),
                "healthy_codes_evaluated": int(healthy_summary["heldout_codes"]),
                "healthy_cells_matched": int(healthy_summary["heldout_cells"]),
                "mm_positive_cells_matched": mm_positive,
                "mm_cells_matched": mm_n,
                "mm_annotated_cells": mm_annotation_n,
                "mm_unmatched_cells": mm_missing,
                "mm_matched_fraction": mm_matched_fraction,
                "mm_all_missing_negative_lower": mm_lower,
                "mm_all_missing_positive_upper": mm_upper,
                "training_threshold": float(mm_summary["threshold_per_10k"]),
                "threshold_state": mm_summary["threshold_state"],
            }
            values.append(row)

    OUT.mkdir(parents=True, exist_ok=True)
    columns = list(values[0])
    with PLOT_DATA.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(values)

    fig, (ax_healthy, ax_mm) = plt.subplots(1, 2, figsize=(13.5, 6.2), gridspec_kw={"width_ratios": [0.93, 1.25]})
    fig.suptitle("Held-out code outcomes under a 10% Healthy reference budget", fontsize=15, y=0.98)

    group_x = {marker: i for i, marker in enumerate(MARKERS)}
    offsets = {"pooled": -0.24, "equal_code": 0.0, "donor_envelope": 0.24}
    for policy in POLICIES:
        selected = [next(r for r in values if r["marker"] == marker and r["policy"] == policy) for marker in MARKERS]
        xs = [group_x[str(row["marker"])] + offsets[policy] for row in selected]
        healthy_y = [int(r["healthy_codes_observed_over_budget"]) for r in selected]
        healthy_hi = [int(r["healthy_codes_max_possible_over_budget"]) for r in selected]
        healthy_err = [hi - low for low, hi in zip(healthy_y, healthy_hi)]
        ax_healthy.bar(xs, healthy_y, color=COLORS[policy], width=0.22, label=POLICY_LABELS[policy], zorder=3)
        ax_healthy.errorbar(xs, healthy_y, yerr=[[0] * len(xs), healthy_err], fmt="none", ecolor="#333333", capsize=3, lw=1.2, zorder=4)
        for x, low, high in zip(xs, healthy_y, healthy_hi):
            ax_healthy.text(x, high + 0.12, f"{low}" if low == high else f"{low}–{high}", ha="center", va="bottom", fontsize=8)
    ax_healthy.set_xticks(list(group_x.values()), list(group_x.keys()), fontsize=8)
    ax_healthy.set_ylabel("Healthy codes above 10% positive")
    ax_healthy.set_title("Held-out Healthy codes (n=9)")
    ax_healthy.set_xlim(-0.55, len(MARKERS) - 0.45)
    ax_healthy.set_ylim(0, 5.1)
    ax_healthy.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax_healthy.grid(axis="y", alpha=0.25, zorder=0)
    ax_healthy.text(
        0.02, -0.23,
        "Bars: exact-matched outcomes. Whiskers: maximum if all unmatched\nannotations were positive; not an interval estimate.",
        transform=ax_healthy.transAxes, fontsize=8, va="top",
    )

    xbase = list(range(len(MARKERS)))
    offsets = {"pooled": -0.19, "equal_code": 0.0, "donor_envelope": 0.19}
    for policy in POLICIES:
        selected = [next(r for r in values if r["marker"] == marker and r["policy"] == policy) for marker in MARKERS]
        xs = [xbase[i] + offsets[policy] for i in range(len(MARKERS))]
        lows = [float(r["mm_all_missing_negative_lower"]) for r in selected]
        highs = [float(r["mm_all_missing_positive_upper"]) for r in selected]
        observed = [float(r["mm_matched_fraction"]) for r in selected]
        for x, low, high in zip(xs, lows, highs):
            ax_mm.vlines(x, low, high, color=COLORS[policy], alpha=0.75, lw=2)
            ax_mm.hlines([low, high], x - 0.045, x + 0.045, color=COLORS[policy], alpha=0.8, lw=1.2)
        ax_mm.scatter(xs, observed, color=COLORS[policy], marker="o", s=36, label=POLICY_LABELS[policy], zorder=3)

    ax_mm.set_xticks(xbase, MARKERS)
    ax_mm.set_ylabel("MM/Malignant positive cells (%)")
    ax_mm.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax_mm.set_ylim(0, 0.42)
    ax_mm.set_title("Held-out MM cells (n=90 codes; 44,916 matched)")
    ax_mm.grid(axis="y", alpha=0.25, zorder=0)
    ax_mm.legend(frameon=False, loc="upper right")
    ax_mm.text(
        0.02, -0.22,
        "Dots: matched-cell fraction. Whiskers: arithmetic bounds if all unmatched\nannotations were negative or positive; they are not confidence intervals.",
        transform=ax_mm.transAxes, fontsize=8, va="top",
    )
    fig.text(
        0.5, 0.02,
        "Same-study MARS assay; source-labelled Healthy/Normal_PC and MM/Malignant groups. No clinical or patient-level inference.",
        ha="center", fontsize=8,
    )
    fig.tight_layout(rect=(0.02, 0.13, 0.99, 0.93), w_pad=2.1)
    fig.savefig(OUT / "gse305971_full_policy_tradeoff.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "gse305971_full_policy_tradeoff.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / 'gse305971_full_policy_tradeoff.png'} and .svg; values={PLOT_DATA}")


if __name__ == "__main__":
    main()
