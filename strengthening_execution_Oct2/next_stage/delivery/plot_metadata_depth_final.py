"""Plot the locked same-cell GSE193531 TNFRSF17 depth-thinning diagnostic."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
NEXT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "strengthening_execution_Oct2" / "dependencies"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter


RESULTS = NEXT / "analysis" / "results"
OUT = NEXT / "delivery" / "figures"
VALUES_OUT = NEXT / "delivery" / "figure_metadata_depth_values.tsv"
STAGE_COLORS = {"NBM": "#38598c", "MM": "#d97835"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    retention = read_tsv(RESULTS / "metadata_depth_retention.tsv")
    paired = [r for r in read_tsv(RESULTS / "metadata_thinning_paired.tsv") if r["gene"] == "TNFRSF17"]
    if len(retention) != 17 or len(paired) != 17:
        raise ValueError(f"Expected 17 retained donors and TNFRSF17 pairs; found {len(retention)} and {len(paired)}")
    order = {"NBM": 0, "MM": 1}
    retention.sort(key=lambda r: (order[r["stage"]], r["sample_id"]))
    paired_by_sample = {r["sample_id"]: r for r in paired}
    rows: list[dict[str, str]] = []
    for r in retention:
        p = paired_by_sample[r["sample_id"]]
        rows.append({
            "stage": r["stage"], "sample_id": r["sample_id"],
            "common_depth_UMI_per_cell": r["common_depth"],
            "eligible_cells": r["eligible_cells"], "retained_cells": r["retained_cells"],
            "retained_fraction": r["retained_fraction"],
            "same_cell_pre_TNFRSF17_detection_fraction": p["same_cells_pre_detection_fraction"],
            "post_thinning_TNFRSF17_detection_MC_mean": p["post_detection_fraction_mc_mean"],
            "post_thinning_MC_q025": p["post_detection_fraction_mc_q025"],
            "post_thinning_MC_q975": p["post_detection_fraction_mc_q975"],
            "MC_draws": p["monte_carlo_draws"],
        })
    OUT.mkdir(parents=True, exist_ok=True)
    with VALUES_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(rows)

    labels = [r["sample_id"] for r in rows]
    colors = [STAGE_COLORS[r["stage"]] for r in rows]
    x = list(range(len(rows)))
    fig, (ax_keep, ax_detect) = plt.subplots(1, 2, figsize=(15, 7.2), sharex=True)
    fig.suptitle("Same-cell TNFRSF17 detection at a fixed 5,189-UMI depth", fontsize=15, y=0.99)

    ax_keep.bar(x, [float(r["retained_fraction"]) for r in rows], color=colors, width=0.72, zorder=3)
    ax_keep.set_title("Cells retained after common-depth thinning")
    ax_keep.set_ylabel("Retained cells / eligible cells")
    ax_keep.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax_keep.set_ylim(0, 1.02)
    ax_keep.grid(axis="y", alpha=0.25, zorder=0)

    for i, r in enumerate(rows):
        pre = float(r["same_cell_pre_TNFRSF17_detection_fraction"])
        post = float(r["post_thinning_TNFRSF17_detection_MC_mean"])
        low = float(r["post_thinning_MC_q025"])
        high = float(r["post_thinning_MC_q975"])
        color = colors[i]
        ax_detect.plot([i - 0.12, i + 0.12], [pre, post], color="#888888", lw=0.8, zorder=1)
        ax_detect.scatter(i - 0.12, pre, marker="s", color=color, s=28, zorder=3)
        ax_detect.errorbar(i + 0.12, post, yerr=[[post - low], [high - post]], fmt="o", color=color,
                           ecolor=color, capsize=2.5, markersize=4.5, lw=1.1, zorder=3)
    ax_detect.set_title("Detected in the same retained cells")
    ax_detect.set_ylabel("Cells with TNFRSF17 molecules / retained cells")
    ax_detect.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax_detect.set_ylim(0, 1.02)
    ax_detect.grid(axis="y", alpha=0.25, zorder=0)

    for axis in (ax_keep, ax_detect):
        axis.set_xticks(x, labels, rotation=45, ha="right", fontsize=8)
        axis.axvline(8.5, color="#555555", lw=0.8, ls="--", alpha=0.6)
    handles = [
        Line2D([0], [0], color=STAGE_COLORS["NBM"], lw=5, label="NBM"),
        Line2D([0], [0], color=STAGE_COLORS["MM"], lw=5, label="MM"),
        Line2D([0], [0], color="#555555", marker="s", lw=0, label="Before thinning"),
        Line2D([0], [0], color="#555555", marker="o", lw=0, label="After thinning (MC mean)"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=4, frameon=False, fontsize=9)
    fig.text(
        0.5, 0.035,
        "Post-thinning whiskers are 2.5th–97.5th percentiles across 100 hypergeometric draws (seed 20261003); they are not donor confidence intervals. Same donor cells before/after; 31 MM-stage Normal-labeled cells excluded.",
        ha="center", fontsize=8,
    )
    fig.tight_layout(rect=(0.02, 0.12, 0.99, 0.89), w_pad=1.7)
    fig.savefig(OUT / "metadata_depth_thinning_final.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "metadata_depth_thinning_final.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / 'metadata_depth_thinning_final.png'} and .svg; values={VALUES_OUT}")


if __name__ == "__main__":
    main()
