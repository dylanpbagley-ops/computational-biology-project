"""Bounded, exploratory metadata and depth diagnostics for GSE193531."""

from __future__ import annotations

import csv
import gzip
import hashlib
import html
import json
import math
import pathlib
import platform
import sys
from collections import Counter, defaultdict

import numpy as np


HERE = pathlib.Path(__file__).resolve().parent
NEXT = HERE.parent
EXECUTION = NEXT.parent
ROOT = EXECUTION.parent
DATA = ROOT / "ev_capture_data" / "healthy_bone_marrow_boiarsky"
GENES = ("TNFRSF17", "GPRC5D", "SLAMF7")
NBM_IDS = tuple(f"NBM-{n}" for n in (1, 2, 3, 4, 6, 7, 8, 10, 11))
MM_IDS = tuple(f"MM-{n}" for n in range(1, 9))
EXPECTED = {"NBM": set(NBM_IDS), "MM": set(MM_IDS)}
SEED = 20261003
N_DRAWS = 100
V2_PROTOCOL = HERE / "diagnostic_protocol_v2.md"
V2_HASH = HERE / "diagnostic_protocol_v2.sha256"
APPEND_ONLY_SNAPSHOT_BACKUPS = {
    "PROJECT_STATUS_Oct2.md": "before/PROJECT_STATUS_Oct2.md",
    "strengthening_execution_Oct2/STATUS.md": "before/execution_STATUS.md",
}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _append_only_snapshot_matches(current_path: pathlib.Path, backup_path: pathlib.Path, expected_sha: str) -> bool:
    """Allow only a byte-preserving suffix append to a backed-up status file."""
    if not current_path.is_file() or not backup_path.is_file() or sha256(backup_path) != expected_sha:
        return False
    return current_path.read_bytes().startswith(backup_path.read_bytes())


def verify_input_snapshot(
    manifest_path: pathlib.Path | None = None,
    root: pathlib.Path | None = None,
    next_stage: pathlib.Path | None = None,
) -> int:
    """Verify protected inputs and the initial snapshot, allowing only backed-up status appends."""
    manifest_path = NEXT / "input_hashes.json" if manifest_path is None else pathlib.Path(manifest_path)
    root = ROOT if root is None else pathlib.Path(root)
    next_stage = NEXT if next_stage is None else pathlib.Path(next_stage)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    mismatches = []
    protected_manifest_rel = manifest.get("protected_manifest_path")
    protected_manifest_sha = manifest.get("protected_manifest_sha256")
    if protected_manifest_rel and protected_manifest_sha:
        path = root / pathlib.PurePosixPath(protected_manifest_rel)
        if not path.is_file() or sha256(path) != protected_manifest_sha:
            mismatches.append(f"changed protected manifest: {protected_manifest_rel}")
    for rel, expected in manifest["protected_inputs"].items():
        path = root / pathlib.PurePosixPath(rel)
        if not path.is_file():
            mismatches.append(f"missing protected input: {rel}")
        elif sha256(path) != expected:
            mismatches.append(f"changed protected input: {rel}")
    for rel, expected in manifest["current_execution_snapshot"].items():
        path = root / pathlib.PurePosixPath(rel)
        if not path.is_file():
            mismatches.append(f"missing starting-state file: {rel}")
        elif sha256(path) != expected:
            backup_rel = APPEND_ONLY_SNAPSHOT_BACKUPS.get(rel)
            backup_path = next_stage / pathlib.PurePosixPath(backup_rel) if backup_rel else None
            if backup_path is None or not _append_only_snapshot_matches(path, backup_path, expected):
                mismatches.append(f"changed starting-state file: {rel}")
    if mismatches:
        raise ValueError("Input preservation check failed: " + "; ".join(mismatches))
    return len(manifest["protected_inputs"])


def source_field_or_unknown(record: dict, field: str) -> str:
    value = record.get(field)
    if value is None or not str(value).strip():
        return "NR"
    return str(value).strip()


def join_metadata_by_cell_id(cell_ids: list[str], metadata_rows: list[dict]) -> list[dict]:
    """Return source rows in expression order, requiring one-to-one IDs."""
    if len(set(cell_ids)) != len(cell_ids):
        raise ValueError("Duplicate expression cell ID")
    by_id = {}
    for row in metadata_rows:
        cell_id = source_field_or_unknown(row, "index")
        if cell_id == "NR":
            raise ValueError("Missing metadata cell ID")
        if cell_id in by_id:
            raise ValueError(f"Duplicate metadata cell ID: {cell_id}")
        by_id[cell_id] = row
    if set(cell_ids) != set(by_id):
        missing = set(cell_ids) - set(by_id)
        extra = set(by_id) - set(cell_ids)
        raise ValueError(f"Cell ID mismatch: missing={len(missing)}, extra={len(extra)}")
    return [by_id[cell_id] for cell_id in cell_ids]


def eligible_sample_groups(metadata_rows: list[dict]) -> dict[str, dict[str, list[int]]]:
    groups: dict[str, dict[str, list[int]]] = {"NBM": defaultdict(list), "MM": defaultdict(list)}
    for index, row in enumerate(metadata_rows):
        stage = source_field_or_unknown(row, "disease_stage")
        label = source_field_or_unknown(row, "normal_or_neoplastic")
        sample_id = source_field_or_unknown(row, "sample_ID")
        if stage == "NBM" and label == "normal":
            if sample_id not in EXPECTED["NBM"]:
                raise ValueError(f"Unknown eligible NBM sample ID: {sample_id}")
            groups["NBM"][sample_id].append(index)
        elif stage == "MM" and label == "neoplastic":
            if sample_id not in EXPECTED["MM"]:
                raise ValueError(f"Unknown eligible MM sample ID: {sample_id}")
            groups["MM"][sample_id].append(index)
    return {stage: dict(samples) for stage, samples in groups.items()}


def integer_count(value, label: str) -> int:
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"Invalid {label}: {value!r}") from error
    if not math.isfinite(number) or number < 0 or not number.is_integer():
        raise ValueError(f"{label} must be a finite nonnegative integer: {value!r}")
    return int(number)


def validate_umi_pair(marker_count, total_count) -> tuple[int, int]:
    marker = integer_count(marker_count, "marker UMI count")
    total = integer_count(total_count, "total UMI count")
    if total <= 0 or marker > total:
        raise ValueError(f"Invalid UMI pair: marker={marker}, total={total}")
    return marker, total


def thin_marker_count(marker_count, total_count, depth: int, rng) -> int:
    marker, total = validate_umi_pair(marker_count, total_count)
    if not isinstance(depth, (int, np.integer)) or depth < 0 or depth > total:
        raise ValueError(f"Invalid common depth {depth} for total {total}")
    if depth == 0:
        return 0
    return int(rng.hypergeometric(marker, total - marker, int(depth)))


def parse_matrix_counts_row(values, line_no: int) -> np.ndarray:
    """Parse one matrix row in NumPy while enforcing integer UMI counts."""
    try:
        counts = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"Invalid UMI matrix values at row {line_no}") from error
    if (counts.ndim != 1 or not np.isfinite(counts).all() or (counts < 0).any()
            or (counts != np.floor(counts)).any() or (counts > np.iinfo(np.int64).max).any()):
        raise ValueError(f"UMI matrix row {line_no} must contain finite nonnegative integer counts")
    return counts.astype(np.int64)


def choose_common_depth(donor_median_depths: list[float]) -> int:
    values = np.asarray(donor_median_depths, dtype=float)
    if values.ndim != 1 or not values.size or not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("Positive finite donor median depths are required")
    return int(math.floor(float(np.quantile(values, 0.10, method="linear"))))


def verify_matrix_column_totals(matrix_path, cell_ids: list[str], expected_totals: list[int]) -> list[int]:
    """Stream the public UMI matrix and require its column sums to match metadata."""
    if len(cell_ids) != len(expected_totals) or len(set(cell_ids)) != len(cell_ids):
        raise ValueError("Expected cell IDs/totals must be one-to-one")
    expected = dict(zip(cell_ids, (integer_count(v, "metadata total UMI count") for v in expected_totals)))
    opener = gzip.open if str(matrix_path).endswith(".gz") else open
    totals = None
    with opener(matrix_path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError("UMI matrix has no cell columns")
        matrix_ids = header[1:]
        if len(set(matrix_ids)) != len(matrix_ids):
            raise ValueError("Duplicate UMI matrix cell ID")
        if set(matrix_ids) != set(cell_ids):
            raise ValueError("UMI matrix cell IDs do not match metadata")
        totals = np.zeros(len(matrix_ids), dtype=np.int64)
        for line_no, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(f"Malformed UMI matrix row {line_no}: {len(row)} columns")
            values = parse_matrix_counts_row(row[1:], line_no)
            if len(values) != len(matrix_ids):
                raise ValueError(f"Malformed UMI matrix row {line_no}: wrong count vector length")
            totals += values
    actual = dict(zip(matrix_ids, map(int, totals)))
    ordered = [actual[cell_id] for cell_id in cell_ids]
    if ordered != [expected[cell_id] for cell_id in cell_ids]:
        mismatch_count = sum(a != b for a, b in zip(ordered, (expected[cell_id] for cell_id in cell_ids)))
        raise ValueError(f"UMI matrix column sums disagree with n_counts for {mismatch_count} cells")
    return ordered


def average_ranks(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2.0
        start = end
    return ranks


def spearman_rho(x, y) -> float | None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 2:
        raise ValueError("Spearman inputs must be paired vectors with at least two values")
    rx, ry = average_ranks(x), average_ranks(y)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return None
    return float(np.corrcoef(rx, ry)[0, 1])


def summarize(values, threshold: float) -> dict:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not values.size or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Summary requires a nonempty finite nonnegative vector")
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "detection_fraction": float(np.mean(values > 0)),
        "p90": float(np.quantile(values, 0.90, method="linear")),
        "tail_fraction": float(np.mean(values >= threshold)),
    }


def quantile_summary(values) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    return (
        float(np.mean(values)),
        float(np.quantile(values, 0.025, method="linear")),
        float(np.quantile(values, 0.975, method="linear")),
    )


def tsv_write(path: pathlib.Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def read_thresholds() -> tuple[dict[str, float], dict[str, float]]:
    path = EXECUTION / "analysis" / "results" / "b_summary.tsv"
    normalized, raw = {}, {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if float(row["budget"]) == 0.1 and row["policy"] == "pooled":
                normalized[row["gene"]] = float(row["threshold"])
    raw_path = EXECUTION / "analysis" / "results" / "b_raw_sensitivity.tsv"
    with raw_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if float(row["budget"]) == 0.1 and row["policy"] == "pooled":
                raw[row["gene"]] = float(row["threshold"])
    if set(normalized) != set(GENES) or set(raw) != set(GENES):
        raise ValueError("Frozen thresholds are incomplete for the three markers")
    return normalized, raw


def q1_median(values) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    return tuple(float(v) for v in np.quantile(values, [0.25, 0.50, 0.75], method="linear"))


def write_depth_tail_svg(path: pathlib.Path, donor_rows: list[dict]) -> None:
    width, height = 960, 360
    left, top, panel_w, panel_h, gap = 72, 42, 260, 250, 32
    colors = {"NBM": "#155e75", "MM": "#b45309"}
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#172033}.axis{stroke:#475569;stroke-width:1}.grid{stroke:#dbe3ec;stroke-width:1}</style>',
             '<text x="480" y="22" text-anchor="middle" font-size="16" font-weight="bold">Donor depth and frozen-threshold tail fraction</text>']
    for j, gene in enumerate(GENES):
        x0 = left + j * (panel_w + gap)
        rows = [r for r in donor_rows if r["gene"] == gene]
        xvals = [float(r["median_n_counts"]) for r in rows]
        ymax = max(0.05, max(float(r["normalized_tail_fraction"]) for r in rows) * 1.12)
        xmin, xmax = min(xvals), max(xvals)
        if xmin == xmax:
            xmax = xmin + 1.0
        tick_rows = []
        for tick_index in range(5):
            fraction = tick_index / 4.0
            y = top + panel_h * (1 - fraction)
            label = f"{ymax * fraction:.0%}"
            depth = xmin + (xmax - xmin) * fraction
            x = x0 + panel_w * fraction
            tick_rows.extend([
                f'<line class="grid" x1="{x0}" y1="{y:.2f}" x2="{x0+panel_w}" y2="{y:.2f}"/>',
                f'<text x="{x0-7}" y="{y+3:.2f}" text-anchor="end" font-size="9">{label}</text>',
                f'<line class="axis" x1="{x:.2f}" y1="{top+panel_h}" x2="{x:.2f}" y2="{top+panel_h+4}"/>',
                f'<text x="{x:.2f}" y="{top+panel_h+15}" text-anchor="middle" font-size="8">{depth:.0f}</text>',
            ])
        parts.extend([
            f'<text x="{x0+panel_w/2}" y="{top-9}" text-anchor="middle" font-size="14" font-weight="bold">{html.escape(gene)}</text>',
            f'<line class="axis" x1="{x0}" y1="{top+panel_h}" x2="{x0+panel_w}" y2="{top+panel_h}"/>',
            f'<line class="axis" x1="{x0}" y1="{top}" x2="{x0}" y2="{top+panel_h}"/>',
            *tick_rows,
            f'<text x="{x0+panel_w/2}" y="{top+panel_h+31}" text-anchor="middle" font-size="11">Median cell UMI depth</text>',
            f'<text x="{x0-56}" y="{top+panel_h/2}" text-anchor="middle" font-size="10" transform="rotate(-90 {x0-56} {top+panel_h/2})">Fraction above pooled 10% threshold</text>',
        ])
        for row in rows:
            x = x0 + panel_w * (float(row["median_n_counts"]) - xmin) / (xmax - xmin)
            y = top + panel_h * (1 - float(row["normalized_tail_fraction"]) / ymax)
            c = colors[row["stage"]]
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.5" fill="{c}" fill-opacity="0.8"><title>{html.escape(row["stage"])} {html.escape(row["sample_id"])}: depth={row["median_n_counts"]:.1f}, tail={row["normalized_tail_fraction"]:.4f}</title></circle>')
    parts.extend(['<circle cx="746" cy="335" r="5" fill="#155e75"/><text x="757" y="339" font-size="11">NBM</text>', '<circle cx="811" cy="335" r="5" fill="#b45309"/><text x="822" y="339" font-size="11">MM</text>', '<text x="480" y="355" text-anchor="middle" font-size="10">Descriptive sample-level association only; no causal adjustment or population inference.</text>', '</svg>'])
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def depth_tail_figure_path() -> pathlib.Path:
    return NEXT / "delivery" / "figures" / "metadata_depth_tail.svg"


def read_local_inputs():
    extracted = json.loads((DATA / "extracted_markers.json").read_text(encoding="utf-8"))
    with gzip.open(DATA / "GSE193531_cell-level-metadata.csv.gz", "rt", newline="", encoding="utf-8") as handle:
        metadata = list(csv.DictReader(handle))
    metadata = join_metadata_by_cell_id(extracted["cell_ids"], metadata)
    groups = eligible_sample_groups(metadata)
    if set(groups["NBM"]) != EXPECTED["NBM"] or set(groups["MM"]) != EXPECTED["MM"]:
        raise ValueError("Eligible donor IDs differ from the locked 9 NBM / 8 MM set")
    totals = np.asarray([integer_count(row["n_counts"], "n_counts") for row in metadata], dtype=np.int64)
    if (totals <= 0).any():
        raise ValueError("Every selected or unselected cell must have positive n_counts")
    raw = {}
    normalized = {}
    for gene in GENES:
        values = np.asarray(extracted["expr"][gene], dtype=float)
        if len(values) != len(metadata) or not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError(f"Marker data are not finite nonnegative integer UMI counts: {gene}")
        if (values > totals).any():
            raise ValueError(f"Marker UMI count exceeds n_counts for {gene}")
        raw[gene] = values.astype(np.int64)
        normalized[gene] = values / totals * 10000.0
    for stage in groups:
        for sample_id, ix in groups[stage].items():
            if not ix:
                raise ValueError(f"No selected cells for {sample_id}")
    return extracted["cell_ids"], metadata, groups, totals, raw, normalized


def main() -> None:
    protected_count = verify_input_snapshot()
    protocol_hash = sha256(V2_PROTOCOL)
    recorded_hash = V2_HASH.read_text(encoding="utf-8").strip()
    if protocol_hash != recorded_hash:
        raise ValueError("Diagnostic protocol v2 hash changed")
    cell_ids, metadata, groups, totals, raw, normalized = read_local_inputs()
    frozen_norm, frozen_raw = read_thresholds()

    donor_rows = []
    donor_cache = {}
    for stage in ("NBM", "MM"):
        for sample_id in sorted(groups[stage]):
            ix = groups[stage][sample_id]
            sample_meta = [metadata[i] for i in ix]
            count_summary = {name: q1_median([float(row[name]) for row in sample_meta]) for name in ("n_counts", "n_genes", "frac_mito")}
            for gene in GENES:
                raw_values = raw[gene][ix]
                norm_values = normalized[gene][ix]
                raw_metrics = summarize(raw_values, frozen_raw[gene])
                norm_metrics = summarize(norm_values, frozen_norm[gene])
                row = {
                    "stage": stage,
                    "sample_id": sample_id,
                    "gene": gene,
                    "eligible_cells": len(ix),
                    "median_n_counts": count_summary["n_counts"][1],
                    "n_counts_q1": count_summary["n_counts"][0],
                    "n_counts_q3": count_summary["n_counts"][2],
                    "median_n_genes": count_summary["n_genes"][1],
                    "n_genes_q1": count_summary["n_genes"][0],
                    "n_genes_q3": count_summary["n_genes"][2],
                    "median_frac_mito": count_summary["frac_mito"][1],
                    "frac_mito_q1": count_summary["frac_mito"][0],
                    "frac_mito_q3": count_summary["frac_mito"][2],
                    "raw_mean": raw_metrics["mean"],
                    "raw_median": raw_metrics["median"],
                    "raw_detection_fraction": raw_metrics["detection_fraction"],
                    "raw_p90": raw_metrics["p90"],
                    "raw_tail_fraction": raw_metrics["tail_fraction"],
                    "normalized_mean": norm_metrics["mean"],
                    "normalized_median": norm_metrics["median"],
                    "normalized_detection_fraction": norm_metrics["detection_fraction"],
                    "normalized_p90": norm_metrics["p90"],
                    "normalized_tail_fraction": norm_metrics["tail_fraction"],
                    "frozen_raw_threshold": frozen_raw[gene],
                    "frozen_normalized_threshold": frozen_norm[gene],
                }
                donor_rows.append(row)
                donor_cache[(stage, sample_id, gene)] = {"indices": ix, "raw": raw_values, "normalized": norm_values, "row": row}
    tsv_write(HERE / "results" / "metadata_donors.tsv", donor_rows)
    print(f"DONOR_DIAGNOSTICS_WRITTEN rows={len(donor_rows)} donors=17 NBM_cells={sum(len(v) for v in groups['NBM'].values())} MM_neoplastic_cells={sum(len(v) for v in groups['MM'].values())}", flush=True)

    association_rows = []
    for stage in ("NBM", "MM"):
        samples = sorted(groups[stage])
        depth = [float(np.median(totals[groups[stage][sample]])) for sample in samples]
        for gene in GENES:
            for outcome in ("normalized_tail_fraction", "normalized_detection_fraction", "raw_tail_fraction"):
                y = [float(donor_cache[(stage, sample, gene)]["row"][outcome]) for sample in samples]
                association_rows.append({
                    "stage": stage,
                    "gene": gene,
                    "depth_summary": "donor_median_n_counts",
                    "marker_summary": outcome,
                    "donor_n": len(samples),
                    "spearman_rho": ("NR_constant_rank" if spearman_rho(depth, y) is None else spearman_rho(depth, y)),
                    "sample_ids": ",".join(samples),
                    "interpretation": "descriptive within-stratum association; not causal",
                })
    tsv_write(HERE / "results" / "metadata_associations.tsv", association_rows)
    write_depth_tail_svg(depth_tail_figure_path(), donor_rows)

    matrix_path = DATA / "GSE193531_umi-count-matrix.csv.gz"
    matrix_status = "verified"
    matrix_error = None
    try:
        verified_totals = verify_matrix_column_totals(matrix_path, cell_ids, totals.tolist())
        if verified_totals != totals.tolist():
            raise ValueError("Returned matrix totals do not match metadata order")
    except (ValueError, StopIteration, OSError) as error:
        matrix_status = "refused"
        matrix_error = str(error)
    thinning_rows = []
    depth_info = {"status": matrix_status, "reason": matrix_error}
    if matrix_status == "verified":
        donor_medians = {stage: {sample: float(np.median(totals[ix])) for sample, ix in groups[stage].items()} for stage in groups}
        common_depth = choose_common_depth([v for stage in donor_medians.values() for v in stage.values()])
        retained = {stage: {sample: [i for i in ix if totals[i] >= common_depth] for sample, ix in group.items()} for stage, group in groups.items()}
        no_cells = [sample for stage in retained.values() for sample, ix in stage.items() if not ix]
        if no_cells:
            depth_info = {"status": "refused", "reason": "No cells meet common depth for donors: " + ",".join(no_cells)}
        else:
            retained_rows = []
            for stage in ("NBM", "MM"):
                for sample in sorted(groups[stage]):
                    retained_rows.append({"stage": stage, "sample_id": sample, "eligible_cells": len(groups[stage][sample]), "retained_cells": len(retained[stage][sample]), "retained_fraction": len(retained[stage][sample]) / len(groups[stage][sample]), "common_depth": common_depth})
            tsv_write(HERE / "results" / "metadata_depth_retention.tsv", retained_rows)
            rng = np.random.default_rng(SEED)
            for gene in GENES:
                for stage in ("NBM", "MM"):
                    for sample in sorted(groups[stage]):
                        indices = retained[stage][sample]
                        n_values = totals[indices]
                        before_values = raw[gene][indices] / n_values * 10000.0
                        threshold = frozen_norm[gene]
                        before = summarize(before_values, threshold)
                        post_by_metric = defaultdict(list)
                        for replicate in range(1, N_DRAWS + 1):
                            sampled = np.fromiter((thin_marker_count(raw[gene][i], totals[i], common_depth, rng) for i in indices), dtype=np.int64, count=len(indices))
                            after_values = sampled / common_depth * 10000.0
                            after = summarize(after_values, threshold)
                            for metric, value in after.items():
                                post_by_metric[metric].append(value)
                        row = {
                            "stage": stage,
                            "sample_id": sample,
                            "gene": gene,
                            "common_depth": common_depth,
                            "retained_cells": len(indices),
                            "monte_carlo_draws": N_DRAWS,
                        }
                        for metric in ("detection_fraction", "median", "p90", "tail_fraction"):
                            mean, low, high = quantile_summary(post_by_metric[metric])
                            row[f"same_cells_pre_{metric}"] = before[metric]
                            row[f"post_{metric}_mc_mean"] = mean
                            row[f"post_{metric}_mc_q025"] = low
                            row[f"post_{metric}_mc_q975"] = high
                        thinning_rows.append(row)
            tsv_write(HERE / "results" / "metadata_thinning_paired.tsv", thinning_rows)
            depth_info = {
                "status": "completed",
                "reason": "All UMI and matrix-total checks passed; all 17 donors retained cells.",
                "common_depth_m": common_depth,
                "m_rule": "floor of linear-interpolated Q10 of 17 donor median n_counts",
                "seed": SEED,
                "monte_carlo_draws": N_DRAWS,
                "same_retained_cells_pre_post": True,
                "marker_draws_are_marginal": True,
                "retained_cell_counts": retained_rows,
            }
    else:
        depth_info["thinning_omitted"] = True
        depth_info["refusal_scope"] = "No count-based thinning or corrected biological interpretation is made."

    runtime = {"python": platform.python_version(), "numpy": np.__version__}
    report = {
        "input_preservation_verified": True,
        "protected_input_count": protected_count,
        "protocol_v2_sha256": protocol_hash,
        "seed": SEED,
        "monte_carlo_draws": N_DRAWS,
        "eligibility": {
            "NBM_samples": sorted(groups["NBM"]),
            "NBM_cells": sum(len(v) for v in groups["NBM"].values()),
            "MM_samples": sorted(groups["MM"]),
            "MM_neoplastic_cells": sum(len(v) for v in groups["MM"].values()),
            "MM_stage_normal_labeled_cells_excluded": sum(1 for r in metadata if r["disease_stage"] == "MM" and r["normal_or_neoplastic"] == "normal"),
        },
        "diagnostic_scope": "descriptive sample-level associations; no causal attribution or clinical endpoint",
        "metadata_covariates_not_in_cell_file": ["age", "sex", "preparation", "batch", "immunoglobulin composition"],
        "count_matrix_verification": depth_info,
        "runtime": runtime,
        "outputs": [
            "results/metadata_donors.tsv",
            "results/metadata_associations.tsv",
            "results/metadata_depth_retention.tsv if thinning completed",
            "results/metadata_thinning_paired.tsv if thinning completed",
            "delivery/figures/metadata_depth_tail.svg",
        ],
    }
    (HERE / "results" / "metadata_diagnostics.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
