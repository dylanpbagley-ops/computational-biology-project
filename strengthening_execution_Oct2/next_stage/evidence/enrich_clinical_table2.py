"""Apply exact event counts from the indexed primary Table 2 to the fixed safety rows."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
INPUT = HERE / "clinical_observations_normalized_v3_pre_table2.tsv"
OUTPUT = HERE / "clinical_observations.tsv"
MANIFEST = HERE / "clinical_observations_normalization.json"
BASE_MANIFEST = HERE / "clinical_observations_normalization_v3_pre_table2.json"
OBSERVATION_STATES = {"event_supported", "explicit_nonobservation", "indeterminate", "ineligible"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    if not INPUT.exists() or not BASE_MANIFEST.exists():
        raise FileNotFoundError("The preserved 18-row base extraction and its manifest are required")
    base_manifest = json.loads(BASE_MANIFEST.read_text(encoding="utf-8"))
    if digest(INPUT) != base_manifest["normalized_observations_sha256"]:
        raise ValueError("Preserved base observation table does not match its recorded SHA-256")
    with INPUT.open("r", encoding="utf-8", newline="") as handle:
        base_rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = list(base_rows[0])
    by_id = {row["observation_id"]: row for row in base_rows}
    if len(base_rows) != 18 or len(by_id) != 18:
        raise ValueError("Expected the preserved 18-row event extraction, with unique IDs")

    replace_ids = {
        "OBS-TAL-SKIN-405-ANY", "OBS-TAL-RASH405-G34", "OBS-TAL-SKIN-800-ANY",
        "OBS-TAL-RASH800-G3", "OBS-TAL-NAIL-405", "OBS-TAL-NAIL-800",
        "OBS-TAL-ORAL-405", "OBS-TAL-ORAL-800",
    }
    if not replace_ids.issubset(by_id):
        raise ValueError("Expected source-normalized talquetamab pilot rows were not found")

    cohort_defs = {
        "405_QW": {
            "cohort_id": "CO-TAL-RP2D-405-QW", "n": 30,
            "denominator_description": "patients treated subcutaneously at 405 µg/kg weekly RP2D",
            "exposure_activity": "Talquetamab 405 µg/kg weekly subcutaneous RP2D.",
            "followup": "11.7",
        },
        "800_Q2W": {
            "cohort_id": "CO-TAL-RP2D-800-Q2W", "n": 44,
            "denominator_description": "patients treated subcutaneously at 800 µg/kg every other week RP2D",
            "exposure_activity": "Talquetamab 800 µg/kg every other week subcutaneous RP2D.",
            "followup": "4.2",
        },
    }
    events = {
        "SKIN": {"name": "skin-related adverse events", "tissue": "skin", "any": (20, 31), "g34": (0, 1), "pct": ("67", "67%", "70", "70%")},
        "NAIL": {"name": "nail-related adverse events", "tissue": "nails", "any": (17, 12), "g34": (0, 1), "pct": ("57", "57%", "27", "27%")},
        "RASH": {"name": "rash-related adverse events", "tissue": "skin", "any": (14, 13), "g34": (0, 7), "pct": ("47", "47%", "30", "30%"), "severe_pct": {"405_QW": ("0", "0%"), "800_Q2W": ("16", "16%")}},
        "DYSGEUSIA": {"name": "dysgeusia adverse events", "tissue": "taste function", "any": (19, 25), "g34": ("NA", "NA"), "pct": None},
        "DRY_MOUTH": {"name": "dry-mouth adverse events", "tissue": "oral mucosa", "any": (9, 25), "g34": (0, 0), "pct": None},
        "DYSPHAGIA": {"name": "dysphagia adverse events", "tissue": "swallowing function", "any": (11, 12), "g34": (0, 0), "pct": None},
    }

    result = [row.copy() for row in base_rows if row["observation_id"] not in replace_ids]
    for cohort_index, (cohort_key, cohort) in enumerate(cohort_defs.items()):
        template = by_id["OBS-TAL-SKIN-405-ANY" if cohort_key == "405_QW" else "OBS-TAL-SKIN-800-ANY"]
        for event_key, event in events.items():
            any_count = event["any"][cohort_index]
            severe_count = event["g34"][cohort_index]
            percent_any = event["pct"]
            if percent_any is None:
                any_pct, any_pct_display, any_pct_state = "NR", "NR", "not extracted"
            else:
                any_pct, any_pct_display = percent_any[cohort_index * 2:cohort_index * 2 + 2]
                any_pct_state = "source-reported rounded percentage"

            for grade_suffix, count, grade, pct, pct_display, pct_state in (
                ("ANY", any_count, "any_grade", any_pct, any_pct_display, any_pct_state),
                ("G34", severe_count, "grade_3_or_4", "NR", "NR", "not extracted"),
            ):
                row = template.copy()
                observation_id = f"OBS-TAL-{event_key}-{cohort_key}-{grade_suffix}"
                row.update({
                    "observation_id": observation_id,
                    "event_group_id": f"TAL-{event_key}",
                    "cohort_id": cohort["cohort_id"],
                    "event_name": event["name"],
                    "normal_tissue": event["tissue"],
                    "event_count": str(count),
                    "event_count_state": (
                        "source_NA" if count == "NA" else
                        "source_explicit_no_reported_events" if count == 0 and grade_suffix == "G34" else
                        "reported_exact_count"
                    ),
                    "event_percentage": pct,
                    "event_percentage_display": pct_display,
                    "percentage_source_state": pct_state,
                    "denominator": str(cohort["n"]),
                    "denominator_unit": "patients",
                    "denominator_description": cohort["denominator_description"],
                    "grade": grade,
                    "event_context": "Talquetamab × CD3 bispecific antibody; phase 1/2 MonumenTAL-1 RP2D cohort. Table 2 gives the exact event count by grade.",
                    "event_window": f"Median follow-up {cohort['followup']} months; event-specific onset/window is not given in extracted Table 2 text.",
                    "attribution_location": "NEJM Table 2, Safety",
                    "exposure_activity": cohort["exposure_activity"],
                    "observation_state": "indeterminate",
                    "adequacy_limits": (
                        "Table 2 records NA for grade 3/4 dysgeusia; unavailable is not a zero. This trial-level event does not establish target-specific injury."
                        if count == "NA" else
                        "Safety events in a single-arm treated population do not isolate GPRC5D-mediated normal-tissue injury; independent causal adjudication was not performed."
                    ),
                    "attribution_context": (
                        "Authors describe nail toxicity as probably related to on-target activity; that is a mechanistic interpretation, not independent patient-level causal adjudication."
                        if event_key == "NAIL" else
                        "Observed during treatment; incidence alone does not establish target-specific injury."
                    ),
                    "review_status": "pending independent review",
                })
                if event_key == "RASH" and grade_suffix == "G34":
                    row["event_percentage"], row["event_percentage_display"] = event.get("severe_pct", {}).get(cohort_key, ("NR", "NR"))
                    row["percentage_source_state"] = "source-reported rounded percentage" if row["event_percentage"] != "NR" else "not extracted"
                    if cohort_key == "800_Q2W":
                        row["supportive_care"] = "Grade 3 rashes responded to oral and topical glucocorticoids."
                        row["recovery"] = "All but one patient with grade 3 rash resumed treatment after resolution."
                result.append(row)

    ids = [row["observation_id"] for row in result]
    if len(result) != 34 or len(ids) != len(set(ids)):
        raise ValueError(f"Expected 34 unique event/grade rows, observed {len(result)} rows")
    if any(row["denominator_unit"] != "patients" for row in result):
        raise ValueError("All denominators must count patients")
    for row in result:
        value = row["event_count"]
        if value not in {"NR", "NA"} and (not value.isdigit() or int(value) < 0):
            raise ValueError(f"Invalid event count in {row['observation_id']}: {value}")
        if row["observation_state"] not in OBSERVATION_STATES:
            raise ValueError(f"Invalid endpoint observation_state in {row['observation_id']}: {row['observation_state']}")
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(result)

    manifest = {
        **base_manifest,
        "normalization_protocol": "18-row fixed-source extraction followed by exact Table 2 enrichment",
        "observation_state_contract": {
            "allowed_values": sorted(OBSERVATION_STATES),
            "assigned_value": "indeterminate",
            "reason": "trial event facts and source NA/zero counts are preserved separately; target attribution and endpoint qualification remain unadjudicated",
            "validated": True,
        },
        "normalization_input_sha256": digest(INPUT),
        "normalized_observations_sha256": digest(OUTPUT),
        "normalized_event_grade_rows": len(result),
        "source_access_note": "Indexed primary NEJM Results and Table 2 text supplied the bounded extracted results; direct article page opens returned HTTP 403. Table 2 exact event counts were retained; supplementary tables were not retrieved. Review remains pending.",
        "table2_enrichment": {
            "source": "NEJM Table 2, Safety, doi:10.1056/NEJMoa2204591",
            "exact_count_events": "skin, nail, rash, dysgeusia, dry mouth, and dysphagia, any grade and grade 3/4, by RP2D cohort",
            "dysgeusia_grade_3_or_4": "NA retained as source NA; never converted to zero",
            "percentages": "only source-reported percentages retained; no percentage-to-count reconstruction",
            "causal_labels_changed": False,
        },
        "labels_predictions_and_clinical_AUC_changed": False,
        "review_status": "pending independent review",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
