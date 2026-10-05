"""Split bounded primary-study safety extracts into source-faithful event/grade rows."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "clinical_observations.tsv"
BACKUP = HERE / "clinical_observations_pre_normalization.tsv"
OUTPUT = HERE / "clinical_observations.tsv"
MANIFEST = HERE / "clinical_observations_normalization.json"
OBSERVATION_STATES = {"event_supported", "explicit_nonobservation", "indeterminate", "ineligible"}


FIELDS = (
    "observation_id", "event_group_id", "study_id", "cohort_id", "source_id", "event_name", "normal_tissue",
    "event_count", "event_count_state", "event_percentage", "event_percentage_display", "percentage_source_state",
    "denominator", "denominator_unit", "denominator_description", "grade", "event_context", "event_window",
    "injury_cessation", "death", "protocol_DLT", "attribution_code", "attribution_context", "attribution_location",
    "depletion", "consequences", "supportive_care", "recovery", "ascertainment", "exposure_activity",
    "missingness_censoring", "observation_state", "adequacy_limits", "source_primary_url", "extractor", "review_status",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def row(
    observation_id: str, event_group_id: str, study_id: str, cohort_id: str, source_id: str,
    event_name: str, normal_tissue: str, event_count: int | str, event_count_state: str,
    event_percentage: int | str, event_percentage_display: str, percentage_source_state: str,
    denominator: int | str, denominator_description: str, grade: str, event_context: str,
    event_window: str, attribution_location: str, adequacy_limits: str,
    *, injury_cessation: str = "NR", death: str = "NR", protocol_DLT: str = "NA",
    attribution_code: str = "unresolved", attribution_context: str = "",
    depletion: str = "NR", consequences: str = "NR", supportive_care: str = "NR",
    recovery: str = "NR", ascertainment: str = "trial adverse-event reporting",
    exposure_activity: str = "", missingness_censoring: str = "event-specific censoring not extracted",
    observation_state: str = "indeterminate", review_status: str = "pending independent review",
) -> dict[str, object]:
    return {
        "observation_id": observation_id, "event_group_id": event_group_id, "study_id": study_id,
        "cohort_id": cohort_id, "source_id": source_id, "event_name": event_name, "normal_tissue": normal_tissue,
        "event_count": event_count, "event_count_state": event_count_state,
        "event_percentage": event_percentage, "event_percentage_display": event_percentage_display,
        "percentage_source_state": percentage_source_state, "denominator": denominator,
        "denominator_unit": "patients", "denominator_description": denominator_description,
        "grade": grade, "event_context": event_context, "event_window": event_window,
        "injury_cessation": injury_cessation, "death": death, "protocol_DLT": protocol_DLT,
        "attribution_code": attribution_code, "attribution_context": attribution_context,
        "attribution_location": attribution_location, "depletion": depletion, "consequences": consequences,
        "supportive_care": supportive_care, "recovery": recovery, "ascertainment": ascertainment,
        "exposure_activity": exposure_activity, "missingness_censoring": missingness_censoring,
        "observation_state": observation_state, "adequacy_limits": adequacy_limits,
        "source_primary_url": (
            "https://www.nejm.org/doi/full/10.1056/NEJMoa2024850" if source_id == "S-TNFRSF17-KARMMA"
            else "https://www.nejm.org/doi/full/10.1056/NEJMoa2204591"
        ),
        "extractor": "assistant; bounded primary-source extraction", "review_status": review_status,
    }


def main() -> None:
    if not BACKUP.exists():
        shutil.copyfile(SOURCE, BACKUP)
    with BACKUP.open("r", encoding="utf-8", newline="") as handle:
        old = list(csv.DictReader(handle, delimiter="\t"))
    if len(old) != 7:
        raise ValueError(f"Expected the preserved seven-row mixed draft, found {len(old)} rows")
    expected_ids = {
        "OBS-KARMMA-HYPOGAMMA", "OBS-KARMMA-INFECTION", "OBS-KARMMA-CYTOPENIA-CONTEXT",
        "OBS-TAL-DLT-RASH", "OBS-TAL-SKIN-405", "OBS-TAL-SKIN-800", "OBS-TAL-ORAL-NAIL",
    }
    if {r["observation_id"] for r in old} != expected_ids:
        raise ValueError("Clinical draft rows changed; inspect rather than normalize unexpected inputs")

    rows: list[dict[str, object]] = []
    kar = "ST-TNFRSF17-KARMMA"
    kar_cohort = "CO-KARMMA-INFUSED"
    kar_src = "S-TNFRSF17-KARMMA"
    kar_context = "Idecabtagene vicleucel (ide-cel) phase 2 KarMMa treated cohort."
    kar_window = (
        "Adverse events were collected through trial follow-up; source notes that most adverse events, except infections and "
        "hypogammaglobulinemia, occurred within the first 8 weeks. Event-specific onset and risk window were not extracted."
    )
    kar_limits = (
        "Treatment-emergent event counts are not a direct count of normal plasma-cell depletion or proof of BCMA-mediated causality; "
        "patient-level trajectories and exact event-specific censoring were not extracted."
    )
    rows.extend([
        row("OBS-KARMMA-HYPOGAMMA-ANY", "KARMMA-HYPOGAMMA", kar, kar_cohort, kar_src,
            "treatment-emergent hypogammaglobulinemia adverse event", "humoral immune function; indirect treatment-emergent endpoint",
            27, "reported_exact_count", 21, "21%", "source-reported rounded percentage", 128,
            "ide-cel infused patients", "any_grade", kar_context, kar_window,
            "NEJM Table 2; Results Safety timing statement", kar_limits,
            attribution_context="Observed during ide-cel plus lymphodepletion; target-specific causality not adjudicated.",
            depletion="Direct normal plasma-cell-depletion count not extracted.",
            supportive_care="Immune globulin use described as common; exact count not extracted.",
            exposure_activity="ide-cel anti-BCMA CAR-T plus protocol lymphodepletion."),
        row("OBS-KARMMA-HYPOGAMMA-G34", "KARMMA-HYPOGAMMA", kar, kar_cohort, kar_src,
            "treatment-emergent hypogammaglobulinemia adverse event", "humoral immune function; indirect treatment-emergent endpoint",
            1, "reported_exact_count", "NR", "<1%", "source-reported less-than percentage", 128,
            "ide-cel infused patients", "grade_3_or_4", kar_context, kar_window,
            "NEJM Table 2; Results Safety timing statement", kar_limits,
            attribution_context="Observed during ide-cel plus lymphodepletion; target-specific causality not adjudicated.",
            depletion="Direct normal plasma-cell-depletion count not extracted.",
            supportive_care="Immune globulin use described as common; exact count not extracted.",
            exposure_activity="ide-cel anti-BCMA CAR-T plus protocol lymphodepletion."),
        row("OBS-KARMMA-INFECTION-ANY", "KARMMA-INFECTION", kar, kar_cohort, kar_src,
            "infection adverse events", "immune function; indirect treatment-emergent endpoint", 88, "reported_exact_count", 69,
            "69%", "source-reported rounded percentage", 128, "ide-cel infused patients", "any_grade", kar_context,
            kar_window, "NEJM Results Safety; Table 2", kar_limits,
            attribution_context="Infections occurred during ide-cel treatment and protocol context; target-specific causality not adjudicated.",
            supportive_care="Antimicrobials and immune globulin described as common; exact counts not extracted.",
            exposure_activity="ide-cel anti-BCMA CAR-T plus protocol lymphodepletion."),
        row("OBS-KARMMA-INFECTION-G34", "KARMMA-INFECTION", kar, kar_cohort, kar_src,
            "infection adverse events", "immune function; indirect treatment-emergent endpoint", 28, "reported_exact_count", 22,
            "22%", "source-reported rounded percentage", 128, "ide-cel infused patients", "grade_3_or_4", kar_context,
            kar_window, "NEJM Results Safety; Table 2", kar_limits,
            attribution_context="Infections occurred during ide-cel treatment and protocol context; target-specific causality not adjudicated.",
            supportive_care="Antimicrobials and immune globulin described as common; exact counts not extracted.",
            exposure_activity="ide-cel anti-BCMA CAR-T plus protocol lymphodepletion."),
        row("OBS-KARMMA-NEUTROPENIA-ANY", "KARMMA-NEUTROPENIA", kar, kar_cohort, kar_src,
            "neutropenia adverse event", "marrow/hematologic function", 117, "reported_exact_count", 91, "91%",
            "source-reported rounded percentage", 128, "ide-cel infused patients", "any_grade", kar_context,
            "Adverse-event follow-up to trial data cutoff; exact event-specific window not extracted.",
            "NEJM Table 2; Results Safety", "Competing treatment-context event; not a plasma-cell-depletion endpoint.",
            attribution_context="Authors state most grade 3/4 events were hematologic and at least partly related to pre-infusion lymphodepleting chemotherapy.",
            supportive_care="Growth-factor use described as common; exact count not extracted.",
            recovery="Among patients with persistent grade 3/4 neutropenia at month 1 (n=52), median time to grade 2 or lower was 1.9 months.",
            exposure_activity="ide-cel anti-BCMA CAR-T plus protocol lymphodepletion."),
        row("OBS-KARMMA-NEUTROPENIA-G34", "KARMMA-NEUTROPENIA", kar, kar_cohort, kar_src,
            "neutropenia adverse event", "marrow/hematologic function", 114, "reported_exact_count", 89, "89%",
            "source-reported rounded percentage", 128, "ide-cel infused patients", "grade_3_or_4", kar_context,
            "Adverse-event follow-up to trial data cutoff; exact event-specific window not extracted.",
            "NEJM Table 2; Results Safety", "Competing treatment-context event; not a plasma-cell-depletion endpoint.",
            attribution_context="Authors state most grade 3/4 events were hematologic and at least partly related to pre-infusion lymphodepleting chemotherapy.",
            supportive_care="Growth-factor use described as common; exact count not extracted.",
            recovery="Among patients with persistent grade 3/4 neutropenia at month 1 (n=52), median time to grade 2 or lower was 1.9 months.",
            exposure_activity="ide-cel anti-BCMA CAR-T plus protocol lymphodepletion."),
    ])

    tal = "ST-GPRC5D-MONUMENTAL1"
    tal_src = "S-GPRC5D-MONUMENTAL1"
    tal_context = "Talquetamab × CD3 bispecific antibody; phase 1/2 MonumenTAL-1. Dose-escalation DLT and RP2D safety cohorts are distinct."
    tal_limits = (
        "Safety events in a single-arm treated population do not isolate GPRC5D-mediated normal-tissue injury; source percentages may be rounded. "
        "No causal safety rate or population risk is inferred."
    )
    rows.append(row("OBS-TAL-DLT-RASH-135-QW", "TAL-DLT-RASH-135-QW", tal, "CO-TAL-DOSE-ESCALATION-135-QW", tal_src,
        "dose-limiting grade 3 rash", "skin", 1, "reported_exact_count", "NR", "NR", "not reported",
        "NR", "dose escalation; denominator not recovered", "grade_3",
        "One grade 3 maculopapular rash dose-limiting toxicity at 135 µg/kg weekly after the second dose; rash began resolving after treatment discontinuation.",
        "Dose-escalation; exact patient-level follow-up not extracted.",
        "NEJM Results Safety, dose-limiting toxic effects paragraph",
        "Dose escalation denominator and full DLT rule/version not retrieved.",
        protocol_DLT="yes; one event at this schedule-specific dose level; denominator NR",
        attribution_context="Protocol DLT classification; investigator causality beyond the DLT label not independently adjudicated.",
        recovery="Source reports rash began resolving after treatment discontinuation; subsequent course not extracted.",
        exposure_activity="Talquetamab 135 µg/kg weekly; dose escalation."))
    for schedule, schedule_text, cohort, dose_count, followup, row_id in (
        ("QW", "800 µg/kg weekly dose-escalation cohort", "CO-TAL-DOSE-ESCALATION-800-QW", "NR", "dose escalation; denominator not recovered", "QW"),
        ("Q2W", "800 µg/kg every-other-week dose-escalation cohort", "CO-TAL-DOSE-ESCALATION-800-Q2W", "NR", "dose escalation; denominator not recovered", "Q2W"),
    ):
        rows.append(row(f"OBS-TAL-DLT-RASH-800-{row_id}", f"TAL-DLT-RASH-800-{schedule}", tal, cohort, tal_src,
            "dose-limiting grade 3 rash", "skin", 1, "reported_exact_count", "NR", "NR", "not reported",
            dose_count, followup, "grade_3", f"One grade 3 rash dose-limiting toxicity after the first {schedule_text}; resolved and treatment resumed after skipped dose(s).",
            "Dose-escalation; event followed first scheduled dose; exact patient-level follow-up not extracted.",
            "NEJM Results Safety, dose-limiting toxic effects paragraph; primary abstract",
            "Dose escalation denominator and full DLT rule/version not retrieved; a later protocol rule excluded an initial glucocorticoid-responsive grade 3 rash resolving within 7 days.",
            protocol_DLT="yes; one event in this schedule-specific dose level; denominator NR",
            attribution_context="Protocol DLT classification; investigator causality beyond the DLT label not independently adjudicated.",
            supportive_care="Glucocorticoids described for rash management; exact dose-level patient care not extracted.",
            recovery="Source reports rash resolution and treatment resumption after skipping dose(s).",
            exposure_activity=f"Talquetamab 800 µg/kg {schedule}; dose escalation."))

    rows.extend([
        row("OBS-TAL-SKIN-405-ANY", "TAL-SKIN-405", tal, "CO-TAL-RP2D-405-QW", tal_src,
            "skin-related adverse events", "skin", "NR", "NR_not_extracted", 67, "67%", "source-reported rounded percentage",
            30, "patients treated subcutaneously at 405 µg/kg weekly RP2D", "any_grade", tal_context,
            "Median follow-up 11.7 months; event-specific onset/window is not given in extracted text.",
            "NEJM indexed primary Results safety text; abstract", tal_limits,
            attribution_context="Skin events were reported during treatment; incidence alone does not establish target-specific injury.",
            supportive_care="Proactive topical and short-course oral glucocorticoids described; exact cohort use not extracted.",
            recovery="Skin/nail timing and resolution referred to supplement; supplement not retrieved.",
            exposure_activity="Talquetamab 405 µg/kg weekly subcutaneous RP2D."),
        row("OBS-TAL-RASH405-G34", "TAL-RASH-405", tal, "CO-TAL-RP2D-405-QW", tal_src,
            "rash-related adverse events", "skin", 0, "source_explicit_no_reported_events", 0, "0%", "source-explicit zero",
            30, "patients treated subcutaneously at 405 µg/kg weekly RP2D", "grade_3_or_4_rash", tal_context,
            "Median follow-up 11.7 months; event-specific onset/window is not given in extracted text.",
            "NEJM indexed primary Results safety text; abstract", tal_limits,
            attribution_context="The source states grade 3 rash-related events were not reported in this cohort.",
            exposure_activity="Talquetamab 405 µg/kg weekly subcutaneous RP2D."),
        row("OBS-TAL-SKIN-800-ANY", "TAL-SKIN-800", tal, "CO-TAL-RP2D-800-Q2W", tal_src,
            "skin-related adverse events", "skin", "NR", "NR_not_extracted", 70, "70%", "source-reported rounded percentage",
            44, "patients treated subcutaneously at 800 µg/kg every other week RP2D", "any_grade", tal_context,
            "Median follow-up 4.2 months; event-specific onset/window is not given in extracted text.",
            "NEJM indexed primary Results safety text; abstract", tal_limits,
            attribution_context="Skin events were reported during treatment; incidence alone does not establish target-specific injury.",
            supportive_care="Proactive topical and short-course oral glucocorticoids described; exact cohort use not extracted.",
            recovery="Skin/nail timing and resolution referred to supplement; supplement not retrieved.",
            exposure_activity="Talquetamab 800 µg/kg every other week subcutaneous RP2D."),
        row("OBS-TAL-RASH800-G3", "TAL-RASH-800", tal, "CO-TAL-RP2D-800-Q2W", tal_src,
            "rash-related adverse events", "skin", "NR", "NR_not_extracted", 16, "16%", "source-reported rounded percentage",
            44, "patients treated subcutaneously at 800 µg/kg every other week RP2D", "grade_3_rash", tal_context,
            "Median follow-up 4.2 months; event-specific onset/window is not given in extracted text.",
            "NEJM indexed primary Results safety text; abstract", tal_limits,
            attribution_context="The source reports grade 3 rash-related events in 16%; exact numerator was not extracted and is not reconstructed.",
            supportive_care="Grade 3 rashes responded to oral and topical glucocorticoids.",
            recovery="All but one patient with grade 3 rash resumed treatment after resolution.",
            exposure_activity="Talquetamab 800 µg/kg every other week subcutaneous RP2D."),
        row("OBS-TAL-NAIL-405", "TAL-NAIL", tal, "CO-TAL-RP2D-405-QW", tal_src,
            "nail-related adverse events", "nails", "NR", "NR_not_extracted", 57, "57%", "source-reported rounded percentage",
            30, "patients treated subcutaneously at 405 µg/kg weekly RP2D", "grade_not_reported", tal_context,
            "Median follow-up 11.7 months; event-specific onset/window is not given in extracted text.",
            "NEJM indexed primary Results safety text", tal_limits,
            attribution_context="Authors describe nail toxicity as probably related to on-target activity; this is an author mechanistic interpretation, not independent patient-level causal adjudication.",
            exposure_activity="Talquetamab 405 µg/kg weekly subcutaneous RP2D."),
        row("OBS-TAL-NAIL-800", "TAL-NAIL", tal, "CO-TAL-RP2D-800-Q2W", tal_src,
            "nail-related adverse events", "nails", "NR", "NR_not_extracted", 27, "27%", "source-reported rounded percentage",
            44, "patients treated subcutaneously at 800 µg/kg every other week RP2D", "grade_not_reported", tal_context,
            "Median follow-up 4.2 months; event-specific onset/window is not given in extracted text.",
            "NEJM indexed primary Results safety text", tal_limits,
            attribution_context="Authors describe nail toxicity as probably related to on-target activity; this is an author mechanistic interpretation, not independent patient-level causal adjudication.",
            exposure_activity="Talquetamab 800 µg/kg every other week subcutaneous RP2D."),
        row("OBS-TAL-ORAL-405", "TAL-ORAL", tal, "CO-TAL-RP2D-405-QW", tal_src,
            "oral-related adverse events (dysgeusia, dry mouth, dysphagia)", "oral mucosa/taste function", "NR", "NR_not_extracted",
            "NR", "NR", "not extracted", 30, "patients treated subcutaneously at 405 µg/kg weekly RP2D", "grade_1_or_2",
            tal_context, "Median follow-up 11.7 months; event-specific onset/window is not given in extracted text.",
            "NEJM indexed primary Results safety text", tal_limits,
            attribution_context="The source describes oral-related events as frequent and grade 1/2; event-specific rates and counts were not extracted.",
            exposure_activity="Talquetamab 405 µg/kg weekly subcutaneous RP2D."),
        row("OBS-TAL-ORAL-800", "TAL-ORAL", tal, "CO-TAL-RP2D-800-Q2W", tal_src,
            "oral-related adverse events (dysgeusia, dry mouth, dysphagia)", "oral mucosa/taste function", "NR", "NR_not_extracted",
            "NR", "NR", "not extracted", 44, "patients treated subcutaneously at 800 µg/kg every other week RP2D", "grade_1_or_2",
            tal_context, "Median follow-up 4.2 months; event-specific onset/window is not given in extracted text.",
            "NEJM indexed primary Results safety text", tal_limits,
            attribution_context="The source describes oral-related events as frequent and grade 1/2; event-specific rates and counts were not extracted.",
            exposure_activity="Talquetamab 800 µg/kg every other week subcutaneous RP2D."),
        row("OBS-TAL-DISCONTINUATION-800", "TAL-DISCONTINUATION", tal, "CO-TAL-RP2D-800-Q2W", tal_src,
            "adverse-event-related treatment discontinuation", "system-level treatment outcome", 1, "reported_exact_count",
            "NR", "NR", "not reported", 44, "patients treated subcutaneously at 800 µg/kg every other week RP2D",
            "any_grade", tal_context, "Through data cutoff; median follow-up 4.2 months; exact discontinuation time not extracted.",
            "NEJM indexed primary Results safety text; abstract", tal_limits,
            attribution_context="Investigators did not consider the single reported discontinuation related to talquetamab.",
            exposure_activity="Talquetamab 800 µg/kg every other week subcutaneous RP2D."),
    ])

    # Validate normalized fields: exact event rows, patient denominators, and no reconstructed counts.
    ids = [str(item["observation_id"]) for item in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate normalized observation IDs")
    for item in rows:
        count = item["event_count"]
        if count != "NR" and (not isinstance(count, int) or count < 0):
            raise ValueError(f"Event count is not a nonnegative integer or NR: {item['observation_id']}")
        if item["denominator_unit"] != "patients":
            raise ValueError(f"Denominator unit must be patients: {item['observation_id']}")
        if item["observation_state"] not in OBSERVATION_STATES:
            raise ValueError(f"Invalid endpoint observation_state in {item['observation_id']}: {item['observation_state']}")
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "source_mixed_draft_sha256": sha256(BACKUP),
        "normalized_observations_sha256": sha256(OUTPUT),
        "source_draft_rows": len(old),
        "normalized_event_grade_rows": len(rows),
        "normalization_rule": "one event, cohort/schedule, and grade row; integer patient count or NR; source percentage in a separate field; denominator_unit is patients; no count inferred from a rounded percentage",
        "observation_state_contract": {
            "allowed_values": sorted(OBSERVATION_STATES),
            "assigned_value": "indeterminate",
            "reason": "reported trial adverse-event counts are preserved descriptively, but endpoint qualification and target attribution have not been independently adjudicated; zero or unavailable grade counts do not establish adequate negative observation",
            "validated": True,
        },
        "dlt_note": "Indexed primary Results text reports four dose-limiting toxicities in the relevant dose-escalation passage: one lipase increase and three grade 3 rash cases. This target-focused table normalizes the three rash cases in separate 135 µg/kg weekly, 800 µg/kg weekly, and 800 µg/kg every-other-week rows; dose-level denominators are NR. The lipase DLT is outside this normal-tissue event table.",
        "source_access_note": "Indexed primary NEJM Results and Table 2 text supplied the bounded extracted results; direct article page opens returned HTTP 403. Table 2 exact event counts were retained; supplementary tables were not retrieved. Review remains pending.",
        "labels_predictions_and_clinical_AUC_changed": False,
        "review_status": "pending independent clinical review",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
