"""Export curated evidence, source-only forms, and auditable coverage. No scoring."""
import csv
import hashlib
import json
from pathlib import Path
import random
import re

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
BASE = PARENT.parent


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    original_manifest = json.loads((PARENT / "before_hashes.json").read_text())
    old_changes = {r["file"].replace("\\", "/") for r in json.loads(
        (PARENT / "document_changes.json").read_text())}
    preserved = []
    for rel, expected in original_manifest.items():
        if rel.replace("\\", "/") not in old_changes:
            assert digest(BASE / rel) == expected, rel
            preserved.append(rel)

    data = json.loads((HERE / "curated_evidence.json").read_text(encoding="utf-8"))
    rows = data["records"]
    assert len(rows) == len({r["id"] for r in rows}) == 15
    assert len({r["pmid"] for r in rows}) == 12
    assert len({r["gene"] for r in rows}) == 9
    allowed = {"direct_support", "author_probable", "author_possible",
               "alternative_attribution", "unresolved", "no_event_observed"}
    index = read_tsv(PARENT / "agent_evidence_index.tsv")
    flags = ["already_in_original_table", "already_in_validation_table",
             "already_in_locked_predictions"]
    for row in rows:
        assert row["attribution"] in allowed
        assert re.fullmatch(r"\d{8}", row["pmid"])
        assert row["source_url"].startswith("https://")
        assert all(row.values())
        target_rows = [r for r in index if r["gene_or_target_string"] == row["gene"]]
        assert target_rows, row["gene"]
        links = [r["evidence_id"] for r in target_rows
                 if row["pmid"] in re.findall(r"\d{8}", r["pmid_or_url"])]
        row["prior_index_ids"] = ";".join(links) or "No exact PMID match; target already in source pool"
        for flag in flags:
            values = {r[flag] for r in target_rows}
            assert len(values) == 1, (row["gene"], flag, values)
            row[flag] = values.pop()
        row["final_endpoint"] = "PENDING_REVIEW"
        row["extraction_date"] = data["date"]
        row["independence"] = "Not established; purposive review with prior outcome knowledge"
    write_tsv(HERE / "clinical_evidence.tsv", rows)

    studies = {}
    for row in rows:
        if row["pmid"] not in studies:
            studies[row["pmid"]] = {
                "study_id": "PMID:" + row["pmid"], "gene": row["gene"],
                "agent": row["agent"], "primary_source": row["source_url"],
                "verified_trial_and_cohort": "", "event_and_severity": "",
                "protocol_DLT_status": "", "normal_tissue_attribution": "",
                "negative_observation_adequacy": "", "antigen_mapping_issue": "",
                "overlapping_publications": "", "source_location": "",
                "endpoint_judgment_and_reason": "", "reviewer_name_date": ""}
    form = list(studies.values())
    random.Random(20261002).shuffle(form)
    write_tsv(HERE / "study_source_form.tsv", form)
    loaded = read_tsv(HERE / "study_source_form.tsv")
    assert len(loaded) == 12
    assert all(not r["reviewer_name_date"] and not r["endpoint_judgment_and_reason"] for r in loaded)
    assert not ({"attribution", "event", "score", "top_share", "historical_outcome"} & set(loaded[0]))

    reviewed_keys = {(r["gene"], r["pmid"]) for r in rows}
    queue = []
    for item in index:
        matched = [(item["gene_or_target_string"], pmid) in reviewed_keys
                   for pmid in re.findall(r"\d{8}", item["pmid_or_url"])]
        queue.append({"evidence_id": item["evidence_id"],
                      "gene": item["gene_or_target_string"],
                      "source": item["pmid_or_url"],
                      "source_review_status": "Partial extraction in pass 2" if any(matched) else "Not extracted in pass 2",
                      "endpoint_status": "PENDING_REVIEW",
                      "next_step": "Resolve full-text gaps and independent judgment" if any(matched) else "Verify source and extract separate event fields"})
    write_tsv(HERE / "source_queue.tsv", queue)
    gene_groups = {flag: sorted({r["gene"] for r in rows if r[flag] == "True"}) for flag in flags}
    result = {
        "evidence_units": len(rows), "publications": len(studies),
        "target_genes": len({r["gene"] for r in rows}),
        "prior_pool_rows_partially_matched": sum(r["source_review_status"].startswith("Partial") for r in queue),
        "prior_pool_rows_total": len(queue), "target_overlap": gene_groups,
        "historical_files_hash_preserved": len(preserved),
        "final_endpoint_judgments": 0, "new_model_scores_or_performance_estimates": 0,
        "source_only_form_rows": len(loaded), "reviewer_fields_blank": True,
        "limits": "Schema, linkage and hash checks only; not independent verification of clinical truth. No pooled patient counts."
    }
    (HERE / "coverage_and_checks.json").write_text(json.dumps(result, indent=2) + "\n")
    files = [p for p in HERE.iterdir() if p.is_file() and p.name != "artifact_hashes.json"]
    (HERE / "artifact_hashes.json").write_text(json.dumps({p.name: digest(p) for p in sorted(files)}, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
