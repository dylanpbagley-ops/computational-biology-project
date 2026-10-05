"""Assemble original-panel source coverage without assigning outcomes or scores."""
import csv
import hashlib
import json
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
BASE = PARENT.parent


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    original = [r for r in read_tsv(PARENT / "operational_labels_44_Oct2.tsv")
                if r["outcome"] in {"safe", "toxic"}]
    genes = {r["gene"] for r in original}
    assert len(original) == len(genes) == 30
    additions = read_tsv(HERE / "extraction.tsv")
    prior = read_tsv(PARENT / "clinical_review_pass2/clinical_evidence.tsv")
    assert len(additions) == len({r["gene"] for r in additions}) == 25
    assert all(set(r) == set(additions[0]) and all(r.values()) for r in additions)
    prior_original = [r for r in prior if r["gene"] in genes]
    assert not ({r["gene"] for r in prior_original} & {r["gene"] for r in additions})
    assert {r["gene"] for r in prior_original + additions} == genes

    factual = {"GPRC5D", "FAP", "FOLR1", "CEACAM5"}
    attribution = {"CA9", "CD44", "CD47", "NECTIN4"}
    context = {"MSLN", "ERBB2", "EGFR", "EPCAM", "FOLH1"}
    sources = []
    form = []
    for gene in sorted(genes):
        new = [r for r in additions if r["gene"] == gene]
        old = [r for r in prior_original if r["gene"] == gene]
        urls = sorted({r["primary_url"] for r in new} | {r["source_url"] for r in old})
        agents = sorted({r["agent"] for r in new + old})
        priority = ("Factual rationale correction" if gene in factual else
                    "Assess event severity and target attribution" if gene in attribution else
                    "Resolve context and negative-observation scope" if gene in context else
                    "Complete mechanism and observation-adequacy extraction")
        sources.append({"gene": gene, "review_priority": priority,
                        "source_pass": "pass3" if new else "pass2",
                        "source_urls": "; ".join(urls),
                        "coverage": "Representative source extraction only; not exhaustive",
                        "human_clinical_eligibility": "Cited evidence is preclinical" if gene == "FAP" else "Human report located; endpoint eligibility pending",
                        "final_outcome": "PENDING_REVIEW"})
        form.append({"gene": gene, "agent_or_agents": "; ".join(agents),
                     "primary_sources": "; ".join(urls), "verified_study_cohorts": "",
                     "human_or_preclinical": "", "event_grade_and_denominator": "",
                     "DLT_and_protocol_definition": "", "normal_tissue_attribution_and_evidence": "",
                     "intentional_depletion_and_consequences": "", "negative_observation_adequacy": "",
                     "antigen_gene_mapping": "", "overlap_with_other_reports": "",
                     "proposed_endpoint_and_reason": "", "reviewer_name_date": ""})
    write_tsv(HERE / "original_panel_coverage.tsv", sources)
    random.Random(20261003).shuffle(form)
    write_tsv(HERE / "source_only_review_form.tsv", form)
    assert all(not r["proposed_endpoint_and_reason"] and not r["reviewer_name_date"] for r in form)
    assert not {"outcome", "score", "top_share", "review_priority"} & set(form[0])

    # A compact view of what is still required, not an automated clinical judgment.
    work = [{"gene": r["gene"], "source": r["primary_url"],
             "remaining_work": r["remaining_work"], "status": "OPEN"} for r in additions]
    work += [{"gene": r["gene"], "source": r["source_urls"],
              "remaining_work": "Resolve pass 2 limitations and independent endpoint judgment",
              "status": "OPEN"} for r in sources if r["source_pass"] == "pass2"]
    write_tsv(HERE / "remaining_work.tsv", sorted(work, key=lambda r: r["gene"]))

    manifest = json.loads((PARENT / "before_hashes.json").read_text())
    modified = {r["file"].replace("\\", "/") for r in json.loads((PARENT / "document_changes.json").read_text())}
    protected = [rel for rel in manifest if rel.replace("\\", "/") not in modified]
    for rel in protected:
        assert digest(BASE / rel) == manifest[rel], rel
    prior_manifest = json.loads((PARENT / "clinical_review_pass2/artifact_hashes.json").read_text())
    for rel, expected in prior_manifest.items():
        assert digest(PARENT / "clinical_review_pass2" / rel) == expected, rel
    checks = {
        "new_target_source_records": len(additions), "original_targets_with_representative_source_coverage": len(sources),
        "original_targets_total": len(genes), "original_targets_with_human_report_located": sum(r["gene"] != "FAP" for r in sources),
        "preclinical_only_source_used_for_original_target": ["FAP"],
        "complete_clinical_adjudications": 0, "new_performance_estimates": 0,
        "review_form_rows": len(form), "all_review_judgments_blank": True,
        "protected_historical_files_unchanged": len(protected),
        "prior_pass_manifest_files_unchanged": len(prior_manifest),
        "scope": "File/schema/coverage checks; coverage is not comprehensive review or clinical validation."
    }
    (HERE / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    artifacts = [p for p in HERE.iterdir() if p.is_file() and p.name != "artifact_hashes.json"]
    (HERE / "artifact_hashes.json").write_text(json.dumps({p.name: digest(p) for p in sorted(artifacts)}, indent=2) + "\n")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
