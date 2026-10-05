# Computational biology research project

Private research workspace for an existing antigen-distribution / cancer-immunotherapy project and its subsequent computational extension investigations. Code, research notes, derived results, figures, and the latest compact processed input set are preserved at their original relative paths.

## Start here

1. Read [`PROJECT_STATUS_Oct2.md`](PROJECT_STATUS_Oct2.md) for the original project's assessment and validity boundaries.
2. Read [`opportunity_search_Oct4/STATUS.md`](opportunity_search_Oct4/STATUS.md) for the newest extension checkpoint.
3. Use [`WHERE-EVERYTHING-IS.md`](WHERE-EVERYTHING-IS.md) to locate original code and results. This older map contains historical claims; newer assessments and decision records qualify them.
4. The original implementation and detailed results are in [`ev_capture_model/`](ev_capture_model/).

The original clinical result is a selected-panel retrospective association; software checks do not establish general clinical prediction or clinical safety. Preserve original labels and locked predictions. No substantial extension has yet passed the project-selection gate.

## Latest experimental-data pilot

[`opportunity_search_Oct4/cycle4_functional/decision.md`](opportunity_search_Oct4/cycle4_functional/decision.md) records a bounded comparison of published mouse tumor CRISPR-selection data. It produced five exploratory cross-study flags. MED16 is the only lead advanced to a conditional independent functional-replication/data-access gate. No new mechanism, causal immune-modality effect, human effect, or clinical result is established.

The next handoff is [`docs/superpowers/plans/2026-10-04-med16-replication-decision-gate.md`](docs/superpowers/plans/2026-10-04-med16-replication-decision-gate.md). It was prepared but has not been executed.

## Contents and data

The repository includes original code/tests, Markdown documentation and plans, derived JSON/CSV/TSV results, figures, and the verified cycle4 processed guide-count/workbook input set (about 21 MB). It excludes multi-gigabyte raw datasets, installed dependency folders, caches, downloaded publication PDFs, large archives, and personal fair forms. Excluded files remain in the local workspace.

This is a research snapshot, not a packaged application. Some original analyses require excluded datasets and additional dependencies. Source inventories and data manifests record provenance; do not assume every historical entry point runs from a fresh clone. Historical preservation validators may reference intentionally excluded files and are not portable installation tests.

For the latest pilot helper checks, use Python 3.12 and [`requirements-pilot.txt`](requirements-pilot.txt):

```console
python -m pip install -r requirements-pilot.txt
python opportunity_search_Oct4/cycle4_functional/pilot/test_transport.py
```

These requirements describe the verified pilot environment only, not dependencies for the entire historical project. The published-input extraction script writes regenerated results; run it in a separate copy if preserving the frozen snapshot. The original protocol's checksum-description inconsistency is documented in [`pilot/hash_definition_note.md`](opportunity_search_Oct4/cycle4_functional/pilot/hash_definition_note.md).

## Research assistance

Use the current status, decision records, source tables, and actual data before proposing work. Keep validity repair separate from substantive contribution. Do not infer novelty from a missing search hit, count guides as independent tumors, reinterpret absent data as zero effects, or rewrite negative findings as successes. The student writes submitted paper, abstract, poster, and research-plan prose; assistance provides evidence, analyses, figures, outlines, and feedback.

Keep this repository private. No collaborator or third-party application access has been granted by creating it. Third-party data and materials retain their original provenance and terms; this repository does not grant a new license over them.
