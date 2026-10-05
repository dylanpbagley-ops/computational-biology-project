# Functional outcome extension selection — status

Last updated: 2026-10-04T07:19:44Z.

## Scope and runtime

Completed the authorized selection of exactly two candidate questions: antibiotic survival and HIV bnAb susceptibility/finite sampling. No third candidate, full project, pilot, raw-data download, old-file edit, or external action was performed. Runtime: Python 3.12.14, NumPy 2.3.5, read-only h5py 3.16.0 from the existing longitudinal dependency tree. Parent execution clock: start 2026-10-04T06:37:16Z; 60-minute eligibility checkpoint 07:37:16Z. Gate review remains within that window.

## Gate decisions

| Candidate | Exact question screened | Novelty gate | Data gate | Stage verdict |
|---|---|---|---|---|
| A — antibiotic survival | Does the pre-exposure fraction of cells that neither grow nor divide predict later ampicillin survivor frequency beyond a simple population model? | **Fail:** Umetani et al. (2025), Fig. 7B–C and Eq. 1, directly relates pre-exposure non-growing fraction to later survivor frequency using f=sP_n, estimating s from the late-stationary condition. | Not advanced after novelty failure. No extra source-data inspection or retrieval. | **reject** |
| B — HIV bnAb susceptibility | Would subsampling baseline single-bnAb susceptibility change a fixed pair-selection decision, and would the selected pair retain measured combination IC80 coverage on viruses withheld from selection? | **Unresolved / insufficiently differentiated:** Wagh et al. already rank combinations and compare Bliss–Hill predictions with measured combinations; Williamson et al. test held-out per-virus susceptibility prediction and synthetic sample-size performance. Neither directly reports finite-panel decision stability, but that remaining calculation is not yet a distinct biological contribution. | **Partial, not a pass:** Wagh Fig. 2 reports measured combination IC80 for 10 regimens on 20 deliberately selected viruses. Their broader 1,622-combination/200-virus results are model predictions. Williamson CATNAP data measure individual antibodies; Fig. 4 combination values are calculated from singles under additive/Bliss–Hill models. Exact row-level joins, donor/lineage grouping and representative independent evaluation were not reconstructed. | **reject for this stage** |

No candidate passed both novelty and data gates. No pilot, protocol, tests or figure was created. This is a bounded rejection of the current two proposals, not a conclusion that no future biological extension is possible.

## Current blocker and next actions

There is no nominated question here with a defensible incremental result and a verified independent measured-outcome evaluation. The next two actions are: (1) parent independently reviews the source/gate and preservation files; (2) Dylan decides whether to make a broader pivot in a later, separately scoped project or present the existing project as completed work without claiming a functional extension. Parent's role is source and integrity review; Dylan authors any submitted prose. No clinician, lab or external collaborator was contacted or assumed.

## Stage files

`selection.tsv` records both questions and gate status; `sources.tsv` gives DOI/URL, source sections and access limits; `data_manifest.json` distinguishes measured records, predictions and simulations; `decision.md` contains the no-selection recommendation and explanation points; `reproduction.md` lists exact verification commands; `preservation/` records baseline hashes and verification; `output_hashes.json` inventories this stage. No new empirical data or dependencies were downloaded. The old myeloma result and prior files remain untouched.
