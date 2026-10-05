# Round 7: searching for the missing data (Sept 29, 2026)

Question from Dylan: could a better rule exist that we cannot find because we lack data? Which data, and can we get it? This note reports the search. **No rule was scored and no label was changed.** All files are in `ev_capture_data/round7_agent_data/` (hashes in `round7_data.hashes`).

## What was found, in the order searched

| # | Data | Result | Rows | Main source | Quality / gaps |
|---|---|---|---|---|---|
| 1 | Drug type, dose and toxicity behind each existing label | **Found for all 30 panel targets and all 12 validation/hedged targets** | 55 approved-label rows, 111 trial-paper rows | FDA labels (openFDA), PubMed and PMC full text | some fields NOT_FOUND; many trial rows are abstract-only (flagged in the notes column); CD44v6 CAR-T has no published clinical results |
| 2 | More labelled drugs | **Found: 77 genes total, 35 new** (not in the panel or validation set); 8 drug formats | included in the row counts above | FDA approval lists, review papers, primary trial papers | Attribution (on-target or not) is stated in about 40% of trial rows and 15% of approved-label rows; 42 development statuses were reset to NOT_FOUND because they were not verified |
| 3 | Graded, organ-level toxicity | **Found**: label dose-modification organs and grade 3-4 rates (approved); DLT event, organ, grade, dose (trials) | same tables | same | Approved-label `max_severity` is almost always "fatal", so use the organ and rate columns instead |
| 4 | Tumor antigen expression | **Found for solid tumors**: TCGA median and quartiles per cancer, tumor and tumor-adjacent normal, all 42 genes | 823 rows | UALCAN pan-cancer page (fetched), HPA (category only), cBioPortal (SCLC) | **Not obtainable:** myeloma, Hodgkin, B-ALL, CLL/FL, CTCL, neuroendocrine, Merkel cell, synovial sarcoma. Hand-transcribed; 2 rows re-checked against the source and matched. Units differ from the project's normal-tissue nCPM (see routes file) |
| 5 | Protein / surface accessibility | **Partly found**: curated membrane topology and extracellular length from UniProt for all 42 genes | 42 rows | UniProt (fetched) | HPA subcellular location is unreliable for this purpose (DLL3 "nucleoplasm", F3 "vesicles"). Tumor protein (IHC in cancer) not obtainable by the tools available |

Blocked routes: ClinicalTrials.gov (robots.txt), GDC expression values (POST only), DepMap (bot wall), most full-text pages behind reCAPTCHA. None was bypassed.

## The main finding: the existing labels are less uniform than a per-target label implies

The per-agent search read the primary papers. Its notes (agent reading; **needs a human to verify against the papers before anything is changed**) flag these panel labels:

| target (label) | what the sources say |
|---|---|
| CEACAM5 (toxic) | The Parkhurst 2011 paper is a **TCR-T**, not a CAR-T, and the authors do attribute severe colitis to CEA. A CEA CAR-T (Zhang 2017) was tolerated to 1e8 cells/kg. |
| ERBB2 (toxic) | The fatal 2010 case is one patient, and the authors only postulate on-target lung damage. A HER2 CAR-T in sarcoma (Ahmed 2015) had no DLT from 1e4 to 1e8 cells/m2. |
| EGFR (toxic) | One CAR-T study calls skin and lung target-mediated effects grade 1-2 only; another (1-3e6 cells/kg) saw no on-target toxicity beyond a grade 1 rash. |
| EPCAM (toxic) | Adecatumumab was well tolerated with no maximum tolerated dose. Catumaxomab toxicity is attributed mainly to immune activation. |
| NCAM1 (toxic) | Lorvotuzumab's neuropathy is not attributed to CD56 by the authors. |
| CD47 (toxic) | On-target anaemia was real but not dose-limiting (no MTD to 45 mg/kg). The sheet's "programme discontinued 2024" was not re-checked. |
| FAP (toxic) | A human intrapleural FAP CAR-T (1e6 cells, low dose) reported no toxicity; the basis of the "toxic" label (marrow stromal damage) was not re-checked. |
| CD44 (toxic) | Rests on the bivatuzumab ADC (fatal skin necrolysis, authors say "most probably" on-target). No human CD44v6 CAR-T results exist. |
| CA9 (toxic) | Consistent: authors explicitly call the liver toxicity on-target. |
| MSLN (safe) | Consistent for mesothelin, but anetumab had dose-limiting eye and nerve toxicity attributed to the ADC class. |
| GPC3, CTAG1B (safe) | Consistent for the agents studied; the same lab's MAGE-A3 TCRs killed patients (attributed to cross-reactivity, so "safe for the target" does not extend to other TCRs). |
| GUCY2C, PSCA, MAGEA12 (validation) | GUCY2C diarrhoea is on-target but not dose-limiting; PSCA cystitis tracked cyclophosphamide dose; MAGEA12 attribution is disputed. |

**Also:** in the approved-drug labels, no organ toxicity is attributed to target expression in normal tissue. Only the CD19 CAR-T labels use the phrase "on-target" (B-cell aplasia). So regulatory text cannot supply on-target attribution; it has to come from trial papers.

**What this means:** the same target can be tolerated with one drug format, route or dose and toxic with another (HER2, EGFR, CEACAM5, DLL3, MSLN). This is the modality and dose effect we suspected. It suggests that part of the remaining error in the frozen rule is label noise, not a missing rule.

## Other findings worth knowing

- **Tumor level is often not higher than normal.** TCGA tumor/adjacent-normal ratios: PSCA stomach 0.05, PSCA prostate 0.39, TACSTD2 breast 0.69, GUCY2C colon 0.74, CLDN18 stomach 0.76, CEACAM5 colon 0.84; above 1: CA9 kidney about 2558 (normal near zero), GPC3 liver 95, MSLN pancreas 108 (only 4 normals), FOLH1 prostate 3.3. TCGA has no normal for ovary, mesothelioma, DLBC or LAML. Several normal groups are tiny (F3 cervix n=3).
- **Membrane structure (UniProt):** 26 single-pass, 7 multi-pass, 6 GPI-anchored, and 3 with no membrane annotation (MAGEA4, MAGEA12, CTAG1B). Those three are intracellular antigens seen as peptide-HLA by TCR-T, so a surface-expression logic does not apply to them the same way.
- **Route recipes** (UALCAN needs legacy gene symbols, for example PVRL4 for NECTIN4) are in `tumor_expression_routes.md`.

## Verification done on this data

- All files re-parsed: consistent column counts, no tabs inside fields.
- 8 randomly chosen trial-row PubMed IDs re-fetched: all 8 titles match the agent named in the row.
- PTK7 (cofetuzumab) row compared with the abstract: dose range, DLT and recommended dose match.
- 2 tumor-expression rows (PSCA in stomach and prostate, normal median and quartiles) re-fetched from UALCAN: match to three decimals.
- **Not verified:** the rest. Treat everything as an agent-collected first pass. Rows marked abstract-only have unchecked fields marked NOT_FOUND.

## Proposed next analyses (NOT run; each would be precommitted with hashes first, and each adds to the count of alternatives, now 74)

1. **Strict-label sensitivity.** Define in advance a "strict" label set: toxic only when the authors of a primary paper attribute a dose-limiting toxicity to the target on a named normal tissue; safe only when tolerated with no such attribution. Score the frozen rule on it. Must be written down before any score is seen, and the frozen labels stay primary. Needs a human (Professor Wideman) to check the label calls.
2. **Agent-level table.** One row per agent with modality as a variable; test whether modality explains any of the errors.
3. **Tumor/normal fold-change form** of the relative-exposure rule, using UALCAN ratios for solid tumors with normals only (about 20 genes). Low power; descriptive.
4. **Scale-up.** Extend the labelled set toward the 150-300 targets needed to detect a 0.04 improvement, always recording modality and dose.

## Caveats

Everything here was gathered by AI search agents from public sources; the project rules (Society for Science: no generative-AI-written submitted text) apply to the paper and poster, not to data gathering, but keep the AI-use log current. Do not change any frozen label on the strength of this note.
