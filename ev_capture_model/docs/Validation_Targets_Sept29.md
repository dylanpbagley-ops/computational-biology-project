# New labelled targets and label-consistency findings — Sept 29, 2026

For Dylan and for outside review (Professor Wideman). Labels are frozen in `target_labels_validation.py` (hashes in `validation_precommitment.md.hashes`). All sources were opened by the searchers; two of them were abstracts only (flagged). Nothing here was chosen using expression data.

## The frozen rule for labelling (same as the main panel)
safe = approval or documented clinical tolerability, no dose-limiting on-target normal-tissue toxicity. toxic = documented severe, dose-limiting or fatal toxicity attributed by the authors to the target on a normal tissue. Excluded: non-protein antigens, early-only evidence, same gene with both an approval and a toxicity record, blockers rather than cell-killers.

## Primary validation set (pre-specified)
| gene | agent | outcome | confidence | why | source |
|---|---|---|---|---|---|
| MAGEA4 | afamitresgene autoleucel (TCR-T) | safe | high | FDA approval 2024; adverse events are lymphodepletion cytopenias and CRS; no treatment-related deaths | PMID 38554725 |
| SSTR2 | 177Lu-dotatate | safe | medium-high | approved 2017/2018; NETTER-1 no renal toxicity, grade 3-4 neutropenia 1%; kidney dose is reabsorption | PMID 28076709 |
| CCR4 | mogamulizumab | safe | medium | approved 2018; rash and Treg depletion manageable, no boxed warning | PMID 30100375 |
| CTAG1B (NY-ESO-1) | NY-ESO-1 TCR-T | safe | medium | well tolerated in several series; not approved. **No measured expression in the reference, so it could only be scored by a post-hoc convention** | PMID 26193344, 29891538, 39836945 |
| IL2RA (CD25) | denileukin diftitox | safe | medium | approved 2024; capillary-leak boxed warning is a toxin-class effect not attributed to the target | PMID 39700456 |
| GUCY2C | IM96 CAR-T; PF-07062119 | toxic | medium-high | 55% grade 3 diarrhoea called a notable on-target off-tumour toxicity; TCE dose-limiting colitis/diarrhoea plus CRS. TCE source is an abstract | PMID 42241671, 42696650 |
| PSCA | PSCA CAR-T | toxic | medium | two dose-limiting grade 3 cystitis events; authors say on-target with a cyclophosphamide contribution | PMID 38867077 |

Result of the frozen rule: separation 0.875 on 4 v 2 scorable targets (exact p = 0.133), 0.900 with the post-hoc absent-gene convention (5 v 2, p = 0.095). Too small to be evidence by itself.

## Hedged toxic (reported separately, never pooled)
TFRC (CX-2029; authors: multifactorial, MMAE overlap), SLC44A4 (ASG-5ME; on- and off-target), MLANA (MART-1 TCR-T melanocyte destruction; no severity data seen; abstract only), MAGEA12 (anti-MAGE-A3 TCR; two deaths with brain necrosis; the culprit gene is disputed; **no measured expression**), EPHA2 (MEDI-547; attribution unresolved). Excluded: PMEL (approval plus on-target toxicity on one gene), SLC6A2 (transporter, product withdrawn), CD276, PRAME (too early).

## Label-consistency issues in the MAIN panel (for the outside reviewer; labels were NOT changed)
- **NECTIN4 is labelled safe but enfortumab vedotin carries a boxed warning for fatal SJS/TEN.** F3 (tisotumab) and FOLR1 (mirvetuximab) are labelled toxic on an ocular boxed warning, and a BCMA ADC (belantamab) carries the same ocular boxed warning while TNFRSF17 is labelled safe from CAR-T/bispecific evidence. Labels are agent-specific verdicts, not gene-level facts; say so in the methods.
- Boxed-warning text does not say whether a toxicity is on-target or payload-driven. Wording taken from web-extracted label text; check the exact wording against the current labels before publishing.
- Sensitivity: flipping NECTIN4 to toxic gives 0.926; F3 and FOLR1 to safe gives 0.884; both 0.890; TNFRSF17 to toxic 0.838. All p < 0.001.
- Not verified this session (from memory only): these panel genes have no approved cytotoxic agent I know of: MSLN, GPC3, CEACAM5, CA9, CD47, CD44, NCAM1, FAP, FCER2, GPNMB, CD40, CD80, HLA-DRA, SDC1, MUC1, MUC16, SLC34A2, CD74, AXL, PTPRC. If any is labelled safe because of an approval, check it.

## Why the set is so small
Review-based searches mostly recycled genes already in the panel; clinical papers almost always hedge on-target attribution (payload, lymphodepletion, CRS, off-target cross-reactivity). Cleanly attributable toxic examples are rare, so a properly powered validation (roughly 150 to 300 labelled targets to detect a 0.04 improvement) is not possible from the literature. The OSF-registered prospective predictions remain the real out-of-sample test.
