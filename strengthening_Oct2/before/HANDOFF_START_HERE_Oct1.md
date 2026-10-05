# HANDOFF: START HERE (written Oct 1, 2026)

For a new AI assistant (ChatGPT) picking up Dylan's ISEF project after a switch from Claude. Read this file first, then do the "Verify before trusting" step (Part 8). It replaces `EV-Capture-CURRENT.md` as the entry point. That older file (Sept 10) is out of date: its "headline" and next steps describe the project as it stood before the Sept 26-29 work.

---

## 1. Who Dylan is and how to work with him
High school student, working solo, ISEF-level computational biology project. Beginner coder (Python basics).
- Plain language. Define a technical term the first time you use it. Concise, but never at the cost of leaving something out.
- Explain concepts in prose, not walls of bullets. Walk him through both the science and any debugging. If he runs code, tell him exactly what to type and what the output should look like.
- Stress-test before executing. Say "this could fail because X" now rather than in January. Do not reassure reflexively.
- He catches errors and expects them acknowledged directly.
- Verify claim by claim against real files and data before trusting any summary, including this one and your own first-pass output.
- If a search or attempt comes back empty or negative, report that plainly and ask how to proceed. Do not quietly work around it.
- If you notice yourself narrowing a broad search to the most exciting thing, stop and say so.
- You cannot see his computer. Ask him to run commands and paste the output. Never invent a file path, a number or a citation. Check the file or say you cannot.

## 2. Hard rule: AI may not write his submitted text
Society for Science (ISEF) bars generative AI from writing the submitted abstract, research plan, paper or poster text. Dylan writes all of that himself. You may give structure, outlines, figures, analysis, checks, and feedback on drafts he wrote. The AI-use disclosure log (`ev_capture_model/docs/AI_Use_Disclosure_Log.md`) must list ChatGPT use as well as the earlier Claude use. It is stale and needs Dylan's own update (his task, Part 6).

## 3. The project in one paragraph
**Claim:** normal-tissue data alone predicts which antigen targets cause dose-limiting on-target toxicity. The frozen rule is `score = top_share + 1[blood is the largest lineage family]`. `top_share` is the top cell type's share of a gene's expression across the 154 Human Protein Atlas (HPA) cell types; blood is the largest of ten lineage families.
**Panel:** 19 safe vs 11 toxic targets. Separation (AUC, the chance a random SAFE target outscores a random toxic one — higher score = safer) = **0.914**; confinement alone = 0.880. A ROC-hull theorem bounds what any monotone therapy can achieve. A lemma (results §40.2) shows `top_share` is the tightest worst-case guarantee on how many other cell types express the gene.
**Honest limits (stay in the paper):** 74 alternative rules were tried; labels are agent-specific verdicts and a few are weakly supported; no tumour data and no drug-type effect yet (future work).
Older Claim 1 (EV capture) and the CAR-T pivot are earlier phases, not the current headline.

## 4. Settled decisions: do not reopen
- In Sept 2026 Dylan decided to narrow the project to one core claim, cut peripheral material, and show only the most important limitations. Do not suggest widening it. Award-worthiness and real-world impact are different things; do not conflate them.
- The frozen rule stays primary. **No further rule search** (count stays 74; promotion needs separation >= 0.944 and family-wise p < 0.05).
- Frozen labels are never edited. The strict-label sensitivity analysis was deliberately skipped.
- Precommit and hash before any new analysis (pattern: `*_precommitment.md` plus `.hashes`).
- OSF predictions are already submitted. The locked list is `prospective_predictions.json` (2026-09-08): 52 targets, 14 tolerated, 18 toxic, 20 no-call. Never revise it.

## 5. State of the work as of Oct 1
Done: theorem and ROC hull (§29, with the corrected wording that confinement is empirical, not "the leading term"); label fragility (§40.1: fragility index 4, weak-label subsets keep p < 0.05 in 97.3%); breadth lemma (§40.2); locked-predictions check (§41: 0 forward verdicts, 5 pre-lock verdicts scored, 2 hits and 3 misses, below the precommitted floor of 8, so no conclusion); 7-target external validation (separation 0.875, p = 0.133, too small to count); paper skeleton rewritten around the one claim.
Verifier `verify_documented_numbers.py`: 221 checks passing. **§41 is not yet wired into the verifier**; `tests/test_locked_check.py` covers it instead.

## 6. What is next
Dylan's own tasks:
1. Update the AI-use log (Claude work through Oct 1, ChatGPT from now on).
2. Write the paper and poster text himself, using `ISEF_Paper_Skeleton.md` (Part 9 is his to-do list).
3. Send `Label_Review_Sheet.md` and `Validation_Targets_Sept29.md` to Prof. Jeremy Wideman (ASU) when reachable, for outside review of the labels.
4. Check the ISEF and fair paperwork deadlines.
Timeline: school fair end of Jan 2027, AzSEF March 2027, then Regeneron ISEF.
Optional, not requested: polish figures `fig_hull`, `fig_ladder`, `fig_forest`; prepare a CEACAM5 worked example; wire §41 into the verifier.

## 7. Where things are
Device root: `C:\Users\dylan.bagley\Desktop\outputs`. Docs in `ev_capture_model\docs\`, code in `ev_capture_model\src\aim13_scoping\`, tests in `ev_capture_model\tests\`.

| need | file |
|---|---|
| every number and method (authoritative) | `Aim13_Results_Two_Comparator_Ceiling.md` (§1-41; §40-41 are newest) |
| paper structure, strength table, numbers table, glossary, judge questions | `ISEF_Paper_Skeleton.md` |
| the Sept 29 decisions and what was cut | `Scope_Narrowing_Plan_Sept29.md` |
| rule definition and label rules | `Validation_Targets_Sept29.md`, `target_labels_frozen.tsv` |
| why no better rule was found | `Round7_Data_Search_Findings_Sept29.md` |
| comparison with field actions | `Framework_vs_Field_Actions_Sept26.md` |
| file map (Sept 9; code sections still valid) | `WHERE-EVERYTHING-IS.md` (root) |
| how the handoff system worked | `PROCESS.md` (root) |

Snapshots, not current (numbers differ): `Status_ISEF_Ratings_and_Expansion_v2.md` (Sept 8, still says "leading term" and quotes 19 v 10 / 0.868), `Full_Assessment_..._Sept8.md`, `Work_Completed_Sept9.md`, `EV-Capture-CURRENT.md`. If a snapshot disagrees with the results doc, the results doc wins.

## 8. Verify before trusting
Ask Dylan to run, in a terminal from `ev_capture_model\src\aim13_scoping`:
`python verify_documented_numbers.py` (expect 221 checks, all passing), then from `ev_capture_model`:
`python -m pytest -q tests\test_fragility_and_bound.py tests\test_locked_check.py` (expect 10 passing).
The full suite (about 195 tests) takes about 2.5 minutes and only needs running after code changes. After editing any script, run its own entry point as well as its tests.

## 9. Known traps
- Labels are agent-specific, not gene-level facts (e.g. NECTIN4 is "safe" although enfortumab vedotin has a fatal-skin boxed warning). Never describe them as gene facts.
- Never state a PMID or a clinical claim you have not opened. Two validation sources were abstracts only.
- Do not use `top_share` described as "the leading term" of the harm formula. That wording was corrected.
- Sensitivity numbers must come from unrounded `top_share` in `wholebody_v2.json`, not the rounded values in `flag_cleanup.json` (rounding creates false ties).

## 10. First message and standing instructions for ChatGPT
See `ChatGPT_First_Prompt.md` (Part A: standing instructions to paste into a ChatGPT Project or custom instructions; Part B: the first message to send).
