# Review of the implemented expansion angles, and two drift findings

**27 September 2026.** Written after implementing Angles 1-4 and then going back over all of
them looking for anything overlooked. Two findings here are more important than the review
itself and are stated first.

---

## FINDING 1 (act on this) -- the results document is stale relative to the 44-target panel

`Aim13_Results_Two_Comparator_Ceiling.md` still reports two numbers computed on the OLD
42-target panel. The September 10 session added ERBB3 and AXL to `target_labels.py` but the
change did not propagate to the document in these two places, and `wholebody_v2.json` was
never regenerated.

| location | document says | correct on the 44-target panel |
| --- | --- | --- |
| §23.4 specificity control | n = 10, separation **0.700**, p = 0.043 | n = 12, separation **0.746**, p = 0.011 |
| §23.7 reference agreement | Spearman ρ = **0.550**, n = 42 | Spearman ρ = **0.592**, n = 44 |

Group means for §23.4 on the current panel: safe 0.370 (n=19), efficacy-fail **0.219**
(n=12), toxic 0.121 (n=11). The document's table says efficacy-fail 0.244 (n=10).

**Why this went undetected:** `verify_documented_numbers.py` is designed to catch exactly
this (prose drifting from code -- the §20 precedent), but `EV-Capture-CURRENT.md` records
that the full verifier was not re-run end-to-end after the panel growth. It was run this
session and caught both immediately.

**Note that `EV-Capture-CURRENT.md` states §23.4 "was rewritten" on September 10.** The
rewrite is not in the document. Whatever happened, the numbers there are from the 42-target
panel. Worth knowing that the handoff and the artifact disagree.

**Not fixed here, deliberately.** §23.4's passage contains interpretive prose explaining
*why* the value moved (it attributes the shift to GPNMB specifically). With n going 10 -> 12
that explanation needs re-deriving, and the careful hedging in that paragraph is Dylan's
analysis to rewrite, not a mechanical find-and-replace. §23.7's is a simpler two-number fix.
Both should be done before anything is written from this document.

**What was regenerated:** `wholebody_v2.json` (now 44 labels). Every safe-vs-toxic headline
number is byte-identical after regeneration -- 0.880, 0.914, 0.751, 0.8125, leave-one-out
0.905 all unchanged, because ERBB3 and AXL are both `efficacy_fail` and no safe-vs-toxic
test touches that class. Only the Spearman moved, and it moved INTO agreement with the value
the paper skeleton already cites.

---

## FINDING 2 -- the novelty claim and the positive control were double-counting the same genes

The screen's exclusion set was widened from 44 genes (the label panel) to **111**:

| tier | n | source |
| --- | --- | --- |
| panel | 44 | `target_labels.py` |
| registry | 72 | De Schepper et al. 2025 ADC trial registry, extracted from the paper's own supplementary data already local in `ev_capture_data/adc_autopsy/` |
| curated | 35 | CAR-T / bispecific / TCR targets in development, one line of justification each |
| **union** | **111** | |

It flags **6 of the 214 "NEW" candidates**: LILRB4, TNFRSF13C, CD27, TNFRSF13B, CD79A,
CLEC12A. Corrected count: **208 genuinely new**.

**Those six are precisely the genes §30.2 cites as its positive control.** The screen was
counting the same genes twice -- once as validation ("it recovered known targets nobody told
it about") and once as discovery (inside the 214). Both framings cannot be true of the same
gene. §30.2's positive control is unaffected and remains valid; what needs correcting is the
candidate count, which should be stated as 208 with the exclusion basis named.

Implementation: `pursued_targets.py`, wired into `screen_rerank_weighted.py`, which now
flags rather than silently drops so the correction is visible.

**Honest limit:** still not exhaustive. There is no single public registry of every antigen
under development across all modalities, and the curated tier is hand-assembled. The
defensible claim is "not in the label panel, the De Schepper registry, or a curated
CAR-T/bispecific list" -- not "novel".

---

## Review of each implemented angle

### Angle 3 -- confidence intervals. Extended this session.

Originally covered three statistics. **Added the solid-tumour-only test**, which is the
smallest and most fragile (8 v 10) and therefore the one where the interval matters most.

| statistic | point | 95% CI | width |
| --- | --- | --- | --- |
| one-factor | 0.880 | [0.732, 0.986] | 0.254 |
| two-factor | 0.914 | [0.789, 0.995] | 0.206 |
| reference B | 0.751 | [0.507, 0.952] | 0.445 |
| solid-only | 0.812 | [0.575, 0.988] | 0.413 |

**A bug was caught while adding the solid-only row.** The first version reconstructed
"solid" by filtering on non-haematopoietic home, giving 0.851 at n=11v11 instead of the
documented 0.8125 at n=8v10. The real definition is an explicit hand-curated 20-gene BLOOD
set in `wholebody_v2.py`. Reconstructing a published subset by inference rather than reusing
its actual definition is precisely how a number drifts; the fix reuses the literal list and
now reproduces 0.812 exactly.

**Still uncovered:** the rapid-autopsy replication (0.750) and the trial-registry control
(0.495) live in different JSONs with a different computation and did not get intervals. Worth
doing; not done here.

### Angle 2 -- ceiling in clinical units. Unchanged, and still the weakest angle.

Machinery verified (reproduces all 20 published AUCs to three decimals via independent
code). The result remains small numbers -- BCMA 1.0% max kill at a 1% harm budget -- and the
decision recorded earlier stands: either report as an honest negative, or test the DLBCL
CyTOF cohort before concluding. **Not tested on DLBCL yet.** That remains the single most
valuable open item for this angle.

The deeper issue found and not yet resolved: applying a UNIFORM harm budget across
compartments contradicts the project's own theorem, which weights compartments. The code now
reports both uniform and non-expendable-only ceilings, but the framing question of which to
present is open.

### Angle 1 -- continuous weights. Split result, both halves locked in tests.

Clinical prediction: **worse** (0.809 vs 0.914 binary), does not clear its own scrambled-weight
null (p = 0.15). Cause identified -- continuous weighting trusts the computed home cell type,
which is an artifact for broadly-expressed targets (EGFR's computed home is placental tissue,
scored fully expendable, so a toxic target collects a large bonus). Binary is robust because
it is coarse. Binary remains primary; the negative is pinned by a test that fails loudly if
a future change reverses it.

Screen ranking: **works as intended.** DCT 135 -> 4, ERVFRD-1 136 -> 2, PSG1 138 -> 3;
AKAP10 134 -> 214, TGFBR1 131 -> 211, CD3D 7 -> 74.

**Weakness worth naming:** the gate (top_share >= 0.25) is a judgement. The threshold was
scanned across 0.00-0.40 and every value reported, but a reader could reasonably ask why 0.25
and not 0.30. The ordering is stable across the upper half of the range, which is the honest
defence, and it should be stated rather than left implicit.

### Angle 4 -- DCT dossier. Complete; see `Candidate_Dossier_DCT_Sept27.md`.

Pre-commitment recorded before investigation. Verdict: survives, with surface accessibility
as the principal unresolved risk and a narrower novelty claim (surface-directed targeting is
unexplored; DCT itself is a known antigen via vaccine/TCR routes).

---

## Infrastructure added this session

**`verify_documented_numbers.py` now covers the new work** -- previously none of it was
wired in, which broke the project's own standing discipline. 35 checks total (was 22),
covering the four confidence intervals, the label-sensitivity range and worst p, the
continuous-vs-binary comparison (as a *pinned negative* -- the check fails if continuous ever
beats binary, forcing re-examination rather than silent adoption), and the re-ranked screen's
top candidates plus the widened exclusion count.

**Test suite: 234 passing** (was 195). New files: `test_auc_confidence_intervals.py` (9),
`test_label_sensitivity.py` (7), `test_ceiling_in_clinical_units.py` (8),
`test_continuous_weights.py` (15).

---

## Open items, ranked

1. **Fix §23.4 and §23.7 in the results document.** Blocking anything written from it.
2. **Decide Angle 2's fate** -- test DLBCL CyTOF, or report the honest negative.
3. Correct the "214 candidates" figure to 208 wherever it appears.
4. Extend confidence intervals to the autopsy and trial-registry statistics.
5. Re-check whether a surface-directed DCT programme has appeared (the review consulted is
   not recent).
