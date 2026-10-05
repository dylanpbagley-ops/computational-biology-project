"""Apply dated source/document corrections against preserved inputs."""
import hashlib,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parent;BASE=ROOT.parent
manifest=json.loads((ROOT/'before_hashes.json').read_text())
changes=[]
def edit(rel,transform):
    p=BASE/rel
    assert hashlib.sha256(p.read_bytes()).hexdigest()==manifest[str(pathlib.Path(rel))], 'File changed since snapshot: '+rel
    old=p.read_text(encoding='utf-8');new=transform(old);assert old!=new
    p.write_text(new,encoding='utf-8');changes.append({'file':rel,'before':manifest[str(pathlib.Path(rel))],'after':hashlib.sha256(p.read_bytes()).hexdigest()})
def theorem(text):
    end=text.index('"""',3)+3
    header='''"""
WEIGHTED MULTI-COMPARTMENT RESPONSE BOUNDS — scope corrected October 2, 2026.

For fixed nonnegative weights, specified compartment prevalences, and a shared
response f, H/W equals E_pW[f] when W>0. This is an exact mixture identity,
not a clinical harm calibration. W=0 is a separate zero-modeled-harm case.

For nondecreasing antigen responses, attainable pairs form the convex hull of
antigen-threshold ROC points. For arbitrary responses, they form the convex hull
of ALL deterministic subsets, with a likelihood-ratio upper boundary. The latter
can exceed the former. AUC does not determine either entire boundary.

The numerical code below is preserved: it checks arbitrary-response upper-bound
containment and linear-utility optimality, NOT monotone-response attainability.
Legacy special-case keys check a single-reference identity, one-hot reference
identities, and share-weighted algebra. They do not establish joint tissue-budget
attainability or derive the empirical confinement-plus-blood score.

See outputs/strengthening_Oct2/mathematical_scope.md and its independently tested
frontiers.py implementation for corrected statements and counterexamples.
The pre-correction source and output are preserved in strengthening_Oct2/before.
"""'''
    text=header+text[end:]
    start=text.index('    """For each random problem, check three things:');end=text.index('    worst_excess',start)
    text=text[:start]+'''    """Check arbitrary-response upper-bound containment and LR utility optimality.

    The loop performs 40 containment checks and one utility check per problem by
    default. It does not check the monotone class or a third single-case identity;
    legacy algebraic identities are evaluated separately in special_cases().
    """
'''+text[end:]
    text=text.replace('"""The three results the theorem is supposed to contain, checked numerically."""','"""Legacy algebraic checks only; see corrected scope note for limitations."""')
    text=text.replace('# 2. worst-compartment: min over one-hot weightings equals the min over compartments of AUC','# 2. One-hot reference identity only; NOT a simultaneous-constraint optimum.')
    text=text.replace("print('\\n  the three special cases the theorem should contain:')","print('\\n  legacy algebraic identities (not clinical or joint-attainability tests):')")
    return text
edit('ev_capture_model/src/aim13_scoping/weighted_ceiling_theorem.py',theorem)
banner='''> **October 2 correction notice:** This is a historical analysis record. For current scientific scope and unresolved issues use `strengthening_Oct2/Strengthening_Report.md` and `strengthening_Oct2/mathematical_scope.md` in the outputs root. Original scores and labels are preserved, not newly adjudicated. The autopsy comparison reuses the project outcome labels; it is not independent clinical replication. Nominal discovery-panel p-values do not account for label, selection, or full research-history uncertainty. Pre-correction copies are retained in `strengthening_Oct2/before/`.

'''
def results(text):
    start=text.index('### 29.1 Setup');end=text.index('### 29.4 Numerical verification',start)
    new='''### 29.1–29.3 Corrected scope (October 2, 2026)

These mathematical statements supersede the earlier formulation; the original is archived. With fixed nonnegative weights and a shared response, H/W = E_pW[f] for W>0. For nondecreasing antigen responses, the feasible region is the convex hull of antigen-threshold ROC points. Arbitrary responses have the larger all-subsets feasible region, whose upper boundary is ordered by likelihood ratio. AUC alone determines neither full frontier. W=0 needs separate treatment. The harm functional is a specified model, not a calibrated measure of clinical injury.

One-hot weighting recovers each compartment's reference, but the minimum of individually optimized ceilings need not be attainable under simultaneous compartment constraints. The expression-share identity does not derive the empirical two-factor score. The corrected proof, assumptions, counterexamples, and verification are in `strengthening_Oct2/mathematical_scope.md` and `frontier_results.json` in the outputs root.

'''
    text=text[:start]+new+text[end:]
    text=text.replace('**It is not. They are the same statement.**','**They are different: a response bound under explicit assumptions and an empirical ranking hypothesis. The bound motivates investigation but does not derive the score.**')
    text=text.replace('**Section 22\'s finding replicates in a solid tumour, in an independent cohort, from an\nindependent laboratory, on a different technology, with a matched-normal design.**','**This supplies independent expression measurements for a comparator, but reuses the original outcome labels. It is not independent clinical replication. The paired AUC difference has an interval spanning zero (October 1 audit).**')
    return banner+text
edit('ev_capture_model/docs/Aim13_Results_Two_Comparator_Ceiling.md',results)
edit('HANDOFF_START_HERE_Oct1.md',lambda text:'''# Current entry point: PROJECT_STATUS_Oct2.md

October 2 update: Dylan authorized strengthening and relevant expansion of the existing biology project. The older restrictions against any expansion and strict-label analysis below are superseded. Frozen historical artifacts and locked predictions remain preserved. Read `PROJECT_STATUS_Oct2.md` first; this October 1 handoff is retained as history and contains superseded scientific claims and test counts.

'''+text)
edit('ev_capture_model/docs/Label_Review_Sheet.md',lambda text:'''# Superseded review sheet — use the October 2 score-free form

This historical sheet exposes the original labels to reviewers and contains known source errors. Use `strengthening_Oct2/reviewer_form_score_free.tsv` together with `Reviewer_Instructions.md` in that directory. Do not send the review key or model scores to the reviewer. The historical content below is preserved for provenance, not current clinical adjudication.

'''+text)
edit('Research_Direction_Reset_Oct2.md',lambda text:'''# Superseded direction decision

Dylan subsequently chose to retain and strengthen the existing biology project, including relevant expansion. The external-field pivot described below is not the active plan. Read `PROJECT_STATUS_Oct2.md`.

'''+text)
(ROOT/'document_changes.json').write_text(json.dumps(changes,indent=2));print(json.dumps(changes,indent=2))
