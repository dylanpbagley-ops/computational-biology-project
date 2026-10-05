"""Candidate hygiene: are the screen's candidates things an antibody or CAR could actually bind?

WHY THIS EXISTS
---------------
The discovery screen's surface filter accepts a gene if HPA lists it as a 'CD marker' OR reports
'plasma membrane' among its observed immunofluorescence (IF) locations. The second criterion is
weak: it fires on ANY plasma-membrane staining pattern in ANY cell line, including antibody
cross-reactivity. Inspecting the top of the re-ranked list showed the consequence -- immunoglobulin
kappa variable genes (clonal B-cell-receptor fragments, not a target class), GAGE-family genes
(multi-copy, intracellular cancer-testis antigens), and secreted proteins (PSG1, SPINK13, LIPI).

RULES
-----
HARD EXCLUSIONS (unambiguous cases only):
  1. ig_or_tcr        protein class includes 'Immunoglobulin genes' or 'T-cell receptor genes'
                      (clonotype-specific variable regions; multi-mapping in RNA-seq).
  2. no_membrane_support  the class carries NEITHER 'Predicted membrane proteins' NOR 'CD markers'
                      NOR a receptor/transporter class, yet DOES carry 'Predicted intracellular
                      proteins' or 'Predicted secreted proteins'. The gene's only surface claim is
                      then a stray immunofluorescence pattern that HPA's own topology prediction
                      contradicts (GAGE-family genes, LIPI, SPINK13).

SOFT FLAGS (reported, never used to exclude):
  secretome_annotated  HPA lists a Secretome location. Shed and cleaved receptors (MET, FAP,
                       PTK7) carry it legitimately, so it cannot be an exclusion criterion.
  if_uncertain         IF reliability is 'Uncertain'.
  tier                 A = 'CD markers' or 'Predicted membrane proteins' in the class;
                       B = surface status rests on the IF plasma-membrane call alone.

REVISION HISTORY (kept because how a filter was arrived at is part of its validity)
----------------------------------------------------------------------------------
v1 also hard-excluded any 'Predicted intracellular' gene without a CD-marker class, and any
gene with a Secretome location but no membrane class. It was checked against the 111 known
pursued targets and falsely excluded 27 of 107 annotated ones (MET, MSLN, DLL3, FOLH1, ROR1,
PTK7 ...): HPA routinely gives real membrane proteins BOTH labels. v1 was abandoned and v2
(above) was built to keep positive membrane evidence decisive. The false-exclusion rate on
known targets is re-measured in false_exclusion_check() and pinned in the tests.

This is quality control on the pipeline, not a new result. Its own validity check is that known
targets must survive it: see positive_control_survival().
"""
import csv
import io
import json
import os
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))
OUT = os.path.join(HERE, 'candidate_hygiene.json')


def load_hpa(genes):
    genes = set(genes)
    z = zipfile.ZipFile(os.path.join(D, 'proteinatlas.tsv.zip'))
    out = {}
    with z.open(z.namelist()[0]) as f:
        rd = csv.DictReader(io.TextIOWrapper(f, 'utf-8', newline=''), delimiter='\t')
        for r in rd:
            g = r['Gene']
            if g in genes:
                out[g] = {'class': (r.get('Protein class') or '').lower(),
                          'secretome': (r.get('Secretome location') or '').strip(),
                          'main': (r.get('Subcellular main location') or '').lower(),
                          'rel_if': (r.get('Reliability (IF)') or '').strip()}
    return out


def classify(rec):
    """Return dict of flags for one HPA record (None -> unannotated)."""
    if rec is None:
        return {'annotated': False, 'hard_exclude': False, 'reasons': [], 'tier': None,
                'if_uncertain': False}
    c = rec['class']
    cd = 'cd markers' in c
    memb = 'predicted membrane proteins' in c
    recept = any(k in c for k in ('g-protein coupled receptors', 'transporters', 'voltage-gated',
                                  'ion channels'))
    reasons = []
    if 'immunoglobulin genes' in c or 't-cell receptor genes' in c:
        reasons.append('ig_or_tcr')
    if (not (cd or memb or recept)) and ('predicted intracellular proteins' in c or
                                         'predicted secreted proteins' in c):
        reasons.append('no_membrane_support')
    return {'annotated': True, 'hard_exclude': bool(reasons), 'reasons': reasons,
            'tier': 'A' if (cd or memb) else 'B', 'if_uncertain': rec['rel_if'] == 'Uncertain',
            'secretome_annotated': bool(rec['secretome'])}


def annotate(genes):
    hpa = load_hpa(genes)
    return {g: classify(hpa.get(g)) for g in genes}


def positive_control_survival():
    """Every known target the screen recovered must pass hygiene, or the filter is destroying
    real signal. Reads the recovery list from the existing positive-control output."""
    pc = json.load(open(os.path.join(HERE, 'target_discovery_positive_control.json')))
    recovered = [g for g, r in pc['ranks'].items() if r is not None]
    ann = annotate(recovered)
    return {'recovered': recovered,
            'fail_hygiene': {g: ann[g]['reasons'] for g in recovered if ann[g]['hard_exclude']},
            'tier_B': [g for g in recovered if ann[g]['tier'] == 'B']}


def false_exclusion_check():
    """Apply the filter to clinically pursued surface antigens (the 44-gene label panel plus the
    ADC registry and curated CAR-T/bispecific lists). Any hard exclusion here is a candidate false
    exclusion -- an UPPER bound, because a few pursued targets are not membrane proteins at all
    (fibronectin FN1 and the secreted kallikrein KLK2 are pursued, but are not cell-surface)."""
    import pursued_targets as P
    known = sorted(P.pursued())
    ann = annotate(known)
    annotated = [g for g in known if ann[g]['annotated']]
    fails = {g: ann[g]['reasons'] for g in annotated if ann[g]['hard_exclude']}
    return {'n_known': len(known), 'n_annotated': len(annotated), 'n_fail': len(fails),
            'false_exclusion_rate_upper_bound': len(fails) / max(1, len(annotated)),
            'fails': fails}


def compute_all():
    cands = json.load(open(os.path.join(HERE, 'target_discovery.json')))['candidates']
    genes = [c['gene'] for c in cands]
    ann = annotate(genes)
    n_excl = sum(a['hard_exclude'] for a in ann.values())
    by_reason = {}
    for a in ann.values():
        for r in a['reasons']:
            by_reason[r] = by_reason.get(r, 0) + 1
    import pursued_targets as P
    widened, panel = P.pursued(), P.panel_symbols()
    pursued_now = {g for g in genes if g in widened and g not in panel}
    passing = [g for g in genes if not ann[g]['hard_exclude']]
    return {'n_candidates': len(genes), 'n_hard_excluded': n_excl,
            'n_pass_and_not_already_pursued': len([g for g in passing if g not in pursued_now]),
            'n_pass': len(genes) - n_excl, 'by_reason': by_reason,
            'tier_counts_among_pass': {t: sum(1 for a in ann.values()
                                              if not a['hard_exclude'] and a['tier'] == t)
                                       for t in ('A', 'B')},
            'positive_control': positive_control_survival(),
            'false_exclusion_check': false_exclusion_check(),
            'annotations': ann}


if __name__ == '__main__':
    r = compute_all()
    print(f"candidates {r['n_candidates']}: hard-excluded {r['n_hard_excluded']} "
          f"({r['by_reason']}), pass {r['n_pass']}; tiers among pass {r['tier_counts_among_pass']}")
    pcs = r['positive_control']
    print('positive control: recovered', len(pcs['recovered']), '| fail hygiene:', pcs['fail_hygiene'],
          '| tier B:', pcs['tier_B'])
    fx = r['false_exclusion_check']
    print(f"false-exclusion check on known targets: {fx['n_fail']}/{fx['n_annotated']} "
          f"({100*fx['false_exclusion_rate_upper_bound']:.0f}%)  fails: {fx['fails']}")
    json.dump(r, open(OUT, 'w'), indent=1)
    print('wrote', OUT)
