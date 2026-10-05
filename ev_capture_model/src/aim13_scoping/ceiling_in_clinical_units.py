"""The ceiling, expressed in units a clinician recognises.

WHY THIS EXISTS
---------------
Every ceiling in this project is currently reported as an AUC -- rank-based, dimensionless,
and abstract. An AUC of 0.79 does not tell a reader what the target actually permits.

The same ROC curves already computed contain a far more legible quantity: for a fixed harm
budget, the maximum fraction of tumour cells any therapy against that target could kill.
That is read directly off the curve -- no new data, no new assumptions, no new model. It is
the same theorem, reported at a point rather than integrated.

    "At most X% of malignant plasma cells can be killed while sparing 99% of <compartment>."

WHY THIS IS NOT BLOCKED BY THE CALIBRATION GAP
----------------------------------------------
Sections 11 and 11a establish that CITE-seq and CyTOF signal is NOT absolutely calibrated to
molecules per cell, so a real drug's operating point cannot be placed on these curves. That
limitation does not apply here. Every quantity below is a FRACTION OF A POPULATION -- what
share of tumour cells sit above a threshold versus what share of normal cells do. Fractions
are invariant to any monotone rescaling of the antigen axis, so an uncalibrated axis changes
where the threshold sits without changing the achievable trade-off. This must be stated
explicitly in the write-up so the result is not mistaken for a calibrated operating point.

METHOD
------
Population gating mirrors `myeloma_ceiling.py` exactly (CLR normalisation, same percentile
gates, same subsampling seed) so the numbers are directly comparable to the published
ceilings. For each protein and each normal compartment the full empirical ROC is built, and
the maximum tumour-kill fraction is read off at each harm budget. A target must spare ALL
normal compartments, so the reported ceiling is the MINIMUM across compartments -- the worst
compartment is the binding constraint, exactly as in the AUC version.

Self-check: the area under the ROC curves built here must reproduce the AUCs in
`myeloma_ceiling_results.json`. The test file asserts this. If the curve were wrong, the
area under it would not match a number computed by entirely separate code.
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ADT = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'myeloma_adt'))
OUT = os.path.join(HERE, 'ceiling_in_clinical_units.json')

HARM_BUDGETS = (0.01, 0.05)
SUB = 8000
SEED = 53

COMPARTMENT_NAMES = {'T': 'T cells', 'MY': 'myeloid cells', 'NK': 'NK cells',
                     'NB': 'normal B cells'}

# The tumour gate is defined as CD138-high OR CD269(BCMA)-high. Evaluating either marker
# as a TARGET against a population it helped define is circular, and measurably so: the
# circular version inflates CD138 from 1.6% to 2.7% and CD269 from 0.8% to 1.6% at a 1%
# harm budget. For these two markers the gate is rebuilt using only the OTHER one.
GATE_MARKERS = ('CD138', 'CD269')

# Compartments the project argues are expendable. Applying a uniform harm budget across
# all compartments contradicts the central claim -- for a target whose binding constraint
# is an expendable lineage, that constraint is precisely the harm clinical practice
# accepts. Both variants are reported.
EXPENDABLE_COMPARTMENTS = frozenset({'NB'})

APPROVED = {'CD269': 'BCMA -- ide-cel / cilta-cel / teclistamab',
            'CD38': 'daratumumab',
            'CD319': 'SLAMF7 -- elotuzumab'}
TRIALS = {'CD307': 'FcRH5 -- cevostamab'}
FAILED = {'CD138': 'indatuximab ravtansine -- discontinued'}


def clr(M):
    X = M + 1.0
    return np.log(X) - np.log(X).mean(axis=1, keepdims=True)


def max_kill_at_harm(tumor, normal, budget):
    """Maximum fraction of `tumor` killed while killing at most `budget` of `normal`.

    Thresholds the antigen value and keeps the best achievable tumour kill subject to the
    harm constraint. Equivalent to reading TPR off the empirical ROC at FPR = budget.
    """
    tumor = np.asarray(tumor, float)
    normal = np.asarray(normal, float)
    # The harm constraint binds at the (1-budget) quantile of the normal population: killing
    # everything above it costs exactly `budget` of that compartment.
    thresh = np.quantile(normal, 1.0 - budget)
    return float((tumor > thresh).mean())


def roc_auc_from_curve(tumor, normal, n_grid=2001):
    """AUC computed by integrating the empirical ROC built the same way max_kill_at_harm
    reads it. Used as a self-check against the independently-computed published AUCs."""
    tumor = np.asarray(tumor, float)
    normal = np.asarray(normal, float)
    fpr_grid = np.linspace(0.0, 1.0, n_grid)
    tpr = np.array([max_kill_at_harm(tumor, normal, max(f, 1e-12)) for f in fpr_grid])
    return float(np.trapezoid(tpr, fpr_grid)) if hasattr(np, 'trapezoid') else float(np.trapz(tpr, fpr_grid))


def load_populations(exclude_gate_marker=None):
    """Gate each sample into tumour and normal compartments. Mirrors myeloma_ceiling.py.

    `exclude_gate_marker` drops one of the two gate-defining markers from the tumour
    definition, so that marker can be evaluated as a target without circularity.
    """
    rng = np.random.default_rng(SEED)
    sub = lambda a: a[rng.choice(len(a), SUB, replace=False)] if len(a) > SUB else a
    files = sorted(glob.glob(os.path.join(ADT, '*.npz')))
    if not files:
        raise FileNotFoundError(f'no myeloma ADT .npz files under {ADT}')
    names = list(np.load(files[0], allow_pickle=True)['names'])
    ix = {n: i for i, n in enumerate(names)}
    samples = []
    for f in files:
        sid = os.path.basename(f)[:-4]
        M = clr(np.load(f, allow_pickle=True)['adt'].astype(float))
        cd3, cd14, cd33, cd56, cd19, cd138, bcma = (
            M[:, ix[k]] for k in ['CD3E', 'CD14', 'CD33', 'CD56', 'CD19', 'CD138', 'CD269'])
        hi = lambda v, p=75: v > np.percentile(v, p)
        T = hi(cd3, 70)
        MY = (~T) & (hi(cd14, 80) | hi(cd33, 80))
        NK = (~T) & (~MY) & hi(cd56, 85)
        NB = (~T) & (~MY) & (~NK) & hi(cd19, 85) & (cd138 <= np.percentile(cd138, 70))
        crit = []
        if exclude_gate_marker != 'CD138':
            crit.append(hi(cd138, 70))
        if exclude_gate_marker != 'CD269':
            crit.append(hi(bcma, 80))
        tumor_crit = crit[0]
        for extra in crit[1:]:
            tumor_crit = tumor_crit | extra
        TU = (~T) & (~MY) & (~NK) & (~NB) & tumor_crit
        if TU.sum() < 100:
            continue
        comps = {}
        for lab, mask in (('T', T), ('MY', MY), ('NK', NK), ('NB', NB)):
            if mask.sum() >= 100:
                comps[lab] = sub(M[mask])
        samples.append({'sid': sid, 'tumor': sub(M[TU]), 'comps': comps})
    return names, ix, samples


def compute_all():
    names, ix, samples = load_populations()
    # separate, non-circular gates for the two markers that define the default gate
    alt = {m: load_populations(exclude_gate_marker=m)[2] for m in GATE_MARKERS}
    out = {}
    for n in names:
        use = alt[n] if n in GATE_MARKERS else samples
        per_budget = {}
        for budget in HARM_BUDGETS:
            per_comp = {}
            for lab in ('T', 'MY', 'NK', 'NB'):
                vals = [max_kill_at_harm(s['tumor'][:, ix[n]], s['comps'][lab][:, ix[n]], budget)
                        for s in use if lab in s['comps']]
                if vals:
                    per_comp[lab] = float(np.median(vals))
            if per_comp:
                worst = min(per_comp, key=per_comp.get)
                nonexp = {k: v for k, v in per_comp.items()
                          if k not in EXPENDABLE_COMPARTMENTS}
                worst_ne = min(nonexp, key=nonexp.get) if nonexp else None
                per_budget[f'{budget:.2f}'] = {
                    'per_compartment': per_comp,
                    'ceiling': per_comp[worst],
                    'binding_compartment': worst,
                    'binding_compartment_name': COMPARTMENT_NAMES[worst],
                    'ceiling_nonexpendable_only': nonexp[worst_ne] if worst_ne else None,
                    'binding_compartment_nonexpendable': worst_ne,
                    'gate_decircularised': n in GATE_MARKERS,
                }
        if per_budget:
            out[n] = per_budget
    out['_meta'] = {
        'harm_budgets': list(HARM_BUDGETS),
        'n_samples_used': len(samples),
        'gate_markers_decircularised': list(GATE_MARKERS),
        'expendable_compartments': sorted(EXPENDABLE_COMPARTMENTS),
        'units_note': ('All values are FRACTIONS OF A POPULATION and are invariant to any '
                       'monotone rescaling of the antigen axis. They are therefore NOT '
                       'affected by the absolute-calibration gap in sections 11/11a, and '
                       'must not be presented as calibrated operating points.'),
    }
    return out


def print_summary():
    r = compute_all()
    print('=' * 86)
    print('THE CEILING IN CLINICAL UNITS -- multiple myeloma CITE-seq')
    print('Max %% of malignant plasma cells killable while sparing 99%% / 95%% of every normal')
    print('compartment. Fractions, not molecules -- calibration-invariant.')
    print('=' * 86)
    for budget in HARM_BUDGETS:
        key = f'{budget:.2f}'
        print(f'\n--- harm budget: at most {budget*100:.0f}% of any normal compartment killed ---')
        print(f"{'PROTEIN':10s} {'MAX KILL':>9s}  {'LIMITED BY':16s} note")
        rows = [(n, v[key]) for n, v in r.items()
                if not n.startswith('_') and key in v]
        for n, v in sorted(rows, key=lambda x: -x[1]['ceiling']):
            note = ''
            if n in APPROVED:
                note = '<-- APPROVED: ' + APPROVED[n]
            elif n in TRIALS:
                note = '<-- in trials: ' + TRIALS[n]
            elif n in FAILED:
                note = '<-- FAILED: ' + FAILED[n]
            if not note and v['ceiling'] < 0.20:
                continue
            print(f"{n:10s} {v['ceiling']*100:8.1f}%  {v['binding_compartment_name']:16s} {note}")
    print()
    print('NOTE:', r['_meta']['units_note'])
    return r


if __name__ == '__main__':
    res = print_summary()
    with open(OUT, 'w') as f:
        json.dump(res, f, indent=1)
    print(f'\nwrote {OUT}')
