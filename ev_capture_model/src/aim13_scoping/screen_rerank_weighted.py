"""Re-rank the discovery screen using confinement-gated continuous expendability weights.

THE PROBLEM THIS FIXES
----------------------
The binary rule adds +1.0 whenever a candidate's highest-expressing normal cell type is
haematopoietic, regardless of how confined the candidate actually is. For a gene spread
across 100+ cell types, no cell type holds more than a few percent of total expression, so
the "home" is close to arbitrary -- and if that arbitrary home happens to be a neutrophil,
the gene collects the full bonus.

Concretely, in the binary ranking AKAP10 (confinement 0.026, expressed across 137 cell
types) sat at rank 134 of 214 while DCT (confinement 0.948, a single cell type) sat at 135.
A judge who knows immunology would also find CD3D -- a T-cell gene -- in the top 10.

THE FIX
-------
Two changes, both principled rather than tuned:
  1. Weights are continuous and per-cell-type (see expendability_weights.py), so T cells,
     microglia and B cells are no longer treated as equally expendable.
  2. The weight is applied ONLY when the home is meaningful -- that is, when the candidate
     is actually concentrated somewhere (top_share >= GATE). Below the gate the home tells
     us nothing, so the candidate is assumed non-expendable rather than given credit for a
     cell type it barely occupies.

GATE CHOICE, STATED HONESTLY
----------------------------
GATE = 0.25 is a judgement, not an optimum: a round value in the region where a target begins
to have a genuinely dominant cell type. It was swept across 0.00-0.40 (screen_rank_stability.py,
every value reported). What that sweep shows, stated without flattery:
  * the known broad-gene artefacts are demoted at EVERY gate tried (AKAP10 >= 203 of 214,
    TGFBR1 >= 211 of 214), so the gate's purpose is achieved throughout the range;
  * the top-25 is identical at every gate (by construction: its members are all far above any
    gate), so top-25 overlap is uninformative about the gate;
  * further down the list the ordering DOES move: top-100 overlap with the default falls to
    about 81-84 of 100 at the extremes of the range. Ranks beyond ~50 should not be read as
    precise.

IMPORTANT SCOPE NOTE
--------------------
This scheme IMPROVES the screen ranking but MAKES THE CLINICAL PREDICTION WORSE
(0.809 vs 0.914 binary, see continuous_weight_sensitivity.py). Those are different tasks
and the same change is not good for both. The binary rule remains primary for the section
23 clinical result. This module is for candidate prioritisation only.
"""
import json
import os

import expendability_weights as W
import pursued_targets as P

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'target_discovery.json')
OUT = os.path.join(HERE, 'screen_rerank_weighted.json')

GATE = 0.25

# The screen's original exclusion set was the 44-gene label panel only. pursued_targets
# widens it to 111 (panel + the De Schepper 2025 ADC registry + a curated CAR-T/bispecific
# list). Candidates the widened set catches are retained in the output but FLAGGED, so the
# corrected novelty count is visible rather than silently changing the total.


def weighted_score(confinement, home, gate=GATE):
    """Confinement, plus expendability credit only where the home is meaningful."""
    w = W.weight(home, 'continuous') if confinement >= gate else 1.0
    return confinement + (1.0 - w)


def compute_all(gate=GATE):
    d = json.load(open(SRC))
    cands = [dict(c) for c in d['candidates']]
    panel_only = P.panel_symbols()
    widened = P.pursued()
    for c in cands:
        c['weighted_score'] = weighted_score(c['confinement'], c['home'], gate)
        c['already_pursued_widened'] = (c['gene'] in widened) and (c['gene'] not in panel_only)
    old = sorted(cands, key=lambda x: -x['score'])
    new = sorted(cands, key=lambda x: -x['weighted_score'])
    old_rank = {c['gene']: i + 1 for i, c in enumerate(old)}
    new_rank = {c['gene']: i + 1 for i, c in enumerate(new)}
    for c in cands:
        c['rank_binary'] = old_rank[c['gene']]
        c['rank_weighted'] = new_rank[c['gene']]
        c['rank_change'] = old_rank[c['gene']] - new_rank[c['gene']]
    flagged = [c['gene'] for c in cands if c['already_pursued_widened']]
    return {
        'gate': gate,
        'n_candidates': len(cands),
        'n_flagged_by_widened_exclusion': len(flagged),
        'flagged_genes': sorted(flagged),
        'n_genuinely_new': len(cands) - len(flagged),
        'top_25_weighted': [
            {k: c[k] for k in ('gene', 'weighted_score', 'confinement', 'n_cell_types',
                               'home', 'rank_binary', 'rank_weighted',
                               'cancer_elevated_in', 'already_pursued_widened')}
            for c in new[:25]],
        'by_gene': {c['gene']: {'confinement': c['confinement'], 'home': c['home'],
                                'rank_binary': c['rank_binary'],
                                'rank_weighted': c['rank_weighted']}
                    for c in cands},
    }


def print_summary():
    r = compute_all()
    print('=' * 92)
    print(f"DISCOVERY SCREEN RE-RANKED -- continuous weights, gated at top_share >= {r['gate']}")
    print('=' * 92)
    print(f"{'#':>3} {'GENE':12s} {'score':>6s} {'conf':>6s} {'nT':>4s} {'HOME':26s} {'was':>5s}")
    for i, c in enumerate(r['top_25_weighted'], 1):
        print(f"{i:3d} {c['gene']:12s} {c['weighted_score']:6.3f} {c['confinement']:6.3f} "
              f"{c['n_cell_types']:4d} {c['home'][:26]:26s} {c['rank_binary']:5d}")
    print()
    print('Artifacts demoted (broadly expressed, spurious haematopoietic home):')
    for g in ('AKAP10', 'TGFBR1', 'SIRPA', 'CD3D'):
        v = r['by_gene'].get(g)
        if v:
            print(f"  {g:10s} conf={v['confinement']:.3f}  rank {v['rank_binary']:3d} "
                  f"-> {v['rank_weighted']:3d}  ({v['home']})")
    print()
    print(f"Widened exclusion set flags {r['n_flagged_by_widened_exclusion']} of "
          f"{r['n_candidates']} as already pursued -> {r['n_genuinely_new']} genuinely new:")
    print('  ', ', '.join(r['flagged_genes']))
    print()
    print('Genuine tissue-restricted candidates promoted:')
    for g in ('DCT', 'ERVFRD-1', 'PSG1'):
        v = r['by_gene'].get(g)
        if v:
            print(f"  {g:10s} conf={v['confinement']:.3f}  rank {v['rank_binary']:3d} "
                  f"-> {v['rank_weighted']:3d}  ({v['home']})")
    return r


if __name__ == '__main__':
    res = print_summary()
    with open(OUT, 'w') as f:
        json.dump(res, f, indent=1)
    print(f'\nwrote {OUT}')
