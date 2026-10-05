"""How much of the discovery-screen shortlist is an artefact of hand-assigned weights?

The continuous expendability weights (expendability_weights.py) are judgements. A shortlist that
only exists at one exact set of weights is not a finding, so this asks two questions and reports
the answers whichever way they fall:

1. WEIGHT UNCERTAINTY. Every cell-type weight is perturbed by independent uniform noise of
   +/- DELTA (clipped to [0, 1]) and the screen is re-ranked, N_DRAWS times. For each candidate:
   its median rank, its 5th-95th percentile rank, and the fraction of draws in which it stays in
   the top 10 / top 25. A candidate that is top-25 in >= 90% of draws survives the judgement
   uncertainty; one that appears in 40% does not.

2. GATE SENSITIVITY. The confinement gate (default 0.25) is also a judgement. It is swept across
   0.00-0.40 and the top-25 list is compared with the default by overlap. The scheme's docstring
   claims the ordering is stable across the upper half of that range; this measures whether that
   claim is true, rather than assuming it.

Candidates already pursued (widened exclusion set) are kept in the ranking but reported
separately, so the 'stable AND genuinely new' shortlist is explicit.
"""
import json
import os
import numpy as np

import expendability_weights as W
import pursued_targets as P
import candidate_hygiene as HY

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'target_discovery.json')
OUT = os.path.join(HERE, 'screen_rank_stability.json')

GATE = 0.25
DELTA = 0.25
N_DRAWS = 2000
SEED = 20260928
GATES = (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)


def _load(hygiene=True):
    """Candidate arrays. With hygiene=True (default) the hard exclusions of candidate_hygiene.py
    (immunoglobulin/TCR genes; genes whose only surface claim is contradicted by HPA's topology
    prediction) are removed BEFORE ranking, so ranks are ranks among plausible surface targets."""
    cands = json.load(open(SRC))['candidates']
    if hygiene:
        ann = HY.annotate([c['gene'] for c in cands])
        cands = [c for c in cands if not ann[c['gene']]['hard_exclude']]
    genes = [c['gene'] for c in cands]
    conf = np.array([c['confinement'] for c in cands], float)
    homes = [c['home'] for c in cands]
    return genes, conf, homes


def _ranks(scores):
    order = np.argsort(-scores, kind='stable')
    r = np.empty(len(scores), int)
    r[order] = np.arange(1, len(scores) + 1)
    return r


def _base_weights(homes):
    return {h: W.weight(h, 'continuous') for h in set(homes)}


def _scores(conf, homes, wmap, gate):
    w = np.array([wmap[h] for h in homes], float)
    return conf + np.where(conf >= gate, 1.0 - w, 0.0)


def weight_perturbation(n_draws=N_DRAWS, delta=DELTA, seed=SEED, gate=GATE, hygiene=True):
    genes, conf, homes = _load(hygiene)
    base = _base_weights(homes)
    rng = np.random.default_rng(seed)
    types = sorted(base)
    ranks = np.empty((n_draws, len(genes)), int)
    for i in range(n_draws):
        wmap = {t: float(np.clip(base[t] + rng.uniform(-delta, delta), 0.0, 1.0)) for t in types}
        ranks[i] = _ranks(_scores(conf, homes, wmap, gate))
    base_rank = _ranks(_scores(conf, homes, base, gate))
    return genes, base_rank, ranks


ARTEFACTS = ('AKAP10', 'TGFBR1', 'SIRPA', 'CD3D')   # broad genes the binary rule wrongly promoted


def gate_sweep():
    """The gate exists to demote broadly-expressed genes whose 'home' cell type is arbitrary. The
    top-25 is insensitive to it BY CONSTRUCTION (top-25 members all have confinement well above
    any gate tried), so top-25 overlap says nothing. What the gate controls is where the known
    artefacts land, and how far down the list the ordering changes. Both are reported, using the
    UNFILTERED candidate list because the artefacts are hygiene-passing broad genes."""
    genes, conf, homes = _load(hygiene=False)
    base = _base_weights(homes)
    ref_rank = _ranks(_scores(conf, homes, base, GATE))
    out = {}
    for g in GATES:
        rk = _ranks(_scores(conf, homes, base, g))
        row = {'artefact_ranks': {a: int(rk[genes.index(a)]) for a in ARTEFACTS if a in genes}}
        for k in (25, 50, 100):
            t = set(np.array(genes)[np.argsort(rk)[:k]])
            r = set(np.array(genes)[np.argsort(ref_rank)[:k]])
            row[f'top{k}_overlap'] = len(t & r)
        out[str(g)] = row
    return out


def compute_all():
    genes, base_rank, ranks = weight_perturbation()
    widened = P.pursued()
    panel = P.panel_symbols()
    order = np.argsort(base_rank)
    rows = []
    for j in order[:40]:
        r = ranks[:, j]
        rows.append({
            'gene': genes[j], 'baseline_rank': int(base_rank[j]),
            'median_rank': float(np.median(r)),
            'p5_rank': float(np.percentile(r, 5)), 'p95_rank': float(np.percentile(r, 95)),
            'frac_top10': float(np.mean(r <= 10)), 'frac_top25': float(np.mean(r <= 25)),
            'already_pursued': bool((genes[j] in widened) and (genes[j] not in panel)),
        })
    stable = [r for r in rows if r['frac_top25'] >= 0.90]
    stable_new = [r['gene'] for r in stable if not r['already_pursued']]
    return {
        'method': {'delta': DELTA, 'n_draws': N_DRAWS, 'seed': SEED, 'gate': GATE},
        'top40_baseline': rows,
        'n_top25_in_>=90pct_of_draws': len(stable),
        'stable_and_not_already_pursued': stable_new,
        'n_ranked_after_hygiene': len(genes),
        'gate_sweep': gate_sweep(),
    }


def print_summary(r):
    print('=' * 92)
    print(f"SCREEN RANK STABILITY  (weights +/-{r['method']['delta']}, {r['method']['n_draws']} draws)")
    print('=' * 92)
    print(f"{'base':>4} {'GENE':12s} {'median':>7s} {'5-95%':>11s} {'top10':>6s} {'top25':>6s}  pursued")
    for x in r['top40_baseline'][:30]:
        print(f"{x['baseline_rank']:4d} {x['gene']:12s} {x['median_rank']:7.1f} "
              f"{x['p5_rank']:5.0f}-{x['p95_rank']:<5.0f} {x['frac_top10']:6.2f} "
              f"{x['frac_top25']:6.2f}  {'YES' if x['already_pursued'] else ''}")
    print(f"\nbaseline top-25 members in the top 25 in >= 90% of draws: "
          f"{r['n_top25_in_>=90pct_of_draws']} (of the top 40 examined)")
    print('stable AND not already pursued:', ', '.join(r['stable_and_not_already_pursued']))
    print('\nGATE SWEEP (unfiltered list of', 214, '): overlap with the default gate, and where the')
    print('known broad-gene artefacts land (rank of 214; higher = more demoted)')
    for g, v in r['gate_sweep'].items():
        print(f"  gate {g:5s} top25 {v['top25_overlap']:2d}/25  top50 {v['top50_overlap']:2d}/50  "
              f"top100 {v['top100_overlap']:3d}/100   artefacts {v['artefact_ranks']}")


if __name__ == '__main__':
    res = compute_all()
    print_summary(res)
    json.dump(res, open(OUT, 'w'), indent=1)
    print('\nwrote', OUT)
