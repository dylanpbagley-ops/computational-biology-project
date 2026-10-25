"""
MAGNITUDE-BASELINE ABLATION -- is confinement (top_share) actually doing the work,
or would any expression-magnitude feature separate safe from toxic just as well?

WHY THIS EXISTS
---------------
Published systematic attempts to predict clinical tolerability from normal-tissue
expression (Thomas et al. 2025; Madan et al. 2023) reported NO signal. This project's
claim is that the signal lives in the CONCENTRATION of expression (top_share: the share
of total expression concentrated in the single most-expressing cell type), not in its
MAGNITUDE (mean or peak level). Until now that claim was asserted, not demonstrated.
This script demonstrates or refutes it head-to-head: same targets, same labels, same
separation definition, same bootstrap -- only the feature changes.

DESIGN (pre-stated before running)
----------------------------------
Features, per gene, all derived from the SAME HPA single-cell consensus
(rna_single_cell_type.tsv.zip), parsed exactly as wholebody_v2.load() parses it:
  top_share  max(nCPM) / sum(nCPM) over the cell-type rows present   (concentration)
  mean       mean(nCPM) over those rows (monotone with total)        (magnitude)
  max        max(nCPM) over those rows (= top_share's numerator)     (peak magnitude)

Separation convention (identical to wholebody_v2 / auc_confidence_intervals):
  sep(safe_values, toxic_values) = P(safe > toxic) (+0.5 for ties), so HIGHER sep
  means the feature ranks SAFE targets higher. If magnitude were the predictive
  feature, mean/max would land BELOW 0.5. The headline number is the paired bootstrap
  delta sep(top_share-based) - sep(magnitude-based) on the SAME resampled targets.

Feature sets:
  one-factor   top_share / mean / max                       (primary comparison; sep is
                                                            rank-based, so no scaling)
  two-factor   top_share + haem                             (the PUBLISHED rule; anchor)
               rank01(mean) + haem, rank01(max) + haem      (scale-matched analogues: the
                                                            published rule adds a 0/1 flag
                                                            to a [0,1] feature, which is
                                                            meaningless for raw nCPM scale)
  exploratory  rank01(mean) + top_share + haem              (does concentration add anything
                                                            beyond magnitude? diagnostic
                                                            only -- NOT a scoring-rule
                                                            proposal; the published rule
                                                            stays frozen)

Post-hoc rows (added 2026-10-25 after the first run; labeled in output and JSON):
  flipped_mean / flipped_max          magnitude used in its BEST direction -- the field's own
                                      heuristic (low normal expression = safer target).
                                      One-factor flips are exact (sep is rank-based, so
                                      negating a feature reverses its orientation); two-factor
                                      flips are scale-matched as (1 - rank01) + haem so the
                                      flag addition stays meaningful.
  Purpose: the "did you try the obvious baseline?" defense -- is concentration still
  ahead when magnitude gets its best shot? Interpretation branches, pre-stated:
    concentration still ahead (delta CI excludes 0)  -> claim survives its strongest challenge
    flipped magnitude ties concentration             -> claim shifts to "concentration captures
                                                        directly what magnitude needs help to see"
    flipped magnitude ahead                          -> feature claim fails; report it

Anchors (hard checks): using the COMMITTED top_share values from wholebody_v2.json, the
one-factor and two-factor toxicity separations must reproduce the published 0.8804 and
0.9139, and the specificity control must reproduce 0.7456. The script hard-fails if any
anchor misses by > 1e-3 (beyond rounding).

Drift guards (hard checks): (1) labels in wholebody_v2.json must equal target_labels.py;
(2) top_share recomputed from the atlas zip must equal the committed values per gene.
Guard (2) is skipped when ABLATION_ALLOW_SYNTHETIC=1 (plumbing tests on synthetic data
only; never publish a run with synthetic_data_mode true).

Secondary contrast (pre-stated): safe vs efficacy_fail (the section 23.4 specificity
control) under the same feature sets.

RUNTIME / REQUIREMENTS
----------------------
Needs ev_capture_data/normal_tissue_atlas/rna_single_cell_type.tsv.zip, which is NOT in
the git clone -- run on the machine that has the atlas. Override its directory with the
ABLATION_DATA_DIR environment variable. Runtime: well under a minute.

OUTPUT
------
  magnitude_baseline_ablation.json   full results (next to this script)
  stdout                             readable tables

This script reads labels and features; it never writes to target_labels.py or any frozen
artifact.
"""
import io
import json
import os
import zipfile
import collections
from datetime import date

import numpy as np

import target_labels as TL
from wholebody_v2 import sep
from auc_confidence_intervals import bootstrap_ci

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'wholebody_v2.json')
OUT = os.path.join(HERE, 'magnitude_baseline_ablation.json')
D = os.environ.get(
    'ABLATION_DATA_DIR',
    os.path.abspath(os.path.join(HERE, '..', '..', '..',
                                 'ev_capture_data', 'normal_tissue_atlas')))
Z = os.path.join(D, 'rna_single_cell_type.tsv.zip')

N_BOOT = 10000
SEED = 0
LEVEL = 95.0
GENES = {t[0] for t in TL.TARGETS}

PUBLISHED_ANCHORS = {
    'toxicity_one_factor_top_share': 0.8804,
    'toxicity_two_factor_top_share': 0.9139,
    'specificity_one_factor_top_share': 0.7456,
}
ANCHOR_TOL = 1e-3


# ----------------------------------------------------------------------------------------
# feature loading (mirrors wholebody_v2.load() parsing exactly)
# ----------------------------------------------------------------------------------------
def load_features(zip_path):
    """Per-gene top_share / mean / max from the single-cell consensus zip.

    mean and max are taken over the cell-type rows PRESENT for that gene -- the same
    rows that enter top_share's denominator -- so all three features see identical data.
    """
    z = zipfile.ZipFile(zip_path)
    name = z.namelist()[0]
    vals = collections.defaultdict(dict)
    with z.open(name) as f:
        fh = io.TextIOWrapper(f, 'utf-8')
        hdr = fh.readline().rstrip('\n').split('\t')
        gi = hdr.index('Gene name')
        ci = next(hdr.index(c) for c in ('Cell type',) if c in hdr)
        vi = next(hdr.index(c) for c in ('nCPM', 'nTPM', 'pTPM') if c in hdr)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if p[gi] in GENES:
                try:
                    vals[p[gi]][p[ci]] = float(p[vi])
                except (ValueError, IndexError):
                    pass
    feats = {}
    for g, d in vals.items():
        v = np.array(list(d.values()), float)
        tot = v.sum()
        if tot <= 0:
            continue
        feats[g] = {'top_share': float(v.max() / tot),
                    'mean': float(v.mean()),
                    'max': float(v.max()),
                    'n_cell_types': int(len(v))}
    return feats


def rank01(values):
    """Average-rank normalization to [0, 1] (ties share the average rank).

    sep() depends only on ranks, so this leaves every ONE-factor result unchanged; it
    exists solely to give mean/max a scale-matched two-factor combination, because the
    published rule adds a 0/1 flag to a feature already on the [0, 1] scale.
    """
    v = np.asarray(values, float)
    n = len(v)
    order = np.argsort(v, kind='mergesort')
    ranks = np.empty(n, float)
    i = 0
    sv = v[order]
    while i < n:
        j = i
        while j + 1 < n and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j)
        i = j + 1
    if n > 1:
        ranks /= (n - 1.0)
    return ranks


# ----------------------------------------------------------------------------------------
# scoring and testing
# ----------------------------------------------------------------------------------------
def split_sep(scores, genes, labels, pos, neg):
    a = [scores[g] for g in genes if labels[g] == pos]
    b = [scores[g] for g in genes if labels[g] == neg]
    return sep(a, b), a, b


def build_feature_sets(genes, labels, haem, committed, feats):
    """All scored feature dictionaries for one contrast, on a common gene list."""
    ts = {g: float(committed[g]['top_share']) for g in genes}
    mn = {g: feats[g]['mean'] for g in genes}
    mx = {g: feats[g]['max'] for g in genes}
    r_mn = dict(zip(genes, rank01([mn[g] for g in genes]).tolist()))
    r_mx = dict(zip(genes, rank01([mx[g] for g in genes]).tolist()))

    def plus_haem(base):
        return {g: base[g] + (1.0 if haem.get(g) else 0.0) for g in genes}

    # post-hoc (2026-10-25): magnitude in its best orientation. One-factor flips negate
    # the feature (exact orientation reversal for a rank-based sep); two-factor flips use
    # 1 - rank01 so the continuous part stays on the [0, 1] scale the flag rule assumes.
    f_mn = {g: -mn[g] for g in genes}
    f_mx = {g: -mx[g] for g in genes}
    fr_mn = {g: 1.0 - r_mn[g] for g in genes}
    fr_mx = {g: 1.0 - r_mx[g] for g in genes}

    return {
        'top_share': ts,                                   # concentration (published feature)
        'mean': mn,                                        # magnitude
        'max': mx,                                         # peak magnitude
        'top_share+haem': plus_haem(ts),                   # published two-factor rule
        'rank_mean+haem': plus_haem(r_mn),                 # scale-matched analogue
        'rank_max+haem': plus_haem(r_mx),                  # scale-matched analogue
        'rank_mean+top_share+haem': {g: r_mn[g] + ts[g] + (1.0 if haem.get(g) else 0.0)
                                     for g in genes},      # exploratory only
        'flipped_mean': f_mn,                              # post-hoc: magnitude, best direction
        'flipped_max': f_mx,                               # post-hoc
        'flipped_rank_mean+haem': plus_haem(fr_mn),        # post-hoc, scale-matched
        'flipped_rank_max+haem': plus_haem(fr_mx),         # post-hoc, scale-matched
    }


def paired_bootstrap(feature_sets, genes, labels, pos, neg,
                     n_boot=N_BOOT, seed=SEED, level=LEVEL):
    """Bootstrap every feature on the SAME resampled targets; return per-feature sep
    draws. Deltas between feature columns are then paired by construction."""
    names = list(feature_sets)
    safe_g = [g for g in genes if labels[g] == pos]
    neg_g = [g for g in genes if labels[g] == neg]
    S = np.array([[feature_sets[n][g] for g in safe_g] for n in names], float)
    T = np.array([[feature_sets[n][g] for g in neg_g] for n in names], float)
    ns, nt = len(safe_g), len(neg_g)
    rng = np.random.default_rng(seed)
    draws = np.empty((n_boot, len(names)), float)
    for i in range(n_boot):
        sa = rng.integers(0, ns, ns)
        to = rng.integers(0, nt, nt)
        A = S[:, sa]                      # F x ns
        B = T[:, to]                      # F x nt
        diff = A[:, :, None] - B[:, None, :]   # F x ns x nt
        gt = (diff > 0).sum(axis=(1, 2))
        eq = (diff == 0).sum(axis=(1, 2))
        draws[i] = (gt + 0.5 * eq) / (ns * nt)
    tail = (100.0 - level) / 2.0
    lo = np.percentile(draws, tail, axis=0)
    hi = np.percentile(draws, 100.0 - tail, axis=0)
    return {'names': names, 'draws': draws,
            'ci_lo': lo.tolist(), 'ci_hi': hi.tolist()}


def delta_stats(draws, names, hi_name, lo_name, point, level=LEVEL):
    """Percentile CI and P(delta <= 0) for sep(hi) - sep(lo) from paired draws.
    `point` is the observed full-sample difference (bootstrap draws give the CI)."""
    i, j = names.index(hi_name), names.index(lo_name)
    d = draws[:, i] - draws[:, j]
    tail = (100.0 - level) / 2.0
    return {'point': float(point),
            'ci_lo': float(np.percentile(d, tail)),
            'ci_hi': float(np.percentile(d, 100.0 - tail)),
            'p_delta_le_0': float((d <= 0).mean())}


# ----------------------------------------------------------------------------------------
def main():
    allow_synthetic = os.environ.get('ABLATION_ALLOW_SYNTHETIC', '') == '1'
    d = json.load(open(SRC))
    labels = d['labels']
    haem = d['haem_home_A']
    committed = d['reference_A_cell_type']

    # drift guard 1: committed labels must equal target_labels.py
    tl_labels = {t[0]: t[3] for t in TL.TARGETS}
    assert labels == tl_labels, 'wholebody_v2.json labels drifted from target_labels.py'

    feats = load_features(Z)
    missing = sorted(GENES - set(feats))
    assert not missing, 'genes missing from the atlas zip: %s' % missing

    # drift guard 2: the zip must reproduce the committed top_share values
    worst_g, worst_v = None, 0.0
    for g in sorted(GENES):
        dev = abs(feats[g]['top_share'] - float(committed[g]['top_share']))
        if dev > worst_v:
            worst_g, worst_v = g, dev
    if allow_synthetic:
        print('SYNTHETIC DATA MODE: top_share drift guard skipped '
              '(worst gene %s dev %.4g -- expected for synthetic data)' % (worst_g, worst_v))
    else:
        assert worst_v < 1e-6, (
            'atlas zip top_share does not match committed values (worst %s dev %.3g) '
            '-- parsing drifted or the wrong file is mounted' % (worst_g, worst_v))

    safe_g = sorted(g for g in labels if labels[g] == 'safe')
    tox_g = sorted(g for g in labels if labels[g] == 'toxic')
    eff_g = sorted(g for g in labels if labels[g] == 'efficacy_fail')
    assert (len(safe_g), len(tox_g), len(eff_g)) == (19, 11, 12), 'outcome groups drifted'

    results = {
        '_meta': {
            'date': date.today().isoformat(),
            'n_boot': N_BOOT, 'seed': SEED, 'level': LEVEL,
            'data_dir': D,
            'synthetic_data_mode': allow_synthetic,
            'method': ('head-to-head feature ablation on the frozen outcome test; '
                       'sep and bootstrap imported from the published pipeline '
                       '(wholebody_v2.sep, auc_confidence_intervals.bootstrap_ci); '
                       'paired deltas share resampled targets'),
            'note': ('Exploratory feature sets are diagnostics, not scoring-rule '
                     'proposals. The published two-factor rule stays frozen.'),
            'posthoc_rows': ('flipped-orientation rows added 2026-10-25 after the first '
                             'run, to test magnitude in its best direction; the '
                             'pre-stated comparisons are unchanged and remain primary'),
        },
        'anchors': {},
        'toxicity_contrast': {},
        'specificity_contrast': {},
        'paired_deltas': {},
        'feature_values': {g: dict(feats[g],
                                   haem_home=bool(haem.get(g)),
                                   label=labels[g]) for g in sorted(GENES)},
    }

    for contrast, pos, neg, groups in (
            ('toxicity_contrast', 'safe', 'toxic', safe_g + tox_g),
            ('specificity_contrast', 'safe', 'efficacy_fail', safe_g + eff_g)):
        fsets = build_feature_sets(groups, labels, haem, committed, feats)
        block = {}
        for name, scores in fsets.items():
            pt, a, b = split_sep(scores, groups, labels, pos, neg)
            # published bootstrap procedure, published seed; bootstrap_ci's own point
            # estimate equals `pt` (same sep on the same full arrays)
            _, lo, hi, deg = bootstrap_ci(a, b)
            block[name] = {
                'point': pt, 'ci_lo': lo, 'ci_hi': hi,
                'n_pos': len(a), 'n_neg': len(b),
                'n_degenerate_resamples': deg,
                'mean_pos': float(np.mean(a)), 'mean_neg': float(np.mean(b)),
            }
        # anchor checks on the published feature
        key = ('toxicity' if contrast == 'toxicity_contrast' else 'specificity')
        for nm, anchor_key in (('top_share', key + '_one_factor_top_share'),
                               ('top_share+haem', key + '_two_factor_top_share')):
            if anchor_key in PUBLISHED_ANCHORS:
                dev = abs(block[nm]['point'] - PUBLISHED_ANCHORS[anchor_key])
                results['anchors'][anchor_key] = {
                    'published': PUBLISHED_ANCHORS[anchor_key],
                    'reproduced': block[nm]['point'], 'abs_dev': dev}
                assert dev < ANCHOR_TOL, (
                    'anchor %s failed: published %.4f vs reproduced %.4f'
                    % (anchor_key, PUBLISHED_ANCHORS[anchor_key], block[nm]['point']))
        pb = paired_bootstrap(fsets, groups, labels, pos, neg)
        results[contrast] = block
        results['paired_deltas'][contrast] = {
            'top_share_minus_mean_one_factor':
                delta_stats(pb['draws'], pb['names'], 'top_share', 'mean',
                            block['top_share']['point'] - block['mean']['point']),
            'top_share_minus_max_one_factor':
                delta_stats(pb['draws'], pb['names'], 'top_share', 'max',
                            block['top_share']['point'] - block['max']['point']),
            'top_share+haem_minus_rank_mean+haem':
                delta_stats(pb['draws'], pb['names'], 'top_share+haem', 'rank_mean+haem',
                            block['top_share+haem']['point'] - block['rank_mean+haem']['point']),
            'top_share+haem_minus_rank_max+haem':
                delta_stats(pb['draws'], pb['names'], 'top_share+haem', 'rank_max+haem',
                            block['top_share+haem']['point'] - block['rank_max+haem']['point']),
            # ---- post-hoc rows (2026-10-25): magnitude at its best orientation ----
            'top_share_minus_flipped_mean_one_factor':
                delta_stats(pb['draws'], pb['names'], 'top_share', 'flipped_mean',
                            block['top_share']['point'] - block['flipped_mean']['point']),
            'top_share_minus_flipped_max_one_factor':
                delta_stats(pb['draws'], pb['names'], 'top_share', 'flipped_max',
                            block['top_share']['point'] - block['flipped_max']['point']),
            'top_share+haem_minus_flipped_rank_mean+haem':
                delta_stats(pb['draws'], pb['names'], 'top_share+haem', 'flipped_rank_mean+haem',
                            block['top_share+haem']['point'] - block['flipped_rank_mean+haem']['point']),
            'top_share+haem_minus_flipped_rank_max+haem':
                delta_stats(pb['draws'], pb['names'], 'top_share+haem', 'flipped_rank_max+haem',
                            block['top_share+haem']['point'] - block['flipped_rank_max+haem']['point']),
        }

    # ------------------------------------------------------------------------------------
    # readable summary
    # ------------------------------------------------------------------------------------
    print('=' * 78)
    print('MAGNITUDE-BASELINE ABLATION  (%d bootstraps, seed %d, %s data)'
          % (N_BOOT, SEED, 'SYNTHETIC' if allow_synthetic else 'atlas'))
    print('=' * 78)
    print('sep(safe, other) = P(safe scores higher); 0.5 = chance. Anchors reproduced:')
    for k, v in results['anchors'].items():
        print('  %-36s published %.4f  reproduced %.4f  (dev %.2g)'
              % (k, v['published'], v['reproduced'], v['abs_dev']))
    for contrast, title in (('toxicity_contrast', 'TOXICITY TEST (safe vs toxic)'),
                            ('specificity_contrast', 'SPECIFICITY CONTROL (safe vs efficacy_fail)')):
        print('\n%s' % title)
        print('  %-26s %7s   %-22s %s' % ('feature', 'sep', '95% CI', 'mean safe / mean other'))
        for name, v in results[contrast].items():
            if 'rank_mean+top_share' in name:
                tag = '   (exploratory)'
            elif name.startswith('flipped'):
                tag = '   (post-hoc flip)'
            else:
                tag = ''
            print('  %-26s %7.4f   [%6.4f, %6.4f]   %.3f / %.3f%s'
                  % (name, v['point'], v['ci_lo'], v['ci_hi'],
                     v['mean_pos'], v['mean_neg'], tag))
        print('  paired deltas (concentration minus magnitude, same resamples):')
        for name, v in results['paired_deltas'][contrast].items():
            print('    %-38s %+.4f  [%+.4f, %+.4f]  P(delta<=0)=%.4f'
                  % (name, v['point'], v['ci_lo'], v['ci_hi'], v['p_delta_le_0']))
    print('\nInterpretation: if the concentration claim is right, top_share rows beat the')
    print('mean/max rows and the paired deltas exclude 0. If magnitude matches top_share,')
    print('the novelty claim must shift from the FEATURE to the panel/label design. Either')
    print('result is informative; report whichever occurs.')
    print('Post-hoc flipped rows (magnitude at its best orientation) are the strongest')
    print('version of the baseline challenge; the pre-stated comparisons remain primary.')
    with open(OUT, 'w') as f:
        json.dump(results, f, indent=1)
    print('\nwrote %s' % OUT)
    return results


if __name__ == '__main__':
    main()
