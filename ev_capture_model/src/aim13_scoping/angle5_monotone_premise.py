"""Angle 5: do Cas9 variants move ALONG one shared curve, or do they change the curve itself?

Implements exactly the tests in angle5_precommitment.md (T1, T2, T3), written before any of them
was run. Data: NAR 2023 (gkad702) Table S10 -- 306 guides x up to 7 designed mismatch sites,
wild-type SpCas9 vs HiFi vs LZ3. See the precommitment for the reasoning and for how each result
will be interpreted.

Isotonic (monotone non-decreasing) regression is implemented directly (pool-adjacent-violators)
rather than imported, so the module has no dependency beyond numpy/pandas and the algorithm is
visible. Cross-validation is grouped by GUIDE so no guide contributes to both fitting and testing.
"""
import json
import os
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'crispr_offtarget', 'nar2023'))
OUT = os.path.join(HERE, 'angle5_monotone_premise.json')
SEED = 20260928
N_FOLDS = 5
N_BOOT = 5000
VARIANTS = ('hf', 'lz')
# (source, destination, precommitted?). WT->HiFi and WT->LZ3 are the tests written into
# angle5_precommitment.md. HiFi->LZ3 was added AFTER the first pass, because the first pass
# showed the WT on-target rate is compressed (saturated), which makes WT a poor x-axis; it is
# reported alongside, clearly labelled as post-hoc, never in place of the precommitted result.
PAIRS = (('wt', 'hf', True), ('wt', 'lz', True), ('hf', 'lz', False))


def load():
    x = pd.read_excel(os.path.join(DATA, 'TableS10(Off-on-ratios_combined).xlsx'), header=None)
    d = x.iloc[2:].copy()
    d.columns = ['sgRNA', 'on', 'off', 'mm', 'wt_on', 'wt_off', 'wt_r', 'hf_on', 'hf_off', 'hf_r',
                 'lz_on', 'lz_off', 'lz_r']
    for c in d.columns[4:]:
        d[c] = pd.to_numeric(d[c], errors='coerce')
    d['mmclass'] = d['mm'].astype(str).str.extract(r'^(M\d)')[0]
    return d.reset_index(drop=True)


def pav(x, y):
    """Isotonic (non-decreasing) fit. Returns (unique sorted x, fitted values)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    order = np.argsort(x, kind='mergesort')
    x, y = x[order], y[order]
    ux, inv = np.unique(x, return_inverse=True)
    w = np.bincount(inv).astype(float)
    m = np.bincount(inv, weights=y) / w
    # pool adjacent violators over the unique-x means
    vals, wts, cnt = [], [], []
    for mi, wi in zip(m, w):
        vals.append(mi); wts.append(wi); cnt.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            v = (vals[-2] * wts[-2] + vals[-1] * wts[-1]) / (wts[-2] + wts[-1])
            wt = wts[-2] + wts[-1]
            c = cnt[-2] + cnt[-1]
            vals[-2:], wts[-2:], cnt[-2:] = [v], [wt], [c]
    fitted = np.repeat(vals, cnt)
    return ux, fitted


def predict(model, xnew):
    ux, fit = model
    return np.interp(xnew, ux, fit)


def r2(y, yhat):
    y, yhat = np.asarray(y, float), np.asarray(yhat, float)
    sst = ((y - y.mean()) ** 2).sum()
    return float(1 - ((y - yhat) ** 2).sum() / sst) if sst > 0 else float('nan')


def folds(guides, seed=SEED, k=N_FOLDS):
    ug = np.array(sorted(set(guides)))
    rng = np.random.default_rng(seed)
    rng.shuffle(ug)
    fid = {g: i % k for i, g in enumerate(ug)}
    return np.array([fid[g] for g in guides])


def cv_isotonic(x, y, guides):
    f = folds(guides)
    pred = np.empty(len(y))
    for k in range(N_FOLDS):
        tr, te = f != k, f == k
        pred[te] = predict(pav(x[tr], y[tr]), x[te])
    return pred


def spearman(a, b):
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])


def t1_on_target(d):
    g = d.groupby('sgRNA').first().reset_index()
    out = {}
    for src, v, pre in PAIRS:
        key = f'{src}->{v}'
        x, y = g[f'{src}_on'].to_numpy(), g[f'{v}_on'].to_numpy()
        ok = np.isfinite(x) & np.isfinite(y)
        x, y, gg = x[ok], y[ok], g['sgRNA'].to_numpy()[ok]
        pred = cv_isotonic(x, y, gg)
        out[key] = {'precommitted': pre, 'n_guides': int(len(x)), 'cv_r2': r2(y, pred),
                    'spearman': spearman(x, y), 'mean_ratio_dst_over_src': float(np.mean(y) / np.mean(x))}
    return out


def t2_off_target(d):
    out = {}
    for src, v, pre in PAIRS:
        block = {'precommitted': pre}
        for name, sel in (('all', np.ones(len(d), bool)),
                          ('M1', (d['mmclass'] == 'M1').to_numpy()),
                          ('M2', (d['mmclass'] == 'M2').to_numpy()),
                          ('M3', (d['mmclass'] == 'M3').to_numpy())):
            s = d[sel]
            x, y = s[f'{src}_off'].to_numpy(), s[f'{v}_off'].to_numpy()
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 30:
                continue
            x, y, gg = x[ok], y[ok], s['sgRNA'].to_numpy()[ok]
            pred = cv_isotonic(x, y, gg)
            block[name] = {'n_pairs': int(len(x)), 'cv_r2': r2(y, pred), 'spearman': spearman(x, y)}
        out[f'{src}->{v}'] = block
    return out


def t3_same_curve(d):
    """Fit the map on OFF-target pairs from training guides; predict the ON-target variant rate of
    held-out guides; look at the residual where WT rates of both site types overlap."""
    g = d.groupby('sgRNA').first().reset_index()
    f_on = dict(zip(g['sgRNA'], folds(g['sgRNA'].to_numpy())))
    out = {}
    rng = np.random.default_rng(SEED)
    for src, v, pre in PAIRS:
        offx, offy = d[f'{src}_off'].to_numpy(), d[f'{v}_off'].to_numpy()
        ok = np.isfinite(offx) & np.isfinite(offy)
        d2 = d[ok]
        gfold_off = d2['sgRNA'].map(f_on).to_numpy()
        pred_on = np.full(len(g), np.nan)
        for k in range(N_FOLDS):
            tr = gfold_off != k
            model = pav(d2[f'{src}_off'].to_numpy()[tr], d2[f'{v}_off'].to_numpy()[tr])
            te = np.array([f_on[s] == k for s in g['sgRNA']])
            pred_on[te] = predict(model, g[f'{src}_on'].to_numpy()[te])
        obs = g[f'{v}_on'].to_numpy()
        wt_on = g[f'{src}_on'].to_numpy()
        lo = max(np.nanpercentile(wt_on, 5), np.nanpercentile(d[f'{src}_off'], 5))
        hi = min(np.nanpercentile(wt_on, 95), np.nanpercentile(d[f'{src}_off'], 95))
        in_rng = (wt_on >= lo) & (wt_on <= hi) & np.isfinite(obs) & np.isfinite(pred_on)
        resid = (obs - pred_on)[in_rng]
        boots = [np.median(rng.choice(resid, len(resid), replace=True)) for _ in range(N_BOOT)] \
            if len(resid) else []
        out[f'{src}->{v}'] = {'precommitted': pre, 'overlap_range_src_rate': [float(lo), float(hi)],
                  'n_guides_in_range': int(in_rng.sum()),
                  'median_signed_residual': float(np.median(resid)) if len(resid) else None,
                  'ci95': [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))] if boots else None,
                  'frac_positive': float(np.mean(resid > 0)) if len(resid) else None}
    return out


def compute_all():
    d = load()
    return {'n_guides': int(d['sgRNA'].nunique()), 'n_pairs': int(len(d)),
            'T1_on_target_one_map': t1_on_target(d),
            'T2_off_target_one_map': t2_off_target(d),
            'T3_on_vs_off_same_curve': t3_same_curve(d),
            'meta': {'seed': SEED, 'folds': N_FOLDS, 'n_boot': N_BOOT,
                     'precommitment': 'angle5_precommitment.md'}}


def print_summary(r):
    print('=' * 92)
    print(f"ANGLE 5 -- do variants move along one curve?  ({r['n_guides']} guides, {r['n_pairs']} pairs)")
    print('=' * 92)
    tag = lambda x: 'precommitted' if x['precommitted'] else 'POST-HOC'
    for k, x in r['T1_on_target_one_map'].items():
        print(f"T1 on-target  {k:7s} [{tag(x)}]: CV R2 {x['cv_r2']:.3f}  Spearman {x['spearman']:.3f}  (n={x['n_guides']})")
    for k, blk in r['T2_off_target_one_map'].items():
        print(f"T2 off-target {k:7s} [{tag(blk)}]: " + '  '.join(
            f"{c}: R2 {x['cv_r2']:.3f} rho {x['spearman']:.3f} (n={x['n_pairs']})"
            for c, x in blk.items() if c != 'precommitted'))
    for k, x in r['T3_on_vs_off_same_curve'].items():
        print(f"T3 same curve {k:7s} [{tag(x)}]: median residual {x['median_signed_residual']:+.4f} "
              f"95% CI [{x['ci95'][0]:+.4f}, {x['ci95'][1]:+.4f}]  frac>0 {x['frac_positive']:.2f}  "
              f"(n={x['n_guides_in_range']}, src-rate range {x['overlap_range_src_rate'][0]:.3f}-{x['overlap_range_src_rate'][1]:.3f})")


if __name__ == '__main__':
    r = compute_all()
    print_summary(r)
    json.dump(r, open(OUT, 'w'), indent=1)
    print('\nwrote', OUT)
