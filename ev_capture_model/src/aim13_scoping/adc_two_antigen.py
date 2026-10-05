"""
THE TWO-ANTIGEN ANALYSIS, IN A SOLID TUMOUR AGAINST MATCHED NORMAL TISSUE.

Sections 17, 21 and 25 established the two-antigen result in blood cancers:
  - the hierarchy SINGLE <= max(AND, NOT) <= LIKELIHOOD RATIO holds with 0 violations;
  - the STATISTICAL gain from a second antigen is at the shuffled-label noise floor;
  - the ENGINEERING gain -- what a second receptor unlocks that no single ACTIVATING receptor
    can reach, because a CAR cannot fire on the absence of a protein -- is ~130x that floor.

This script repeats all of it in metastatic breast cancer, against normal tissue taken from the
same autopsies (De Schepper et al., Nature Communications 2025). 72 ADC targets, 2,556 antigen
pairs, tumour and normal from the same bodies.

RECONSTRUCTION NOTE. The per-sample matrix is rebuilt from the long-format table, which is
gene-major with a consistent replicate count in all 315 (patient, kind, tissue) blocks. The
reconstruction is verified against the paper's own wide table before anything is computed; the
check is a hard assertion, not a printed warning.

SCOPE NOTE, stated first. These are BULK samples. The separation is between tumour SAMPLES and
normal SAMPLES, not tumour CELLS and normal CELLS, so these numbers are not comparable to the
CyTOF results in sections 17-21 and are reported as their own quantity. There are also only 96
normal samples, so the held-out halves are small; the shuffled-label null is what makes the
numbers readable, and it is run here too.
"""
import os, csv, json, itertools, collections
import numpy as np
import two_antigen_ceiling as TA

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'adc_autopsy'))
LONG = os.path.join(D, 'supfig8_long.tsv')
rng = np.random.default_rng(5)


def build():
    blocks = collections.defaultdict(lambda: collections.defaultdict(list))
    for row in csv.DictReader(open(LONG, encoding='utf-8'), delimiter='\t'):
        try:
            v = float(row['value'])
        except (ValueError, TypeError):
            v = np.nan
        blocks[(row['Patient_ID'], row['samplekind'], row['tissue'])][row['variable']].append(v)
    genes = sorted(next(iter(blocks.values())))
    rows, meta = [], []
    for key, d in blocks.items():
        n = len(d[genes[0]])
        assert all(len(d[g]) == n for g in genes), key
        for k in range(n):
            rows.append([d[g][k] for g in genes])
            meta.append(key)
    return np.array(rows, float), genes, meta


def verify(X, genes, meta):
    """Hard check against the paper's own wide table."""
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(D, '41467_2025_67840_MOESM13_ESM.xlsx'), read_only=True)
    ws = wb['Sup fig 3']
    it = ws.iter_rows(values_only=True)
    hdr = next(it)
    gcols = {g: i for i, g in enumerate(hdr) if g in genes}
    ours = collections.defaultdict(list)
    for r, m in zip(X, meta):
        if m[1] == 'Metastatic':
            ours[(m[0], m[2])].append(r)
    checked = matched = 0
    for r in it:
        if r[0] is None:
            continue
        key = (r[0], r[1])
        if key not in ours:
            continue
        for g, ci in list(gcols.items())[:8]:
            try:
                want = float(r[ci])
            except (TypeError, ValueError):
                continue
            gi = genes.index(g)
            checked += 1
            if any(abs(v[gi] - want) < 1e-6 for v in ours[key]):
                matched += 1
    wb.close()
    assert checked > 200 and matched / checked > 0.98, \
        'reconstruction check failed: %d of %d values matched' % (matched, checked)
    print('  reconstruction verified against the paper\'s own wide table: %d/%d values matched'
          % (matched, checked))


def split(a):
    p = rng.permutation(len(a))
    h = max(1, len(a) // 2)
    return a[p[:h]], a[p[h:]]


def run(T, N, genes, pairs):
    tr_t, te_t = split(T); tr_n, te_n = split(N)
    E = [TA.edges(np.concatenate([tr_t[:, c], tr_n[:, c]])) for c in range(len(genes))]
    j1, _ = TA.J_single(tr_t, tr_n, te_t, te_n)
    viol = 0
    stat, eng = [], []
    for i, j in pairs:
        r = TA.gates_and_lr(tr_t, tr_n, te_t, te_n, i, j, E[i], E[j])
        # the theorem's ordering is an IN-SAMPLE statement, so it is checked in-sample
        if max(r['AND_train'], r['NOT_train']) > r['LR_train'] + 1e-9:
            viol += 1
        if r['SINGLE_train'] > max(r['AND_train'], r['NOT_train']) + 1e-9:
            viol += 1
        best2 = max(r['AND'], r['NOT'])
        stat.append(best2 - r['SINGLE_in_pair'])
        eng.append(best2 - r['SINGLEHI'])
    return j1, np.array(stat), np.array(eng), viol


def main():
    X, genes, meta = build()
    print('reconstructed matrix: %d samples x %d genes' % X.shape)
    verify(X, genes, meta)
    keep = ~np.isnan(X).any(1)
    X, meta = X[keep], [m for m, k in zip(meta, keep) if k]
    kind = np.array([m[1] for m in meta])
    T, N = X[kind == 'Metastatic'], X[kind == 'Normal']
    print('  usable: %d metastatic samples, %d normal samples' % (len(T), len(N)))

    all_pairs = list(itertools.combinations(range(len(genes)), 2))
    print('  antigen pairs: %d' % len(all_pairs))

    j1, stat, eng, viol = run(T, N, genes, all_pairs)
    print('\n' + '=' * 78)
    print('REAL DATA')
    print('=' * 78)
    print('  theorem violations (in-sample, both orderings): %d of %d checks'
          % (viol, 2 * len(all_pairs)))
    print('  best single antigen, held out:            J = %.3f' % j1)
    print('  statistical gain from a 2nd antigen:      median %+.4f  (90th pct %+.4f)'
          % (np.median(stat), np.percentile(stat, 90)))
    print('  engineering gain from a 2nd antigen:      median %+.4f  (90th pct %+.4f)'
          % (np.median(eng), np.percentile(eng, 90)))

    print('\n' + '=' * 78)
    print('SHUFFLED-LABEL NULL (same samples, same code, tumour/normal labels destroyed)')
    print('=' * 78)
    pooled = np.concatenate([T, N])
    nT = len(T)
    ns, ne, nj = [], [], []
    sub = [all_pairs[i] for i in rng.choice(len(all_pairs), 400, replace=False)]
    for s in range(3):
        r2 = np.random.default_rng(200 + s)
        perm = r2.permutation(len(pooled))
        a, b, c, viol2 = run(pooled[perm[:nT]], pooled[perm[nT:]], genes, sub)
        nj.append(a); ns.append(np.median(b)); ne.append(np.median(c))
    print('  best single antigen:                      J = %.3f   (real %.3f)' % (np.median(nj), j1))
    print('  statistical gain:                         median %+.4f   (real %+.4f)'
          % (np.median(ns), np.median(stat)))
    print('  engineering gain:                         median %+.4f   (real %+.4f)'
          % (np.median(ne), np.median(eng)))

    print("""
READING. The blood-cancer result was: no new information from a second antigen (statistical
gain at the noise floor), but a large amount of information unlocked that a single ACTIVATING
receptor cannot use (engineering gain far above it). Whether the same split appears here is the
question, and the null column is the only thing that makes either number interpretable at 96
normal samples.""")

    json.dump({'n_tumour': int(len(T)), 'n_normal': int(len(N)), 'n_pairs': len(all_pairs),
               'violations': viol, 'single_J': j1,
               'stat_median': float(np.median(stat)), 'eng_median': float(np.median(eng)),
               'null_single': float(np.median(nj)), 'null_stat': float(np.median(ns)),
               'null_eng': float(np.median(ne))},
              open(os.path.join(HERE, 'adc_two_antigen.json'), 'w'), indent=1)
    print('\nwrote adc_two_antigen.json')


if __name__ == '__main__':
    main()


def top_pairs(k=12):
    """Which pairs carry the engineering gain, and what would you actually build?"""
    X, genes, meta = build()
    keep = ~np.isnan(X).any(1)
    X, meta = X[keep], [m for m, kk in zip(meta, keep) if kk]
    kind = np.array([m[1] for m in meta])
    T, N = X[kind == 'Metastatic'], X[kind == 'Normal']
    tr_t, te_t = split(T); tr_n, te_n = split(N)
    E = [TA.edges(np.concatenate([tr_t[:, c], tr_n[:, c]])) for c in range(len(genes))]
    out = []
    for i, j in itertools.combinations(range(len(genes)), 2):
        r = TA.gates_and_lr(tr_t, tr_n, te_t, te_n, i, j, E[i], E[j])
        best2 = max(r['AND'], r['NOT'])
        out.append((best2 - r['SINGLEHI'], best2, r['SINGLEHI'], r['SINGLE_in_pair'],
                    genes[i], genes[j], 'NOT' if r['NOT'] >= r['AND'] else 'AND'))
    out.sort(reverse=True)
    print('\nthe %d pairs with the largest ENGINEERING gain (held out):' % k)
    print('  %-9s %-9s %5s %8s %8s %8s' % ('antigen A', 'antigen B', 'gate', 'gain', '2-receptor', '1-receptor'))
    for g, b2, hi, anyd, a, bb, gate in out[:k]:
        print('  %-9s %-9s %5s %+8.3f %8.3f %8.3f' % (a, bb, gate, g, b2, hi))
    frac = sum(1 for o in out if o[0] > 0.05) / len(out)
    print('\n  %.1f%% of the %d pairs show an engineering gain above 0.05.' % (100 * frac, len(out)))
    json.dump([{'A': a, 'B': b, 'gate': gate, 'gain': g, 'two_receptor': b2, 'one_receptor': hi}
               for g, b2, hi, anyd, a, b, gate in out[:40]],
              open(os.path.join(HERE, 'adc_two_antigen_top_pairs.json'), 'w'), indent=1)


if __name__ == '__main__' and os.environ.get('TOP'):
    top_pairs()
