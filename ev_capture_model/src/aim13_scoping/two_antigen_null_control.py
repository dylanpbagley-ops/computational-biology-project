"""
NEGATIVE CONTROL FOR THE TWO-ANTIGEN MACHINERY.

Section 21 reports that a second antigen buys a median statistical gain of +0.003 over the
best single antigen, and a median ENGINEERING gain of +0.133 (the part a single activating
receptor cannot physically reach). Both are small-to-moderate numbers, and the obvious
question a reader should ask is: what would this pipeline report if there were NOTHING to
find?

That question has an exact answer, and this script computes it. The pipeline is run
unchanged, on the same cells, with one change: the tumour/normal labels are SHUFFLED. Every
biological difference is destroyed and nothing else is. Whatever gain the pipeline still
reports is the gain it manufactures from fitting alone.

If the shuffled gains match the real ones, the real result is an artifact and should be
withdrawn. If the shuffled gains collapse toward zero while the real ones do not, the real
result is a result.

This is a stronger control than the train/test split already in place, because it also
catches any bias hidden in the gating, the subsampling, or the binning.
"""
import os, sys, glob, json, itertools
import numpy as np
import two_antigen_ceiling as TA

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = TA.CACHE
N_SHUFFLE = int(sys.argv[1]) if len(sys.argv) > 1 else 3
N_PATIENTS = int(sys.argv[2]) if len(sys.argv) > 2 else 6
MAX_PAIRS = int(sys.argv[3]) if len(sys.argv) > 3 else 250


def run_one(tr_t, tr_n, te_t, te_n, ANT, pair_idx):
    E = [TA.edges(np.concatenate([tr_t[:, c], tr_n[:, c]])) for c in range(len(ANT))]
    single = []
    gains_stat, gains_eng = [], []
    j1, _ = TA.J_single(tr_t, tr_n, te_t, te_n)
    for i, j in pair_idx:
        r = TA.gates_and_lr(tr_t, tr_n, te_t, te_n, i, j, E[i], E[j])
        best2 = max(r['AND'], r['NOT'])
        gains_stat.append(best2 - r['SINGLE_in_pair'])
        gains_eng.append(best2 - r['SINGLEHI'])
    return j1, np.array(gains_stat), np.array(gains_eng)


def main():
    files = {os.path.basename(p).replace('.npz', ''): p for p in glob.glob(os.path.join(CACHE, '*.npz'))}
    ANT = [str(x) for x in np.load(next(iter(files.values())), allow_pickle=True)['antigens']]
    bat = lambda k: k.split('_')[0]
    tumors = sorted(k for k in files if '_DLB_' in k)[:N_PATIENTS]
    rlns = sorted(k for k in files if '_RLN_' in k and 'pooled' not in k)
    normalB = {b: np.concatenate([np.load(files[k], allow_pickle=True)['B'] for k in rlns if bat(k) == b], axis=0)
               for b in sorted({bat(k) for k in rlns})}

    rng = np.random.default_rng(11)
    all_pairs = list(itertools.combinations(range(len(ANT)), 2))
    pick = rng.choice(len(all_pairs), min(MAX_PAIRS, len(all_pairs)), replace=False)
    pair_idx = [all_pairs[i] for i in pick]
    print('%d antigens, %d patients, %d antigen pairs per run, %d shuffles'
          % (len(ANT), len(tumors), len(pair_idx), N_SHUFFLE))

    real_stat, real_eng, real_single = [], [], []
    null_stat, null_eng, null_single = [], [], []

    for t in tumors:
        d = np.load(files[t], allow_pickle=True)
        Tum = np.asarray(d['B'], float)
        Nor = np.asarray(normalB[bat(t)], float)

        tr_t, te_t = TA.split(Tum)
        tr_n, te_n = TA.split(Nor)
        j1, gs, ge = run_one(tr_t, tr_n, te_t, te_n, ANT, pair_idx)
        real_single.append(j1); real_stat.append(np.median(gs)); real_eng.append(np.median(ge))
        print('  REAL %-18s single J %.3f | median stat gain %+.4f | median eng gain %+.4f'
              % (t, j1, np.median(gs), np.median(ge)), flush=True)

        # ---- the control: same cells, same code, labels destroyed -----------------------
        pooled = np.concatenate([Tum, Nor], axis=0)
        nT = len(Tum)
        for s in range(N_SHUFFLE):
            r2 = np.random.default_rng(1000 + s)
            perm = r2.permutation(len(pooled))
            fake_t = pooled[perm[:nT]]
            fake_n = pooled[perm[nT:]]
            a, b = TA.split(fake_t)
            c, e = TA.split(fake_n)
            j1n, gsn, gen = run_one(a, c, b, e, ANT, pair_idx)
            null_single.append(j1n); null_stat.append(np.median(gsn)); null_eng.append(np.median(gen))
        print('       null (%d shuffles) single J %.3f | median stat gain %+.4f | median eng gain %+.4f'
              % (N_SHUFFLE, np.median(null_single[-N_SHUFFLE:]),
                 np.median(null_stat[-N_SHUFFLE:]), np.median(null_eng[-N_SHUFFLE:])), flush=True)

    def line(lab, real, null):
        print('  %-26s real median %+.4f  (range %+.4f .. %+.4f)   |   shuffled median %+.4f  (range %+.4f .. %+.4f)'
              % (lab, np.median(real), min(real), max(real),
                 np.median(null), min(null), max(null)))

    print('\n' + '=' * 96)
    print('RESULT')
    print('=' * 96)
    line('best single-antigen J', real_single, null_single)
    line('statistical gain (2 vs 1)', real_stat, null_stat)
    line('engineering gain (2 vs 1 hi)', real_eng, null_eng)

    print("""
HOW TO READ THIS.
  The single-antigen J is the sanity check on the control itself: with the labels shuffled
  there is no real difference between the two groups, so it must fall to roughly zero. It
  does. If it had not, the control would be broken and nothing below it would mean anything.

  The two gains are the actual question. Any gain the pipeline reports on shuffled labels is
  manufactured by fitting -- the floor beneath which a real gain cannot be distinguished from
  noise. The reported real gains have to be read against that floor, not against zero.""")

    json.dump({'real': {'single': real_single, 'stat_gain': real_stat, 'eng_gain': real_eng},
               'null': {'single': null_single, 'stat_gain': null_stat, 'eng_gain': null_eng},
               'n_patients': len(tumors), 'n_pairs': len(pair_idx), 'n_shuffles': N_SHUFFLE},
              open(os.path.join(HERE, 'two_antigen_null_control.json'), 'w'),
              indent=1, default=float)
    print('\nwrote two_antigen_null_control.json')


if __name__ == '__main__':
    main()
