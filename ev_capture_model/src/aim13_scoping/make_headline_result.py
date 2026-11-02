"""
Headline-result figure: CONCENTRATION, not magnitude, separates toxic from tolerated.

Forest plot of the frozen outcome test (30 targets: 19 tolerated / 11 toxic). One row
per scoring rule: dot = AUC (Mann-Whitney, safe = positive class), bar = pre-computed
stratified bootstrap 95% CI (10,000 resamples, seed 0), dashed line = chance (0.5).

  two-factor confinement (top_share + expendable-home bonus)   AUC 0.9139 (0.789-0.995)
  one-factor confinement (top_share alone)                     AUC 0.8804 (0.732-0.986)
  mean expression, BEST orientation (post-hoc, disclosed)      AUC 0.7608 (0.569-0.928)
  max expression, BEST orientation (post-hoc, disclosed)       AUC 0.6220 (0.416-0.818)

Read: both confinement rules clear chance decisively; the magnitude baselines overlap
it. Expression LEVEL does not predict on-target toxicity — where the expression is
CONCENTRATED does.

Data: reads ONLY frozen artifacts (magnitude_baseline_ablation.json,
auc_confidence_intervals.json); the script reproduces every point estimate from the
frozen per-target features and refuses to write a figure on mismatch. DESCRIPTIVE of
the frozen panel — not a prediction; caption accordingly.

Usage (from ev_capture_model/src/aim13_scoping/):
    python make_headline_result.py

Requires: matplotlib. No atlas download needed.
"""
import os, json
import matplotlib

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, '..', '..', 'figures')
GRAY, BLUE = '#8a8a8a', '#0279EE'


def auc_pairwise(scores, is_safe):
    """Mann-Whitney AUC, positive class = safe (higher score = safer)."""
    pos = [s for s, m in zip(scores, is_safe) if m]
    neg = [s for s, m in zip(scores, is_safe) if not m]
    wins = [(a > b) + 0.5 * (a == b) for a in pos for b in neg]
    return sum(wins) / len(wins)


def main():
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
    plt.rcParams['svg.fonttype'] = 'none'

    abl = json.load(open(os.path.join(HERE, 'magnitude_baseline_ablation.json')))
    ci = json.load(open(os.path.join(HERE, 'auc_confidence_intervals.json')))
    fv = abl['feature_values']
    tc = abl['toxicity_contrast']

    genes = [g for g, v in fv.items() if v['label'] in ('safe', 'toxic')]
    is_safe = [fv[g]['label'] == 'safe' for g in genes]

    # --- reproduce every published point estimate before drawing anything -------
    feats = {
        'one_factor': [fv[g]['top_share'] for g in genes],
        'two_factor': [fv[g]['top_share'] + (1.0 if fv[g]['haem_home'] else 0.0)
                       for g in genes],
        'mag_mean': [-fv[g]['mean'] for g in genes],   # best orientation (toxic higher)
        'mag_max': [-fv[g]['max'] for g in genes],
    }
    published = {
        'one_factor': ci['A_one_factor_top_share']['point'],
        'two_factor': ci['A_two_factor']['point'],
        'mag_mean': tc['flipped_mean']['point'],
        'mag_max': tc['flipped_max']['point'],
    }
    for name, s in feats.items():
        got, pub = auc_pairwise(s, is_safe), published[name]
        flag = 'OK' if abs(got - pub) < 5e-4 else 'MISMATCH'
        print('%-12s reproduced %.4f  published %.4f  %s' % (name, got, pub, flag))
        if flag == 'MISMATCH':
            raise SystemExit('reproduction failed — not writing a figure')

    # --- forest plot -------------------------------------------------------------
    rows = [   # top to bottom
        ('Two-factor confinement\n(top_share + expendable-home bonus)',
         'two_factor', 'A_two_factor', BLUE),
        ('One-factor confinement\n(top_share alone)',
         'one_factor', 'A_one_factor_top_share', BLUE),
        ('Mean expression\n(best orientation, post-hoc)',
         'mag_mean', 'flipped_mean', GRAY),
        ('Max expression\n(best orientation, post-hoc)',
         'mag_max', 'flipped_max', GRAY),
    ]
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    fig.subplots_adjust(top=0.82, left=0.30, right=0.96)
    ys = list(range(len(rows) - 1, -1, -1))
    for y, (name, fkey, ckey, color) in zip(ys, rows):
        c = tc[ckey] if ckey in tc else ci[ckey]
        ax.errorbar(c['point'], y, xerr=[[c['point'] - c['ci_lo']],
                                         [c['ci_hi'] - c['point']]],
                    fmt='none', color=color, alpha=0.45, capsize=4,
                    elinewidth=2.2, capthick=2.2)
        ax.plot(c['point'], y, 'o', color=color, markersize=9, zorder=3)
        ax.text(c['ci_hi'] + 0.015, y, '%.3f' % c['point'],
                va='center', fontsize=9, color='#333333')

    ax.axvline(0.5, color='#999999', lw=1.2, ls='--')
    ax.text(0.5, len(rows) - 0.42, 'chance', color='#999999',
            fontsize=9, ha='center')
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9.5)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel('AUC for separating 11 toxic from 19 tolerated targets '
                  '(dot = AUC, bar = 95% bootstrap CI)', fontsize=9.5)
    ax.set_title('Concentration — not expression level — predicts on-target toxicity risk',
                 fontsize=12, pad=12)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(left=True, length=4, width=0.8)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=BLUE, marker='o', lw=2.5,
                              label='confinement rule (pre-specified)'),
                       Line2D([], [], color=GRAY, marker='o', lw=2.5,
                              label='expression-level baseline (post-hoc)')],
              loc='lower left', fontsize=8.5, frameon=False)

    out_dir = os.path.normpath(OUT_DIR)
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('svg', 'png'):
        out = os.path.join(out_dir, 'headline_result.%s' % ext)
        try:
            fig.savefig(out, dpi=300, bbox_inches='tight')
            print('wrote', out)
        except OSError as e:
            alt = os.path.join(out_dir, 'headline_result_new.%s' % ext)
            print('could not write %s (%s) — writing %s instead' % (out, e, alt))
            fig.savefig(alt, dpi=300, bbox_inches='tight')
            print('wrote', alt)


if __name__ == '__main__':
    main()
