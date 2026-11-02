"""
Confinement-ranking figure: the project's core quantity for the whole frozen panel.

One horizontal bar per scored target (19 tolerated / 11 toxic), sorted by top_share —
the largest single normal cell type's share of the gene's total normal expression, the
exact quantity frozen in confinement_table.tsv. The diffuse-vs-confined contrast that
the confinement-map pair figure tries to show appears here as the panel-wide pattern:
toxic targets cluster near zero, tolerated targets reach high. CD47 and MS4A1 are
annotated as the anchor examples.

Data: reads ONLY the two frozen artifacts (confinement_table.tsv,
target_labels_frozen.tsv) — no atlas download needed, runs in under a second.
Efficacy-fail and excluded targets are omitted (efficacy is a different endpoint);
note this in the caption.

Figure is DESCRIPTIVE of the frozen panel — it is not a prediction and must be
captioned accordingly. top_share is one component of the confinement score, not the
full two-factor score.

Usage (from ev_capture_model/src/aim13_scoping/):
    python make_confinement_ranking.py

Requires: matplotlib (pip install matplotlib). No other dependencies.
"""
import os, csv
import matplotlib

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, '..', '..', 'figures')
GRAY, BLUE = '#8a8a8a', '#0279EE'


def main():
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    plt.rcParams['font.family'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
    plt.rcParams['svg.fonttype'] = 'none'

    conf = {}
    with open(os.path.join(HERE, 'confinement_table.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            conf[r['gene']] = r
    labels = {}
    with open(os.path.join(HERE, 'target_labels_frozen.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            labels[r['gene']] = r['outcome']

    rows, missing = [], []
    for g, oc in labels.items():
        if oc in ('safe', 'toxic'):
            (rows if g in conf else missing).append(
                (g, oc, float(conf[g]['top_share']), conf[g]['home_cell_type'])
                if g in conf else g)
    if missing:
        print('WARNING — no confinement_table row for:', ', '.join(missing))
    rows.sort(key=lambda x: x[2])
    n_toxic = sum(1 for r in rows if r[1] == 'toxic')
    print('plotted %d targets (%d toxic / %d tolerated)'
          % (len(rows), n_toxic, len(rows) - n_toxic))

    fig, ax = plt.subplots(figsize=(8.5, 9.5))
    fig.subplots_adjust(top=0.94, left=0.16, right=0.97)
    ys = list(range(len(rows)))
    vals = [r[2] * 100 for r in rows]
    cols = [BLUE if r[1] == 'toxic' else GRAY for r in rows]
    ax.barh(ys, vals, color=cols, height=0.68)
    for y, r, v in zip(ys, rows, vals):
        if r[0] in ('CD47', 'MS4A1'):   # anchor examples of the confinement story
            ax.text(v + 1.2, y, '%s — %.1f%% (%s)' % (r[0], v, r[3]),
                    va='center', fontsize=9, color='#333333', fontweight='bold')
        else:
            ax.text(v + 1.2, y, '%.1f' % v if v >= 1 else '%.2f' % v,
                    va='center', fontsize=8, color='#666666')
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_xlabel("top_share: largest single cell type's share of the gene's "
                  "total normal expression (%)", fontsize=9.5)
    ax.set_title('Target confinement: 30 frozen antibody-target candidates',
                 fontsize=12.5, pad=12)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(handles=[Patch(color=BLUE, label='toxic (%d)' % n_toxic),
                       Patch(color=GRAY, label='tolerated (%d)' % (len(rows) - n_toxic))],
              loc='lower right', fontsize=9, frameon=False)

    out_dir = os.path.normpath(OUT_DIR)
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('svg', 'png'):
        out = os.path.join(out_dir, 'confinement_ranking.%s' % ext)
        try:
            fig.savefig(out, dpi=300, bbox_inches='tight')
            print('wrote', out)
        except OSError as e:
            alt = os.path.join(out_dir, 'confinement_ranking_new.%s' % ext)
            print('could not write %s (%s) — writing %s instead' % (out, e, alt))
            fig.savefig(alt, dpi=300, bbox_inches='tight')
            print('wrote', alt)


if __name__ == '__main__':
    main()
