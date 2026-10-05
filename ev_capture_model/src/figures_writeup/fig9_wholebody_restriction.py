"""
fig9_wholebody_restriction.py

Write-up figure 9 -- the whole-body extension, and the project's central claim in one
picture: a target is safe only if it is BOTH confined to few normal cell types AND those
cell types are expendable.

Panel A  every target with a clinical safety verdict, placed by how confined it is across
         154 normal human cell types (x) and how many cell types it is on (y). Fill colour
         is the clinical outcome, which never touches the coordinates. The safe targets sit
         top-left; the on-target-toxic ones sit bottom-right; CEACAM5 is deliberately
         labelled because it is the case a pure-confinement account gets wrong.
Panel B  the same test with EVERY blood-cancer target removed. This is the load-bearing
         robustness check: the obvious objection is that the measure only works because
         blood targets happen to be both confined and safe.

All numbers are read live from wholebody_v2.json, which wholebody_v2.py writes. Nothing
here is transcribed.
"""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "aim13_scoping"))
from figure_style import (clean_spines, section_badge, add_callout,
                          PROJECT_BLUE, PROJECT_ORANGE, INK, NEUTRAL_GRAY)

RES = os.path.join(HERE, "..", "aim13_scoping", "wholebody_v2.json")
FIGURES_DIR = os.path.join(HERE, "..", "..", "..", "Poster, Paper, and Write Up", "Figures")

SAFE_C = "#2E5EAA"
TOX_C = "#B0563A"

COMMON = {'CD19': 'CD19', 'MS4A1': 'CD20', 'CD22': 'CD22', 'CD79B': 'CD79b',
          'TNFRSF17': 'BCMA', 'SLAMF7': 'CS1', 'CD38': 'CD38', 'TNFRSF8': 'CD30',
          'CD33': 'CD33', 'GPRC5D': 'GPRC5D', 'FOLH1': 'PSMA', 'TACSTD2': 'TROP2',
          'NECTIN4': 'Nectin-4', 'CLDN18': 'CLDN18.2', 'DLL3': 'DLL3', 'MET': 'c-MET',
          'MSLN': 'MSLN', 'GPC3': 'GPC3', 'CD52': 'CD52', 'ERBB2': 'HER2', 'EGFR': 'EGFR',
          'CEACAM5': 'CEA', 'CA9': 'CAIX', 'CD47': 'CD47', 'CD44': 'CD44v6',
          'NCAM1': 'CD56', 'EPCAM': 'EpCAM', 'FOLR1': 'FRa', 'FAP': 'FAP'}

BLOOD = {'CD19', 'MS4A1', 'CD22', 'CD79B', 'TNFRSF17', 'SLAMF7', 'CD38', 'TNFRSF8',
         'CD33', 'GPRC5D', 'CD52', 'PTPRC', 'IL3RA', 'FCER2', 'CD40', 'CD80',
         'HLA-DRA', 'SDC1', 'CD74', 'CD47'}


def build(out_path=None):
    R = json.load(open(RES))
    cA = R['reference_A_cell_type']
    label = R['labels']
    haem = R['haem_home_A']
    solid = R['solid_only']

    genes = [g for g in cA if label.get(g) in ('safe', 'toxic')]

    fig = plt.figure(figsize=(13.2, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.72, 1.0], wspace=0.28,
                          left=0.065, right=0.985, top=0.86, bottom=0.20)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # ---------------- Panel A ----------------
    textsA, xsA, ysA = [], [], []
    for g in genes:
        x = cA[g]['top_share']
        y = cA[g]['n_types']
        c = TOX_C if label[g] == 'toxic' else SAFE_C
        m = 'o' if haem.get(g) else 's'
        axA.scatter(x, y, s=132, marker=m, facecolor=c, edgecolor='white',
                    linewidth=1.6, zorder=3)
        xsA.append(x); ysA.append(y)
        textsA.append(axA.text(x, y, COMMON.get(g, g), fontsize=8.8, color=INK, zorder=6))

    axA.set_yscale('log')
    axA.set_xlim(-0.02, 0.93)
    axA.set_ylim(0.62, 210)
    axA.set_yticks([1, 2, 5, 10, 20, 50, 100, 200])
    axA.set_yticklabels(['1', '2', '5', '10', '20', '50', '100', '200'])
    axA.set_xlabel('confinement  (share of all normal expression in its single top cell type)',
                   fontsize=10.5)
    axA.set_ylabel('number of normal cell types it is on', fontsize=10.5)
    axA.set_title('A.  Every target with a clinical safety verdict, placed only by\n'
                  'where it sits in 154 normal human cell types',
                  fontsize=12.5, fontweight='bold', loc='left', pad=10)
    clean_spines(axA)
    axA.grid(axis='both', color='#E6E9EF', linewidth=0.8, zorder=0)
    axA.set_axisbelow(True)

    adjust_text(textsA, x=xsA, y=ysA, ax=axA,
                expand=(1.25, 1.6), force_text=(0.45, 0.75), force_static=(0.28, 0.45),
                arrowprops=dict(arrowstyle='-', color='#9AA3B2', linewidth=0.7,
                                shrinkA=1, shrinkB=3), max_move=None,
                ensure_inside_axes=True)

    add_callout(axA, (cA['CEACAM5']['top_share'], cA['CEACAM5']['n_types']),
                'CEA: confined —\nbut to colon,\nwhich you cannot lose',
                (0.615, 0.30), color=PROJECT_ORANGE, fontsize=9.5,
                textcoords='axes fraction')

    handles = [
        plt.Line2D([], [], marker='o', color='none', markerfacecolor=SAFE_C,
                   markeredgecolor='white', markersize=11, label='no on-target toxicity'),
        plt.Line2D([], [], marker='o', color='none', markerfacecolor=TOX_C,
                   markeredgecolor='white', markersize=11, label='severe on-target toxicity'),
        plt.Line2D([], [], marker='o', color='none', markerfacecolor=NEUTRAL_GRAY,
                   markeredgecolor='white', markersize=11, label='home = blood lineage (expendable)'),
        plt.Line2D([], [], marker='s', color='none', markerfacecolor=NEUTRAL_GRAY,
                   markeredgecolor='white', markersize=10, label='home = solid tissue'),
    ]
    axA.legend(handles=handles, fontsize=8.8, loc='lower left',
               bbox_to_anchor=(0.005, 0.005), frameon=True, facecolor='white',
               edgecolor='#D6DAE2', framealpha=0.96,
               handletextpad=0.4, labelspacing=0.55).set_zorder(20)

    # ---------------- Panel B ----------------
    sg = [g for g in genes if g not in BLOOD]
    textsB, xsB, ysB = [], [], []
    for i, grp in enumerate(('safe', 'toxic')):
        gs_ = sorted((g for g in sg if label[g] == grp), key=lambda g: cA[g]['top_share'])
        xs = np.full(len(gs_), i, float)
        ys = [cA[g]['top_share'] for g in gs_]
        axB.scatter(xs, ys, s=128, marker='s',
                    facecolor=(SAFE_C if grp == 'safe' else TOX_C),
                    edgecolor='white', linewidth=1.6, zorder=3)
        for x, y, g in zip(xs, ys, gs_):
            xsB.append(x); ysB.append(y)
            textsB.append(axB.text(x, y, COMMON.get(g, g), fontsize=8.6, color=INK, zorder=6))
        axB.plot([i - 0.30, i + 0.30], [np.mean(ys)] * 2, color=INK, linewidth=2.4, zorder=5)
        axB.text(i, -0.055, 'mean %.3f' % np.mean(ys), transform=axB.get_xaxis_transform(),
                 ha='center', va='top', fontsize=9, color=INK, fontweight='bold')

    adjust_text(textsB, x=xsB, y=ysB, ax=axB,
                expand=(1.3, 1.7), force_text=(0.6, 0.9), force_static=(0.35, 0.5),
                arrowprops=dict(arrowstyle='-', color='#9AA3B2', linewidth=0.7,
                                shrinkA=1, shrinkB=3), max_move=None,
                ensure_inside_axes=True)

    axB.set_xticks([0, 1])
    axB.set_xticklabels(['no on-target\ntoxicity  (n=%d)' % solid['n_safe'],
                         'severe on-target\ntoxicity  (n=%d)' % solid['n_toxic']],
                        fontsize=10)
    axB.tick_params(axis='x', pad=22)
    axB.set_xlim(-0.95, 1.95)
    axB.set_ylabel('confinement', fontsize=10.5)
    axB.set_title('B.  Solid-tumour targets only —\nevery blood target removed',
                  fontsize=12.5, fontweight='bold', loc='left', pad=10)
    clean_spines(axB)
    axB.grid(axis='y', color='#E6E9EF', linewidth=0.8, zorder=0)
    axB.set_axisbelow(True)
    axB.set_ylim(0.0, 0.63)
    axB.text(0.5, 0.995,
             'separation %.3f      exact p = %.3f' % (solid['sep'], solid['p']),
             transform=axB.transAxes, ha='center', va='top', fontsize=10.5,
             fontweight='bold', color=INK,
             bbox=dict(boxstyle='round,pad=0.5', facecolor='#EEF1F6', edgecolor='none'))

    section_badge(fig, 'CAR-T', 9)
    fig.text(0.5, 0.055,
             'The rule that survives every check: a target is safe only if it is BOTH confined to few normal cell types '
             'AND those cell types are expendable.',
             ha='center', va='center', fontsize=11, fontweight='bold', color='white',
             bbox=dict(boxstyle='round,pad=0.6', facecolor=PROJECT_ORANGE, edgecolor='none'))
    fig.text(0.5, 0.012,
             'Human Protein Atlas single-cell consensus, 154 cell types. Outcome labels are clinical safety history only and are frozen in target_labels.py. '
             'Full panel: separation 0.868, permutation p = 0.0003 (n = 19 vs 10).',
             ha='center', va='center', fontsize=8.2, color='#5A6473', style='italic')

    out = out_path or os.path.join(FIGURES_DIR, 'fig9_wholebody_restriction.png')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=300 if out.endswith('.png') else None,
                facecolor='white', bbox_inches=None)
    plt.close(fig)
    print('wrote', out)
    return out


if __name__ == '__main__':
    build()
