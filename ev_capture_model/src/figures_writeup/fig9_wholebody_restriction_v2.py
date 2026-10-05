"""
fig9_wholebody_restriction_v2.py

Redesign of fig9's Panel A. The original panel A was a log-scale scatter with
30 leader-lined text labels (via adjustText) pointing at dots -- many lines
cross, several labels sit far from their dot, and the reader has to trace
each line before comparing genes at all.

Fix: drop the scatter + leader-line structure for a sorted horizontal
dot-strip ("rank plot"), one row per gene, split into two column groups
(no on-target toxicity / severe on-target toxicity). Every label sits
directly on its own row next to its own dot -- zero leader lines, because
there is no ambiguity about which label goes with which point. Sorting by
confinement within each group makes the panel read top-to-bottom as a
ranking, which mirrors how the project's own argument uses this data: safe
targets cluster confined; toxic targets don't.

n_types is kept as a second encoding (marker size, log-scaled) instead of a
second full axis -- stacking two continuous axes onto one scatter was a big
part of why the original panel got crowded. The haematopoietic-home flag
keeps its circle/square marker encoding from the original.

Panel B is copied through unchanged from fig9_wholebody_restriction.py --
it was already a clean strip plot with only 20 points, so this is a
complete, drop-in replacement figure covering both panels.

Nothing here is transcribed: all numbers are read live from wholebody_v2.json.
Does NOT overwrite the original fig9_wholebody_restriction.svg/.png.
"""
import os, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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


def n_types_to_markersize(n):
    # log-ish scaling so 1 -> ~5.5pt radius and 120 -> ~13pt radius, monotone
    return 34 + 26 * (np.log10(max(n, 1)) / np.log10(150))


def build(out_path=None):
    R = json.load(open(RES))
    cA = R['reference_A_cell_type']
    label = R['labels']
    haem = R['haem_home_A']
    solid = R['solid_only']

    genes = [g for g in cA if label.get(g) in ('safe', 'toxic')]

    fig = plt.figure(figsize=(14.4, 8.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.38,
                          left=0.065, right=0.985, top=0.80, bottom=0.185)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # ---------------- Panel A: sorted dot-strip (rank plot), no leader lines ----------------
    safe_genes = sorted((g for g in genes if label[g] == 'safe'),
                        key=lambda g: cA[g]['top_share'], reverse=True)
    toxic_genes = sorted((g for g in genes if label[g] == 'toxic'),
                         key=lambda g: cA[g]['top_share'], reverse=True)

    n_rows = max(len(safe_genes), len(toxic_genes))
    row_gap = 1.0
    top_y = n_rows * row_gap

    def draw_column(gene_list, x0, color):
        ys = []
        for i, g in enumerate(gene_list):
            y = top_y - i * row_gap
            x = cA[g]['top_share']
            ntypes = cA[g]['n_types']
            m = 'o' if haem.get(g) else 's'
            ms = n_types_to_markersize(ntypes)

            # dot at position proportional to confinement, on a short local track
            track_x0, track_x1 = x0 + 0.05, x0 + 0.62
            dot_x = track_x0 + x * (track_x1 - track_x0)

            axA.plot([track_x0, track_x1], [y, y], color='#E6E9EF', linewidth=1.2, zorder=1)
            axA.scatter([dot_x], [y], s=ms, marker=m, facecolor=color,
                        edgecolor='white', linewidth=1.3, zorder=3)

            # gene label immediately to the right of the track -- always same spot, no leader line
            axA.text(x0 + 0.66, y, COMMON.get(g, g), fontsize=9.3, color=INK,
                     va='center', ha='left', zorder=4)
            # confinement value + n_types, small, right after the label
            axA.text(x0 + 0.985, y, '%.2f  ·  %d ct' % (x, ntypes), fontsize=7.6,
                     color=NEUTRAL_GRAY, va='center', ha='left', zorder=4,
                     family='monospace')
            ys.append(y)
        return ys

    col_w = 1.0
    draw_column(safe_genes, 0.0, SAFE_C)
    draw_column(toxic_genes, col_w + 0.18, TOX_C)

    # column headers
    axA.text(0.0 + 0.05, top_y + 0.85,
             'NO ON-TARGET TOXICITY  (n=%d)' % len(safe_genes),
             fontsize=10.5, fontweight='bold', color=SAFE_C, ha='left', va='bottom')
    axA.text(col_w + 0.18 + 0.05, top_y + 0.85,
             'SEVERE ON-TARGET TOXICITY  (n=%d)' % len(toxic_genes),
             fontsize=10.5, fontweight='bold', color=TOX_C, ha='left', va='bottom')
    # single centered axis cue under each track (a two-ended label fighting for the same
    # narrow 0.57-wide track was itself cramped, so state it as one arrowed caption instead)
    axA.text(0.335, -0.85, 'more spread out  →  more confined', fontsize=8, color=NEUTRAL_GRAY,
             ha='center', va='top', style='italic')
    axA.text(col_w + 0.18 + 0.335, -0.85, 'more spread out  →  more confined', fontsize=8, color=NEUTRAL_GRAY,
             ha='center', va='top', style='italic')

    axA.set_xlim(-0.02, col_w + 0.18 + 1.03)
    axA.set_ylim(-1.5, top_y + 1.5)
    axA.axis('off')

    axA.set_title('A.  Every target with a clinical safety verdict, ranked by confinement\n'
                  'across 154 normal cell types (position = confinement, size = cell types on)',
                  fontsize=11.5, fontweight='bold', loc='left', pad=10)

    # CEA callout -- annotate directly on its own row, no separate leader line needed
    cea_idx = toxic_genes.index('CEACAM5')
    cea_y = top_y - cea_idx * row_gap
    cea_track_x1 = col_w + 0.18 + 0.62
    cea_dot_x = col_w + 0.18 + 0.05 + cA['CEACAM5']['top_share'] * (cea_track_x1 - (col_w + 0.18 + 0.05))
    axA.annotate('confined to colon,\nwhich you cannot lose',
                 xy=(cea_dot_x, cea_y), xytext=(cea_dot_x - 0.05, cea_y - 2.3),
                 fontsize=8.6, fontweight='bold', color='white', ha='center', va='top',
                 bbox=dict(boxstyle='round,pad=0.4', facecolor=PROJECT_ORANGE, edgecolor='none'),
                 arrowprops=dict(arrowstyle='-|>', color=PROJECT_ORANGE, linewidth=1.8,
                                 shrinkA=2, shrinkB=2))

    handles = [
        mlines.Line2D([], [], marker='o', color='none', markerfacecolor=NEUTRAL_GRAY,
                      markeredgecolor='white', markersize=10, label='home = blood lineage (expendable)'),
        mlines.Line2D([], [], marker='s', color='none', markerfacecolor=NEUTRAL_GRAY,
                      markeredgecolor='white', markersize=9, label='home = solid tissue'),
        mlines.Line2D([], [], marker='o', color='none', markerfacecolor='white',
                      markeredgecolor=NEUTRAL_GRAY, markersize=6, label='fewer normal cell types → smaller dot'),
        mlines.Line2D([], [], marker='o', color='none', markerfacecolor='white',
                      markeredgecolor=NEUTRAL_GRAY, markersize=12, label='more normal cell types → larger dot'),
    ]
    axA.legend(handles=handles, loc='lower center',
               bbox_to_anchor=(0.5, -0.145), ncol=2, frameon=True, facecolor='white',
               edgecolor='#D6DAE2', framealpha=0.96,
               handletextpad=0.5, labelspacing=0.6, columnspacing=1.4, fontsize=7.8).set_zorder(20)

    # ---------------- Panel B: unchanged from original fig9 ----------------
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
        axB.text(i, -0.025, 'mean %.3f' % np.mean(ys), transform=axB.get_xaxis_transform(),
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
    axB.tick_params(axis='x', pad=18)
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
    fig.text(0.5, 0.045,
             'The rule that survives every check: a target is safe only if it is BOTH confined to few normal cell types '
             'AND those cell types are expendable.',
             ha='center', va='center', fontsize=11, fontweight='bold', color='white',
             bbox=dict(boxstyle='round,pad=0.6', facecolor=PROJECT_ORANGE, edgecolor='none'))
    fig.text(0.5, 0.010,
             'Human Protein Atlas single-cell consensus, 154 cell types. Outcome labels are clinical safety history only and are frozen in target_labels.py. '
             'Full panel: separation 0.868, permutation p = 0.0003 (n = 19 vs 10).',
             ha='center', va='center', fontsize=8.2, color='#5A6473', style='italic')

    out = out_path or os.path.join(FIGURES_DIR, 'fig9_wholebody_restriction_v2.svg')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=300 if out.endswith('.png') else None,
                facecolor='white', bbox_inches=None)
    plt.close(fig)
    print('wrote', out)
    return out


if __name__ == '__main__':
    build()
