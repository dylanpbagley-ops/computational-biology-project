"""
fig10_pancreatic_thread.py

Write-up figure 10 -- the pancreatic thread, end to end, in one figure.

A  the real CLDN18.2 dose-response calibration on pancreatic and gastric lines, with the
   project's own binding model fitted to it (Aim 10).
B  CPTAC pancreatic ductal adenocarcinoma, n = 135: mesothelin protein level vs survival.
   This panel is an HONEST NEGATIVE -- the curves do not separate (log-rank p = 0.35) -- and
   it is included precisely because a figure set that only shows wins is not evidence.
C  the whole-body normal atlas ranking of 14 real PDAC surface targets, with the clinical
   verdict shown as colour. The atlas contains no tumour data and no outcome data.

All numbers read live from pancreatic_thread.json.
"""
import os, sys, json, pickle
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from figure_style import clean_spines, section_badge, PROJECT_BLUE, PROJECT_ORANGE, INK, NEUTRAL_GRAY

SRC = os.path.abspath(os.path.join(HERE, '..'))
RES = os.path.join(SRC, 'aim13_scoping', 'pancreatic_thread.json')
FIGURES_DIR = os.path.join(HERE, "..", "..", "..", "Poster, Paper, and Write Up", "Figures")

SAFE_C = "#2E5EAA"
TOX_C = "#B0563A"
NEU_C = "#8A93A6"

# clinical verdict per PDAC target, for panel C colour only (never used in the ranking)
VERDICT = {'claudin-18.2': 'approved', 'mesothelin': 'tolerable', 'CEA': 'on-target toxicity',
           'MUC1': 'failed', 'EGFR': 'on-target toxicity', 'HER2': 'on-target toxicity',
           'PSCA': 'in trials', 'TROP2': 'approved elsewhere', 'glypican-1': 'diagnostic marker',
           'CEACAM6': 'preclinical', 'CA125/MUC16': 'failed', 'B7-H3': 'in trials',
           'integrin b6': 'early', 'CDH17': 'early'}
VCOLOR = {'approved': SAFE_C, 'approved elsewhere': SAFE_C, 'tolerable': SAFE_C,
          'on-target toxicity': TOX_C, 'failed': TOX_C,
          'in trials': NEU_C, 'preclinical': NEU_C, 'early': NEU_C, 'diagnostic marker': NEU_C}


def kaplan_meier(t, e):
    order = np.argsort(t); t, e = t[order], e[order]
    times, surv, n = [0.0], [1.0], len(t)
    s = 1.0
    for i, ti in enumerate(np.unique(t)):
        at_risk = (t >= ti).sum()
        d = ((t == ti) & (e == 1)).sum()
        if d and at_risk:
            s *= (1 - d / at_risk)
            times.append(ti); surv.append(s)
    return np.array(times), np.array(surv)


def build(out_path=None):
    R = json.load(open(RES))
    c = pickle.load(open(os.path.join(SRC, 'pancreatic_calibration', 'cldn182_fit_results.pkl'), 'rb'))
    rows = json.load(open(os.path.join(SRC, 'cptac_msln_outcomes', 'cptac_msln_merged.json')))

    fig = plt.figure(figsize=(15.0, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.28], wspace=0.32,
                          left=0.055, right=0.985, top=0.84, bottom=0.235)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))

    # ---------- A: CLDN18.2 calibration ----------
    data = c['data']
    abc = np.array([d[1] for d in data], float)
    cyt = np.array([d[2] for d in data], float)
    panc = np.array([d[3] for d in data], bool)
    k, nh = c['k_intra_fit'], c['n_half_fit']
    grid = np.logspace(2, np.log10(6e5), 300)
    # the project's saturating occupancy form, normalised to the fitted plateau
    pred = 100 * (1 - np.exp(-k * grid)) / (1 - np.exp(-k * 6e5))
    plateau = cyt[abc > 0].max()
    pred = pred / pred.max() * plateau
    axA.plot(grid, pred, color=INK, linewidth=2.0, zorder=2, label='fitted binding model')
    axA.scatter(np.clip(abc[panc], 300, None), cyt[panc], s=150, marker='o',
                facecolor=PROJECT_ORANGE, edgecolor='white', linewidth=1.6, zorder=4,
                label='pancreatic line')
    axA.scatter(abc[~panc], cyt[~panc], s=120, marker='s', facecolor=NEU_C,
                edgecolor='white', linewidth=1.5, zorder=3, label='gastric line')
    axA.set_xscale('log')
    axA.set_xlim(250, 7e5)
    axA.set_ylim(-4, 100)
    axA.set_xlabel('CLDN18.2 molecules per cell (antibody binding capacity)', fontsize=10)
    axA.set_ylabel('% tumour cells killed at 48 h', fontsize=10)
    axA.set_title('A.  A real pancreatic dose-response,\nfitted by this project\'s model',
                  fontsize=12, fontweight='bold', loc='left', pad=9)
    axA.legend(fontsize=8.6, loc='center right', bbox_to_anchor=(1.0, 0.42), frameon=False)
    clean_spines(axA); axA.grid(color='#E6E9EF', linewidth=0.8, zorder=0); axA.set_axisbelow(True)
    axA.text(0.03, 0.955, 'Spearman r = %.3f\np = %.4f     RMSE %.2f'
             % (c['spearman_r'], c['spearman_p'], c['rmse']),
             transform=axA.transAxes, va='top', fontsize=9.5, fontweight='bold', color=INK,
             bbox=dict(boxstyle='round,pad=0.45', facecolor='#EEF1F6', edgecolor='none'))
    axA.annotate('the KO line (zero antigen, no killing)\nis plotted at the axis floor',
                 xy=(320, cyt[abc == 0][0]), xytext=(0.33, 0.10), textcoords='axes fraction',
                 fontsize=8.2, color='#5A6473',
                 arrowprops=dict(arrowstyle='-', color='#9AA3B2', linewidth=0.8))

    # ---------- B: CPTAC survival ----------
    x = np.array([r['msln_protein_log2'] for r in rows], float)
    t = np.array([r['follow_up_days'] for r in rows], float)
    e = np.array([r['event'] for r in rows], int)
    med = np.median(x)
    for grp, col, lab in ((x > med, PROJECT_ORANGE, 'high mesothelin'),
                          (x <= med, PROJECT_BLUE, 'low mesothelin')):
        tt, ss = kaplan_meier(t[grp], e[grp])
        axB.step(tt, ss, where='post', color=col, linewidth=2.4,
                 label='%s  (n=%d)' % (lab, grp.sum()))
    axB.set_xlabel('days of follow-up', fontsize=10)
    axB.set_ylabel('fraction surviving', fontsize=10)
    axB.set_ylim(0, 1.02)
    axB.set_title('B.  CPTAC pancreatic cancer, n = 135:\nan honest negative result',
                  fontsize=12, fontweight='bold', loc='left', pad=9)
    axB.legend(fontsize=9, loc='upper right', frameon=False)
    clean_spines(axB); axB.grid(color='#E6E9EF', linewidth=0.8, zorder=0); axB.set_axisbelow(True)
    axB.text(0.03, 0.06, 'log-rank p = %.2f — the curves do NOT separate.\n'
                         'Target level does not predict survival here,\nand this is reported as such.'
             % R['cptac']['logrank_p'],
             transform=axB.transAxes, va='bottom', fontsize=9, color=INK,
             bbox=dict(boxstyle='round,pad=0.45', facecolor='#F3EFE6', edgecolor='#DED5C4'))

    # ---------- C: PDAC ranking ----------
    rank = R['ranking']
    names = [r[1] for r in rank][::-1]
    vals = [r[2] for r in rank][::-1]
    cols = [VCOLOR[VERDICT[n]] for n in names]
    ypos = np.arange(len(names))
    axC.barh(ypos, vals, color=cols, edgecolor='white', linewidth=1.0, height=0.72, zorder=3)
    axC.set_yticks(ypos)
    axC.set_yticklabels(names, fontsize=9.2)
    axC.set_xlabel('confinement across 154 normal human cell types', fontsize=10)
    axC.set_xlim(0, 0.46)
    axC.set_title('C.  14 real pancreatic targets, ranked by a normal-tissue\n'
                  'atlas holding no tumour and no outcome data',
                  fontsize=11.5, fontweight='bold', loc='left', pad=9)
    clean_spines(axC, keep=('left', 'bottom'))
    axC.grid(axis='x', color='#E6E9EF', linewidth=0.8, zorder=0); axC.set_axisbelow(True)
    for y, n, v in zip(ypos, names, vals):
        axC.text(v + 0.006, y, VERDICT[n], va='center', fontsize=8.0, color='#5A6473')
    handles = [plt.Line2D([], [], marker='s', color='none', markerfacecolor=SAFE_C,
                          markersize=11, label='approved / tolerable'),
               plt.Line2D([], [], marker='s', color='none', markerfacecolor=TOX_C,
                          markersize=11, label='failed or on-target toxicity'),
               plt.Line2D([], [], marker='s', color='none', markerfacecolor=NEU_C,
                          markersize=11, label='no verdict yet')]
    axC.legend(handles=handles, fontsize=8.6, loc='lower right', frameon=False)

    section_badge(fig, 'CAR-T', 10)
    fig.text(0.5, 0.085,
             'The pancreatic thread runs from the project\'s first result to its last — and the parts that did not work are shown too.',
             ha='center', va='center', fontsize=11.5, fontweight='bold', color='white',
             bbox=dict(boxstyle='round,pad=0.6', facecolor=PROJECT_ORANGE, edgecolor='none'))
    fig.text(0.5, 0.022,
             'A: Barrett et al. 2024 Fig 5D, values read from the published figure. B: CPTAC PDAC proteomics, median split, log-rank. '
             'C: Human Protein Atlas single-cell consensus, 154 cell types. Clinical verdicts are shown for context only and never enter the ranking.',
             ha='center', va='center', fontsize=8.0, color='#5A6473', style='italic')

    out = out_path or os.path.join(FIGURES_DIR, 'fig10_pancreatic_thread.png')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=300 if out.endswith('.png') else None, facecolor='white')
    plt.close(fig)
    print('wrote', out)
    return out


if __name__ == '__main__':
    build()
