from __future__ import annotations
import csv, hashlib, json, os, sys
from pathlib import Path

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parent
os.environ.setdefault('MPLCONFIGDIR', str(STAGE / 'dependencies' / 'matplotlib-config'))
sys.path.insert(0, str(ROOT / 'strengthening_execution_Oct2' / 'dependencies'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

source = STAGE / 'analysis' / 'results' / 'figure_joint_category_deltas.tsv'
with source.open(encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
if len(rows) != 36:
    raise RuntimeError(f'Expected 36 predeclared participant×pair×category values, got {len(rows)}')
patients = ['PM01','PM07','PM10']
pairs = [('TNFRSF17','GPRC5D'),('TNFRSF17','SLAMF7'),('GPRC5D','SLAMF7')]
categories = ['both','A_only','B_only','neither']
colors = {'PM01':'#0072B2','PM07':'#D55E00','PM10':'#009E73'}
markers = {'PM01':'o','PM07':'s','PM10':'^'}
lookup = {(r['participant'],r['marker_A'],r['marker_B'],r['category']): float(r['post_minus_pre_fraction']) for r in rows}
max_abs = max(abs(v) for v in lookup.values())
lim = max(0.025, max_abs * 1.15)
fig, axes = plt.subplots(1,3,figsize=(12.5,4.7),sharey=True)
x = range(4)
offsets = {'PM01':-0.12,'PM07':0.0,'PM10':0.12}
for ax,(a,b) in zip(axes,pairs):
    for patient in patients:
        ys=[lookup[(patient,a,b,cat)] for cat in categories]
        ax.scatter([i+offsets[patient] for i in x], ys, s=54, marker=markers[patient],
                   color=colors[patient], edgecolors='white', linewidths=0.5, label=patient, zorder=3)
    ax.axhline(0,color='#555555',linewidth=0.9,zorder=1)
    ax.set_title(f'{a} / {b}',fontsize=11,pad=9)
    ax.set_xticks(list(x),['Both','A only','B only','Neither'])
    ax.set_xlim(-0.4,3.4)
    ax.set_ylim(-lim,lim)
    ax.grid(axis='y',color='#D9DEE5',linewidth=0.65,zorder=0)
    ax.tick_params(axis='x',labelrotation=0,labelsize=8)
    ax.tick_params(axis='y',labelsize=8)
axes[0].set_ylabel('Post − pre change in category fraction')
fig.suptitle('Paired CD138-selected RNA detection categories\nRaw deposited counts; detection = count > 0; not depth-adjusted',fontsize=13,y=1.04)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles[:3],labels[:3],loc='upper center',bbox_to_anchor=(0.5,0.91),ncol=3,frameon=False,title='Participant')
fig.text(0.5,0.015,'All filtered-matrix barcodes retained; categories use fixed target order. No per-cell malignant labels, protein readout, or lineage tracking.',ha='center',fontsize=8,color='#444444')
fig.tight_layout(rect=(0,0.07,1,0.80))
outdir=STAGE/'analysis'/'results'
svg=outdir/'paired_joint_category_change.svg'
png=outdir/'paired_joint_category_change.png'
fig.savefig(svg,format='svg',bbox_inches='tight',metadata={'Title':'Paired joint-category fractions before and after venetoclax','Creator':'plot_pilot_figure.py'})
fig.savefig(png,format='png',dpi=300,bbox_inches='tight',metadata={'Title':'Paired joint-category fractions before and after venetoclax'})
plt.close(fig)
for p in (svg,png):
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    print(f'{p.relative_to(STAGE)}\t{p.stat().st_size}\t{h}')

# Companion view: raw versus the mean of the locked 100 replicate joint-thinning changes.
results = STAGE / 'analysis' / 'results'
depth_D = json.loads((results / 'pilot_run_manifest.json').read_text(encoding='utf-8'))['depth_sensitivity']['D']
with (results / 'participant_paired_deltas.tsv').open(encoding='utf-8-sig', newline='') as f:
    raw_rows = list(csv.DictReader(f, delimiter='\t'))
raw = {(r['participant'],r['marker_A'],r['marker_B']): float(r['post_minus_pre'])
       for r in raw_rows if r['metric'] == 'n00_fraction'}
replicate_values = {}
with (results / 'depth_sensitivity_paired_deltas.tsv').open(encoding='utf-8-sig', newline='') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['metric'] == 'n00_fraction':
            key=(r['participant'],r['marker_A'],r['marker_B'])
            replicate_values.setdefault(key,[]).append(float(r['post_minus_pre']))
if set(raw) != set(replicate_values) or any(len(v) != 100 for v in replicate_values.values()):
    raise RuntimeError('Expected 9 paired rows with 100 joint-thinning replicates each')
d_rows=[]
for pair in pairs:
    for patient in patients:
        key=(patient,pair[0],pair[1])
        vals=replicate_values[key]
        d_rows.append({
            'participant':patient,'marker_A':pair[0],'marker_B':pair[1],
            'raw_post_minus_pre_double_undetected':raw[key],
            'mean_100_joint_depth_sensitivity_replicates':sum(vals)/len(vals),
            'replicate_count':len(vals),'D':depth_D,
            'interpretation':'within-specimen thinning sensitivity mean; not participant uncertainty'
        })
values_path=results/'figure_depth_sensitivity_values.tsv'
with values_path.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(d_rows[0]),delimiter='\t',lineterminator='\n')
    w.writeheader();w.writerows(d_rows)
fig,axes=plt.subplots(1,3,figsize=(11.5,4.6),sharey=True)
shape={'raw':'o','depth':'s'}
for ax,pair in zip(axes,pairs):
    for i,patient in enumerate(patients):
        r=next(r for r in d_rows if r['participant']==patient and r['marker_A']==pair[0] and r['marker_B']==pair[1])
        raw_y=r['raw_post_minus_pre_double_undetected']
        depth_y=r['mean_100_joint_depth_sensitivity_replicates']
        ax.plot([i-0.13,i+0.13],[raw_y,depth_y],color=colors[patient],linewidth=1.0,alpha=.7,zorder=2)
        ax.scatter(i-0.13,raw_y,s=56,marker=shape['raw'],color=colors[patient],edgecolors='white',linewidths=.5,zorder=3)
        ax.scatter(i+0.13,depth_y,s=56,marker=shape['depth'],color=colors[patient],edgecolors='white',linewidths=.5,zorder=3)
    ax.axhline(0,color='#555555',linewidth=.9,zorder=1)
    ax.set_title(f'{pair[0]} / {pair[1]}',fontsize=10,pad=8)
    ax.set_xticks(range(3),patients)
    ax.set_xlim(-.5,2.5)
    ax.grid(axis='y',color='#D9DEE5',linewidth=.65,zorder=0)
    ax.tick_params(axis='x',labelsize=8)
    ax.tick_params(axis='y',labelsize=8)
axes[0].set_ylabel('Post − pre fraction with neither transcript detected')
fig.suptitle(f'Raw and depth-sensitivity changes in double-undetected fractions\nDepth sensitivity = mean of 100 joint-thinning replicates; D={depth_D:,}',fontsize=12,y=1.05)
from matplotlib.lines import Line2D
legend=[Line2D([0],[0],marker='o',color='#555555',linestyle='None',label='Raw deposited counts'),
        Line2D([0],[0],marker='s',color='#555555',linestyle='None',label='Mean joint-thinning sensitivity')]
fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.91),ncol=2,frameon=False)
fig.text(.5,.015,'100 deterministic replicate draws per specimen. These are sampling-sensitivity summaries, not confidence intervals or participant uncertainty.',ha='center',fontsize=8,color='#444444')
fig.tight_layout(rect=(0,.07,1,.79))
svg=results/'paired_double_undetected_raw_vs_depth_sensitivity.svg'
png=results/'paired_double_undetected_raw_vs_depth_sensitivity.png'
fig.savefig(svg,format='svg',bbox_inches='tight',metadata={'Title':'Raw and depth-sensitivity change in double-undetected fraction','Creator':'plot_pilot_figure.py'})
fig.savefig(png,format='png',dpi=300,bbox_inches='tight',metadata={'Title':'Raw and depth-sensitivity change in double-undetected fraction'})
plt.close(fig)
for p in (svg,png):
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    print(f'{p.relative_to(STAGE)}\t{p.stat().st_size}\t{h}')
