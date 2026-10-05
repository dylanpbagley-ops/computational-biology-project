"""Scientific figures from saved outputs; no scoring or source mutation."""
import pathlib,sys,json
HERE=pathlib.Path(__file__).resolve().parent; OUT=HERE.parent
sys.path.insert(0,str(OUT/'dependencies'));sys.path.insert(0,str(OUT))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from startup import read
from run_b_pilot import load,GENES
from thresholds import fraction_above

F=OUT/'delivery/figures'
COLORS={'pooled':'#41658a','equal_donor':'#6c757d','donor_envelope':'#b26342'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})

def save(fig,name):
    F.mkdir(parents=True,exist_ok=True)
    fig.savefig(F/(name+'.png'),dpi=300,bbox_inches='tight')
    fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig)

def main():
    groups=load();summary=read(HERE/'results/b_summary.tsv')
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for j,g in enumerate(GENES):
        ax=axes[0,j]; names=list(groups['NBM'])
        for policy in ['pooled','donor_envelope']:
            r=next(r for r in summary if r['gene']==g and r['policy']==policy and float(r['budget'])==.1)
            t=float(r['threshold']);v=[100*fraction_above(groups['NBM'][n][g],t) for n in names]
            ax.plot(range(len(names)),v,'o-',color=COLORS[policy],label=policy.replace('_',' '))
        ax.axhline(10,color='#333',linestyle='--',linewidth=1);ax.set_title(g);ax.set_ylim(0,32)
        ax.set_xticks(range(len(names)),names,rotation=55,ha='right');ax.set_ylabel('Measured normal cells above threshold (%)')
        ax=axes[1,j]
        for k,policy in enumerate(['pooled','donor_envelope']):
            r=next(r for r in summary if r['gene']==g and r['policy']==policy and float(r['budget'])==.1)
            val=100*float(r['tumor_pooled']);ax.bar(k,val,color=COLORS[policy]);ax.text(k,val+.5,f'{val:.1f}%',ha='center')
        ax.set_xticks([0,1],['Pooled policy','Donor envelope']);ax.set_ylim(0,32);ax.set_ylabel('Pooled measured MM cells above threshold (%)')
    axes[0,0].legend(frameon=False,loc='upper left',fontsize=8)
    fig.suptitle('A pooled RNA constraint does not protect every measured donor\n9 healthy and 8 myeloma samples; illustrative 10% normal-cell fraction budget',fontsize=13)
    save(fig,'b_observed_donor_tradeoff')
    rows=read(HERE/'results/b_held_out_donors.tsv')
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for j,g in enumerate(GENES):
        for i,stage in enumerate(['NBM','MM']):
            ax=axes[i,j]
            for k,policy in enumerate(['pooled','equal_donor','donor_envelope']):
                rs=[r for r in rows if r['gene']==g and r['policy']==policy and r['stage']==stage]
                y=[100*float(r['mean_held_out_fraction']) for r in rs]
                ax.scatter(np.repeat(k,len(y))+np.linspace(-.10,.10,len(y)),y,color=COLORS[policy],s=28)
            if stage=='NBM':ax.axhline(10,color='#333',linestyle='--',linewidth=1)
            ax.set_ylim(0,62 if stage=='MM' else 40)
            ax.set_xticks([0,1,2],['Pooled','Equal donor','Envelope'],rotation=20)
            ax.set_ylabel(('Held-out normal' if stage=='NBM' else 'Held-out MM')+' fraction (%)')
            if i==0:ax.set_title(g)
    fig.suptitle('Protection of observed donors does not guarantee a withheld donor\nEach dot is one donor, averaged across partner folds; 72 folds are dependent',fontsize=13)
    save(fig,'b_held_out_donors')
    rows=read(OUT/'feasibility/candidate_a/comparability.tsv');fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    y=np.arange(len(rows));reg=np.array([int(r['regulatory_rows']) for r in rows]);allr=np.array([int(r['source_rows']) for r in rows])
    ax.barh(y,allr-reg,label='Publication leads / missing-report lead',color='#41658a');ax.barh(y,reg,left=allr-reg,label='Regulatory descriptions',color='#b3b9c0')
    ax.set_yticks(y,[r['gene'] for r in rows]);ax.set_xlabel('Source rows — not independent outcomes');ax.legend(frameon=False)
    ax.set_title('Clinical feasibility workload: 35 rows, 12 target strings\nNo independently adjudicated comparable binary pair established')
    save(fig,'diagnostic_a_source_types')
    old=json.loads((HERE/'results/historical_reproduction.json').read_text());fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
    labels=['Frozen rule','Specificity / tau','Blood flag'];v=[old['metrics'][k] for k in ['frozen_rule','top_share','blood_flag']]
    ax.bar(labels,v,color=['#41658a','#809bb1','#b3b9c0']);ax.set_ylim(0,1);ax.axhline(.5,color='#333',linestyle='--');ax.set_ylabel('AUC on historical labels')
    for k,z in enumerate(v):ax.text(k,z+.02,f'{z:.3f}',ha='center')
    ax.set_title('Historical selected panel reproduced\n19 favorable / 11 unfavorable; labels remain unadjudicated')
    save(fig,'historical_baselines')
    print('Saved four figures as PNG (300 dpi) and SVG; matplotlib '+matplotlib.__version__)

if __name__=='__main__':main()
