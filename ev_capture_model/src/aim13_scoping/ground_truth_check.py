import json, os, numpy as np
d = json.load(open(os.path.expanduser('~/scratch/pilot_roc.json')))
pids = sorted(d); ants=[k for k in d[pids[0]] if isinstance(d[pids[0]][k],dict)]
med_auc={a: float(np.median([d[p][a]['auc'] for p in pids])) for a in ants}
med_mean={a: float(np.median([d[p][a]['mean_B'] for p in pids])) for a in ants}

# GROUND TRUTH, fixed BEFORE looking at ranks (stated in comments):
# Clinically validated B-NHL surface immunotherapy targets = an approved drug
# (CAR-T, bispecific, mAb, or ADC) directed at that antigen in B-cell lymphoma.
VALIDATED = ['CD19','CD20','CD22','CD79B']   # axi/tisa/liso-cel & blina; rituximab/mosun/glofit/epcorit; inotuzumab & CD22 CAR-T; polatuzumab vedotin
# Antigens known to be broadly expressed on normal leukocytes -> known/expected on-target off-tumor problem
BROAD     = ['CD44','HLA-ABC','CD48','CD27','CD43','CD38','CD84','HLA-E','CD25','FasR','CD45RA']

def ranks(score):
    order = sorted(ants, key=lambda a: -score[a])
    return {a:i+1 for i,a in enumerate(order)}
ra, rm = ranks(med_auc), ranks(med_mean)
N=len(ants)
print('N antigens =',N)
print('%-9s %10s %10s %9s %9s' % ('ANTIGEN','rank_AUC','rank_MEAN','medAUC','medMean'))
for a in VALIDATED: print('%-9s %10d %10d %9.3f %9.1f   <- clinically validated target' % (a,ra[a],rm[a],med_auc[a],med_mean[a]))
print()
for a in BROAD: print('%-9s %10d %10d %9.3f %9.1f   <- broadly expressed / known off-tumor risk' % (a,ra[a],rm[a],med_auc[a],med_mean[a]))

def mw(group, r):
    g=[r[a] for a in group]; o=[r[a] for a in ants if a not in group]
    # P(group ranks better i.e. smaller) 
    wins=sum(1 for x in g for y in o if x<y)+0.5*sum(1 for x in g for y in o if x==y)
    return wins/(len(g)*len(o)), float(np.mean(g))
for name,grp in [('VALIDATED (want near top)',VALIDATED),('BROAD (want near bottom)',BROAD)]:
    pa,ma=mw(grp,ra); pm,mm=mw(grp,rm)
    print('\n%s' % name)
    print('  by AUC  : mean rank %5.1f/%d   P(beats a random other antigen)=%.3f' % (ma,N,pa))
    print('  by MEAN : mean rank %5.1f/%d   P(beats a random other antigen)=%.3f' % (mm,N,pm))

# direct separation test: do the 4 validated targets rank above the 11 broad ones?
def sep(r):
    v=[r[a] for a in VALIDATED]; b=[r[a] for a in BROAD]
    w=sum(1 for x in v for y in b if x<y)
    return w/(len(v)*len(b))
print('\nSeparation of VALIDATED above BROAD:')
print('  overlap ranking (AUC)      : %.3f  (1.0 = perfect)' % sep(ra))
print('  mean-expression ranking    : %.3f' % sep(rm))
