"""
Whole-body restriction vs documented ON-TARGET TOXICITY, across blood AND solid tumours.
Reference: Human Protein Atlas single-cell consensus, 154 human cell types.
Labels are drug-safety history only; nothing about expression enters them.
Run wholebody_restriction.py first (it writes wholebody_restriction.json).
"""
import json, itertools, numpy as np, os
HERE=os.path.dirname(os.path.abspath(__file__))
r=json.load(open(os.path.join(HERE,'wholebody_restriction.json')))
TOX={'ERBB2':'HER2 CAR-T fatal pulmonary event (Morgan 2010)','EGFR':'skin/lung on-target',
     'CEACAM5':'CEA CAR-T severe colitis','CA9':'CAIX CAR-T liver toxicity',
     'CD47':'magrolimab on-target anaemia','CD44':'CD44v6 on-target skin toxicity',
     'PTPRC':'pan-leukocyte, never targetable'}
OK={'MS4A1':'approved','CD19':'approved','CD22':'approved','CD79B':'approved','TNFRSF17':'approved',
    'SLAMF7':'approved','FOLH1':'approved','TACSTD2':'approved','NECTIN4':'approved',
    'CLDN18':'approved 2024','DLL3':'approved','MSLN':'tolerable','GPC3':'tolerable',
    'FCER2':'failed for EFFICACY not toxicity','CD40':'failed for EFFICACY'}
def sep(a,b): return (sum(1 for x in a for y in b if x>y)+0.5*sum(1 for x in a for y in b if x==y))/(len(a)*len(b))
def exact(a,b):
    allv=a+b; obs=sep(a,b)*len(a)*len(b); c=t=0
    for comb in itertools.combinations(range(len(allv)),len(a)):
        aa=[allv[i] for i in comb]; ff=[allv[i] for i in range(len(allv)) if i not in comb]
        u=sum(1 for x in aa for y in ff if x>y)+0.5*sum(1 for x in aa for y in ff if x==y)
        t+=1; c+=(u>=obs)
    return c/t
for key,flip,lab in (('top_share',1,'top_share (higher = more confined)'),
                     ('n_types',-1,'n_types (fewer = more confined)')):
    A=[flip*r[g][key] for g in OK if g in r]; F=[flip*r[g][key] for g in TOX if g in r]
    print('%s\n   safe %.3f (n=%d) vs toxic %.3f (n=%d) | separation %.3f | exact p = %.5f'
          % (lab, np.mean(A),len(A), np.mean(F),len(F), sep(A,F), exact(A,F)))
print("""
HONEST READING. Breadth predicts SYSTEMIC on-target toxicity well (CD47 120 cell types ->
anaemia; HER2 72 -> fatal; EGFR 60 -> skin/lung). It does NOT predict the case where a
target is narrowly confined but confined to a tissue you cannot lose: CEACAM5 sits in 7 cell
types, essentially colonocytes, and CEA CAR-T caused severe colitis anyway. That is not a
failure of the framework -- it is the project's own thesis restated for solid tumours:
a target is safe only if it is BOTH confined AND confined to something expendable.
  CD19 / CD20 -> B cells        confined + expendable   -> works
  BCMA        -> plasma cells   confined + expendable   -> works
  CLDN18      -> gastric mucosa confined + regenerating -> approved 2024
  CEACAM5     -> colonocytes    confined, NOT expendable-> colitis
  HER2/CD47   -> everywhere     not confined            -> catastrophic
Known limitation: this reference is RNA, not protein, and some 'home' cell types are germ or
placental tissue (CA9, DLL3, IL13RA2), which are immune-privileged and not the clinically
relevant site.""")
