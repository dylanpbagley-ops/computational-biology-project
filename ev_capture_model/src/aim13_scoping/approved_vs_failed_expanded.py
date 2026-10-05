"""
The outcome test, expanded and generalised.

GENERALISATION: the predictor is no longer "B-lineage restriction" but "restriction to the
lineage the cancer came from" -- B cells for B-cell malignancies, myeloid cells for AML.
That is the same idea (is the collateral confined to one, expendable lineage?) stated in a
form that covers more than one disease, and it roughly doubles the number of labelled
targets.

INCLUSION RULE, fixed before computing anything: a target is included only if
  (a) a real drug against it was given to patients, and its outcome (approved / never
      approved) is documented, and
  (b) the malignancy it was aimed at has a clear normal counterpart lineage present in
      peripheral blood.
(b) excludes myeloma targets (plasma cells are not in blood -- already shown to make them
unscoreable here) and Hodgkin targets such as CD30 and CD25, whose Reed-Sternberg origin
has no clean PBMC counterpart. Excluding them is a judgement made on the disease, not on
any number.

TWO PREDICTORS, both fixed in advance:
  primary   restriction  = AUC(disease lineage vs every other healthy cell)
  secondary margin       = restriction MINUS the highest restriction this antigen achieves
                           for any OTHER lineage. This penalises antigens that are also high
                           somewhere else, which is the HLA-DR failure mode noted earlier.
                           It is specified here BEFORE being computed, because inventing it
                           after seeing that it fixes HLA-DR would be fitting to the answer.
"""
import os, json, math, itertools, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','healthy_pbmc_citeseq'))
rng=np.random.default_rng(313); SUB=8000
z=np.load(os.path.join(D,'healthy_adt.npz'), allow_pickle=True)
NAMES=[str(x) for x in z['names']]; L1=np.array([str(x) for x in z['celltype_l1']])
X=np.log(z['adt'].astype(np.float32)+1.0); X=X-X.mean(axis=1,keepdims=True)
def auc(x,y):
    n1,n2=len(x),len(y); zz=np.concatenate([x,y]); o=np.argsort(zz,kind='mergesort')
    r=np.empty(len(zz)); r[o]=np.arange(1,len(zz)+1)
    return (r[:n1].sum()-n1*(n1+1)/2)/(n1*n2)
def sub(i): return rng.choice(i,SUB,replace=False) if len(i)>SUB else i
LIN=[l for l in sorted(set(L1)) if l!='other']
idx={l:np.flatnonzero(L1==l) for l in LIN}
oth={l:np.flatnonzero((L1!=l)&(L1!='other')) for l in LIN}
S={}
for k,n in enumerate(NAMES):
    S[n]={l: float(auc(X[sub(idx[l]),k], X[sub(oth[l]),k])) for l in LIN}

# antigen -> (healthy-panel name, disease lineage, outcome, drug + evidence)
T=[
 ('CD19','B','APPROVED','axi-cel/tisa-cel/liso-cel, blinatumomab'),
 ('CD20','B','APPROVED','rituximab, mosunetuzumab, glofitamab, epcoritamab'),
 ('CD22','B','APPROVED','inotuzumab ozogamicin'),
 ('CD79b','B','APPROVED','polatuzumab vedotin'),
 ('CD52','B','APPROVED','alemtuzumab (CLL)'),
 ('CD23','B','FAILED','lumiliximab - LUCID ph2/3 n=627 stopped early, doi 10.1111/bjh.13061'),
 ('CD80','B','FAILED','galiximab - CALGB 50402, never approved'),
 ('CD40','B','FAILED','dacetuzumab/SGN-40, never approved'),
 ('HLA-DR','B','FAILED','apolizumab/Hu1D10, never approved'),
 ('CD200','B','FAILED','samalizumab, never approved'),
 ('CD33','Mono','APPROVED','gemtuzumab ozogamicin (AML)'),
 ('CD123','Mono','FAILED','talacotuzumab (AML), PMID 32203138'),
]
# SENSITIVITY ADDITIONS, flagged because they are UNBALANCED: two further failures and no
# further approvals. That is not cherry-picking -- the 228-antibody panel simply contains no
# other approved blood-cancer surface target meeting the inclusion rule (CD23 and CD33 are
# absent from the panel; CD30/CD25 are Hodgkin; CD38/BCMA/SLAMF7 are myeloma). Adding only
# failures can only move the separation in one direction, so this is reported SEPARATELY
# from the primary result and never as the headline.
T_SENS=[
 ('CD47','Mono','FAILED','magrolimab - ENHANCE phase 3 programme discontinued (PMIDs 36888930, 40009500)'),
 ('CD70','Mono','FAILED','cusatuzumab (AML), never approved'),
]
rows=[]
print('%-8s %-6s %-9s %11s %9s  %s' % ('ANTIGEN','lineage','outcome','restriction','margin','drug'))
for n,lin,out,drug in T:
    if n not in S: print('%-8s  NOT IN PANEL' % n); continue
    r=S[n][lin]; other=max(v for l,v in S[n].items() if l!=lin)
    rows.append((out,n,lin,r,r-other))
    print('%-8s %-6s %-9s %11.3f %9.3f  %s' % (n,lin,out,r,r-other,drug))

def sep(A,F):
    return (sum(1 for a in A for f in F if a>f)+0.5*sum(1 for a in A for f in F if a==f))/(len(A)*len(F))
def exact_p(A,F):
    allv=A+F; obs=sep(A,F)*len(A)*len(F); c=t=0
    for comb in itertools.combinations(range(len(allv)),len(A)):
        aa=[allv[i] for i in comb]; ff=[allv[i] for i in range(len(allv)) if i not in comb]
        u=sum(1 for a in aa for f in ff if a>f)+0.5*sum(1 for a in aa for f in ff if a==f)
        t+=1; c+= (u>=obs)
    return c/t
print()
for lab,ix in (('PRIMARY: restriction to the disease lineage',3),
               ('SECONDARY (pre-specified): margin over the next-best lineage',4)):
    A=[r[ix] for r in rows if r[0]=='APPROVED']; F=[r[ix] for r in rows if r[0]=='FAILED']
    print('%s' % lab)
    print('   approved mean %.3f (n=%d) | failed mean %.3f (n=%d) | separation %.3f | exact p = %.4f'
          % (np.mean(A),len(A),np.mean(F),len(F),sep(A,F),exact_p(A,F)))
    mis=[(r[1],r[ix]) for r in rows if r[0]=='FAILED' and r[ix]>min(a for a in A)]
    print('   failed targets scoring above the weakest approved one: %s' % (mis if mis else 'none'))
print()
print('='*72)
print('SENSITIVITY ANALYSIS -- unbalanced (adds 2 failures, 0 approvals; see note in code)')
rows2=list(rows)
for n,lin,out,drug in T_SENS:
    if n not in S: print('%-8s NOT IN PANEL' % n); continue
    r=S[n][lin]; other=max(v for l,v in S[n].items() if l!=lin)
    rows2.append((out,n,lin,r,r-other))
    print('%-8s %-6s %-9s %11.3f %9.3f  %s' % (n,lin,out,r,r-other,drug))
for lab,ix in (('restriction',3),('margin',4)):
    A=[r[ix] for r in rows2 if r[0]=='APPROVED']; F=[r[ix] for r in rows2 if r[0]=='FAILED']
    print('   %-12s approved %.3f (n=%d) vs failed %.3f (n=%d) | separation %.3f | exact p = %.4f'
          % (lab,np.mean(A),len(A),np.mean(F),len(F),sep(A,F),exact_p(A,F)))
print()
print('CEILING ON THIS TEST: the number of APPROVED surface targets in blood cancer is small')
print('in the world, not just in this panel. n cannot be raised much by more work; this test')
print('will stay underpowered no matter how carefully it is done, and should be reported as')
print('suggestive evidence rather than a significant result.')
json.dump({'rows':[list(r) for r in rows],'rows_sensitivity':[list(r) for r in rows2]},
          open(os.path.join(HERE,'approved_vs_failed_expanded.json'),'w'), indent=1)
