"""
THE OUTCOME FLIP, in its strongest form.

Earlier tests used "approved drug" versus "broadly expressed antigen" -- but the second
group was labelled from expression biology, which is close to what the model measures, so
the test is partly circular. This version fixes that: BOTH groups are antigens against which
a real drug was built and taken into patients. The only difference is what happened.

  APPROVED : the drug reached regulatory approval.
  FAILED   : the drug entered clinical trials and never reached approval.

Nothing about expression enters the labels. This is a clinical-outcome test.

Two candidate predictors are compared:
  (a) the tumor-vs-normal CEILING measured in DLBCL (3 compartments), and
  (b) B-LINEAGE RESTRICTION measured in 161,764 cells from 8 HEALTHY donors --
      specifically AUC(B cells vs every other healthy cell), NOT the max over lineages,
      because a target homed to monocytes is useless for a B-cell cancer no matter how
      restricted it is.
"""
import os, sys, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','healthy_pbmc_citeseq'))
rng=np.random.default_rng(131); SUB=8000

z=np.load(os.path.join(D,'healthy_adt.npz'), allow_pickle=True)
M=z['adt'].astype(np.float32); NAMES=[str(x) for x in z['names']]
L1=np.array([str(x) for x in z['celltype_l1']])
X=np.log(M+1.0); X=X-X.mean(axis=1,keepdims=True); del M
def auc(x,y):
    n1,n2=len(x),len(y); zz=np.concatenate([x,y]); o=np.argsort(zz,kind='mergesort'); zs=zz[o]
    r=np.empty(len(zz)); i=0
    while i<len(zz):
        j=i
        while j+1<len(zz) and zs[j+1]==zs[i]: j+=1
        r[i:j+1]=(i+j)/2.0+1; i=j+1
    rr=np.empty(len(zz)); rr[o]=r
    return (rr[:n1].sum()-n1*(n1+1)/2.0)/(n1*n2)
def sub(idx): return rng.choice(idx,SUB,replace=False) if len(idx)>SUB else idx
bi=sub(np.flatnonzero(L1=='B')); oi=sub(np.flatnonzero((L1!='B')&(L1!='other')))
brest={n: float(auc(X[bi,k],X[oi,k])) for k,n in enumerate(NAMES)}
ISO=['Rat-IgG1-1','Rat-IgG1-2','Rat-IgG2b','Rag-IgG2c']
floor=max(abs(brest[n]-0.5) for n in ISO)+0.5
print('isotype-control B-restriction values: %s' % {n:round(brest[n],3) for n in ISO})
print('noise floor (max isotype deviation from 0.5, mirrored): %.3f\n' % floor)

ceil3=json.load(open(os.path.join(HERE,'three_compartment_results.json')))
CEIL={a: ceil3[a]['joint3'] for a in ceil3}

# ---- the labelled set: every B-cell-malignancy surface target in this panel that had a
# real clinical program. Labels are drug outcomes only, assigned from the literature.
APPROVED={
 'CD19':  ('CD19',  'axi-cel / tisa-cel / liso-cel / blinatumomab'),
 'CD20':  ('CD20',  'rituximab / mosunetuzumab / glofitamab / epcoritamab'),
 'CD22':  ('CD22',  'inotuzumab ozogamicin (ADC)'),
 'CD79b': ('CD79B', 'polatuzumab vedotin (ADC)'),
}
FAILED={
 'CD23':  ('CD23',  'lumiliximab - LUCID ph2/3 n=627 STOPPED EARLY for lack of efficacy (Awan 2014)'),
 'CD80':  ('CD80',  'galiximab - CALGB 50402 ph2, never approved'),
 'CD40':  ('CD40',  'dacetuzumab / SGN-40, never approved'),
 'HLA-DR':('HLA-DR','apolizumab / Hu1D10, never approved'),
}
# EXPANSION, added by a fixed rule and NOT after inspecting results: every remaining
# antigen in the healthy panel that had a real clinical antibody programme in a lymphoid
# malignancy. CD52 is included even though it is expected to HURT the separation
# (alemtuzumab is approved for CLL yet is a pan-lymphocyte antigen causing profound
# lymphopenia) -- excluding an inconvenient case would be exactly the error this test exists
# to avoid.
APPROVED_EXT={'CD52': ('CD52','alemtuzumab - approved in CLL (pan-lymphocyte, severe lymphopenia)')}
FAILED_EXT={'CD25':  ('CD25','camidanlumab tesirine - never approved'),
            'CD200': ('CD200','samalizumab - never approved')}
print('%-8s %-9s %14s %12s   %s' % ('ANTIGEN','outcome','B-restriction','DLBCL ceiling','drug'))
rows=[]
for grp,dd in (('APPROVED',APPROVED),('FAILED',FAILED)):
    for hname,(cname,drug) in dd.items():
        b=brest.get(hname); c=CEIL.get(cname)
        rows.append((grp,hname,b,c))
        fb=('%.3f'%b) if b is not None else 'not in panel'
        fc=('%.3f'%c) if c is not None else 'not in panel'
        print('%-8s %-9s %14s %12s   %s' % (hname,grp,fb,fc,drug))

print()
rows_ext=list(rows)
for grp,dd in (('APPROVED',APPROVED_EXT),('FAILED',FAILED_EXT)):
    for hname,(cname,drug) in dd.items():
        b=brest.get(hname); c=CEIL.get(cname)
        rows_ext.append((grp,hname,b,c))
        fb=('%.3f'%b) if b is not None else 'not in panel'
        fc=('%.3f'%c) if c is not None else 'not in panel'
        print('%-8s %-9s %14s %12s   %s  [EXPANSION]' % (hname,grp,fb,fc,drug))

def sep(vals_a, vals_b):
    """P(a random APPROVED scores above a random FAILED)."""
    w=sum(1 for a in vals_a for b in vals_b if a>b)+0.5*sum(1 for a in vals_a for b in vals_b if a==b)
    return w/(len(vals_a)*len(vals_b))
for label,idx in (('B-lineage restriction (healthy donors)',2), ('tumor-vs-normal ceiling (DLBCL)',3)):
    A=[r[idx] for r in rows if r[0]=='APPROVED' and r[idx] is not None]
    F=[r[idx] for r in rows if r[0]=='FAILED'   and r[idx] is not None]
    s=sep(A,F)
    print('\n%s' % label)
    print('   approved mean %.3f   failed mean %.3f' % (np.mean(A), np.mean(F)))
    print('   P(approved ranks above failed) = %.3f  %s' % (s, '<-- perfect separation' if s==1.0 else ''))
    import math
    if s==1.0:
        print('   exact p under random ordering = 1/C(%d,%d) = %.4f'
              % (len(A)+len(F), len(A), 1/math.comb(len(A)+len(F), len(A))))
    print('   (n = %d approved, %d failed with a value)' % (len(A), len(F)))
print('\n' + '='*70)
print('EXPANDED SET (pre-specified rule, includes the inconvenient CD52 case)')
import math
for label,idx in (('B-lineage restriction (healthy donors)',2), ('tumor-vs-normal ceiling (DLBCL)',3)):
    A=[r[idx] for r in rows_ext if r[0]=='APPROVED' and r[idx] is not None]
    F=[r[idx] for r in rows_ext if r[0]=='FAILED'   and r[idx] is not None]
    s2=sep(A,F)
    print('\n%s' % label)
    print('   approved mean %.3f (n=%d)   failed mean %.3f (n=%d)' % (np.mean(A),len(A),np.mean(F),len(F)))
    print('   P(approved ranks above failed) = %.3f' % s2)
    print('   best possible exact p at this sample size = %.4f' % (1/math.comb(len(A)+len(F),len(A))))
json.dump({'b_restriction':brest,'floor':floor,'rows':[list(r) for r in rows_ext]},
          open(os.path.join(HERE,'approved_vs_failed.json'),'w'), indent=1)
