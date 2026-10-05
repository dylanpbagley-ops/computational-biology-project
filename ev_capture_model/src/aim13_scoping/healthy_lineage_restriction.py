"""
Testing the project's central thesis directly, in healthy people.

The thesis (section 4 of the results doc): the antigens that became approved drugs are not
selective for tumor -- they are restricted to a cell lineage the body can do without.

That claim is about NORMAL biology, so it can be tested entirely inside a healthy dataset,
with no cross-technology comparison. For each of 228 surface proteins measured on 161,764
PBMCs from 8 healthy donors (GSE164378, Hao et al. Cell 2021, with the authors' own
cell-type labels), compute:

  restriction = max over healthy lineages L of AUC(L vs every other healthy cell)
                -> 1.0 means the protein marks exactly one lineage; 0.5 means it is
                   spread evenly across all of them.
  home        = which lineage L achieves that maximum.
  breadth     = the number of lineages whose median exceeds the whole-sample median,
                a second, independent read on how widely it is expressed.

A target is safe-by-lineage only if it is BOTH highly restricted AND its home lineage is
expendable. This makes "dispensable lineage" a measured quantity rather than an argument.

Controls built into the panel: 4 isotype antibodies (Rat-IgG1-1, Rat-IgG1-2, Rat-IgG2b,
Rag-IgG2c) that bind nothing and must score ~0.5, and 11 antigens measured with two
independent clones, which must agree with each other.
"""
import os, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','healthy_pbmc_citeseq'))
rng=np.random.default_rng(71); SUB=6000

z=np.load(os.path.join(D,'healthy_adt.npz'), allow_pickle=True)
M=z['adt'].astype(np.float32); NAMES=[str(x) for x in z['names']]
L1=np.array([str(x) for x in z['celltype_l1']]); DON=np.array([str(x) for x in z['donor']])
X=M+1.0; X=np.log(X); X=X-X.mean(axis=1,keepdims=True)     # CLR per cell
del M

def auc(x,y):
    n1,n2=len(x),len(y)
    zz=np.concatenate([x,y]); o=np.argsort(zz,kind='mergesort'); zs=zz[o]
    r=np.empty(len(zz)); i=0
    while i<len(zz):
        j=i
        while j+1<len(zz) and zs[j+1]==zs[i]: j+=1
        r[i:j+1]=(i+j)/2.0+1; i=j+1
    rr=np.empty(len(zz)); rr[o]=r
    return (rr[:n1].sum()-n1*(n1+1)/2.0)/(n1*n2)
def sub(idx): return rng.choice(idx,SUB,replace=False) if len(idx)>SUB else idx

LIN=[l for l in sorted(set(L1)) if l!='other']
idx={l: np.flatnonzero(L1==l) for l in LIN}
print('healthy lineages:', {l:len(idx[l]) for l in LIN})
overall_med=np.median(X,axis=0)

res={}
for k,n in enumerate(NAMES):
    col=X[:,k]
    best=(0.5,None)
    for l in LIN:
        a=sub(idx[l]); b=sub(np.flatnonzero((L1!=l)&(L1!='other')))
        v=auc(col[a],col[b])
        if v>best[0]: best=(v,l)
    breadth=sum(1 for l in LIN if np.median(col[idx[l]])>overall_med[k])
    res[n]={'restriction':float(best[0]),'home':best[1],'breadth':int(breadth)}
json.dump(res, open(os.path.join(HERE,'healthy_lineage_restriction.json'),'w'), indent=1)

ISO=['Rat-IgG1-1','Rat-IgG1-2','Rat-IgG2b','Rag-IgG2c']
print('\n=== CONTROL 1: isotype antibodies (bind nothing; must sit near 0.5) ===')
for n in ISO: print('   %-12s restriction %.3f  home %-8s breadth %d' % (n,res[n]['restriction'],res[n]['home'],res[n]['breadth']))
print('   mean isotype restriction: %.3f' % np.mean([res[n]['restriction'] for n in ISO]))

DUP=[('CD4-1','CD4-2'),('CD38-1','CD38-2'),('CD138-1','CD138-2'),('CD45-1','CD45-2'),
     ('CD56-1','CD56-2'),('CD44-1','CD44-2'),('CD26-1','CD26-2'),('CD11b-1','CD11b-2'),
     ('CD133-1','CD133-2'),('CD275-1','CD275-2'),('CD3-1','CD3-2')]
print('\n=== CONTROL 2: same antigen, two independent antibody clones (must agree) ===')
d=[]
for a,b in DUP:
    d.append(abs(res[a]['restriction']-res[b]['restriction']))
    print('   %-9s %.3f (%-6s)  vs  %-9s %.3f (%-6s)   diff %.3f' %
          (a,res[a]['restriction'],res[a]['home'],b,res[b]['restriction'],res[b]['home'],d[-1]))
print('   median |difference| between clone pairs: %.3f' % np.median(d))

TARGETS={
 'CD19':'DLBCL - approved (axi-cel, tisa-cel, liso-cel, blinatumomab)',
 'CD20':'DLBCL - approved (rituximab, mosunetuzumab, glofitamab, epcoritamab)',
 'CD22':'DLBCL - approved (inotuzumab ozogamicin)',
 'CD79b':'DLBCL - approved (polatuzumab vedotin)',
 'CD269':'Myeloma - approved (ide-cel, cilta-cel, teclistamab) [BCMA]',
 'CD38-1':'Myeloma - approved (daratumumab, isatuximab)',
 'CD319':'Myeloma - approved (elotuzumab) [SLAMF7]',
 'CD307d':'Myeloma - in trials (cevostamab) [FcRH5]',
 'CD138-1':'Myeloma - DISCONTINUED (indatuximab ravtansine)',
 'CD30':'Lymphoma - approved (brentuximab vedotin)',
 'CD33':'AML - approved (gemtuzumab ozogamicin)' if 'CD33' in res else None,
 'CD117':'AML/mast - not successfully targeted',
 'CD44-1':'broadly expressed - never successfully targeted',
 'CD45-1':'pan-leukocyte - not targetable',
 'HLA-DR':'broadly expressed on antigen-presenting cells',
 'CD48':'broadly expressed on leukocytes',
 'CD27':'T/B memory marker',
}
print('\n=== real immunotherapy targets, measured in healthy people ===')
print('%-9s %12s %-9s %8s  %s' % ('PROTEIN','restriction','home','breadth','clinical status'))
for n,note in TARGETS.items():
    if note is None or n not in res: continue
    r=res[n]
    print('%-9s %12.3f %-9s %8d  %s' % (n,r['restriction'],r['home'],r['breadth'],note))
