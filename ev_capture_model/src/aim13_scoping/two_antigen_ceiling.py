"""
THE TWO-ANTIGEN EXTENSION.

The single-antigen theorem says: for any agent whose killing rises with one antigen's
density, the achievable (tumor killed, normal spared) points lie under the ROC of that
antigen's tumor-vs-normal distributions. With TWO antigens each cell is a point on a
plane rather than a position on a line, and two populations that overlap completely along
either axis alone can still be separated on the plane. So the single-antigen ceiling stops
applying, and a new hierarchy takes its place.

THE HIERARCHY (each family contains the one above it, so the chain is forced):

  1. SINGLE      kill if  x > a                      -- one antigen, one threshold
  2. AND gate    kill if  x > a AND y > b            -- contains (1) at b = -inf
  3. NOT gate    kill if  x > a AND y < b            -- contains (1) at b = +inf
  4. LIKELIHOOD  kill if  P_tumor(x,y)/P_normal(x,y) > c
                                                     -- optimal over ALL rules, by
                                                        Neyman-Pearson; contains 2 and 3

Therefore  J1 <= max(J2, J3) <= J4  for Youden's J = sensitivity + specificity - 1,
and the gap J4 - J1 is exactly what combining two antigens can buy for that pair.

THE METHODOLOGICAL TRAP, AND THE GUARD AGAINST IT:
a two-antigen rule has more free parameters than a one-antigen rule, so on the data used
to fit it, it will ALWAYS look at least as good -- whether or not the gain is real. Every
number here is therefore fitted on a random half of the cells and scored on the held-out
half. A gain that survives the split is real; one that does not is overfitting, and is
reported as such.
"""
import os, sys, glob, json, itertools, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
CACHE=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','gated_cache'))
rng=np.random.default_rng(97)
NBIN=16; SUB=6000

def split(a):
    a=np.asarray(a,float)
    if len(a)>2*SUB: a=a[rng.choice(len(a),2*SUB,replace=False)]
    p=rng.permutation(len(a)); h=len(a)//2
    return a[p[:h]], a[p[h:]]

def edges(v, nb=NBIN):
    """Bin edges for the 2-D histograms.

    A CORRECTION TO AN EARLIER COMMENT IN THIS FILE, kept rather than deleted because the
    mistake is instructive. This function once carried a comment calling the +/- infinity
    outer bins "the third bug", and claiming the fix moved the best single-antigen J from
    0.618 down to 0.410. That diagnosis was WRONG. The change itself is correct and is kept
    -- a test cell outside the train range should not be silently dropped -- but it did NOT
    move the number, because dropped cells were never the cause.

    The real explanation is median-of-maximum versus maximum-of-median. 0.618 is the median
    over patients of the best single antigen chosen PER PATIENT; 0.410 is the best value for
    any single antigen held fixed across all patients. Both are correct and they answer
    different questions. Verified directly on CD27+FasR: median-of-max 0.616 against
    per-antigen medians of 0.384 and 0.412."""
    q=np.unique(np.percentile(v, np.linspace(0,100,nb+1)))
    if len(q)<3: q=np.linspace(float(np.min(v)), float(np.max(v))+1e-9, 4)
    q=q.astype(float); q[0]=-np.inf; q[-1]=np.inf
    return q

def J_single(tr_t, tr_n, te_t, te_n):
    """best single-antigen threshold, fitted on train, scored on test."""
    best=(-9,None)
    for col in range(tr_t.shape[1]):
        cand=np.unique(np.percentile(np.concatenate([tr_t[:,col],tr_n[:,col]]), np.linspace(2,98,40)))
        for a in cand:
            j=np.mean(tr_t[:,col]>a)-np.mean(tr_n[:,col]>a)
            if j>best[0]: best=(j,(col,a))
    col,a=best[1]
    return float(np.mean(te_t[:,col]>a)-np.mean(te_n[:,col]>a)), best[1]

def _hist2(t, n, i, j, ex, ey):
    Ht,_,_=np.histogram2d(t[:,i],t[:,j],bins=[ex,ey])
    Hn,_,_=np.histogram2d(n[:,i],n[:,j],bins=[ex,ey])
    return Ht/max(Ht.sum(),1), Hn/max(Hn.sum(),1)

def _rev_cum2(H):
    """P(x > edge_a, y > edge_b) for every (a,b) at once."""
    return np.cumsum(np.cumsum(H[::-1,::-1],axis=0),axis=1)[::-1,::-1]

def _rev_fwd_cum2(H):
    """P(x > edge_a, y < edge_b) for every (a,b) at once."""
    return np.cumsum(np.cumsum(H[::-1,:],axis=0)[::-1,:],axis=1)

def gates_and_lr(tr_t,tr_n,te_t,te_n,i,j,ex,ey):
    """Fit every rule family on train; report BOTH the train (in-sample) and the held-out
    score for each. Returns Youden's J.

    TWO BUGS FOUND AND FIXED (September 7, 2026) -- both were inflating the apparent value
    of combining antigens, and both were caught by the myeloma replication:

    BUG 1, the single-antigen baseline only tested "kill where the marker is HIGH".
    Plenty of real markers are LOW on tumor (CD3E here, IgD in the lymphoma data), for which
    the correct one-antigen rule is "kill where it is LOW". Omitting that made a NOT gate
    look like a two-antigen gain when it was really a one-antigen effect with the opposite
    sign. The baseline now tries both directions.

    BUG 2, the likelihood-ratio rule is the optimum in the POPULATION but is not estimable
    from this much data. It selects an arbitrary subset of 16x16 = 256 histogram bins --
    hundreds of effective free parameters against a gate's two thresholds -- so it wins
    in-sample by construction and generalises worse. Both numbers are now reported: the
    in-sample value verifies the theorem's ordering, the held-out value says what actually
    works.
    """
    Ht,Hn=_hist2(tr_t,tr_n,i,j,ex,ey)
    Et,En=_hist2(te_t,te_n,i,j,ex,ey)
    out={}
    def score(sel_tr):
        return float(Ht[sel_tr].sum()-Hn[sel_tr].sum()), float(Et[sel_tr].sum()-En[sel_tr].sum())
    # AND: x above a, y above b
    A_tr=_rev_cum2(Ht)-_rev_cum2(Hn)
    k=np.unravel_index(np.argmax(A_tr),A_tr.shape)
    selA=np.zeros(Ht.shape,bool); selA[k[0]:,k[1]:]=True
    out['AND_train'],out['AND']=score(selA)
    # NOT: x above a, y below b  (both orientations)
    best=None
    for flip in (False,True):
        H1,H2=(Ht,Hn) if not flip else (Ht.T,Hn.T)
        N_tr=_rev_fwd_cum2(H1)-_rev_fwd_cum2(H2)
        kk=np.unravel_index(np.argmax(N_tr),N_tr.shape)
        sel=np.zeros(H1.shape,bool); sel[kk[0]:,:kk[1]+1]=True
        if flip: sel=sel.T
        v=score(sel)
        if best is None or v[0]>best[0]: best=v
    out['NOT_train'],out['NOT']=best
    # SINGLE, two versions, because they answer different questions and BOTH are needed:
    #
    #   SINGLE_hi   kill only where the antigen is HIGH. This is what ONE receptor can
    #               physically do -- a CAR kills cells that HAVE its antigen; you cannot
    #               build a receptor that fires on the ABSENCE of a protein. This is the
    #               honest ENGINEERING baseline: 1 receptor vs 2.
    #   SINGLE_any  kill where the antigen is high OR where it is low. Not buildable with
    #               one receptor (killing on absence needs an inhibitory receptor plus an
    #               activating one, i.e. two), but it is the right STATISTICAL control for
    #               "is the two-antigen gain just an inverted marker with the sign flipped?"
    #
    # Reporting only the first overstates the value of combining; reporting only the second
    # understates it. Both are computed.
    bt=bs=-9; btr=-9; hi_tr=-9; hi_te=-9
    for axis in (1,0):
        H1t=Ht.sum(axis); H1n=Hn.sum(axis); E1t=Et.sum(axis); E1n=En.sum(axis)
        for direction in ('hi','lo'):
            if direction=='hi':
                ct=np.cumsum(H1t[::-1])[::-1]; cn=np.cumsum(H1n[::-1])[::-1]
                et=np.cumsum(E1t[::-1])[::-1]; en=np.cumsum(E1n[::-1])[::-1]
            else:
                ct=np.cumsum(H1t); cn=np.cumsum(H1n)
                et=np.cumsum(E1t); en=np.cumsum(E1n)
            k1=int(np.argmax(ct-cn))
            if ct[k1]-cn[k1]>btr: btr=float(ct[k1]-cn[k1]); bs=float(et[k1]-en[k1])
            if direction=='hi' and ct[k1]-cn[k1]>hi_tr:
                hi_tr=float(ct[k1]-cn[k1]); hi_te=float(et[k1]-en[k1])
    out['SINGLE_train']=btr; out['SINGLE_in_pair']=bs
    out['SINGLEHI_train']=hi_tr; out['SINGLEHI']=hi_te
    # Likelihood ratio: population optimum, finite-sample overfitter
    LR=(Ht+1e-6)/(Hn+1e-6)
    order=np.argsort(LR.ravel())[::-1]
    ct=np.cumsum(Ht.ravel()[order]); cn=np.cumsum(Hn.ravel()[order])
    b=int(np.argmax(ct-cn))
    sel=np.zeros(LR.size,bool); sel[order[:b+1]]=True; sel=sel.reshape(LR.shape)
    out['LR_train'],out['LR']=score(sel)
    return out

def main():
    files={os.path.basename(p).replace('.npz',''):p for p in glob.glob(os.path.join(CACHE,'*.npz'))}
    ANT=[str(x) for x in np.load(next(iter(files.values())),allow_pickle=True)['antigens']]
    bat=lambda k:k.split('_')[0]
    tumors=sorted(k for k in files if '_DLB_' in k)
    rlns=sorted(k for k in files if '_RLN_' in k and 'pooled' not in k)
    normalB={b: np.concatenate([np.load(files[k],allow_pickle=True)['B'] for k in rlns if bat(k)==b],axis=0)
             for b in sorted({bat(k) for k in rlns})}
    STORE=os.path.join(HERE,'two_antigen_perpatient.json')
    done=json.load(open(STORE)) if os.path.exists(STORE) else {}
    npat=int(sys.argv[1]) if len(sys.argv)>1 else 8
    use=[t for t in tumors if t not in done][:npat]
    print('cohort %d, already done %d, processing %d now' % (len(tumors),len(done),len(use)))
    # COMPARATOR = batch-matched RLN NORMAL B CELLS ONLY.
    #
    # Two reasons, and the first is a correctness fix. The B gate is
    # "(CD19 high OR CD20 high) AND CD3 low", so CD19 and CD20 separate B cells from T and
    # myeloid cells BY CONSTRUCTION -- any pair containing them was partly circular when the
    # comparator pooled all three compartments. Normal B cells are gated by exactly the same
    # rule as tumor B cells, so the gate cannot manufacture a difference between them and no
    # marker in the panel is circular here.
    # Second, sections 3-4 established that normal B cells are the BINDING compartment for
    # every B-lineage target anyway, so this is also the comparator that matters clinically.
    acc={}
    for t in use:
        d=np.load(files[t],allow_pickle=True)
        T=np.asarray(normalB[bat(t)])
        tr_t,te_t=split(d['B']); tr_n,te_n=split(T)
        j1,which=J_single(tr_t,tr_n,te_t,te_n)
        acc.setdefault('single',[]).append(j1)
        E=[edges(np.concatenate([tr_t[:,c],tr_n[:,c]])) for c in range(len(ANT))]
        for i,jx in itertools.combinations(range(len(ANT)),2):
            r=gates_and_lr(tr_t,tr_n,te_t,te_n,i,jx,E[i],E[jx])
            acc.setdefault('%s|%s'%(ANT[i],ANT[jx]),[]).append((r['AND'],r['NOT'],r['LR'],r['SINGLE_in_pair'],r['AND_train'],r['NOT_train'],r['LR_train'],r['SINGLE_train'],r['SINGLEHI'],r['SINGLEHI_train']))
        done[t]={k:v[-1] for k,v in acc.items() if k!='single'}
        json.dump(done, open(STORE,'w'))
        print('  processed %s (%d/%d)' % (t,len(done),len(tumors)), flush=True)
    # rebuild from the on-disk store so every patient ever processed is included
    pairs={}
    for t,rec in done.items():
        for k,v in rec.items(): pairs.setdefault(k,[]).append(v)
    out={'n_patients': len(done)}
    rows=[]
    for k,v in pairs.items():
        k=tuple(k.split('|')); v=np.array(v)
        rows.append({'pair':list(k),'AND':float(np.median(v[:,0])),'NOT':float(np.median(v[:,1])),
                     'LR':float(np.median(v[:,2])),'SINGLE':float(np.median(v[:,3])),
                     'AND_tr':float(np.median(v[:,4])),'NOT_tr':float(np.median(v[:,5])),
                     'LR_tr':float(np.median(v[:,6])),'SINGLE_tr':float(np.median(v[:,7])),
                     'SINGLEHI':float(np.median(v[:,8])),'SINGLEHI_tr':float(np.median(v[:,9]))})
    out['pairs']=rows
    json.dump(out, open(os.path.join(HERE,'two_antigen_results.json'),'w'))
    print('\npairs evaluated: %d, patients accumulated: %d of %d' % (len(rows), len(done), len(tumors)))

if __name__=='__main__': main()
