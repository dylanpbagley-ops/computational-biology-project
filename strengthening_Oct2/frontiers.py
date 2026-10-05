"""Finite-distribution response bounds, not a clinical outcome model.

Nondecreasing responses are convex combinations of ordered threshold responses.
Arbitrary responses obey a different bound given by fractional knapsack on p/q.
The weighted mixture is an exact identity only for a shared response and specified
nonnegative weights. The empirical expression score is not derived from it.
"""
import json,pathlib
import numpy as np
def probability(x):
    x=np.asarray(x,dtype=float)
    if x.ndim!=1 or len(x)==0 or not np.isfinite(x).all() or (x<0).any() or not np.isclose(x.sum(),1,atol=1e-12,rtol=0):raise ValueError('Expected finite nonnegative probability vector summing to one')
    return x
def upper_at(p,q,budget,kind):
    p=probability(p);q=probability(q)
    if p.shape!=q.shape or not np.isfinite(budget) or not 0<=budget<=1:raise ValueError('Incompatible distributions or budget')
    if kind=='arbitrary':
        value=float(p[q==0].sum());remaining=float(budget)
        positive=np.flatnonzero(q>0)
        for i in positive[np.argsort(-(p[positive]/q[positive]))]:
            amount=min(remaining,float(q[i]));value+=amount*p[i]/q[i];remaining-=amount
        return float(value)
    if kind!='monotone':raise ValueError('Unknown response class')
    # Every feasible monotone response mixes upper-tail thresholds (including constants).
    x=np.r_[0,np.cumsum(q[::-1])];y=np.r_[0,np.cumsum(p[::-1])]
    best=max(float(v) for u,v in zip(x,y) if u<=budget+1e-12)
    # A linear objective with a single budget constraint requires at most two vertices.
    for i in range(len(x)):
        for j in range(i+1,len(x)):
            if x[i]<=budget<=x[j] and x[j]>x[i]:
                t=(budget-x[i])/(x[j]-x[i]);best=max(best,float((1-t)*y[i]+t*y[j]))
    return best
def weighted_reference(compartments,prevalence,weights):
    pc=np.asarray([probability(p) for p in compartments]);q=probability(prevalence);w=np.asarray(weights,dtype=float)
    if pc.ndim!=2 or len(pc)!=len(q) or w.shape!=q.shape or not np.isfinite(w).all() or (w<0).any():raise ValueError('Invalid compartment weights')
    mass=float(q@w)
    return (None,0.) if mass==0 else (((q*w)[:,None]*pc).sum(axis=0)/mass,mass)
def tau_from_share(share,n):
    if n<2 or not 1/n-1e-12<=share<=1:raise ValueError('Invalid maximum share or type count')
    return (n-1/share)/(n-1)
def main():
    p=[.45,.1,.45];q=[.1,.8,.1]
    out={'nonmonotone_example':{'budget':.2,'arbitrary':upper_at(p,q,.2,'arbitrary'),'monotone':upper_at(p,q,.2,'monotone')},'multiple_compartments':{'p':[.5,.5],'normal_1':[1,0],'normal_2':[0,1],'each_budget':.5,'min_of_individual_ceilings':min(upper_at([.5,.5],n,.5,'arbitrary') for n in [[1,0],[0,1]]),'joint_optimum':.5,'joint_reason':'Both constraints require f[0] <= 0.5 and f[1] <= 0.5; f=(0.5,0.5) attains 0.5.'},'scope':'Abstract response geometry; no therapeutic construction or clinical validation.'}
    path=pathlib.Path(__file__).with_name('frontier_results.json');path.write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
