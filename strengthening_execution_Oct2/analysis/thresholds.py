"""Exact empirical deterministic-threshold policies; no clinical calibration."""
import math
import numpy as np

def clean(values):
    a=np.asarray(values,dtype=float)
    if a.ndim!=1 or not a.size or not np.isfinite(a).all() or (a<0).any():
        raise ValueError('Nonempty finite nonnegative measurement vector required')
    return np.sort(a)

def fraction_above(values,threshold):
    a=clean(values)
    return float((len(a)-np.searchsorted(a,threshold,side='left'))/len(a))

def auc(favorable,unfavorable):
    a=clean(favorable); b=clean(unfavorable)
    return float(np.mean([(np.searchsorted(b,x,'left')+.5*(np.searchsorted(b,x,'right')-np.searchsorted(b,x,'left')))/len(b) for x in a]))

def validate_count(count,denominator,unit):
    if not isinstance(count,int) or not isinstance(denominator,int) or count<0 or denominator<=0:
        raise ValueError('Explicit nonnegative integer count and positive denominator required')
    if unit not in ('patients','events'): raise ValueError('Specify patients or events')
    if unit=='patients' and count>denominator: raise ValueError('Patient count exceeds denominator')
    return count

def select_threshold(tumor,normal,budget,policy):
    if not 0<=budget<=1 or policy not in ('pooled','equal_donor','donor_envelope'):
        raise ValueError('Invalid budget or policy')
    if not tumor or not normal: raise ValueError('Both donor groups required')
    ts=[clean(x) for x in tumor]; ns=[clean(x) for x in normal]
    thresholds=np.r_[0,np.unique(np.concatenate(ts+ns)),np.inf]
    nt=np.array([(len(x)-np.searchsorted(x,thresholds,'left'))/len(x) for x in ns])
    tt=np.array([(len(x)-np.searchsorted(x,thresholds,'left'))/len(x) for x in ts])
    if policy=='pooled': criterion=np.average(nt,axis=0,weights=[len(x) for x in ns])
    elif policy=='equal_donor': criterion=nt.mean(axis=0)
    else: criterion=nt.max(axis=0)
    # All tails decrease with t: smallest feasible t optimizes every tumor tail.
    i=int(np.flatnonzero(criterion<=budget+1e-14)[0])
    t=float(thresholds[i])
    return dict(threshold=t,normal_mean=float(nt[:,i].mean()),normal_worst=float(nt[:,i].max()),
                tumor_mean=float(tt[:,i].mean()),tumor_worst=float(tt[:,i].min()),
                normal_pooled=float(np.average(nt[:,i],weights=[len(x) for x in ns])),
                tumor_pooled=float(np.average(tt[:,i],weights=[len(x) for x in ts])))
