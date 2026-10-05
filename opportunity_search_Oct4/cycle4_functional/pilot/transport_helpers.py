import numpy as np

def count_only_positions(full_header):
    names=[name for name in full_header if name.startswith('MUS_')]
    return names,{name:i for i,name in enumerate(names)}

def normalized_log_frequency(counts):
    arr=np.asarray(counts,dtype=np.float64)
    if arr.ndim not in (1,2) or arr.shape[0]==0: raise ValueError('counts must be nonempty guide-by-sample array')
    if not np.isfinite(arr).all() or (arr<0).any(): raise ValueError('counts must be finite and nonnegative')
    totals=arr.sum(axis=0); n_guides=arr.shape[0]
    return np.log2(1e6*(arr+0.5)/(totals+0.5*n_guides))

def summarize_envelope(values):
    vals=np.asarray(values,dtype=np.float64)
    if vals.shape!=(9,) or not np.isfinite(vals).all(): raise ValueError('expected exactly nine finite culture-pair values')
    med=float(np.median(vals)); lo=float(np.min(vals)); hi=float(np.max(vals))
    if hi<0: cls='sign_stable_sensitizing'
    elif lo>0: cls='sign_stable_resistant'
    else: cls='mixed'
    return med,lo,hi,cls

def guide_filter_reasons(guide_id,gene_symbol,library_gene,t0_count,control_counts,duplicate_id=False,missing_t0=False,missing_control=False):
    reasons=[]
    if not guide_id: reasons.append('missing_guide_id')
    if duplicate_id: reasons.append('duplicate_guide_id')
    if not gene_symbol: reasons.append('missing_gene_symbol')
    # Gene.symbol in the measured count matrix is the frozen aggregation source; library annotation is QC only.
    if missing_t0 or t0_count is None: reasons.append('missing_T0_count')
    elif t0_count<20: reasons.append('T0_below_20')
    if missing_control or len(control_counts)!=3: reasons.append('missing_primary_control_count')
    else:
        for i,value in enumerate(control_counts,1):
            if value is None: reasons.append(f'control_rep_{i}_missing')
            elif value<20: reasons.append(f'control_rep_{i}_below_20')
    return reasons

def exact_source_rows(rows,source):
    return [r for r in rows if r.get('source')==source]
