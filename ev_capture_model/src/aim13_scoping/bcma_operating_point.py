"""
Placing the Bluhm et al. 2018 BCMA CAR-T activation threshold on the myeloma ROC --
the same exercise section 11 already ran for CD19 (CLN-978 / blinatumomab), now
extended to the second disease with a third real, independently calibrated agent.

Bluhm et al. (Mol Ther 2018, doi:10.1016/j.ymthe.2018.06.012, PMID 30078440) quantified
BCMA density with QuantiBRITE PE calibration beads -- the same absolute-molecule standard
used for CLN-978/blinatumomab -- and reported that BCMA CAR T cells lysed target cells with
an activation threshold "in the range of 100 BCMA molecules" (JeKo-1, the lowest-density
responding line, carried <=~100 molecules/cell; NCI-H929 carried ~25,000).

One honest methods difference from CLN-978/blinatumomab, stated plainly: those two were fit
to a continuous dose-response titration, giving a precise k_intra. Bluhm's ~100-molecule
number is read off a small set of discretely bead-quantified cell lines correlated with
lysis/IFN-gamma output, not a fitted curve. It is real, absolute-molecule, peer-reviewed
data -- just a coarser derivation. Treated here exactly as section 11 treats its own three
agents: an operating-point placement, not a claim of clinical performance.

Data: this project's own myeloma CITE-seq ADT panel (11 samples, ev_capture_data/myeloma_adt),
raw (non-CLR) CD269/BCMA antibody-derived-tag counts. Populations follow myeloma_ceiling.py's
gating exactly (tumor = malignant plasma cells; normal compartments = T, myeloid, NK, normal B).

Expected honest result, stated up front: like section 11, this is very likely ANOTHER
absolute-calibration gap, not a working operating point -- CITE-seq ADT counts (raw UMI/read
counts, median ~7, max ~221 for CD269 in this data) are not molecule-calibrated, so n*=100
cannot be compared to them on a common scale without a bridging assumption. That gap, found a
second time in a second disease and a second technology, is itself the reportable result.
"""
import os, glob, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
ADT=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','myeloma_adt'))

N_STAR_BCMA = 100.0   # Bluhm et al. 2018, QuantiBRITE-calibrated activation threshold

files=sorted(glob.glob(os.path.join(ADT,'*.npz')))
NAMES=list(np.load(files[0],allow_pickle=True)['names'])
ix={n:i for i,n in enumerate(NAMES)}

def gate(M):
    cd3,cd14,cd33,cd56,cd19,cd138,bcma = (M[:,ix[k]].astype(float) for k in
        ['CD3E','CD14','CD33','CD56','CD19','CD138','CD269'])
    hi=lambda v,p=75: v > np.percentile(v,p)
    T  = hi(cd3,70)
    MY = (~T) & (hi(cd14,80) | hi(cd33,80))
    NK = (~T) & (~MY) & hi(cd56,85)
    NB = (~T) & (~MY) & (~NK) & hi(cd19,85) & (cd138<=np.percentile(cd138,70))
    TU = (~T) & (~MY) & (~NK) & (~NB) & (hi(cd138,70) | hi(bcma,80))
    return {'T':T,'MY':MY,'NK':NK,'NB':NB,'TU':TU}

print('=== RAW (non-CLR) CD269/BCMA ADT counts across the 11-sample myeloma cohort ===')
print('(raw counts because Bluhm\'s n* is in absolute molecules; CLR log-ratios have no')
print(' molecule scale at all, so raw counts are the closer -- if still not equivalent -- unit)\n')

per_sample=[]
for f in files:
    sid=os.path.basename(f)[:-4]
    M=np.load(f,allow_pickle=True)['adt'].astype(float)
    bcma_raw = M[:,ix['CD269']]
    pops = gate(M)
    tu_n = int(pops['TU'].sum())
    if tu_n < 100:
        print('  skip %s: too few tumor cells (%d)' % (sid, tu_n)); continue
    row={'sample':sid,'n_tumor':tu_n}
    row['tumor_median']=float(np.median(bcma_raw[pops['TU']]))
    row['tumor_p90']=float(np.percentile(bcma_raw[pops['TU']],90))
    row['tumor_max']=float(bcma_raw[pops['TU']].max())
    for lab in ('T','MY','NK','NB'):
        mask=pops[lab]
        if mask.sum()<50: row[lab+'_median']=np.nan; continue
        row[lab+'_median']=float(np.median(bcma_raw[mask]))
    per_sample.append(row)

print('%-6s %6s %10s %8s %8s  %6s %6s %6s %6s' % ('sample','n_tu','tu_median','tu_p90','tu_max','T_med','MY_med','NK_med','NB_med'))
for r in per_sample:
    print('%-6s %6d %10.1f %8.1f %8.1f  %6.1f %6.1f %6.1f %6.1f' % (
        r['sample'], r['n_tumor'], r['tumor_median'], r['tumor_p90'], r['tumor_max'],
        r['T_median'], r['MY_median'], r['NK_median'], r['NB_median']))

all_tumor_max = max(r['tumor_max'] for r in per_sample)
all_tumor_p90 = np.median([r['tumor_p90'] for r in per_sample])
print('\nHighest single-cell raw CD269 ADT count observed anywhere in the cohort: %.0f' % all_tumor_max)
print('Median of per-sample tumor-compartment 90th-percentile raw CD269 ADT count: %.1f' % all_tumor_p90)
print('Bluhm et al. n* (absolute molecules/cell): %.0f' % N_STAR_BCMA)
print()
print('=== VERDICT ===')
print('The raw ADT scale in this dataset tops out in the hundreds (max %.0f across %d samples);' % (all_tumor_max, len(per_sample)))
print('Bluhm\'s n*=100 is stated in absolute molecules/cell, calibrated by QuantiBRITE beads.')
print('These are not the same unit and cannot be equated without an independent bridging')
print('calibration (e.g. running QuantiBRITE beads on cells from THIS cohort, which this')
print('project does not have access to do). This is the same absolute-calibration gap')
print('section 11 already reported for CD19/CLN-978/blinatumomab -- found again, independently,')
print('in a second disease, second technology (CITE-seq vs mass cytometry), second lab.')
print('That consistency is itself informative: it is not a one-off data problem with the')
print('DLBCL CyTOF panel, it is a structural feature of how these datasets are generated.')
