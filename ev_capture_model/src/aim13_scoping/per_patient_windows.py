"""
Per-patient joint therapeutic window, benchmarked against the healthy donor-to-donor
noise floor measured in joint_window_and_controls.py step (A).
"""
import os, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
d=json.load(open(os.path.join(HERE,'two_comparator_results.json')))
pp=d['per_patient']; summ=d['summary']
ANT=list(summ.keys())
NOISE_MED, NOISE_P90 = 0.089, 0.220   # |AUC-0.5| between independent healthy donors

rows={}
for pid, rec in pp.items():
    joint={a: min(rec[a]['vsT'], rec[a]['vsNormalB']) for a in ANT}
    best=max(joint, key=joint.get)
    rows[pid]={'best':best,'best_joint':joint[best],
               'cd19':joint['CD19'],'cd20':joint['CD20'],'cd22':joint['CD22'],'cd79b':joint['CD79B'],
               'n_above_floor':sum(1 for a in ANT if joint[a]-0.5 > NOISE_P90)}
print('%-22s %-9s %10s %8s %8s %8s' % ('PATIENT','best','bestJoint','CD19','CD20','CD79B'))
for pid in sorted(rows, key=lambda p:-rows[p]['best_joint']):
    r=rows[pid]
    print('%-22s %-9s %10.3f %8.3f %8.3f %8.3f' % (pid.split('_DLB_')[-1]+' ('+pid.split('_')[0]+')', r['best'], r['best_joint'], r['cd19'], r['cd20'], r['cd79b']))

bj=np.array([rows[p]['best_joint'] for p in rows])
cd19=np.array([rows[p]['cd19'] for p in rows])
print('\nBest-available joint window per patient: median %.3f, range %.3f-%.3f' % (np.median(bj), bj.min(), bj.max()))
print('CD19 joint window per patient:            median %.3f, range %.3f-%.3f' % (np.median(cd19), cd19.min(), cd19.max()))
print('Patients where CD19 IS the best available target: %d/%d' % (sum(1 for p in rows if rows[p]['best']=='CD19'), len(rows)))
from collections import Counter
print('Best-target identity across the cohort:', Counter(rows[p]['best'] for p in rows).most_common())
above=[p for p in rows if rows[p]['best_joint']-0.5 > NOISE_P90]
print('\nPatients with ANY antigen whose joint window exceeds the healthy donor-variation')
print('90th percentile (|AUC-0.5| > %.3f): %d/%d' % (NOISE_P90, len(above), len(rows)))
above_med=[p for p in rows if rows[p]['best_joint']-0.5 > NOISE_MED]
print('...exceeds the MEDIAN healthy donor variation (%.3f): %d/%d' % (NOISE_MED, len(above_med), len(rows)))
json.dump(rows, open(os.path.join(HERE,'per_patient_windows.json'),'w'), indent=1)
