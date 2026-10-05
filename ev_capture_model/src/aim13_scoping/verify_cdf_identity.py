import sys, os, numpy as np
SRC = os.path.expanduser('~/mnt/outputs/ev_capture_model/src')
sys.path.insert(0, SRC)
import aim12_single_cell_cytof_reordering as a12
from aim10_car_t_kill_probability_model import escape_fraction, kill_probability

RAW = os.path.expanduser('~/mnt/outputs/ev_capture_data/raw_cytof')
f = os.path.join(RAW, 'CytB2_DLB_P5904.fcs')

cd19 = a12.load_marker_values_from_fcs_file(f, 'CD19')
print('CD19 n=', len(cd19), 'mean=%.3f median=%.3f zerofrac=%.4f max=%.3f' % (
    cd19.mean(), np.median(cd19), float((cd19<=0).mean()), cd19.max()))
print('EXPECTED from doc: n=47568 mean=62.701 median=16.895 zero=0.1469 max=2041.712')

# ---- 1. verify the n* identity ----
def n_star(car_avidity, k_intra, z, tau):
    return np.log(1.0 + tau/(z*(1.0-tau))) / (car_avidity*k_intra)

print('\n--- identity check: escape_fraction == fraction(n < n*) ---')
ok = True
for k_intra in [1e-5, 1e-4, 1.2e-3, 0.01]:
    for tau in [0.1,0.3,0.5,0.7,0.9]:
        for z in [0.5, 1.0, 2.0]:
            ef = escape_fraction(cd19, 1.0, k_intra, z, kill_threshold=tau)['escape_fraction']
            ns = n_star(1.0, k_intra, z, tau)
            cdf = float(np.mean(cd19 < ns))
            if abs(ef-cdf) > 1e-12:
                ok = False
                print('MISMATCH', k_intra, tau, z, ef, cdf, ns)
print('identity holds exactly for every combination tested:', ok)

# show n* values for the 3 real aim10 fits
print('\n--- n* for the real calibrations (car_avidity=1, z=1) ---')
import aim10_car_t_kill_probability_model as a10
fits = a10.fit_all_real_thresholds()
print(type(fits))
try:
    for k,v in fits.items():
        print(' ', k, v)
except Exception as e:
    print(fits)
