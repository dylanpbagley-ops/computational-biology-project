"""
CLDN18.2 (claudin 18.2) real dose-response calibration, extracted from
Barrett et al. 2024, Clinical Cancer Research 30(23):5413-5429,
DOI: 10.1158/1078-0432.CCR-24-1853, Figure 5D.

Values read directly off the zoomed figure image (broken-axis bar chart,
NOT log scale -- bottom segment 0-4,000 linear, break, then 4,000-500,000
linear on top segment). These are careful visual estimates from the
published figure, not exact extracted numbers (no data table was found in
Supplementary Data1) -- explicitly labeled as such throughout.

ABC = Antibody Binding Capacity, quantitated via Quantum Simply Cellular
beads (Bangs Laboratories) -- the same real calibration-bead methodology
used in the mesothelin/Chu et al. 2024 paper this project already uses.
This is a real, same-assay-type molecules/cell proxy.

Cell line identity note: NUGC4 and SNU601 are gastric cancer lines (not
pancreatic) -- included in the fit for statistical power since they are
real points on the SAME real dose-response curve, but flagged so pancreatic-
specific claims can be made separately using only the 3 pancreatic lines
(AsPC1+CLDN18.2, PaTu8988s HS, BxPC3+CLDN18.2) if a stricter pancreatic-only
claim is wanted.
"""

import numpy as np
from scipy import stats, optimize
import pickle

# (cell_line, ABC molecules/cell [estimated from Fig 5D], % cytolysis at 48h [estimated], is_pancreatic)
CLDN182_FIGURE_5D_DATA = [
    ("PaTu8988s_18.2_KO", 0,        1.5,  True),   # negative control, pancreatic
    ("DAN-G",             850,      15,   False),  # gastric
    ("IM-95",             2750,     65,   False),  # gastric
    ("NUGC4",             45000,    68,   False),  # gastric
    ("ASPC1_CLDN18.2",    45000,    68,   True),   # pancreatic (engineered overexpression)
    ("SNU601",            50000,    77,   False),  # gastric
    ("PaTu8988s_HS",      185000,   88,   True),   # pancreatic (native high expresser)
    ("BxPC3_CLDN18.2",    350000,   82,   True),   # pancreatic (engineered overexpression)
]

print("=" * 78)
print("CLDN18.2 REAL DOSE-RESPONSE DATA (estimated from Fig 5D, Barrett et al. 2024)")
print("=" * 78)
for name, abc, cytolysis, is_panc in CLDN182_FIGURE_5D_DATA:
    tag = "PANCREATIC" if is_panc else "gastric"
    print(f"  {name:22s}  ABC~{abc:>7,}  cytolysis~{cytolysis:>5.1f}%   [{tag}]")

abc_vals = np.array([d[1] for d in CLDN182_FIGURE_5D_DATA], dtype=float)
cytolysis_vals = np.array([d[2] for d in CLDN182_FIGURE_5D_DATA], dtype=float)
is_panc = np.array([d[3] for d in CLDN182_FIGURE_5D_DATA])

# Exclude the true-zero KO point for correlation/log-fit purposes (log(0) undefined;
# it's still a valid "no antigen -> no lysis" sanity check, reported separately)
nonzero_mask = abc_vals > 0
abc_nz = abc_vals[nonzero_mask]
cyto_nz = cytolysis_vals[nonzero_mask]

# Spearman correlation across all non-zero points (mirrors the mesothelin paper's own
# Spearman r=0.623 statistic, for direct comparability)
rho, pval = stats.spearmanr(abc_nz, cyto_nz)
print(f"\nSpearman correlation (ABC vs %cytolysis, n={len(abc_nz)} non-KO points): "
      f"r={rho:.3f}, p={pval:.4f}")

# Pancreatic-only subset
panc_mask = is_panc & nonzero_mask
abc_panc = abc_vals[panc_mask]
cyto_panc = cytolysis_vals[panc_mask]
print(f"\nPancreatic-only subset (n={len(abc_panc)}): "
      f"{list(zip([d[0] for d in CLDN182_FIGURE_5D_DATA if d[3] and d[1]>0], abc_panc, cyto_panc))}")

# ---------------------------------------------------------------------------
# Fit k_intra using this project's own capture_probability functional form:
# theta(n, sigma=1, k_intra, z=1) = zqb/(1+zqb), qb = exp(n*k_intra) - 1
# Fit k_intra (single free parameter, sigma=z=1 fixed per this project's own
# convention) via least-squares against the REAL (ABC, %cytolysis/100) pairs.
# ---------------------------------------------------------------------------
import sys
sys.path.insert(0, "/sessions/intelligent-practical-cannon/mnt/Desktop--outputs/ev_capture_model/src")
from binding_model import capture_probability

def residuals(log_k_intra):
    k_intra = np.exp(log_k_intra[0])
    theta_pred = capture_probability(abc_nz, 1.0, k_intra, 1.0)
    return theta_pred - (cyto_nz / 100.0)

# Reasonable starting guess: k_intra ~ ln(2)/(typical mid-range ABC ~40,000)
x0 = [np.log(np.log(2) / 40000.0)]
result = optimize.least_squares(residuals, x0=x0, bounds=([-30], [0]))
k_intra_fit = float(np.exp(result.x[0]))
n_half_fit = np.log(2) / k_intra_fit  # the antigen count at which theta=0.5

print(f"\n" + "=" * 78)
print("FITTED k_intra (least-squares against real Fig 5D data)")
print("=" * 78)
print(f"  k_intra = {k_intra_fit:.6e}")
print(f"  Implied half-kill threshold (n_half, theta=0.5 at sigma=z=1): "
      f"{n_half_fit:,.0f} molecules/cell")

# Report fit quality
theta_pred_final = capture_probability(abc_nz, 1.0, k_intra_fit, 1.0)
residual_final = theta_pred_final - (cyto_nz/100.0)
rmse = np.sqrt(np.mean(residual_final**2))
print(f"  RMSE (predicted vs actual cytolysis fraction): {rmse:.4f}")
print(f"\n  Predicted vs actual, per cell line:")
names_nz = [d[0] for d in CLDN182_FIGURE_5D_DATA if d[1] > 0]
for name, abc, actual, pred in zip(names_nz, abc_nz, cyto_nz, theta_pred_final*100):
    print(f"    {name:22s} ABC~{abc:>7,.0f}  actual~{actual:5.1f}%  predicted={pred:5.1f}%")

with open("/sessions/intelligent-practical-cannon/mnt/outputs/cldn182_calibration/cldn182_fit_results.pkl", "wb") as f:
    pickle.dump({
        "data": CLDN182_FIGURE_5D_DATA,
        "spearman_r": rho, "spearman_p": pval,
        "k_intra_fit": k_intra_fit, "n_half_fit": n_half_fit, "rmse": rmse,
    }, f)
print("\nSaved fit results.")
