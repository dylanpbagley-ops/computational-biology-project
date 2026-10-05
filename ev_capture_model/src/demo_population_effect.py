"""
demo_population_effect.py

A demonstration (not a formal test) that shows WHY the population-level
model matters: it compares the population-average yield for two
populations that have the SAME MEAN receptor count but different
distribution SHAPES (Poisson vs. negative binomial / overdispersed).

If the "fold-change" approach (which only looks at the mean) were a good
enough approximation, these two populations -- despite having very
different shapes -- should give similar yields. If the model shows they
give meaningfully different yields, that's a concrete, visual demonstration
of the project's central claim: shape matters, not just the mean.

This produces a plot comparing yield vs. ligand density (sigma) for:
    1. A Poisson-distributed population (mean_n receptors, "textbook" shape)
    2. A negative-binomial population with the SAME mean but heavy
       zero-inflation/overdispersion (closer to real EV data per the
       handoff)

Run with: python demo_population_effect.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from population_model import (
    poisson_distribution,
    negative_binomial_distribution,
    population_yield,
)

MEAN_N = 20.0          # same mean receptor count for both populations
K_INTRA = 0.05
Z = 1.0
DISPERSION_R = 0.5     # small r = strongly overdispersed/zero-inflated
N_MAX = 400

sigma_range = np.logspace(-2, 1.5, 40)  # from 0.01 to about 31.6

poisson_n, poisson_p = poisson_distribution(MEAN_N, N_MAX)
nb_n, nb_p = negative_binomial_distribution(MEAN_N, DISPERSION_R, N_MAX)

# Sanity check: report what fraction of each population has ZERO receptors,
# since that's the "zero-inflation" the handoff specifically flags as a
# real feature of EV data.
frac_zero_poisson = poisson_p[0]
frac_zero_nb = nb_p[0]
print(f"Fraction with 0 receptors -- Poisson: {frac_zero_poisson:.1%}, "
      f"Negative binomial (r={DISPERSION_R}): {frac_zero_nb:.1%}")
print(f"Both distributions have the same mean receptor count: {MEAN_N}")

yields_poisson = [population_yield(poisson_n, poisson_p, s, K_INTRA, Z) for s in sigma_range]
yields_nb = [population_yield(nb_n, nb_p, s, K_INTRA, Z) for s in sigma_range]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left panel: the two receptor-count distributions themselves
ax = axes[0]
plot_n_max = 80
ax.bar(poisson_n[:plot_n_max], poisson_p[:plot_n_max], alpha=0.6, label="Poisson", width=0.9)
ax.bar(nb_n[:plot_n_max], nb_p[:plot_n_max], alpha=0.6, label=f"Neg. binomial (r={DISPERSION_R})", width=0.9)
ax.set_xlabel("Receptor count (n)")
ax.set_ylabel("Probability")
ax.set_title(f"Same mean (n={MEAN_N}), different shapes")
ax.legend()
ax.grid(alpha=0.3)

# Right panel: resulting yield curves
ax = axes[1]
ax.plot(sigma_range, yields_poisson, "-o", markersize=3, label="Poisson population")
ax.plot(sigma_range, yields_nb, "-o", markersize=3, label=f"Neg. binomial population (r={DISPERSION_R})")
ax.set_xscale("log")
ax.set_xlabel("Ligand surface density (sigma)")
ax.set_ylabel("Population-average yield")
ax.set_title("Same mean receptor count -> different yield")
ax.legend()
ax.grid(alpha=0.3)

fig.tight_layout()
save_path = "../figures/population_shape_effect.png"
fig.savefig(save_path, dpi=150)
print(f"\nSaved figure to: {save_path}")

# Quantify the gap at a representative sigma
mid_idx = len(sigma_range) // 2
gap = yields_poisson[mid_idx] - yields_nb[mid_idx]
print(f"\nAt sigma={sigma_range[mid_idx]:.3f}: Poisson yield={yields_poisson[mid_idx]:.4f}, "
      f"Neg. binomial yield={yields_nb[mid_idx]:.4f}, difference={gap:.4f} "
      f"({gap/yields_poisson[mid_idx]:.1%} relative)")
