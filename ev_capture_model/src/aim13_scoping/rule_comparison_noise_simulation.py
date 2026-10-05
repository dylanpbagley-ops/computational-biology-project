"""How often does a truly better rule LOSE to a coarser one on a small labelled panel?
Binormal model: each rule's score ~ N(mu_r * y, 1), noise correlated rho between the two rules.
AUC_r = Phi(mu_r / sqrt 2). Paired DeLong test for the difference."""
import numpy as np
from scipy.stats import norm
rng = np.random.default_rng(20260929)

def auc_and_placements(x, y):
    pos = x[y == 1]; neg = x[y == 0]
    cmp = (pos[:, None] > neg[None, :]) + 0.5 * (pos[:, None] == neg[None, :])
    return cmp.mean(), cmp.mean(1), cmp.mean(0)

def delong_p(sa, sb, y):
    a, va10, va01 = auc_and_placements(sa, y); b, vb10, vb01 = auc_and_placements(sb, y)
    m, n = len(va10), len(va01)
    s10 = np.cov(np.vstack([va10, vb10])); s01 = np.cov(np.vstack([va01, vb01]))
    var = (s10[0, 0] + s10[1, 1] - 2 * s10[0, 1]) / m + (s01[0, 0] + s01[1, 1] - 2 * s01[0, 1]) / n
    z = (a - b) / np.sqrt(var) if var > 0 else 0.0
    return a - b, z

def run(n_safe, n_tox, auc_b, delta, rho, reps=3000):
    mu = lambda A: np.sqrt(2) * norm.ppf(A)
    mb, mt = mu(auc_b), mu(auc_b + delta)
    y = np.r_[np.zeros(n_safe), np.ones(n_tox)].astype(int)
    lose = ties = sig_win = sig_lose = 0
    for _ in range(reps):
        z1 = rng.standard_normal(len(y)); z2 = rho * z1 + np.sqrt(1 - rho ** 2) * rng.standard_normal(len(y))
        sb = mb * y + z1; st = mt * y + z2            # 'true' rule = second
        d, z = delong_p(st, sb, y)
        lose += d < 0; ties += d == 0
        sig_win += z > 1.96; sig_lose += z < -1.96
    return lose / reps, ties / reps, sig_win / reps, sig_lose / reps

print('Panel of n_safe/n_toxic; truer rule has population AUC = 0.90 + delta; coarse rule 0.90; noise correlation rho')
print('%-9s %-6s %-5s | %-24s %-20s' % ('panel', 'delta', 'rho', 'truer rule scores LOWER', 'detected as better (p<.05)'))
for (ns, nt) in ((19, 11), (38, 22), (95, 55), (190, 110), (475, 275)):
    for delta in (0.02, 0.04):
        for rho in (0.8,):
            lo, ti, sw, sl = run(ns, nt, 0.90 - delta / 2 * 0, delta, rho, reps=1500)
            print('%-9s %-6.2f %-5.1f | %-24.0f%% %-20.0f%%' % ('%d/%d' % (ns, nt), delta, rho, 100 * lo, 100 * sw))
E=1
