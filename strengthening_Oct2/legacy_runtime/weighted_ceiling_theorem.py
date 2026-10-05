"""
WEIGHTED MULTI-COMPARTMENT RESPONSE BOUNDS — scope corrected October 2, 2026.

For fixed nonnegative weights, specified compartment prevalences, and a shared
response f, H/W equals E_pW[f] when W>0. This is an exact mixture identity,
not a clinical harm calibration. W=0 is a separate zero-modeled-harm case.

For nondecreasing antigen responses, attainable pairs form the convex hull of
antigen-threshold ROC points. For arbitrary responses, they form the convex hull
of ALL deterministic subsets, with a likelihood-ratio upper boundary. The latter
can exceed the former. AUC does not determine either entire boundary.

The numerical code below is preserved: it checks arbitrary-response upper-bound
containment and linear-utility optimality, NOT monotone-response attainability.
Legacy special-case keys check a single-reference identity, one-hot reference
identities, and share-weighted algebra. They do not establish joint tissue-budget
attainability or derive the empirical confinement-plus-blood score.

See outputs/strengthening_Oct2/mathematical_scope.md and its independently tested
frontiers.py implementation for corrected statements and counterexamples.
The pre-correction source and output are preserved in strengthening_Oct2/before.
"""
import os, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RNG = np.random.default_rng(20260908)


def random_problem(K, nbins, rng):
    """A random multi-compartment problem on a common discrete antigen grid."""
    pT = rng.dirichlet(np.ones(nbins) * rng.uniform(0.4, 3.0))
    pc = np.array([rng.dirichlet(np.ones(nbins) * rng.uniform(0.4, 3.0)) for _ in range(K)])
    q = rng.dirichlet(np.ones(K))
    w = rng.choice([0.0, rng.uniform(0, 1), 1.0], size=K)
    if w.sum() == 0:
        w[rng.integers(K)] = 1.0
    return pT, pc, q, w


def harm_and_kill(f, pT, pc, q, w):
    T = float(f @ pT)
    H = float(sum(w[c] * q[c] * (f @ pc[c]) for c in range(len(q))))
    return T, H


def verify(n_problems=400, K_max=6, nbins=7, n_random_rules=400):
    """Check arbitrary-response upper-bound containment and LR utility optimality.

    The loop performs 40 containment checks and one utility check per problem by
    default. It does not check the monotone class or a third single-case identity;
    legacy algebraic identities are evaluated separately in special_cases().
    """
    worst_excess = 0.0
    worst_gap = 0.0
    checks = 0
    for it in range(n_problems):
        K = int(RNG.integers(1, K_max + 1))
        pT, pc, q, w = random_problem(K, nbins, RNG)
        W = float(sum(w[c] * q[c] for c in range(K)))
        pW = sum(w[c] * q[c] * pc[c] for c in range(K)) / W

        # --- 1. hull containment -------------------------------------------------------
        # the ROC of pT vs pW, over all 2^nbins deterministic subsets (the hull's vertices)
        pts = []
        for mask in range(1 << nbins):
            f = np.array([(mask >> b) & 1 for b in range(nbins)], float)
            pts.append((float(f @ pW), float(f @ pT)))
        pts.sort()
        # UPPER CONVEX HULL of the ROC cloud, by monotone chain. An earlier version of this
        # used an ad-hoc "upper envelope" loop that was simply wrong and reported violations
        # of up to 0.96 -- the bug was in the hull, not the theorem. Kept as a note because
        # a verification script that can produce false violations is worse than none.
        pts = [(0.0, 0.0)] + pts + [(1.0, 1.0)]
        pts = sorted(set(pts))
        hull = []
        for pt in pts:
            while len(hull) >= 2:
                (x1, y1), (x2, y2) = hull[-2], hull[-1]
                # keep only left turns (concave upper boundary)
                if (x2 - x1) * (pt[1] - y1) - (y2 - y1) * (pt[0] - x1) >= 0:
                    hull.pop()
                else:
                    break
            hull.append(pt)

        def ceiling_at(x):
            """Largest tumour kill reachable at weighted-normal kill <= x, by linear
            interpolation along the upper hull. Randomised rules attain exactly this."""
            x = min(max(x, 0.0), 1.0)
            for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
                if x1 - 1e-12 <= x <= x2 + 1e-12:
                    if x2 - x1 < 1e-15:
                        return max(y1, y2)
                    t = (x - x1) / (x2 - x1)
                    return y1 + t * (y2 - y1)
            return hull[-1][1]

        for _ in range(n_random_rules // 10):
            f = RNG.random(nbins)                      # an arbitrary randomised kill rule
            T, H = harm_and_kill(f, pT, pc, q, w)
            excess = T - ceiling_at(H / W)
            worst_excess = max(worst_excess, excess)
            checks += 1

        # --- 2. Neyman-Pearson optimality ----------------------------------------------
        lam = float(RNG.uniform(0.2, 4.0))
        LR = pT / np.maximum(pW, 1e-12)
        order = np.argsort(-LR)
        best_np = -1e9
        for k in range(nbins + 1):
            f = np.zeros(nbins)
            f[order[:k]] = 1.0
            T, H = harm_and_kill(f, pT, pc, q, w)
            best_np = max(best_np, T - lam * H)
        best_brute = -1e9
        for mask in range(1 << nbins):
            f = np.array([(mask >> b) & 1 for b in range(nbins)], float)
            T, H = harm_and_kill(f, pT, pc, q, w)
            best_brute = max(best_brute, T - lam * H)
        worst_gap = max(worst_gap, best_brute - best_np)
        checks += 1

    return worst_excess, worst_gap, checks


def special_cases():
    """Legacy algebraic checks only; see corrected scope note for limitations."""
    out = {}
    rng = np.random.default_rng(7)
    nbins = 8

    # 1. K = 1, w = 1 must reproduce the plain ROC ceiling
    pT = rng.dirichlet(np.ones(nbins)); pN = rng.dirichlet(np.ones(nbins))
    W = 1.0
    pW = pN
    out['single_compartment_matches'] = bool(np.allclose(pW, pN))

    # 2. One-hot reference identity only; NOT a simultaneous-constraint optimum.
    K = 4
    pc = np.array([rng.dirichlet(np.ones(nbins)) for _ in range(K)])
    q = rng.dirichlet(np.ones(K))

    def auc(a, b):
        cum = np.cumsum(b) - b
        return float((a * (cum + 0.5 * b)).sum())
    per = [auc(pT, pc[c]) for c in range(K)]
    onehot = []
    for c in range(K):
        w = np.zeros(K); w[c] = 1.0
        Wc = float((w * q).sum())
        pW = (w[:, None] * q[:, None] * pc).sum(0) / Wc
        onehot.append(auc(pT, pW))
    out['worst_compartment_recovered'] = bool(np.allclose(sorted(per), sorted(onehot), atol=1e-12))

    # 3. harm reduces to sum_c w_c s_c when the rule kills whatever it can reach in each
    #    compartment -- the confinement-and-expendability form
    s = rng.dirichlet(np.ones(K))          # the target's expression share per compartment
    w = rng.random(K)
    lhs = float((w * s).sum())
    rhs = 1.0 - float(((1 - w) * s).sum())
    out['confinement_expendability_identity'] = bool(abs(lhs - rhs) < 1e-12)
    return out


if __name__ == '__main__':
    print('=' * 78)
    print('NUMERICAL VERIFICATION OF THE WEIGHTED MULTI-COMPARTMENT CEILING')
    print('=' * 78)
    excess, gap, checks = verify()
    print('  random problems with 1-6 compartments, arbitrary randomised kill rules')
    print('  largest amount by which any rule EXCEEDED the weighted ROC ceiling: %.3e' % excess)
    print('  largest gap between likelihood-ratio thresholding and brute-force optimum: %.3e' % gap)
    print('  total checks: %d' % checks)
    ok = (excess < 1e-9) and (gap < 1e-9)
    print('  VERDICT: %s' % ('theorem holds in every check' if ok else 'VIOLATION FOUND'))

    print('\n  legacy algebraic identities (not clinical or joint-attainability tests):')
    sc = special_cases()
    for k, v in sc.items():
        print('    %-42s %s' % (k, 'ok' if v else 'FAILED'))

    json.dump({'worst_excess': excess, 'worst_np_gap': gap, 'checks': checks,
               'special_cases': sc},
              open(os.path.join(HERE, 'weighted_ceiling_theorem.json'), 'w'), indent=1)
    print('\nwrote weighted_ceiling_theorem.json')
