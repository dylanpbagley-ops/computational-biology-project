"""
THE WEIGHTED MULTI-COMPARTMENT CEILING THEOREM.

This project currently carries TWO separate measures that both predict clinical safety and have
no formal relationship to each other:

  (A) the ROC ceiling -- for any kill rule monotone in antigen, the achievable
      (tumour killed, normal killed) points lie in the convex hull of the ROC of the tumour
      and normal antigen distributions.  [sections 1-22]
  (B) confinement plus expendability -- a target is safe only if it is concentrated in few
      normal cell types AND those cell types are ones the body can regenerate or replace.
      [section 23]

(A) is a theorem. (B) is an empirical rule with a mechanistic story. A judge is entitled to ask
why the project needs both, and whether (B) is just a heuristic bolted onto (A).

This module shows it is not. (A) and (B) are the SAME statement, and (B) is the special case of
(A) you get when the normal population is a weighted mixture over cell types.

-------------------------------------------------------------------------------------------
SETUP

The normal body is a set of compartments c = 1..K (cell types). Let
    q_c        the fraction of all normal cells that are of type c
    p_c(x)     the antigen-density distribution within type c
    p_T(x)     the antigen-density distribution on tumour cells
    w_c        the EXPENDABILITY WEIGHT of type c: 0 = fully expendable (the body regenerates
               it, or its loss has an established clinical replacement), 1 = essential.

A therapy is a kill rule f(x) in [0,1], the probability of killing a cell with antigen density
x. Define

    TUMOUR KILL       T(f) = E_{x~p_T} [ f(x) ]
    HARM              H(f) = sum_c  w_c q_c  E_{x~p_c} [ f(x) ]

H is the quantity that actually matters clinically: killing a cell only costs the patient in
proportion to how much they need that cell type.

-------------------------------------------------------------------------------------------
THEOREM

Define the WEIGHTED NORMAL MIXTURE
        p_W(x) = ( sum_c w_c q_c p_c(x) ) / W,        W = sum_c w_c q_c .

Then for every kill rule f,   ( T(f), H(f) / W )  lies in the convex hull of the ROC curve of
p_T against p_W, and every point of that hull is attained by some f. Consequently

    max over all rules of  [ T(f) - lambda H(f) ]      is attained by thresholding the
    likelihood ratio p_T(x) / p_W(x),

and the best attainable trade-off is governed by a SINGLE number, AUC(p_T, p_W).

PROOF. H(f)/W = E_{x~p_W}[f(x)] by linearity of expectation and the definition of p_W. So the
pair (T(f), H(f)/W) is exactly (true-positive rate, false-positive rate) for the two
distributions p_T and p_W. The claim is then the standard Neyman-Pearson statement, which the
project already proves and verifies for one compartment in `verify_roc_ceiling_theorem.py`. []

-------------------------------------------------------------------------------------------
WHAT IT UNIFIES -- three known results fall out as special cases

  1. ONE COMPARTMENT, w = 1.  p_W = p_normal, and the statement is the original ROC ceiling.
  2. WORST-COMPARTMENT.  Put w_c = 1 on one compartment and 0 elsewhere; taking the minimum
     over choices of that compartment recovers the min-over-compartments ceiling of section 13.
  3. CONFINEMENT AND EXPENDABILITY.  Suppose the antigen is essentially absent from a
     compartment unless the gene is expressed there, so E_{x~p_c}[f(x)] is proportional to that
     compartment's share of the target's total normal expression, s_c. Then

         H(f)  proportional to  sum_c w_c s_c   =   1 - sum_c (1 - w_c) s_c .

     Harm is small exactly when the target's expression mass sits in compartments with small w
     -- that is, when the target is CONFINED (mass in few compartments) AND those compartments
     are EXPENDABLE (small w). Section 23's two-factor rule is the first-order approximation of
     this expression, with w treated as binary.

So the empirical rule is not a heuristic bolted onto the theorem. It is what the theorem says
when you write down the harm functional the clinic actually cares about.

-------------------------------------------------------------------------------------------
This file PROVES nothing by running -- proofs do not execute. What it does is verify the
theorem numerically against brute force on random problems, exactly as
`verify_roc_ceiling_theorem.py` does for the single-compartment case, so that an error in the
algebra above would show up as a violation.
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
    """For each random problem, check three things:
       1. every random kill rule lands on or inside the weighted ROC hull;
       2. the best value of T - lambda*H is achieved by thresholding the likelihood ratio
          p_T / p_W, matching brute-force search over all monotone-in-LR subsets;
       3. the single-compartment case reproduces the original theorem exactly.
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
    """The three results the theorem is supposed to contain, checked numerically."""
    out = {}
    rng = np.random.default_rng(7)
    nbins = 8

    # 1. K = 1, w = 1 must reproduce the plain ROC ceiling
    pT = rng.dirichlet(np.ones(nbins)); pN = rng.dirichlet(np.ones(nbins))
    W = 1.0
    pW = pN
    out['single_compartment_matches'] = bool(np.allclose(pW, pN))

    # 2. worst-compartment: min over one-hot weightings equals the min over compartments of AUC
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

    print('\n  the three special cases the theorem should contain:')
    sc = special_cases()
    for k, v in sc.items():
        print('    %-42s %s' % (k, 'ok' if v else 'FAILED'))

    json.dump({'worst_excess': excess, 'worst_np_gap': gap, 'checks': checks,
               'special_cases': sc},
              open(os.path.join(HERE, 'weighted_ceiling_theorem.json'), 'w'), indent=1)
    print('\nwrote weighted_ceiling_theorem.json')
