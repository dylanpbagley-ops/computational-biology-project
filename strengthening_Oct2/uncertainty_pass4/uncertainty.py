"""Exact rank-based diagnostics under hypothetical binary outcome changes."""
import math


def _validate(scores, labels):
    if len(scores) != len(labels) or not scores:
        raise ValueError("Nonempty paired scores and labels required")
    if not all(math.isfinite(v) for v in scores):
        raise ValueError("Finite scores required")
    if not all(isinstance(v, bool) for v in labels):
        raise ValueError("Boolean labels required")
    if not any(labels) or all(labels):
        raise ValueError("Both groups required")


def _ranks(scores):
    return [sum(v < x for v in scores) + (sum(v == x for v in scores) + 1)/2 for x in scores]


def auc(scores, labels):
    """Higher score predicts True; ties contribute one-half."""
    _validate(scores, labels)
    m = sum(labels)
    return (sum(r for r, y in zip(_ranks(scores), labels) if y) - m*(m+1)/2)/(m*(len(labels)-m))


def flip_bounds(scores, labels, k):
    """Sharp extrema after exactly k flips, with explicit attaining witnesses.

    For fixed counts a=True->False and b=False->True, group sizes are fixed.
    Minimizing AUC therefore removes the a largest True ranks and adds the
    b smallest False ranks; maximizing reverses those choices. Enumerate a+b=k.
    This is an adversarial mathematical check, not a clinical relabeling model.
    """
    _validate(scores, labels)
    if not isinstance(k, int) or isinstance(k, bool) or not 0 <= k <= len(labels):
        raise ValueError("k must be an integer in [0,n]")
    ranks = _ranks(scores)
    positive = sorted([i for i, y in enumerate(labels) if y], key=lambda i: ranks[i])
    negative = sorted([i for i, y in enumerate(labels) if not y], key=lambda i: ranks[i])
    baseline_sum = sum(ranks[i] for i in positive)
    result = {}
    for a in range(k+1):
        b = k-a
        if a > len(positive) or b > len(negative):
            continue
        m = len(positive)-a+b
        if m == 0 or m == len(labels):
            continue
        for side, remove_order, add_order in [
            ('min', positive[::-1], negative), ('max', positive, negative[::-1])
        ]:
            remove = remove_order[:a]
            add = add_order[:b]
            rank_sum = baseline_sum - sum(ranks[i] for i in remove) + sum(ranks[i] for i in add)
            value = (rank_sum - m*(m+1)/2)/(m*(len(labels)-m))
            if side not in result or (value < result[side]['auc'] if side == 'min' else value > result[side]['auc']):
                result[side] = {'auc': value, 'flipped_indices': sorted(remove+add), 'favorable_count': m}
    if not result:
        raise ValueError("No assignment retaining both groups")
    return result


def zero_event_upper(n, alpha=.05):
    """Exact one-sided binomial upper bound, conditional on ideal iid detection."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 1 or not 0 < alpha < 1:
        raise ValueError("Positive integer n and alpha in (0,1) required")
    return -math.expm1(math.log(alpha)/n)
