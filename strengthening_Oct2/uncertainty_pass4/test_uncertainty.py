import itertools
import random
import unittest

from uncertainty import auc, flip_bounds, zero_event_upper


def direct_auc(scores, labels):
    pairs = [(a, b) for a, ya in zip(scores, labels) for b, yb in zip(scores, labels) if ya and not yb]
    if not pairs:
        raise ValueError("Both groups required")
    return sum((a > b) + .5 * (a == b) for a, b in pairs) / len(pairs)


class UncertaintyTests(unittest.TestCase):
    def test_auc_direction_and_ties(self):
        self.assertEqual(auc([1, 2, 3, 4], [False, False, True, True]), 1)
        self.assertEqual(auc([1, 2, 3, 4], [True, True, False, False]), 0)
        self.assertEqual(auc([1, 1, 1], [True, False, False]), .5)
        self.assertAlmostEqual(auc([1, 2, 2, 4], [True, True, False, False]), .125)

    def test_auc_matches_pairwise(self):
        rng = random.Random(412)
        for _ in range(100):
            scores = [rng.randrange(4) for _ in range(8)]
            labels = [True] * 3 + [False] * 5
            rng.shuffle(labels)
            self.assertAlmostEqual(auc(scores, labels), direct_auc(scores, labels))

    def test_exact_flip_bounds_and_witnesses_match_bruteforce(self):
        rng = random.Random(419)
        for _ in range(12):
            scores = [rng.randrange(4) for _ in range(7)]
            labels = [True] * 3 + [False] * 4
            rng.shuffle(labels)
            for k in range(8):
                values = []
                for subset in itertools.combinations(range(7), k):
                    y = [not v if i in subset else v for i, v in enumerate(labels)]
                    if any(y) and not all(y):
                        values.append(direct_auc(scores, y))
                result = flip_bounds(scores, labels, k)
                self.assertAlmostEqual(result['min']['auc'], min(values))
                self.assertAlmostEqual(result['max']['auc'], max(values))
                for side in ['min', 'max']:
                    changed = result[side]['flipped_indices']
                    self.assertEqual(len(changed), k)
                    self.assertEqual(len(set(changed)), k)
                    y = [not v if i in changed else v for i, v in enumerate(labels)]
                    self.assertAlmostEqual(result[side]['auc'], direct_auc(scores, y))

    def test_zero_event_bound_inverts_binomial(self):
        for n in [1, 4, 19, 27, 40, 100]:
            bound = zero_event_upper(n)
            self.assertAlmostEqual((1-bound)**n, .05)
        self.assertGreater(zero_event_upper(4), zero_event_upper(27))

    def test_invalid_inputs(self):
        for scores, labels in [([], []), ([1, 2], [True]), ([1, 2], [True, True]), ([float('nan'), 2], [True, False])]:
            with self.assertRaises(ValueError):
                auc(scores, labels)
        for k in [-1, 3, .5]:
            with self.assertRaises(ValueError):
                flip_bounds([1, 2], [True, False], k)
        for n in [0, -1, .5]:
            with self.assertRaises(ValueError):
                zero_event_upper(n)


if __name__ == '__main__':
    unittest.main()
