"""Boundary tests for the public GSE305971 threshold-policy pilot."""

from __future__ import annotations

import sys
import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_external_pilot as pilot  # noqa: E402


class ThresholdPolicyTests(unittest.TestCase):
    def test_pooled_and_equal_code_constraints_use_different_weights(self):
        training = {
            "A": np.asarray([0.0, 1.0, 1.0]),
            "B": np.asarray([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        }

        pooled, pooled_rate = pilot.fit_thresholds("pooled", training, 0.2)
        equal_code, equal_rate = pilot.fit_thresholds("equal_code", training, 0.2)

        self.assertEqual(pooled, np.nextafter(0.0, np.inf))
        self.assertEqual(pooled_rate, 0.2)
        self.assertEqual(equal_code, np.nextafter(1.0, np.inf))
        self.assertEqual(equal_rate, 0.0)

    def test_shared_candidate_grid_contains_both_sides_of_each_tie(self):
        training = {
            "A": np.asarray([0.0, 1.0, 2.0]),
            "B": np.asarray([0.0, 1.0, 3.0]),
        }

        candidates = pilot.candidate_thresholds(training)

        self.assertEqual(
            candidates.tolist(),
            [0.0, np.nextafter(0.0, np.inf), 1.0, np.nextafter(1.0, np.inf),
             2.0, np.nextafter(2.0, np.inf), 3.0, np.nextafter(3.0, np.inf)],
        )

    def test_envelope_can_meet_both_constraints_just_above_a_tied_value(self):
        training = {
            "A": np.asarray([0.0, 1.0, 2.0, 3.0]),
            "B": np.asarray([0.0, 1.0, 1.5, 3.0]),
        }

        threshold, constraint_rate = pilot.fit_thresholds("donor_envelope", training, 0.5)

        self.assertEqual(threshold, np.nextafter(1.0, np.inf))
        self.assertLessEqual(np.mean(training["A"] >= threshold), 0.5)
        self.assertLessEqual(np.mean(training["B"] >= threshold), 0.5)
        self.assertEqual(constraint_rate, 0.5)

    def test_all_zero_reference_can_still_detect_a_positive_tumor_value(self):
        training = {"A": np.asarray([0.0, 0.0]), "B": np.asarray([0.0, 0.0])}

        threshold, constraint_rate = pilot.fit_thresholds("donor_envelope", training, 0.1)

        self.assertEqual(threshold, np.nextafter(0.0, np.inf))
        self.assertEqual(pilot.threshold_positive(np.asarray([0.0, 0.0]), threshold), 0)
        self.assertEqual(pilot.threshold_positive(np.asarray([2.0]), threshold), 1)
        self.assertEqual(constraint_rate, 0.0)
        self.assertEqual(pilot.classify_threshold(threshold, max(training["A"])), "ABOVE_MAX")
        self.assertEqual(float(pilot.format_threshold(threshold)), threshold)


class PairedContrastTests(unittest.TestCase):
    @staticmethod
    def rows(values: list[float]) -> list[dict[str, float]]:
        return [
            {
                "TNFRSF17_per_10k": value,
                "GPRC5D_per_10k": value,
                "SLAMF7_per_10k": value,
            }
            for value in values
        ]

    def test_zero_group_mean_keeps_mean_difference_defined(self):
        grouped = {
            ("mm_malignant", "P1"): self.rows([1.0, 1.0]),
            ("mm_normal_pc_prior_art_reference", "P1"): self.rows([0.0, 0.0]),
        }

        results = pilot.paired_contrasts(grouped, {"P1"})
        tnfrsf17 = next(row for row in results if row["marker"] == "TNFRSF17")

        self.assertEqual(tnfrsf17["paired_mean_difference_malignant_minus_Normal_PC"], 1.0)
        self.assertEqual(tnfrsf17["standardized_mean_contrast"], "NR")
        self.assertEqual(tnfrsf17["standardized_contrast_state"], "undefined_constant_Normal_PC_expression")

    def test_empty_group_is_undefined_but_not_zero(self):
        grouped = {("mm_malignant", "P1"): self.rows([0.0, 0.0])}

        results = pilot.paired_contrasts(grouped, {"P1"})
        tnfrsf17 = next(row for row in results if row["marker"] == "TNFRSF17")

        self.assertEqual(tnfrsf17["paired_mean_difference_malignant_minus_Normal_PC"], "NR")
        self.assertEqual(tnfrsf17["standardized_contrast_state"], "missing_or_empty_group")


class CodeEligibilityTests(unittest.TestCase):
    def test_code_metrics_label_50_cells_as_sensitivity(self):
        cell = {"analysis_group": "healthy_normal_pc", "Hospital.Code": "H1", "total_mRNA_molecule_counts": 10}
        for marker in pilot.MARKERS:
            cell[f"{marker}_molecules"] = 0
            cell[f"{marker}_per_10k"] = 0.0

        code_rows, _ = pilot.metric_rows([cell], {"H1"})

        self.assertFalse(code_rows[0]["meets_50_cell_sensitivity"])
        self.assertNotIn("evaluable_50_cell_tail", code_rows[0])

    def test_nonempty_small_codes_remain_in_primary_heldout_results(self):
        healthy_codes = {"H1", "H2", "H3", "H4"}
        mm_codes = {"M1", "M2", "M3", "M4"}
        grouped: dict[tuple[str, str], list[dict[str, float]]] = {}
        for code in healthy_codes:
            grouped[("healthy_normal_pc", code)] = PairedContrastTests.rows([0.0, 0.0])
        for code in mm_codes:
            grouped[("mm_malignant", code)] = PairedContrastTests.rows([0.0, 2.0])

        result = pilot.evaluate_policy_set("all_nonempty_primary", healthy_codes, mm_codes, grouped)
        code_rows = result["cv_code_rows"]
        heldout_healthy = [
            row for row in code_rows
            if row["marker"] == "TNFRSF17" and row["budget"] == 0.1
            and row["policy"] == "pooled" and row["holdout_group"] == "heldout_healthy_normal_pc"
        ]

        self.assertEqual({row["Hospital.Code"] for row in heldout_healthy}, healthy_codes)
        self.assertTrue(all(row["n_cells"] == 2 for row in heldout_healthy))

    def test_heldout_healthy_and_mm_perturbations_do_not_change_fitted_thresholds(self):
        healthy_codes = {"H1", "H2", "H3", "H4"}
        mm_codes = {"M1", "M2", "M3", "M4"}
        grouped: dict[tuple[str, str], list[dict[str, float]]] = {}
        for code in healthy_codes:
            grouped[("healthy_normal_pc", code)] = PairedContrastTests.rows([0.0, 0.0])
        for code in mm_codes:
            grouped[("mm_malignant", code)] = PairedContrastTests.rows([0.0, 2.0])
        baseline = pilot.evaluate_policy_set("all_nonempty_primary", healthy_codes, mm_codes, grouped)
        heldout_fold = baseline["folds"]["healthy"]["H1"]

        perturbed = {key: [dict(row) for row in values] for key, values in grouped.items()}
        heldout_codes = {
            ("healthy_normal_pc", code)
            for code in healthy_codes
            if baseline["folds"]["healthy"][code] == heldout_fold
        } | {
            ("mm_malignant", code)
            for code in mm_codes
            if baseline["folds"]["mm"][code] == heldout_fold
        }
        for key in heldout_codes:
            for index, row in enumerate(perturbed[key]):
                row["TNFRSF17_per_10k"] = 1_000_000.0 + index
        changed = pilot.evaluate_policy_set("all_nonempty_primary", healthy_codes, mm_codes, perturbed)

        def thresholds_at(result: dict[str, object]) -> dict[tuple[float, str, str], float]:
            return {
                (row["budget"], row["policy"], row["heldout_group"]): float(row["threshold_per_10k"])
                for row in result["cv_fold_rows"]
                if row["analysis_set"] == "all_nonempty_primary" and row["marker"] == "TNFRSF17"
                and row["fold"] == heldout_fold
            }

        self.assertEqual(thresholds_at(baseline), thresholds_at(changed))
        baseline_positive_counts = {
            (row["budget"], row["policy"], row["holdout_group"], row["Hospital.Code"]): row["positive_cells"]
            for row in baseline["cv_code_rows"]
            if row["analysis_set"] == "all_nonempty_primary" and row["marker"] == "TNFRSF17"
            and row["fold"] == heldout_fold
        }
        changed_positive_counts = {
            (row["budget"], row["policy"], row["holdout_group"], row["Hospital.Code"]): row["positive_cells"]
            for row in changed["cv_code_rows"]
            if row["analysis_set"] == "all_nonempty_primary" and row["marker"] == "TNFRSF17"
            and row["fold"] == heldout_fold
        }
        self.assertNotEqual(baseline_positive_counts, changed_positive_counts)


class CountParserTests(unittest.TestCase):
    @staticmethod
    def write_table(folder: Path, other_row: str) -> Path:
        path = folder / "counts.tsv.gz"
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            handle.write("cell-a\tcell-b\n")
            handle.write(other_row + "\n")
            for marker in pilot.MARKERS:
                handle.write(f"{marker}\t0\t0\n")
        return path

    def test_parser_refuses_fractional_count_token(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.write_table(Path(temp), "OTHER\t1\t2.5")

            with self.assertRaises(ValueError):
                pilot.parse_count_file(path, set())

    def test_parser_refuses_trailing_malformed_count_token(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.write_table(Path(temp), "OTHER\t1\t2x")

            with self.assertRaises(ValueError):
                pilot.parse_count_file(path, set())

    def test_parser_refuses_integer_sum_overflow(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.write_table(Path(temp), "OTHER\t9223372036854775807\t0")
            with gzip.open(path, "at", encoding="utf-8", newline="") as handle:
                handle.write("OTHER2\t1\t0\n")

            with self.assertRaisesRegex(ValueError, "overflow"):
                pilot.parse_count_file(path, set())

    def test_exact_cell_join_follows_ids_when_count_columns_are_reordered(self):
        annotations = {
            "cell-a": {"Amp.Batch.ID": "batch-1", "Hospital.Code": "code-a"},
            "cell-b": {"Amp.Batch.ID": "batch-1", "Hospital.Code": "code-b"},
            "cell-c": {"Amp.Batch.ID": "batch-2", "Hospital.Code": "code-c"},
        }
        ids = ["cell-b", "not-annotated", "cell-a", "cell-c"]

        matches = pilot.exact_batch_matches(ids, annotations, "batch-1")

        self.assertEqual([(pos, cell_id, row["Hospital.Code"]) for pos, cell_id, row in matches], [
            (0, "cell-b", "code-b"),
            (2, "cell-a", "code-a"),
        ])


if __name__ == "__main__":
    unittest.main()
