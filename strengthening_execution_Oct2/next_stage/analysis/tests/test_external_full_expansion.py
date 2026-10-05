"""Tests for the exposure-aware fixed-code expansion analysis."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_external_full_expansion as full  # noqa: E402
import run_external_pilot as pilot  # noqa: E402


class FixedTrainingExpansionTests(unittest.TestCase):
    @staticmethod
    def rows(n: int, values: dict[str, float]) -> list[dict[str, float]]:
        result = []
        for _ in range(n):
            row: dict[str, float] = {}
            for marker in pilot.MARKERS:
                row[f"{marker}_per_10k"] = values.get(marker, 0.0)
            result.append(row)
        return result

    @staticmethod
    def grouped() -> dict[tuple[str, str], list[dict[str, float]]]:
        grouped = {}
        for code, values in {
            "H1": [0.0, 0.0], "H2": [0.0, 0.5], "H3": [0.0, 1.0], "H4": [0.0, 1.5],
        }.items():
            grouped[("healthy_normal_pc", code)] = [
                FixedTrainingExpansionTests.rows(1, {"TNFRSF17": value})[0] for value in values
            ]
        for code in ("H5", "H6"):
            grouped[("healthy_normal_pc", code)] = FixedTrainingExpansionTests.rows(2, {"TNFRSF17": 0.0})
        for code in ("M5", "M6"):
            grouped[("mm_malignant", code)] = FixedTrainingExpansionTests.rows(2, {"TNFRSF17": 2.0})
        grouped[("healthy_ambiguous_malignant_label", "H5")] = FixedTrainingExpansionTests.rows(
            1, {"TNFRSF17": 0.5}
        )
        grouped[("healthy_ambiguous_malignant_label", "H1")] = FixedTrainingExpansionTests.rows(
            2, {"TNFRSF17": 5.0}
        )
        return grouped

    def test_heldout_marker_changes_do_not_change_any_fitted_threshold(self):
        grouped = self.grouped()
        baseline = full.evaluate_fixed_training({"H1", "H2", "H3", "H4"}, {"H5", "H6"}, {"M5", "M6"}, grouped)
        perturbed = {key: [dict(row) for row in rows] for key, rows in grouped.items()}
        for key in (("healthy_normal_pc", "H5"), ("healthy_normal_pc", "H6"), ("mm_malignant", "M5"), ("mm_malignant", "M6")):
            for index, row in enumerate(perturbed[key]):
                row["TNFRSF17_per_10k"] = 1_000_000.0 + index
        changed = full.evaluate_fixed_training({"H1", "H2", "H3", "H4"}, {"H5", "H6"}, {"M5", "M6"}, perturbed)
        threshold = lambda result: {
            (row["marker"], row["budget"], row["policy"], row["fit_reference"]): float(row["threshold_per_10k"])
            for row in result["summary_rows"]
        }
        self.assertEqual(threshold(baseline), threshold(changed))
        self.assertNotEqual(baseline["code_rows"], changed["code_rows"])

    def test_each_marker_uses_its_own_training_values(self):
        grouped = self.grouped()
        for code in ("H1", "H2", "H3", "H4"):
            grouped[("healthy_normal_pc", code)] = [
                {"TNFRSF17_per_10k": 0.0, "GPRC5D_per_10k": 1.0, "SLAMF7_per_10k": 2.0},
                {"TNFRSF17_per_10k": 0.0, "GPRC5D_per_10k": 1.0, "SLAMF7_per_10k": 2.0},
            ]
        result = full.evaluate_fixed_training({"H1", "H2", "H3", "H4"}, {"H5", "H6"}, {"M5", "M6"}, grouped)
        thresholds = {
            row["marker"]: float(row["threshold_per_10k"])
            for row in result["summary_rows"]
            if row["budget"] == 0.1 and row["policy"] == "pooled" and row["heldout_group"] == "heldout_new_healthy_normal_pc"
        }
        self.assertEqual(thresholds["TNFRSF17"], float(np.nextafter(0.0, np.inf)))
        self.assertEqual(thresholds["GPRC5D"], float(np.nextafter(1.0, np.inf)))
        self.assertEqual(thresholds["SLAMF7"], float(np.nextafter(2.0, np.inf)))

    def test_fixed_training_refuses_code_overlap(self):
        with self.assertRaisesRegex(ValueError, "disjoint"):
            full.evaluate_fixed_training({"H1", "H2"}, {"H1", "H3"}, {"M5"}, self.grouped())

    def test_empty_heldout_code_is_not_counted_as_a_zero_rate(self):
        grouped = self.grouped()
        grouped[("healthy_normal_pc", "H6")] = []
        result = full.evaluate_fixed_training({"H1", "H2", "H3", "H4"}, {"H5", "H6"}, {"M5", "M6"}, grouped)
        rows = [row for row in result["summary_rows"] if row["marker"] == "TNFRSF17" and row["budget"] == 0.1
                and row["policy"] == "pooled" and row["heldout_group"] == "heldout_new_healthy_normal_pc"]
        self.assertEqual(rows[0]["heldout_codes"], 1)
        self.assertEqual(rows[0]["heldout_cells"], 2)
        self.assertEqual({row["Hospital.Code"] for row in result["code_rows"] if row["heldout_group"] == "heldout_new_healthy_normal_pc"}, {"H5"})

    def test_uncertain_reference_sensitivity_is_separately_labeled(self):
        grouped = self.grouped()
        baseline = full.evaluate_fixed_training({"H1", "H2", "H3", "H4"}, {"H5", "H6"}, {"M5", "M6"}, grouped)
        augmented = full.evaluate_fixed_training(
            {"H1", "H2", "H3", "H4"}, {"H5", "H6"}, {"M5", "M6"}, grouped,
            include_uncertain_reference=True,
        )
        baseline_threshold = next(row["threshold_per_10k"] for row in baseline["summary_rows"]
                                  if row["marker"] == "TNFRSF17" and row["budget"] == 0.1 and row["policy"] == "pooled")
        augmented_threshold = next(row["threshold_per_10k"] for row in augmented["summary_rows"]
                                   if row["marker"] == "TNFRSF17" and row["budget"] == 0.1 and row["policy"] == "pooled")
        self.assertNotEqual(baseline_threshold, augmented_threshold)
        self.assertEqual({row["fit_reference"] for row in augmented["summary_rows"]}, {"healthy_normal_pc_plus_uncertain_label_cells"})

    def test_locked_split_has_no_pilot_mm_overlap_in_primary_evaluation(self):
        lock_path = Path(__file__).resolve().parents[1] / "results" / "external_full_selection_lock.json"
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
        universe = lock["code_universe"]
        self.assertEqual(len(universe["pilot_healthy_training_codes_seen_in_pilot"]), 8)
        self.assertEqual(len(universe["new_healthy_evaluation_codes_not_seen_in_pilot"]), 9)
        self.assertEqual(len(universe["pilot_mm_codes_excluded_from_confirmatory_evaluation"]), 16)
        self.assertEqual(len(universe["new_mm_evaluation_codes_not_seen_in_pilot"]), 90)
        self.assertFalse(set(universe["pilot_mm_codes_excluded_from_confirmatory_evaluation"]) &
                         set(universe["new_mm_evaluation_codes_not_seen_in_pilot"]))
        self.assertNotIn("tlv53", universe["secondary_full_code_cv_mm_codes_excluding_preexposed"])


if __name__ == "__main__":
    unittest.main()
