import gzip
import hashlib
import importlib.util
import json
import pathlib
import tempfile
import unittest


SOURCE = pathlib.Path(__file__).resolve().parents[1] / "metadata_diagnostics.py"


def api(name):
    if not SOURCE.is_file():
        raise AssertionError(f"metadata diagnostic API {name} is missing")
    spec = importlib.util.spec_from_file_location("metadata_diagnostics_under_test", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, name):
        raise AssertionError(f"metadata diagnostic API {name} is missing")
    return getattr(module, name)


class MetadataDiagnosticTests(unittest.TestCase):
    def test_cell_metadata_joins_by_key_and_reorders_explicitly(self):
        join = api("join_metadata_by_cell_id")
        rows = [{"index": "cell-b", "sample_ID": "NBM-1"}, {"index": "cell-a", "sample_ID": "NBM-2"}]
        result = join(["cell-a", "cell-b"], rows)
        self.assertEqual([row["sample_ID"] for row in result], ["NBM-2", "NBM-1"])

    def test_duplicate_metadata_cell_id_is_rejected_even_when_conflicting(self):
        join = api("join_metadata_by_cell_id")
        rows = [{"index": "cell-a", "sample_ID": "NBM-1"}, {"index": "cell-a", "sample_ID": "MM-1"}]
        with self.assertRaises(ValueError):
            join(["cell-a"], rows)

    def test_unknown_cell_id_is_rejected(self):
        join = api("join_metadata_by_cell_id")
        with self.assertRaises(ValueError):
            join(["not-in-metadata"], [{"index": "cell-a", "sample_ID": "NBM-1"}])

    def test_absent_covariate_remains_explicitly_unknown(self):
        field = api("source_field_or_unknown")
        self.assertEqual(field({}, "batch"), "NR")
        self.assertEqual(field({"batch": ""}, "batch"), "NR")
        self.assertEqual(field({"batch": "Batch-2"}, "batch"), "Batch-2")

    def test_unknown_selected_sample_is_rejected(self):
        select = api("eligible_sample_groups")
        rows = [
            {"index": "a", "sample_ID": "NBM-1", "disease_stage": "NBM", "normal_or_neoplastic": "normal"},
            {"index": "b", "sample_ID": "MM-999", "disease_stage": "MM", "normal_or_neoplastic": "neoplastic"},
        ]
        with self.assertRaises(ValueError):
            select(rows)

    def test_grouping_keeps_all_cells_from_each_donor_together(self):
        select = api("eligible_sample_groups")
        rows = [
            {"index": f"n{i}", "sample_ID": "NBM-1", "disease_stage": "NBM", "normal_or_neoplastic": "normal"}
            for i in range(3)
        ] + [
            {"index": f"m{i}", "sample_ID": "MM-1", "disease_stage": "MM", "normal_or_neoplastic": "neoplastic"}
            for i in range(2)
        ]
        groups = select(rows)
        self.assertEqual(groups["NBM"]["NBM-1"], [0, 1, 2])
        self.assertEqual(groups["MM"]["MM-1"], [3, 4])

    def test_count_validation_rejects_fractional_or_impossible_counts(self):
        validate = api("validate_umi_pair")
        with self.assertRaises(ValueError):
            validate(1.5, 10)
        with self.assertRaises(ValueError):
            validate(11, 10)
        self.assertEqual(validate("3.0", "10.0"), (3, 10))

    def test_hypergeometric_thinning_respects_boundaries_and_depth(self):
        thin = api("thin_marker_count")
        import numpy as np
        rng = np.random.default_rng(123)
        self.assertEqual(thin(0, 8, 4, rng), 0)
        self.assertEqual(thin(8, 8, 4, rng), 4)
        self.assertEqual(thin(3, 8, 8, rng), 3)
        self.assertGreaterEqual(thin(4, 8, 4, rng), 0)
        self.assertLessEqual(thin(4, 8, 4, rng), 4)

    def test_common_depth_uses_only_donor_depth_summaries(self):
        choose = api("choose_common_depth")
        self.assertEqual(choose([100, 300, 500, 700, 900]), 180)

    def test_thinning_replicates_are_reproducible_from_fixed_seed(self):
        import numpy as np
        thin = api("thin_marker_count")
        def sequence():
            rng = np.random.default_rng(20261003)
            return [thin(4, 10, 5, rng) for _ in range(20)]
        self.assertEqual(sequence(), sequence())

    def test_depth_figure_is_written_under_next_stage_delivery(self):
        figure_path = api("depth_tail_figure_path")
        self.assertEqual(
            figure_path(),
            SOURCE.resolve().parents[1] / "delivery" / "figures" / "metadata_depth_tail.svg",
        )

    def test_append_only_status_update_uses_unchanged_backed_up_prefix(self):
        verify = api("verify_input_snapshot")
        with tempfile.TemporaryDirectory() as temp:
            temp_path = pathlib.Path(temp)
            root = temp_path / "root"
            next_stage = temp_path / "next_stage"
            status = root / "strengthening_execution_Oct2" / "STATUS.md"
            backup = next_stage / "before" / "execution_STATUS.md"
            status.parent.mkdir(parents=True)
            backup.parent.mkdir(parents=True)
            original = b"# Existing status\r\nold history\r\n"
            backup.write_bytes(original)
            status.write_bytes(original + b"\r\n## New pointer\r\nAdded at end.\r\n")
            manifest_path = next_stage / "input_hashes.json"
            manifest_path.write_text(json.dumps({
                "protected_inputs": {},
                "current_execution_snapshot": {
                    "strengthening_execution_Oct2/STATUS.md": hashlib.sha256(original).hexdigest()
                },
            }), encoding="utf-8")
            self.assertEqual(verify(manifest_path, root, next_stage), 0)

            status.write_bytes(b"# Altered status\r\n" + original)
            with self.assertRaisesRegex(ValueError, "changed starting-state file"):
                verify(manifest_path, root, next_stage)

    def test_root_project_status_append_uses_verified_backup_and_rejects_edited_prefix(self):
        verify = api("verify_input_snapshot")
        with tempfile.TemporaryDirectory() as temp:
            temp_path = pathlib.Path(temp)
            root = temp_path / "root"
            next_stage = temp_path / "next_stage"
            status = root / "PROJECT_STATUS_Oct2.md"
            backup = next_stage / "before" / "PROJECT_STATUS_Oct2.md"
            root.mkdir(parents=True)
            backup.parent.mkdir(parents=True)
            original = b"# Project status\r\nprotected history\r\n"
            backup.write_bytes(original)
            status.write_bytes(original + b"\r\n## Next-stage pointer\r\nSee strengthening_execution_Oct2/next_stage.\r\n")
            manifest_path = next_stage / "input_hashes.json"
            manifest_path.write_text(json.dumps({
                "protected_inputs": {},
                "current_execution_snapshot": {
                    "PROJECT_STATUS_Oct2.md": hashlib.sha256(original).hexdigest()
                },
            }), encoding="utf-8")
            self.assertEqual(verify(manifest_path, root, next_stage), 0)

            status.write_bytes(b"# Rewritten project status\r\n" + original)
            with self.assertRaisesRegex(ValueError, "changed starting-state file: PROJECT_STATUS_Oct2.md"):
                verify(manifest_path, root, next_stage)

    def test_spearman_matches_independent_rank_and_tie_example(self):
        spearman = api("spearman_rho")
        # Average ranks are [1, 2.5, 2.5, 4] and [4, 1.5, 1.5, 3],
        # whose Pearson correlation is exactly -1/3.
        self.assertAlmostEqual(spearman([1, 2, 2, 4], [4, 1, 1, 2]), -1 / 3, places=12)

    def test_thinning_refuses_unverified_total_counts(self):
        verify = api("verify_matrix_column_totals")
        with tempfile.TemporaryDirectory() as temp:
            path = pathlib.Path(temp) / "counts.csv.gz"
            with gzip.open(path, "wt", newline="") as f:
                f.write(",cell-a,cell-b\nGENE1,1,0\nGENE2,1,3\n")
            self.assertEqual(verify(path, ["cell-a", "cell-b"], [2, 3]), [2, 3])
            with self.assertRaises(ValueError):
                verify(path, ["cell-a", "cell-b"], [2, 4])

    def test_matrix_row_parser_vectorizes_and_rejects_invalid_umi_counts(self):
        parse_row = api("parse_matrix_counts_row")
        self.assertEqual(parse_row(["1.0", "0", 3], 2).tolist(), [1, 0, 3])
        for invalid in (["1", "1.5"], ["1", "-1"], ["1", "inf"], ["1", "bad"]):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                parse_row(invalid, 2)

    def test_duplicate_matrix_cells_cannot_pass_total_verification(self):
        verify = api("verify_matrix_column_totals")
        with tempfile.TemporaryDirectory() as temp:
            path = pathlib.Path(temp) / "counts.csv.gz"
            with gzip.open(path, "wt", newline="") as f:
                f.write(",cell-a,cell-a\nGENE1,1,1\n")
            with self.assertRaises(ValueError):
                verify(path, ["cell-a"], [1])


if __name__ == "__main__":
    unittest.main()
