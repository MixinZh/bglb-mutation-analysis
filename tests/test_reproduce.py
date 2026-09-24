"""Check real results against the verified comparison made before packaging."""

import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from bglb_analysis.analysis import read_csv, write_csv
from bglb_analysis.reproduce import reproduce, validate


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/bglb"


class RealDataTests(unittest.TestCase):
    def test_real_data_reproduces_previous_results_and_accounts_for_every_record(self):
        expected = json.loads((ROOT / "tests/expected_results.json").read_text())
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "results"
            result = reproduce(DATA, ROOT / "selection_policy.json", output)
            self.assertEqual(result["selection"], expected["selection"])
            for subset, scores in expected["metrics"].items():
                for metric, value in scores.items():
                    actual = result["metrics"][subset][metric]
                    if isinstance(value, float):
                        self.assertAlmostEqual(actual, value, places=12)
                    else:
                        self.assertEqual(actual, value)
            comparison = read_csv(output / "comparison.csv", [])
            self.assertEqual(len(comparison), len(expected["mutations"]))
            numeric = {"observed_log10", "predicted_log10", "observed_fold_vs_wt", "predicted_fold_vs_wt"}
            for actual, previous in zip(comparison, expected["mutations"]):
                for key, value in previous.items():
                    if key in numeric:
                        self.assertAlmostEqual(float(actual[key]), float(value), places=12)
                    else:
                        self.assertEqual(actual[key], value)
            records = read_csv(output / "record_selection.csv", [])
            self.assertEqual(len(records), 1192)
            self.assertEqual(sum(r["included_in_primary_comparison"] == "True" for r in records), 63)
            self.assertEqual(sum(bool(r["original_exclusion_reasons"]) for r in records), 224)
            self.assertFalse(any(r["row_exclusion_reason"] == "excluded_in_earlier_raw_data_audit" for r in records))
            n220f = next(r for r in comparison if r["mutation_id"] == "N220F")
            self.assertEqual(n220f["direct_wt_records"], "10")
            self.assertGreater(float(n220f["record_min_fold_vs_wt"]), 1)
            self.assertLess(float(n220f["predicted_fold_vs_wt"]), 1)

    def changed_input(self, name, change, update_hash=True):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        data = Path(folder.name) / "data"
        shutil.copytree(DATA, data)
        path = data / name
        rows = read_csv(path, [])
        change(rows)
        write_csv(path, rows)
        if update_hash:
            manifest_path = data / "manifest.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["files"][name]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
        return data

    def test_changed_input_fails_before_calculation(self):
        data = self.changed_input("measurements.csv", lambda rows: rows[0].update(delta_log10_eff="9"), False)
        with self.assertRaisesRegex(ValueError, "Input hash changed"):
            validate(data)

    def test_effect_must_follow_measured_kinetics_even_if_manifest_updated(self):
        data = self.changed_input("measurements.csv", lambda rows: rows[0].update(delta_log10_eff="9"))
        with self.assertRaisesRegex(ValueError, "effect for record"):
            validate(data)

    def test_wrong_wt_record_cannot_be_used(self):
        data = self.changed_input("measurements.csv", lambda rows: rows[0].update(baseline_source_record_ids="1940"))
        with self.assertRaisesRegex(ValueError, "WT source mismatch"):
            validate(data)

    def test_catpred_wt_normalization_is_checked(self):
        data = self.changed_input("predictions.csv", lambda rows: rows[0].update(catpred_wt_efficiency_ratio="1"))
        with self.assertRaisesRegex(ValueError, "CatPred normalization"):
            validate(data)

    def test_an_earlier_exclusion_requires_a_reason(self):
        def remove_reason(rows):
            next(r for r in rows if r["row_audit_status"] == "QUARANTINE")["quarantine_reason_codes"] = ""
        data = self.changed_input("audit.csv", remove_reason)
        with self.assertRaisesRegex(ValueError, "Missing exclusion reason"):
            validate(data)

    def test_known_groups_are_preserved_without_names(self):
        for row in read_csv(DATA / "audit.csv", []):
            self.assertRegex(row["context_key"], r"^(group_\d{4}|UNKNOWN)$")
        for row in read_csv(DATA / "measurements.csv", []):
            self.assertRegex(row["context_key"], r"^(group_\d{4}|UNVERIFIED::\d+)$")


if __name__ == "__main__":
    unittest.main()
