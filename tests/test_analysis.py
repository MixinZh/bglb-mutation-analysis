import copy
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

from bglb_analysis.analysis import compare, load_policy, metrics, read_csv, run, select


ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "examples/synthetic"


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.policy = load_policy(ROOT / "selection_policy.json")
        self.consensus = read_csv(DEMO / "consensus.csv", [])
        self.measurements = read_csv(DEMO / "measurements.csv", [])
        self.predictions = read_csv(DEMO / "predictions.csv", [])

    def selection(self):
        return select(self.consensus, self.measurements, self.policy)

    def test_four_records_and_existing_quality_both_required(self):
        ledger = self.selection()
        self.assertEqual([r['mutation_id'] for r in ledger if r['primary']], ['A10V', 'G20S'])
        self.assertEqual(sum(r['all_repeated'] for r in ledger), 3)
        self.assertIn('fewer_than_4', ledger[3]['exclusion_reasons'])

    def test_changing_predictions_cannot_change_selection(self):
        ledger = self.selection()
        before = copy.deepcopy(ledger)
        for row in self.predictions:
            row['prediction'] = '2.0'
        compare(ledger, self.consensus, self.measurements, self.predictions, 'demo_model', self.policy)
        self.assertEqual(ledger, before)

    def test_duplicate_record_ids_are_rejected(self):
        self.measurements.append(self.measurements[0])
        with self.assertRaisesRegex(ValueError, 'duplicate record_id'):
            self.selection()

    def test_claimed_four_records_requires_four_source_ids(self):
        self.measurements.pop(0)
        with self.assertRaisesRegex(ValueError, 'count does not reconcile'):
            self.selection()

    def test_source_id_mismatch_is_rejected(self):
        self.consensus[0]['target_source_record_ids'] = 'x|demo02|demo03|demo04'
        with self.assertRaisesRegex(ValueError, 'Source IDs'):
            self.selection()

    def test_opposing_effects_cannot_be_mislabeled_core(self):
        self.consensus[2]['quality_tier'] = 'CORE'
        self.consensus[2]['quality_reason_codes'] = ''
        selected = next(r for r in self.selection() if r['mutation_id']=='L30A')
        self.assertFalse(selected['primary'])
        self.assertIn('opposing', selected['exclusion_reasons'])

    def test_core_does_not_allow_global_wt(self):
        self.measurements[0]['baseline_type'] = 'global_wt'
        selected = next(r for r in self.selection() if r['mutation_id']=='A10V')
        self.assertEqual(selected['retained_records'], 3)
        self.assertFalse(selected['primary'])

    def test_missing_prediction_stops_instead_of_changing_subset(self):
        with self.assertRaisesRegex(ValueError, 'Missing prediction'):
            compare(self.selection(), self.consensus, self.measurements, self.predictions[1:], 'demo_model', self.policy)

    def test_duplicate_selected_predictions_stop(self):
        self.predictions.append(self.predictions[0])
        with self.assertRaisesRegex(ValueError, 'duplicate mutation_id'):
            compare(self.selection(), self.consensus, self.measurements, self.predictions, 'demo_model', self.policy)

    def test_mismatched_targets_stop(self):
        self.predictions[0]['observed_delta_log10_eff'] = '0.5'
        with self.assertRaisesRegex(ValueError, 'different measured target'):
            compare(self.selection(), self.consensus, self.measurements, self.predictions, 'demo_model', self.policy)

    def test_nonfinite_values_stop(self):
        self.measurements[0]['delta_log10_eff'] = 'nan'
        with self.assertRaisesRegex(ValueError, 'finite'):
            self.selection()

    def test_mutation_position_mismatch_stops(self):
        self.consensus[0]['position'] = '99'
        with self.assertRaisesRegex(ValueError, 'position disagree'):
            self.selection()

    def test_high_spread_cannot_be_mislabeled_core(self):
        self.measurements[0]['delta_log10_eff'] = '-3.0'
        selected = next(r for r in self.selection() if r['mutation_id']=='A10V')
        self.assertFalse(selected['primary'])
        self.assertIn('spread', selected['exclusion_reasons'])

    def test_one_excluded_source_record_does_not_veto_valid_repeats(self):
        self.consensus[1].update(quality_tier='EXTENDED', quality_reason_codes='EXPRESSION_CONTRADICTION',
                                 n_raw_rows='5', quarantined_record_ids='unusable01',
                                 source_record_ids='demo05|demo06|demo07|demo08|unusable01')
        selected = next(r for r in self.selection() if r['mutation_id']=='G20S')
        self.assertTrue(selected['primary'])
        self.assertEqual(selected['retained_records'], 4)
        self.assertEqual(selected['total_database_records'], 5)
        self.assertEqual(selected['earlier_excluded_record_ids'], 'unusable01')

    def test_individual_row_flag_excludes_only_that_row(self):
        self.measurements[0]['phase1_downweight_reason_codes'] = 'EFFICIENCY_COMPONENT_CONFLICT'
        selected = next(r for r in self.selection() if r['mutation_id']=='A10V')
        self.assertEqual(selected['retained_records'], 3)
        self.assertEqual(selected['newly_excluded_usable_record_ids'], 'demo01')

    def test_rank_by_retained_records_not_mutation_position(self):
        new = dict(self.measurements[4], record_id='extra01', context_key='extra_group', delta_log10_eff='0.3')
        self.measurements.append(new)
        self.consensus[1]['n_rows']='5'
        self.consensus[1]['n_verified_contexts']='5'
        self.consensus[1]['target_source_record_ids']+='|extra01'
        self.assertEqual(self.selection()[0]['mutation_id'], 'G20S')

    def test_constant_predictions_have_no_rank_correlation(self):
        result = metrics([{'observed_log10': -1., 'predicted_log10': 0.},
                          {'observed_log10': 1., 'predicted_log10': 0.}], 0.2)
        self.assertIsNone(result['spearman'])
        self.assertEqual(result['direction_correct'], 0)

    def test_no_beneficial_observations_are_counted_as_zero_examples(self):
        result = metrics([{'observed_log10': -1., 'predicted_log10': -0.1}], 0.2)
        self.assertEqual(result['beneficial_measured'], 0)
        self.assertIsNone(result['spearman'])

    def test_context_first_summary_does_not_overweight_same_group(self):
        rows = self.measurements[:4]
        for row, value in zip(rows, [-1.0, -1.0, -1.0, -0.5]):
            row['delta_log10_eff'] = str(value)
        for row in rows[:3]:
            row['context_key'] = 'same_group'
        self.consensus = self.consensus[:1]
        self.consensus[0]['n_verified_contexts'] = '2'
        self.consensus[0]['target_delta_log10_eff'] = '-0.75'
        self.measurements = rows
        selected = self.selection()[0]
        self.assertFalse(selected['primary'])
        self.assertEqual(selected['retained_records'], 4)
        self.assertEqual(selected['retained_observed_log10'], -.75)

    def test_end_to_end_and_existing_output_preservation(self):
        with tempfile.TemporaryDirectory() as scratch:
            args = Namespace(consensus=DEMO/'consensus.csv', measurements=DEMO/'measurements.csv',
                             predictions=DEMO/'predictions.csv', policy=ROOT/'selection_policy.json',
                             out=Path(scratch)/'result', model='demo_model', dataset_label='SYNTHETIC', plot=False)
            summary = run(args)
            self.assertEqual(summary['selection']['primary'], 2)
            self.assertEqual(summary['metrics']['primary']['direction_correct'], 0)
            self.assertEqual(summary['metrics']['primary']['spearman'], -1)
            self.assertEqual(json.loads((args.out/'summary.json').read_text()), summary)
            before = (args.out/'summary.json').read_bytes()
            with self.assertRaisesRegex(ValueError, 'existing results'):
                run(args)
            self.assertEqual((args.out/'summary.json').read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
