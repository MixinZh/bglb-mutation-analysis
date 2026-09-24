import copy
import unittest

from bglb_analysis.d2d_source import AUDIT_KINETICS, FIELDS, check_audit, check_records, export_row


class PublicSourceTests(unittest.TestCase):
    def setUp(self):
        self.public = [{key: None for key in FIELDS.values()}]
        self.public[0].update(id=1, resid='A', resnum=10, resmut='V', expressed=True,
                              curated=True, submitted_for_curation=True,
                              KM_avg=2.0, kcat_avg=4.0, kcat_over_KM=2.0)
        self.local = [export_row(self.public[0])]

    def test_equivalent_decimal_formats_match_without_tolerance(self):
        self.local[0]['KM (mM)'] = '2'
        result = check_records(self.local, self.public)
        self.assertEqual(result['local_records_matching_all_columns'], 1)
        self.public[0]['KM_avg'] = 2.000001
        result = check_records(self.local, self.public)
        self.assertEqual(result['local_records_matching_kinetic_columns'], 0)

    def test_curation_change_is_distinct_from_a_kinetic_change(self):
        self.public[0]['curated'] = False
        result = check_records(self.local, self.public)
        self.assertEqual(result['local_records_matching_kinetic_columns'], 1)
        self.assertEqual(result['not_currently_curated_ids'], ['1'])
        self.assertEqual(result['not_currently_displayed_ids'], [])
        self.public[0]['submitted_for_curation'] = False
        self.assertEqual(check_records(self.local, self.public)['not_currently_displayed_ids'], ['1'])

    def test_missing_id_cannot_count_as_a_match(self):
        self.public[0]['id'] = 2
        result = check_records(self.local, self.public)
        self.assertEqual(result['missing_ids'], ['1'])
        self.assertEqual(result['local_records_matching_kinetic_columns'], 0)

    def test_duplicate_public_id_stops_check(self):
        self.public.append(copy.deepcopy(self.public[0]))
        with self.assertRaisesRegex(ValueError, 'duplicate id'):
            check_records(self.local, self.public)

    def test_unrelated_incomplete_public_entry_does_not_block_local_check(self):
        self.public.append(dict(id=2, expressed=None, curated=False, submitted_for_curation=False))
        result = check_records(self.local, self.public)
        self.assertEqual(result['local_records_matching_all_columns'], 1)
        self.assertEqual(result['public_endpoint_records'], 2)

    def test_selected_records_must_belong_to_the_local_mutation(self):
        comparison = [dict(mutation_id='A10V', primary='True', highest_repeat='False',
                           all_repeated='True', source_record_ids='missing')]
        with self.assertRaisesRegex(ValueError, 'outside the local mutations'):
            check_records(self.local, self.public, comparison)

    def test_subset_join_uses_audit_mutation_case_without_changing_source_values(self):
        self.public[0]['resmut'] = 'v'
        self.local = [export_row(self.public[0])]
        comparison = [dict(mutation_id='A10V', primary='True', highest_repeat='False',
                           all_repeated='True', source_record_ids='1')]
        result = check_records(self.local, self.public, comparison)
        self.assertEqual(result['comparison_sets']['primary']['total_local_records'], 1)
        self.assertEqual(result['local_records_matching_all_columns'], 1)
        self.assertEqual(self.local[0]['Variant'], 'A10v')

    def test_changed_contributor_is_detected_without_exporting_the_name(self):
        self.public[0]['creator'] = 'invented_example_contributor'
        result = check_records(self.local, self.public)
        self.assertEqual(result['differences'][0]['fields'], ['Created by'])
        self.assertNotIn('invented_example_contributor', str(result))

    def test_formatter_keeps_yield_zero_but_matches_kinetic_zero_export_blank(self):
        self.public[0].update(yield_avg=0, KM_avg=0)
        row = export_row(self.public[0])
        self.assertEqual(row['Yield (mg/mL)'], '0')
        self.assertEqual(row['KM (mM)'], '')
        self.assertEqual(self.public[0]['KM_avg'], 0)

    def audit_row(self):
        return {"record_id": "1", "variant_norm": "A10V", "expressed_yes": "true",
                **{field: str(self.public[0][key]) if self.public[0][key] else ""
                   for field, key in AUDIT_KINETICS.items()}}

    def test_numerical_extract_check_finds_changed_values(self):
        local = [self.audit_row()]
        self.assertEqual(check_audit(local, self.public)['local_records_matching_kinetic_columns'], 1)
        self.public[0]['kcat_avg'] = 8
        result = check_audit(local, self.public)
        self.assertEqual(result['local_records_matching_kinetic_columns'], 0)
        self.assertEqual(result['differences'][0]['fields'], ['kcat_value'])

    def test_numerical_extract_checks_identity_and_expression_separately(self):
        local = [self.audit_row()]
        self.public[0].update(resmut='L', expressed=False, curated=False, submitted_for_curation=False)
        result = check_audit(local, self.public)
        self.assertEqual(result['local_records_matching_kinetic_columns'], 1)
        self.assertEqual(result['differences'][0]['fields'], ['variant_norm', 'expressed_yes'])
        self.assertEqual(result['not_currently_displayed_ids'], ['1'])

    def test_numerical_extract_missing_id_is_not_a_match(self):
        local = [self.audit_row()]
        self.public[0]['id'] = 2
        result = check_audit(local, self.public)
        self.assertEqual(result['missing_ids'], ['1'])
        self.assertEqual(result['local_records_matching_kinetic_columns'], 0)


if __name__ == '__main__':
    unittest.main()
