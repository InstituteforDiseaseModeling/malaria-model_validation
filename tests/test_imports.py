import unittest
from BaseTest import BaseTest


class ImportTest(BaseTest):
    def setUp(self) -> None:
        super(BaseTest, self).setUp()
        self.expected_items = None
        self.found_items = None
        pass

    def verify_expected_items_present(self, namespace):
        self.found_items = dir(namespace)
        for item in self.expected_items:
            self.assertIn(
                item,
                self.found_items
            )

    def tearDown(self) -> None:
        pass

    def test_simulations_helpers_import(self):
        self.expected_items = [
            'update_sim_random_seed',
            'mAb_vs_EIR',
            'update_mab',
            'set_param_fn',
            'set_simulation_scenario',
            'build_standard_campaign_object',
            'build_camp',
            'add_hfca_hs',
            'add_hs_from_file',
            'add_nmf_hs',
            'add_nmf_hs_from_file',
            'ptr_config_builder',
            'add_broadcasting_survey',
            'build_demog',
            'get_comps_id_filename',
            'load_coordinator_df',
            'get_suite_id'
        ]
        import simulations.helpers as helpers
        self.verify_expected_items_present(namespace=helpers)
        pass

    def test_wait_for_experiment_import(self):
        self.expected_items = [
            'check_experiment'
        ]
        import simulations.wait_for_experiment as wait_for_experiment
        self.verify_expected_items_present(namespace=wait_for_experiment)
        pass

    def test_InsetChartAnalyzer_import(self):
        self.expected_items = [
            'InsetChartAnalyzer'
        ]
        import simulations.analyzers.InsetChartAnalyzer as InsetChartAnalyzer
        self.verify_expected_items_present(namespace=InsetChartAnalyzer)
        pass

    def test_ParDensAgeAnalyzer_import(self):
        self.expected_items = [
            'ParDensAgeAnalyzer'
        ]
        import simulations.analyzers.ParDensAgeAnalyzer as ParDensAgeAnalyzer
        self.verify_expected_items_present(namespace=ParDensAgeAnalyzer)
        pass

    def test_PatientReportAnalyzer_import(self):
        self.expected_items = [
            'PatientAnalyzer'
        ]
        import simulations.analyzers.PatientReportAnalyzer as PatientReportAnalyzer
        self.verify_expected_items_present(namespace=PatientReportAnalyzer)
        pass

    def test_PatientReportAnalyzer_laterDays_import(self):
        self.expected_items = [
            'PatientAnalyzer'
        ]
        import simulations.analyzers.PatientReportAnalyzer_laterDays as PatientReportAnalyzer_laterDays
        self.verify_expected_items_present(namespace=PatientReportAnalyzer_laterDays)
        pass

    def test_AnnualSummaryReportAnalyzer_import(self):
        self.expected_items = [
            'AnnualSummaryReportAnalyzer'
        ]
        import simulations.analyzers.AnnualSummaryReportAnalyzer as AnnualSummaryReportAnalyzer
        self.verify_expected_items_present(namespace=AnnualSummaryReportAnalyzer)
        pass

    def test_run_sims_import(self):
        self.expected_items = [
            'submit_sim'
        ]
        import simulations.run_sims as run_sims
        self.verify_expected_items_present(namespace=run_sims)
        pass

    def test_add_suite_import(self):
        self.expected_items = [
            'add_suite'
        ]
        import simulations.add_suite as add_suite
        self.verify_expected_items_present(namespace=add_suite)
        pass


if __name__ == '__main__':
    unittest.main()
