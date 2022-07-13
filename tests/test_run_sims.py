import unittest
from BaseTest import BaseTest
from simulations.run_sims import create_exp
import emod_malaria.bootstrap as dtk

import simulations.manifest as manifest
from simulations.get_eradication import get_eradication
import pathlib
import json


class RunSimsTest(BaseTest):
    @classmethod
    def setUpClass(cls) -> None:
        super(BaseTest, cls).setUpClass()
        era_folder = pathlib.Path(__file__).resolve().parent / 'stash'
        dtk.setup(era_folder)
        manifest.eradication_path = era_folder / 'Eradication'
        manifest.schema_file = era_folder / 'schema.json'
        get_eradication(use_local_eradication=False, my_manifest=manifest)
        cls.maxDiff = None

    def setUp(self) -> None:
        super(BaseTest, self).setUp()
        self.site = 'test_site'
        self.nSims = 2
        self.characteristic = False
        pass

    def test_submit_sim_config(self):
        exp = create_exp(characteristic=self.characteristic, nSims=self.nSims, site=self.site, my_manifest=manifest,
                         not_use_singularity=True)

        self.assertEqual(exp.simulation_count, self.nSims)

        for i, sim in enumerate(exp.simulations):
            # test config
            # config = sim.task.config
            # config.parameters.finalize()
            # with open("config.json", "w") as config_file:
            #     json.dump(config, config_file, indent=4, sort_keys=True)
            task = sim.task
            task.gather_transient_assets()
            config = task.transient_assets.assets[1]
            with open("inputs/new_config.json", "w") as config_file:
                # for debugging
                config_file.write(config.content)
            new_config_dict = json.loads(config.content)['parameters']

            with open("inputs/old_my_config.json", "r") as config_file:
                old_config_dict = json.load(config_file)['parameters']

            self.assertEqual(i, new_config_dict["Run_Number"])

            # this is not the final config file
            # todo: need to test the actual config file that we submitted to comps
            for key in ["Run_Number", "Campaign_Filename", "Custom_Reports_Filename",
                        "Enable_Interventions", "Enable_Termination_On_Zero_Total_Infectivity"]:
                new_config_dict.pop(key, None)
                old_config_dict.pop(key, None)
            self.assertDictEqual(new_config_dict, old_config_dict)
        pass

    def test_submit_sim_demog(self):
        exp = create_exp(characteristic=self.characteristic, nSims=self.nSims, site=self.site, my_manifest=manifest,
                         not_use_singularity=True)
        for i, sim in enumerate(exp.simulations):
            task = sim.task
            task.gather_transient_assets()
            demog = task.transient_assets.assets[0]
            print(demog.absolute_path)
            with open(demog.absolute_path, 'r') as demog_file:
                new_demog_dict = json.load(demog_file)

            with open("inputs/old_demog.json", "r") as demog_file:
                old_demog_dict = json.load(demog_file)
            self.assertDictEqual(new_demog_dict, old_demog_dict)
        pass

    def test_submit_sim_campaign(self):
        exp = create_exp(characteristic=self.characteristic, nSims=self.nSims, site=self.site, my_manifest=manifest,
                         not_use_singularity=True)
        for i, sim in enumerate(exp.simulations):
            task = sim.task
            task.gather_transient_assets()
            campaign = task.transient_assets.assets[2]
            with open("inputs/new_campaign.json", "w") as campaign_file:
                # for debugging
                campaign_file.write(campaign.content)
            new_campaign_dict = json.loads(campaign.content)

            with open("inputs/old_campaign.json", "r") as campaign_file:
                old_campaign_dict = json.load(campaign_file)
            self.assertDictEqual(new_campaign_dict, old_campaign_dict)
        pass

    def test_submit_sim_custom_report(self):
        exp = create_exp(characteristic=self.characteristic, nSims=self.nSims, site=self.site, my_manifest=manifest,
                         not_use_singularity=True)
        for i, sim in enumerate(exp.simulations):
            task = sim.task
            task.gather_transient_assets()
            custom_reports = task.transient_assets.assets[3]
            with open("inputs/new_custom_reports.json", "w") as custom_reports_file:
                # for debugging
                custom_reports_file.write(custom_reports.content)
            new_custom_reports_dict = json.loads(custom_reports.content)

            with open("inputs/old_custom_reports.json", "r") as custom_reports_file:
                old_custom_reports_dict = json.load(custom_reports_file)
            self.assertDictEqual(new_custom_reports_dict, old_custom_reports_dict)
        pass


if __name__ == '__main__':
    unittest.main()
