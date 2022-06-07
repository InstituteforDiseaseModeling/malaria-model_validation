from BaseTest import BaseTest
import unittest
from simulations.generate_site_rules import generate_rule, run
import os
import simulations.manifest as manifest


class GenerateRuleTest(BaseTest):
    def test_generate_rule(self):
        site = "test_site_name"
        n = 99
        script_name = "test.py"
        test_string = generate_rule(site=site, n=n, script_name=script_name)
        rules = [f'rule {site}_run_sim', f'rule {site}_analyzer', f'rule {site}_download']
        for rule in rules:
            self.assertIn(rule, test_string)

    def test_run(self):
        snakefile_bak = manifest.CURRENT_DIR / 'snakefile_bak'
        snakefile = manifest.CURRENT_DIR / 'snakefile_test'
        if os.path.isfile(snakefile):
            os.remove(snakefile)
        run(snakefile_bak=snakefile_bak, snakefile=snakefile)
        self.assertTrue(os.path.isfile(snakefile))


if __name__ == '__main__':
    unittest.main()
