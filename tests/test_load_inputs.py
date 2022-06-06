from BaseTest import BaseTest
import unittest
from simulations.load_inputs import load_sites


class LoadInputsTest(BaseTest):
    def test_load_sites(self):
        sites, nSims, script_names = load_sites()
        self.assertEqual(len(sites), len(nSims))
        self.assertEqual(len(sites), len(script_names))


if __name__ == '__main__':
    unittest.main()
