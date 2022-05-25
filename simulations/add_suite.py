from COMPS.Data import Experiment
from idmtools.core.platform_factory import Platform
from idmtools.entities import Suite
from simulations.load_inputs import load_sites
import simulations.manifest as manifest
from simulations.helpers import get_comps_id_filename


def add_suite(sites: list, suite_name: str = 'Malaria Model Validation Suite') -> object:
    """
    Add all experiments to the suite of the 1st experiment, if the 1st experiment doesn't belong to any suite, add them
    to a new suite.
    Args:
        sites (): list of site names
        suite_name (): str for suite name

    Returns:
        suite_id
    """
    platform = Platform(manifest.platform_name)
    first_exp_found = True
    for site in sites:
        exp_id_file = get_comps_id_filename(site, level=0)
        with open(exp_id_file, "r") as f:
            exp_id = f.readline().rstrip()
            exp = Experiment.get(exp_id)
            if first_exp_found:
                first_exp_found = False
                if exp.suite_id:
                    suite_id = exp.suite_id
                    print(f"Found existing suite {suite_id} for the 1st experiment, will save to this suite instead.")
                    continue
                else:
                    suite = Suite(name=suite_name)
                    suite.update_tags({'name': suite_name})
                    platform.create_items([suite])
                    suite_id = suite.id
                    print(f"Try to add all experiments to suite.id={suite_id}:")
            exp.suite_id = suite_id
            exp.save()
    with open(manifest.suite_id_file, "w") as file:
        file.write(str(suite_id))
    return suite_id


if __name__ == '__main__':
    all_sites, _, _ = load_sites()
    add_suite(all_sites)
