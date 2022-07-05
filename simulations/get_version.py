import argparse
import re
from importlib.metadata import version

import simulations.manifest as manifest

from idmtools.core.platform_factory import Platform
from idmtools.core import ItemType


def re_search(regex, line):
    match = re.search(regex, line)
    if match:
        return match.group(1)
    else:
        raise LookupError


def get_eradication_version(line):
    regex = 'EMOD Disease Transmission Kernel ' + r"(\d*\.*\d*\.*\d*\.*\d*)"
    return re_search(regex, line)


def get_eradication_branch(line):
    regex = 'from ' + r"(\w*\-\w*\(\d*\w*\d*\w*\d*\w*\))"
    return re_search(regex, line)


def get_eradication_info(exp_id=None, experiment=None):
    if not experiment and not exp_id:
        raise ValueError("Please provide an experiment object or an experiment id.")
    if not experiment:
        platform = Platform(manifest.platform_name)
        experiment = platform.get_item(item_id=exp_id, item_type=ItemType.EXPERIMENT)
    sim = experiment.simulations[0]

    # Get stdout.txt from Comps
    report = experiment.platform.get_files(sim, ['stdout.txt'])
    stdout = report['stdout.txt'].decode("utf-8")
    version = get_eradication_version(stdout)
    branch = get_eradication_branch(stdout)
    return version, branch


def get_emodpy_malaria_version():
    return version("emodpy_malaria")


def write_to_file(exp_id=None, experiment=None):
    if not experiment and not exp_id:
        raise ValueError("Please provide an experiment object or an experiment id.")
    if not experiment:
        era_version, era_branch = get_eradication_info(exp_id=exp_id)
    else:
        era_version, era_branch = get_eradication_info(experiment=experiment)
    emodpy_malaria_version = get_emodpy_malaria_version()
    version_file = manifest.version_file
    with open(version_file, 'w') as file:
        file.write(f"Eradication version: {era_version}\n")
        file.write(f"branch: {era_branch}\n")
        file.write(f"emodpy_malaria version: {emodpy_malaria_version}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='experiment id')
    parser.add_argument('--expid', '-e', type=str, help='experiment id', default='8a008695-bee5-ec11-a9f9-b88303911bc1')  # not sure if we want to make this required argument
    args = parser.parse_args()
    write_to_file(exp_id=args.expid)

