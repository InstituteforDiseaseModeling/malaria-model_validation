#!/usr/bin/env python3
import argparse
from idmtools.core.platform_factory import Platform
from idmtools.core import ItemType

import simulations.params as params
import simulations.manifest as manifest
from simulations.helpers import get_comps_id_filename


def check_experiment(site, platform=None):
    compid_file = get_comps_id_filename(site=site)
    with open(compid_file, 'r') as id_file:
        exp_id = id_file.readline()

    if not platform:
        platform = Platform(manifest.platform_name)
    experiment = platform.get_item(item_id=exp_id, item_type=ItemType.EXPERIMENT)

    # Wait for the experiment to be done in Comps.
    experiment.wait(wait_on_done_progress=True, refresh_interval=30)

    # Check result
    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
    else:
        print(f"Experiment {experiment.uid} succeeded.")

    return experiment.succeeded


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process site name')
    parser.add_argument('--site', '-s', type=str, help='site name', default=params.sites[0]) # not sure if we want to make this required argument
    args = parser.parse_args()
    check_experiment(args.site)

