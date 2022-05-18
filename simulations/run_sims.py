#!/usr/bin/env python3
import argparse

# idmtools
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment

# emodpy
from emodpy.emod_task import EMODTask
from emodpy_malaria.reporters.builtin import add_report_intervention_pop_avg

from simulations.helpers import *
import simulations.params as params
from simulations import manifest as manifest


def general_sim(site=None, nSims=1, characteristic=False, priority=manifest.priority):
    """
    This function is designed to be a parameterized version of the sequence of things we do 
    every time we run an emod experiment. 
    """

    # Create a platform
    # Show how to dynamically set priority and node_group
    platform = Platform(manifest.platform_name, priority=priority, node_group=manifest.node_group)

    # create EMODTask 
    print("Creating EMODTask (from files)...")
    
    task = EMODTask.from_default2(
            config_path="my_config.json",
            eradication_path=str(manifest.eradication_path),
            ep4_custom_cb=None,
            campaign_builder=None,
            schema_path=str(manifest.schema_file),
            param_custom_cb=set_param_fn,
            demog_builder=None,
        )

    # add html intervention-visualizer asset to COMPS
    add_inter_visualizer = False
    if add_inter_visualizer:
        task.common_assets.add_asset(manifest.intervention_visualizer_path)
        add_report_intervention_pop_avg(task, manifest)

    # Create simulation sweep with builder
    builder = SimulationBuilder()


    exp_name = "validation_" + site

    # Sweep run number
    builder.add_sweep_definition(update_sim_random_seed, range(nSims))

    # Sweep sites and seeds - based on values in simulation_coordinator csv
    # builder.add_sweep_definition(set_simulation_scenario, [site])

    if characteristic:
        builder.add_sweep_definition(set_simulation_scenario_for_characteristic_site, [site])
    else:
        builder.add_sweep_definition(set_simulation_scenario_for_matched_site, [site])

    # create experiment from builder
    print( f"Prompting for COMPS creds if necessary..." )
    experiment = Experiment.from_builder(builder, task, name=exp_name)

    # The last step is to call run() on the ExperimentManager to run the simulations.
    experiment.run(wait_until_done=False, platform=platform)

    # Save experiment id to file
    with open(f"COMPS_ID/{site}_COMPS_ID_submit", "w") as fd:
        fd.write(experiment.uid.hex)
    print()
    print(experiment.uid.hex)
    return experiment.uid.hex


if __name__ == "__main__":
    # TBD: user should be allowed to specify (override default) erad_path and input_path from command line 
    # plan = EradicationBambooBuilds.MALARIA_LINUX
    # print("Retrieving Eradication and schema.json from Bamboo...")
    # get_model_files( plan, manifest )
    # print("...done.")
    
    parser = argparse.ArgumentParser(description='Process site name')
    parser.add_argument('--site', '-s', type=str, help='site name', default="dapelogo_2007")#params.sites[0]) # todo: not sure if we want to make this required argument
    parser.add_argument('--nSims', '-n', type=int, help='number of simulations', default=params.nSims)
    parser.add_argument('--characteristic', '-c', action='store_true', help='site-characteristic sweeps')
    parser.add_argument('--priority', '-p', type=str,
                        choices=['Lowest', 'BelowNormal', 'Normal', 'AboveNormal', 'Highest'],
                        help='Comps priority', default=manifest.priority)

    args = parser.parse_args()

    general_sim(site=args.site, nSims=args.nSims, characteristic=args.characteristic, priority=args.priority)
