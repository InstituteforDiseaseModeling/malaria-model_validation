#!/usr/bin/env python3

import pathlib # for a join

# idmtools ...
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment

# emodpy
from emodpy.emod_task import EMODTask
from emodpy_malaria.reporters.builtin import MalariaSummaryReport

from ewan_sim_match.helpers import *
import ewan_sim_match.params as params
import ewan_sim_match.manifest as manifest


def general_sim():
    """
    This function is designed to be a parameterized version of the sequence of things we do 
    every time we run an emod experiment. 
    """

    # Create a platform
    # Show how to dynamically set priority and node_group
    platform = Platform("SLURM")

    # create EMODTask 
    print("Creating EMODTask (from files)...")
    
    task = EMODTask.from_default2(
            config_path="my_config.json",
            eradication_path=manifest.eradication_path,
            ep4_custom_cb=None,
            campaign_builder=build_camp,
            schema_path=manifest.schema_file,
            param_custom_cb=set_param_fn,
            demog_builder=None,
        )

    # Create simulation sweep with builder
    builder = SimulationBuilder()

    # Add asset
    task.common_assets.add_asset("C:\\Users\\pselvaraj\\Github\\malaria-model_validation\\ewan_sim_match\\input_files\\demographics.json")

    # Sweep run number
    builder.add_sweep_definition(update_sim_random_seed, range(params.nSims))

    exp_name = params.exp_name

    # Add reporter
    reporter = MalariaSummaryReport()  # Create the reporter
    reporter.config(msr_config_builder, manifest)  # Config the reporter
    task.reporters.add_reporter(reporter)  # Add thre reporter

    # create experiment from builder
    print( f"Prompting for COMPS creds if necessary..." )
    experiment = Experiment.from_builder(builder, task, name=exp_name)

    # The last step is to call run() on the ExperimentManager to run the simulations.
    experiment.run(wait_until_done=True, platform=platform)

    # Check result
    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
        exit()

    print(f"Experiment {experiment.uid} succeeded.")

    # Save experiment id to file
    with open("COMPS_ID", "w") as fd:
        fd.write(experiment.uid.hex)
    print()
    print(experiment.uid.hex)


if __name__ == "__main__":
    # TBD: user should be allowed to specify (override default) erad_path and input_path from command line 
    # plan = EradicationBambooBuilds.MALARIA_LINUX
    # print("Retrieving Eradication and schema.json from Bamboo...")
    # get_model_files( plan, manifest )
    # print("...done.")

    general_sim()
