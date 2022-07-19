
# This is a user-modifiable Python file designed to be a set of simple input file and directory settings that you can choose and change.
from pathlib import Path

# ========================================================
"""
This section contains the most important parameters that you may want to change before running the snakemake workflow.
use_local_eradication: set to 0 to run workflow using the Eradication and schema currently under download_folder.
                       set to 1 to get Eradicaiton and schema from emod_malaria and unzip them to download_folder.
singularity_id: the AC id for singularity image that the simulation will be run with.
                Set it to None if you don't want to run with any singularity image.
"""
use_local_eradication = 0
singularity_id = "8df53802-53f3-ec11-a9f9-b88303911bc1"
# ========================================================

CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = CURRENT_DIR.parent

# The script is going to use this to store the downloaded schema file. Create 'download' directory or change to your preferred (existing) location.
download_folder = CURRENT_DIR / "download"
schema_file = download_folder / "schema.json"
# The script is going to use this to store the downloaded Eradication binary. Create 'download' directory or change to your preferred (existing) location.
eradication_path = download_folder / "Eradication"
# Create 'Assets' directory or change to a path you prefer. idmtools will upload files found here.
assets_input_dir = CURRENT_DIR / "Assets"
plugins_folder = download_folder / "reporter_plugins"
analyzed_ouptut_path = PROJECT_DIR / "EMOD_validation_recalibration/simulation_output"
comps_id_folder = "COMPS_ID/"
suite_id_file = comps_id_folder + 'Suite'
eradication_found = comps_id_folder + 'eradication_found'
sif_id = comps_id_folder + 'sif.id'

# TODO: remove following lines
input_files_path = PROJECT_DIR / "simulation_inputs"
asset_path = input_files_path / "demographics_files/demographics_vital_1000.json"
simulation_coordinator_path = input_files_path / "simulation_coordinator.csv"
sweep_sim_coordinator_path = input_files_path / "sweep_sim_coordinator.csv"
intervention_visualizer_path = CURRENT_DIR / "download/index.html"

my_ep4_assets = None
requirements = PROJECT_DIR / "requirements.txt"

# Define Comps platform
platform_name = "Calculon"
priority = 'Normal'
node_group_private = 'idm_48cores'
node_group = 'idm_abcd'
