
# This is a user-modifiable Python file designed to be a set of simple input file and directory settings that you can choose and change.
from pathlib import Path
CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = CURRENT_DIR.parent

# The script is going to use this to store the downloaded schema file. Create 'download' directory or change to your preferred (existing) location.
schema_file= CURRENT_DIR / "download/schema.json"
# The script is going to use this to store the downloaded Eradication binary. Create 'download' directory or change to your preferred (existing) location.
#eradication_path="download/Eradication"
eradication_path=CURRENT_DIR / "download/Eradication"
# Create 'Assets' directory or change to a path you prefer. idmtools will upload files found here.
assets_input_dir= CURRENT_DIR / "Assets"
plugins_folder = CURRENT_DIR / "download/reporter_plugins"
analyzed_ouptut_path = PROJECT_DIR / 'EMOD_validation_recalibration/simulation_output'

# TODO: remove following lines
input_files_path = PROJECT_DIR / "simulation_inputs"
asset_path = input_files_path / "demographics_files/demographics_vital_1000.json"
# todo: add path for the csv file for site charatistic sweeps
simulation_coordinator_path = input_files_path / "simulation_coordinator.csv"

my_ep4_assets = None
requirements = PROJECT_DIR / "requirements.txt"

# Define Comps platform
platform_name = "Calculon"
priority = 'Normal'
node_group_private = 'idm_48cores'
node_group = 'idm_abcd'
