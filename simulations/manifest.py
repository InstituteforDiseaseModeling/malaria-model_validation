
# This is a user-modifiable Python file designed to be a set of simple input file and directory settings that you can choose and change.
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = CURRENT_DIR.parent
DOWNLOAD_DIR = CURRENT_DIR / "download"

# The script is going to use this to store the downloaded schema file. Create 'download' directory or change to your preferred (existing) location.
schema_file = DOWNLOAD_DIR / "schema.json"
# The script is going to use this to store the downloaded Eradication binary. Create 'download' directory or change to your preferred (existing) location.
eradication_path = DOWNLOAD_DIR / "Eradication"
plugins_folder = DOWNLOAD_DIR / "reporter_plugins"

# Create 'Assets' directory or change to a path you prefer. idmtools will upload files found here.
assets_input_dir = CURRENT_DIR / "Assets"
# analyzed_ouptut_path = PROJECT_DIR / "EMOD_validation_recalibration" / "simulation_output"
comps_id_folder = "COMPS_ID/"
suite_id_file = comps_id_folder + 'Suite'
version_file = comps_id_folder + "version.txt"

simulation_output_filepath = CURRENT_DIR / "output"
benchmark_simulation_filepath = CURRENT_DIR / "output"
input_files_path = PROJECT_DIR / "simulation_inputs"
base_script_plot_filepath = PROJECT_DIR / "create_plots"
base_reference_filepath = PROJECT_DIR / "reference_datasets"
plot_output_filepath = PROJECT_DIR / "report" / "_plots"

# TODO: remove following lines
asset_path = input_files_path / "demographics_files" / "demographics_vital_1000.json"
simulation_coordinator_path = input_files_path / "simulation_coordinator.csv"
sweep_sim_coordinator_path = input_files_path / "sweep_sim_coordinator.csv"
intervention_visualizer_path = DOWNLOAD_DIR / "index.html"

my_ep4_assets = None
requirements = PROJECT_DIR / "requirements.txt"

# Define Comps platform
platform_name = "Calculon"
priority = 'Normal'
node_group_private = 'idm_48cores'
node_group = 'idm_abcd'
