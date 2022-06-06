# malaria-model_validation

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [User Installation](#user-installation)
  - [Developer Installation](#developer-installation)
  - [Regular Installation](#regular-installation)
  - [Pre-requisites](#pre-requisites)
- [Run Simulations](#run-simulations)
  - [Login to Comps](#login-to-comps)
  - [Option 1_Run One Site with Python Scripts](#option-1_run-one-site-with-python-scripts)
  - [Option 2_Run all Sites with Snakemake (Recommended)](#option-2_run-all-sites-with-snakemake-recommended)


<!-- END doctoc generated TOC please keep comment here to allow auto update -->


# User Installation

## Developer Installation
If you want to install this project in a editable mode('develop' mode), run the following command: 
```bash
pip install -e . -r requirements.txt
```
or you can run 
```bash
pip install -e . -i https://packages.idmod.org/api/pypi/pypi-production/simple
```

## Regular Installation
```bash
pip install wheel
python3 setup.py bdist_wheel
pip install dist/{update-this-with-the-wheel-file-name}.whl --index-url=https://packages.idmod.org/api/pypi/pypi-production/simple
```

## Pre-requisites
- Python 3.9 x64

If you see the following error during installation:
![alt text](./datrie_error.png?raw=true)
You may need to install C++ Build Tools from this link: https://visualstudio.microsoft.com/visual-cpp-build-tools/

# Run Simulations

## Login to Comps
If you haven't login to Comps for a while, you can run the following script to login and cache your credential:
```bash
cd simulations
python3 create_auth_token_args.py -u {your-comps-username}
```
When you see "Password:" in the terminal, enter your comps password and hit enter key. 
![alt text](./comps_login.PNG?raw=true)


## Option 1_Run One Site with Python Scripts
```bash
cd simulations
python3 run_sims.py -s {site_name} -n {nSims}
python3 run_analyzers.py -s {site_name}
python3 download_wi.py -s {site_name}
```

Run Plotting and reportting scripts with site(s) that you ran:
```bash
Rscript create_plots\run_generate_validation_comparisons_site.R
python3 report\create_pdf_report_3.py
```

## Option 2_Run all Sites with Snakemake (Recommended)
Run the whole pipeline with all sites in simulation_coordinator.csv:
```bash
snakemake -j
```

- If you want to re-run the whole workflow, clean up your environment with "snakemake clean -j" and run "snakemake -j" again:
```bash
snakemake clean -j
snakemake -j
```

- If you make change locally in simulation_coordinator.csv, run the generate_site_rules.py to regenerate snakemake rules(run this script every time you update your simulation_coordinator.csv):
```bash
python3 generate_site_rules.py
snakemake -j
```

- If you want to re-run the analyzers steps with previous experiments you ran, you can delete the analyzer id files and run:
```bash
snakemake clean_ana clean_download -j
snakemake -j
```

- Simular to previous scenario, if you want to run only the download and plotting steps:
```bash
snakemake clean_download -j
snakemake -j
```

- If you want to re-run simulations for certain sites, delete COMPS ID files for those sites that you want to -re-run(/simulations/COMPS_ID/{site_name}_COMPS_ID_submit and _done files) and run "snakemake -j" again.

- If you want to re-run the analyzers and plotting steps for certain sites, delete the {site_name}_COMPS_ID_done files only and re-run "snakemake -j".

