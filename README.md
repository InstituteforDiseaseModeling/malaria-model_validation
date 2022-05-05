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
python3 run_sims.py -s {site_name} -n {nSims}
python3 wait_for_experiment.py -s {site_name}
```

## Option 2_Run all Sites with Snakemake (Recommended)
Run the generate_site_rules.py to generate snakemake rules based on information in simulation_coordinator.csv(run this script every time you update your simulation_coordinator.csv):
```bash
python3 generate_site_rules.py
```

Run "snakemake clean" to clean up your environment:
```bash
snakemake clean -j
```

Run the whole pipeline with all sites in simulation_coordinator.csv:
```bash
snakemake -j
```