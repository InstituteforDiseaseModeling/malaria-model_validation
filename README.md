# malaria-model_validation

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- .[User Installation](#user-installation)
  - .[Developer Installation](#developer-installation)
  - .[Regular Installation](#regular-installation)
  - .[Pre-requisites](#pre-requisites)
- .[Run Simulations](#run-simulations)
  - .[Login to Comps](#login-to-comps)
  - .[Option 1_Run One Site with Python Scripts](#option-1_run-one-site-with-python-scripts)
  - .[Option 2_Run all Sites with Snakemake (Recommended)](option-2_run-all-sites-with-snakemake-recommended)


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

# Run Simulations

## Login to Comps
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
```bash
snakemake clean -j
snakemake -j
```