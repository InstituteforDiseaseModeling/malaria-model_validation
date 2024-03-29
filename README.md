# malaria-model_validation

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [User Installation](#user-installation)
  - [Installation](#installation)
  - [Pre-requisites](#pre-requisites)
- [Run workflow](#run-workflow)
  - [Login to Comps](#login-to-comps)
  - [Option 1: Run all Sites with Snakemake](#option-1-run-all-sites-with-snakemake)
  - [Option 2: Run Sites in Certain Subset(s) with Snakemake](#option-2-run-sites-in-certain-subsets-with-snakemake)
    - [How to Change Default Setting](#how-to-change-default-setting)
    - [Snakemake Tips](#snakemake-tips)
- [Check Plots and Final Report](#check-plots-and-final-report)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


# User Installation
Note: we recemend to upgrade pip to the latest version before installation:
```bash
pip install --upgrade pip
```

## Installation
Run the following command: ( it will install this project in editable/'develop' mode): 
```bash
pip install -e . -r requirements.txt
```
or you can run 
```bash
pip install -e . -i https://packages.idmod.org/api/pypi/pypi-production/simple
```

## Pre-requisites
- Python 3.9 x64

If you see the following error during installation:
![alt text](./datrie_error.png?raw=true)
You may need to install C++ Build Tools from this link: https://visualstudio.microsoft.com/visual-cpp-build-tools/

If you see error about conflicting dependencies on Pandas:
```bash
ERROR: Cannot install idmtools-cli and malaria-model-validation because these package versions have conflicting dependencies.
The conflict is caused by:
    datar 0.0.0 depends on pandas<2.0 and >=1.2
    plotnine 0.1.0 depends on pandas>=0.19.0
    idmtools 1.6.6 depends on pandas<1.2 and >=1.1.4
```
Idmtools has dependencies with an old pandas version, we have opened a ticket to update this dependency on idmtools and are waiting for a patch release: https://github.com/InstituteforDiseaseModeling/idmtools/issues/1774.
For now, you can workarround this error by removing 'datar' from our requirements.txt and install it after installating our package with "pip install datar".

# Run Workflow

## Login to Comps
If you haven't login to Comps for a while, you can run the following script to login and cache your credential:
```bash
cd simulations
python3 create_auth_token_args.py -u {your-comps-username}
```
When you see "Password:" in the terminal, enter your comps password and hit enter key. 
![alt text](./comps_login.PNG?raw=true)


## Option 1: Run all Sites with Snakemake
Run the snakemake pipeline with all sites in simulation_coordinator.csv:
```bash
snakemake -j
```

## Option 2: Run Sites in Certain Subset(s) with Snakemake
Run the snakemake pipeline with sites in one or multiple subsets in simulation_coordinator.csv:
```bash
snakemake --config -s="core_relationship" -j
```
or 
```bash
snakemake --config -s="core_relationship, infection_duration" -j
```


### How to Change Default Setting
- Some details about our default setting:
  - By default, the workflow will copy the Eradication and schema files from corresponding emod_malaria package to 
    "simulations\download" folder and run experiments with them. If you would like to run against a certain copy of 
    Eradication, you can set "use_local_eradication = 1" in "\simulations\manifest.py" and copy the Eradication and 
    schema files to "simulations\download" folder by yourself.
  - By default, the workflow will run simulations in singularity container in Comps. The singularity image used in the 
    workflow is defined as "singularity_id" in "\simulations\manifest.py". Set this id to None will disable running with
    singularity container and run directly in Comps environment(CentOS + Python 3.6 in Calculon).
 
### Snakemake Tips
- Some snakemake tips about running the workflow:
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

  - If you want to re-run the analyzers steps with previous experiments you ran, you can delete the analyzer id files and all following the analyzers steps and run:
  ```bash
  snakemake clean_ana -j
  snakemake -j
  ```

  - Simular to previous scenario, if you want to only re-run starting from the download steps:
  ```bash
  snakemake clean_download -j
  snakemake -j
  ```
  
  - Simular to previous scenario, if you want to only re-run starting from the plotting steps:
  ```bash
  snakemake clean_plot -j
  snakemake -j
  ```
  - Simular to previous scenario, if you want to only re-run starting from the reporting steps:
  ```bash
  snakemake clean_report -j
  snakemake -j
  ```

    - If you want to re-run simulations for certain sites, delete COMPS ID files for those sites that you want to -re-run 
      (/simulations/COMPS_ID/{site_name}_COMPS_ID_exp_submit, _analyzers and _download files) and run
      "snakemake clean_plot -j" and then "snakemake -j" again.

    - If you want to re-run the analyzers and plotting steps for certain sites, delete the {site_name}_COMPS_ID_analyzers 
      and _download files only and run "snakemake clean_plot -j" and then "snakemake -j".

- Other snakemake tips:
  - Snakemake log can be found in '\malaria-model_validation\simulations\.snakemake\log'.
  - If your working directory is locked by snakemake by accident (for example, your machine is powered off with a 
    running Snakemake instance will cause a stale lock), you can run "snakemake --unlock -j" to remove the stale lock.
  - If you have incompleted file from previous snakemake instance, you can run "snakemake --cleanup-metadata 
    <filenames>" to clean up the metadata. Or you can delete the '.snakemake' folder to clean up the metadata manually.

# Check Plots and Final Report
You can check the plots in \malaria-model_validation\report\_plots folder and final report:
\malaria-model_validation\report\Malaria_model_validation_output_{date}}({time}).pdf.