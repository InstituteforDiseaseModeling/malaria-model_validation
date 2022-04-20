# malaria-model_validation

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- .[User Installation](#user-installation)
  - .[Developer Installation](#developer-installation)
  - .[Regular Installation](#regular-installation)
  - .[Pre-requisites](#pre-requisites)
- .[Run Simulations](#run-simulations)


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
```bash
cd simulations
python3 run_sims.py
```