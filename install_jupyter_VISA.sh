#!/bin/bash
set -e
set -o pipefail
echo "Installing the necessary packages and softwares as subdirectories of the current directory"

# first clone the instrument database
git clone  --recurse-submodules https://github.com/PaNOSC-ViNYL/instrument_database.git
cd instrument_database/

# create a python environment
python3 -m venv --system-site-packages --symlinks python_packages

# activate the environment and upgrade it
source python_packages/bin/activate
pip install --upgrade pip

# install the dependencies
pip install -r instrumentDataBaseAPI/requirements.txt
pip install -e instrumentDataBaseAPI/

# clone some python packages from GIT since they are in development mode
git clone --depth 1 -b custom_component_dir https://github.com/PaNOSC-ViNYL/McStasScript.git
pip install -e McStasScript

git clone --depth 1 https://github.com/PaNOSC-ViNYL/libpyvinyl.git
pip install -e libpyvinyl/

# creating an environment for jupyter lab
python -m ipykernel install --user --name=simulation

sed -i '2 i "env":{"MCSTAS":"/usr/share/mcstas/3.1"},' ~/.local/share/jupyter/kernels/simulation/kernel.json 
python -c "
MCSTAS_PATH = '/usr/share/mcstas/3.1'
from mcstasscript.interface import functions
my_configurator = functions.Configurator()
my_configurator.set_mcstas_path(MCSTAS_PATH)
my_configurator.set_mcrun_path(MCSTAS_PATH + '/bin/')
"

