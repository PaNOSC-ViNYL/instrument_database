git clone --recurse-submodules https://github.com/PaNOSC-ViNYL/instrument_database.git
cd instrument_database
git checkout pixi
pixi install
pixi run pip install --no-deps --force-reinstall git+https://github.com/PaNOSC-ViNYL/McStasScript.git@custom_component_dir#egg=mcstasscript

pixi add ipykernel
pixi run python -m ipykernel install --user --name InstDB --display-name "Python (Instrument Database)"
