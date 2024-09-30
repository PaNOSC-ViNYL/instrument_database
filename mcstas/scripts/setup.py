"""
This script sets up the configuration for mcstasscript for the current environment
"""

import os

if "MCSTAS" in os.environ:
    MCSTAS_PATH = os.environ["MCSTAS"]
    MCSTAS_BIN_PATH = os.environ["MCSTAS"]
elif "CONDA_PREFIX" in os.environ:
    MCSTAS_PATH = os.environ["CONDA_PREFIX"] + "/share/mcstas/resources"
    MCSTAS_BIN_PATH = os.environ["CONDA_PREFIX"]
else:
    raise RuntimeError(
        "CONDA not configured and MCSTAS environment variable not defined"
    )

from mcstasscript.interface import functions

my_configurator = functions.Configurator()
my_configurator.set_mcstas_path(MCSTAS_PATH)
my_configurator.set_mcrun_path(MCSTAS_BIN_PATH + "/bin/")
