"""
This script sets up the configuration for mcstasscript for the current environment
"""

import os

MCSTAS_PATH = os.environ["MCSTAS"]
from mcstasscript.interface import functions

my_configurator = functions.Configurator()
my_configurator.set_mcstas_path(MCSTAS_PATH)
my_configurator.set_mcrun_path(MCSTAS_PATH + "/bin/")
