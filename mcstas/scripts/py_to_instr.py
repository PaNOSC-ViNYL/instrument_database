import os
MCSTAS_PATH = os.environ['MCSTAS']
import sys

instrument_file = os.path.abspath(sys.argv[1])
# parent directory of the facility:
facilitydir = os.path.abspath(os.path.dirname(instrument_file)+"/../../../")
sys.path.append(facilitydir)

(rest, filename) = os.path.split(instrument_file)
filename = filename[:-3] # remove .py
(rest, instrument_version) = os.path.split(rest)
(rest, instrument_name) = os.path.split(rest)


print("Reading instrument description from file: ", instrument_file)
print(" Instrument name: ", instrument_name, "\n",
      "Instrument version: ", instrument_version, "\n",
      "Filename: ", filename+".py")

from mcstasscript.interface import functions
my_configurator = functions.Configurator()
my_configurator.set_mcrun_path("/usr/bin/")
my_configurator.set_mcstas_path(MCSTAS_PATH)
print("McStas path: ",MCSTAS_PATH)

import importlib

my_instrument_module = importlib.import_module("instruments."+instrument_name+"."+instrument_version+"."+filename)
#print(vars())

#instrument = locals()[instrument_name]
instrument = my_instrument_module.def_instrument()
instrument.write_full_instrument()
