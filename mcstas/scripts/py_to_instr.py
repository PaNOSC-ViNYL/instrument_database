import sys
instrument_file = sys.argv[1]
instrument_name = sys.argv[2]

print("Reading instrument: ", instrument_name, " from file: ", instrument_file)

from mcstasscript.interface import functions
my_configurator = functions.Configurator()
my_configurator.set_mcrun_path("/usr/bin/")
my_configurator.set_mcstas_path("/usr/share/mcstas/2.7/")

exec(open(instrument_file).read())
#D22_quick.write_full_instrument()
#print(vars())

instrument = locals()[instrument_name]
instrument.write_full_instrument()