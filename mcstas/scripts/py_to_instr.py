import sys
import os
from instrumentdatabaseapi import instrumentdatabaseapi as API
from mcstasscript.interface import instr

institute = sys.argv[1]
instrument = sys.argv[2]
version = sys.argv[3]
software = "mcstas"
print(sys.argv)
if len(sys.argv) >= 5:
    flavour = sys.argv[4]
else:
    flavour = ""

sample_name = sys.argv[5]
sample_environment_name = sys.argv[6]

repo = API.Repository(local_repo=os.path.dirname(__file__) + "/../../")

myinstrument = repo.load(institute, instrument, version, software, flavour, dep=False)

myinstrument.set_sample_by_name(sample_name)
myinstrument.set_sample_environment_by_name(sample_environment_name)

calc_with_sample = myinstrument._calculator_with_sample

for calcname in myinstrument.calculators:
    calc = myinstrument.calculators[calcname]
    # calc.show_components()
    if calcname == calc_with_sample.name:
        calc.name = (
            calcname + "_" + flavour + "_" + sample_name + "_" + sample_environment_name
        )
    else:
        calc.name = calcname + "_" + flavour

    if isinstance(calc, instr.McStas_instr):
        calc.write_full_instrument()
