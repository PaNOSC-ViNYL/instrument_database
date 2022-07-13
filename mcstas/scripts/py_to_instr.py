import sys
import os
from instrumentdatabaseapi import instrumentdatabaseapi as API
from mcstasscript.interface import instr

institute =sys.argv[1]
instrument=sys.argv[2]
version=sys.argv[3]
software="mcstas"
print(sys.argv)
if len(sys.argv)==5:
    flavour=sys.argv[4]
else:
    flavour=""

repo = API.Repository(local_repo=os.path.dirname(__file__)+"/../../")

myinstrument = repo.load(institute, instrument, version, software, flavour, dep=False)
for calcname in myinstrument.calculators:
    calc=myinstrument.calculators[calcname]
    #calc.show_components()
    if isinstance(calc, instr.McStas_instr):
        calc.write_full_instrument()
