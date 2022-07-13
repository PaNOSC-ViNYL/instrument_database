from instrumentdatabaseapi import instrumentdatabaseapi as API

# import os
# print(os.getenv("MCSTAS"))
repo = API.Repository(local_repo=".")
instrument_name = "ThALES"
myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", dep=False)
myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", "merge", dep=False)
import sys

sys.exit(0)
# print(myinstrument)
myThALES = myinstrument.calculators[instrument_name]
# myThALES.show_components()

import pint

ureg = pint.get_application_registry()

a2 = myinstrument.parameters["ThALES"]["a2"]
a2.energy = 4.98 * ureg.meV

basedir = "/tmp/ThALES_scan/"
myinstrument.set_instrument_base_dir(basedir)

a4 = myinstrument.parameters["ThALES"]["a4"]
a3 = myinstrument.parameters["ThALES"]["a3"]
a3.value = 0 * ureg.degree
a4.value = 50 * ureg.degree

a6 = myinstrument.parameters["ThALES"]["a6"]

dE = 0.05

energy = 4.98
# myThALES.settings(force_compile=False)

myThALES.settings(ncount=10000000)
a6.energy = a2.energy
a4.value = 60 * ureg.degree
myThALES.calculator_base_dir = "/tmp/test/"
# myThALES.show_diagram()
myThALES.show_components()
# import sys

# sys.exit(0)

# myThALES.show_instrument()
print(myinstrument)
myinstrument.run()


sys.exit(0)
energy = 4.48
for i in range(0, 21):
    if i == 0:
        myThALES.settings(force_compile=True)
    else:
        myThALES.settings(force_compile=False)

    a6.energy = energy * ureg.meV
    print("Scan Energy: %.2f" % energy)
    print(f"Scan Energy: {a6.energy}")
    print(f"Scan angle: {a6.value}")
    print(f"Generate energy: {a2.energy}")
    print(f"Generate angle: {a2.value}")
    myThALES.calculator_base_dir = "{0:.2f}/".format(energy)
    myinstrument.run()
    energy += dE

import os

# gnuplot  -c mysim.gpl /tmp/ThALES_scan/4.48_5/
# for dir in /tmp/ThALES_scan/*; set -l d (basename $dir); echo -en "$d\t";tail -1 $dir/detector_all.dat; end | sort > /tmp/scan.dat'
# gnuplot p '/tmp/scan.dat' u 1:2:3 w yerr


# cd "/tmp/ThALES_scan/4.48"
# set xrange [4:5]
# p  'H5.E',    'H53_D.E', 'H53_7.E','slit_mono.E',  'ThALES_mono_in.E',  'ThALES_mono_out.E'
