from instrumentdatabaseapi import instrumentdatabaseapi as API

import sys
import os

# print(os.getenv("MCSTAS"))
repo = API.Repository(local_repo=".")
instrument_name = "ThALES"

myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", dep=False)
# myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", "merge", dep=False)


# import the units
import pint

ureg = pint.get_application_registry()

# setting the base directory for the simulation output
basedir = "/tmp/ThALES_scan/"
myinstrument.set_instrument_base_dir(basedir)

# generation energy (monochromator)
myinstrument.master["a2"] = myinstrument.energy_to_angle(4.98 * ureg.meV)
myinstrument.master["a4"] = 60 * ureg.degree
myinstrument.master["a6"] = myinstrument.master["a2"].pint_value

# myinstrument.calculators[myinstrument._calculator_name].settings(force_compile=False)
myinstrument.sim_neutrons(10000000)
myinstrument.set_instrument_base_dir(basedir + "generation/")
# myinstrument.run()


for rh in [2.0, 2.5, 3.00, 3.50, 4.00]:
    for rv in [1.0, 1.25, 1.5, 1.75, 2.0, 2.5]:
        bdir = basedir + "/rv_{0:0.2f}-rh_{1:0.2f}".format(rv, rh)
        if os.path.exists(bdir):
            continue
        os.mkdir(bdir)
        myinstrument.set_instrument_base_dir(bdir)
        mono = myinstrument.calculators[myinstrument._calculator_name].get_component(
            "Monochromator"
        )
        mono.RV = rv
        mono.RH = rh
        print("Simulating (RV,RH)= (", rv, rh, ")")
        myinstrument.run()


# python institutes/ILL/instruments/ThALES/HEAD/mcstas/other/mysim_rv.py
# gnuplot -c institutes/ILL/instruments/ThALES/HEAD/mcstas/other/mysim_rv.gpl
# zathura /tmp/ThALES_scan/optim_rv.pdf
