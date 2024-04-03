from instrumentdatabaseapi import instrumentdatabaseapi as API

import sys
import os

# print(os.getenv("MCSTAS"))
repo = API.Repository(local_repo=".")
instrument_name = "ThALES"
instrument_name = "Panther"
instrument_name = "D11"

repo.ls_flavours("ILL", instrument_name, "HEAD", "mcstas")
flavour = "full"
flavour = "nosection"
myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", flavour, dep=False)
# myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", "merge", dep=False)


# import the units
import pint

ureg = pint.get_application_registry()

# setting the base directory for the simulation output
basedir = "/tmp/" + instrument_name
myinstrument.set_instrument_base_dir(basedir)

# generation energy (monochromator)
# myinstrument.master["a2"] = myinstrument.energy_to_angle(4.98 * ureg.meV)
# myinstrument.master["a4"] = 60 * ureg.degree
# myinstrument.master["a6"] = myinstrument.master["a2"].pint_value
# print(myinstrument.get_total_SPLIT())
# myinstrument.set_sample_by_name("vanadium")
# myinstrument.set_sample_by_name("H2O")
# myinstrument.sample_cylinder_shape(0.005, 0.01)
# print(myinstrument)
myinstrument.sim_neutrons(500000)
myinstrument.set_seed(654321)


test_number = 0  # None 1 or 2


def set_tests(myinstrument, test_number):
    instrument_name = myinstrument.name

    if instrument_name == "ThALES":
        myinstrument.master["a2"] = myinstrument.energy_to_angle(4.98 * ureg.meV)
        myinstrument.master["a4"] = 60 * ureg.degree
        myinstrument.master["a6"] = myinstrument.master["a2"].pint_value
    if instrument_name == "Panther":
        myinstrument.master["energy"] = 19 * ureg.meV
    if instrument_name == "D11":
        if test_number is not None:

            myinstrument.master["detpos"] = 2 * ureg.m
            myinstrument.master["attenuator_index"] = 0
            myinstrument.master["collimation"] = 8 * ureg.m
            myinstrument.sample_holder(
                material="quartz", shape="box", w=0.02, h=0.03, d=0.0135, th=0.00125
            )
            myinstrument.sample_shape("holder")
            if test_number == 0:
                myinstrument.set_sample_by_name("None")
                myinstrument.sample_holder(None, None)
                myinstrument.master["attenuator_index"] = 6
            elif test_number == 1:
                myinstrument.set_sample_by_name("None")
            elif test_number == 2:
                myinstrument.set_sample_by_name("qSq")
                myinstrument.master[
                    "qSq_file"
                ] = '"/users/nourbakhsh/digitaltwin/instrument_database/institutes/ILL/instruments/D11/HEAD/mcstas/data/simul_5711.sq"'


def run_test(myinstrument, test_number):
    calcname = "OriginCalc"
    calcname_data = calcname + "_data"

    set_tests(myinstrument, test_number)

    myinstrument.run()
    data = myinstrument.output
    detectors = data[calcname_data].get_data()["data"]
    for detector in detectors:
        if detector.name == "detector_central":
            return detector.Intensity


Intensities = []


for itest in range(2, 3):
    Intensities.append(run_test(myinstrument, itest))

for itest in range(2, 3):
    print(itest, " : ", Intensities[itest][1])

sys.exit(0)

calcname = "OriginCalc"
attenuation_values = [
    8.325,  # attenuator 1
    26.21,  # attenuator 2
    72.23,  # attenuator 3
    216.5,  # attenuator 1+2
    594.6,  # attenuator 1+3
    1702,  # attenuator 2+3
    13480,  # attenuator 1+2+3
]
for att in range(0, len(attenuation_values)):
    myinstrument.calculators["OriginCalc"].parameters["attenuator_index"] = att
    #    myinstrument.run()
    data = myinstrument.output
    calcname_data = calcname + "_data"
    detectors = data[calcname_data].get_data()["data"]
    for detector in detectors:
        print(detector.name)
        if detector.name == "PSD_attenuator":
            Intensities.append(detector.Intensity)
        # print(detector, "\n\t I = ",detector.Intensity, "\n\t E =", detector.Error, "\n\t N =", detector.Ncount)

for att in range(0, len(attenuation_values)):
    print(
        att, " : ", Intensities[att][1], Intensities[att][1] * attenuation_values[att]
    )
# print("Ltof: "+str(myinstrument.calcLtof(myinstrument.calculators["OriginCalc"], "Chopper0", "Chopper1", True)))
# sys.exit(0)
# diagnostics
# for calc in myinstrument.calculators.values():
#    calc.show_diagram(analysis=True)
# mycalc = myinstrument.calculators["OriginCalc"]
# mycalc.show_diagram(analysis=True)
import mcstasscript as ms

for att in range(0, len(attenuation_values)):
    myinstrument.calculators["OriginCalc"].parameters["attenuator_index"] = att
    print(
        att,
        attenuation_values[att],
        myinstrument.calculators["OriginCalc"].parameters["attenuator_index"],
    )
# myinstrument.run()
sys.exit(0)
np = 21
np = (np - 1) / 2
dEI = 0.05
import numpy

for energy in numpy.arange(4.98 - np * dEI, 4.98 + np * dEI, dEI):
    angle = myinstrument.energy_to_angle(energy * ureg.meV)
    print(f"{energy}\t{angle}")

myinstrument.run()
sys.exit(0)

# myinstrument.calculators[myinstrument._calculator_name].settings(force_compile=False)
myinstrument.sim_neutrons(10000000)
myinstrument.set_instrument_base_dir(basedir + "generation/")

# the calculator has the same name as the instrument just for ThALES because it has been created like this. It is not true in general
# myThALES = myinstrument.calculators[instrument_name]
# myThALES.show_diagram()
# myThALES.show_components()

# myinstrument.set_sample_by_name("sqw")
print(myinstrument)
print("------------------------------")
print("Implemented samples: ", myinstrument.samples)
print("Current sample name: ", myinstrument.sample_name)
print("Current sample object: \n", myinstrument.sample)

print("------------------------------")
print("Implemented sample environments: ", myinstrument.sample_environments)
print("Current sample environment name: ", myinstrument.sample_environment_name)
print("Current sample environment object: \n", myinstrument.sample_environment)

# sys.exit(0)
s = myinstrument.sample
# myinstrument.run()

# myinstrument.sim_neutrons(5000)
# myinstrument.run()
# print((vars(s).values()))
# print("hash: ", (hash(vars(s).values())))
# myinstrument.sample.thickness = 0.002
# print("hash after:", (frozenset(vars(s))))
# myinstrument.run()

for rv in [1.0, 1.5, 2.0, 2.5, 3.0]:
    bdir = basedir + "/rv_{0:.2}".format(rv)
    os.mkdir(bdir)
    myinstrument.set_instrument_base_dir(bdir)
    mono = myinstrument.calculators[myinstrument._calculator_name].get_component(
        "Monochromator"
    )
    mono.RV = rv
    myinstrument.run()


sys.exit(0)

myinstrument = repo.load(
    "ILL", instrument_name, "HEAD", "mcstas", "from_sample", dep=False
)
myinstrument.set_instrument_base_dir(basedir)
# myThALES = myinstrument.calculators[instrument_name]
# myThALES.parameters["filelist"] = '"/tmp/ThALES_scan/generation/sSAMPLE.mcpl.gz"'
# print(myinstrument)


sys.exit(0)
import os

for r in [0.002, 0.005, 0.01, 0.02]:
    for h in [0.005, 0.01, 0.02]:
        os.mkdir(basedir + "r{0:.3f}_h{1:.3f}/".format(r, h))
        myinstrument.master["sample_size_r"] = r
        myinstrument.master["sample_size_y"] = h
        energy = 4.48  # starting energy of the scan
        dE = 0.05  # scan steps
        myThALES.show_components()
        # sys.exit(0)

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
            myThALES.calculator_base_dir = "r{0:.3f}_h{1:.3f}/{2:.2f}/".format(
                r, h, energy
            )
            myinstrument.run()
            energy += dE

# for r in /tmp/ThALES_scan/r0.0*/;  for dir in $r/{4,5}*; set -l d (basename $dir); echo -en "$d\t";tail -1 $dir/detector_all.dat; end | sort > $r/scan.dat;end

#  for file in /tmp/ThALES_scan/r*/scan.dat; gnuplot -persistent -c scan.gpl  $file;end
# gnuplot  -c mysim.gpl /tmp/ThALES_scan/4.48_5/
# for dir in /tmp/ThALES_scan/*; set -l d (basename $dir); echo -en "$d\t";tail -1 $dir/detector_all.dat; end | sort > /tmp/scan.dat'
# gnuplot p '/tmp/scan.dat' u 1:2:3 w yerr
# A=0.1
# sigma=0.05
# mu=5
# gauss(x) = A*exp(-(x-mu)*(x-mu)/(sigma*sigma))
# fit gauss(x) '/tmp/scan.dat' u 1:2 via A,sigma,mu

# set fit errorvariables
# set samples 500

# cd "/tmp/ThALES_scan/4.48"
# set xrange [4:5]
# p  'H5.E',    'H53_D.E', 'H53_7.E','slit_mono.E',  'ThALES_mono_in.E',  'ThALES_mono_out.E'
