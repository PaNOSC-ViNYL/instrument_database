# #python institutes/ILL/instruments/D11/HEAD/mcstas/test.py | tail -12 >  institutes/ILL/instruments/D11/HEAD/mcstas/test/test.log; git diff institutes/ILL/instruments/D11/HEAD/mcstas/test/test.log

from instrumentdatabaseapi import instrumentdatabaseapi as API

import sys
import os

# import the units
import pint
import numpy as np
import h5py

# import matplotlib.pyplot as plt
from mcstasscript.interface.functions import load_metadata, load_monitor

repo = API.Repository(local_repo=".")
instrument_name = "D11"

repo.ls_flavours("ILL", instrument_name, "HEAD", "mcstas")
flavour = "full"
# flavour = "nosection"
myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", flavour, dep=False)


ureg = pint.get_application_registry()

# setting the base directory for the simulation output
basedir = "/tmp/" + instrument_name
myinstrument.set_instrument_base_dir(basedir)

simulation_dir = None
if len(sys.argv) == 2:
    simulation_dir = sys.argv[1]


myinstrument.sim_neutrons(80000000)
myinstrument.sim_neutrons(8000000)
myinstrument.sim_neutrons(800000)
myinstrument.set_seed(654321)


acquisition_time = 60  # seconds
detector_names = ["detector_central", "detector_left", "detector_right"]


def center_of_mass(data, nx_min, nx_max, ny_min, ny_max):
    cx = 0.0
    cy = 0.0
    nx = 0.0
    ny = 0.0
    for i in range(nx_min, nx_max + 1):
        cx += np.sum(data[i, :]) * i
        nx += np.sum(data[i, :])
    for i in range(ny_min, ny_max + 1):
        cy += np.sum(data[:, i]) * i
        ny += np.sum(data[:, i])
    return cx / nx + 1, cy / ny + 1


def set_tests(myinstrument, test_number):
    if test_number is not None:

        myinstrument.master["lambda"] = 6 * ureg.angstrom
        myinstrument.master["detpos"] = 2 * ureg.m
        myinstrument.master["attenuator_index"] = 0
        myinstrument.master["collimation"] = 8 * ureg.m
        myinstrument.master["bs_index"] = 0
        myinstrument.sample_holder(
            material="quartz", shape="box", w=0.02, h=0.03, d=0.0135, th=0.00125
        )
        myinstrument.sample_shape("holder")
        if test_number == 0:  # direct attenuated beam
            myinstrument.set_sample_by_name("None")
            myinstrument.sample_holder(None, None)
            myinstrument.master["attenuator_index"] = 6
        elif test_number == 1:  # direct beam with empty sample holder
            myinstrument.set_sample_by_name("None")
        elif test_number == 2:  # with sample
            myinstrument.set_sample_by_name("qSq")
            myinstrument.master[
                "qSq_file"
            ] = '"./institutes/ILL/instruments/D11/HEAD/mcstas/data/simul_5711.sq"'
            #                '"simul_5711.sq"'
        elif test_number == -1:  # direct beam no beamstop
            myinstrument.set_sample_by_name("None")
            myinstrument.sample_holder(None, None)
            myinstrument.master["attenuator_index"] = 6
            myinstrument.master["bs_index"] = -1
        else:
            raise RuntimeError(f"Test number {test_number} out of range")


def read_test(myinstrument, test_number, acquisition_time):
    calcname = "OriginCalc"
    calcname_data = calcname + "_data"

    set_tests(myinstrument, test_number)
    metadata_list = load_metadata(simulation_dir)
    # print(left_detector_data.shape, central_detector_data.shape, right_detector_data.shape)
    # print(metadata_list)

    detectors_simulation = {}
    detectors_trueMC = {}
    for detector in metadata_list:
        #        if detector.name in detector_names or detector.name in ["PSD_sample"]:
        #            detectors_data[detector.name] = detector.Intensity * acquisition_time
        #            detectors_trueMC[detector.name] = detector.Ncount
        #    return detectors_data, detectors_trueMC
        if detector.component_name in detector_names or detector.component_name in [
            "PSD_sample"
        ]:
            monitor = load_monitor(detector, simulation_dir)
            detectors_simulation[monitor.name] = monitor.Intensity * acquisition_time
            detectors_trueMC[monitor.name] = monitor.Ncount
        # print("### MONITOR:")
        # help(monitor)

    return detectors_simulation, detectors_trueMC


def run_test(myinstrument, test_number, acquisition_time):
    calcname = "OriginCalc"
    calcname_data = calcname + "_data"

    set_tests(myinstrument, test_number)

    myinstrument.run()
    data = myinstrument.output
    # detectors = data[calcname_data].get_data()["data"]
    detectors = data.get_data()["data"]
    print(detectors)
    detectors_data = {}
    detectors_trueMC = {}
    for detector in detectors:
        if detector.name in detector_names or detector.name in ["PSD_sample"]:
            detectors_data[detector.name] = detector.Intensity * acquisition_time
            detectors_trueMC[detector.name] = detector.Ncount
    return detectors_data, detectors_trueMC


def data_test(test_number):
    file = ""
    if test_number == 0 or test_number == -1:  # direct attenuated beam
        file = "institutes/ILL/instruments/D11/HEAD/mcstas/data/005708.nxs"
    elif test_number == 1:  # direct beam with empty sample holder
        file = "institutes/ILL/instruments/D11/HEAD/mcstas/data/005721.nxs"
    elif test_number == 2:
        file = "institutes/ILL/instruments/D11/HEAD/mcstas/data/005711.nxs"
    else:
        raise RuntimeError(f"Test number {test_number} out of range")

    # Read NeXus
    f = h5py.File(file, "r")

    detectors_data = {
        "detector_left": f["entry0"]["D11"]["Detector 2"]["data"][:, :, 0],
        "detector_central": f["entry0"]["D11"]["Detector 1"]["data"][:, :, 0],
        "detector_right": f["entry0"]["D11"]["Detector 3"]["data"][:, :, 0],
    }
    return detectors_data


def compare(d, s, mc):
    print(
        "{:18}: {:>10} {:>10} | {:>10} | {:>10}".format(
            "", "data", "sim", "data/sim", "MC count"
        )
    )
    for key in detector_names:
        dsum = np.sum(d[key])
        ssum = np.sum(s[key])
        mcsum = np.sum(mc[key])
        ratio = dsum / ssum
        print(
            "{:18}: {:10.0f} {:10.0f} | {:10.2f} | {:10.0f}".format(
                key, dsum, ssum, ratio, mcsum
            )
        )


if simulation_dir is not None:
    intensity, count = read_test(myinstrument, 0, acquisition_time)
    data = data_test(0)
    compare(data, intensity, count)
    sys.exit(0)
darr = []
sarr = []
mcarr = []
ntests = 0
for itest in range(-1, ntests):
    darr.append(data_test(itest))
    intensity, count = run_test(myinstrument, itest, acquisition_time)
    sarr.append(intensity)
    mcarr.append(count)

    # print(sarr[0]["PSD_sample"])
    # sys.exit(0)

for itest in range(-1, ntests):
    d = darr[itest]
    s = sarr[itest]
    mc = mcarr[itest]
    compare(d, s, mc)

# print("PSD_sample: ", center_of_mass(sarr[0]["PSD_sample"], 0, 100, 0, 100))
sys.exit(0)
print(
    detectors_simulation["left"].Intensity.transpose().shape,
    detectors_simulation["central"].Intensity.transpose().shape,
    detectors_simulation["right"].Intensity.transpose().shape,
)
ratio = {}
ratio["left"] = (
    detectors_simulation["left"].Intensity.transpose() + detectors_data["left"]
)


Intensities = []
import numpy as np

a = run_test(myinstrument, 0)
print(np.sum(a))
sys.exit(0)
for itest in range(0, 3):
    Intensities.append(run_test(myinstrument, itest))

for itest in range(0, 3):
    print(itest, " : ", Intensities[itest])

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
