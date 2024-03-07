# #python institutes/ILL/instruments/D11/HEAD/mcstas/test.py | tail -12 >  institutes/ILL/instruments/D11/HEAD/mcstas/test/test.log; git diff institutes/ILL/instruments/D11/HEAD/mcstas/test/test.log

from instrumentdatabaseapi import instrumentdatabaseapi as API
from D11 import def_instrument
import sys, os


import pytest
import pytest_functions as pyf


# import the units
import pint
import numpy as np
import h5py

import matplotlib.pyplot as plt
from mcstasscript.interface.functions import load_metadata, load_monitor


repo = API.Repository(local_repo=".")
instrument_name = "D11"

repo.ls_flavours("ILL", instrument_name, "HEAD", "mcstas")
flavour = "full"
# flavour = "nosection"
flavour = "simplefull"
# myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", flavour, dep=False)
myinstrument = def_instrument(flavour)

ureg = pint.get_application_registry()

# setting the base directory for the simulation output
basedir = "/tmp/" + instrument_name
myinstrument.set_instrument_base_dir(basedir)

myinstrument.set_seed(654321)

# for calc in myinstrument.calculators:
#    myinstrument.calculators[calc].settings(mpi=8)

detector_names = ["detector_central", "detector_left", "detector_right"]


def set_test_dir(instrument_name, itest, base="/tmp"):
    return base + "/{}/test_{:d}".format(instrument_name, itest)


def get_detector_data(file):
    # file = myinstrument.test_datafile(test_number)
    # Read NeXus
    f = h5py.File(file, "r")

    detectors_data = {
        "detector_left": f["entry0"]["D11"]["Detector 2"]["data"][:, :, 0],
        "detector_central": f["entry0"]["D11"]["Detector 1"]["data"][:, :, 0],
        "detector_right": f["entry0"]["D11"]["Detector 3"]["data"][:, :, 0],
    }
    return detectors_data


electronic_noise = pyf.get_electronic_noise(
    file="institutes/ILL/instruments/D11/HEAD/mcstas/data/005708.nxs",
    acquisition_time=60,
    detector_names=["detector_left", "detector_right"],
    get_detector_data_fun=get_detector_data,
)


def compare_data(file, acquisition_time, mc_detectors, electronic_noise):
    # now check the compatibility with the data withing the background noise levels
    detector_data = get_detector_data(file)
    #
    assert np.sum(mc_detectors["detector_central"].Intensity) == pytest.approx(
        np.sum(detector_data["detector_central"]) / acquisition_time, rel=0.05
    )
    #
    for d in ["detector_left", "detector_right"]:
        assert np.sum(mc_detectors[d].Intensity) == pytest.approx(
            np.sum(detector_data[d]) / acquisition_time,
            abs=electronic_noise
            * detector_data[d].shape[0]
            * detector_data[d].shape[1]
            * 2,  # 2sigma of the noise
        ), d


def plot_data_sim(data, sim, tmp_path):
    fig = plt.figure(layout="constrained")
    axs = fig.subplot_mosaic(
        """
        lcr
        LCR
        """,
        # set the height ratios between the rows
        # height_ratios=[1, 1],
        # set the width ratios between the columns
        # width_ratios=[1, 2, 1],
    )
    for ax in axs:
        # print(ax)
        axs[ax].set_aspect("equal", adjustable="box")
    fig.suptitle("Vertically stacked subplots")
    if data is not None:
        axs["l"].imshow(
            data["detector_left"].transpose(),
            aspect="auto",
            cmap="seismic",
            origin="lower",
        )
        axs["c"].imshow(
            data["detector_central"].transpose(),
            aspect="auto",
            cmap="seismic",
            origin="lower",
        )
        axs["r"].imshow(
            data["detector_right"].transpose(),
            aspect="auto",
            cmap="seismic",
            origin="lower",
        )
    if sim is not None:
        axs["L"].imshow(
            sim["detector_left"].Ncount,
            aspect="auto",
            cmap="seismic",
            origin="lower",
        )
        axs["C"].imshow(
            sim["detector_central"].Ncount,
            aspect="auto",
            cmap="seismic",
            origin="lower",
        )
        axs["R"].imshow(
            sim["detector_right"].Ncount,
            aspect="auto",
            cmap="seismic",
            origin="lower",
        )

    # plt.show()
    plt.savefig(os.path.join(tmp_path, "plot.pdf"))
    print(f"Check file: {tmp_path}/plot.pdf")


def test_direct_beam(tmp_path):
    """Direct attenuated beam"""
    myinstrument.set_instrument_base_dir(str(tmp_path))
    myinstrument.sim_neutrons(1e6)
    myinstrument.master["lambda"] = 6 * ureg.angstrom
    myinstrument.master["detpos"] = 2 * ureg.m
    myinstrument.master["attenuator_index"] = 0
    myinstrument.master["collimation"] = 8 * ureg.m
    myinstrument.master["bs_index"] = 0
    myinstrument.sample_holder(
        material="quartz", shape="box", w=0.02, h=0.03, d=0.0135, th=0.00125
    )
    myinstrument.sample_shape("holder")
    myinstrument.set_sample_by_name("None")
    myinstrument.sample_holder(None, None)

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    for d in detector_names:
        assert np.sum(detector[d].Intensity) == 0, d
        assert np.sum(detector[d].Ncount) == 0, d


def test_direct_beam_noBS(tmp_path):
    myinstrument.set_instrument_base_dir(str(tmp_path))
    myinstrument.sim_neutrons(1e7)
    myinstrument.master["lambda"] = 6 * ureg.angstrom
    myinstrument.master["detpos"] = 2 * ureg.m
    myinstrument.master["attenuator_index"] = 6
    myinstrument.master["collimation"] = 8 * ureg.m
    myinstrument.master["bs_index"] = -1
    myinstrument.sample_holder(
        material="quartz", shape="box", w=0.02, h=0.03, d=0.0135, th=0.00125
    )
    myinstrument.sample_shape("holder")
    myinstrument.set_sample_by_name("None")
    myinstrument.sample_holder(None, None)

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    data_file = "institutes/ILL/instruments/D11/HEAD/mcstas/data/005708.nxs"
    plot_data_sim(get_detector_data(data_file), detector, tmp_path)

    # first check impact of any change in the simulation
    assert np.sum(detector["detector_central"].Intensity) == pytest.approx(
        8991, rel=1e-2
    )
    assert np.sum(detector["detector_central"].Ncount) == pytest.approx(
        40801, abs=1
    )  # pytest.approx(3977, rel=1e-2)

    for d in ["detector_left", "detector_right"]:
        assert np.sum(detector[d].Intensity) == 0, d
        assert np.sum(detector[d].Ncount) == 0, d

    # now check the compatibility with the data withing the background noise levels
    compare_data(
        file=data_file,
        acquisition_time=60,
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


def assert_detector(detector, dict_vals, rel=0.05):
    """
    dict_vals : { "detector_name": [intensity, ncount]}
    """

    for name, value in dict_vals.items():
        assert np.sum(detector[name].Intensity) == pytest.approx(
            value[0], abs=np.sum(detector[name].Error)
        ), name
        # pytest.approx(value[0], rel), name
        assert np.sum(detector[name].Ncount) == value[1], name


def test_empty_holder(tmp_path):
    myinstrument.set_instrument_base_dir(str(tmp_path))
    myinstrument.sim_neutrons(1e8)
    for calc in myinstrument.calculators:
        myinstrument.calculators[calc].settings(mpi=3)

    myinstrument.master["lambda"] = 6 * ureg.angstrom
    myinstrument.master["detpos"] = 2 * ureg.m
    myinstrument.master["attenuator_index"] = 0
    myinstrument.master["collimation"] = 8 * ureg.m
    myinstrument.master["bs_index"] = 0
    myinstrument.sample_holder(
        material="quartz", shape="box", w=0.02, h=0.03, d=0.0135, th=0.00125
    )
    myinstrument.sample_shape("holder")
    myinstrument.set_sample_by_name("None")

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    data_file = "institutes/ILL/instruments/D11/HEAD/mcstas/data/005721.nxs"
    plot_data_sim(get_detector_data(data_file), detector, tmp_path)

    # assert np.sum(detector["detector_central"].Intensity) == pytest.approx(989, 1)
    # assert np.sum(detector["detector_central"].Ncount) == 915  # with p_interact=1

    assert_detector(
        detector,
        {
            "detector_central": [1147, 915],
            "detector_left": [117, 119],
            "detector_right": [117, 144],
        },
    )

    #    assert np.sum(detector["detector_left"].Intensity) ==
    #    assert np.sum(detector["detector_left"].Ncount) == 119  # 70  # 7

    #    assert np.sum(detector["detector_right"].Intensity) == pytest.approx(179, 1)
    #    assert np.sum(detector["detector_right"].Ncount) == 144  # 10

    # now check the compatibility with the data withing the background noise levels

    compare_data(
        file=data_file,
        acquisition_time=60,
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


def test_sample1(tmp_path):
    myinstrument.set_instrument_base_dir(str(tmp_path))
    myinstrument.sim_neutrons(1e7)
    myinstrument.master["lambda"] = 6 * ureg.angstrom
    myinstrument.master["detpos"] = 2 * ureg.m
    myinstrument.master["attenuator_index"] = 0
    myinstrument.master["collimation"] = 8 * ureg.m
    myinstrument.master["bs_index"] = 0
    myinstrument.sample_holder(
        material="quartz", shape="box", w=0.02, h=0.03, d=0.0135, th=0.00125
    )
    myinstrument.sample_shape("holder")

    myinstrument.set_sample_by_name("qSq")
    myinstrument.master["sqw_file"] = (
        '"'
        + str(
            os.path.relpath(
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)), "data/simul_5711.sq"
                )
            )
        )
        + '"'
    )

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    assert np.sum(detector["detector_central"].Intensity) == pytest.approx(
        38165, rel=0.05
    )
    assert np.sum(detector["detector_central"].Ncount) == 152909

    assert np.sum(detector["detector_left"].Intensity) == pytest.approx(2036, rel=0.05)
    assert np.sum(detector["detector_left"].Ncount) == 9062

    assert np.sum(detector["detector_right"].Intensity) == pytest.approx(2036, rel=0.05)
    assert np.sum(detector["detector_right"].Ncount) == 9004

    # now check the compatibility with the data withing the background noise levels
    compare_data(
        file="institutes/ILL/instruments/D11/HEAD/mcstas/data/005711.nxs",
        acquisition_time=60,
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


"""
sys.exit(0)


def read_test(myinstrument, test_number, acquisition_time):
    simulation_dir = args.simulation_dir
    if simulation_dir is None:
        simulation_dir = set_test_dir(myinstrument.name, test_number) + "/OriginCalc"
    print("simulation dir: " + simulation_dir)
    myinstrument.set_test(test_number)
    metadata_list = load_metadata(simulation_dir)
    detectors_simulation = {}
    detectors_trueMC = {}
    for detector in metadata_list:
        if detector.component_name in detector_names or detector.component_name in [
            "PSD_sample"
        ]:
            monitor = load_monitor(detector, simulation_dir)
            detectors_simulation[monitor.name] = monitor.Intensity * acquisition_time
            detectors_trueMC[monitor.name] = monitor.Ncount

    return detectors_simulation, detectors_trueMC


def run_test(myinstrument, test_number, acquisition_time):

    myinstrument.set_test(test_number)

    myinstrument.run()
    data = myinstrument.output
    detectors = data.get_data()["data"]
    print(detectors)
    detectors_data = {}
    detectors_trueMC = {}
    for detector in detectors:
        if detector.name in detector_names or detector.name in ["PSD_sample"]:
            detectors_data[detector.name] = detector.Intensity * acquisition_time
            detectors_trueMC[detector.name] = detector.Ncount
    return detectors_data, detectors_trueMC




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
        strans = s[key].transpose()
        com_s = center_of_mass(strans, 0, strans.shape[0] - 1, 0, strans.shape[1] - 1)
        # del com_s
        # com_s = []
        # com_s.append(np.average(strans, 0))
        # com_s.append(np.average(strans, 1))
        com_d = center_of_mass(d[key], 0, d[key].shape[0] - 1, 0, d[key].shape[1] - 1)
        print(
            "{:18}: {:10.0f} {:10.0f} | {:10.2f} | {:10.0f} | {} | {}".format(
                key, dsum, ssum, ratio, mcsum, com_d, com_s
            )
        )


if args.simulation_dir is not None or args.plot:
    intensity, count = read_test(myinstrument, args.tests[-1], acquisition_time)
    data = data_test(myinstrument, args.tests[-1])
    compare(data, intensity, count)
    print(data["detector_central"].transpose().shape)
    fig, axs = plt.subplots(2, 3)
    fig.suptitle("Vertically stacked subplots")
    axs[0][0].imshow(
        data["detector_left"].transpose(),
        aspect="auto",
        cmap="seismic",
        origin="lower",
    )
    axs[0][1].imshow(
        data["detector_central"].transpose(),
        aspect="auto",
        cmap="seismic",
        origin="lower",
    )
    axs[0][2].imshow(
        data["detector_right"].transpose(),
        aspect="auto",
        cmap="seismic",
        origin="lower",
    )

    axs[1][0].imshow(
        intensity["detector_left"],
        aspect="auto",
        cmap="seismic",
        origin="lower",
    )
    axs[1][1].imshow(
        intensity["detector_central"],
        aspect="auto",
        cmap="seismic",
        origin="lower",
    )
    axs[1][2].imshow(
        intensity["detector_right"],
        aspect="auto",
        cmap="seismic",
        origin="lower",
    )
    plt.show()

    sys.exit(0)


darr = []
sarr = []
mcarr = []


for itest in args.tests:
    basedir = set_test_dir(instrument_name, itest)
    os.makedirs(basedir)
    myinstrument.set_instrument_base_dir(basedir)
    # if itest > 0:
    #    myinstrument.force_compile(True)
    darr.append(data_test(myinstrument, itest))
    intensity, count = run_test(myinstrument, itest, acquisition_time)
    sarr.append(intensity)
    mcarr.append(count)

    # print(sarr[0]["PSD_sample"])
    # sys.exit(0)

for itest in range(0, len(darr)):
    d = darr[itest]
    s = sarr[itest]
    mc = mcarr[itest]
    compare(d, s, mc)

# print("PSD_sample: ", center_of_mass(sarr[0]["PSD_sample"], 0, 100, 0, 100))
sys.exit(0)


axs[1][0].imshow(
    detectors_simulation["left"].Intensity,
    aspect="auto",
    cmap="seismic",
    origin="lower",
)
axs[1][1].imshow(
    detectors_simulation["central"].Intensity,
    aspect="auto",
    cmap="seismic",
    origin="lower",
)
axs[1][2].imshow(
    detectors_simulation["right"].Intensity,
    aspect="auto",
    cmap="seismic",
    origin="lower",
)
# cbar = fig.colorbar(detectors_simulation["central"])
# cbar.set_label("ZLabel", loc="top")
# cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# plt.colorbar(cax=cax)
plt.show()
f.close()

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
"""
