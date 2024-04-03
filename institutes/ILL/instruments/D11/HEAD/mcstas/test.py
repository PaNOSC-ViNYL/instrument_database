# pytest  -k holder --comparison_type=stat
# pytest  -k holder --comparison_type=exact
# python institutes/ILL/instruments/D11/HEAD/mcstas/test.py | tail -12 >  institutes/ILL/instruments/D11/HEAD/mcstas/test/test.log; git diff institutes/ILL/instruments/D11/HEAD/mcstas/test/test.log

from D11 import def_instrument
import sys, os

import pytest
import pytest_functions as pyf


# import the units
import pint
import numpy as np
import h5py

import matplotlib.pyplot as plt

# from mcstasscript.interface.functions import load_metadata, load_monitor
@pytest.fixture
def flavour():
    """
    Default flavour to test
    """
    return "simplefull"


@pytest.fixture
def set_instr(flavour, mpi, tmp_path):
    """Setup the instrument with common settings"""
    myinstrument = def_instrument(flavour)
    myinstrument.set_seed(654321)
    myinstrument.set_instrument_base_dir(str(tmp_path))
    for calc in myinstrument.calculators:
        myinstrument.calculators[calc].settings(mpi=mpi)

    return myinstrument


@pytest.fixture
def detector_names():
    """
    Instrument detectors
    """
    return ["detector_central", "detector_left", "detector_right"]


@pytest.fixture
def config0(set_instr):
    myinstrument = set_instr
    myinstrument.master["lambda"] = 6 * ureg.angstrom
    myinstrument.master["detpos"] = 2 * ureg.m
    myinstrument.master["attenuator_index"] = 0
    myinstrument.master["collimation"] = 8 * ureg.m
    myinstrument.master["bs_index"] = 0
    myinstrument.set_sample_by_name("None")
    myinstrument.sample_holder(None, None)
    return myinstrument


ureg = pint.get_application_registry()

# setting the base directory for the simulation output


def test_comparison_type(comparison_type):
    assert comparison_type in ["exact", "phys"]


def test_name(set_instr):
    myinstrument = set_instr
    assert myinstrument.name == "D11"


def get_detector_data(file):
    """This can be a member function of the class"""
    # file = myinstrument.test_datafile(test_number)
    # Read NeXus
    f = h5py.File(file, "r")

    detectors_data = {
        "detector_left": f["entry0"]["D11"]["Detector 2"]["data"][:, :, 0],
        "detector_central": f["entry0"]["D11"]["Detector 1"]["data"][:, :, 0],
        "detector_right": f["entry0"]["D11"]["Detector 3"]["data"][:, :, 0],
    }
    return detectors_data


@pytest.fixture
def data_direct_attenuated():
    acquisition_time = 60
    return get_detector_data(
        "institutes/ILL/instruments/D11/HEAD/mcstas/data/005708.nxs"
    )


@pytest.fixture
def data_empty_holder():
    return get_detector_data(
        "institutes/ILL/instruments/D11/HEAD/mcstas/data/005721.nxs"
    )


@pytest.fixture
def data_sample1():
    return get_detector_data(
        "institutes/ILL/instruments/D11/HEAD/mcstas/data/005711.nxs"
    )


@pytest.fixture
def electronic_noise(data_direct_attenuated):
    detector_names = ["detector_left", "detector_right"]
    detector_data = data_direct_attenuated
    acquisition_time = 60
    I = 0
    n = 0
    for name in detector_names:
        d = detector_data[name]
        I = I + np.sum(d)
        n = n + d.shape[0] * d.shape[1]

    # I should use a Poissonian instead of a Gaus distribution :-)
    return I / n / acquisition_time


# this can be generalized and put in the general functions
def compare_data(detector_data, acquisition_time, mc_detectors, electronic_noise):
    # now check the compatibility with the data withing the background noise levels
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


def test_direct_beam(config0, detector_names):
    """Direct attenuated beam"""
    myinstrument = config0
    myinstrument.sim_neutrons(1e6)

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    for d in detector_names:
        assert np.sum(detector[d].Intensity) == 0, d
        assert np.sum(detector[d].Ncount) == 0, d


import math


def assert_detector(detector, dict_vals, comparison_type):
    """
    dict_vals : { "detector_name": [intensity, ncount]}
    """

    for name, value in dict_vals.items():
        if comparison_type == "exact":
            assert np.sum(detector[name].Intensity) == pytest.approx(
                value[0], abs=1
            ), name
            assert np.sum(detector[name].Ncount) == value[1], name
        elif comparison_type == "stat":
            assert np.sum(detector[name].Intensity) == pytest.approx(
                value[0], abs=np.sum(detector[name].Error)
            ), name
            assert np.sum(detector[name].Ncount) == pytest.approx(
                value[1], abs=math.sqrt(value[1])
            ), name


def test_direct_beam_noBS(
    config0,
    detector_names,
    tmp_path,
    comparison_type,
    data_direct_attenuated,
    electronic_noise,
):
    myinstrument = config0
    myinstrument.sim_neutrons(1e7)
    myinstrument.master["attenuator_index"] = 6
    myinstrument.master["bs_index"] = -1

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    plot_data_sim(data_direct_attenuated, detector, tmp_path)

    # first check impact of any change in the simulation
    assert_detector(
        detector,
        {
            "detector_central": [8929, 40588],
            "detector_left": [0, 0],
            "detector_right": [0, 0],
        },
        comparison_type,
    )

    # now check the compatibility with the data withing the background noise levels
    compare_data(
        detector_data=data_direct_attenuated,
        acquisition_time=60,
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


def test_attenuation(config0, detector_names, comparison_type):
    myinstrument = config0
    myinstrument.sim_neutrons(1e6)
    myinstrument.master["attenuator_index"] = 6
    myinstrument.master["bs_index"] = -1

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    myinstrument.master["attenuator_index"] = 6
    myinstrument.run()
    attdetectors = myinstrument.output.get_data()["data"]
    attdetector = {d.name: d for d in detectors if d.name in detector_names}


@pytest.mark.xfail  # (True, reason="need to debug the difference between data and simulation", True)
def test_empty_holder(
    config0, detector_names, tmp_path, comparison_type, data_empty_holder
):
    myinstrument = config0
    myinstrument.sim_neutrons(1e8)

    myinstrument.sample_holder(
        material="quartz", shape="box", w=0.02, h=0.03, d=0.0135, th=0.00125
    )
    myinstrument.sample_shape("holder")
    myinstrument.set_sample_by_name("None")

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    plot_data_sim(data_empty_holder, detector, tmp_path)

    # assert np.sum(detector["detector_central"].Intensity) == pytest.approx(989, 1)
    # assert np.sum(detector["detector_central"].Ncount) == 915  # with p_interact=1

    assert_detector(
        detector,
        {
            "detector_central": [1135, 893],
            "detector_left": [205, 142],
            "detector_right": [174, 146],
        },
        comparison_type,
    )

    # now check the compatibility with the data withing the background noise levels
    compare_data(
        data_empty_holder,
        acquisition_time=60,
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


def test_sample1(config0, detector_names, tmp_path, comparison_type, data_sample1):
    myinstrument = config0
    myinstrument.sim_neutrons(1e9)
    myinstrument.master["lambda"] = 6 * ureg.angstrom
    myinstrument.master["detpos"] = 2 * ureg.m
    myinstrument.master["attenuator_index"] = 0
    myinstrument.master["collimation"] = 8 * ureg.m
    myinstrument.master["bs_index"] = 3
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

    plot_data_sim(data_sample1, detector, tmp_path)

    assert_detector(
        detector,
        {
            "detector_central": [37638, 151577],
            "detector_left": [2004, 8807],
            "detector_right": [2079, 9036],
        },
        comparison_type,
    )

    # now check the compatibility with the data withing the background noise levels
    compare_data(
        data_sample1,
        acquisition_time=60,
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


"""
attenuation_values = [
    8.325,  # attenuator 1
    26.21,  # attenuator 2
    72.23,  # attenuator 3
    216.5,  # attenuator 1+2
    594.6,  # attenuator 1+3
    1702,  # attenuator 2+3
    13480,  # attenuator 1+2+3
]
"""
