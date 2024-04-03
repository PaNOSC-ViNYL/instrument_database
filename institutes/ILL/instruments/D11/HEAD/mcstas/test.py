# Author: Shervin Nourbakhsh
# pytest  -k holder --comparison_type=stat
# pytest  -k holder --comparison_type=exact

from conftest import ureg
import sys, os, pytest, math

# import the units
import h5py
import numpy as np
import matplotlib.pyplot as plt

from D11 import def_instrument

# ------------------------------ helper functions
def get_parameters_from_NX(file):
    """Read the instrument parameters from the Nexus file"""
    # Read NeXus
    f = h5py.File(file, "r")
    # create a dict with all the detectors
    pars = {
        #        "energy": f["entry0"]["Ei"][0],
        #        "a2": f["entry0"]["a2"]["target_value"][0],
        #        "chopper_rpm": f["entry0"]["BC1"]["rotation_speed"][0],
        #        "chopper_ratio": f["entry0"]["FC"]["ratio"][0],
        #        "Efoc": f["entry0"]["Detector"]["elasticpeak"][0],
    }


def get_NXdetectors(file):
    """Function that returns the NXdetector object"""
    # Read NeXus
    f = h5py.File(file, "r")
    # create a dict with all the detectors
    detectors = {
        "detector_left": f["entry0"]["D11"]["Detector 2"],
        "detector_central": f["entry0"]["D11"]["Detector 1"],
        "detector_right": f["entry0"]["D11"]["Detector 3"],
    }
    acquisition_time = f["entry0"]["duration"][0]
    return detectors, acquisition_time


def get_NXdetector_data(file):
    """Function that returns a dict of detectors with their data as ndarray"""
    d, time = get_NXdetectors(file)

    detectors_data = {
        "detector_left": d["detector_left"]["data"][:, :, 0],
        "detector_central": d["detector_central"]["data"][:, :, 0],
        "detector_right": d["detector_right"]["data"][:, :, 0],
    }
    return detectors_data, time


def compare_data(detector_data, acquisition_time, mc_detectors, electronic_noise):
    # now check the compatibility with the data withing the background noise levels
    #
    #
    for d in ["detector_left", "detector_right"]:
        assert np.sum(mc_detectors[d].Intensity) == pytest.approx(
            np.sum(detector_data[d]) / acquisition_time,
            abs=electronic_noise
            * detector_data[d].shape[0]
            * detector_data[d].shape[1]
            * 2,  # 2sigma of the noise
        ), d


def stat(data):
    if len(data.shape) != 2:
        raise RuntimeError("wrong number of dimensions")

    def vv(data, axis):
        # array with the bin indexes
        x = np.array(range(0, data.shape[axis]))

        # cumulative weights along one axis
        w = np.sum(data, axis=axis)

        sumw = np.sum(w)
        sumwx = np.sum(x * w)
        sumwx2 = np.sum(x * x * w)

        m = sumwx / sumw
        std = sumwx2 / sumw - m * m
        return (m, np.sqrt(std))

    mx, stdx = vv(data, 0)
    my, stdy = vv(data, 1)

    return (mx, my, stdx, stdy)


def plot_sample_monitor(sim, plot_settings):
    tmp_path = plot_settings[1]
    settings = plot_settings[0]

    fig = plt.figure(layout="constrained", figsize=(15, 10))
    axs = fig.subplot_mosaic("""IN""")

    immap = {}
    immap["I"] = axs["I"].imshow(sim.Intensity.T, **settings["sample_monitor"])
    immap["N"] = axs["N"].imshow(sim.Ncount.T, **settings["sample_monitor"])
    fig.colorbar(ax=axs["I"], mappable=immap["I"])
    fig.colorbar(ax=axs["N"], mappable=immap["N"])

    plt.savefig(os.path.join(tmp_path, "sample_monitor.pdf"))


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


def plot_data_sim(data, sim, plot_settings):
    tmp_path = plot_settings[1]
    settings = plot_settings[0]
    fig = plt.figure(layout="constrained", dpi=92, figsize=(15, 15))
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
    # fig.suptitle("Vertically stacked subplots")

    immap = {}
    if data is not None:
        imshow = lambda ax, name: ax.imshow(data[name].transpose(), **settings[name])
        imshow(axs["l"], "detector_left")
        immap["data"] = imshow(axs["c"], "detector_central")
        imshow(axs["r"], "detector_right")

    if sim is not None:
        imshow = lambda ax, name: ax.imshow(sim[name].Intensity, **settings[name])
        imshow(axs["L"], "detector_left")
        immap["sim"] = imshow(axs["C"], "detector_central")
        imshow(axs["R"], "detector_right")

    fig.colorbar(ax=axs["r"], mappable=immap["data"])
    fig.colorbar(ax=axs["R"], mappable=immap["sim"])

    # plt.show()
    plt.savefig(os.path.join(tmp_path, "plot.pdf"))
    print(f"Check file: {tmp_path}/plot.pdf")

    # def test_plot(data_direct_attenuated, tmp_path, detector_sizes):
    #    plot_data_sim(
    #        data_direct_attenuated, data_direct_attenuated, tmp_path, detector_sizes
    #    )
    return
    fig = plt.figure(layout="constrained", figsize=(15, 15))
    axs = fig.subplot_mosaic(
        """
        xyt
        XYT
        """
    )

    axs["x"].scatter(
        np.linspace(0, sim[name].I.shape[0], sim[name].I.shape[0]),
        np.sum(data[name], axis=(1, 2)),
    )
    axs["y"].scatter(
        np.linspace(0, sim[name].I.shape[1], sim[name].I.shape[1]),
        np.sum(data[name], axis=(0, 2)),
    )
    axs["t"].scatter(
        np.linspace(0, sim[name].I.shape[2], sim[name].I.shape[2]),
        np.sum(data[name], axis=(0, 1)),
    )

    axs["X"].scatter(
        np.linspace(0, sim[name].I.shape[0], sim[name].I.shape[0]),
        np.sum(sim[name].I, axis=(1, 2)),
    )
    axs["Y"].scatter(
        np.linspace(0, sim[name].I.shape[1], sim[name].I.shape[1]),
        np.sum(sim[name].I, axis=(0, 2)),
    )
    axs["T"].scatter(
        np.linspace(0, sim[name].I.shape[2], sim[name].I.shape[2]),
        np.sum(sim[name].I, axis=(0, 1)),
    )
    plt.savefig(os.path.join(tmp_path, "plot2.pdf"))
    print(f"Check file: {tmp_path}/plot.pdf")


def assert_detector(detector, dict_vals, comparison_type):
    """
    dict_vals : { "detector_name": [intensity, ncount]}
    """

    for name, value in dict_vals.items():
        if comparison_type == "exact":
            assert np.sum(detector[name].Intensity) == pytest.approx(
                value[0], abs=1
            ), "SIM {} intensity: {} != {}".format(
                name, detector[name].Intensity, value[0]
            )
            assert np.sum(detector[name].Ncount) == value[1], name
        elif comparison_type == "stat":
            assert np.sum(detector[name].Intensity) == pytest.approx(
                value[0], abs=np.sum(detector[name].Error)
            ), name
            assert np.sum(detector[name].Ncount) == pytest.approx(
                value[1], abs=math.sqrt(value[1])
            ), name


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


# ------------------------------ fixtures
@pytest.fixture
def flavour():
    """
    Default flavour to test
    """
    # return "nosection"
    return "simplefull"


@pytest.fixture
def mpi():
    """Setting default mpi"""
    return 8


@pytest.fixture
def set_instr(flavour, mpi, tmp_path):
    """Setup the instrument with common settings"""
    myinstrument = def_instrument(flavour)
    myinstrument.set_instrument_base_dir(str(tmp_path))
    myinstrument.settings(run_path=(str(tmp_path)), mpi=mpi, seed=654321)
    return myinstrument


@pytest.fixture
def detector_names():
    """
    Instrument detectors
    """
    return ["detector_central", "detector_left", "detector_right"]


@pytest.fixture
def config0(set_instr):
    """Base configuration for tests: direct beam with beam stop"""
    myinstrument = set_instr
    myinstrument.master["lambda"] = 6 * ureg.angstrom
    myinstrument.master["detpos"] = 2 * ureg.m
    myinstrument.master["attenuator_index"] = 0
    myinstrument.master["collimation"] = 8 * ureg.m
    myinstrument.master["bs_index"] = 0
    myinstrument.set_sample_by_name("None")
    myinstrument.sample_holder(None, None)
    return myinstrument


@pytest.fixture
def data_direct_attenuated():
    return get_NXdetector_data(
        os.path.abspath("institutes/ILL/instruments/D11/HEAD/mcstas/data/005708.nxs")
    )


@pytest.fixture
def data_empty_holder():
    return get_NXdetector_data(
        os.path.abspath("institutes/ILL/instruments/D11/HEAD/mcstas/data/005721.nxs")
    )


@pytest.fixture
def data_sample1():
    return get_NXdetector_data(
        os.path.abspath("institutes/ILL/instruments/D11/HEAD/mcstas/data/005711.nxs")
    )


@pytest.fixture
def electronic_noise(data_direct_attenuated):
    detector_names = ["detector_left", "detector_right"]
    detector_data, acquisition_time = data_direct_attenuated
    I = 0
    n = 0
    for name in detector_names:
        d = detector_data[name]
        I = I + np.sum(d)
        n = n + d.shape[0] * d.shape[1]

    # I should use a Poissonian instead of a Gaus distribution :-)
    return I / n / acquisition_time


@pytest.fixture
def plot_settings(detector_names, tmp_path):
    file = "institutes/ILL/instruments/D11/HEAD/mcstas/data/005708.nxs"
    # Read NeXus
    detectors, time = get_NXdetectors(file)

    settings = {}
    for d in detector_names:
        settings[d] = {
            "cmap": "hsv",
            "origin": "lower",
            "extent": [
                0,
                detectors[d]["pixel_size_x"][0]
                * detectors[d]["data"][:, :, 0].shape[0],
                0,
                detectors[d]["pixel_size_y"][0]
                * detectors[d]["data"][:, :, 0].shape[1],
            ],
        }
    return settings, tmp_path


# ------------------------------ tests
def test_comparison_type(comparison_type):
    assert comparison_type in ["exact", "stat"]


def test_plot_settings(plot_settings):
    print(plot_settings)
    assert True == True


def test_name(set_instr):
    myinstrument = set_instr
    assert myinstrument.name == "D11"
    # print(myinstrument.calculators["OriginCalc"].get_component("VCS").dist)
    # assert myinstrument.calculators["OriginCalc"].get_component("VCS").dist == 61.97


def test_detector_size(data_direct_attenuated):
    detectors, time = data_direct_attenuated
    assert detectors["detector_central"].shape == (256, 192)
    assert detectors["detector_left"].shape == (32, 256)
    assert detectors["detector_right"].shape == (32, 256)


def test_master_parameters(config0):
    myinstr = config0

    # print(help(myinstr.master))
    master_parameters = [
        "lambda",
        "collimation",
        "disk1_index",
        "disk2_index",
        "disk3_index",
        "disk4_index",
        "disk5_index",
        "detpos",
        "attenuator_index",
        "bs_x",
        "bs_y",
        "bs_index",
    ]
    # print(myinstr.master)
    for par in master_parameters:
        assert par in myinstr.master


def test_direct_beam_noBS(
    config0,
    detector_names,
    comparison_type,
    data_direct_attenuated,
    electronic_noise,
    plot_settings,
):
    myinstrument = config0
    myinstrument.settings(ncount=1e7)
    myinstrument.master["attenuator_index"] = 6
    myinstrument.master["bs_index"] = -1

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    plot_data_sim(data_direct_attenuated[0], detector, plot_settings)

    # first check impact of any change in the simulation
    assert_detector(
        detector,
        {
            "detector_central": [9304, 48421],
            "detector_left": [0, 0],
            "detector_right": [0, 0],
        },
        comparison_type,
    )

    # now check the compatibility with the data withing the background noise levels
    compare_data(
        detector_data=data_direct_attenuated[0],
        acquisition_time=data_direct_attenuated[1],
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


def test_direct_beam_withBS(
    config0, detector_names, data_direct_attenuated, plot_settings
):  # ):

    """Direct attenuated beam"""
    myinstrument = config0
    myinstrument.settings(ncount=1e6)

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    for d in detector_names:
        assert np.sum(detector[d].Intensity) == 0, d
        assert np.sum(detector[d].Ncount) == 0, d

    plot_data_sim(data_direct_attenuated[0], detector, plot_settings)


# def test_attenuation(config0, detector_names, comparison_type):
#     myinstrument = config0
#     myinstrument.settings(ncount=1e6)
#     myinstrument.master["attenuator_index"] = 6
#     myinstrument.master["bs_index"] = -1

#     myinstrument.run()

#     detectors = myinstrument.output.get_data()["data"]
#     detector = {d.name: d for d in detectors if d.name in detector_names}

#     myinstrument.master["attenuator_index"] = 6
#     myinstrument.run()
#     attdetectors = myinstrument.output.get_data()["data"]
#     attdetector = {d.name: d for d in detectors if d.name in detector_names}


@pytest.mark.xfail  # (True, reason="need to debug the difference between data and simulation", True)
def test_empty_holder(
    config0,
    detector_names,
    comparison_type,
    data_empty_holder,
    electronic_noise,
    plot_settings,
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

    plot_data_sim(data_empty_holder[0], detector, plot_settings)

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
        data_empty_holder[0],
        acquisition_time=data_empty_holder[1],
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


@pytest.mark.xfail  # (True, reason="need to debug the difference between data and simulation", True)
def test_sample1(
    config0,
    detector_names,
    comparison_type,
    data_sample1,
    plot_settings,
    electronic_noise,
):
    myinstrument = config0
    myinstrument.settings(ncount=1e9)
    myinstrument.master["bs_index"] = 1
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

    plot_data_sim(data_sample1[0], detector, plot_settings)

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
        data_sample1[0],
        acquisition_time=data_sample1[1],
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


def test_diagram(config0):
    myinstrument = config0
    myinstrument.settings(ncount=1e5)

    mycalc = myinstrument.calculators["OriginCalc"]
    # myinstrument.run()
    mycalc.show_diagram(analysis=True, variable=None)
    assert True == True


def test_diagnostics(config0, tmp_path):
    myinstrument = config0

    mycalc = myinstrument.calculators["OriginCalc"]
    # mycalc.get_component("VCS").focus_yh = 0.04
    import mcstasscript as ms

    # fermi = mycalc.get_component("fermi_chopper")

    diag = ms.Diagnostics(mycalc)
    diag.settings(ncount=1e6, suppress_output=False, mpi=1)
    diag.show_settings()
    diag.clear_points()

    for component in diag.instr.component_list:
        if component.component_name in ["Progress_bar", "Arm"]:
            continue
        if component.component_name in ["Beamstop"]:
            diag.add_point(before=component.name)
            continue
        if component.name in ["VCS", "HCS"]:
            diag.add_point(after=component.name)
            continue
        if component.component_name in [
            "Disk_Chopper",
            "Slit",
            "Monochromator_curved",
            "Fermi_chopper",
        ]:
            diag.add_point(before=component.name)
            diag.add_point(after=component.name)
        else:
            # print("#"+component.name+"#")
            diag.add_point(before=component.name)

    if mycalc.name == "SampleCalc":
        diag.clear_points()
        diag.add_point(after="Vin")
        diag.add_point(before="sqw")
        diag.add_point(after="Sample_Out")
        diag.add_point(before="after_sample_slit")
        diag.add_point(after="after_sample_slit")

    diag.run()

    diag.clear_views()
    diag.add_view("x", bins=50)  # , limits=[-0.04,0.04])
    diag.add_view("y", bins=50)  # , same_scale=False, limits=[-0.07,0.07])
    # diag.add_view("x","y",bins=[30,30])
    # diag.add_view("e", same_scale=False)
    # diag.add_view("l", same_scale=False)
    diag.add_view("t", same_scale=False)
    diag.add_view("t", "x", same_scale=False, log=True)
    diag.add_view("t", "y", same_scale=False)
    diag.add_view("x", "y", same_scale=False, log=True)
    # diag.add_view(
    #    "dx", "dy", same_scale=False, log=True
    # )  # , limits=[-50,-10],left_min=-45,right_lim=-20,bins=100)
    # diag.add_view("dx","x",same_scale=False,log=True) #, limits=[-50,-10],left_min=-45,right_lim=-20,bins=100)
    # diag.add_view("dx","y",same_scale=False,log=True) #, limits=[-50,-10],left_min=-45,right_lim=-20,bins=100)

    diag.add_view("vz", same_scale=True)
    print(diag)
    fig = diag.plot()
    fig.savefig(os.path.join(tmp_path, "diagnostics.pdf"))

    assert True == True


def test_optim_source_focus(config0, tmp_path):
    """Make sure that the distribution of lambda si wider than that coming out of the velocity selector
    and that the x and y distribution are almost flat"""
    import mcstasscript as ms

    myinstrument = config0
    myinstrument.set_sample_by_name("monitor")
    myinstrument.sample.set_parameters(
        xwidth=0.04,
        yheight=0.04,
        zdepth=0.01,
    )

    mycalc = myinstrument.calculators["OriginCalc"]

    last = mycalc.get_component("sg30")
    ilast = mycalc.component_list.index(last)
    components_to_remove = reversed(mycalc.component_list[ilast + 1 :])
    for component in components_to_remove:
        mycalc.remove_component(component)
    # mycalc.remove_component("detector")

    def create_diag(mycalc):
        diag = ms.Diagnostics(mycalc)
        diag.settings(ncount=1e6, suppress_output=False, mpi=1)
        # diag.show_settings()
        diag.clear_points()

        diag.clear_views()
        diag.add_view("x", bins=100)  # , limits=[-0.04,0.04])
        diag.add_view("y", bins=100)  # , same_scale=False, limits=[-0.07,0.07])
        # diag.add_view("e", same_scale=False)
        diag.add_view("l", same_scale=False)
        # diag.add_view("t", same_scale=False)
        # diag.add_view("t", "x", same_scale=False, log=True)
        # diag.add_view("t", "y", same_scale=False)
        diag.add_view("x", "y", same_scale=False, log=True)
        diag.add_view("vz", same_scale=True)

        for component in mycalc.component_list:
            if component.name in ["sample_monitor", "Dolores", "sg30"]:
                diag.add_point(before=component.name)

        return diag

    source = mycalc.get_component("VCS")
    results = []
    X = np.arange(0.02, 0.045, 0.005)
    # X = [HCS.focus_xw]
    Y = X
    L = np.arange(0.04, 0.06, 0.01)
    for l in L:
        mycalc.parameters["dlambda"] = l
        diag = create_diag(mycalc)
        diag.run()
        # print(diag.instr.output.get_data()["data"][0])
        # print(diag.instr.output.get_data()["data"][0].make_1d("l"))
        fig = diag.plot()
        fig.savefig(os.path.join(tmp_path, "optim-l_{}.pdf".format(l)))

    for x in X:
        y = 0
        source.focus_xw = x
        diag = create_diag(mycalc)
        diag.run()
        fig = diag.plot()
        fig.savefig(os.path.join(tmp_path, "optim-x_{}.pdf".format(x)))

    for y in Y:
        x = 0
        source.focus_yh = y
        diag = create_diag(mycalc)
        diag.run()
        #########################################################################
        # detectors = diag.instr.output.get_data()["data"]                      #
        # sample_monitor = [d for d in detectors if d.name == "sample_monitor"] #
        #                                                                       #
        # results.append(                                                       #
        #     [                                                                 #
        #         x,                                                            #
        #         y,                                                            #
        #         np.sum(sample_monitor[0].Ncount),                             #
        #         np.sum(sample_monitor[0].Intensity),                          #
        #         stat(sample_monitor[0].Ncount),                               #
        #         stat(sample_monitor[0].Intensity),                            #
        #     ]                                                                 #
        # )                                                                     #
        #########################################################################
        # print(diag)
        fig = diag.plot()
        fig.savefig(os.path.join(tmp_path, "optim-y_{}.pdf".format(y)))
    return
    #    assert np.sum(sample_monitor[0].Ncount) == 8373
    #    assert np.sum(sample_monitor[0].Intensity) == pytest.approx(4152070)
    # print(results)
    for r in results:
        print(
            "{:.2f},{:.2f}\t{:.0f} - {:.0f}\t{}\t{}".format(
                r[0], r[1], r[2], r[3], r[4][2:], r[5][2:]
            )
        )
    assert True == True
