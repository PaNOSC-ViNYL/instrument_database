# Author: Shervin Nourbakhsh
# pytest  -k holder --comparison_type=stat
# pytest  -k holder --comparison_type=exact
# TODO:
# - [ ] add a test to verify attenuation values
# - [X] add noise estimate from data
#       def electronic_noise(data_direct_attenuated):
#
from pytest_check import check
from conftest import ureg, PlotHelper, ComparePlots, PdfPages, plt
import sys, os, pytest, math

# import the units
import h5py
import numpy as np

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
    detectors, time = get_NXdetectors(file)

    detectors_data = {}
    for d in detectors:
        detectors_data[d] = detectors[d]["data"][:, :, 0]

    return detectors_data, time


def stat(data):
    """This function calculates the mean and standard deviation of the values along x and y axes on 2D rectangular detectors.
    data   : numpy 2d array
    return : (mean_x, mean_y, std_x, std_y)

    McStasScript detector data should be transposed first
    """

    if len(data.shape) != 2:
        raise RuntimeError("wrong number of dimensions")

    def vv(data, axis):
        # array with the bin indexes
        x = np.array(range(0, data.T.shape[axis]))

        # cumulative weights along one axis
        w = np.sum(data, axis=axis)

        sumw = np.sum(w)
        if sumw != 0:
            sumwx = np.sum(x * w)
            sumwx2 = np.sum(x * x * w)
            m = sumwx / sumw
            std = sumwx2 / sumw - m * m
        else:
            m = 0
            std = 0
        return (m, np.sqrt(std))

    mx, stdx = vv(data, 0)
    my, stdy = vv(data, 1)

    return (mx, my, stdx, stdy)


# this can be generalized and put in the general functions
def compare_data(detector_data, acquisition_time, mc_detectors, electronic_noise):
    # now check the compatibility with the data withing the background noise levels

    with check:
        # Compare total counts
        # assert np.sum(mc_detectors["detector_central"].Intensity) == pytest.approx(
        #     np.sum(detector_data["detector_central"]) / acquisition_time, rel=0.05
        # )
        # #
        # for d in ["detector_left", "detector_right"]:
        #     assert np.sum(mc_detectors[d].Intensity) == pytest.approx(
        #         np.sum(detector_data[d]) / acquisition_time,
        #         abs=electronic_noise
        #         * detector_data[d].shape[0]
        #         * detector_data[d].shape[1]
        #         * 2,  # 2sigma of the noise
        #     ), d

        # Compare beam shape
        # just central since the rest is zero in this test
        detector = "detector_central"
        sd = stat(detector_data[detector])
        sm = stat(mc_detectors[detector].Intensity.T)
        assert sd[0] == pytest.approx(sm[0], 4), "Beam center - x axis"
        assert sd[1] == pytest.approx(sm[1], 4), "Beam center - y axis"
        assert sd[2] == sm[2], "Beam width - x axis"
        assert sd[3] == sm[3], "Beam width - y axis"


def plot_data_sim(data, sim, time, plot_settings):
    """
    data : ntuple (dict: {detector_name: numpy 2D array}, time (float))
    sim  : dict: {detector_name: McStasBinned}
    plot_settings : [other_settings, outpath]
    """
    pixel_sizes = plot_settings[0]
    tmp_path = plot_settings[1]

    ph1 = {}
    ph2 = {}
    for d in sim:
        ph1[d] = PlotHelper(sim[d].Intensity.T * time)
        ph1[d].set_pixel_size(pixel_sizes[d]["px"], pixel_sizes[d]["py"])
        if data is not None:
            ph2[d] = PlotHelper(data[d])
            ph2[d].set_pixel_size(pixel_sizes[d]["px"], pixel_sizes[d]["py"])

    compare_plots = ComparePlots(ph1, ph2, "Simulation", "Data")

    plot_positions = {
        "detector_left": "l",
        "detector_central": "c",
        "detector_right": "r",
    }

    fig1D = plt.figure(layout="constrained", dpi=92, figsize=(15, 15))
    axs = fig1D.subplot_mosaic([["lx", "cx", "rx"], ["ly", "cy", "ry"]])
    compare_plots.plot_marginals(plot_positions, axs)

    fig2D = plt.figure(layout="constrained", dpi=92, figsize=(15, 15))
    axs = fig2D.subplot_mosaic(
        [["l_data", "c_data", "r_data"], ["l_sim", "c_sim", "r_sim"]]
    )
    compare_plots.plot_2Ds(fig2D, plot_positions, axs, "r")

    with PdfPages(f"{tmp_path}/plot.pdf") as pdf:
        pdf.savefig(fig1D)
        pdf.savefig(fig2D)

    print(f"Check file: {tmp_path}/plot.pdf")


def assert_detector(detector, dict_vals, comparison_type):
    """
    dict_vals : { "detector_name": [intensity, ncount]}
    """
    with check:
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


@pytest.fixture(autouse=False)
def detector_names():
    """
    Name of the detectors
    """
    detector_names = ["detector_central", "detector_left", "detector_right"]
    return detector_names


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
def config_sample_monitor(config0):
    """
    Instrument configuration config0 with sample monitor of 2x2 cm2 at sample position
    """
    myinstrument = config0
    myinstrument.set_sample_by_name("monitor")
    myinstrument.sample.set_parameters(
        xwidth=0.02,
        yheight=0.02,
        zdepth=0.01,
    )

    # remove all simulation components after the sample monitor
    *_, calc = myinstrument.calculators

    while calc != myinstrument._calculator_with_sample.name:
        myinstrument.remove_calculator(calc)

    calc = myinstrument._calculator_with_sample
    comp = calc.get_last_component()

    while comp.name != myinstrument.sample.name:
        calc.remove_component(comp)
        comp = calc.get_last_component()

    return myinstrument


@pytest.fixture
def data_direct_attenuated():
    """Data with direct attenuated beam, no sample, no sample holder, no beamstop"""
    return get_NXdetector_data(
        os.path.abspath("institutes/ILL/instruments/D11/HEAD/mcstas/data/005708.nxs")
    )


@pytest.fixture
def data_empty_holder():
    """Data acquired with direct beam, with beam stop and empty sample holder"""
    return get_NXdetector_data(
        os.path.abspath("institutes/ILL/instruments/D11/HEAD/mcstas/data/005721.nxs")
    )


@pytest.fixture
def data_sample1():
    """Data acquired with quartz sample holder, and sample"""
    return get_NXdetector_data(
        os.path.abspath("institutes/ILL/instruments/D11/HEAD/mcstas/data/005711.nxs")
    )


@pytest.fixture
def electronic_noise(data_direct_attenuated, detector_names):
    """Return the electronic noise per pixel per second"""
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
            "px": detectors[d]["pixel_size_x"][0],
            "py": detectors[d]["pixel_size_y"][0],
        }

    return settings, tmp_path


# ------------------------------ tests
def test_plot_settings(plot_settings):
    """Trivial test printing the plot settings"""
    print(plot_settings)
    assert True == True


def test_name(set_instr):
    """Trivial test checking the instrument name"""
    myinstrument = set_instr
    assert myinstrument.name == "D11"


def test_detector_size(data_direct_attenuated):
    """Check the size of the detectors in a data file"""
    detectors, time, pixel_sizes = data_direct_attenuated
    assert detectors["detector_central"].shape == (256, 192)
    assert detectors["detector_left"].shape == (32, 256)
    assert detectors["detector_right"].shape == (32, 256)


def test_master_parameters(config0):
    """Check if the master parameters are all defined"""
    myinstr = config0

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

    # check that master parameters here defined are present in the simulation
    for par in master_parameters:
        assert par in myinstr.master

    # check that the test accounts for all the master parameters defined in the simulation
    for par in myinstr.master:
        assert par.name in master_parameters


def test_intensity_at_sample(config_sample_monitor, plot_settings):
    myinstrument = config_sample_monitor

    myinstrument.run()
    detectors = myinstrument.output.get_data()["data"]
    detector = [d for d in detectors if d.name in "sample_monitor"][0]

    # check simulation vs previous version
    assert np.sum(detector.Ncount) == 5092
    assert np.sum(detector.Intensity) == pytest.approx(4062845)
    ph = PlotHelper(detector.Intensity.T)
    ph.stats[0].mean = pytest.approx(49.68, abs=0.02)
    ph.stats[1].mean = pytest.approx(49.58, abs=0.02)
    ph.stats[0].std = pytest.approx(22.52, abs=0.02)
    ph.stats[1].std = pytest.approx(22.40, abs=0.02)

    detector = {d.name: d for d in detectors if d.name in "sample_monitor"}
    plot_sample_monitor(detector, plot_settings)

    # assert stat(detector.Intensity) == pytest.approx(
    #    (49.68, 49.58, 22.52, 22.40), abs=0.02
    # )


# @pytest.mark.xfail(
#    True,
#    reason="The beam width should be tuned",
# )
def test_direct_beam_noBS(
    config0,
    data_direct_attenuated,
    electronic_noise,
    plot_settings,
    detector_names,
    comparison_type,
):
    data, time = data_direct_attenuated
    myinstrument = config0
    myinstrument.settings(ncount=1e7)
    myinstrument.master["attenuator_index"] = 6
    myinstrument.master["bs_index"] = -1

    mycalc = myinstrument.calculators["OriginCalc"]
    source = mycalc.get_component("VCS")

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    plot_data_sim(data, detector, time, plot_settings)

    # first check impact of any change in the simulation
    # assert_detector(
    #     detector,
    #     {
    #         "detector_central": [9304, 48421],
    #         "detector_left": [0, 0],
    #         "detector_right": [0, 0],
    #     },
    #     comparison_type,
    # )

    # now check the compatibility with the data withing the background noise levels
    compare_data(
        detector_data=data_direct_attenuated[0],
        acquisition_time=data_direct_attenuated[1],
        electronic_noise=electronic_noise,
        mc_detectors=detector,
    )


def test_direct_beam_withBS(config0, data_direct_attenuated, plot_settings):  # ):

    """Direct attenuated beam"""
    myinstrument = config0
    myinstrument.settings(ncount=1e6)

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    detector = {d.name: d for d in detectors if d.name in detector_names}

    for d in detector_names:
        assert np.sum(detector[d].Intensity) == 0, d
        assert np.sum(detector[d].Ncount) == 0, d

    plot_data_sim(data_direct_attenuated, detector, plot_settings)


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

    plot_data_sim(data_empty_holder, detector, plot_settings)

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

    plot_data_sim(data_sample1, detector, plot_settings)

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


def test_diagram(config0, tmp_path):
    myinstrument = config0
    myinstrument.settings(ncount=1e5)

    mycalc = myinstrument.calculators["OriginCalc"]
    # myinstrument.run()
    fig = mycalc.show_diagram(analysis=True, variable=None)
    fig.set_size_inches(10, 30)
    fig.savefig(os.path.join(tmp_path, "diagram.pdf"))
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

    # remove components after VS
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
