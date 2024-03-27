# pytest  -k holder --comparison_type=stat
# pytest  -k holder --comparison_type=exact
#
# - [ ] flux @ sample (2cm x 2cm) = 5e5 n/cm/s
# - [ ] vanadium sample: r=5mm h=40mm th=1mm
# - [ ] plot of 3D data for comparison
# - [ ] profile of 3D data for comparison
from conftest import ureg
import sys, os, pytest, math

# import the units
import h5py
import numpy as np
import matplotlib.pyplot as plt

from Panther import def_instrument


def get_distance(a, b, instr):
    a = instr.get_component(a)
    b = instr.get_component(b)
    ia = instr.component_list.index(a)
    ib = instr.component_list.index(b)
    d = None
    z = 0
    for c in reversed(instr.component_list[ia:ib]):
        if d is not None and c.name != d:
            continue
        print(c)
        z = z + float(c.AT_data[2])
        d = c.AT_reference
    return z


# ------------------------------ helper functions
def get_parameters_from_NX(file):
    """Read the instrument parameters from the Nexus file"""
    # Read NeXus
    f = h5py.File(file, "r")
    # create a dict with all the detectors
    pars = {
        "energy": f["entry0"]["Ei"][0],
        "a2": f["entry0"]["a2"]["target_value"][0],
        "chopper_rpm": f["entry0"]["BC1"]["rotation_speed"][0],
        "chopper_ratio": f["entry0"]["FC"]["ratio"][0],
        "Efoc": f["entry0"]["Detector"]["elasticpeak"][0],
    }


def get_NXdetectors(file):
    """Function that returns the NXdetector object"""
    # Read NeXus
    f = h5py.File(file, "r")
    # create a dict with all the detectors
    detectors = {
        "detector": f["entry0"]["instrument"]["Detector"],
    }
    acquisition_time = f["entry0"]["duration"][0]
    return detectors, acquisition_time


def get_NXdetector_data(file):
    """Function that returns a dict of detectors with their data as ndarray"""
    d, time = get_NXdetectors(file)

    detectors_data = {"detector": d["detector"]["data"]}
    return detectors_data, time


def compare_data(detector_data, acquisition_time, mc_detectors, electronic_noise):
    # now check the compatibility with the data withing the background noise levels
    #
    #
    for d in ["detector"]:
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


def plot_data_sim(data, sim, plot_settings):
    tmp_path = plot_settings[1]
    settings = plot_settings[0]

    fig = plt.figure(layout="constrained", dpi=92, figsize=(15, 15))
    axs = fig.subplot_mosaic(
        """
        c
        C
        """,
        # set the height ratios between the rows
        # height_ratios=[1, 1],
        # set the width ratios between the columns
        # width_ratios=[1, 2, 1],
    )
    # fig.suptitle("Vertically stacked subplots")
    name = next(iter(data.keys()))  # first, and unique, detector name
    immap = {}
    if data is not None:
        # print(stat(data[name][:, :, 100].T))
        imshow = lambda ax, t: ax.imshow(data[name][:, :, t].T, **settings[name])
        immap["data"] = imshow(axs["c"], 100)

    if sim is not None:
        # print(stat(np.flip(np.sum(sim[name].I, axis=0), 0)))
        imshow = lambda ax, t: ax.imshow(
            np.sum(sim[name].I, axis=-1).T, **settings[name]
        )
        immap["sim"] = imshow(axs["C"], 100)

    fig.colorbar(ax=axs["c"], mappable=immap["data"])
    fig.colorbar(ax=axs["C"], mappable=immap["sim"])

    # plt.show()
    plt.savefig(os.path.join(tmp_path, "plot.pdf"))

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
def detector_names():
    """
    Instrument detectors
    """
    return ["detector"]


@pytest.fixture
def flavour():
    """
    Default flavour to test
    """
    return "nosection"
    # return "quicknosection"


@pytest.fixture
def mpi():
    """Setting default mpi"""
    return 7


@pytest.fixture
def set_instr(flavour, mpi, tmp_path):
    """Setup the instrument with common settings"""
    myinstrument = def_instrument(flavour)
    # myinstrument.set_seed(654321)
    myinstrument.set_instrument_base_dir(str(tmp_path))
    myinstrument.settings(run_path=(str(tmp_path)), mpi=mpi, seed=654321)
    # os.chdir(str(tmp_path))
    return myinstrument


@pytest.fixture
def config0(set_instr):
    """Base configuration for tests: direct beam with beam stop"""
    myinstrument = set_instr

    myinstrument.calculators["OriginCalc"].remove_component("collimator")

    myinstrument.master["energy"] = 110 * ureg.meV
    myinstrument.set_sample_by_name("None")
    myinstrument.sample_holder(None, None)
    return myinstrument


@pytest.fixture
def config_vanadium(config0):
    """Test with vanadium sample"""
    myinstrument = config0
    myinstrument.set_sample_by_name("vanadium")
    myinstrument.sample_shape("cylinder", r=0.01, h=0.04, th=0.001)
    myinstrument.sample.set_SPLIT(
        1
    )  # disabling the SPLIT due to the focusing on detector
    myinstrument.master["chopper_ratio"] = 3
    myinstrument.master["Efoc"] = 100
    return myinstrument


@pytest.fixture
def data_vanadium():
    return get_NXdetector_data(
        os.path.abspath(
            "institutes/ILL/instruments/Panther/HEAD/mcstas/data/032156.nxs"
        )
    )[0]


# @pytest.fixture
# def data_empty_holder():
#     return get_detector_data(
#         "institutes/ILL/instruments/D11/HEAD/mcstas/data/005721.nxs"
#     )


# @pytest.fixture
# def data_sample1():
#     return get_detector_data(
#         "institutes/ILL/instruments/D11/HEAD/mcstas/data/005711.nxs"
#     )


@pytest.fixture
def plot_settings(detector_names, tmp_path):
    file = os.path.abspath(
        "institutes/ILL/instruments/Panther/HEAD/mcstas/data/032156.nxs"
    )
    # Read NeXus
    detectors, time = get_NXdetectors(file)

    settings = {}
    for d in detector_names:
        settings[d] = {
            "cmap": "hsv",
            "origin": "lower",
            "norm": "log",
            "extent": [
                0,
                detectors[d]["pixel_size_x"][0]
                * detectors[d]["data"][:, :, 0].shape[0],
                0,
                detectors[d]["pixel_size_y"][0]
                * detectors[d]["data"][:, :, 0].shape[1],
            ],
        }

    settings["sample_monitor"] = {
        "cmap": "hsv",
        "origin": "lower",
    }
    return settings, tmp_path


"""
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
"""
# ------------------------------ tests
def test_comparison_type(comparison_type):
    assert comparison_type in ["exact", "stat"]


def test_name(set_instr):
    myinstrument = set_instr
    assert myinstrument.name == "Panther"


def test_detector_size(data_vanadium):
    detectors = data_vanadium
    assert detectors["detector"].shape == (288, 256, 512)


def test_plot_settings(plot_settings):
    print(plot_settings)
    assert True == True


def test_master_parameters(config0):
    myinstr = config0

    # print(help(myinstr.master))
    master_parameters = ["chopper_rpm", "chopper_ratio", "energy", "Efoc"]

    for par in master_parameters:
        assert par in myinstr.master


class PSD_TOF:
    I = None
    E = None
    N = None

    def __init__(self, detectors):
        i = []
        e = []
        n = []
        for d in detectors:
            if d.name == "detector":
                i.append(np.swapaxes(d.Intensity, 0, 1))
                e.append(np.swapaxes(d.Error, 0, 1))
                n.append(np.swapaxes(d.Ncount, 0, 1))
        self.I = np.stack(i)
        self.E = np.stack(e)
        self.N = np.stack(n)

        return


def test_direct_beam(config0, plot_settings, data_vanadium, detector_names):
    myinstrument = config0
    myinstrument.settings(ncount=1e5, mpi=1)
    myinstrument.set_sample_by_name("monitor")
    myinstrument.sample.set_parameters(
        xwidth=0.02,
        yheight=0.02,
        zdepth=0.01,
    )
    mycalc = myinstrument.calculators["OriginCalc"]
    mycalc.get_component("Monochromator").verbose = 1

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    sample_monitor = [d for d in detectors if d.name == "sample_monitor"]
    # print(sample_monitor)

    TOF = PSD_TOF(detectors)

    plot_sample_monitor(sample_monitor[0], plot_settings)

    # each tube is saved in a different file, need to compose them back
    plot_data_sim(data_vanadium, {"detector": TOF}, plot_settings)

    assert np.sum(sample_monitor[0].Ncount) == 8373
    assert np.sum(sample_monitor[0].Intensity) == pytest.approx(4152070)

    assert TOF.I.shape == (288, 256, 512)
    print(data_vanadium["detector"].shape)

    assert True == True


def test_vanadium(config_vanadium, plot_settings, data_vanadium, detector_names):
    myinstrument = config_vanadium
    myinstrument.sim_neutrons(1e5)

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]

    TOF = PSD_TOF(detectors)

    assert TOF.I.shape == (288, 256, 512)
    print(data_vanadium["detector"].shape)

    # each tube is saved in a different file, need to compose them back
    plot_data_sim(data_vanadium, {"detector": TOF}, plot_settings)

    assert True == True


def test_diagram(config0):
    myinstrument = config0
    myinstrument.sim_neutrons(1e5)

    mycalc = myinstrument.calculators["OriginCalc"]
    # myinstrument.run()
    mycalc.show_diagram(analysis=True, variable=None)
    assert True == True


def test_optim_HCS_focus(config0, tmp_path):
    import mcstasscript as ms

    myinstrument = config0
    myinstrument.sim_neutrons(1e8)
    myinstrument.set_sample_by_name("monitor")
    myinstrument.sample.set_parameters(
        xwidth=0.04,
        yheight=0.04,
        zdepth=0.01,
    )

    mycalc = myinstrument.calculators["OriginCalc"]
    mycalc.remove_component("detector")

    def create_diag(mycalc):
        diag = ms.Diagnostics(mycalc)
        diag.settings(ncount=1e6, suppress_output=True, mpi=1)
        # diag.show_settings()
        diag.clear_points()

        diag.clear_views()
        diag.add_view("x", bins=100)  # , limits=[-0.04,0.04])
        diag.add_view("y", bins=100)  # , same_scale=False, limits=[-0.07,0.07])
        # diag.add_view("e", same_scale=False)
        # diag.add_view("t", same_scale=False)
        # diag.add_view("t", "x", same_scale=False, log=True)
        # diag.add_view("t", "y", same_scale=False)
        diag.add_view("x", "y", same_scale=False, log=True)
        # diag.add_view("vz", same_scale=True)

        for component in mycalc.component_list:
            if component.name in ["BS", "slit_fermi"]:
                diag.add_point(before=component.name)

        return diag

    HCS = mycalc.get_component("HCS")
    results = []
    X = np.arange(0.05, 0.15, 0.01)
    # X = [HCS.focus_xw]
    Y = np.arange(0.02, 0.06, 0.01)

    X = np.arange(0.01, 0.08, 0.01)
    Y = np.arange(0.01, 0.3, 0.05)
    for x in X:
        y = 0
        HCS.focus_xw = x
        diag = create_diag(mycalc)
        diag.run()
        detectors = diag.instr.output.get_data()["data"]
        sample_monitor = [d for d in detectors if d.name == "sample_monitor"]

        results.append(
            [
                x,
                y,
                np.sum(sample_monitor[0].Ncount),
                np.sum(sample_monitor[0].Intensity),
                stat(sample_monitor[0].Ncount),
                stat(sample_monitor[0].Intensity),
            ]
        )
        # print(diag)
        fig = diag.plot()
        fig.savefig(os.path.join(tmp_path, "optim-{}_{}.pdf".format(x, y)))

    for y in Y:
        x = 0
        HCS.focus_yh = y
        diag = create_diag(mycalc)
        diag.run()
        detectors = diag.instr.output.get_data()["data"]
        sample_monitor = [d for d in detectors if d.name == "sample_monitor"]

        results.append(
            [
                x,
                y,
                np.sum(sample_monitor[0].Ncount),
                np.sum(sample_monitor[0].Intensity),
                stat(sample_monitor[0].Ncount),
                stat(sample_monitor[0].Intensity),
            ]
        )
        # print(diag)
        fig = diag.plot()
        fig.savefig(os.path.join(tmp_path, "optim-{}_{}.pdf".format(x, y)))

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


def test_diagnostics(config_vanadium, tmp_path):
    myinstrument = config_vanadium
    myinstrument.sim_neutrons(1e5)

    mycalc = myinstrument.calculators["OriginCalc"]
    mycalc.remove_component("detector")
    mycalc.get_component("HCS").focus_yh = 0.04
    import mcstasscript as ms

    fermi = mycalc.get_component("fermi_chopper")
    # fermi.delay = 0.0038
    fermi.phase = -10
    fermi.zero_time = 1
    print(fermi)

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
    diag.add_view("e", same_scale=False)
    # diag.add_view("l",same_scale=False)
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
