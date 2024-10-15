# pytest  -k holder --comparison_type=stat
# pytest  -k holder --comparison_type=exact
#
# - [ ] flux @ sample (2cm x 2cm) = 5e5 n/cm/s
# - [ ] vanadium sample: r=5mm h=40mm th=1mm
# - [ ] plot of 3D data for comparison
# - [ ] profile of 3D data for comparison
from pytest_check import check
from conftest import ureg, PlotHelper, ComparePlots, PdfPages, plt
import sys, os, pytest, math
import h5py
import numpy as np

from Panther import def_instrument
from conftest import plot_sample_monitor


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


def plot_data_sim(data, sim, time, plot_settings):
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
        "detector": "c",
    }

    fig = plt.figure(layout="constrained", dpi=92, figsize=(15, 15))
    axs = fig.subplot_mosaic([["cx"], ["cy"]])
    return
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

    myinstrument.master["energy"] = 10 * ureg.meV
    myinstrument.master["mono_index"] = '"pg002"'
    myinstrument.master["a2"] = 50.46 * ureg.degree
    myinstrument.master["chopper_rpm"] = 7998

    myinstrument.master["a2"] = 0 * ureg.degree
    myinstrument.master["mono_index"] = '"pg"'
    myinstrument.master["chopper_rpm"] = 0

    # myinstrument.calculators["OriginCalc"].parameters["mono_index"] = "0"
    myinstrument.set_sample_by_name("None")
    myinstrument.sample_holder(None, None)
    return myinstrument


@pytest.fixture
def config_before_fermi(config0):
    """
    Instrument configuration config0 with MCPL output before the fermi
    """
    myinstrument = config0
    myinstrument.master["a2"] = 0 * ureg.degree
    myinstrument.master["mono_index"] = "0"
    myinstrument.master["chopper_rpm"] = 0

    myinstrument.do_section(True)

    mycalc = myinstrument.calculators["OriginCalc"]
    for component in reversed(mycalc.component_list):
        if component.name == "monitor":
            break
        mycalc.remove_component(component)
    oldfermi = mycalc.get_component("fermi_chopper")
    comp = mycalc.get_last_component()

    mycalculator, mono_out = myinstrument.add_new_section(
        "FermiChopperCalc", output_arm="Monochromator_Out"
    )

    fermi = mycalculator.add_component(oldfermi.name, oldfermi.component_name)
    fermi.set_parameters(oldfermi.__dict__)
    fermi.set_RELATIVE(mono_out)

    monitor = mycalculator.add_component(comp.name, comp.component_name)
    monitor.set_parameters(comp.__dict__)
    # print(fermi)

    chopper_rpm = mycalculator.add_parameter(
        "double",
        "chopper_rpm",
        comment="Fermi chopper speed",
        unit="",
        value=0,
    )
    myinstrument.add_parameter_to_master(chopper_rpm.name, mycalculator, chopper_rpm)

    ei = mycalculator.add_parameter("double", "Ei", unit="meV")
    myinstrument.add_parameter_to_master("energy", mycalculator, ei)
    mycalculator.add_declare_var("double", "lambda")
    mycalculator.append_initialize("lambda = sqrt(81.80421036/Ei);")
    mycalculator.add_declare_var("double", "neutron_velocity")
    mycalculator.append_initialize("neutron_velocity = 3956.034012/lambda;")
    mycalculator.append_initialize('printf("nv = %2f\\n", neutron_velocity);')
    mycalculator.append_initialize('printf("lambda = %.2f\\n", lambda);')

    # print(fermi)
    mycalc.remove_component(oldfermi)
    mycalc.remove_component(comp)
    # print(myinstrument)
    # print(mycalculator)
    return myinstrument


@pytest.fixture
def config_sample_monitor(config0):
    """
    Instrument configuration config0 with sample monitor of 2x2 cm2 at sample position
    """
    myinstrument = config0
    myinstrument.master["energy"] = 110 * ureg.meV
    myinstrument.set_sample_by_name("monitor")
    myinstrument.sample.set_parameters(
        xwidth=0.04,
        yheight=0.08,
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
def config_vanadium(config0):
    """Test with vanadium sample"""
    myinstrument = config0
    myinstrument.set_sample_by_name("vanadium")
    myinstrument.sample_shape("cylinder", r=0.01, h=0.04, th=0.001)
    myinstrument.sample.set_SPLIT(
        1
    )  # disabling the SPLIT due to the focusing on detector
    myinstrument.master["chopper_ratio"] = 1
    # myinstrument.master["Efoc"] = 100
    return myinstrument


@pytest.fixture
def data_vanadium():
    return get_NXdetector_data(
        os.path.abspath(
            "institutes/ILL/instruments/Panther/HEAD/mcstas/data/032156.nxs"
        )
    )


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
            "px": detectors[d]["pixel_size_x"][0],
            "py": detectors[d]["pixel_size_y"][0],
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
    detectors = data_vanadium[0]
    assert detectors["detector"].shape == (288, 256, 512)


def test_plot_settings(plot_settings):
    print(plot_settings)
    assert True == True


def test_master_parameters(config0):
    myinstr = config0

    # print(help(myinstr.master))
    master_parameters = [
        "mono_index",
        "chopper_rpm",
        "chopper_ratio",
        "energy",
        "Efoc",
        "tdelay",
        "twidth",
    ]

    # check that master parameters here defined are present in the simulation
    for par in master_parameters:
        assert par in myinstr.master

    # check that the test accounts for all the master parameters defined in the simulation
    for par in myinstr.master:
        assert par.name in master_parameters


def test_mono_index(config0):
    myinstrument = config0
    myinstrument.master["chopper_rpm"] = 0
    myinstrument.settings(ncount=10, mpi=1)
    *_, calc = myinstrument.calculators

    while calc != "OriginCalc":
        myinstrument.remove_calculator(calc)

    mycalc = myinstrument.calculators["OriginCalc"]
    for component in reversed(mycalc.component_list):
        if component.name == "HCS":
            break
        mycalc.remove_component(component)

    myinstrument.run()


def test_intensity_at_sample(config_sample_monitor, plot_settings):
    myinstrument = config_sample_monitor
    myinstrument.settings(ncount=2e8)

    try:
        myinstrument.run()
    except:
        print("ciao")
        print(myinstrument.output)
    detectors = myinstrument.output.get_data()["data"]
    print(detectors)
    detector = [d for d in detectors if d.name in ["sample_monitor"]][0]

    # check simulation vs previous version
    assert np.sum(detector.Ncount) == 32575
    assert np.sum(detector.Intensity) == pytest.approx(9580, abs=1)
    ph = PlotHelper(detector.Intensity.T)
    ph.set_limits(detector.metadata.limits)
    ph.stats[0].mean = pytest.approx(49.68, abs=0.02)
    ph.stats[1].mean = pytest.approx(49.58, abs=0.02)
    ph.stats[0].std = pytest.approx(22.52, abs=0.02)
    ph.stats[1].std = pytest.approx(22.40, abs=0.02)

    detector = {d.name: d for d in detectors if d.name in ["sample_monitor"]}
    plot_sample_monitor(detector, plot_settings)


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
        self.Intensity = np.stack(i)
        self.E = np.stack(e)
        self.N = np.stack(n)

        return


def test_direct_beam(config0, plot_settings, data_vanadium, detector_names):
    data, time = data_vanadium
    myinstrument = config0
    myinstrument.settings(ncount=1e5, mpi=1)
    mycalc = myinstrument.calculators["OriginCalc"]
    mycalc.get_component("Monochromator").verbose = 1

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]
    sample_monitor = [d for d in detectors if d.name == "sample_monitor"]
    # print(sample_monitor)

    TOF = PSD_TOF(detectors)

    # each tube is saved in a different file, need to compose them back
    plot_data_sim(data, {"detector": TOF}, time, plot_settings)

    assert np.sum(sample_monitor[0].Ncount) == 8373
    assert np.sum(sample_monitor[0].Intensity) == pytest.approx(4152070)

    assert TOF.I.shape == (288, 256, 512)
    print(data_vanadium["detector"].shape)

    assert True == True


def test_vanadium(config_vanadium, plot_settings, data_vanadium, detector_names):
    myinstrument = config_vanadium
    myinstrument.settings(ncount=1e5)

    myinstrument.run()

    detectors = myinstrument.output.get_data()["data"]

    TOF = PSD_TOF(detectors)

    assert TOF.I.shape == (288, 256, 512)
    print(data_vanadium["detector"].shape)

    # each tube is saved in a different file, need to compose them back
    plot_data_sim(data_vanadium, {"detector": TOF}, plot_settings)

    assert True == True


def test_diagram(config0, tmp_path):
    myinstrument = config0
    myinstrument.settings(ncount=1e5)

    mycalc = myinstrument.calculators["OriginCalc"]
    # myinstrument.run()
    mycalc.show_diagram(analysis=True, variable=None)
    plt.savefig(os.path.join(tmp_path, "diagram.pdf"))
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


def test_BCs(config0, tmp_path):
    myinstrument = config0
    myinstrument.settings(ncount=1e5)

    # remove all simulation components after the OriginCalc calculator
    *_, calc = myinstrument.calculators

    while calc != "OriginCalc":
        myinstrument.remove_calculator(calc)

    # remove all components after BC5
    import mcstasscript as ms

    mycalc = myinstrument.calculators["OriginCalc"]

    for component in reversed(mycalc.component_list):
        if component.name == "BC5":
            break
        mycalc.remove_component(component)

    monitor = mycalc.add_component(
        "monitor", "Monitor_nD", AT=2.001, RELATIVE="diaphragm1"
    )
    monitor.set_parameters(
        xwidth=0.150,
        yheight=0.500,
        zdepth=0.01,
        filename='"monitor.dat"',
        restore_neutron=1,
        options='"box intensity, bins=1"',
        # options='"x bins={} y bins={} file={}"'.format(1, 1, "counter.dat"
    )

    def create_diag(mycalc):
        diag = ms.Diagnostics(mycalc)
        diag.settings(ncount=1e5, suppress_output=False, mpi=1)
        diag.show_settings()
        diag.clear_points()

        for component in diag.instr.component_list:
            if component.component_name in [
                "DiskChopper",
                "Monochromator_curved",
                "FermiChopper",
            ]:
                diag.add_point(before=component.name)
                diag.add_point(after=component.name)
            if component.name in ["monitor"]:
                diag.add_point(before=component.name)

        return diag

    counts = []
    intensity = []
    phases = [30]
    distances = myinstrument.distances()
    print(distances)
    for energy in [10, 50, 100]:
        myinstrument.master["energy"] = energy * ureg.meV

        for phase in phases:
            for ibc in [2, 3, 4, 5]:
                bc = mycalc.get_component("BC{}".format(ibc))
                # bc.phase = "{}+{}".format(bc.phase, phase)
                # del bc.phase
                d = myinstrument.distances()["L1{}".format(ibc)]
                bc.delay = (
                    "{dist}/neutron_velocity + {phase_init}/({omega})/360".format(
                        dist=d,
                        phase_init=phase,
                        omega=bc.nu,
                    )
                )

                print(bc)
            mycalc.append_initialize(
                'printf("BC2 delay = %.6e\\n", {phase_init}/({omega})/360);\n'.format(
                    phase_init=phase,
                    omega=bc.nu,
                )
            )
            diag = create_diag(mycalc)
            diag.run()

            diag.clear_views()
            diag.add_view("e", same_scale=False)
            diag.add_view(
                "t",
                same_scale=True,
                left_min=-0.003,
            )
            diag.add_view("t", "x", same_scale=False, log=True)
            diag.add_view("t", "y", same_scale=False)
            diag.add_view("x", "y", same_scale=False, log=True)

            print(diag)
            fig = diag.plot()
            fig.savefig(
                os.path.join(
                    tmp_path, "diagnostics_bc1_e{}_ph{}.pdf".format(phase, energy)
                )
            )

            detectors = diag.instr.output.get_data()["data"]
            sample_monitor = [d for d in detectors if d.name == "monitor"]
            print(sample_monitor)
            counts.append(sample_monitor[0].Ncount)
            intensity.append(sample_monitor[0].Intensity)
            # assert np.sum(sample_monitor[0].Ncount) == 8373
            # assert np.sum(sample_monitor[0].Intensity) == pytest.approx(4152070)

    print(phases)
    print(counts)
    print(intensity)

    assert True == True


def test_fermi_phase(config_before_fermi, tmp_path):
    myinstrument = config_before_fermi
    myinstrument.settings(ncount=1e8, mpi=8)
    print(myinstrument.master)

    mycalc = myinstrument.calculators["FermiChopperCalc"]
    # print(mycalc.parameters)
    # print(mycalc.component_list)

    import mcstasscript as ms

    def create_diag(mycalc):
        diag = ms.Diagnostics(mycalc)
        diag.settings(ncount=1e7, suppress_output=False, mpi=8)
        # diag.show_settings()
        diag.clear_points()

        for component in diag.instr.component_list:
            if component.component_name in [
                "Monochromator_curved",
                "FermiChopper",
            ]:
                diag.add_point(before=component.name)
                diag.add_point(after=component.name)
            if component.name in ["monitor"]:
                diag.add_point(before=component.name)

        diag.clear_views()
        diag.add_view("e", same_scale=False)
        diag.add_view(
            "t",
            same_scale=True,
            # left_lim=-0.003,
            right_lim=0.01,
        )
        diag.add_view("t", "x", same_scale=False, log=True)
        diag.add_view("t", "y", same_scale=False)
        diag.add_view("x", "y", same_scale=False, log=True)
        diag.add_view("x", "z", same_scale=False, log=True)

        return diag

    counts = []
    intensity = []
    phases = range(50, 100, 10)
    for energy in [10, 50, 100]:
        myinstrument.master["energy"] = energy * ureg.meV
        myinstrument.run()

        for phase in phases:
            bc = mycalc.get_component("fermi_chopper")
            d = myinstrument.distances()["L1c"]
            bc.delay = "{dist}/neutron_velocity + {phase_init}/({omega})/ 360".format(
                dist=myinstrument.distances()["L1c"],
                phase_init=phase,
                omega="chopper_rpm/60",
            )
            # bc.verbose = 2
            # bc.phase = -90
            # print(bc)

            diag = create_diag(mycalc)
            diag.run()

            # print(diag)
            fig = diag.plot()
            fig.savefig(
                os.path.join(
                    tmp_path, "diagnostics_fermi_e{}_ph{}.pdf".format(energy, phase)
                )
            )

            detectors = diag.instr.output.get_data()["data"]
            sample_monitor = [d for d in detectors if d.name == "monitor"]
            print(sample_monitor)
            counts.append(sample_monitor[0].Ncount)
            intensity.append(sample_monitor[0].Intensity)
            # assert np.sum(sample_monitor[0].Ncount) == 8373
            # assert np.sum(sample_monitor[0].Intensity) == pytest.approx(4152070)

    print(phases)
    print(counts)
    print(intensity)

    assert True == True


def test_BC1(config0, tmp_path):
    """used to test the phase_init for the first BC"""
    myinstrument = config0
    myinstrument.settings(ncount=1e5)

    # remove all simulation components after the OriginCalc calculator
    *_, calc = myinstrument.calculators

    while calc != "OriginCalc":
        myinstrument.remove_calculator(calc)

    # remove all components after BC5
    import mcstasscript as ms

    mycalc = myinstrument.calculators["OriginCalc"]

    for component in reversed(mycalc.component_list):
        if component.name == "diaphragm1":
            break
        mycalc.remove_component(component)

    monitor = mycalc.add_component(
        "monitor", "Monitor_nD", AT=2.001, RELATIVE="diaphragm1"
    )
    monitor.set_parameters(
        xwidth=0.150,
        yheight=0.500,
        zdepth=0.01,
        filename='"monitor.dat"',
        restore_neutron=1,
        options='"box intensity, bins=1"',
        # options='"x bins={} y bins={} file={}"'.format(1, 1, "counter.dat"
    )

    def create_diag(mycalc):
        diag = ms.Diagnostics(mycalc)
        diag.settings(ncount=1e5, suppress_output=False, mpi=1)
        diag.show_settings()
        diag.clear_points()

        for component in diag.instr.component_list:
            if component.component_name in [
                "DiskChopper",
                "Slit",
                "Monochromator_curved",
                "FermiChopper",
            ]:
                diag.add_point(before=component.name)
                diag.add_point(after=component.name)
            if component.name in ["monitor"]:
                diag.add_point(before=component.name)

        return diag

    counts = []
    intensity = []
    phases = [20, 29, 30, 31]
    for energy in [10, 50, 100]:
        myinstrument.master["energy"] = energy * ureg.meV

        for phase in phases:
            bc = mycalc.get_component("BC1")
            bc.phase = phase
            print(bc)

            diag = create_diag(mycalc)
            diag.run()

            diag.clear_views()
            diag.add_view("e", same_scale=False)
            diag.add_view("t", same_scale=True, left_min=-0.003, right_lim=0.003)
            diag.add_view("t", "x", same_scale=False, log=True)
            diag.add_view("t", "y", same_scale=False)
            diag.add_view("x", "y", same_scale=False, log=True)

            print(diag)
            fig = diag.plot()
            fig.savefig(
                os.path.join(
                    tmp_path, "diagnostics_bc1_ph{}_e{}.pdf".format(phase, energy)
                )
            )

            detectors = diag.instr.output.get_data()["data"]
            sample_monitor = [d for d in detectors if d.name == "monitor"]
            print(sample_monitor)
            counts.append(sample_monitor[0].Ncount)
            intensity.append(sample_monitor[0].Intensity)
            # assert np.sum(sample_monitor[0].Ncount) == 8373
            # assert np.sum(sample_monitor[0].Intensity) == pytest.approx(4152070)

    print(phases)
    print(counts)
    print(intensity)

    assert True == True


def test_diagnostics(config_vanadium, tmp_path):
    myinstrument = config_vanadium
    myinstrument.settings(ncount=1e5)

    myinstrument.master["energy"] = 50 * ureg.meV

    mycalc = myinstrument.calculators["OriginCalc"]
    mycalc.remove_component("detector")
    # mycalc.get_component("HCS").focus_yh = 0.04
    import mcstasscript as ms

    phase = 30
    fermi = mycalc.get_component("fermi_chopper")
    fermi.verbose = 0
    d = myinstrument.distances()["L1c"]
    # fermi.delay = "{dist}/neutron_velocity + {phase_init}/({omega})/ 360".format(
    #    dist=myinstrument.distances()["L1c"],
    #    phase_init=phase,
    #    omega="chopper_rpm/60",
    # )
    #    fermi.phase = -90
    print(fermi)

    BC1 = mycalc.get_component("BC1")
    # BC1.nu = "chopper_rpm/chopper_ratio/60*2"
    for iBC in range(2, 5):
        nBC = "BC{}".format(iBC)
        dist = myinstrument.distances()["L1{}".format(iBC)]
        BC = mycalc.get_component(nBC)
        # BC.nu = BC1.nu
        # BC.delay = "{dist}/neutron_velocity + {phase_init}/({omega})/360".format(
        #    dist=dist, phase_init=24, omega=BC1.nu
        # )
        print(BC)

    diag = ms.Diagnostics(mycalc)
    diag.settings(ncount=1e6, suppress_output=False, mpi=1)
    diag.show_settings()
    diag.clear_points()

    for component in diag.instr.component_list:
        if component.component_name in ["Progress_bar", "Arm", "Monochromator_curved"]:
            continue
        if component.component_name in ["Beamstop"]:
            diag.add_point(before=component.name)
            continue
        if component.name in ["VCS", "HCS"]:
            diag.add_point(after=component.name)
            continue
        if component.component_name in [
            "DiskChopper",
            "Slit",
            "FermiChopper",
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
    diag.add_view("e", same_scale=True)
    # diag.add_view("l",same_scale=False)
    diag.add_view("t", same_scale=True)
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
