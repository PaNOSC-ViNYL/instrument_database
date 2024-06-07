import pytest
import pint

ureg = pint.get_application_registry()

from mcstas.helper_functions import PlotHelper, ComparePlots
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def pytest_addoption(parser):
    parser.addoption(
        "--comparison_type",
        action="store",
        default="exact",
        help="comparison: exact or stat",
    )


@pytest.fixture(autouse=False)
def comparison_type(request):
    """
    sets the type of comparison performed:
    - exact: checks the values to be identical (used to verify changes in the code)
    - stat:  used to verify if quantities are statistically compatible
    """
    comparison_type = request.config.getoption("--comparison_type")
    return comparison_type


@pytest.fixture
def mpi():
    """Setting default mpi"""
    return 2


def test_comparison_type(request):
    assert comparison_type in ["exact", "stat"]


def plot_sample_monitor(sim, plot_settings):
    pixel_sizes = plot_settings[0]
    tmp_path = plot_settings[1]
    time = 1

    ph1 = {}
    for d in sim:
        print("Limits: ", sim[d].metadata.limits)
        ph1[d] = PlotHelper(sim[d].Intensity.T * time)
        ph1[d].set_limits(sim[d].metadata.limits)
        # ph1[d].set_pixel_size(0.0002, 0.0002)

    compare_plots = ComparePlots(ph1, None, "Simulation", "Data")

    plot_positions = {
        "sample_monitor": "c",
    }

    fig1D = plt.figure(layout="constrained", dpi=92, figsize=(10, 15))
    axs = fig1D.subplot_mosaic([["cx"], ["cy"]])
    compare_plots.plot_marginals(plot_positions, axs)

    fig2D = plt.figure(layout="constrained", dpi=92, figsize=(10, 15))
    axs = fig2D.subplot_mosaic([["c_sim"]])
    compare_plots.plot_2Ds(fig2D, plot_positions, axs, "c")

    with PdfPages(f"{tmp_path}/plot.pdf") as pdf:
        pdf.savefig(fig1D)
        pdf.savefig(fig2D)

    print(f"Check file: {tmp_path}/plot.pdf")
