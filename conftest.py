import pytest
import pint

ureg = pint.get_application_registry()


def pytest_addoption(parser):
    parser.addoption(
        "--comparison_type",
        action="store",
        default="exact",
        help="comparison: exact or stat",
    )


@pytest.fixture
def comparison_type(request):
    """
    sets the type of comparison performed:
    - exact: checks the values to be identical (used to verify changes in the code)
    - stat:  used to verify if quantities are statistically compatible
    """
    return request.config.getoption("--comparison_type")


@pytest.fixture
def mpi():
    """Setting default mpi"""
    return 2


def test_comparison_type(comparison_type):
    assert comparison_type in ["exact", "stat"]
