from mcstasscript.interface.instr import McCode_instr
from mcstasscript.helper.mcstas_objects import Component


def HCS_source(mcstas_instrument: McCode_instr) -> Component:
    """
    Simple Gaussian source with mean energy and standard deviation to be set according to
    energy selection parameters according to the instrument.

    :param mcstas_instrument: previously declared McStas instrument to which add the source component

    The source should be added after the origin component is declared
    """
    Ei = mcstas_instrument.add_parameter(
        "double", "Ei", comment="Initial neutron energy", value=0, unit="meV"
    )
    dE = mcstas_instrument.add_parameter(
        "double", "dE", comment="Initial neutron energy dispertion", value=0, unit="meV"
    )

    slambda = mcstas_instrument.add_parameter(
        "double", "lambda", value=0, unit="angstrom"
    )
    sdlambda = mcstas_instrument.add_parameter(
        "double", "dlambda", value=0, unit="angstrom"
    )

    HCS = mcstas_instrument.add_component("HCS", "Source_simple")
    HCS.radius = 0.21 / 2
    #    HCS.dist = 2.155
    #    HCS.focus_xw = 0.170
    #    HCS.focus_yh = 0.120
    # HCS.target_index = +1
    # HCS.xwidth = 0.06
    # HCS.yheight = 0.12
    HCS.E0 = "Ei"
    HCS.dE = "dE"
    HCS.flux = 1e13
    #    HCS.lambda0 = "lambda"
    #    HCS.dlambda = "dlambda"
    #    HCS.Lmin = 0
    #    HCS.Lmax = 0

    HCS.gauss = 1
    HCS.set_AT(["0", " 0", " 0"], RELATIVE="Origin")

    # mcstas_instrument.add_parameter("lambda", comment="wavelength in angstroms")
    return HCS
