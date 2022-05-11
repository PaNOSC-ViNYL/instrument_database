from mcstasscript.interface.instr import McCode_instr
from mcstasscript.helper.mcstas_objects import Component


def HCS_source(mcstas_instrument: McCode_instr) -> Component:
    """
    Description of the Hot source of the ILL

    :param mcstas_instrument: previously declared McStas instrument to which add the source component

    The source should be added after the origin component is declared
    """
    HCS = mcstas_instrument.add_component("HCS", "Source_gen")
    HCS.radius = 0.21 / 2
    HCS.dist = 2.155
    HCS.focus_xw = 0.170
    HCS.focus_yh = 0.120
    HCS.E0 = "Ei"
    HCS.dE = "ThALES_dE"
    HCS.verbose = 0
    HCS.I1 = 2.78e13
    HCS.T1 = 40.1
    HCS.T2 = 145.8
    HCS.I2 = 3.44e13
    HCS.T3 = 413.5
    HCS.I3 = 10.22e12
    HCS.zdepth = 0.15
    HCS.set_AT(["0", " 0", " 0"], RELATIVE="Origin")

    mcstas_instrument.add_parameter("lambda", comment="wavelength in angstroms")
    return HCS


def add_source_parameters():
    QuickSource.lambda0 = "lambda"
    QuickSource.dlambda = "dlambda"
