from mcstasscript.interface.instr import McCode_instr
from mcstasscript.helper.mcstas_objects import Component


def VCS_source(mcstas_instrument: McCode_instr) -> Component:
    """
    Description of the VCS source of the ILL

    :param mcstas_instrument: previously declared McStas instrument to which add the source component

    The source should be added after the origin component is declared
    """
    slambda = mcstas_instrument.add_parameter(
        "double", "lambda", value=0, unit="angstrom"
    )
    sdlambda = mcstas_instrument.add_parameter(
        "double", "dlambda", value=0, unit="angstrom"
    )

    VCS = mcstas_instrument.add_component(
        "VCS", "Source_gen", AT=[0, 0, 0], RELATIVE="Origin"
    )

    VCS.set_parameters(
        yheight=0.22,
        xwidth=0.14,
        focus_xw=0.038,
        focus_yh=0.2,
        lambda0="lambda",
        dlambda="dlambda",
        T1=216.8,
        I1=1.24e13,  # /* VCS parameters */
        T2=33.9,
        I2=1.02e13,
        T3=16.7,
        I3=3.0423e12,
        verbose=1,
    )


def HCS_source(mcstas_instrument: McCode_instr) -> Component:
    """
    Description of the Hot source of the ILL

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

    HCS = mcstas_instrument.add_component("HCS", "Source_gen")
    HCS.radius = 0.21 / 2
    HCS.dist = 2.155
    HCS.focus_xw = 0.170
    HCS.focus_yh = 0.120
    # HCS.target_index = +1
    HCS.E0 = "Ei"
    HCS.dE = "dE"
    #    HCS.lambda0 = "lambda"
    #    HCS.dlambda = "dlambda"
    #    HCS.Lmin = 0
    #    HCS.Lmax = 0
    HCS.verbose = 1
    HCS.I1 = 2.78e13
    HCS.T1 = 40.1
    HCS.T2 = 145.8
    HCS.I2 = 3.44e13
    HCS.T3 = 413.5
    HCS.I3 = 10.22e12
    HCS.zdepth = 0.15
    HCS.set_AT(["0", " 0", " 0"], RELATIVE="Origin")

    # mcstas_instrument.add_parameter("lambda", comment="wavelength in angstroms")
    return HCS
