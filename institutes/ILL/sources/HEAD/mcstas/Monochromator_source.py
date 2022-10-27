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

    #    slambda = mcstas_instrument.add_parameter(
    #        "double", "lambda", value=0, unit="angstrom"
    #    )
    sdlambda = mcstas_instrument.add_parameter(
        "double", "dlambda", value=0, unit="angstrom"
    )

    last_component = mcstas_instrument.get_last_component()
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
    HCS.set_AT(["0", " 0", " 0"], RELATIVE=last_component)

    HCS.target_index = 1

    # override the value of lambda by the value of the angle
    mcstas_instrument.add_declare_var("double", "lambda")

    mcstas_instrument.append_initialize(
        "lambda =2*sin(0.5*a2*DEG2RAD)*monochromator_d;"
    )
    mcstas_instrument.append_initialize('printf("lambda: %.2f\\n", lambda);')
    mcstas_instrument.append_initialize("dlambda=0.01*lambda;")
    # the following python string is an alternative way of doing:
    mcstas_instrument.add_declare_var("double", "energy")
    mcstas_instrument.append_initialize("Ei = 81.80421/(lambda*lambda);")
    mcstas_instrument.add_declare_var("double", "denergy")
    mcstas_instrument.append_initialize("dE = 2*Ei*dlambda/lambda;")
    mcstas_instrument.append_initialize(
        'printf("lambda: %.2f +/- %.2f\\n", lambda, dlambda);'
    )
    mcstas_instrument.append_initialize('printf("energy: %.2f +/- %.2f\\n", Ei, dE);')

    # mcstas_instrument.add_parameter("lambda", comment="wavelength in angstroms")
    return HCS
