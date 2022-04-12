from mcstasscript.interface.instr import McCode_instr
from mcstasscript.helper.mcstas_objects import Component


def D22_quick(instrument: McCode_instr) -> Component:
    QuickSource = instrument.add_component("QuickSource", "Source_simple")
    QuickSource.yheight = 0.04
    QuickSource.xwidth = 0.04
    QuickSource.focus_xw = "sample_size_r"
    QuickSource.focus_yh = "sample_size_y"
    QuickSource.lambda0 = "lambda"
    QuickSource.dlambda = "dlambda"
    QuickSource.flux = 1e13
    QuickSource.gauss = 0
    QuickSource.target_index = +5
    QuickSource.set_AT(["0", " 0", " 0"], RELATIVE="H512_Before_VS")
    return QuickSource
