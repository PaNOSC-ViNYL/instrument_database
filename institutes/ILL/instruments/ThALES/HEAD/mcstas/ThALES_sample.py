"""
"""
from .ThALES import def_instrument as def_thales


def def_instrument():
    instr = def_thales()
    instr.sample = "empty"
    for calcname in instr.calculators:

        calc = instr.calculators[calcname]
        calc.name = calcname + "_sample"
        calc.settings(checks=False)
        remove_comp = []
        first_component_to_remove = "sample_mcpl_arm"
        first_component_to_remove_found = False
        for comp in calc.component_list[::-1]:
            if comp.name == first_component_to_remove:
                first_component_to_remove_found = True
                comp.set_AT([0, 0, 0], RELATIVE="ABSOLUTE")
                continue
            if first_component_to_remove_found is False:
                continue
            if comp.component_name == "Progress_bar":
                break
            calc.remove_component(comp)

        comp = calc.component_list[2]
        if comp.component_name != "MCPL_output":
            print(calc.show_components())
            raise RuntimeError("MCPL output component is not immediately after the Arm")
        calc.remove_component(comp)

        vin = calc.add_component(
            "Vin",
            "MCPL_input",
            AT=[0, 0, 0],
            #            RELATIVE="detector_arm",
            after="Origin",
        )
        vin.filename = "filelist"

    #

    return instr
