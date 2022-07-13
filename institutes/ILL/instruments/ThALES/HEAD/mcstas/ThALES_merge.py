"""
"""
from .ThALES import def_instrument as def_thales


def def_instrument():
    merge_instr = def_thales()
    for calcname in merge_instr.calculators:

        calc = merge_instr.calculators[calcname]
        calc.settings(checks=None)
        remove_comp = []
        for comp in calc.component_list[::-1]:
            if comp.component_name == "Progress_bar":
                break
            if comp.name == "detector_arm":
                comp.set_AT([0, 0, 0], RELATIVE="Origin")
                continue
            if comp.name == "detector_all":
                continue
            calc.remove_component(comp)

        calc.add_component(
            "Vin",
            "MCPL_input",
            AT=[0, 0, 0],
            RELATIVE="detector_arm",
            after="detector_arm",
        )
    #        print(calc.show_components())

    return merge_instr
