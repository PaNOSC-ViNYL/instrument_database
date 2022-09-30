# ------------------------------ For McStasscript instruments
import mcstasscript as ms
from mcstasscript.interface import instr

# ------------------------------ Mandatory classes to use
from libpyvinyl.Instrument import Instrument
from libpyvinyl.Parameters import Parameter

# ------------------------------ Extras
# import os  # to add the path of custom mcstas components


# list here all the common parts to be imported
from typing import List, Optional, Any


class McStasInstrumentBase(Instrument):
    """:class: BaseClass to be used for creating a good instrument"""

    def __init__(self, name, do_section=True):

        super().__init__(name, instrument_base_dir=".")

        self.__do_section = do_section
        self._temp_directory = "/dev/shm/mcstasscript/"

        self._single_calculator = False
        # this is specific for McStasscript instruments:
        # the components of the position for the sample and sample environment
        self._sample_environment_arm = None
        self._sample_arm = None

        self._calculator_with_sample = None

        self._sample_hash = None

    def add_sample_arms(self, mycalculator, previous_component):
        """The set_AT and set_ROTATE should be called afterwards.
        Setting at the same position as the previous component by default"""

        self._calculator_with_sample = mycalculator

        self._sample_environment_arm = mycalculator.add_component(
            "sample_environment_arm",
            "Arm",
            AT=[0, 0, 0],
            RELATIVE=previous_component,
        )

        # do we want a sample rotation parameter by default?
        self._sample_arm = mycalculator.add_component(
            "sample_arm",
            "Arm",
            AT=[0, 0, 0],
            ROTATED=[0, "sample_rotation", 0],
            RELATIVE=previous_component,
        )

        a3 = mycalculator.add_parameter(
            "double",
            "sample_rotation",
            comment="sample table rotation angle",
            unit="degree",
            value=0,
        )

    def add_new_section(
        self, new_calcname, output_arm=None, section_name=None, hasSample=False
    ):
        """
        Method to divide the instrument into sections
        each storing the neutrons in an MCPL file
        that is used as input of the following section.

        If output_arm is not provided or None, there is no bridge being created since there is not previous section

        Always returns the origin of the neutrons and sample and sample environments at that same position. They need to be displaced afterwards with a
        self._sample_environment_arm.set_AT([x,y,z], RELATIVE=vin)
        """

        # This is to obtain the same instrument as a single McStas instrument

        calculatorname = list(self.calculators.keys())[-1]
        mycalculator = self.calculators[calculatorname]
        vin = output_arm

        if output_arm is not None and self.__do_section:
            oldcalculator = mycalculator
            output = oldcalculator.add_component(
                section_name, "MCPL_output", AT=[0, 0, 0], RELATIVE=output_arm
            )
            output.filename = '"' + section_name + '"'  #'"sSAMPLE"'

        # ------------------------------------------------------------
        if output_arm is None or self.__do_section:
            mycalculator = instr.McStas_instr(
                new_calcname, input_path=self._temp_directory
            )
            self.add_calculator(mycalculator)

            Origin = mycalculator.add_component("Origin", "Progress_bar")
            Origin.set_AT(["0", "0", "0"], RELATIVE="ABSOLUTE")

            vin = Origin
            if output_arm is not None:
                mycalculator.add_parameter(
                    "string",
                    "vin_filename",
                    unit="",
                    comment="",
                    value='"none"',  # + oldcalculator.output["mcpl"] + '"',
                )
                vin = mycalculator.add_component(
                    "Vin",
                    "MCPL_input",
                    AT=[0, 0, 0],
                    after="Origin",
                )
                vin.filename = "vin_filename"

        if hasSample:
            add_sample_arms(mycalculator, vin)

        return (mycalculator, vin)

    def calculator_before_run(self, calcname) -> None:
        """
        Operations to performe before calling the backengine
        """
        calc = self.calculators[calcname]
        previous_calc = None
        keys = list(self.calculators.keys())
        index = keys.index(calcname)
        if index > 0:
            previous_calc = self.calculators[keys[index - 1]]

        vin = None
        try:
            vin = calc.parameters["vin_filename"]
        except KeyError:
            return

        if vin is not None and previous_calc is not None:
            if "mcpl" not in previous_calc.output_keys:
                raise RuntimeError(
                    f"Something is wrong with the instrument implementation:\ncalculator {calcname} has an MCPL_input component but the previous does not have an output"
                )
            # update the parameter value
            # this is needed because the output is generated after calling the backengine
            # so the filename is not set before
            calc.parameters["vin_filename"].value = (
                '"' + previous_calc.output["mcpl"].filename + '"'
            )

    # ------------------------------ utility methods made available for the users
    def sim_neutrons(self, number) -> None:
        """Method to set the number of neutrons to be simulated"""
        for calc in self.calculators:
            mycalc = self.calculators[calc]
            mycalc.settings(ncount=number)
