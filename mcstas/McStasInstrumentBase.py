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
        self._temp_directory = "./"  # "/dev/shm/mcstasscript/"

        self._single_calculator = False
        # this is specific for McStasscript instruments:
        # the components of the position for the sample and sample environment
        self._sample_environment_arm = None
        self._sample_arm = None

        self._calculator_with_sample = None

        self._sample_hash = None

        self._add_monitors = True

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
        self, new_calcname, output_arm=None, hasSample=False, section_name=None
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
        if section_name is None:
            section_name = new_calcname

        if output_arm is not None and self.__do_section:
            calculatorname = list(self.calculators.keys())[-1]
            oldcalculator = self.calculators[calculatorname]
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
                # this parameter is just to have some flexibility on the compiled instrument
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

                mycalculator.input = oldcalculator.output
            else:
                vin = output_arm
                calculatorname = list(self.calculators.keys())[-1]
                mycalculator = self.calculators[calculatorname]

        if hasSample:
            self.add_sample_arms(mycalculator, vin)

        return (mycalculator, vin)

    def add_monitor(
        self,
        calculator,
        name: str,
        position=[0, 0, 0],
        width=0.10,
        height=0.15,
        bins=100,
        radius=0.30,
        energy=5,
        denergy=4.9,
    ):
        if self._add_monitors is False:
            return

        instr = calculator
        monitor = instr.add_component(name + "_DEBUG", "Monitor_nD")
        monitor.filename = '"' + name + '_DEBUG"'
        monitor.xwidth = width
        monitor.yheight = height
        monitor.bins = bins
        monitor.options = '"energy limits=[{0:.2f} {1:.2f}] x y, parallel"'.format(
            #        energy - denergy, energy + denergy
            0,
            25,
        )

        # monitor.options = '"energy limits=[Ei-dE, Ei+dE] x y, borders, parallel"'
        monitor.set_AT(position, RELATIVE="PREVIOUS")

        psd = instr.add_component(name + "_psd_DEBUG", "Monitor_nD")
        psd.filename = '"' + name + '_DEBUG"'
        psd.xwidth = width
        psd.yheight = height
        psd.bins = bins
        psd.options = '"x y, parallel"'
        psd.set_AT([0, 0, 0], RELATIVE="PREVIOUS")

        psd = instr.add_component(name + "_psdcyl_DEBUG", "Monitor_nD")
        psd.filename = '"' + name + '_DEBUG"'
        psd.xwidth = 2 * radius
        psd.yheight = height
        psd.options = '"theta bins=360, y bins=20, cylinder, parallel"'
        psd.set_AT([0, 0, 0], RELATIVE="PREVIOUS")

        return psd

    # ------------------------------
    # this implements what is foreseen in the libpyvinyl.Instrument class
    def set_sample_by_name(self, name: str) -> None:
        """Set the sample component
        Always put a sample relative to the _sample_arm and after the _sample_arm component
        In case of an Sqw sample, a parameter Sqw_file is added"""
        # print(f"Setting sample to: {name}")
        mycalculator = self._calculator_with_sample

        if self.sample is not None:
            mycalculator.remove_component(self.sample)
        if name in ["empty", "Empty", "None", "none"]:
            self.sample = None
        elif name in ["v_sample", "vanadium"]:
            self.sample_name = "vanadium"
            self.sample = mycalculator.add_component(
                self.sample_name,
                "V_sample",
                AT=[0, 0, 0],
                # ROTATED=[0, "a4", 0],
                RELATIVE=self._sample_arm,
                after=self._sample_arm,
            )
            sample = self.sample
            # sample_radius, sample_height and sample_thickness added to the
            sample.radius = 0.02  # "sample_radius"
            sample.yheight = 0.02  # "sample_height"
            sample.thickness = 0.001  # "sample_thickness"

            # calculator parameters
            sample.focus_xw = 0.04
            sample.focus_yh = 0.12
            sample.target_z = 0.25
            #    v_sample.append_EXTEND("if(flag==SCATTERED) ABSORB;")
            # Absorption fraction           =0.0425179
            # Single   scattering intensity =1.65546e+07 (coh=1.65473e+07 inc=7331.45)
            # Multiple scattering intensity =276313
        elif name in ["sqw", "H2O", "D2O", "quartz"]:
            self.sample_name = "sqw"
            self.sample = mycalculator.add_component(
                self.sample_name,
                "Isotropic_Sqw",
                after=self._sample_arm,
                AT=[0, 0, 0],
                # ROTATED=[0, "a4", 0],
                RELATIVE=self._sample_arm,
            )
            s = self.sample
            s.Sqw_coh = 0
            s.Sqw_inc = 0
            s.sigma_coh = -1
            s.sigma_inc = -1
            # s.xwidth = 0.01
            s.radius = 0.01
            s.yheight = 0.01
            s.thickness = 0.001
            s.set_SPLIT(20)

            if name == "H2O":
                s.Sqw_coh = '"H2O_liq.qSq"'
            elif name == "D2O":
                s.Sqw_coh = '"D2O_liq.qSq"'
            elif name == "quartz":
                s.Sqw_coh = '"SiO2_liq.qSq"'
                s.Sqw_inc = '"SiO2_liq.qSq"'

            else:
                self.calculators[self._calculator_name].add_parameter(
                    "string", "sqw_file", comment="File of the Sqw in McStas convention"
                )
                s.Sqw_coh = "sqw_file"
                self.add_master_parameter(
                    "sqw_file",
                    # here I would need to get the name of the calculator in which the sample is defined
                    {self.calculators[self._calculator_name].name: "sqw_file"},
                )

        # quartz_sample.radius = 0.005
        # quartz_sample.yheight = 0.05
        # # quartz_sample.classical = 0
        # quartz_sample.p_interact = 0.5
        # quartz_sample.Sqw_coh = '"SiO2_liq.qSq"'
        # quartz_sample.Sqw_inc = '"SiO2_liq.qSq"'
        # quartz_sample.sigma_coh = 10.6
        # quartz_sample.sigma_inc = 0.0056
        # quartz_sample.sigma_abs = 0.17
        # quartz_sample.density = 2.2
        # quartz_sample.weight = 60.08
        # # quartz_sample.powder_format = '"qSq"'
        # quartz_sample.verbose = 1
        # quartz_sample.T = 273.21
        # #    quartz_sample.target_index = 3
        # # quartz_sample.append_EXTEND("if(!SCATTERED) ABSORB;")
        # #    Monochromator.

        else:
            raise NameError(f"Sample with name {name} not implemented")

        #        self.__sample_hash = hash(frozenset(vars(self.sample)))
        return self.sample

    # ------------------------------ utility methods made available for the users
    def sim_neutrons(self, number) -> None:
        """Method to set the number of neutrons to be simulated"""
        for calc in self.calculators:
            mycalc = self.calculators[calc]
            mycalc.settings(ncount=number)
