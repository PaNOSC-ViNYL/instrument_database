# ------------------------------ For McStasscript instruments
import mcstasscript as ms
from mcstasscript.interface import instr
from mcstasscript.helper.mcstas_objects import Component

# ------------------------------ Mandatory classes to use
from libpyvinyl.Instrument import Instrument, BaseCalculator
from libpyvinyl.Parameters import Parameter

# ------------------------------ Extras
import os  # to add the path of custom mcstas components

import math

# list here all the common parts to be imported
from typing import List, Optional, Any


class McStasInstrumentBase(Instrument):
    """:class: BaseClass to be used for creating a good instrument"""

    def __init__(self, name, do_section=True):

        super().__init__(name, instrument_base_dir=".")

        self.__do_section = do_section
        self._temp_directory = "./"  # "/dev/shm/mcstasscript/"
        self._custom_component_dirs = [
            os.path.join(os.path.dirname(__file__), "components")
        ]
        # print(self._custom_component_dirs)
        self._single_calculator = False
        # this is specific for McStasscript instruments:
        # the components of the position for the sample and sample environment
        self._sample_environment_arm = None
        self._sample_arm = None
        self.sample = None
        self._calculator_with_sample = None

        self._sample_hash = None

        self._add_monitors = True
        self.samples = ["None", "vanadium", "H2O", "D2O", "sqw"]

        self.focus_xw = None
        self.focus_yh = None
        self.target_z = None

        self.sample_environments = ["None"]

    def calcLtof(self, fcomp: Component, lcomp: Component):
        L = 0
        for calc in self.calculators:
            L = L + calcLtof(calc, fcomp.name, lcomp.name)
        return L

    def calcLtof(
        self, mycalculator: BaseCalculator, fcomp: str, lcomp: str, debug=False
    ):
        complist = mycalculator.component_list.copy()

        firstfound = False
        L = 0
        for comp in complist[::-1]:
            if debug:
                print(comp.name)
            if comp.name == lcomp:
                firstfound = True
                if debug:
                    print("first found")
            if comp.name == fcomp:
                if not firstfound:
                    raise Exception("Inverted order of arguments")
                break
            if not firstfound:
                continue
            if comp.AT_relative == "ABSOLUTE":
                print(
                    "WARNING: component "
                    + comp.name
                    + " is being ignored because positioned in ABSOLUTE coordinate system"
                )
                continue
            if debug:
                print("Calc TOF for " + comp.name)
            nextcomp = comp.AT_reference
            L = L + math.sqrt(
                comp.AT_data[0] * comp.AT_data[0]
                + comp.AT_data[1] * comp.AT_data[1]
                + comp.AT_data[2] * comp.AT_data[2]
            )
            if debug:
                print("L = " + str(L))
        return L

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
            ROTATED=[0, "sample_y_rotation", 0],
            RELATIVE=previous_component,
        )

        a3 = mycalculator.add_parameter(
            "double",
            "sample_y_rotation",
            comment="sample table rotation angle",
            unit="degree",
            value=0,
        )

        mycalculator.add_parameter(
            "double",
            "sample_width",
            comment="width of the sample if shape is box",
            unit="m",
            value=0.0,
        )
        mycalculator.add_parameter(
            "double",
            "sample_height",
            comment="height of the sample if shape is box or cylinder",
            unit="m",
            value=0.0,
        )
        mycalculator.add_parameter(
            "double",
            "sample_depth",
            comment="depth of the sample if shape is box",
            unit="m",
            value=0.0,
        )
        mycalculator.add_parameter(
            "double",
            "sample_radius",
            comment="radius of the sample if shape is sphere or cylinder",
            unit="m",
            value=0.0,
        )
        mycalculator.add_parameter(
            "double",
            "sample_thickness",
            comment="thickness of the sample if it is hollow",
            unit="m",
            value=0.0,
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
            for d in self._custom_component_dirs:
                mycalculator.add_component_dir(d)

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

        monitor = calculator.add_component(name + "_DEBUG", "Monitor_nD")
        monitor.filename = '"' + name + '_DEBUG"'
        monitor.xwidth = width
        monitor.yheight = height
        monitor.bins = bins
        monitor.options = '"energy limits=[{0:.2f} {1:.2f}] x y time auto, parallel"'.format(
            #        energy - denergy, energy + denergy
            0,
            25,
        )

        # monitor.options = '"energy limits=[Ei-dE, Ei+dE] x y, borders, parallel"'
        monitor.set_AT(position, RELATIVE="PREVIOUS")

        psd = calculator.add_component(name + "_psd_DEBUG", "Monitor_nD")
        psd.filename = '"' + name + '_DEBUG"'
        psd.xwidth = width
        psd.yheight = height
        psd.bins = bins
        psd.options = '"x y, parallel"'
        psd.set_AT([0, 0, 0], RELATIVE="PREVIOUS")

        psd = calculator.add_component(name + "_psdcyl_DEBUG", "Monitor_nD")
        psd.filename = '"' + name + '_DEBUG"'
        psd.xwidth = 2 * radius
        psd.yheight = height
        psd.options = '"theta bins=360, y bins=20, cylinder, parallel"'
        psd.set_AT([0, 0, 0], RELATIVE="PREVIOUS")

        return psd

    def set_sample_focus(self, xwidth, yheight, zdistance):
        self.focus_xw = xwidth
        self.focus_yh = yheight
        self.target_z = zdistance

    def focus_angle(self, h, z):
        return 2 * math.atan(h / 2 / z)

    def _set_sample_shape(
        self, radius: float, width: float, height: float, depth: float, thickness: float
    ) -> None:
        mycalculator = self._calculator_with_sample
        mycalculator.parameters["sample_radius"].value = radius
        mycalculator.parameters["sample_width"].value = width
        mycalculator.parameters["sample_height"].value = height
        mycalculator.parameters["sample_depth"].value = depth
        mycalculator.parameters["sample_thickness"].value = thickness

    def sample_sphere_shape(self, radius: float, thickness: float = 0) -> None:
        if radius <= 0:
            raise RuntimeError("Radius of spheric sample should be >0")
        self._set_sample_shape(radius, 0, 0, 0, thickness)

    def sample_cylinder_shape(
        self, radius: float, height: float, thickness: float = 0
    ) -> None:
        if radius <= 0:
            raise RuntimeError("Radius of cylindric sample should be >0")
        if height <= 0:
            raise RuntimeError("Height of cylindric sample should be >0")
        self._set_sample_shape(radius, 0, height, 0, thickness)

    def sample_box_shape(
        self, width: float, height: float, depth: float, thickness: float = 0
    ) -> None:
        if width <= 0:
            raise RuntimeError("Width of box sample should be >0")
        if height <= 0:
            raise RuntimeError("Height of box sample should be >0")
        if depth <= 0:
            raise RuntimeError("Depth of box sample should be >0")

        self._set_sample_shape(0, width, height, depth, thickness)

    def _check_sample_shape(self):
        mycalculator = self._calculator_with_sample

        r = mycalculator.parameters["sample_radius"].value
        w = mycalculator.parameters["sample_width"].value
        # h = mycalculator.parameters["sample_height"].value
        # d = mycalculator.parameters["sample_depth"].value
        # t = mycalculator.parameters["sample_thickness"].value

        if r <= 0 and w <= 0:
            raise RuntimeError(
                "Sample shape should be defined by: sample_sphere_shape or sample_cylinder_shape or sample_box_shape methods"
            )

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
        elif name in ["v_sample"]:
            self.sample_name = "vanadium_old"
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
        elif name in ["sqw", "H2O", "D2O", "quartz", "Vanadium", "vanadium"]:
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
            # s.append_EXTEND("if(!SCATTERED) ABSORB;")
            s.Sqw_coh = 0
            s.Sqw_inc = 0
            s.sigma_coh = -1
            s.sigma_inc = -1
            s.xwidth = "sample_width"
            s.radius = "sample_radius"  # 0.005
            s.yheight = "sample_height"
            s.thickness = "sample_thickness"  # 0  # 0.002
            s.verbose = 2
            s.p_interact = 1
            s.d_phi = 180 / math.pi * self.focus_angle(self.focus_yh, self.target_z)

            s.set_SPLIT(
                2 * round(2 * math.pi / self.focus_angle(self.focus_xw, self.target_z))
            )

            if name == "H2O":
                s.Sqw_coh = '"H2O_liq_290_coh.sqw"'
                s.Sqw_inc = '"H2O_liq_290_inc.sqw"'
                s.sigma_coh = 7.75
                s.sigma_inc = 161
            elif name == "D2O":
                s.Sqw_coh = '"D2O_liq_290_coh.sqw"'
                s.Sqw_inc = '"D2O_liq_290_inc.sqw"'
                s.sigma_coh = 15.411
                s.sigma_inc = 4.1008

            elif name == "quartz":
                s.Sqw_coh = '"SiO2_liq.qSq"'
                s.Sqw_inc = '"SiO2_liq.qSq"'
            elif name in ["Vanadium", "vanadium"]:
                s.rho = 1 / 13.827
                s.sigma_abs = 5.08
                s.sigma_inc = 4.935
                s.sigma_coh = 0
            else:
                s.sigma_abs = 0
                s.sigma_coh = 0
                s.sigma_inc = 0
                mycalculator.add_parameter(
                    "string", "sqw_file", comment="File of the Sqw in McStas convention"
                )
                mycalculator.add_parameter(
                    "string",
                    "sqw_inc",
                    comment="File with incoherent Sqw in McStas convention (disabled by default)",
                    value='"0"',
                )
                s.Sqw_coh = "sqw_file"
                s.Sqw_inc = "sqw_inc"
                self.add_master_parameter(
                    "sqw_file",
                    # here I would need to get the name of the calculator in which the sample is defined
                    {mycalculator.name: "sqw_file"},
                    comment=mycalculator.parameters["sqw_file"].comment,
                )
                self.add_master_parameter(
                    "sqw_inc",
                    {mycalculator.name: "sqw_inc"},
                    comment=mycalculator.parameters["sqw_file"].comment,
                )
                self.master["sqw_inc"] = '"0"'

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

    def run(self):
        self._check_sample_shape()
        self.custom_flags("-I mcstas/components")

        return super().run()

    # ------------------------------ utility methods made available for the users
    def sim_neutrons(self, number) -> None:
        """Method to set the number of neutrons to be simulated"""
        for mycalc in self.calculators.values():
            mycalc.settings(ncount=number)

    def set_seed(self, number) -> None:
        """Setting the simulation seed"""
        for mycalc in self.calculators.values():
            mycalc.settings(seed=number)

    def force_compile(self, force_compile: bool) -> None:
        """Setting the force compile on all the calculators"""
        for mycalc in self.calculators.values():
            mycalc.settings(force_compile=force_compile)

    def custom_flags(self, flags: str) -> None:
        """Additional custom flags for McRun"""
        for mycalc in self.calculators.values():
            mycalc.settings(custom_flags=flags)

    def get_total_SPLIT(self):
        split = 1
        for mycalc in self.calculators.values():
            # print("Calculator: ", calc)
            # mycalc = self.calculators[calc]
            for component in mycalc.component_list:
                # print(component.name, component.SPLIT)

                if component.SPLIT > 1:
                    split = split * component.SPLIT
        return split
