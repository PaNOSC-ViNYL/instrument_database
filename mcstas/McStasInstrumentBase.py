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
        self._sample_holder = None
        self._sample_shape = None
        self._calculator_with_sample = None

        self._sample_hash = None

        self._add_monitors = True
        self.samples = ["None", "vanadium", "H2O", "D2O", "sqw", "monitor"]

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

        mycalculator.add_parameter(
            "double",
            "sample_holder_width",
            comment="width of the sample if shape is box",
            unit="m",
            value=0.0,
        )
        mycalculator.add_parameter(
            "double",
            "sample_holder_height",
            comment="height of the sample if shape is box or cylinder",
            unit="m",
            value=0.0,
        )
        mycalculator.add_parameter(
            "double",
            "sample_holder_depth",
            comment="depth of the sample if shape is box",
            unit="m",
            value=0.0,
        )
        mycalculator.add_parameter(
            "double",
            "sample_holder_radius",
            comment="radius of the sample if shape is sphere or cylinder",
            unit="m",
            value=0.0,
        )
        mycalculator.add_parameter(
            "double",
            "sample_holder_thickness",
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

    def add_multislit(
        self,
        mycalculator: BaseCalculator,
        name: str,
        forms,
        at=0,
        relative="PREVIOUS",
        after=None,
        align=None,
    ):
        """
        Function that simplifies the addition of multiform/multisize collimation holes

        :param mycalculator: the current calculator
        :param name:         base name of the new components, an index will be appended to each
        :param forms:        list of dict with 'x', 'y', and 'r' values for the slits
        :param at:           position according to McStasScript AT convention
        :param relative:     position according to McStasScript RELATIVE convention
        :param after:        as in McstasScript
        :param align:        None or 'b'. 'b' indicates that the rectangular openings are all aligned at the bottom

        Example:
        forms = [
            {"x": 0.035, "y": 0.035, "r": None},
            {"x": 0.020, "y": 0.020, "r": None},
            {"x": 0.015, "y": 0.015, "r": None},
            {"x": None, "y": None, "r": 0.020},
            {"x": None, "y": None, "r": 0.015},
            {"x": None, "y": None, "r": 0.010},
            {"x": 0.025, "y": 0.040, "r": None},
            {"x": 0.030, "y": 0.030, "r": None},
        ]

        """

        disk_index = mycalculator.add_parameter(
            "int",
            "{}_index".format(name),
            comment="Position of {}".format(name),
            value=0,
        )
        allowed_values = []
        allowed_values.extend(range(0, len(forms)))
        disk_index.add_option(allowed_values, True)

        ywidthmin = None
        if align is not None:
            if align != "b":
                raise RuntimeError(
                    "parameter 'align' for add_multislit function not allowed: only None or 'b'"
                )

            for i in range(0, len(forms)):
                size = forms[i]
                if size["r"] is not None:
                    if ywidthmin is None or size["y"] < ywidthmin:
                        ywidthmin = size["y"]

        for i in range(0, len(forms)):
            size = forms[i]
            diaph = mycalculator.add_component(
                "{}_{:d}".format(name, i),
                "Slit",
                AT=at,
                RELATIVE=relative,
                after=after,
                WHEN="{}=={:d}".format(disk_index.name, i),
            )

            if (
                size["y"] is not None
                and ywidthmin is not None
                and size["y"] > ywidthmin
            ):
                ymin = -size["y"] / 2 + ywidthmin
                ymax = -size["y"] / 2 + ywidthmin
                diaph.set_parameters(
                    xwidth=size["x"], ymin=ymin, ymax=ymax, radius=size["r"]
                )
            else:
                diaph.set_parameters(
                    xwidth=size["x"], yheight=size["y"], radius=size["r"]
                )

        return disk_index

    def add_parameter_to_master(
        self, mastername: str, calc: BaseCalculator, par: Parameter
    ) -> None:
        """
        Simplify the addition of a parameter to the master parameters

        :param mastername: name of the master parameter
        :param calc:       current calculator containing the parameter
        :param par:        parameter object defined in the current calculator

        TODO: move to libpyvinyl library. It seems to me general enough and handy
        """

        if mastername not in self.master:
            self.add_master_parameter(
                mastername,
                {calc.name: par.name},
                unit=par.unit,
                comment=par.comment,
            )
            self.master[mastername] = par.value
        else:
            links = self.master[mastername].links
            if calc.name in links:
                print(
                    f"WARNING: trying to add again parameter {par.name} from calculator {calc.name} to master parameter {mastername}"
                )
                raise RuntimeError(
                    f"WARNING: trying to add again parameter {par.name} from calculator {calc.name} to master parameter {mastername}"
                )
                return
            links[calc.name] = par.name
            self.master[mastername].add_links(links)
        self.master[mastername] = self.master[mastername].value

    def sample_focus(self, xwidth, yheight, zdistance):
        # if xwidth is None or yheight is None or zdistance is None:
        if not (isinstance(xwidth, int) or isinstance(xwidth, float)):
            raise RuntimeError(
                "sample_focus xwidth must be int or float instead of {}".format(
                    type(xwidth)
                )
            )

        if not (isinstance(yheight, int) or isinstance(yheight, float)):
            raise RuntimeError(
                "sample_focus yheight must be int or float instead of {}".format(
                    type(yheight)
                )
            )

        if not (
            isinstance(zdistance, int)
            or isinstance(zdistance, float)
            or isinstance(zdistance, Parameter)
        ):
            raise RuntimeError(
                "sample_focus zdistance must be int or float or libpyvinyl.Parameter instead of {}".format(
                    type(zdistance)
                )
            )

        self.focus_xw = xwidth
        self.focus_yh = yheight
        self.target_z = zdistance

    # ------------------------------
    def set_sample_focus(self):
        focus_angle = lambda h, z: 2 * math.atan(h / 2 / z)

        if self.sample is not None and self.sample.component_name == "Isotropic_Sqw":
            dphi = self._calculator_with_sample.add_declare_var(
                "double",
                "sample_focus_dphi",
                comment="focusing dphi for isotropic sample",
            )
            sample_split = self._calculator_with_sample.add_declare_var(
                "double", "sample_split"
            )
            z = self.target_z
            if isinstance(self.target_z, Parameter):
                z = self.target_z.name

            self._calculator_with_sample.append_initialize(
                "sample_focus_dphi = RAD2DEG * 2 * atan(0.5*{h}/{z});".format(
                    h=self.focus_yh, z=z
                )
            )

            self._calculator_with_sample.append_initialize(
                "sample_split = 2 * floor(2*PI/(0.5*{w}/{z}));".format(
                    w=self.focus_xw, z=z
                )
            )

            if self.focus_yh is not None and self.target_z is not None:
                self.sample.d_phi = "sample_focus_dphi"
                # (
                #    180 / math.pi * self.focus_angle(self.focus_yh, self.target_z)
                # )

            if self.focus_xw is not None and self.target_z is not None:
                self.sample.set_SPLIT("sample_split")
                # 2 * round(2 * math.pi / focus_angle(self.focus_xw, self.target_z))
                # )

    # def focus_angle(self, h, z):
    #    return 2 * math.atan(h / 2 / z)

    def _set_sample_shape(
        self, radius: float, width: float, height: float, depth: float, thickness: float
    ) -> None:
        mycalculator = self._calculator_with_sample

        mycalculator.parameters["sample_radius"].value = radius
        mycalculator.parameters["sample_width"].value = width
        mycalculator.parameters["sample_height"].value = height
        mycalculator.parameters["sample_depth"].value = depth
        mycalculator.parameters["sample_thickness"].value = thickness

    def _set_sample_shape_as_holder(self) -> None:
        """
        Helper method to set the sample parameter values to match the shape of the sample holder.

        This is meant to be used in case of liquids and gases.
        """
        mycalculator = self._calculator_with_sample
        if mycalculator.parameters["sample_holder_thickness"].value > 0:
            for v in ["radius", "height", "depth", "width"]:
                if mycalculator.parameters[f"sample_holder_{v}"].value > 0:
                    mycalculator.parameters[f"sample_{v}"].value = (
                        mycalculator.parameters[f"sample_holder_{v}"].value
                        - 2 * mycalculator.parameters["sample_holder_thickness"].value
                    )
                else:
                    mycalculator.parameters[
                        f"sample_{v}"
                    ].value = mycalculator.parameters[f"sample_holder_{v}"].value

        mycalculator.parameters["sample_thickness"].value = 0

    def sample_shape(self, shape: str, r=None, w=None, h=None, d=None, th=0) -> None:
        """
        Sets the sample's shape parameters

         :param shape: accepted values are shere, cylinder, box, holder
         :param r:     None or float value of the radius for shere or cylinder shapes
         :param w:     None or width of box shape
         :param h:     None or height of cylinder and box shapes
         :param d:     None or depth of box shape
         :param th:    thickness of hollow shapes

        Exceptions:
           RuntimeError in case given values are not defined for the selected shape.

        Only parameters related to the selected shapes are used, the others are discarded.
        """
        self._sample_shape = shape
        if shape in ["sphere", "SPHERE"]:
            if r is None or r <= 0:
                raise RuntimeError("Radius of spheric sample (r) should be set > 0")
            self._set_sample_shape(r, 0, 0, 0, th)
        elif shape in ["cylinder", "CYLINDER"]:
            if r is None or r <= 0:
                raise RuntimeError("Radius of cylindric sample (r) should be >0")
            if h is None or h <= 0:
                raise RuntimeError("Height of cylindric sample (h) should be >0")
            self._set_sample_shape(r, 0, h, 0, th)
        elif shape in ["box", "BOX"]:
            if w is None or w <= 0:
                raise RuntimeError("Width of box sample (w) should be >0")
            if h is None or h <= 0:
                raise RuntimeError("Height of box sample (h) should be >0")
            if d is None or d <= 0:
                raise RuntimeError("Depth of box sample (d) should be >0")
            self._set_sample_shape(0, w, h, d, th)
        elif shape in ["holder", "HOLDER"]:
            self._sample_shape = "holder"
            self._set_sample_shape_as_holder()
        else:
            raise RuntimeError(
                "Sample shape not among allowed ones: sphere, box, cylinder, holder"
            )

    def sample_holder(
        self,
        material: str,
        shape: str,
        r: float = None,
        w: float = None,
        h: float = None,
        d: float = None,
        th: float = 0,
    ) -> None:
        """Method to add a sample holder"""

        mycalculator = self._calculator_with_sample
        if material is None or shape is None:
            if self._sample_holder is not None:
                mycalculator.remove_component("sample_holder_in")
                mycalculator.remove_component("sample_holder_out")
                self._sample_holder = None
            return

        s_out = None
        s_in = None

        if self._sample_holder is None:
            self._sample_holder = mycalculator.add_component(
                "sample_holder_in",
                "Isotropic_Sqw",
                AT=[0, 0, 0],
                # ROTATED=[0, "a4", 0],
                RELATIVE=self._sample_arm,
            )
            s_in = self._sample_holder
            # s.append_EXTEND("if(!SCATTERED) ABSORB;")
            s_in.Sqw_coh = 0
            s_in.Sqw_inc = 0
            # s_in.sigma_coh = -1
            # s_in.sigma_inc = -1
            s_in.xwidth = "sample_holder_width"
            s_in.radius = "sample_holder_radius"
            s_in.yheight = "sample_holder_height"
            s_in.zdepth = "sample_holder_depth"
            s_in.thickness = "sample_holder_thickness"
            s_in.verbose = 1
            s_in.concentric = 1
            s_in.p_interact = 1

            s_out = mycalculator.copy_component(
                "sample_holder_out",
                self._sample_holder,
                AT=0,
                RELATIVE=self._sample_arm,
            )
            s_out.concentric = 0
        else:
            s_out = mycalculator.get_component("sample_holder_out")
            s_in = mycalculator.get_component("sample_holder_in")

        if self.sample is None:
            mycalculator.move_component(s_in, after=self._sample_arm)
            mycalculator.move_component(s_out, after=s_in)
        else:
            mycalculator.move_component(s_in, before=self.sample)
            mycalculator.move_component(s_out, after=self.sample)

        # ------------------------------ material
        if material in ["quartz", "QUARTZ"]:
            s_in.Sqw_coh = '"SiO2_liq.qSq"'
            s_out.Sqw_coh = s_in.Sqw_coh

        else:
            raise RuntimeError(
                "Sample holder material " + material + " not implemented"
            )

        # ------------------------------ shape
        if th <= 0:
            raise RuntimeError("Sample holder thickness should be >0: hollow geometry")
        mycalculator.parameters["sample_holder_thickness"].value = th

        if shape in ["shere", "SPHERE"]:
            if r is None or r <= 0:
                raise RuntimeError(
                    "Radius of spheric sample holder (r) should be set > 0"
                )

            mycalculator.parameters["sample_holder_radius"].value = r
            mycalculator.parameters["sample_holder_width"].value = 0
            mycalculator.parameters["sample_holder_height"].value = 0
            mycalculator.parameters["sample_holder_depth"].value = 0

        elif shape in ["cylinder", "CYLINDER"]:
            if r is None or r <= 0:
                raise RuntimeError("Radius of cylindric sample (r) should be >0")
            if h is None or h <= 0:
                raise RuntimeError("Height of cylindric sample (h) should be >0")
            mycalculator.parameters["sample_holder_radius"].value = r
            mycalculator.parameters["sample_holder_width"].value = 0
            mycalculator.parameters["sample_holder_height"].value = h
            mycalculator.parameters["sample_holder_depth"].value = 0
        elif shape in ["box", "BOX"]:
            if w is None or w <= 0:
                raise RuntimeError("Width of box sample (w) should be >0")
            if h is None or h <= 0:
                raise RuntimeError("Height of box sample (h) should be >0")
            if d is None or d <= 0:
                raise RuntimeError("Depth of box sample (d) should be >0")
            mycalculator.parameters["sample_holder_radius"].value = 0
            mycalculator.parameters["sample_holder_width"].value = w
            mycalculator.parameters["sample_holder_height"].value = h
            mycalculator.parameters["sample_holder_depth"].value = d
        else:
            raise RuntimeError(
                "Sample shape not among allowed ones: shere, box, cylinder"
            )

        if self._sample_shape == "holder":
            self._set_sample_shape_as_holder()

            #        print(s_in)
            #        print(s_out)

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
            for p in ["sqw_file", "sqw_inc"]:
                if p in mycalculator.parameters:
                    del mycalculator.parameters[p]
                    del self.master[p]

        if name in ["empty", "Empty", "None", "none"]:
            self.sample = None
            self._sample_shape = None
        elif name in ["monitor"]:
            self.sample_name = "sample_monitor"
            self.sample = mycalculator.add_component(
                self.sample_name,
                "Monitor_nD",
                AT=0,
                RELATIVE=self._sample_arm,
            )
            self.sample.set_parameters(
                xwidth=0.02,
                yheight=0.04,
                zdepth=0.01,
                filename='"sample_monitor.dat"',
                restore_neutron=1,
                options='"box intensity per cm2, bins=1 pressure=0.001, parallel"',
                # options='"x bins={} y bins={} file={}"'.format(1, 1, "counter.dat"
            )

        elif name in ["v_sample"]:
            self.sample_name = "vanadium_old"
            self.sample = mycalculator.add_component(
                self.sample_name,
                "V_sample",
                AT=[0, 0, 0],
                # ROTATED=[0, "a4", 0],
                RELATIVE=self._sample_arm,
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
        elif name in ["sqw", "H2O", "D2O", "quartz", "Vanadium", "vanadium", "qSq"]:
            self.sample_name = "sqw"
            self.sample = mycalculator.add_component(
                self.sample_name,
                "Isotropic_Sqw",
                AT=[0, 0, 0],
                # ROTATED=[0, "a4", 0],
                RELATIVE=self._sample_arm,
            )
            s = self.sample
            # s.append_EXTEND("if(!SCATTERED) ABSORB;")
            s.Sqw_coh = 0
            s.Sqw_inc = 0
            # s.sigma_coh = -1
            # s.sigma_inc = -1
            s.xwidth = "sample_width"
            s.radius = "sample_radius"  # 0.005
            s.yheight = "sample_height"
            s.zdepth = "sample_depth"
            s.thickness = "sample_thickness"  # 0  # 0.002
            s.verbose = 1
            s.p_interact = 1
            s.set_SPLIT(10)
            self.set_sample_focus()

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
            elif name in ["qSq"]:
                s.powder_format = {-1, 0, 0, 0, 0, 0, 1, 0, 0}
                # sqwfile = mycalculator.add_parameter(
                #    "string", "qSq_file", comment="File of the qSqb in McStas convention"
                # )
                sqwfile = mycalculator.add_parameter(
                    "string", "sqw_file", comment="File of the Sqw in McStas convention"
                )
                s.Sqw_coh = "qSq_file"
                s.Sqw_coh = "sqw_file"
                # self.add_parameter_to_master("qSq_file", mycalculator, sqwfile)
                self.add_parameter_to_master("sqw_file", mycalculator, sqwfile)

            else:
                s.sigma_abs = 0
                s.sigma_coh = 0
                s.sigma_inc = 0
                sqwfile = mycalculator.add_parameter(
                    "string", "sqw_file", comment="File of the Sqw in McStas convention"
                )
                sqwinc = mycalculator.add_parameter(
                    "string",
                    "sqw_inc",
                    comment="File with incoherent Sqw in McStas convention (disabled by default)",
                    value='"0"',
                )
                s.Sqw_coh = "sqw_file"
                s.Sqw_inc = "sqw_inc"
                self.add_parameter_to_master("sqw_file", mycalculator, sqwfile)
                self.add_parameter_to_master("sqw_inc", mycalculator, sqwinc)
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

        if self.sample is not None:
            if self._sample_holder is None:
                mycalculator.move_component(self.sample, after=self._sample_arm)
            else:
                mycalculator.move_component(self.sample, after=self._sample_holder)
            # self.__sample_hash = hash(frozenset(vars(self.sample)))
            self.__sample_hash = hash(name)
        else:
            self.__sample_hash = 0
        return self.sample

    def set_sample_environment_by_name(self, name: str) -> None:
        """Adding a sample environment to the simulation"""
        if self._sample_environment_arm is None:
            raise Exception("No sample environment arm defined in the instrument")
        # if self.sample_environment is not None:
        # self.__remove_sample_environment()

        mycalculator = self._calculator_with_sample
        if name in ["empty", "Empty", "None", "none"]:
            self.sample_environment = None
            return

        #### to implement the rest
        mycryo = None
        exit = None
        if name == "10T":
            mycryo, exit = cryo10T(
                mycalculator,
                "10T",
                [0, 0, 0],
                self._sample_environment_arm,
                2,
            )
        else:
            raise NameError("Sample environment name not recognized or not implemented")

        self.sample_environment_name = name

        exit.set_parameters(
            radius=__sample.radius + 1e-6,
            yheight=__sample.yheight + 1e-6,
            priority=100000,
            material_string='"Exit"',
        )

        union_master_after_sample = mycalculator.add_component(
            "master_after_sample",
            "Union_master",
            after=self.__sample,
            AT=[0, 0, 0],
            RELATIVE=mycryo.name,
        )
        union_master_after_sample.allow_inside_start = 1

    def __update_hash(self, calc):
        h = hash(calc)
        if calc == self._calculator_with_sample and self._sample_hash is not None:
            h = h + self._sample_hash + hash(self._sample_shape)
            if self.sample_holder is not None:
                h = h + hash(self.sample_holder.Sqw_coh)
        return hash(h)

    def run(self) -> None:
        # self._check_sample_shape()
        self.custom_flags("-I mcstas/components")
        return super().run()

        for calculator in self.calculators.values():

            self.__update_hash(calculator)
            if hasattr(calculator, "hash"):
                if calculator.hash != hash(calculator):
                    print("HASH: ", calculator.name, hash(calculator))
                    calculator.backengine()
            else:
                calculator.hash = hash(calculator)
                calculator.backengine()
            # print("HASH parameters: ", calculator.name, hash(calculator.parameters))
            # if calculator.name == "OriginCalc":
            #     try:
            #         print(str(calculator.output))
            #     except TypeError as te:
            #         calculator.backengine()
            # else:
            #    calculator.backengine()
        return
        # return super().run()

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
