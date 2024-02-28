"""

 - The monochromator do not allow transmitted neutrons

TODO: 
 - [ ] Check RV RH for monochromator and analyzer
 - [ ] Understand the A3_offset:
        ThALES.add_parameter("double", "q_x_elastic", value=1.3139)
        ThALES.add_parameter("double", "q_z_elastic", value=0.146)
        ThALES.append_initialize("  A3_offset=atan(q_z_elastic/q_x_elastic)*RAD2DEG; ")
 - [ ] Hashing of the sample object to check differences and trigger a recompilation
"""

# ------------------------------ Mandatory classes to use
from libpyvinyl.Instrument import Instrument
from libpyvinyl.Parameters import Parameter

# ------------------------------ For McStasscript instruments
import mcstasscript as ms
from mcstasscript.interface import functions
from mcstasscript.interface import instr

from mcstas.McStasInstrumentBase import McStasInstrumentBase

# this is needed to get the location of McStas executables and libraries
my_configurator = functions.Configurator()

# ------------------------------ Importing sources
# from institutes.ILL.sources.HEAD.mcstas import Full as source
# from institutes.ILL.sources.HEAD.mcstas import Gauss as source
from institutes.ILL.sources.HEAD.mcstas import Monochromator_source as source


# ------------------------------ Extras
import os  # to add the path of custom mcstas components

# for operations
import math

# list here all the common parts to be imported
from typing import List, Optional, Any

# for unit conversions
import pint
from pint import set_application_registry

ureg = pint.get_application_registry()


############## Mandatory method
def get_flavours():
    return ["full", "nosection", "quick", "quicknosection"]


############## Mandatory method
def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    if flavour in [None, "None", "", "full"]:
        return ThALES()
    if flavour == "nosection":
        return ThALES(False)
    if flavour in ["quick"]:
        return ThALES(True, True)
    if flavour in ["quicknosection"]:
        return ThALES(False, True)
    else:
        raise RuntimeError(f"Flavour {flavour} not implement")


class ThALES(McStasInstrumentBase):
    """:class: Instrument class defining the ThALES instrument at ILL"""

    # ------------------------------ utility methods made available for the users

    def wavelength_to_angle(self, value: [float, pint.Quantity]) -> pint.Quantity:
        """Conversion from wavelength to angle
        for Bragg law with lattice parameter
        equal to the monochromator lattice"""
        d_lattice = self.parameters["OriginCalc"]["monochromator_d"].pint_value
        return (math.asin(value / 2.0 / d_lattice) * 2 * ureg.radians).to("degrees")

    def energy_to_angle(self, value: [float, pint.Quantity]) -> pint.Quantity:
        """Conversion from energy to angle
        for Bragg law with lattice parameter
        equal to the monochromator lattice"""
        d_lattice = self.parameters["OriginCalc"]["monochromator_d"].pint_value
        wl = math.sqrt(81.80421 * ureg.meV / value) * ureg.angstrom
        return self.wavelength_to_angle(wl).to("degrees")

    # ------------------------------ Specialized methods required by the Instrument class [mandatory]

    #    def run(self):
    #        force_compile = self.calculators[self._calculator_name]._run_settings[
    #            "force_compile"
    #        ]
    #        print("Presettings: ")
    #        self.calculators[self._calculator_name].show_settings()
    #        print(self.__sample_hash)

    #        if hash(frozenset(vars(self.sample))) != self.__sample_hash:
    #            self.calculators[self._calculator_name].settings(force_compile=True)
    # super().run()
    #        print("Updated settings: ")
    #        self.calculators[self._calculator_name].show_settings()
    #        self.calculators[self._calculator_name].settings(force_compile=force_compile)

    #        print("Restored settings:")
    #        self.calculators[self._calculator_name].show_settings()

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

    # ------------------------------ Internal methods (not available to users)
    # def __remove_sample_environment(self) -> None:
    #    """Remove any previously defined sample environment"""
    #    if self.__sample_environment is not None:
    #        mycalculator = self.calculators[-1]  # there is only 1
    #        mycalculator.remove_component(self.__sample_environment)

    #    def set_instrument_base_dir(self, dirname) -> None:
    #        """ """
    #        for calc in self.calculators:
    #            c = self.calculators[calc]
    #            c.input_path = dirname
    #        super().set_instrument_base_dir(dirname)

    # ------------------------------ The instrument definition goes in the __init__
    def __init__(self, do_section=True, _start_from_virtual_source=False):
        """Here the real definition of the instrument is performed"""

        super().__init__("ThALES", do_section)

        self.sample_environments = ["None", "10T", "Orange"]

        #
        # self._calculator_name = "Origin"
        # ------------------------------ some local variables
        myinstr = self

        gR0 = 1
        gQc = 0.0216
        gAlpha = 4.07
        gW = 1.0 / 300.0
        Al_Thickness = 0.001
        gGap = 0.001

        ThALES_L = 2.000  # distance between monochromator and sample
        virtsource_mono_distance = 7
        dist_mono_sample = 2.000
        dist_sample_ana = 1.260  # distance between sample and analyzer
        dist_ana_det = 0.640  # distance between analyzer and detector

        # ------------------------------------------------------------
        # Start with a first section and declaring its parameters
        mycalculator, Origin = self.add_new_section("OriginCalc")
        OriginCalc = mycalculator

        a2 = mycalculator.add_parameter(
            "double",
            "a2",
            comment="Angle between beam reflected by monochromator and incident beam",
            unit="degree",
            value=33,
        )
        a2.add_interval(33, 128, True)
        self.add_parameter_to_master("a2", mycalculator, a2)

        monochromator_d = mycalculator.add_parameter(
            "double",
            "monochromator_d",
            comment="Monochromator lattice parameter",
            value=3.355,
            unit="angstrom",
        )

        RHmono = mycalculator.add_parameter(
            "double",
            "RHmono",
            comment="Monochromator horizontal focusing. Calculated by default",
            value=-1,
            unit="",
        )
        # setting default value for RHmono if not provided
        mycalculator.append_initialize(
            f"if(RHmono<0) RHmono= (1/( ( 1/{dist_mono_sample} + 1/{virtsource_mono_distance} ) *sin(DEG2RAD*a2/2) / 2 ));"
        )

        RVmono = mycalculator.add_parameter(
            "double",
            "RVmono",
            comment="Monochromator horizontal focusing. Calculated by default",
            value=-1,
            unit="",
        )
        # setting default value for RVmono if not provided
        mycalculator.append_initialize(
            f"if(RVmono<0) RVmono= (1/( ( 1/{dist_mono_sample} + 1/{virtsource_mono_distance} ) / (2*sin(DEG2RAD*a2/2) )));"
        )
        mycalculator.append_initialize(
            'printf("(RHmono,RVmono) = (%.2f,%.2f)\\n", RHmono, RVmono);'
        )

        # imported source requires the following variables to be defined for the value of lambda:
        # "lambda =2*sin(0.5*a2*DEG2RAD)*monochromator_d;"
        # - a2
        # - monochromator_d
        #
        HCS = source.HCS_source(mycalculator)
        HCS.flux = 16e10

        mycalculator.parameters["dlambda"] = 0.08  # imported from HCS

        # ------------------------------------------------------------
        if not _start_from_virtual_source:
            HCS.dist = 2.155
            HCS.focus_xw = 0.06
            HCS.focus_yh = 0.12

            # start of the guide
            H5 = mycalculator.add_component("H5", "Arm")
            H5.set_AT([0, 0, 2.155], RELATIVE=HCS)

            # adds monitors at the same position of the previous component

            H53_1a = mycalculator.add_component("H5_1", "Guide_gravity")
            # material: Al 5083
            H53_1a.w1 = 0.060
            H53_1a.h1 = 0.120
            H53_1a.l = 1.0
            #    H53_1a.R0 = gR0
            #    H53_1a.Qc = gQc
            #    H53_1a.alpha = gAlpha
            #    H53_1a.m = 3
            #    H53_1a.W = gW
            H53_1a.reflect = '"supermirror_m3.rfl"'
            H53_1a.set_AT([0, 0, 0], RELATIVE=H5)

            H53_1b = mycalculator.copy_component("H5_1b", H53_1a)
            H53_1b.l = 1.775
            H53_1b.set_AT([0, 0, H53_1a.l], RELATIVE=H53_1a)

            # membrane
            # gap 25 mm
            H53_A = mycalculator.add_component("H53_A", "Arm")
            H53_A.set_AT([0, 0, H53_1b.l + 0.025 + 0.006], RELATIVE=H53_1b)

            H53_2a = mycalculator.copy_component("H53_2a", H53_1a)
            H53_2a.l = 2.317
            H53_2a.set_AT([0, 0, 0], RELATIVE=H53_A)

            H53_2b = mycalculator.copy_component("H53_2b", H53_1a)
            H53_2b.l = 0.757
            H53_2b.set_AT([0, 0, H53_2a.l], RELATIVE=H53_2a)

            H53_B = mycalculator.add_component("H53_B", "Arm")
            H53_B.set_AT([0, 0, H53_2b.l + 0.060 + 0.027], RELATIVE=H53_2b)

            H53_3 = mycalculator.copy_component("H53_3", H53_1a)
            H53_3.l = 0.501
            H53_3.set_AT([0, 0, 0], RELATIVE=H53_B)

            H53_4 = mycalculator.copy_component("H53_4", H53_1a)
            H53_4.l = 0.501
            H53_4.set_AT([0, 0, H53_3.l], RELATIVE=H53_3)

            H53_5 = mycalculator.copy_component("H53_5", H53_1a)
            H53_5.l = 5.809
            H53_5.set_AT([0, 0, H53_4.l], RELATIVE=H53_4)

            H53_C = mycalculator.add_component("H53_C", "Arm")
            H53_C.set_AT([0, 0, H53_5.l + 0.0078 + 0.024 + 0.110], RELATIVE=H53_5)

            H53_6 = mycalculator.copy_component("H53_6", H53_1a)
            H53_6.l = 2.420
            H53_6.set_AT([0, 0, 0], RELATIVE=H53_C)

            H53_D = mycalculator.add_component(
                "H53_D", "Arm", AT=[0, 0, H53_6.l + 0.0195 + 0.002], RELATIVE=H53_6
            )

            # adds monitors at the same position of the previous component

            H53_7 = mycalculator.add_component("H53_7_parabolic", "Guide_tapering")
            H53_7.option = '"elliptical"'
            H53_7.w1 = 0.06
            H53_7.h1 = 0.12
            H53_7.l = 2.4915
            H53_7_ah = 2.737445  # major semi-axis
            H53_7_bh = H53_7.w1 / 2.0  # minor semi-axis
            H53_7.linw = math.sqrt(H53_7_ah * H53_7_ah - H53_7_bh * H53_7_bh)
            H53_7.loutw = (
                H53_7.linw - H53_7.l
            )  # the length of the guide must be smaller than the coordinate of the focus
            #    H53_7_linw = 0
            #    H53_7_linh = 0
            H53_7.linh = 0  # plane mirrors (vertical)
            H53_7.louth = 0  #  plane mirrors (vertical)
            H53_7.R0 = gR0
            H53_7.Qcx = gQc
            H53_7.Qcy = gQc
            H53_7.alphax = gAlpha
            H53_7.alphay = gAlpha
            H53_7.W = gW
            H53_7.mx = 3
            H53_7.my = 3
            H53_7.segno = 1  # is this the number of segments?
            H53_7.set_AT([0, 0, 0], RELATIVE=H53_D)

            H53_7out = mycalculator.add_component("H53_7out", "Arm")
            H53_7out.set_AT([0, 0, H53_7.l], RELATIVE=H53_7)

        else:
            HCS.dist = 20
            HCS.focus_xw = 0.04
            HCS.focus_yh = 0.12
            H53_7out = mycalculator.add_component(
                "H53_7out", "Arm", AT=20, RELATIVE=HCS
            )

        # adds monitors at the same position of the previous component

        slit_A = mycalculator.add_component(
            "slit_A", "Slit", AT=0.3, RELATIVE="PREVIOUS"
        )
        slit_A.xwidth = 0.04
        slit_A.yheight = 0.12

        # monitor

        Monochromator_Arm = mycalculator.add_component(
            "Monochromator_Arm", "Arm", AT=2, RELATIVE=slit_A
        )

        # double check the probability of reflection
        Monochromator = mycalculator.add_component(
            "Monochromator", "Monochromator_curved", AT=0, RELATIVE=Monochromator_Arm
        )
        #   Monochromator.reflect = '"HOPG.rfl"'
        #   Monochromator.transmit = '"HOPG.trm"'

        Monochromator.set_parameters(
            gap=0.0005,
            NH=13,
            NV=13,
            mosaich=30,
            mosaicv=30,
            r0=1,
            t0=0,  # remove transmitted neutrons
            RV="RVmono",
            RH="RHmono",
            DM=monochromator_d,
            width=0.25,
            height=0.2,
            verbose=1,
        )
        # Monochromator.append_EXTEND("if(flag!=SCATTERED) ABSORB;")
        Monochromator.set_ROTATED([0, "a2/2", 0], RELATIVE=Monochromator_Arm)

        beamstop = mycalculator.add_component(
            "BS", "Beamstop", AT=[0, 0, 2], RELATIVE=Monochromator_Arm
        )
        beamstop.set_parameters(xwidth=1, yheight=1)

        Monochromator_Out = mycalculator.add_component("Monochromator_Out", "Arm")
        Monochromator_Out.set_AT([0, 0, 0], RELATIVE=Monochromator_Arm)
        Monochromator_Out.set_ROTATED([0, "a2", 0], RELATIVE=Monochromator_Arm)

        counter = mycalculator.add_component(
            "counter", "Monitor_nD", AT=[0.0, 0.0, 0.8], RELATIVE=Monochromator_Out
        )
        counter.set_parameters(
            xwidth=0.01,
            yheight=0.01,
            zdepth=0.01,
            filename='"counter.dat"',
            restore_neutron=1,
            options='"box intensity, bins=1 pressure=0.001"',
        )

        before_sample_slit = mycalculator.add_component("slit_before_sample", "Slit")
        before_sample_slit.xwidth = 0.03
        before_sample_slit.yheight = 0.028
        before_sample_slit.set_AT([0, 0, ThALES_L - 0.250], RELATIVE=Monochromator_Out)

        sample_mcpl_arm = mycalculator.add_component(
            "sample_mcpl_arm",
            "Arm",
            AT=[0, 0, 0],
            RELATIVE=before_sample_slit,
        )

        # ------------------------------------------------------------
        # this new section contains the sample and the sample environment
        SampleCalc, sample_mcpl_arm = self.add_new_section(
            "SampleCalc", sample_mcpl_arm, True
        )
        # ------------------------------------------------------------

        a4 = SampleCalc.add_parameter(
            "double",
            "a4",
            comment="Angle between reflected by sample and incident beams",
            unit="degree",
            value=0,
        )
        a4.add_interval(-128, 128, True)
        self.add_parameter_to_master("a4", SampleCalc, a4)

        self._sample_arm.set_AT(0.250, RELATIVE=sample_mcpl_arm)
        self._sample_environment_arm.set_AT(0.250, RELATIVE=sample_mcpl_arm)

        # default sample
        self.sample_focus(0.03, 0.04, 0.250)

        sample = self.set_sample_by_name("Vanadium")

        Sample_Out = SampleCalc.add_component(
            "Sample_Out", "Arm", AT=0, RELATIVE=self._sample_arm
        )
        Sample_Out.set_ROTATED([0, "a4", 0], RELATIVE=self._sample_arm)

        slit_distance = 0.250
        after_sample_slit = SampleCalc.add_component(
            "slit_after_sample", "Slit", AT=slit_distance, RELATIVE=Sample_Out
        )
        after_sample_slit.xwidth = 0.03
        after_sample_slit.yheight = 0.04

        # ------------------------------------------------------------
        AnalyzerCalc, after_sample_slit = self.add_new_section(
            "AnalyzerCalc", after_sample_slit
        )
        # ------------------------------------------------------------

        a6 = AnalyzerCalc.add_parameter(
            "double",
            "a6",
            comment="Angle between reflected by analyzer and incident beams",
            unit="degree",
            value=33,
        )
        self.add_parameter_to_master("a6", AnalyzerCalc, a6)

        RHanalyzer = AnalyzerCalc.add_parameter(
            "double",
            "RHanalyzer",
            comment="Monochromator horizontal focusing. Calculated by default",
            value=-1,
            unit="",
        )
        # setting default value for RHmono if not provided
        ana_focus = 1.0 / (1.0 / dist_sample_ana + 1.0 / dist_ana_det)
        AnalyzerCalc.append_initialize(
            f"if(RHanalyzer<0) RHanalyzer= 2 * {dist_sample_ana} / sin(DEG2RAD*a6/2);"
        )

        RVanalyzer = AnalyzerCalc.add_parameter(
            "double",
            "RVanalyzer",
            comment="Monochromator horizontal focusing. Calculated by default",
            value=-1,
            unit="",
        )
        # setting default value for RVmono if not provided
        AnalyzerCalc.append_initialize(
            f"if(RVanalyzer<0) RVanalyzer= 2 * {ana_focus} * sin(DEG2RAD*a6/2) ;"
        )
        AnalyzerCalc.append_initialize(
            'printf("(RHanalyzer,RVanalyzer) = (%.2f,%.2f)\\n", RHanalyzer, RVanalyzer);'
        )

        Ana_Cradle = AnalyzerCalc.add_component("Ana_Cradle", "Arm")
        Ana_Cradle.set_AT(
            [0, 0, dist_sample_ana - slit_distance], RELATIVE=after_sample_slit
        )

        analyzer = AnalyzerCalc.add_component("analyzer", "Monochromator_curved")
        analyzer.gap = 0.0005
        analyzer.NH = 11
        analyzer.NV = 9
        analyzer.mosaich = 30
        analyzer.mosaicv = 30
        analyzer.r0 = 0.7
        analyzer.RV = RVanalyzer
        analyzer.RH = RHanalyzer
        analyzer.DM = 3.355  # PG 002
        analyzer.width = 0.17
        analyzer.height = 0.13
        Monochromator.verbose = 1
        analyzer.set_AT([0, 0, 0], RELATIVE=Ana_Cradle)
        analyzer.set_ROTATED([0, "a6*0.5", 0], RELATIVE=Ana_Cradle)

        Ana_Out = AnalyzerCalc.add_component(
            "Ana_Out", "Arm", AT=[0, 0, 0], ROTATED=[0, "a6", 0], RELATIVE=Ana_Cradle
        )

        slit_distance = 0.340
        slit = AnalyzerCalc.add_component("slit", "Slit")
        slit.xwidth = 0.03
        slit.yheight = 0.08
        slit.set_AT([0, 0, slit_distance], RELATIVE=Ana_Out)

        # ------------------------------------------------------------
        DetectorCalc, slit = self.add_new_section("DetectorCalc", slit)
        # ------------------------------------------------------------
        detector_arm = DetectorCalc.add_component(
            "detector_arm",
            "Arm",
            AT=[0, 0, dist_ana_det - slit_distance],
            RELATIVE=slit,
        )

        detector_all = DetectorCalc.add_component(
            "detector_all", "Monitor_nD", AT=[0, 0, 0], RELATIVE=detector_arm
        )
        detector_all.xwidth = 0.05
        detector_all.yheight = 0.12
        detector_all.restore_neutron = 0
        detector_all.options = '"intensity, square bins=1, file=detector_all.dat"'

        # ------------------------------ instrument parameters
        self.add_parameter_to_master(
            "a3", SampleCalc, SampleCalc.parameters["sample_y_rotation"]
        )

        myinstr.master["a2"] = 79.10 * ureg.degree
        myinstr.master["a3"] = 0 * ureg.degree
        myinstr.master["a4"] = 60 * ureg.degree
        myinstr.master["a6"] = 79.10 * ureg.degree
        # ------------------------------ sample parameters
        # Do not add sample parameters. They should be modified externally retrieving
        # sample with .sample
        # this obviously will require the instrument to be recompiled


# ------------------------------ Helper functions


def cryoOrange(
    instrument,
    name="cryoOrange",
    position=[0, 0, 0],
    relative="PREVIOUS",
    number_of_activations=1,
):

    if addCryostat is False:
        mycryo = instrument.add_component(name, "Arm")
        mycryo.set_AT(position, RELATIVE=relative, after=relative)
        return mycryo, None

    mycryo = ms.Cryostat(name, instrument)
    mycryo.set_AT(position, RELATIVE=relative, after=relative)

    mycryo.add_layer(
        inner_radius=150e-3 / 2.0,  # diameter: 150 mm
        thickness=1.0e-3,  # 1 mm
        origin_to_top=150e-3 / 2,  # guessed
        top_thickness=3e-3,
        origin_to_bottom=150e-3 / 2,  # guessed
        bottom_thickness=3e-3,
        p_interact=0.2,
    )

    #    mycryo.last_layer.add_window(
    #        inner_radius=mycryo.last_layer.inner_radius,
    #        thickness=1.0 / 1000,
    #        origin_to_top=11.0 / 1000,  # guessed
    #        origin_to_bottom=11.0 / 1000,  # guessed
    #    )

    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_layer(
        inner_radius=159e-3 / 2.0,  # diameter: 150 mm
        thickness=0.5e-3,  # 1 mm
        origin_to_top=150e-3 / 2,  # guessed
        top_thickness=3e-3,
        origin_to_bottom=150e-3 / 2,  # guessed
        bottom_thickness=3e-3,
        p_interact=0.2,
    )
    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_layer(
        inner_radius=170e-3 / 2.0,  # diameter: 150 mm
        thickness=1e-3,  # 1 mm
        origin_to_top=150e-3 / 2,  # guessed
        top_thickness=3e-3,
        origin_to_bottom=150e-3 / 2,  # guessed
        bottom_thickness=3e-3,
        p_interact=0.2,
    )
    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_layer(
        inner_radius=280e-3 / 2.0,  # diameter: 150 mm
        thickness=1e-3,  # 1 mm
        origin_to_top=150e-3 / 2,  # guessed
        top_thickness=3e-3,
        origin_to_bottom=150e-3 / 2,  # guessed
        bottom_thickness=3e-3,
        p_interact=0.2,
    )
    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_spatial_loggers()
    mycryo.build(include_master=False)

    exit = instrument.add_component(
        "exit_volume", "Union_cylinder", AT=[0, 0, 0], RELATIVE=name
    )

    master_before_sample = instrument.add_component(
        "master_before_sample", "Union_master", RELATIVE=name
    )

    return mycryo, exit


def cryo10T(
    instrument,
    name="10T",
    position=[0, 0, 0],
    relative="PREVIOUS",
    number_of_activations=1,
):
    """Create the cryostat at the sample_arm position previously defined"""
    if addCryostat is False:
        mycryo = instrument.add_component(name, "Arm")
        mycryo.set_AT(position, RELATIVE=relative)
        return mycryo, None

    mycryo = ms.Cryostat(name, instrument)
    mycryo.set_AT(position, RELATIVE=relative)

    mycryo.add_layer(
        inner_radius=50.8 / 2.0 / 1000,  # diameter: 50.8mm
        thickness=2.0 / 1000,  # 1 mm
        origin_to_top=150.0 / 1000,  # guessed
        top_thickness=0.001,
        origin_to_bottom=150.0 / 1000,  # guessed
        bottom_thickness=0.001,
        p_interact=0.2,
    )

    mycryo.last_layer.add_window(
        inner_radius=mycryo.last_layer.inner_radius,
        thickness=1.0 / 1000,
        origin_to_top=11.0 / 1000,  # guessed
        origin_to_bottom=11.0 / 1000,  # guessed
    )

    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_spatial_loggers()
    mycryo.build(include_master=False)

    exit = instrument.add_component("exit_volume", "Union_cylinder")
    exit.set_AT([0, 0, 0], RELATIVE=name)

    master_before_sample = instrument.add_component(
        "master_before_sample", "Union_master", RELATIVE=name
    )

    return mycryo, exit

    ################### ignored for the moment
    mycryo.add_layer(
        inner_radius=592.0 / 2.0 / 1000,  # diameter: 592 mm
        thickness=3.0 / 1000,  # 3 mm
        origin_to_top=140.0 / 1000,  # guessed
        top_thickness=0.001,
        origin_to_bottom=140.0 / 1000,  # guessed
        bottom_thickness=0.001,
        p_interact=0.2,
        material="Al",
    )
    mycryo.last_layer.add_window(
        inner_radius=mycryo.last_layer.inner_radius,
        thickness=1.0 / 1000,
        origin_to_top=51.0 / 1000,  # guessed
        origin_to_bottom=51.0 / 1000,  # guessed
    )
    #    mycryo.add_spatial_loggers()
    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.build(include_master=False)

    return mycryo
