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

# ------------------------------ For McStasscript instruments
import mcstasscript as ms
from mcstasscript.interface import functions
from mcstasscript.interface import instr

# this is needed to get the location of McStas executables and libraries
my_configurator = functions.Configurator()

# ------------------------------ Importing sources
# from institutes.ILL.sources.HEAD.mcstas import Full as source
from institutes.ILL.sources.HEAD.mcstas import Gauss as source

# ------------------------------ Mandatory classes to use
from libpyvinyl.Instrument import Instrument
from libpyvinyl.Parameters import Parameter
from institutes.ILL.instruments.ThALES.HEAD.mcstas.McStasInstrumentBase import (
    McStasInstrumentBase,
)

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


## Global variables for DEBUG to be removed!
addMonitors = True
addCryostat = True


############## Mandatory method
def get_flavours():
    return ["full", "from_sample", "merge", "nosection"]


############## Mandatory method
def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    if flavour in [None, "None", "", "full"]:
        return ThALES()
    elif flavour == "from_sample":
        return ThALES_from_sample()
    elif flavour == "merge":
        return ThALES_merge()
    elif flavour == "nosection":
        return ThALES(False)
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
                ROTATED=[0, "a4", 0],
                RELATIVE=self._sample_arm,
                after=self._sample_arm,
            )
            v_sample = self.sample
            v_sample.radius = 0.01
            v_sample.yheight = 0.01
            v_sample.thickness = 0.001
            v_sample.focus_xw = 0.04
            v_sample.focus_yh = 0.12
            v_sample.target_z = 0.25
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
                ROTATED=[0, "a4", 0],
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
    def __init__(self, do_section=True):
        """Here the real definition of the instrument is performed"""

        super().__init__("ThALES")

        self.samples = ["None", "vanadium", "H2O", "D2O", "sqw"]
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
        dist_sample_ana = 1.260  # distance between sample and analyzer
        dist_ana_det = 0.640  # distance between analyzer and detector

        # ------------------------------------------------------------
        # Start with a first section and declaring its parameters
        mycalculator, Origin = self.add_new_section("OriginCalc")

        a2 = mycalculator.add_parameter(
            "double",
            "a2",
            comment="Angle between beam reflected by monochromator and incident beam",
            unit="degree",
            value=33,
        )
        a2.add_interval(33, 128, True)

        monochromator_d = mycalculator.add_parameter(
            "double",
            "monochromator_d",
            comment="Monochromator lattice parameter",
            value=3.355,
            unit="angstrom",
        )

        # imported source and associated parameters: check the source file!
        # - lambda
        # - dlambda
        # - Ei
        # - dE
        HCS = source.HCS_source(mycalculator)
        HCS.E0 = "Ei"
        HCS.target_index = 1
        HCS.flux = 1e8

        # override the value of lambda by the value of the angle
        mycalculator.append_initialize("lambda =2*sin(0.5*a2*DEG2RAD)*monochromator_d;")
        mycalculator.append_initialize('printf("lambda: %.2f\\n", lambda);')
        mycalculator.append_initialize("dlambda=0.01*lambda;")
        # the following python string is an alternative way of doing:
        mycalculator.add_declare_var("double", "energy")
        mycalculator.append_initialize("Ei = 81.80421/(lambda*lambda);")
        mycalculator.add_declare_var("double", "denergy")
        mycalculator.append_initialize("dE = 2*Ei*dlambda/lambda;")
        mycalculator.append_initialize(
            'printf("lambda: %.2f +/- %.2f\\n", lambda, dlambda);'
        )
        mycalculator.append_initialize('printf("energy: %.2f +/- %.2f\\n", Ei, dE);')

        # ------------------------------------------------------------
        # start of the guide
        H5 = mycalculator.add_component("H5", "Arm")
        H5.set_AT([0, 0, 2.155], RELATIVE=HCS)

        # adds monitors at the same position of the previous component
        addMonitor(mycalculator, "H5")

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

        addMonitor(mycalculator, "H5_1a", H53_1a.l)

        H53_1b = mycalculator.copy_component("H5_1b", H53_1a)
        H53_1b.l = 1.775
        H53_1b.set_AT([0, 0, H53_1a.l], RELATIVE=H53_1a)

        addMonitor(mycalculator, "H5_1b", [0, 0, H53_1b.l])

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
        addMonitor(mycalculator, "H53_D")

        H53_7 = mycalculator.add_component("H53_7", "Guide_tapering")
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

        # adds monitors at the same position of the previous component
        addMonitor(mycalculator, "H53_7")

        slit_A = mycalculator.add_component("slit_A", "Slit")
        slit_A.xwidth = 0.04
        slit_A.yheight = 0.12
        slit_A.set_AT([0, 0, H53_7.l + 0.3], RELATIVE=H53_7)

        # adds monitors at the same position of the previous component
        addMonitor(mycalculator, "slit_A")

        Monochromator_Arm = mycalculator.add_component("Monochromator_Arm", "Arm")
        Monochromator_Arm.set_AT([0, 0, 2], RELATIVE=slit_A)

        addMonitor(mycalculator, "mono_in")

        # double check the probability of reflection
        Monochromator = mycalculator.add_component(
            "Monochromator", "Monochromator_curved"
        )
        #   Monochromator.reflect = '"HOPG.rfl"'
        #   Monochromator.transmit = '"HOPG.trm"'
        Monochromator.gap = 0.0005
        Monochromator.NH = 13
        Monochromator.NV = 13
        Monochromator.mosaich = 30
        Monochromator.mosaicv = 30
        Monochromator.r0 = 1
        Monochromator.t0 = 0  # remove transmitted neutrons
        Monochromator.RV = "2*" + str(ThALES_L) + "*sin(DEG2RAD*a2/2)"
        Monochromator.RV = (
            "(1/( ( 1/" + str(ThALES_L) + " + 1/7.0 ) / (2*sin(DEG2RAD*a2/2))))"
        )
        Monochromator.RH = "2*" + str(ThALES_L) + "/sin(DEG2RAD*a2/2)"
        Monochromator.RH = (
            "(1/( ( 1/" + str(ThALES_L) + " + 1/2.0 ) *sin(DEG2RAD*a2/2) / 2 ))"
        )
        Monochromator.DM = monochromator_d
        Monochromator.width = 0.25
        Monochromator.height = 0.2
        Monochromator.verbose = 1
        #    Monochromator.append_EXTEND("if(flag!=SCATTERED) ABSORB;")
        Monochromator.set_AT([0, 0, 0], RELATIVE=Monochromator_Arm)
        Monochromator.set_ROTATED([0, "a2/2", 0], RELATIVE=Monochromator_Arm)

        addMonitor(mycalculator, "mono_out")

        Monochromator_Out = mycalculator.add_component("Monochromator_Out", "Arm")
        Monochromator_Out.set_AT([0, 0, 0], RELATIVE=Monochromator_Arm)
        Monochromator_Out.set_ROTATED([0, "a2", 0], RELATIVE=Monochromator_Arm)

        addMonitor(mycalculator, "mono_out_rot", [0, 0, 0.1])

        before_sample_slit = mycalculator.add_component("before_sample_slit", "Slit")
        before_sample_slit.xwidth = 0.03
        before_sample_slit.yheight = 0.028
        before_sample_slit.set_AT([0, 0, ThALES_L - 0.250], RELATIVE=Monochromator_Out)

        # addMonitor(mycalculator, "sample_slit_IN", 0.200)

        sample_mcpl_arm = mycalculator.add_component(
            "sample_mcpl_arm",
            "Arm",
            AT=[0, 0, 0],
            RELATIVE=before_sample_slit,
        )

        addMonitor(mycalculator, "sample_IN")

        # ------------------------------------------------------------
        # this new section contains the sample and the sample environment
        mycalculator, sample_mcpl_arm = self.add_new_section(
            "SampleCalc", sample_mcpl_arm, "sSAMPLE", True
        )
        # ------------------------------------------------------------

        a4 = mycalculator.add_parameter(
            "double",
            "a4",
            comment="Angle between reflected by sample and incident beams",
            unit="degree",
            value=0,
        )
        a4.add_interval(-128, 128, True)

        self._sample_environment_arm.set_AT([0, 0, 0.250], RELATIVE=sample_mcpl_arm)

        # default sample
        sample = self.set_sample_by_name("vanadium")

        Sample_Out = mycalculator.add_component("Sample_Out", "Arm")
        Sample_Out.set_AT([0, 0, 0], RELATIVE=self._sample_arm)
        Sample_Out.set_ROTATED([0, "a4", 0], RELATIVE=self._sample_arm)

        addMonitor(mycalculator, "sample_out", 0)

        after_sample_slit = mycalculator.add_component("after_sample_slit", "Slit")
        after_sample_slit.xwidth = 0.03
        after_sample_slit.yheight = 0.04
        after_sample_slit.set_AT([0, 0, 0.250], RELATIVE=Sample_Out)

        addMonitor(mycalculator, "slit_B")

        Ana_Cradle = mycalculator.add_component("Ana_Cradle", "Arm")
        Ana_Cradle.set_AT([0, 0, dist_sample_ana], RELATIVE=Sample_Out)

        # ------------------------------------------------------------
        mycalculator, Ana_Cradle = self.add_new_section(
            "AnalyzerCalc", Ana_Cradle, "sANALYZER"
        )
        # ------------------------------------------------------------

        a6 = mycalculator.add_parameter(
            "double",
            "a6",
            comment="Angle between reflected by analyzer and incident beams",
            unit="degree",
            value=33,
        )

        addMonitor(mycalculator, "analyzer_IN")

        analyzer = mycalculator.add_component("analyzer", "Monochromator_curved")
        analyzer.gap = 0.0005
        analyzer.NH = 11
        analyzer.NV = 9
        analyzer.mosaich = 30
        analyzer.mosaicv = 30
        analyzer.r0 = 0.7
        analyzer.RV = "2*" + str(dist_ana_det) + "*sin(DEG2RAD*a6/2)"
        analyzer.RH = "2*" + str(dist_ana_det) + "/sin(DEG2RAD*a6/2)"
        analyzer.DM = 3.355  # PG 002
        analyzer.width = 0.17
        analyzer.height = 0.13
        analyzer.set_AT([0, 0, 0], RELATIVE=Ana_Cradle)
        analyzer.set_ROTATED([0, "a6*0.5", 0], RELATIVE=Ana_Cradle)

        Ana_Out = mycalculator.add_component(
            "Ana_Out", "Arm", AT=[0, 0, 0], ROTATED=[0, "a6", 0], RELATIVE=Ana_Cradle
        )

        addMonitor(mycalculator, "analyzer_slit_IN")

        slit = mycalculator.add_component("slit", "Slit")
        slit.xwidth = 0.03
        slit.yheight = 0.08
        slit.set_AT([0, 0, 0.340], RELATIVE=Ana_Out)

        #    res_monitor = mycalculator.add_component("res_monitor", "Res_monitor")
        #    res_monitor.res_sample_comp = "res_sample"
        #    res_monitor.filename = '"res_monitor"'
        #    res_monitor.xwidth = 0.05
        #    res_monitor.yheight = 0.12
        #    res_monitor.set_WHEN("SAMPLE==0")
        #    res_monitor.set_AT([0, 0, dist_ana_det], RELATIVE="Ana_Out")

        detector_arm = mycalculator.add_component(
            "detector_arm", "Arm", AT=[0, 0, dist_ana_det], RELATIVE="Ana_Out"
        )
        addMonitor(mycalculator, "detector", [0, 0, 0])

        # ------------------------------------------------------------
        mycalculator, detector_arm = self.add_new_section(
            "DetectorCalc", detector_arm, "sDETECTOR"
        )
        # ------------------------------------------------------------

        detector_all = mycalculator.add_component(
            "detector_all", "Monitor", AT=[0, 0, 0.001], RELATIVE=detector_arm
        )
        detector_all.xwidth = 0.05
        detector_all.yheight = 0.12
        detector_all.restore_neutron = 0

        # detector = mycalculator.add_component("detector_final", "Monitor_nD")
        # detector.xwidth = 0.05
        # detector.yheight = 0.12
        # detector.bins = 1
        # detector.options = '"p n energy"'
        # detector.set_AT([0, 0, dist_ana_det + 0.002], RELATIVE="Ana_Out")

        # ------------------------------ instrument parameters
        myinstr.add_master_parameter("a2", {"OriginCalc": "a2"}, unit="degree")
        myinstr.add_master_parameter("a3", {"SampleCalc": "a3"}, unit="degree")
        myinstr.add_master_parameter("a4", {"SampleCalc": "a4"}, unit="degree")
        myinstr.add_master_parameter("a6", {"AnalyzerCalc": "a6"}, unit="degree")
        myinstr.master["a2"] = 79.10 * ureg.degree
        myinstr.master["a3"] = 0 * ureg.degree
        myinstr.master["a4"] = 60 * ureg.degree
        myinstr.master["a6"] = 79.10 * ureg.degree
        # ------------------------------ sample parameters
        # Do not add sample parameters. They should be modified externally retrieving
        # sample with .sample
        # this obviously will require the instrument to be recompiled


class ThALES_merge(ThALES):
    def set_sample_by_name(self, name):
        pass

    def set_sample_environment_by_name(self, name):
        pass

    def __init__(self):
        super().__init__()

        instr = self
        instr.name = "ThALES_merge"
        instr.set_sample_by_name("None")
        instr.set_sample_environment_by_name("None")

        for calcname in instr.calculators:

            calc = instr.calculators[calcname]
            calc.name = calcname + "_merge"
            calc.settings(checks=False)
            remove_comp = []
            for comp in calc.component_list[::-1]:
                if comp.component_name == "Progress_bar":
                    break
                if comp.name == "detector_arm":
                    comp.set_AT([0, 0, 0], RELATIVE="ABSOLUTE")
                    continue
                if comp.name == "detector_all":
                    continue
                calc.remove_component(comp)

            #        calc.add_component(
            #            "Origin", "Progress_bar", AT=[0, 0, 0], before="detector_arm"
            #        )

            vin = calc.add_component(
                "Vin",
                "MCPL_input",
                AT=[0, 0, 0],
                #            RELATIVE="detector_arm",
                after="Origin",
            )
            vin.filelist = "filelist"

            #
        instr.calculators[instr._calculator_name].add_parameter(
            "string", "filelist", comment="name of MCPL input file", value='"none"'
        )

        instr.add_master_parameter(
            "input_file",
            {instr._calculator_name: "filelist"},
            comment="MCPL file generated before the sample",
        )
        instr.master["input_file"] = '""'

        del instr.master["a2"]
        del instr.master["a3"]
        del instr.master["a4"]
        del instr.master["a6"]


def ThALES_from_sample():
    instr = ThALES()

    instr.name = "ThALESfromsample"
    # Removing any sample
    instr.set_sample_by_name("None")

    # removing all components before the sample arm and adding an MCPL input at the place of the sample_arm before the sample environment
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
    instr.calculators[instr._calculator_name].add_parameter(
        "string", "filelist", comment="name of MCPL input file", value='"none"'
    )

    instr.add_master_parameter(
        "input_file",
        {instr._calculator_name: "filelist"},
        comment="MCPL file generated before the sample",
    )
    instr.master["input_file"] = '""'

    del instr.master["a2"]
    return instr


# ------------------------------ Helper functions
def addMonitor(
    instr,
    name: str,
    position=[0, 0, 0],
    width=0.10,
    height=0.15,
    bins=100,
    radius=0.30,
    energy=5,
    denergy=4.9,
):
    if addMonitors is False:
        return

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
