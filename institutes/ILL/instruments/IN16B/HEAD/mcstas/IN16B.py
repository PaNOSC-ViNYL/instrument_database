"""
 IN16B instrument description.
 Authors: 
    - Markus Appel <masolomaster3000@googlemail.com>
    - Shervin Nourbakhsh <shervin86@posteo.net>

 Documentation:
  - https://www.ill.eu/users/instruments/instruments-list/in16b/description/instrument-layout

 TODO:
 - [ ] check h112_start
 - [ ] check if range in python include the last element 
"""

# ------------------------------ For McStasscript instruments
import mcstasscript as ms
from mcstasscript.interface import functions
from mcstasscript.interface import instr

# this is needed to get the location of McStas executables and libraries
my_configurator = functions.Configurator()

# ------------------------------ Importing sources
# from institutes.ILL.sources.HEAD.mcstas import Full as source
# from institutes.ILL.sources.HEAD.mcstas import Gauss as source
from institutes.ILL.sources.HEAD.mcstas import Gauss as source

# ------------------------------ Mandatory classes to use
from libpyvinyl.Instrument import Instrument
from libpyvinyl.Parameters import Parameter
from mcstas.McStasInstrumentBase import McStasInstrumentBase

# ------------------------------ Extras

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
    return ["None", "nosection"]


############## Mandatory method
def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    if flavour in [None, "None", "", "full"]:
        return IN16B()
    if flavour == "nosection":
        return IN16B(False)
    else:
        raise RuntimeError(f"Flavour {flavour} not implement")


class IN16B(McStasInstrumentBase):
    """:class: Instrument class defining the Panther instrument at ILL"""

    # ------------------------------ utility methods made available for the users

    # def wavelength_to_angle(self, value: [float, pint.Quantity]) -> pint.Quantity:
    #     """Conversion from wavelength to angle
    #     for Bragg law with lattice parameter
    #     equal to the monochromator lattice"""
    #     d_lattice = self.parameters["OriginCalc"]["monochromator_d"].pint_value
    #     return (math.asin(value / 2.0 / d_lattice) * 2 * ureg.radians).to("degrees")

    # def energy_to_angle(self, value: [float, pint.Quantity]) -> pint.Quantity:
    #     """Conversion from energy to angle
    #     for Bragg law with lattice parameter
    #     equal to the monochromator lattice"""
    #     wl = math.sqrt(81.80421 * ureg.meV / value) * ureg.angstrom
    #     return self.wavelength_to_angle(wl).to("degrees")

    # ------------------------------ Internal methods (not available to users)

    # ------------------------------ The instrument definition goes in the __init__
    def __init__(self, do_section=True):
        """Here the real definition of the instrument is performed"""

        super().__init__("IN16B", do_section)

        # ------------------------------ some local variables
        myinstr = self

        # ------------------------------------------------------------
        # Start with a first section and declaring its parameters
        mycalculator, Origin = self.add_new_section("OriginCalc")

        # ------------------------------------------------------------
        # imported source and associated parameters: check the source file!
        # - Ei
        # - dE
        mycalculator.append_initialize("dE = 0.10 * Ei;")
        HCS = source.HCS_source(mycalculator)
        HCS.E0 = "Ei"
        HCS.target_index = 2
        HCS.flux = 2.5e10
        HCS.radius = 0.100 / 2

        Ei = mycalculator.parameters["Ei"]
        Ei.value = 15 * ureg.meV
        Ei.add_interval(7.5, 130, True)
        del mycalculator.parameters["lambda"]
        del mycalculator.parameters["dlambda"]
        mycalculator.add_declare_var("double", "lambda")
        mycalculator.append_initialize("lambda = sqrt(81.80421036/Ei);")
        mycalculator.add_declare_var("double", "neutron_velocity")
        mycalculator.append_initialize("neutron_velocity = 3956.034012/lambda;")
        mycalculator.append_initialize('printf("lambda = %.2f\\n", lambda);')
        # mycalculator.append_initialize(
        #    "if(a2<=0) a2 = asin(lambda/2/mono_d)*RAD2DEG;"
        # )  # put a warning if A2 does not match

        # ------------------------------------------------------------ inpile
        R0_h112 = 0
        Qc_h112 = 0
        alpha_h112 = 0
        m_h112 = 0
        W_h112 = 0
        h112_kink = 0

        # section 1
        H112 = mycalculator.add_component(
            "H112", "Guide_gravity", AT=2.525, RELATIVE=HCS
        )
        H112.set_parameters(
            w1=0.066,
            h1=0.120,
            w2=0.06235,
            h2=0.120,
            l=1.996,
            R0=R0_h112,
            Qc=Qc_h112,
            alpha=alpha_h112,
            m=m_h112,
            W=W_h112,
        )

        # section 2
        carter_pink = mycalculator.copy_component(
            "carter_pink", H112, AT=2, RELATIVE=H112
        )
        carter_pink.set_parameters(w1=0.06235, w2=0.060, l=1.29158)

        # section A
        obt = mycalculator.copy_component("obt", H112, AT=1.292, RELATIVE=carter_pink)
        obt.set_parameters(w1=0.060, w2=0.060, l=0.2325)

        # Section 3
        guide3 = mycalculator.copy_component(
            "guide3", H112, AT=0.2325 + 0.003, RELATIVE=obt
        )
        guide3.set_parameters(w1=0.060, w2=0.075, l=6.0)

        # Section B
        # Here is the security valve with a gap of 148mm
        inpile_end = mycalculator.add_component(
            "inpile_end", "Arm", AT=6.0 + 0.148, RELATIVE=guide3
        )

        # This part begins with the vertical splitter
        # and contains only the lower guide H112A
        # negative y offset of 15.5mm
        h112_start = mycalculator.add_component(
            "h112_start", "Arm", AT=[0, -0.0155, 0], RELATIVE="PREVIOUS"
        )

        # Section 8
        # lower guide part of vertical splitter
        # straight and diverging, split in two parts
        gdiverge8a = mycalculator.copy_component(
            "gdiverge8a", "H112", AT=0.0001, RELATIVE=h112_start
        )
        gdiverge8a.set_parameters(l=1, w1=0.07537, h1=0.089, w2=0.077805, h2=0.089)

        # Section 10
        # curvature starts here, R=2000m to the right
        # guide is still diverging vertically
        # use 10 x 500mm elements, 5m in total

        gcurvediverge = mycalculator.copy_component(
            "gcurvediverge",
            H112,
            AT=4.0001,
            RELATIVE=gdiverge8a,
            ROTATED=[0, h112_kink, 0],
        )
        gcurvediverge.set_parameters(
            l=0.500,
            w1=0.090,
            h1=0.10050,
            w2=0.090,
            h2=0.10195,
        )

        def _gcurvediverge(index, h2, previous):
            gc = mycalculator.copy_component(
                f"gcurvediverge_{index}",
                gcurvediverge,
                AT=0.5001,
                RELATIVE=previous,
                ROTATED=[0, h112_kink, 0],
            )
            gc.set_parameters(h1=previous.h2, h2=h2)
            return gc

        gc = gcurvediverge
        index = 0
        for h in [
            0.10340,
            0.10485,
            0.10630,
            0.10775,
            0.10920,
            0.11065,
            0.11210,
            0.11355,
            0.11500,
        ]:
            gc = _gcurvediverge(index, h2=h, previous=gc)
            index = index + 1

        gcurvediverge_end = mycalculator.add_component(
            "gcurvediverge_end", "Arm", AT=0.5001, RELATIVE=gc
        )

        ##############################################################################
        #         // location of the second flux measurement (23m)                   #
        # COMPONENT PSD_h112_23m = Monitor_nD(options="x y",bins=PSD_bins,           #
        #                              xwidth=0.090,yheight=0.115,restore_neutron=1) #
        # AT (0, 0, 0) RELATIVE gcurvediverge_end
        # COMPONENT L_h112_23m = Monitor_nD(options="lambda per cm2,bins=300,limits=[0.5:20.5]",
        #                              xwidth=0.090,yheight=0.115,restore_neutron=1)
        # AT (0, 0, 0) RELATIVE gcurvediverge_end

        # COMPONENT cap_h112_23m = Monitor_nD(options="capture per cm2",
        #                              xwidth=0.090,yheight=0.115,restore_neutron=1)
        # AT (0, 0, 0) RELATIVE gcurvediverge_end
        #
        ##############################################################################

        # Section 11
        # curved guide with fixed cross section of 90x115 mm^2 R=2000m to the right
        # iterate 11 x 500mm elements, 5.5m in total
        # ***************************************************
        # COMPONENT gcurved11_iteration_start = Arm()
        # AT (0,0,0) RELATIVE gcurvediverge_end

        def _gcurved(i, name):
            gc = mycalculator.copy_component(
                f"{name}_{i}",
                H112,
                AT=0.5001,
                RELATIVE="PREVIOUS",
                ROTATED=[0, h112_kink, 0],
            )

            gc.set_parameters(
                l=0.500,
                w1=0.090,
                h1=0.115,
                w2=0.090,
                h2=0.115,
            )
            return gc

        for i in range(0, 11):
            gcurved11 = _gcurved(i, "gcurved11")

        gcurved11_iteration_end = mycalculator.add_component(
            "gcurved11_iteration_end",
            "Arm",
            AT=gcurved11.l + 0.0001,
            RELATIVE=gcurved11,
        )

        # Section C
        # Out of reactor building
        # 2mm Al window
        membrane = mycalculator.add_component(
            "membrane", "Al_window", AT=0.002, RELATIVE=gcurved11_iteration_end
        )
        membrane.thickness = 0.002

        # Section 12
        # same as Section 11 (copy)
        # curved guide with fixed cross section of 90x115 mm^2 R=2000m to the right
        # iterate 11 x 500mm elements, 5.5m in total
        # ***************************************************
        gcurved12_iteration_start = mycalculator.add_component(
            "gcurved12_iteration_start",
            "Arm",
            AT=0.022,
            RELATIVE=gcurved11_iteration_end,
        )

        for i in range(0, 11):
            gcurved12 = _gcurved(i, "gcurved12")

        gcurved12_iteration_end = mycalculator.add_component(
            "gcurved12_iteration_end", "Arm", AT=0.5001, RELATIVE=gcurved12
        )

        # Section D
        # OT H112 (beamline shutter)
        # gap of 109mm

        OT_H112 = mycalculator.add_component(
            "OT_H112", "Arm", AT=0.109, RELATIVE=gcurved12_iteration_end
        )

        # Section 13
        # looong curved guide part with fixed cross section of 90x115 mm^2 R=2000m to the right
        # iterate 128 x 500mm elements, 64m in total
        # ***************************************************
        for i in range(0, 98):
            gcurved13a = _gcurved(i, "gcurved13a")

        gcurved13a_end = mycalculator.add_component(
            "gcurved13a_end", "Arm", AT=0.5001, RELATIVE=gcurved13a
        )

        flux_83m = mycalculator.add_component(
            "flux83m", "Arm", AT=0, RELATIVE=gcurved13a_end
        )

        # location of the third flux measurement (83m)
        ##########################################################################################
        # COMPONENT PSD_h112_83m = Monitor_nD(options="x y",bins=PSD_bins,                       #
        #                              xwidth=0.090,yheight=0.115,restore_neutron=1)             #
        # AT (0, 0, 0) RELATIVE flux_83m                                                         #
        #                                                                                        #
        # COMPONENT L_h112_83m = Monitor_nD(options="lambda per cm2,bins=300,limits=[0.5:20.5]", #
        #                              xwidth=0.090,yheight=0.115,restore_neutron=1)             #
        # AT (0, 0, 0) RELATIVE flux_83m                                                         #
        #                                                                                        #
        # COMPONENT cap_h112_83m = Monitor_nD(options="capture per cm2",                         #
        #                              xwidth=0.090,yheight=0.115,restore_neutron=1)             #
        # AT (0, 0, 0) RELATIVE flux_83m                                                         #
        ##########################################################################################

        for i in range(0, 25):
            gcurved13b = _gcurved(i, "gcurved13b")

        gcurved13b_end = mycalculator.add_component(
            "gcurved13b_end", "Arm", AT=0.5001, RELATIVE=gcurved13b
        )

        flux_95m = mycalculator.add_component(
            "flux95m", "Arm", AT=0, RELATIVE=gcurved13b_end
        )

        for i in range(0, 5):
            gcurved13c = _gcurved(i, "gcurved13c")

        gcurved13c_end = mycalculator.add_component(
            "gcurved13c_end", "Arm", AT=0.5001, RELATIVE=gcurved13c
        )

        # Section E
        # OS IN16B (instrument shutter)
        # gap of 123mm
        OS_IN16B = mycalculator.add_component(
            "OS_IN16B", "Arm", AT=0.123, RELATIVE=gcurved13c_end
        )

        # ------------------------------
        sample_mcpl_arm = mycalculator.add_component(
            "sample_mcpl_arm",
            "Arm",
            AT=0,
            RELATIVE=before_sample_diaphragm,
        )

        # ------------------------------------------------------------
        # this new section contains the sample and the sample environment
        mycalculator, sample_mcpl_arm = self.add_new_section(
            "SampleCalc", sample_mcpl_arm, True
        )
        # ------------------------------------------------------------
        self._sample_arm.set_AT(Lcs - Lbsd, RELATIVE=sample_mcpl_arm)
        self._sample_environment_arm.set_AT(Lcs - Lbsd, RELATIVE=sample_mcpl_arm)

        # default sample
        self.set_sample_focus(Lsd, 2, Lsd)
        sample = self.set_sample_by_name("vanadium")

        Sample_Out = mycalculator.add_component(
            "Sample_Out", "Arm", AT=0, RELATIVE=self._sample_arm
        )

        # adding a shielding to avoid saving neutrons outside acceptance
        #        bs = mycalculator.add_component(
        #            "acceptance", "Beamstop", AT=[-2, 0, 0], RELATIVE=Sample_Out
        #        )
        #        bs.set_parameters(
        #            xwidth=100,
        #            yheight=100,
        #        )
        # Sample_Out.set_ROTATED([0, "a4", 0], RELATIVE=self._sample_arm)

        # ------------------------------------------------------------
        mycalculator, detector_arm = self.add_new_section("DetectorCalc", Sample_Out)
        # ------------------------------------------------------------

        # collimator = mycalculator.add_component(
        #     "collimator", "Collimator_radial", AT=[0, 0, 0], RELATIVE=detector_arm
        # )
        # collimator.set_parameters(
        #     xwidth=0,
        #     yheight=2,
        #     length=1,
        #     theta_min=5,
        #     theta_max=136,
        #     nchans=9 * 32,
        #     radius=1,
        #     nslit=1,
        #     roc=0.1,
        #     verbose=1,
        # )

        time_channels = 512
        tube_width = 0.022
        theta_bins = 9 * 32 + 8
        theta_min = -5
        angle_increment = math.asin(tube_width / 2.0 / Lsd) * 2  # * 180 / math.pi
        theta_max = theta_bins * angle_increment * 180 / math.pi + theta_min

        y_channels = 512

        theta_max = 180
        theta_min = -180
        theta_bins = 180

        detector = mycalculator.add_component(
            "PSD_TOF",
            "Monitor_nD",
            AT=[0, -0.4, 0],
            RELATIVE=detector_arm,
        )

        detector.set_parameters(
            xwidth=Lsd,  # 2.580 m
            yheight=2.0,
            # detector.restore_neutron = 1
            filename='"detector_TOF.dat"',
            options=(  # 9*32+8(spacers)
                f'" theta bins={theta_bins} limits=[{theta_min}:{theta_max}], y bins={y_channels}, time bins={time_channels} auto all list"'  # , 3He_pressure=10"'
            ),
        )

        """
        angle = 0
        angle_increment = math.asin(tube_width / 2.0 / Lsd) * 2  # * 180 / math.pi
        itube = 0

        for ibank in range(0, 9):  # 9
            for i in range(0, 32):  # 32
                # the beam is at 1.4m from the ground and the detector height is 2m,
                # so need to shift the detector 40cm lower
                # arm = mycalculator.add_component(
                #     f"arm_{itube}",
                #     "Arm",
                #     AT=[0, 0, 0],
                #     ROTATED=[0, angle, 0],
                #     RELATIVE=detector_arm,
                # )
                # print(f"angle: {angle} (+{angle_increment})")

                detector = mycalculator.add_component(
                    f"PSD_{itube}",
                    "Monitor_nD",
                    AT=[Lsd * math.sin(angle), -0.4, Lsd * math.cos(angle)],
                    RELATIVE=detector_arm,
                )

                #        image
                detector.xwidth = tube_width  # Lsd  # 2.580 m
                detector.yheight = 2.0
                # detector.restore_neutron = 1
                detector.filename = f'"detector_{itube}.dat"'
                detector.options = (
                    f'"  y bins=512, time bins=512 parallel 3He_pressure=10"'
                )
   
   
                angle = angle + angle_increment
                itube = itube + 1
            # at the end of each bank there is a Cd spacer of the size of one tube
            angle = angle + angle_increment
        # 32 tubes in 9 banks with space between banks equal to one tube
        """
        # ------------------------------ instrument parameters

        OriginCalc = myinstr.calculators["OriginCalc"]

        # myinstr.add_master_parameter(
        #     "mono_index",
        #     {OriginCalc.name: "mono_index"},
        #     unit=OriginCalc.parameters["mono_index"].unit,
        #     comment=OriginCalc.parameters["mono_index"].comment,
        # )

        myinstr.add_master_parameter(
            "energy",
            {OriginCalc.name: "Ei"},
            unit=OriginCalc.parameters["Ei"].unit,
            comment=OriginCalc.parameters["Ei"].comment,
        )


#        myinstr.add_master_parameter(
#            "a3", {"SampleCalc": "sample_rotation"}, unit="degree"
#        )
#        myinstr.add_master_parameter("a4", {"SampleCalc": "a4"}, unit="degree")
#        myinstr.add_master_parameter("a6", {"AnalyzerCalc": "a6"}, unit="degree")
#        myinstr.master["a2"] = 79.10 * ureg.degree
#        myinstr.master["a3"] = 0 * ureg.degree
#        myinstr.master["a4"] = 60 * ureg.degree
#        myinstr.master["a6"] = 79.10 * ureg.degree
# ------------------------------ sample parameters
# Do not add sample parameters. They should be modified externally retrieving
# sample with .sample
# this obviously will require the instrument to be recompiled
