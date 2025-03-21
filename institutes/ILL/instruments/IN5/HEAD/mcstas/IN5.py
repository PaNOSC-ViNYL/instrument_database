"""
 IN5 instrument description.
 Author: Shervin Nourbakhsh based on the mcstas example IN5 written by
         E. Farhi, J. Ollivier, Celia Castan Guerrero
        
 TODO:
"""

# ------------------------------ For McStasscript instruments
import mcstasscript as ms
from mcstasscript.interface import functions
from mcstasscript.interface import instr

# this is needed to get the location of McStas executables and libraries
my_configurator = functions.Configurator()

# ------------------------------ Importing sources
from institutes.ILL.sources.HEAD.mcstas import Full as source

# from institutes.ILL.sources.HEAD.mcstas import Gauss as source
# from institutes.ILL.sources.HEAD.mcstas import Gauss as source

# from institutes.ILL.samples.vanadium import set_vanadium_sample

# ------------------------------ Mandatory classes to use
from libpyvinyl.Instrument import Instrument
from libpyvinyl.Parameters import Parameter
from mcstas.McStasInstrumentBase import McStasInstrumentBase

# ------------------------------ Extras
# import os  # to add the path of custom mcstas components

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
    return ["None", "full", "nosection"]


############## Mandatory method
def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    if flavour in [None, "None", "", "full"]:
        return IN5()
    if flavour == "nosection":
        return IN5(False)
    else:
        raise RuntimeError(f"Flavour {flavour} not implement")


class IN5(McStasInstrumentBase):
    """:class: Instrument class defining the IN5 instrument at ILL"""

    # ------------------------------ utility methods made available for the users

    # ------------------------------ Internal methods (not available to users)

    # ------------------------------ The instrument definition goes in the __init__
    def __init__(self, do_section=True):
        """Here the real definition of the instrument is performed"""

        super().__init__("IN5", do_section)

        # ------------------------------ some local variables

        # ------------------------------------------------------------
        # Start with a first section and declaring its parameters
        mycalculator, Origin = self.add_new_section("OriginCalc")

        nu = mycalculator.add_parameter(
            "double",  # "int"
            "speed",
            comment="Rotation frequency of disk choppers RPM",
            value=8500,
            unit="",
        )
        nu.add_interval(2000, 17000, True)  # must be positive

        ratio = mycalculator.add_parameter("double", "ratio", comment="", value=0.5)

        # coh = mycalculator.add_parameter(
        #    "string", "coh", comment="", value="Y3Fe5O12_YIG.laz"
        # )
        # inc = mycalculator.add_parameter("string", "inc", comment="", value="NULL")

        # thickness=0, height=0.025, radius=0.005, order=0)
        # ================ Distances
        L_gap = 0.2130  # gap VTE+OT-H16
        L_Guide1 = 4.3900  # for gerade Guide1
        L_Guide21 = 0.6950  # for gerade Guide21
        L_Guide22 = 0.1300  # for gerade Guide22
        L_Guide23 = 0.69500  # for gerade Guide23
        disk_gap = 0.02  # full gap at choppers
        L_Guide3 = 5.5125  # for gerade Guide3
        L_Guide41 = 0.7425  # for gerade Guide41
        L_Guide42 = 0.0350  # for gerade Guide42
        L_Guide43 = 0.7500  # for gerade Guide43
        L_Guide44 = 0.0350  # for gerade Guide44
        L_Guide45 = 0.7900  # for gerade Guide45
        mono_gap = 0.0300  # gap for the 1st monitor
        L_Collimator = 0.1300  # for gerade Collimator
        L_CollSample = 0.2400 - 0.025  # the sample chamber size & keep

        # Alt Guide coating
        alt_Guide_Qc = 0.021745  # for m=1 alpha and W aren't used.
        alt_Guide_Ro = 0.995
        alt_Guide_alpha = 6.07
        alt_Guide_W = 0.0023
        # New Guide and super-mirors
        Guide_Qc = 0.02275
        Guide_Ro = 0.996
        Guide_alpha = 5.75
        Guide_W = 0.00125

        # ==========================================================================
        #                 Choppers
        # ==========================================================================

        mycalculator.append_initialize(
            "if ("
            + nu.name
            + '==0){ printf("FATAL ERROR: Chopper nu = 0 !"); exit(-1); } '
        )

        # ----------------------------------------------------------------
        # Compute the phases of each choppers
        # -------------------------------------
        # Zero time at chopper 0
        # 1st compute the distance from the zero time position for each chopper
        # 2nd compute the phase as distance/velocity, it means, the time delay
        # ----------------------------------------------------------------

        # ========================================
        #   Actual sample and detector
        # ========================================
        # thickness = 0.0125
        # radius = 0.015
        # height = 0.06

        ang_ini = -11.9175
        # angular range of de detector in degrees
        ang_fin = 134.8172
        #
        det_angle = abs(ang_fin - ang_ini) / 2.0 + ang_ini

        mysource = source.VCS_source(mycalculator)

        # mysource.target_index = 2
        #        HCS.flux = 2.5e10

        lambda0 = mycalculator.parameters["lambda"]
        lambda0.value = 4.5

        # Ei = mycalculator.parameters["Ei"]
        # Ei.value = 15 * ureg.meV
        # Ei.add_interval(7.5, 130, True)
        # del mycalculator.parameters["Ei"]
        # del mycalculator.parameters["dE"]
        # del mycalculator.parameters["lambda"]
        # del mycalculator.parameters["dlambda"]
        # mycalculator.add_declare_var("double", "lambda")
        # mycalculator.append_initialize("lambda = sqrt(81.80421036/Ei);")

        mycalculator.add_declare_var("double", "neutron_velocity")
        mycalculator.append_initialize("neutron_velocity = 3956.034012/lambda;")
        mycalculator.append_initialize('printf("nv = %2f\\n", neutron_velocity);')
        mycalculator.append_initialize('printf("lambda = %.2f\\n", lambda);')
        mycalculator.append_initialize(
            'dlambda = 0.1*lambda;printf("dlambda = %.2f\\n", dlambda);'
        )

        def tofdelay(fcomp, lcomp, delay=0):
            if fcomp.name != "Chopper0":
                raise RuntimeError(
                    "First component should be Chopper 0 for detay calculation"
                )
            if lcomp.name == "Chopper1":
                L = 5.3083
            elif lcomp.name == "Chopper2":
                L = 5.4583
            elif lcomp.name == "Chopper3":
                L = 12.4289
            elif lcomp.name == "Chopper4":
                L = 12.4839
            elif lcomp.name == "Chopper5":
                L = 13.2539
            elif lcomp.name == "Chopper6":
                L = 13.3089
            else:
                RuntimeError("Last component should be one of the choppers")

            # L = self.calcLtof(mycalculator, fcomp.name, lcomp.name, debug=True)
            print(f"{fcomp.name} -> {lcomp.name} : L = {L}")
            return str(L) + "/neutron_velocity + " + str(delay)

        SourceTarget = mycalculator.add_component(
            "sourcetarget", "Arm", AT=2.55, RELATIVE="PREVIOUS"
        )
        mysource.dist = SourceTarget.AT_data[2]

        ## CHOPPER TIME-RESET##########################/
        Chopper0 = mycalculator.add_component(
            "Chopper0", "DiskChopper", AT=0.2, RELATIVE="PREVIOUS"
        )
        Chopper0.set_parameters(
            theta_0=20.222,
            radius=0.285,
            yheight=0.2,
            nu="(speed/60)",
            nslit=2,
            delay=0.01,
            isfirst=1,
        )

        Guide1 = mycalculator.add_component(
            "Guide1", "Guide_channeled", AT=L_gap, RELATIVE=Chopper0
        )
        Guide1.set_parameters(
            w1=0.03000,
            h1=0.20000,
            w2=0.03000,
            h2=0.17415,
            l=L_Guide1,
            R0=Guide_Ro,
            Qcx=Guide_Qc,
            Qcy=Guide_Qc,
            alphax=Guide_alpha,
            alphay=Guide_alpha,
            mx=1,
            my=2,
            W=Guide_W,
        )

        Guide21 = mycalculator.copy_component(
            "Guide21", Guide1, AT=L_Guide1 + 0.0003, RELATIVE=Guide1
        )
        Guide21.set_parameters(
            h1=Guide1.h2,
            h2=0.17000,
            l=L_Guide21,
        )

        # P1
        Chopper1 = mycalculator.add_component(
            "Chopper1", "DiskChopper", AT=L_Guide21 + disk_gap / 2, RELATIVE=Guide21
        )
        Chopper1.set_parameters(
            theta_0=0.17,
            radius=0.285,
            yheight=0.17,
            nu="(speed/60)",
            nslit=2,
            delay=tofdelay(Chopper0, Chopper1, Chopper0.delay),  # Ch_phase[1]
        )

        ###GUIDE TO CHOPPER2#######################
        Guide22 = mycalculator.copy_component(
            "Guide22", Guide21, AT=L_Guide21 + disk_gap, RELATIVE=Guide21
        )
        Guide22.set_parameters(h1=Guide21.h2, h2=0.16813, l=L_Guide22)

        # P2
        Chopper2 = mycalculator.add_component(
            "Chopper2", "DiskChopper", AT=L_Guide22 + disk_gap / 2, RELATIVE=Guide22
        )
        Chopper2.set_parameters(
            theta_0=9.0,
            radius=0.285,
            yheight=0.16813,
            nu="(speed/60)",
            nslit=2,
            delay=tofdelay(Chopper0, Chopper2, Chopper0.delay),  # Ch_phase[2]
        )

        # COMPONENT M1 = Monitor_nD(xwidth=0.03, yheight=0.17,
        #  options="auto time")
        # AT (0,0, disk_gap/4+0.002) RELATIVE Chopper2

        Guide23 = mycalculator.copy_component(
            "Guide23", Guide22, AT=L_Guide22 + disk_gap, RELATIVE=Guide22
        )
        Guide23.set_parameters(
            h1=Guide22.h2, w2=0.02856, h2=0.15931, l=L_Guide23, mx=2, my=3
        )

        Guide3 = mycalculator.copy_component(
            "Guide3", Guide23, AT=L_Guide23 + 0.0003, RELATIVE=Guide23
        )
        Guide3.set_parameters(
            w1=Guide23.w2, h1=Guide23.h2, w2=0.01733, h2=0.09041, l=L_Guide3
        )

        Guide41 = mycalculator.copy_component(
            "Guide41", Guide3, AT=L_Guide3 + 0.0003, RELATIVE=Guide3
        )
        Guide41.set_parameters(
            w1=Guide3.w2, h1=Guide3.h2, w2=0.01579, h2=0.08100, l=L_Guide41
        )

        Chopper3 = mycalculator.add_component(
            "Chopper3", "DiskChopper", AT=L_Guide41 + disk_gap / 2, RELATIVE=Guide41
        )
        Chopper3.set_parameters(
            theta_0=9.5,
            radius=0.299,
            yheight=0.081,
            nu="(speed/60.0 * ratio)",
            nslit=2,
            delay=tofdelay(Chopper0, Chopper3, Chopper0.delay),  # Ch_phase[3]
        )

        Guide42 = mycalculator.copy_component(
            "Guide42", Guide41, AT=L_Guide41 + disk_gap, RELATIVE=Guide41
        )
        Guide42.set_parameters(
            w1=0.01577, h1=0.08088, w2=0.01568, h2=0.08031, l=L_Guide42
        )

        Chopper4 = mycalculator.add_component(
            "Chopper4", "DiskChopper", AT=L_Guide42 + disk_gap / 2, RELATIVE=Guide42
        )
        Chopper4.set_parameters(
            theta_0=9.5,
            radius=0.299,
            yheight=0.08031,
            nu="(speed/60.0)",
            nslit=2,
            delay=tofdelay(Chopper0, Chopper4, Chopper0.delay),  # Ch_phase[4]
        )

        Guide43 = mycalculator.copy_component(
            "Guide43", Guide42, AT=L_Guide42 + disk_gap, RELATIVE=Guide42
        )
        Guide43.set_parameters(
            w1=0.01566,
            h1=0.08019,
            w2=0.01411,
            h2=0.07069,
            l=L_Guide43,
        )

        Chopper5 = mycalculator.add_component(
            "Chopper5", "DiskChopper", AT=L_Guide43 + disk_gap / 2, RELATIVE=Guide43
        )
        Chopper5.set_parameters(
            theta_0=3.25,
            radius=0.304,
            yheight=0.07069,
            nu="(speed/60.0)",
            nslit=2,
            delay=tofdelay(Chopper0, Chopper5, Chopper0.delay),  # Ch_phase[5]
        )

        Guide44 = mycalculator.copy_component(
            "Guide44", Guide43, AT=L_Guide43 + disk_gap, RELATIVE=Guide43
        )
        Guide44.set_parameters(
            w1=0.01413,
            h1=0.07081,
            w2=0.01400,
            h2=0.0700,
            l=L_Guide44,
        )

        Chopper6 = mycalculator.add_component(
            "Chopper6", "DiskChopper", AT=L_Guide44 + disk_gap / 2, RELATIVE=Guide44
        )
        Chopper6.set_parameters(
            theta_0=3.25,
            radius=0.304,
            yheight=0.0700,
            nu="(speed/60.0)",
            nslit=2,
            delay=tofdelay(Chopper0, Chopper6, Chopper0.delay),  # Ch_phase[6]
        )

        Guide45 = mycalculator.copy_component(
            "Guide45", Guide44, AT=L_Guide44 + disk_gap, RELATIVE=Guide44
        )
        Guide45.set_parameters(
            w1=Guide44.w2,
            h1=0.06983,
            w2=0.01400,
            h2=0.05663,
            l=L_Guide45,
        )

        Collimator = mycalculator.copy_component(
            "Collimator", Guide45, AT=L_Guide45 + mono_gap, RELATIVE=Guide45
        )
        Collimator.set_parameters(
            w1=Guide45.w2,
            h1=0.05617,
            w2=Guide45.w2,
            h2=0.05400,
            l=L_Collimator,
        )

        Det_sample_t = mycalculator.add_component(
            "Detector", "Monitor_nD", AT=L_Collimator + 0.0002, RELATIVE=Collimator
        )
        Det_sample_t.set_parameters(
            xwidth=0.014, yheight=0.054, options='"auto t bins=20"', restore_neutron=1
        )

        # ------------------------------
        sample_mcpl_arm = mycalculator.add_component(
            "sample_mcpl_arm",
            "Arm",
            AT=0,
            RELATIVE=Det_sample_t,
        )

        # ------------------------------------------------------------
        # this new section contains the sample and the sample environment
        mycalculator, sample_mcpl_arm = self.add_new_section(
            "SampleCalc", sample_mcpl_arm, True
        )
        # ------------------------------------------------------------
        self._sample_arm.set_AT(L_CollSample + 0.025 - 0.0002, RELATIVE=sample_mcpl_arm)
        self._sample_arm.set_ROTATED([0, det_angle, 0])
        self._sample_environment_arm.set_AT(
            L_CollSample + 0.025 - 0.0002, RELATIVE=sample_mcpl_arm
        )
        # self._sample_environment.set_ROTATED([0, det_angle, 0])

        # default sample
        self.set_sample_focus(8, 3, 8)  # FIXME
        sample = self.set_sample_by_name("vanadium")

        Sample_Out = mycalculator.add_component(
            "Sample_Out", "Arm", AT=0, RELATIVE=self._sample_arm
        )

        # arm2 = self._sample_arm

        # COMPONENT SAMPLE = Isotropic_Sqw(
        #  radius = radius, thickness=thickness, yheight = height,
        #  Sqw_coh=coh, Sqw_inc=inc, p_interact=0.9,
        #  order = order, d_phi = 180/PI*atan(1.5/4)*2, verbose=1)
        # AT (0,0,0) RELATIVE arm2
        # ROTATED (0,0,0) RELATIVE arm2
        # EXTEND
        #%{
        #   if(!SCATTERED) ABSORB;
        #%}

        mycalculator, center_det = self.add_new_section("DetectorCalc", Sample_Out)
        nt = mycalculator.add_parameter(
            "double", "nt", comment="Number of time channels", value=512  # int
        )
        ny = mycalculator.add_parameter(
            "double",
            "ny",
            comment="Number of vertical position channels",
            value=256,  # "int"
        )
        epchannel = mycalculator.add_parameter(
            "int",
            "epchannel",
            comment="Elastic peak position in number of channels",
            value=295,
        )
        housing = mycalculator.add_parameter(
            "string", "housing", comment="", value='"Fe.laz"'
        )

        # --------------- DETECTOR IDEAL ----------------------------------------

        detideal = mycalculator.add_component(
            "Det_ideal_ay", "Monitor_nD", AT=[0, 0, 0], RELATIVE=center_det
        )
        detideal.set_parameters(
            xwidth=(4.0 - 0.0005 - 0.00002) * 2,
            yheight=3,
            options='"banana, theta limits=[-73.36735 73.36765] bins=100, y bins=100"',
        )

        # ------------ Fe HOUSING------------------------------------------------

        hous = mycalculator.add_component(
            "hous", "PowderN", AT=[0, 0, 0], RELATIVE=center_det
        )
        hous.set_parameters(
            reflections=housing,
            radius=4.0 - 0.00001,
            thickness=0.0005,
            yheight=3.0,
            p_transmit=0.8,
        )

        # ------------ PSD Detector ---------------------------------------------
        detector = mycalculator.add_component(
            "detector", "Cyl_TOF", AT=[0, 0, 0], RELATIVE=center_det
        )
        detector.set_parameters(
            nphi=384,
            ny=ny,
            nt=nt,
            yheight=3.0,
            radius=4.0,
            phimin=-11,
            phimax=134,
            tmin=0,
            tmax=1,
        )

        Det_PSD = mycalculator.add_component(
            "Det_PSD", "PSD_Detector", AT=[0, 0, 0], RELATIVE=center_det
        )

        Det_PSD.set_parameters(
            yheight=3.0,
            radius=4.0,
            zdepth=0.02600,
            awidth=(ang_fin - ang_ini) * math.pi / 180 * 4.0,
            nx=384,
            ny=128,  # type = "events",
            PressureConv=4.75,
            PressureStop=1.25,
            threshold=100,
            borderx=-1,
            bordery=-1,
            LensOn=1,
            filename='"in5det.dat"',
            FN_Conv='"Gas_tables/He3inHe.table"',
            FN_Stop='"Gas_tables/He3inCF4.table"',
        )

        in5_t = mycalculator.add_component(
            "in5_t", "Monitor_nD", AT=[0, 0, 0], RELATIVE=center_det
        )
        in5_t.set_parameters(
            options='"banana, t limits=[0.0206 0.0216] bins=41, parallel, previous"'
        )

        # ------------------------------ instrument parameters

        OriginCalc = self.calculators["OriginCalc"]
        DetectorCalc = None
        if do_section:
            DetectorCalc = self.calculators["DetectorCalc"]
        else:
            DetectorCalc = OriginCalc

        self.add_master_parameter(
            "speed",
            {OriginCalc.name: "speed"},
            unit=OriginCalc.parameters["speed"].unit,
            comment=OriginCalc.parameters["speed"].comment,
        )

        self.add_master_parameter(
            "ratio",
            {OriginCalc.name: "ratio"},
            unit=OriginCalc.parameters["ratio"].unit,
            comment=OriginCalc.parameters["ratio"].comment,
        )

        self.add_master_parameter(
            "lambda",
            {OriginCalc.name: "lambda"},
            unit=OriginCalc.parameters["lambda"].unit,
            comment=OriginCalc.parameters["lambda"].comment,
        )

        self.add_master_parameter(
            "nt",
            {DetectorCalc.name: "nt"},
            unit=DetectorCalc.parameters["nt"].unit,
            comment=DetectorCalc.parameters["nt"].comment,
        )

        self.add_master_parameter(
            "epchannel",
            {DetectorCalc.name: "epchannel"},
            unit=DetectorCalc.parameters["epchannel"].unit,
            comment=DetectorCalc.parameters["epchannel"].comment,
        )

        for c in [
            "Chopper0",
            "Chopper1",
            "Chopper2",
            "Chopper3",
            "Chopper4",
            "Chopper5",
            "Chopper6",
        ]:
            print(c + " L=" + str(self.calcLtof(mycalculator, "Chopper0", c)))

        self.master["ratio"] = 0.5
        self.master["speed"] = 8500
        self.master["lambda"] = 4.5
        self.master["nt"] = 512
        self.master["epchannel"] = 295
        #        myinstr.add_master_parameter("a4", {"SampleCalc": "a4"}, unit="degree")
        #        myinstr.add_master_parameter("a6", {"AnalyzerCalc": "a6"}, unit="degree")
        #        myinstr.master["a2"] = 79.10 * ureg.degree
        #        myinstr.master["a3"] = 0 * ureg.degree
        #        myinstr.master["a4"] = 60 * ureg.degree
        #        myinstr.master["a6"] = 79.10 * ureg.degree
        # ------------------------------ sample parameters
