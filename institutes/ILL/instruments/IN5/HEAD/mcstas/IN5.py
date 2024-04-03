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
# from institutes.ILL.sources.HEAD.mcstas import Full as source
# from institutes.ILL.sources.HEAD.mcstas import Gauss as source
from institutes.ILL.sources.HEAD.mcstas import Gauss as source

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
    return ["None", "nosection"]


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
        myinstr = self

        # ------------------------------------------------------------
        # Start with a first section and declaring its parameters
        mycalculator, Origin = self.add_new_section("OriginCalc")

        speed = mycalculator.add_parameter("int", "DCnu", comment="Rotation speed of disk choppers", value=8500, unit="RPM")
        speed.set_interval(0,inf,True) # must be positive
        ratio = mycalculator.add_parameter("double", "ratio", comment="", value=0.5)
        housing = mycalculator.add_parameter(
            "string", "housing", comment="", value="Fe.laz"
        )
        #coh = mycalculator.add_parameter(
        #    "string", "coh", comment="", value="Y3Fe5O12_YIG.laz"
        #)
        #inc = mycalculator.add_parameter("string", "inc", comment="", value="NULL")

        # thickness=0, height=0.025, radius=0.005, order=0)
        #================ Distances
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

        mycalculator.append_initialize('if ('+speed.name+'==0){ printf("FATAL ERROR: Chopper speed = 0 !"); exit(-1); } ')
                                 
        #----------------------------------------------------------------
        # Compute the phases of each choppers
        #-------------------------------------
        # Zero time at chopper 0
        # 1st compute the distance from the zero time position for each chopper
        # 2nd compute the phase as distance/velocity, it means, the time delay
        #----------------------------------------------------------------
        
        
        #========================================
        #   Actual sample and detector
        #========================================
        thickness    = 0.0125;
        radius    = 0.015;
        height     = 0.06;
        
        ang_ini = -11.9175;            #angular range of de detector in degrees
        ang_fin = 134.8172;            #
        det_angle = abs(ang_fin-ang_ini)/2.0 + ang_ini;




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

        def calcLtof(fcomp, lcomp, debug=True):
            complist = mycalculator.component_list[fcomp:lcomp]
            for component in complist:
                
                if component.component_name != fcomp:
                    if component.AT[0] == 0 and component.AT[1]==0:
                        L= L+component.AT[3]
                    else:
                        L= L+math.sqrt(component.AT[0]*component.AT[0]+component.AT[1]*component.AT[1]+component.AT[3]+component.AT[3])
                
                if debug:
                    print("@"+component.name+": L = "+L)
            return L

        def tofdelay(fcomp, lcomp):
            return str(calcLtof(Chopper0,Chopper1))+"/neutron_velocity"
        
        ## CHOPPER TIME-RESET##########################/
        Chopper0 = mycalculator.add_component("DC0", "DiskChopper", AT=0.23, RELATIVE="PREVIOUS")
        Chopper0.set_parameters(
            theta_0 = 20.222,
            radius = 0.285,
            yheight = 0.2,
            nu  = speed/60,
            nslit = 2,
            delay = Ch_phase[0],
            isfirst=1
        )
        
        
        """
        GERADE Guide
        """
        Guide1 = mycalculator.add_component("Guide1", "Guide_channeled", AT=L_gap, RELATIVE=Chopper0)
        Guide1.set_parameters(
            w1 = 0.03000,
            h1 = 0.20000,
            w2 = 0.03000,
            h2 = 0.17415,
            l = L_Guide1,
            R0 = Guide_Ro,
            Qcx = Guide_Qc,
            Qcy = Guide_Qc,
            alphax = Guide_alpha,
            alphay = Guide_alpha,
            mx = 1, my=2 ,
            W = Guide_W
        )


        Guide21 = mycalculator.copy_component("Guide21", Guide1, AT=L_Guide1+0.0003, RELATIVE=Guide1)
        Guide21  = mycalculator.add_component("Guide21", "Guide_channeled",
                                              AT=L_Guide1+0.0003, RELATIVE=Guide1)
        Guide21.set_parameters(
            h1 = Guide1.h2,
            h2 = 0.17000,
            l = L_Guide21,
        )
        
        Chopper1  = mycalculator.add_component("Chopper1","DiskChopper",AT=L_Guide21+disk_gap/2,RELATIVE=Guide21)
        Chopper1.set_parameters(
            theta_0 = 0.17,
            radius = 0.285,
            yheight = Ch_height[1],
            nu  = speed/60, nslit = 2, delay = tofdelay(Chopper0,Chopper1) #Ch_phase[1]
        )


        ###GUIDE TO CHOPPER2#######################
        Guide22  = mycalculator.copy_component("Guide22",Guide21,AT=L_Guide21+disk_gap,RELATIVE=Guide21)
        Guide22.set_parameters(
            h1 = Guide21.h2,
            h2 = 0.16813,
            l = L_Guide22
        )

        Chopper2  = mycalculator.add_component("Chopper2","DiskChopper",AT=L_Guide22+disk_gap/2,RELATIVE=Guide22)
        Chopper2.set_parameters(
            theta_0 = 9.0,
            radius = 0.285,
            yheight = 0.16813,
            nu   = speed/60, nslit = 2, delay = tofdelay(Chopper0,Chopper2) #Ch_phase[2]
        )
        

        #COMPONENT M1 = Monitor_nD(xwidth=0.03, yheight=0.17,
        #  options="auto time")
        #AT (0,0, disk_gap/4+0.002) RELATIVE Chopper2
        
        Guide23  = mycalculator.copy_component("Guide23", Guide22,AT=L_Guide22+disk_gap,RELATIVE=Guide22)
        Guide23.set_parameters(
            h1 = Guide22.h2,
            w2 = 0.02856,
            h2 = 0.15931,
            l = L_Guide23,
            mx = 2, my = 3
        )


        Guide3  = mycalculator.copy_component("Guide3",Guide23,AT=L_Guide23+0.0003,RELATIVE=Guide23)
        Guide3.set_parameters(
            w1 = Guide23.w2,
            h1 = Guide23.h2,
            w2 = 0.01733,
            h2 = 0.09041,
            l = L_Guide3
        )

        Guide41  = mycalculator.copy_component("Guide41",Guide3,AT=L_Guide3+0.0003),RELATIVE=Guide3)
        Guide41.set_parameters(
            w1 = Guide3.w2,
            h1 = Guide3.h2,
            w2 = 0.01579,
            h2 = 0.08100,
            l = L_Guide41
        )


        Chopper3  = mycalculator.add_component("Chopper3","DiskChopper",AT=L_Guide41+disk_gap/2,RELATIVE=Guide41)
        Chopper3.set_parameters(
            theta_0 = 9.5, radius = 0.299, yheight =  0.081,
            nu   = speed/60*ratio, nslit = 2, delay = tofdelay(Chopper0,Chopper3) #Ch_phase[3]
        )

        Guide42 = mycalculator.copy_component("Guide42", Guide41,AT=L_Guide41+disk_gap,RELATIVE=Guide41)
        Guide42.set_parameters(
            w1 = 0.01577,
            h1 = 0.08088,
            w2 = 0.01568,
            h2 = 0.08031,
            l = L_Guide42
        )

        Chopper4  = mycalculator.add_component("Chopper4","DiskChopper",AT=L_Guide42+disk_gap/2,RELATIVE=Guide42)
        Chopper4.set_parameters(
            theta_0 = 9.5, radius = 0.299, yheight = 0.08031,
            nu   = speed/60, nslit = 2, delay = tofdelay(Chopper0, Chooper4) # Ch_phase[4]
        )


        Guide43 = mycalculator.copy_component("Guide43",Guide42,AT=L_Guide42+disk_gap,RELATIVE=Guide42)
        Guide43.set_parameters(
            w1 = 0.01566,
            h1 = 0.08019,
            w2 = 0.01411,
            h2 = 0.07069,
            l = L_Guide43,
        )


        Chopper5  = mycalculator.add_component("Chopper5","DiskChopper",AT=L_Guide43+disk_gap/2,RELATIVE=Guide43)
        Chopper5.set_parameters(
            theta_0 = 3.25, radius = 0.304, yheight =  0.07069,
            nu = speed/60, nslit = 2, delay = tofdelay(Chopper0, Chopper5) #Ch_phase[5]
        )

        Guide44 = mycalculator.copy_component("Guide44",Guide43,AT=L_Guide43+disk_gap,RELATIVE=Guide43)
        Guide44.set_parameters(
            w1 = 0.01413,
            h1 = 0.07081,
            w2 = 0.01400,
            h2 = 0.0700,
            l = L_Guide44,
        )
        
        Chopper6  = mycalculator.add_component("Chopper6","DiskChopper",AT=L_Guide44+disk_gap/2,RELATIVE=Guide44)
        Chopper6.set_parameters(
            theta_0 = 3.25,
            radius = 0.304, yheight =  0.0700,
            nu = speed/60,
            nslit = 2,
            delay = tofdelay(Chopper0, Chopper6) #Ch_phase[6]
        )

        Guide45 = mycalculator.copy_component("Guide45",Guide44,AT=L_Guide44+disk_gap,RELATIVE=Guide44)
        Guide45.set_parameters(
            w1 = Guide44.w2,
            h1 = 0.06983,
            w2 = 0.01400,
            h2 = 0.05663,
            l = L_Guide45,
        )
        
        Collimator = mycalculator.copy_component("Collimator", Guide45,AT=L_Guide45+mono_gap,RELATIVE=Guide45)
        Collimator.set_parameters(
            w1 = Guide45.w2, h1 = 0.05617,
            w2 = Guide45.w2, h2 = 0.05400,
            l = L_Collimator,
        )

        Det_sample_t = mycalculator.add_component("Detector", "Monitor_nD",
                                                  AT=L_Collimator+0.0002,RELATIVE=Collimator)
        Det_sample_t.set_parameters(xwidth=0.014, yheight=0.054,
                                    options="auto t bins=20", restore_neutron=1)
        
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
        self._sample_arm.set_AT(L_CollSample+0.025-0.0002, RELATIVE=sample_mcpl_arm)
        self._sample_arm.set_ROTATE([0,det_angle,0])
        self._sample_environment_arm.set_AT(L_CollSample+0.025-0.0002, RELATIVE=sample_mcpl_arm)
        self._sample_environment.set_ROTATE([0,det_angle,0])
        
        # default sample
        #self.set_sample_focus(Lsd, 2, Lsd) # FIXME
        sample = self.set_sample_by_name("vanadium")

        Sample_Out = mycalculator.add_component(
            "Sample_Out", "Arm", AT=0, RELATIVE=self._sample_arm
        )

        arm2 = self._sample_arm
        
        #COMPONENT SAMPLE = Isotropic_Sqw(
        #  radius = radius, thickness=thickness, yheight = height,
        #  Sqw_coh=coh, Sqw_inc=inc, p_interact=0.9,
        #  order = order, d_phi = 180/PI*atan(1.5/4)*2, verbose=1)
        #AT (0,0,0) RELATIVE arm2
        #ROTATED (0,0,0) RELATIVE arm2
        #EXTEND
        #%{
        #   if(!SCATTERED) ABSORB;
        #%}

        center_det = Sample_Out
        
        #--------------- DETECTOR IDEAL ----------------------------------------
        
        detideal = mycalculator.add_component("Det_ideal_ay", "Monitor_nD", AT=[0,0,0], RELATIVE=center_det)
        detideal.set_parameters(xwidth=(4.0-0.0005-0.00002)*2, yheight=3,
                                options="banana, theta limits=[-73.36735 73.36765] bins=100, y bins=100"
                                )
        
        
        #------------ Fe HOUSING------------------------------------------------
        
        hous = mycalculator.add_component("hous", "PowderN", AT=[0,0,0],RELATIVE= center_det)
        hous.set_parameters(
            reflections=housing, radius = 4.0-0.00001, thickness = 0.0005,
        yheight = 3.0,p_transmit=0.8)
        
        
        #------------ PSD Detector ---------------------------------------------
        
        Det_PSD = mycalculator.add_component("Det_PSD", "PSD_Detector",AT=[0,0,0], RELATIVE= center_det)
        
        Det_PSD.set_parameters(
            yheight = 3.0, radius = 4.0, zdepth = 0.02600, awidth=(ang_fin-ang_ini)*PI/180*4.0,
            nx = 384, ny = 128, #type = "events",
            PressureConv = 4.75, PressureStop = 1.25, threshold=100,
            borderx=-1, bordery=-1, LensOn = 1, filename = "in5det.dat",
            FN_Conv="Gas_tables/He3inHe.table", FN_Stop="Gas_tables/He3inCF4.table")
        
        
        in5_t = mycalculator.add_component("in5_t","Monitor_nD",AT=[0,0,0], RELATIVE= center_det)
        in5_t.set_parameters(
            options="banana, t limits=[0.0206 0.0216] bins=41, parallel, previous"
        )
        
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
        

        # Ch_Ltot = L_gap+L_Guide1+0.0003+L_Guide21+disk_gap/2.0;
        # Ch_Ltot[2] = Ch_Ltot[1]+disk_gap+L_Guide22;
        # Ch_Ltot[3] = Ch_Ltot[2]+disk_gap+L_Guide23+L_Guide3+L_Guide41+2*0.0003;
        # Ch_Ltot[4] = Ch_Ltot[3]+disk_gap+L_Guide42;
        # Ch_Ltot[5] = Ch_Ltot[4]+disk_gap+L_Guide43;
        # Ch_Ltot[6] = Ch_Ltot[5]+disk_gap+L_Guide44;
        
        
        # for (i=0;i<=6;i++)
        # {
        #     Ch_phase[i]  =  Ch_Ltot[i]/neutron_velocity;
        # }
