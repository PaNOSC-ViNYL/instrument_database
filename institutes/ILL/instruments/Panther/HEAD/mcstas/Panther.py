"""
 Panther instrument description.
 Author: Shervin Nourbakhsh

 INFO:
  - https://www.ill.eu/users/instruments/instruments-list/panther/characteristics

 TODO:
 - [ ] sapphire xwidth and yheight unknown
 - [ ] diaphragm1 not present in the Panther schematic
 - [ ] parameters of the Monochromator
 - [ ] parameters of the Fermi Chopper
 - [ ] verify the formulat for the RV RH of the monochromator in case the parameter value is 0
 - [ ] Add the phase for the chopper
 - [ ] check what happens if chopper_rpm==0 and Efoc==0
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
    return ["None", "nosection", "quick", "quicknosection"]


def def_tests(flavour: Optional[str] = None):
    myinstrument = def_instrument(flavour)


############## Mandatory method
def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    if flavour in [None, "None", ""]:
        return Panther()
    if flavour in ["nosection", "full"]:
        return Panther(False)
    if flavour in ["quick"]:
        return Panther(True, True)
    if flavour in ["quicknosection"]:
        return Panther(False, True)
    else:
        raise RuntimeError(f"Flavour {flavour} not implement")


class Panther(McStasInstrumentBase):
    """:class: Instrument class defining the Panther instrument at ILL"""

    # ------------------------------ utility methods made available for the users

    # ------------------------------ Specialized methods required by the Instrument class [mandatory]

    # ------------------------------ Internal methods (not available to users)

    # ------------------------------ The instrument definition goes in the __init__
    def __init__(self, do_section=True, _start_from_Fermi=False):
        """Here the real definition of the instrument is performed"""

        super().__init__("Panther", do_section)

        def tofdelay(fcomp, lcomp, delay=None):
            if fcomp.name != "BC1":
                raise RuntimeError(
                    "First component should be Chopper 0 for delay calculation"
                )
            dist = 0.800
            chopper_distances = {
                "BC1": dist,
                "BC2": dist * 2,
                "BC3": dist * 3,
                "BC4": dist * 4,
                "BC5": dist * 5,
            }
            if lcomp.name not in chopper_distances:
                RuntimeError("Last component should be one of the choppers")

            L = chopper_distances[lcomp.name]
            # L = self.calcLtof(mycalculator, fcomp.name, lcomp.name, debug=True)
            if delay == None:
                delay = fcomp.delay

            print(f"{fcomp.name} -> {lcomp.name} : L = {L}; delay = {delay}")
            return str(L) + "/neutron_velocity + " + str(delay)

        # ------------------------------ some local variables
        myinstr = self

        monochromators = [
            # name, lattice parameter,
            ["PG002", 3.3550],
            ["PG004", 1.6775],
            ["Cu220", 1.2763],
            ["PG006", 1.1183],
            ["Cu331", 0.8282],
        ]

        # ------------------------------------------------------------
        # Start with a first section and declaring its parameters
        mycalculator, Origin = self.add_new_section("OriginCalc")

        mono_index = mycalculator.add_parameter(
            "int",
            "mono_index",
            comment="Monochromator index. -1 = auto: based on the energy",
            value=-1,
        )
        mono_index.add_option(-1, True)  # auto
        mono_index.add_option(list(range(0, len(monochromators))), True)

        # gnuplot> lambda(x) = sqrt(81.80421/x)
        # gnuplot> angle(x,d) = 2*asin(lambda(x)/2/d)*180/pi
        # gnuplot> p for[ d in "3.3550 1.6775 1.2763 1.1183 0.8282"]  angle(x,d) t d, 35 lt -1 t "35", 70 lt -1 t "70"

        # put in the binary the lattice parameters
        # list of only the monochromator lattices
        mycalculator.add_declare_var(
            "double",
            "mono_ds",
            array=len(monochromators),
            value=[x[1] for x in monochromators],
            comment="monochromator lattice parameter",
        )

        a2 = mycalculator.add_parameter(
            "double",
            "a2",
            comment="Angle between beam reflected by monochromator and incident beam",
            unit="degree",
            value=0,
        )
        a2.add_option(0, True)  # for automatic calculation
        a2.add_interval(36.00, 58.97, True)

        chopper_rpm = mycalculator.add_parameter(
            "double",
            "chopper_rpm",
            comment="Fermi chopper speed",
            unit="",
            value=0,
        )
        chopper_rpm.add_option(0, True)  # default value for automatic setting
        chopper_rpm.add_interval(6000, 30000, True)

        chopper_ratio = mycalculator.add_parameter(
            "double", "chopper_ratio", comment="", unit="", value=1
        )
        chopper_ratio.add_option([1, 2], True)

        if not _start_from_Fermi:
            mono_rv = mycalculator.add_parameter(
                "double",
                "mono_rv",
                comment="Monochromator vertical focusing: internally calculated if value is negative",
                value=-1,  # flat by default -> calculated internally
                unit="",
            )
            mono_rh = mycalculator.add_parameter(
                "double",
                "mono_rh",
                comment="Monochromator horizontal focusing: internally calculated if value is negative",
                value=-1,  # flat by default -> calculated internally
                unit="",
            )

        # ------------------------------------------------------------
        # imported source and associated parameters: check the source file!
        # - Ei
        # - dE
        mycalculator.append_initialize("dE = dE * Ei;")
        HCS = source.HCS_source(mycalculator)

        HCS.E0 = "Ei"
        # HCS.target_index = 2
        HCS.dist = 4.6975 + 1.4199 + 0.8 + 0.8 + 0.8
        HCS.focus_xw = 0.3
        HCS.focus_yh = 0.3

        HCS.flux = 2.5e10
        HCS.radius = 0.100 / 2

        Ei = mycalculator.parameters["Ei"]
        self.add_parameter_to_master("energy", mycalculator, Ei)

        Ei.value = 15 * ureg.meV
        Ei.add_interval(7.5, 130, True)

        # The energy resolution varies between 3.8 and 5.7% of the i of
        # the incoming energy for the pyrolytic graphite
        # monochromator, depending on the monochromator take-off
        # angle.
        mycalculator.parameters["dE"] = 0.06

        del mycalculator.parameters["lambda"]
        del mycalculator.parameters["dlambda"]
        mycalculator.add_declare_var("double", "lambda")
        mycalculator.append_initialize("lambda = sqrt(81.80421036/Ei);")
        mycalculator.add_declare_var("double", "neutron_velocity")
        mycalculator.append_initialize("neutron_velocity = 3956.034012/lambda;")
        mycalculator.append_initialize('printf("nv = %2f\\n", neutron_velocity);')
        mycalculator.append_initialize('printf("lambda = %.2f\\n", lambda);')
        # mycalculator.append_initialize("time_frame=chopper_ratio/2./chopper_rpm/60.;")
        # mycalculator.append_initialize('printf("time_frame = %.2e\\n", time_frame);')

        # optimal time focusing:
        Lcs = 0.8  # [m] chopper-sample distance
        Lsd = 2.5  # [m] sample-detector distance
        Lmc = 1.7  # [m] monochromator-chopper distance
        Lms = Lmc + Lcs  # [m] monochromator-sample distance
        Lhm = 4.123  # [m] distance between the horizontal virtual source and the monochromator

        a2_interval = a2.get_intervals()[0]
        mycalculator.add_declare_var("double", "mono_d")
        mycalculator.append_initialize('printf("%ld\\n", mono_index);')
        mycalculator.append_initialize(  # TODO: remove the cast to int (int) in the comparison. Now only needed for McStas3.3 because mono_index is long and the compiler treats it as unsigned
            "if((int)mono_index<0){\n"
            + '  printf("Selecting monochromator...");\n'
            + "  for(mono_index=0; mono_index < {nindex} && !(a2 > {amin!s} && a2 < {amax!s} ); ++mono_index)".format(
                nindex=len(monochromators),
                amin=a2_interval[0].m,
                amax=a2_interval[1].m,
            )
            + "{\n"
            + "    mono_d = mono_ds[mono_index];\n"
            + "    a2 = 2*asin(lambda/2/mono_d)*RAD2DEG;\n"
            + 'printf("mono_d = %.4f\\n",mono_d);\n'
            + 'printf("a2 = %.2f\\n",a2);\n'
            + "  }\n"
            + "}else{\n"
            + "  mono_d = mono_ds[mono_index];\n"
            + "  a2 = 2*asin(lambda/2/mono_d)*RAD2DEG;\n"
            + "}\n"
        )
        mycalculator.append_initialize(
            'printf("mono_d = %.4f\\n",mono_d);'
        )  # get the value corresponing to the selected monochromator
        mycalculator.append_initialize(
            'printf("a2 = %.2f\\n",a2);'
        )  # put a warning if A2 does not match

        mycalculator.append_initialize('printf("mono_index = %ld\\n", mono_index);')
        if not _start_from_Fermi:
            # ------------------------------------------------------------
            H12 = mycalculator.add_component("H12", "Arm", AT=[0, 0, 0], RELATIVE=HCS)

            H12_bouchon = mycalculator.add_component(
                "B", "Slit", AT=4.6975, RELATIVE=H12
            )
            H12_bouchon.set_parameters(radius=0.100)

            sapphire_arm = mycalculator.add_component(
                "sapphire_arm", "Arm", AT=0.5, RELATIVE=H12_bouchon
            )
            sapphire = mycalculator.add_component(
                "sapphire", "Filter_gen", AT=0, RELATIVE=sapphire_arm
            )
            # Thickness: 0.081 [m] (3x2.7 cm plates, c-axis along beam).
            sapphire.set_parameters(
                options='"multiply"',
                filename='"Al2O3_sapphire.trm"',
                xwidth=0.3,
                yheight=0.3,
                thickness=1.11,
            )

            # guide = self._add_chopper_guide(  index , position,  guide_template, guide_lenth)
            BC1 = mycalculator.add_component(
                "BC1",
                "DiskChopper",
                AT=1.4199,
                RELATIVE=H12_bouchon,
            )
            bc_nslits = 6
            BC1.set_parameters(
                nu="chopper_rpm/60 * 2 / chopper_ratio /" + str(bc_nslits),
                radius=0.300,
                yheight=0.175,
                nslit=bc_nslits,
                isfirst=1,
                abs_out=1,
                xwidth=0.150,
                delay=0.002,
            )

            diaphragm = mycalculator.add_component(
                "diaphragm1", "Slit", AT=0.5825, RELATIVE=BC1
            )
            diaphragm.set_parameters(xwidth=0.150, yheight=0.500)

            BC2 = mycalculator.copy_component("BC2", BC1, AT=0.800, RELATIVE=BC1)
            BC2.set_parameters(
                delay=tofdelay(BC1, BC2),
                isfirst=0,
                #            nu="-" + BC1.nu,
            )
            BC3 = mycalculator.copy_component("BC3", BC2, AT=0.800, RELATIVE=BC2)
            BC3.delay = tofdelay(BC1, BC3)
            BC4 = mycalculator.copy_component("BC4", BC2, AT=0.800, RELATIVE=BC3)
            BC4.delay = tofdelay(BC1, BC4)
            BC5 = mycalculator.copy_component("BC5", BC2, AT=0.800, RELATIVE=BC4)
            BC5.delay = tofdelay(BC1, BC5)
            BC1.verbose = 1

            L_bc5_fermi = (
                BC2.AT_data[2] + BC3.AT_data[2] + BC4.AT_data[2] + BC5.AT_data[2]
            )
            print("L_bc5_fermi before BC3 = ", BC2.AT_data[2] + BC3.AT_data[2])
            # ------------------------------

            DCH = mycalculator.add_parameter(
                "double",
                "dch",
                comment="Heavy slit horizontal aperture: if zero it is automatically calculated",
                units="m",
                value=0,
            )
            # automatically calculate dhc if not given as input parameter
            mycalculator.append_initialize(
                'if(dch==0) dch=0.0378 + 0.136*sin(DEG2RAD*a2/2);printf("dch = %.3f\\n", dch);'
            )
            diaphragm = mycalculator.add_component(  # AT=9.8185, RELATIVE=HCS
                "heavy_diaphragm", "Slit", AT=0.5011, RELATIVE=BC5
            )
            diaphragm.set_parameters(xwidth=DCH, yheight=0.150)

            # ------------------------------
            # adds monitors at the same position of the previous component

            # AT=12.1485, RELATIVE=HCS)
            Monochromator_Arm = mycalculator.add_component(
                "Monochromator_Arm", "Arm", AT=2.8311, RELATIVE=BC5
            )
            L_bc5_fermi = L_bc5_fermi + Monochromator_Arm.AT_data[2]

            Monochromator = mycalculator.add_component(
                "Monochromator", "Monochromator_curved"
            )
            Monochromator.set_parameters(  # width=0.300, height=0.220
                NH=11,
                NV=15,
                zwidth=0.019,
                yheight=0.019,
                gap=0.0005,
                t0=0,  # remove transmitted neutrons that would be absorbed by the Beamstop
                order=0,  # higher order neutrons would not pass the diaphragm
                r0=1,
            )

            #   Monochromator.reflect = '"HOPG.rfl"'
            #   Monochromator.transmit = '"HOPG.trm"'
            Monochromator.mosaic = "mono_mosaic"  # depends if PG or Cu
            #        Monochromator.r0 = 1
            Monochromator.t0 = 0  # remove transmitted neutrons
            Monochromator.RV = mono_rv
            Monochromator.RH = mono_rh
            Monochromator.DM = "mono_d"
            Monochromator.verbose = 1
            # Monochromator.append_EXTEND("if(flag!=SCATTERED) ABSORB;")
            Monochromator.set_AT([0, 0, 0], RELATIVE=Monochromator_Arm)
            Monochromator.set_ROTATED([0, "a2/2", 0], RELATIVE=Monochromator_Arm)

            mycalculator.add_declare_var(
                "double", "RMV_u", value=12.7, comment="default value for PG"
            )
            mycalculator.add_declare_var(
                "double", "RMV_w", value=0.449, comment="default value for PG"
            )
            mycalculator.add_declare_var(
                "double", "RMH_g", value=45, comment="default value for PG"
            )
            mycalculator.add_declare_var(
                "double", "RMH_h", value=221, comment="default value for PG"
            )
            mycalculator.append_initialize(  # different values for Cu
                "if(mono_index == 2 || mono_index == 4 ){ RMV_u=12.7;  RMV_w=0.449; RMH_g=45; RMH_h=221;}"
            )

            mycalculator.append_initialize(
                "if(mono_rv<0) mono_rv = (RMV_u + acos(1- RMV_w/sin(DEG2RAD*a2/2)))/1000.;"
            )
            mycalculator.append_initialize(
                # "if(mono_rh<0) mono_rh = (RMH_g + RMH_h * sin(DEG2RAD*a2/2))/1000.;"
                f"if(mono_rh<0) mono_rh = 2*(1./(1./{Lhm}+1./{Lms})  )/ sin(DEG2RAD*a2/2);"
            )

            mycalculator.add_declare_var("float", "mono_mosaic")
            mycalculator.append_initialize(
                "mono_mosaic = (mono_index==2 || mono_index==4) ? 30 : 50;"  # 0.3: 0.5
            )

            Monochromator_Out = mycalculator.add_component("Monochromator_Out", "Arm")
            Monochromator_Out.set_AT([0, 0, 0], RELATIVE=Monochromator_Arm)
            Monochromator_Out.set_ROTATED([0, "a2", 0], RELATIVE=Monochromator_Arm)

            # AT=13.8485, RELATIVE=HCS)

            beamstop = mycalculator.add_component(
                "BS", "Beamstop", AT=1, RELATIVE=Monochromator_Arm
            )
            beamstop.set_parameters(xwidth=0.5, yheight=0.5)
        else:
            Monochromator_Out = mycalculator.add_component(
                "Monochromator_Out", "Arm", AT=0, RELATIVE=Origin
            )

        slit_fermi = mycalculator.add_component(
            "slit_fermi", "Slit", AT=1.600, RELATIVE=Monochromator_Out
        )
        fermi = mycalculator.add_component(
            "fermi_chopper", "FermiChopper", AT=1.700, RELATIVE=Monochromator_Out
        )
        # L_bc5_fermi = L_bc5_fermi + fermi.AT_data[2]
        # print(L_bc5_fermi)

        fermi.set_parameters(
            radius=math.sqrt((0.0006 * 45) ** 2 + 0.023 ** 2),
            nslit=45,
            length=0.023,
            w=0.0006,
            yheight=0.113,
            nu="chopper_rpm/60",  # 1/60 to convert RPM into Hz
            verbose=2,
            eff=0.86 * 0.8,
            m=0,
            R0=0,  # no super mirror
            zero_time=2,
            # delay=0.007,
            # delay=str(L_bc5_fermi) + " / neutron_velocity",
            #            delay=0.0065,
            # phase=-23.1209 * 2,
        )
        slit_fermi.set_parameters(
            xwidth=fermi.nslit * fermi.w * 1.4, yheight=fermi.yheight * 1
        )

        if _start_from_Fermi:
            fermi.zero_time = 2
            # focusing on entry of the Fermi Chopper
            HCS.focus_yh = fermi.yheight * 0.9
            HCS.focus_xw = fermi.nslit * fermi.w * 1.2
            HCS.dist = 1.7
            # focusing on the sample area
            HCS.focus_yh = 0.04
            HCS.focus_xw = 0.02
            HCS.dist = 2

            HCS.flux = 2e9
        mycalculator.add_parameter(
            "double", "Efoc", comment="Focusing energy", unit="meV", value=0
        )
        #        mycalculator.append_initialize("if(Efoc==0) Efoc=Ei;")
        mycalculator.append_initialize(
            "if(chopper_rpm==0) chopper_rpm = 60*neutron_velocity * fabs(tan(DEG2RAD*a2 / 2)) / PI / (({Lcs} + {Lsd}*pow((1 - Efoc / Ei),(-3/2))) * (1 - {Lms} / {Lhm}));".format(
                Lcs=Lcs, Lsd=Lsd, Lms=Lms, Lhm=Lhm
            )
        )

        mycalculator.append_initialize(
            'printf("Chopper_rpm = %.2f\\nNeutron velocity: %.2f\\nE_focus: %.2f\\n", chopper_rpm, neutron_velocity, Efoc );'
        )

        # ------------------------------
        bsw = mycalculator.add_parameter(
            "double",
            "before_sample_slit_width",
            comment="Horizontal width of the slit opening",
            unit="m",
            value=0.20,
        )
        bsh = mycalculator.add_parameter(
            "double",
            "before_sample_slit_height",
            comment="Vertical height of the slit opening",
            unit="m",
            value=0.20,
        )

        Lbsd = fermi.radius + 0.04

        monitor = mycalculator.add_component(
            "monitor", "Monitor_nD", AT=Lbsd - 0.01, RELATIVE=fermi
        )
        monitor.set_parameters(
            xwidth=0.01,
            yheight=0.01,
            zdepth=0.01,
            filename='"monitor.dat"',
            restore_neutron=1,
            options='"box intensity, bins=1 pressure=0.001"',
            # options='"x bins={} y bins={} file={}"'.format(1, 1, "counter.dat"
        )

        before_sample_diaphragm = mycalculator.add_component(
            "before_sample_diaphragm", "Slit", AT=Lbsd, RELATIVE=fermi
        )
        before_sample_diaphragm.set_parameters(xwidth=bsw, yheight=bsh)

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
        self.sample_focus(Lsd, 2, Lsd)
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
        """
        the is the possibility to simulate the transmission exactly or by approximation depending if nblades is given or not
        """
        detector_arm = mycalculator.add_component(
            "detector_arm", "Arm", AT=[0, -0.4, 0], RELATIVE=Sample_Out
        )

        collimator = mycalculator.add_component(
            "collimator", "Collimator_radial", AT=[0, 0, 0], RELATIVE=detector_arm
        )
        collimator.set_parameters(  # as reported in https://doi.org/10.1051/epjconf/202227202001
            radius=0.434,
            length=0.2,
            nchan=9 * 32 + 8,
            xwidth=0,
            yheight=0.587,  # 2,
            theta_min=-5,
            theta_max=136,
            nslit=1,
            roc=0.02,
            verbose=1,
        )

        # ------------------------------------------------------------
        mycalculator, detector_arm = self.add_new_section("DetectorCalc", detector_arm)

        tube_width = 0.022
        theta_bins = 9 * 32 + 8
        theta_min = -5
        angle_increment = math.asin(tube_width / 2.0 / Lsd) * 2  # * 180 / math.pi
        theta_max = theta_bins * angle_increment * 180 / math.pi + theta_min

        ny = mycalculator.add_parameter(
            "int",
            "ny",
            value=256,
            comment="Number of channels along the tube (y direction)",
        )

        nt = mycalculator.add_parameter(
            "int", "nt", value=512, comment="Number of time channels"
        )
        # if chopper_ratio is 2: nt=1024
        detector = mycalculator.add_component(
            "detector",
            "Cyl_TOF",
            AT=0,
            RELATIVE=detector_arm,
        )

        chopper_rpm = mycalculator.add_parameter(
            "double",
            "chopper_rpm",
            comment="Fermi chopper speed",
            unit="",
            value=0,
        )

        chopper_ratio = mycalculator.add_parameter(
            "double", "chopper_ratio", comment="", unit="", value=1
        )

        time_frame = mycalculator.add_declare_var("double", "time_frame")

        mycalculator.append_initialize(
            "time_frame = chopper_ratio/2.0/chopper_rpm/60.0;"
        )
        detector.set_parameters(
            ny=ny,
            nt=nt,
            yheight=2.0,
            radius=Lsd,
            phimin=-17.328,
            # phimax=136,
            tmin=1e-3,
            tmax=time_frame,
            nphi_groups=9,
            nphi_pergroup=32,
            phi_groupgap=0.518,
            phi_binwidth=0.518,
            # saveingap= 0|1
        )

        # ------------------------------ instrument parameters

    def set_test(self, test_number: Optional[int] = None):
        myinstrument = self
        myinstrument.sample_shape("holder")
        if test_number == 0:  # flux at sample position
            myinstrument.sample_holder(None, None)
            myinstrument.set_sample_by_name("monitor")
            myinstrument.master["energy"] = 19 * ureg.meV
            myinstrument.sample.xwidth = 0.02
            myinstrument.sample.yheight = 0.04

        elif test_number == 1:  # direct beam with empty sample holder
            myinstrument.set_sample_by_name("None")
        elif test_number == 2:  # with sample
            myinstrument.set_sample_by_name("qSq")
            myinstrument.master[
                "sqw_file"
            ] = '"./institutes/ILL/instruments/D11/HEAD/mcstas/data/simul_5711.sq"'
        elif test_number == -1:  # direct beam no beamstop
            myinstrument.set_sample_by_name("None")
            myinstrument.sample_holder(None, None)
            myinstrument.master["attenuator_index"] = 6
            myinstrument.master["bs_index"] = -1
        else:
            raise RuntimeError(f"Test number {test_number} out of range")

    def test_datafile(self, test_number: Optional[int] = None):
        file = ""
        if test_number == 0 or test_number == -1:  # direct attenuated beam
            file = "institutes/ILL/instruments/D11/HEAD/mcstas/data/005708.nxs"
        elif test_number == 1:  # direct beam with empty sample holder
            file = "institutes/ILL/instruments/D11/HEAD/mcstas/data/005721.nxs"
        elif test_number >= 2:
            file = "institutes/ILL/instruments/D11/HEAD/mcstas/data/005711.nxs"
        else:
            raise RuntimeError(f"Test number {test_number} out of range")
        return file


# ------------------------------ sample parameters
# Do not add sample parameters. They should be modified externally retrieving
# sample with .sample
# this obviously will require the instrument to be recompiled
