"""
 Panther instrument description.
 Author: Shervin Nourbakhsh

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
    return ["None", "nosection"]


############## Mandatory method
def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    if flavour in [None, "None", "", "full"]:
        return Panther()
    if flavour == "nosection":
        return Panther(False)
    else:
        raise RuntimeError(f"Flavour {flavour} not implement")


class Panther(McStasInstrumentBase):
    """:class: Instrument class defining the Panther instrument at ILL"""

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
        wl = math.sqrt(81.80421 * ureg.meV / value) * ureg.angstrom
        return self.wavelength_to_angle(wl).to("degrees")

    # ------------------------------ Specialized methods required by the Instrument class [mandatory]

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

    # ------------------------------ The instrument definition goes in the __init__
    def __init__(self, do_section=True):
        """Here the real definition of the instrument is performed"""

        super().__init__("Panther", do_section)

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
        mono_lattices = []
        for d in monochromators:
            mono_lattices.append(d[1])

        # ------------------------------------------------------------
        # Start with a first section and declaring its parameters
        mycalculator, Origin = self.add_new_section("OriginCalc")

        mono_index = mycalculator.add_parameter(
            "int",
            "mono_index",
            comment="Monochromator index, automatically based on the energy if set to -1",
            value=-1,
        )
        mono_index.add_interval(-1, len(monochromators), True)

        # gnuplot> lambda(x) = sqrt(81.80421/x)
        # gnuplot> angle(x,d) = 2*asin(lambda(x)/2/d)*180/pi
        # gnuplot> p for[ d in "3.3550 1.6775 1.2763 1.1183 0.8282"]  angle(x,d) t d, 35 lt -1 t "35", 70 lt -1 t "70"

        # put in the binary the lattice parameters
        mycalculator.add_declare_var(
            "double",
            "mono_ds",
            array=len(monochromators),
            value=mono_lattices,
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

        mono_rv = mycalculator.add_parameter(
            "double",
            "mono_rv",
            comment="Monochromator vertical focusing: internally calculated if value is zero",
            value=-1,  # flat by default -> calculated internally
            unit="",
        )
        mono_rh = mycalculator.add_parameter(
            "double",
            "mono_rh",
            comment="Monochromator horizontal focusing: internally calculated if value is zero",
            value=-1,  # flat by default -> calculated internally
            unit="",
        )

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
        chopper_ratio.add_interval(1, 2, True)

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
        a2_interval = a2.get_intervals()[0]
        mycalculator.add_declare_var("double", "mono_d")
        mycalculator.append_initialize('printf("%d\\n", mono_index);')
        mycalculator.append_initialize( # TODO: remove the cast to int (int) in the comparison. Now only needed for McStas3.3 because mono_index is long and the compiler treats it as unsigned
            "if((int)mono_index<0){\n"
            + '  printf("Selecting monochromator...");\n'
            + "  for(mono_index=0; mono_index < "
            + str(len(monochromators))
            + " && !(a2 > "
            + str(a2_interval[0].m)
            + " && a2<"
            + str(a2_interval[1].m)
            + "); ++mono_index){\n"
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

        mycalculator.append_initialize('printf("mono_index = %d\\n", mono_index);')

        # ------------------------------------------------------------
        H12 = mycalculator.add_component("H12", "Arm", AT=[0, 0, 0], RELATIVE=HCS)

        H12_bouchon = mycalculator.add_component("B", "Slit", AT=4.6975, RELATIVE=H12)
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
            phase=-41.3346 / 2,
        )

        diaphragm = mycalculator.add_component(
            "diaphragm1", "Slit", AT=0.5825, RELATIVE=BC1
        )
        diaphragm.set_parameters(xwidth=0.150, yheight=0.500)

        BC2 = mycalculator.copy_component("BC2", BC1, AT=0.800, RELATIVE=BC1)
        tdelay = "0.800 / neutron_velocity"
        BC2.delay = tdelay
        BC2.isfirst = 0
        # BC2.phase = 41.3346
        BC2.nu = "-" + BC1.nu
        BC3 = mycalculator.copy_component("BC3", BC2, AT=0.800, RELATIVE=BC2)
        BC3.nu = BC1.nu
        BC3.delay = BC2.delay + " + " + tdelay
        BC4 = mycalculator.copy_component("BC4", BC2, AT=0.800, RELATIVE=BC3)
        BC4.delay = BC3.delay + " + " + tdelay
        BC5 = mycalculator.copy_component("BC5", BC2, AT=0.800, RELATIVE=BC4)
        BC5.delay = BC4.delay + " + " + tdelay
        BC1.verbose = 1

        L_bc5_fermi = BC2.AT_data[2] + BC3.AT_data[2] + BC4.AT_data[2] + BC5.AT_data[2]
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
            'if(dch==0) dch=(37.8 + 136.0*sin(DEG2RAD*a2/2))/1000;printf("dch = %.2f\\n", dch);'
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
        )

        #   Monochromator.reflect = '"HOPG.rfl"'
        #   Monochromator.transmit = '"HOPG.trm"'
        Monochromator.gap = 0.0005
        Monochromator.mosaic = "mono_mosaic"  # depends if PG or Cu
        Monochromator.r0 = 1
        Monochromator.t0 = 1  # remove transmitted neutrons
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
        # optimal time focusing:
        Lcs = 0.8  # [m] chopper-sample distance
        Lsd = 2.5  # [m] sample-detector distance
        Lmc = 1.7  # [m] monochromator-chopper distance
        Lms = Lmc + Lcs  # [m] monochromator-sample distance
        Lhm = 4.123  # [m] distance between the horizontal virtual source and the monochromator

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

        fermi = mycalculator.add_component(
            "fermi_chopper", "FermiChopper", AT=1.700, RELATIVE=Monochromator_Out
        )
        L_bc5_fermi = L_bc5_fermi + fermi.AT_data[2]
        print(L_bc5_fermi)

        fermi.set_parameters(
            radius=0.0006 * 45,
            nslit=45,
            length=0.023,
            w=0.0006,
            yheight=0.113,
            nu="chopper_rpm/60",  # 1/60 to convert RPM into Hz
            verbose=1,
            eff=0.86 * 0.8,
            m=0,
            R0=0,  # no super mirror
            zero_time=2,
            # delay=0.007,
            # delay=str(L_bc5_fermi) + " / neutron_velocity",
            #            delay=0.0065,
            # phase=-23.1209 * 2,
        )

        mycalculator.add_parameter(
            "double", "Efoc", comment="Focusing energy", unit="meV", value=0
        )
        #        mycalculator.append_initialize("if(Efoc==0) Efoc=Ei;")
        mycalculator.append_initialize(
            "if(chopper_rpm==0) chopper_rpm = 60*neutron_velocity * tan(DEG2RAD*a2 / 2) / PI / (({Lcs} + {Lsd}*pow((1 - Efoc / Ei),(-3/2))) * (1 - {Lms} / {Lhm}));".format(
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
