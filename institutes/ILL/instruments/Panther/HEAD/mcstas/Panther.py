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
 - [ ] How is determined the tmin for the detector acquisition?
 - [ ] fix the time_frame
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
# ureg = pint.UnitRegistry()

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
    # public because needed in the tests
    def distances(self):
        """
        From B. Fak
        Puts Panther-2 distances in a dictionary
        Versions:
        13-Dec-2022: add (trivial) L11=  0
        28-Nov-2022: For Panther-2, replaces old panlen() for Panther-1
        """
        # Known distances
        Lob = 4.6975 # H12 (source) to H12 bouchon 
        Lb1 = 1.4199 # H12 bouchon to BC1
        L11 = 0.0  # BC1 to BC1 (somewhat trivial, but needed)
        L12 = 0.723  # BC1 to BC2
        L23 = 0.972  # BC2 to BC3
        L34 = 0.918  # BC3 to BC4
        L45 = 0.919  # BC4 to BC5
        Lcs = 0.800  # Fermi chopper to sample
        Lms = 2.500  # Mono to sample
        Lsd = 2.500  # Sample to detector
        Lcd = Lcs + Lsd  # 3.200 Fermi chopper to detector
        # Guessed distances
        L1m = 6.407  # BC1 to mono (this distance is somewhat uncertain)
        Lhm = 4.123  # Horizontal virtual slit to monochromator
        # Calculted distances
        L13 = L12 + L23  #  1.695
        L14 = L13 + L34  #  2.613
        L15 = L14 + L45  #  3.532
        L1s = L1m + Lms  #  8.907
        L1c = L1s - Lcs  #  8.107
        L1d = L1s + Lsd  # 11.407
        L5m = L1m - L15  # 2.875 shouldn't it be 2.8311 ?
        """
        # Calculated distances between mono and the choppers
        L1m= L15 + L5m # 6.407
        L2m= L1m - L12 # 5.684
        L3m= L1m - L13 # 4.712
        L4m= L1m - L14 # 3.794
        Lmc= Lms - Lcs # 1.700 Mono to Fermi
        # BCi to Fermi distances
        L1c= L1m + Lmc # 8.107
        L2c= L2m + Lmc # 7.384
        L3c= L3m + Lmc # 6.412
        L4c= L4m + Lmc # 5.494
        L5c= L5m + Lmc # 4.575
        """
        return dict(
            Lom=Lob+Lb1+L1m,
            Lob=Lob,
            Lb1=Lb1,
            L11=L11,
            L12=L12,
            L13=L13,
            L14=L14,
            L15=L15,
            L23=L23,
            L34=L34,
            L45=L45,
            L1m=L1m,
            L1s=L1s,
            L1c=L1c,
            L1d=L1d,
            Lcs=Lcs,
            Lms=Lms,
            Lsd=Lsd,
            Lcd=Lcd,
            Lhm=Lhm,
            L5m=L5m,
        )


    # ------------------------------ The instrument definition goes in the __init__
    def __init__(self, do_section=True, _start_from_Fermi=False):
        """Here the real definition of the instrument is performed"""

        super().__init__("Panther", do_section)

        distances = self.distances()

        def get_settings_table():
        """return predefined table of chopper_rpm w.r.t. energy"""
        eis = [
            [7.5, "", 8],
            [8.75, "", 8],
            [10, "", 8],
            [12.5, "", 8],
            [15, "", 8],
            [19, "", 9.2],
            [30, "", 16],
            [35, "", 16],
            [40, "", 16],
            [50, "", 16],
            [60, "", 16],
            [76.11, "pg006", 18.4],
            [67.5, "pg006", 19],
            [78.75, "pg006", 19],
            [90, "pg006", 20],
            [112.5, "pg006", 20],
            [52, "cu220", 19],
            [60, "cu220", 19],
            [75, "cu220", 19],
            [90, "cu220", 19],
            [110, "cu220", 19],
            [130, "cu220", 19],
            [130, "cu220", 20],
            [123, "cu331", 19],
            [150, "cu331", 19],
        ]
        return eis

        def init_chopper_rpm(mycalculator):
            distances = self.distances()
            eis = get_settings_table()
            mycalculator.append_initialize(
                "if(chopper_rpm==0){\n" + '    printf("Ei = %.2f\\n", Ei);\n'
            )
            for i in eis:
                mycalculator.append_initialize(
                    "    if(Ei>={}) chopper_rpm = {};".format(i[0], int(i[2] * 1e3))
                )
            mycalculator.append_initialize("}")

        # ------------------------------ some local variables
        monochromators = [
            # name, lattice parameter,
            ["pg002", 3.3550],
            ["pg004", 1.6775],
            ["cu220", 1.2763],
            ["pg006", 1.1183],
            ["cu331", 0.8282],
        ]

        # ------------------------------------------------------------
        # Start with a first section and declaring its parameters
        mycalculator, Origin = self.add_new_section("OriginCalc")

        mono_index = mycalculator.add_parameter(
            "string",
            "mono_index",
            comment="Monochromator name. empty = auto: based on the energy",
            value="0",
        )
        mono_index.add_option("0", True)  # auto
        mono_index.add_option("", True)  # auto
        mono_index.add_option('"pg"', True)
        mono_index.add_option('"cu"', True)
        mono_index.add_option(['"{}"'.format(m[0]) for m in monochromators], True)
        self.add_parameter_to_master(mono_index.name, mycalculator, mono_index)

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

        # does not work in McStasScript
        # mycalculator.add_declare_var(
        #     "string",
        #     "mono_names",
        #     array=len(monochromators),
        #     value=[str(x[0]) for x in monochromators],
        #     comment="monochromator name",
        # )

        a2 = mycalculator.add_parameter(
            "double",
            "a2",
            comment="Angle between beam reflected by monochromator and incident beam: internally calculated if a2=0",
            unit="degree",
            value=0,
        )
        a2.add_option(0, True)  # for automatic calculation
        a2.add_interval(36.00, 62.0, True)  # 58.97, True)
        #self.add_parameter_to_master(a2.name, mycalculator, a2)

        chopper_rpm = mycalculator.add_parameter(
            "double",
            "chopper_rpm",
            comment="Fermi chopper speed: internally calculated if value =0",
            unit="",
            value=0,
        )
        chopper_rpm.add_option(0, True)  # default value for automatic setting
        chopper_rpm.add_interval(6000, 30000, True)
        self.add_parameter_to_master(chopper_rpm.name, mycalculator, chopper_rpm)

        chopper_ratio = mycalculator.add_parameter(
            "double", "chopper_ratio", comment="", unit="", value=1
        )
        chopper_ratio.add_option([1, 2, 3], True)
        self.add_parameter_to_master(chopper_ratio.name, mycalculator, chopper_ratio)

        if not _start_from_Fermi:
            mono_rv = mycalculator.add_parameter(
                "double",
                "mono_rv",
                comment="Monochromator vertical focusing: internally calculated if value is negative",
                value=-1,
                unit="",
            )
            mono_rh = mycalculator.add_parameter(
                "double",
                "mono_rh",
                comment="Monochromator horizontal focusing: internally calculated if value is negative",
                value=-1,
                unit="",
            )

        # ------------------------------------------------------------
        # imported source and associated parameters: check the source file!
        # - Ei
        # - dE
        mycalculator.append_initialize("dE = dE * Ei;")
        mycalculator.append_initialize('printf("dE = %.2f\\n", dE);')
        mycalculator.add_declare_var(
            "double",
            "HCS_focus_xw",
            comment="focusing of the source based on rotation of monochromator",
        )  # append initialize after monochromator definition
        HCS = source.HCS_source(mycalculator)

        HCS.E0 = "Ei"
        # HCS.target_index = 2
        HCS.dist = distances["Lom"]
        HCS.focus_xw = "HCS_focus_xw"  # 0.07  # optimized @ 110 meV
        HCS.focus_yh = 0.05 # reset after monochromator is defined to fit its height

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
        mycalculator.parameters["dE"] = 0.06 / 3
        # the 3 is to optimize the simulation a bit, to be checked
        # from the test_diagnostics:
        # dE =

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

        # Lhm = 4.123  # [m] distance between the horizontal virtual source and the monochromator

        a2_interval = a2.get_intervals()[0]
        mycalculator.add_declare_var("double", "mono_d")
        mycalculator.append_initialize('printf("mono_index = %s\\n", mono_index);')
        mycalculator.append_initialize(
            "char mono_names[{}][6] = ".format(len(monochromators)) + "{"
        )
        for i in monochromators:
            mycalculator.append_initialize('"{}",'.format(i[0]))
        mycalculator.append_initialize("};\n")
        mycalculator.append_initialize(
            "mono_d = -1;\n"
            + 'printf("Selecting monochromator...\\n");\n'
            + "int mindex=0;\n"
            + "int set_mono=(1==1);\n"
            + "int set_a2 = (a2<=0);\n"
            + "for(mindex=0; mindex < {nindex} && set_mono; ++mindex)".format(
                nindex=len(monochromators)
            )
            + "{\n"
            + '    if(strcmp(mono_index,"0")==0 || strlen(mono_index)==0 || (strncmp(mono_index,mono_names[mindex],2)==0 && strcmp(mono_index,mono_names[mindex])<=0)){\n'
            + "        //mono_index = mono_names[mindex];\n"
            + "        mono_d = mono_ds[mindex];\n"
            + '        printf("   - trying mono: %s\\n", mono_names[mindex]);\n'
            + "        if(set_a2){\n"
            + "            a2 = 2*asin(lambda/2/mono_d)*RAD2DEG;\n"
            + '            printf("   - tring a2 = %.6f\\n",a2);\n'
            + "        }\n"
            + "        if(a2 >= {amin!s}-1e-2 && a2 <= {amax!s} ) set_mono=0;\n".format(
                amin=a2_interval[0].m, amax=a2_interval[1].m
            )
            + "    }\n"
            + "}\n"
            + "if(mindex=="
            + "{nindex}".format(nindex=len(monochromators))
            + '){printf("ERROR: mono_d not set, maybe mono_index not acceptable\\n");return 1;}\n'
            + 'if(strcmp(mono_index,"0")==0 || strlen(mono_index)==0) mono_index = mono_names[--mindex];\n'
        )

        mycalculator.append_initialize(
            'printf("mono_d = %.4f\\n",mono_d);'
        )  # get the value corresponing to the selected monochromator
        mycalculator.append_initialize(
            'printf("a2 = %.2f\\n",a2);'
        )  # put a warning if A2 does not match

        mycalculator.append_initialize('printf("mono_index = %s\\n", mono_index);')

        if not _start_from_Fermi:
            # ------------------------------------------------------------
            H12 = mycalculator.add_component("H12", "Arm", AT=[0, 0, 0], RELATIVE=HCS)

            H12_bouchon = mycalculator.add_component(
                "B", "Slit", AT=distances["Lob"], RELATIVE=H12
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
                AT=distances["Lb1"],
                RELATIVE=H12_bouchon,
            )
            BC1.set_parameters(
                # nu="chopper_rpm / chopper_ratio /60 /3", #this is how it is in Panther
                nu="chopper_rpm / chopper_ratio /60 *2",  # using this to improve simulation performance
                radius=0.300,
                yheight=0.175,
                nslit=6,
                isfirst=1,
                abs_out=1,
                xwidth=0.150,
                phase=30.0,  # test_BC1 shows that at 30 the first neutrons have time=0, it should be improved
            )

            diaphragm = mycalculator.add_component(
                # "diaphragm1", "Slit", AT=0.5825, RELATIVE=BC1
                "diaphragm1",
                "Slit",
                AT=0.0005,
                RELATIVE=BC1,
            )
            diaphragm.set_parameters(xwidth=0.150, yheight=0.500)

            phase_init = 24
            lastBC = BC1
            for iBC in [2, 3, 4, 5]:
                dist = distances["L1{}".format(iBC)]
                BC = mycalculator.copy_component(
                    "BC{}".format(iBC), BC1, AT=dist, RELATIVE=BC1
                )
                BC.set_parameters(
                    delay="{dist}/neutron_velocity + {phase_init}/({omega})/360".format(
                        dist=dist, phase_init=phase_init, omega=BC1.nu
                    ),
                    isfirst=0,
                )
                del BC.phase
                lastBC = BC

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
                "heavy_diaphragm", "Slit", AT=0.5011, RELATIVE=lastBC
            )
            diaphragm.set_parameters(xwidth=DCH, yheight=0.150)

            # ------------------------------
            # adds monitors at the same position of the previous component

            # AT=12.1485, RELATIVE=HCS)
            Monochromator_Arm = mycalculator.add_component(
                "Monochromator_Arm", "Arm", AT=distances["L5m"], RELATIVE=lastBC
            )

            Monochromator = mycalculator.add_component(
                "Monochromator", "Monochromator_curved", AT=0, RELATIVE=Monochromator_Arm
            )
            Monochromator.set_ROTATED([0, "a2/2", 0], RELATIVE=Monochromator_Arm)
            Monochromator.set_parameters(  # width=0.300, height=0.220
                NH=11,
                NV=15,
                zwidth=0.019,
                yheight=0.019,
                gap=0.0005,
                t0=0,  # remove transmitted neutrons that would be absorbed by the Beamstop
                order=0,  # higher order neutrons would not pass the diaphragm
                r0=1,
                mosaic = "mono_mosaic",  # depends if PG or Cu
                RV = mono_rv,
                RH = mono_rh,
                DM = "mono_d",
            )

            #   Monochromator.reflect = '"HOPG.rfl"'
            #   Monochromator.transmit = '"HOPG.trm"'
            # Monochromator.append_EXTEND("if(flag!=SCATTERED) ABSORB;")

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
                'if(strncmp(mono_index,"cu",2) == 0){ RMV_u=12.7;  RMV_w=0.449; RMH_g=45; RMH_h=221;}'
            )

            ana_focus = 1.0 / (1.0 / distances["Lhm"] + 1.0 / distances["Lms"])
            mycalculator.append_initialize(
                # "if(mono_rv<0) mono_rv = (RMV_u + acos(1- RMV_w/sin(DEG2RAD*a2/2)))/1000.;"
                "if(mono_rv<0) mono_rv = 2 * ({ana_focus}) * sin(DEG2RAD*a2/2) ;".format(
                    ana_focus=ana_focus
                )
            )
            mycalculator.append_initialize(
                # "if(mono_rh<0) mono_rh = (RMH_g + RMH_h * sin(DEG2RAD*a2/2))/1000.;"
                "if(mono_rh<0) mono_rh = 2 * ({ana_focus}) / sin(DEG2RAD*a2/2);".format(
                    ana_focus=ana_focus
                )
            )

            mycalculator.add_declare_var("float", "mono_mosaic")
            mycalculator.append_initialize(
                'mono_mosaic = (strncmp(mono_index,"cu",2) == 0) ? 30 : 50;'  # 0.3: 0.5
            )

            Monochromator_Out = mycalculator.add_component("Monochromator_Out", "Arm")
            Monochromator_Out.set_AT([0, 0, 0], RELATIVE=Monochromator_Arm)
            Monochromator_Out.set_ROTATED([0, "a2", 0], RELATIVE=Monochromator_Arm)

            mycalculator.append_initialize(
                "HCS_focus_xw = {}*sin(a2/2*DEG2RAD)/3;".format(
                    Monochromator.zwidth * Monochromator.NH
                )
            )
            mycalculator.append_initialize(
                'printf("HCS focus xw = %.4f\\n", HCS_focus_xw);'
            )
            HCS.focus_yh = Monochromator.NV * Monochromator.yheight
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
            "fermi_chopper",
            "FermiChopper",
            AT=distances["Lms"] - distances["Lcs"],
            RELATIVE=Monochromator_Out,
        )

        fermi.set_parameters(
            radius=math.sqrt((0.0006 * 45) ** 2 + 0.023 ** 2),
            nslit=45,
            length=0.023,
            w=0.0006,
            yheight=0.113,
            nu="chopper_rpm/60",  # 1/60 to convert RPM into Hz
            verbose=0,
            eff=0.86 * 0.8,
            m=0,
            R0=0,  # no super mirror
            zero_time=0,
            delay="{dist}/neutron_velocity + {phase_init}/({omega})/ 360".format(
                dist=distances["L1c"], phase_init=12, omega="chopper_rpm/60"
            ),
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

        init_chopper_rpm(mycalculator)


        # ------------------------------
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
        before_sample_diaphragm = mycalculator.add_component(
            "before_sample_diaphragm", "Slit", AT=Lbsd, RELATIVE=fermi
        )
        before_sample_diaphragm.set_parameters(xwidth=bsw, yheight=bsh)

        # ------- printing infos
        mycalculator.append_initialize(
            'printf("Chopper_rpm = %.2f\\nchopper_ratio = %.2f\\nNeutron velocity: %.2f\\n", chopper_rpm, chopper_ratio, neutron_velocity);'
        )
        for iBC in [2, 3, 4, 5]:
            bc = mycalculator.get_component("BC{}".format(iBC))
            mycalculator.append_initialize(
                'printf("delay {} = %.2e\\n", {});\n'.format(bc.name, bc.delay)
            )
        mycalculator.append_initialize(
            'printf("delay {} = %.2e\\n", {});\n'.format(fermi.name, fermi.delay)
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
        self._sample_arm.set_AT(distances["Lcs"] - Lbsd, RELATIVE=sample_mcpl_arm)
        self._sample_environment_arm.set_AT(
            distances["Lcs"] - Lbsd, RELATIVE=sample_mcpl_arm
        )

        # default sample
        self.sample_focus(distances["Lsd"], 2, distances["Lsd"])
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
            "detector_arm", "Arm", AT=[0, 0.4, 0], RELATIVE=Sample_Out
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
            verbose=0,
        )

        # ------------------------------------------------------------
        mycalculator, detector_arm = self.add_new_section("DetectorCalc", detector_arm)

        tube_width = 0.022
        theta_bins = 9 * 32 + 8
        theta_min = -5
        angle_increment = (
            math.asin(tube_width / 2.0 / distances["Lsd"]) * 2
        )  # * 180 / math.pi
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

        tmin = mycalculator.add_parameter(
            "double",
            "tdelay",
            value=0,
            comment="Time min in the histogram (us)",
            unit="mus",
        )
        self.add_parameter_to_master(tmin.name, mycalculator, tmin)
        twidth = mycalculator.add_parameter(
            "double", "twidth", value=0, comment="Width of the time bin"
        )
        self.add_parameter_to_master(twidth.name, mycalculator, twidth)

        # if chopper_ratio is 2: nt=1024
        detector = mycalculator.add_component(
            "detector",
            "Cyl_TOF",
            AT=0,
            RELATIVE=detector_arm,
        )

        # if this is a separate calculator
        if not chopper_rpm.name in mycalculator.parameters:
            chopper_rpm = mycalculator.add_parameter(
                "double",
                "chopper_rpm",
                comment="Fermi chopper speed",
                unit="",
                value=0,
            )
            self.add_parameter_to_master(chopper_rpm.name, mycalculator, chopper_rpm)

            chopper_ratio = mycalculator.add_parameter(
                "double", "chopper_ratio", comment="", unit="", value=1
            )
            self.add_parameter_to_master(
                chopper_ratio.name, mycalculator, chopper_ratio
            )

            init_chopper_rpm(mycalculator)

        mycalculator.append_initialize(
            "if(twidth<=0) twidth = chopper_ratio/2.0/chopper_rpm/60.0;"
        )
        detector.set_parameters(
            filename='"{}"'.format(detector.name),
            ny=ny,
            nt=nt,
            yheight=2.0,
            radius=distances["Lsd"],
            phimin=-17.328,
            # phimax=136,
            tmin="{} * 1e-6".format(
                tmin.name
            ),  # "({Lcs}+{Lsd})/neutron_velocity".format(Lcs=Lcs, Lsd=Lsd),
            tmax="({}+{}) * 1e-6".format(tmin.name, twidth.name),  # time_frame,
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
