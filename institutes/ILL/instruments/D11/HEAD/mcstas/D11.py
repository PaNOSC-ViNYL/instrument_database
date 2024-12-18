"""
 Written by: FARHI Emmanuel (farhi@ill.fr)
 Date: 2012
 Origin:ILL
 Release: McStas 2.5
 Version: $Revision: 1.0 $
 %INSTRUMENT_SITE: ILL

 Adapted and modified by Shervin NOURBAKHSH for McStasscript
        
 TODO:
  - [ ] controllare il comportamento del MCPL output:
        da verificare che i neutroni nel file abbiamo una posizione non nulla che e' impostata rispetto all'Arm definito per l'MCPL
  - [X] addmultislit: da debuggare
  - [ ] implementare il nuovo detector installato nel 2021
    - [ ] tube_length ?
"""

# ------------------------------ For McStasscript instruments
import mcstasscript as ms
from mcstasscript.interface import functions
from mcstasscript.interface import instr

# this is needed to get the location of McStas executables and libraries
my_configurator = functions.Configurator()

# ------------------------------ Importing sources
from institutes.ILL.sources.HEAD.mcstas import Full as source

from institutes.ILL.sources.HEAD.mcstas import Gauss as sourcesimple

# from institutes.ILL.sources.HEAD.mcstas import Gauss_div as sourcesimple

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
    return [
        "None",
        "full",
        "nosection",
        "Borkron_1972",
        "Borkron_2003",
        "Borofloat_2001",
        "simple",
        "simplefull",
    ]


############## Mandatory method
def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    movable_guide_config = {
        "Borkron_1972": {
            "l": [20, 0, 0, 0, 0, 4, 3, 3, 2.5, 0, 2.5, 0, 1.5, 0, 1.5, 0, 0],
            "n": [40, 0, 0, 0, 0, 8, 6, 6, 5, 0, 5, 0, 3, 0, 3, 0, 0],
            "g": [0, 0],
            "chamfers": 0.0002,
            "waviness": 2.5e-5,
        },
        "Borkron_2003": {
            "l": [6.5, 0, 6, 7.5, 0, 4, 3, 3, 2.5, 0, 2.5, 0, 1.5, 0, 1.5, 0, 1],
            "n": [13, 0, 12, 15, 0, 8, 6, 6, 5, 0, 5, 0, 3, 0, 3, 0, 2],
            "g": [0.002, 0.002],
            "chamfers": 0.0002,
            "waviness": 1e-4,
        },
        "Borofloat_2001": {
            "l": [0.5, 6, 6, 0.5, 7, 4, 3, 3, 2, 0.5, 2, 0.5, 0.5, 1, 0.5, 1, 1],
            "n": [1, 6, 6, 1, 7, 4, 3, 3, 2, 1, 2, 1, 1, 1, 1, 1, 1],
            "g": [0.002, 0.002],
            "chamfers": 0.0008,
            "waviness": 8e-4,
        },
        # 2e-4                Waviness [rad]
        # 0.2                 Chamfers [mm]
    }

    if flavour in [None, "None", "", "full", "Borkron_1972"]:
        return D11(movable_guide_config["Borofloat_2001"], do_section=True)
    if flavour in movable_guide_config:
        return D11(movable_guide_config[flavour], do_section=True)
    if flavour == "nosection":
        return D11(movable_guide_config["Borofloat_2001"], False)
    if flavour == "simple":
        return D11(
            movable_guide_config["Borofloat_2001"], do_section=True, remove_H15=True
        )
    if flavour == "simplefull":
        return D11(
            movable_guide_config["Borofloat_2001"], do_section=False, remove_H15=True
        )
    else:
        raise RuntimeError(f"Flavour {flavour} not implement")


def H15(mycalculator, mysource, SourceTarget):
    """
    Description of the H15 guide
      return: calculator, lastcomponent
    """
    gElementGap = 0.004
    PinkCarter = mycalculator.add_component(
        "pinkcarter", "Guide_gravity", AT=0, RELATIVE=SourceTarget
    )
    PinkCarter.set_parameters(
        w1=0.038,
        h1=0.2,
        w2=0.032,
        h2=0.2,
        l=3.170,
        R0=0.995,
        Qc=0.0218,
        alpha=4.07,
        m=1,
        W=1.0 / 300.0,
    )
    SourceTarget = PinkCarter
    # mysource.focus_xw = SourceTarget.w1
    # mysource.focus_yh = SourceTarget.h1

    AlWindow2 = mycalculator.add_component(
        "Alw2", "Al_window", AT=3.522 + 0.001, RELATIVE=SourceTarget
    )
    AlWindow2.set_parameters(thickness=0.002)

    # /* Followed by a 228 mm guide (for H1) in Obturator, after at 10 mm gap, ends at 5.952 m from source after 80 mm gap */

    LeadShutter = mycalculator.copy_component(
        "LeadShutter", PinkCarter, AT=3.522 + 0.01, RELATIVE=SourceTarget
    )

    LeadShutter.set_parameters(w1=PinkCarter.w2, l=0.228)

    AlWindow3 = mycalculator.copy_component(
        "Alw3", AlWindow2, AT=0.228 + 0.005, RELATIVE=LeadShutter
    )  # TODO: controllare lo spessore

    # /*-------------------------*/
    # /*  Curved Guide  ("MAN")  */
    # /*-------------------------*/

    # /* curvature radius is 2700 m , start at XR=5.913 m from source */

    CurvedGuideStart = mycalculator.add_component(
        "CurvedGuideStart", "Arm", AT=0.228 + 0.08, RELATIVE=LeadShutter
    )

    # /* curved part 1, 25 elements : should be 25.5 m long, rho=2700 */
    #  double gGuideWidth        = 0.03;
    gLength1 = 25.5  # 25 pieces
    gElmtLength1 = gLength1 / 25
    gCurvatureRadius = 2700.0
    gElmtRot1 = gElmtLength1 / gCurvatureRadius * 180 / math.pi

    cg1 = mycalculator.copy_component(
        "cg1",
        PinkCarter,
        AT=[0, 0, 0],
        RELATIVE=CurvedGuideStart,
        ROTATED=[0, gElmtRot1, 0],
    )
    cg1.set_parameters(w1=0.03, w2=0.03, l=gElmtLength1 - 0.004, m=1)

    for i in range(2, 26):
        cg_next = mycalculator.copy_component(
            "cg" + str(i),
            cg1,
            AT=[0, 0, gElmtLength1],
            RELATIVE="PREVIOUS",
            ROTATED=[0, gElmtRot1, 0],
        )

    # /* gap V.T.E. 0.260 at 27408 */

    AlWindow4 = mycalculator.copy_component(
        "Alw4", AlWindow3, AT=gElmtLength1 + 0.001, RELATIVE="PREVIOUS"
    )

    PSD_VTE = mycalculator.add_component(
        "PSD_VTE", "Monitor_nD", AT=0.13, RELATIVE=AlWindow4
    )
    PSD_VTE.set_parameters(xwidth=cg1.w1, yheight=cg1.h1, options='"xy"')

    AlWindow5 = mycalculator.copy_component(
        "Alw5", AlWindow4, AT=0.26 - 0.003, RELATIVE=AlWindow4
    )

    VTEtoIN6GuideStart = mycalculator.add_component(
        "VTEtoIN6GuideStart", "Arm", AT=0.003, RELATIVE=AlWindow5
    )

    gLength2 = 22.284  # 22 pieces
    gElmtLength2 = gLength2 / 22

    gLength3 = 5.4  # 5 pieces
    gElmtLength3 = gLength3 / 5

    # /* Ni guide -> IN6: 22.284 m long at 27.588, 22 elements */
    sg1 = mycalculator.copy_component(
        "sg1", cg1, AT=[0, 0, 0], RELATIVE=VTEtoIN6GuideStart
    )
    sg1.set_parameters(
        l=gElmtLength2 - gElementGap,
    )
    for i in range(2, 23):
        sg_next = mycalculator.copy_component(
            "sg" + str(i), sg1, AT=gElmtLength2, RELATIVE="PREVIOUS"
        )

    # /* gap 0.3 m  OT, OS IN6 H15 */

    AlWindow6 = mycalculator.copy_component(
        "Alw6", AlWindow2, AT=gElmtLength2 + 0.001, RELATIVE="PREVIOUS"
    )

    PSD_IN6 = mycalculator.copy_component(
        "PSD_IN6", PSD_VTE, AT=0.15, RELATIVE=AlWindow6
    )

    AlWindow7 = mycalculator.copy_component(
        "Alw7", AlWindow2, AT=0.29, RELATIVE=AlWindow6
    )

    IN6toD7GuideStart = mycalculator.add_component(
        "IN6toD7GuideStart", "Arm", AT=0.30, RELATIVE=AlWindow6
    )

    # /* Ni guide -> D7 Carter Man 2 L=5.4 */

    sg23 = mycalculator.copy_component(
        "sg23", sg1, AT=[0, 0, 0], RELATIVE=IN6toD7GuideStart
    )

    sg23.l = gElmtLength3 - gElementGap

    for i in range(24, 28):
        sg_next = mycalculator.copy_component(
            "sg" + str(i), sg1, AT=gElmtLength3, RELATIVE="PREVIOUS"
        )

    # /* gap 0.3 m OS D7/D11 at 55572 */

    AlWindow8 = mycalculator.copy_component(
        "Alw8", AlWindow2, AT=gElmtLength3 + 0.001, RELATIVE="PREVIOUS"
    )

    PSD_D7 = mycalculator.copy_component("PSD_D7", PSD_VTE, AT=0.15, RELATIVE=AlWindow8)

    AlWindow9 = mycalculator.copy_component(
        "Alw9", AlWindow2, AT=0.29, RELATIVE=AlWindow8
    )

    D7toD11GuideStart = mycalculator.add_component(
        "D7toD11GuideStart", "Arm", AT=[0, -0.065, 0.30], RELATIVE=AlWindow8
    )

    # /* Glass guide SPRI 1.25 m h=0.05 AT (0,-0.065,0.3) */

    sg28 = mycalculator.copy_component(
        "sg28", sg1, AT=[0, 0, 0], RELATIVE=D7toD11GuideStart
    )

    sg28.set_parameters(h1=0.05, h2=0.05, l=1.25 - gElementGap, m=0.65)
    # /* Glass guide  0.68 h=0.05 */

    sg29 = mycalculator.copy_component("sg29", sg28, AT=1.25, RELATIVE=sg28)
    sg29.l = 0.5 - gElementGap

    # /* Velocity selector 0.3. Path in atm is 7 cm */

    AlWindow10 = mycalculator.copy_component(
        "Alw10", AlWindow2, AT=0.5 + 0.001, RELATIVE=sg29
    )
    AlWindow10.thickness = 0.004

    return mycalculator, AlWindow10


class D11(McStasInstrumentBase):
    """:class: Instrument class defining the D11 instrument at ILL"""

    # ------------------------------ utility methods made available for the users

    # ------------------------------ Internal methods (not available to users)
    gElementGap = 0.004
    # attenuators
    attenuation_values = [
        1,  # no attenuator (attenuator out)
        8.325,  # attenuator 1
        26.21,  # attenuator 2
        72.23,  # attenuator 3
        216.5,  # attenuator 1+2
        594.6,  # attenuator 1+3
        1702,  # attenuator 2+3
        13480,  # attenuator 1+2+3
    ]
    collimation_options = [
        40.5,
        # 37,
        34,
        # 31,
        28,
        # 24,
        20.5,
        16.5,
        13.5,
        10.5,
        8,
        5.5,
        4,
        2.5,
        1.5,  # non presente in Nomad
    ]

    # ------------------------------ The instrument definition goes in the __init__
    def __init__(self, movable_guide_config, do_section=True, remove_H15=False):
        """Here the real definition of the instrument is performed"""

        super().__init__("D11", do_section)

        # ------------------------------ some local variables

        # ------------------------------------------------------------
        # Start with a first section and declaring its parameters
        mycalculator, Origin = self.add_new_section("OriginCalc")

        # ================ Distances
        if remove_H15:
            mysource = sourcesimple.VCS_source(mycalculator)
            mysource.set_parameters(
                xwidth=0.01,
                yheight=0.05,
                # zdepth=0.1,
                focus_xw=0.03,
                focus_yh=0.05,
                dist=2.3,
                # flux=2.22e12,
                flux=16.0e12,
            )
        else:
            mysource = source.VCS_source(mycalculator)
            mysource.set_parameters(
                xwidth=0.10,
                zdepth=0.1,
            )

        lambda0 = mycalculator.parameters["lambda"]
        lambda0.value = 6 * ureg.angstrom
        self.add_parameter_to_master("lambda", mycalculator, lambda0)
        self.master["lambda"].add_interval(0.12, 12, True)

        # del mycalculator.parameters["Ei"]
        # del mycalculator.parameters["dE"]
        # del mycalculator.parameters["lambda"]
        mycalculator.parameters["dlambda"].value = 0.05

        # mycalculator.add_declare_var("double", "lambda")
        # mycalculator.append_initialize("lambda = sqrt(81.80421036/Ei);")

        mycalculator.add_declare_var("double", "neutron_velocity")
        mycalculator.append_initialize("neutron_velocity = 3956.034012/lambda;")
        mycalculator.append_initialize('printf("nv = %2f\\n", neutron_velocity);')

        mycalculator.append_initialize('printf("lambda = %.2f\\n", lambda);')
        mycalculator.append_initialize(
            'dlambda = dlambda*lambda;printf("dlambda = %.2f\\n", dlambda);'
        )

        AlWindow1 = mycalculator.add_component(
            "Alw1", "Al_window", AT=2.33, RELATIVE=mysource
        )
        AlWindow1.set_parameters(thickness=0.002)

        SourceTarget = mycalculator.add_component(
            "SourceTarget", "Arm", AT=AlWindow1.thickness, RELATIVE=AlWindow1
        )
        # mysource.dist = AlWindow1.AT_data[2]

        if remove_H15 is False:
            mycalculator, lastcomponent = H15(mycalculator, mysource, SourceTarget)
        # ------------------------------
        velocity_selector_mcpl_arm = mycalculator.add_component(
            "velocity_selector_mcpl_arm",
            "Arm",
            AT=0,
            RELATIVE="PREVIOUS",
        )

        # ------------------------------------------------------------
        if remove_H15 is False:
            mycalculator, velocity_selector_arm = self.add_new_section(
                "VelocityCalc", velocity_selector_mcpl_arm
            )
        else:
            velocity_selector_arm = velocity_selector_mcpl_arm

        if do_section is True and remove_H15 is False:
            lambda1 = mycalculator.add_parameter(
                "double", lambda0.name, unit=lambda0.unit, comment=lambda0.comment
            )
            self.add_parameter_to_master(lambda0.name, mycalculator, lambda1)

        Vrpm = mycalculator.add_parameter(
            "double", "rpm", comment="velocity selector RPM", value=0
        )
        Vrpm.add_option(0, True)
        Vrpm.add_interval(None, 3100, False)
        Vrpm.add_interval(28300, None, False)
        Vrpm.add_interval(9000, 11900, False)  # resonance speed for D11

        Dolores = mycalculator.add_component(
            "Dolores", "V_selector", AT=0.025 + 0.30 / 2, RELATIVE=velocity_selector_arm
        )
        Dolores.set_parameters(
            xwidth=0.03,
            yheight=0.05,
            zdepth=0.30,
            radius=0.12,
            alpha=48.298,
            length=0.250,
            d=0.0004,
            nu="{:s}/60".format(Vrpm.name),
            nslit=72,
        )
        mycalculator.append_initialize(
            "if({}==0) {} = 60*3956*{:.6f}*DEG2RAD/2/PI/{}/{:.6f};".format(
                Vrpm.name, Vrpm.name, Dolores.alpha, lambda0.name, Dolores.length
            )
        )
        mycalculator.append_initialize(
            'printf("VS rpm = %.2f\\n", {});'.format(Vrpm.name)
        )
        # if remove_H15:
        # mysource.dist = (
        #    61.97  # calculated as the distance with H15 between Dolores and source
        # )
        #            mysource.dist = Dolores.AT_data[2]
        # mysource.focus_xw = Dolores.xwidth
        # mysource.focus_yh = Dolores.yheight

        AlWindow11 = mycalculator.add_component(
            "Alw11", "Al_window", AT=0.15 + 0.01, RELATIVE=Dolores
        )
        AlWindow11.thickness = 0.004  # AlWindow10.thickness

        # /* Glass guide  0.50 h=0.05 */

        sg30 = mycalculator.add_component(
            "sg30", "Guide_gravity", AT=0.15 + 0.02, RELATIVE=Dolores
        )
        sg30.set_parameters(
            w1=0.03, w2=0.03, h1=0.05, h2=0.05, l=0.5 - self.gElementGap, m=0.65
        )

        # /* Gap 16 cm, start of movable guide */
        AlWindow12 = mycalculator.copy_component(
            "Alw12", AlWindow11, AT=0.5 + 0.001, RELATIVE=sg30
        )
        AlWindow12.thickness = AlWindow1.thickness

        # COMPONENT Mon_D11_Out = Monitor_nD(xwidth=gGuideWidth, yheight=gGuideHeight2,
        #  options=xlmonopts)
        # AT (0,0,0.5+0.01) RELATIVE sg30

        AlWindow13 = mycalculator.copy_component(
            "Alw13", AlWindow12, AT=0.5 + 0.03 - 0.003, RELATIVE=sg30
        )

        # ----------------------------------------
        CollimationCalc = mycalculator
        collimation = mycalculator.add_parameter(
            "double",
            "collimation",
            comment="Collimation length: free path between end of the guide and sample",
            unit="m",
            value=1.5,
        )
        self.add_parameter_to_master(collimation.name, mycalculator, collimation)
        collimation.add_option(self.collimation_options, True)

        mycalculator.add_declare_var(
            "double",
            "collimation_options",
            array=len(self.collimation_options),
            value=self.collimation_options,
            comment="accepted values for collimation",
        )

        icollimation = mycalculator.add_declare_var(
            "int",
            "icollimation",
            comment="index of the chosen collimation within the array of accepted values",
            value=0,
        )
        mycalculator.append_initialize(
            "while(collimation_options[icollimation]>collimation)icollimation++;"
        )
        mycalculator.append_initialize(
            "if(collimation!=collimation_options[icollimation]){"
        )
        mycalculator.append_initialize(
            'printf("[ERROR] chosen collimation not within accepted values, exiting\\n");exit(EXIT_FAILURE);}'
        )
        # ----------------------------------------

        disk_index = self.add_multislit(
            mycalculator,
            "disk6",
            [
                {"x": None, "y": None, "r": 0.010},  # 0°
                {"x": 0.035, "y": 0.055, "r": None},  # 45°
            ],
            0.65,
            sg30,
        )
        disk_index.value = 1
        # self.add_parameter_to_master(disk_index.name, mycalculator, disk_index)
        # self.master[disk_index.name] = 1

        MovableGuideStart = mycalculator.add_component(
            "MovableGuideStart", "Arm", AT=0.66, RELATIVE=sg30
        )

        # /* D11 Movable guide start */
        microGap = 0.0001

        def inactive_coll(i, collimation_length, movable_guide_config):
            # print(
            #     i,
            #     "("
            #     + str(collimation_length)
            #     + " - collimation)>0 ? "
            #     + str(movable_guide_config["n"][i])
            #     + " : 0",
            # )

            return (
                "("
                + str(collimation_length)
                + " - collimation)>0 ? "
                + str(movable_guide_config["n"][i])
                + " : 0"
            )

        collimation_length = self.collimation_options[0]
        mg0 = mycalculator.copy_component("mg0", sg30, AT=0, RELATIVE=MovableGuideStart)
        mg0.set_parameters(
            l=movable_guide_config["l"][0],
            chamfers_tb=movable_guide_config["chamfers"],
            chamfers_z=movable_guide_config["chamfers"],
            nelements=inactive_coll(0, collimation_length, movable_guide_config),
            wavy=movable_guide_config["waviness"],
        )
        collimation_length = collimation_length - movable_guide_config["l"][0]

        mg1 = mycalculator.copy_component(
            "mg1", mg0, AT=microGap + movable_guide_config["l"][0], RELATIVE="PREVIOUS"
        )
        mg1.set_parameters(
            l=movable_guide_config["l"][1],
            nelements=inactive_coll(1, collimation_length, movable_guide_config),
        )
        collimation_length = collimation_length - movable_guide_config["l"][1]

        for i in range(2, 4):
            mg_next = mycalculator.copy_component(
                "mg" + str(i),
                mg1,
                AT=microGap
                + movable_guide_config["l"][i - 1]
                + movable_guide_config["g"][i - 2],
                RELATIVE="PREVIOUS",
            )
            mg_next.set_parameters(
                l=movable_guide_config["l"][i],
                nelements=inactive_coll(i, collimation_length, movable_guide_config),
            )
            collimation_length = collimation_length - movable_guide_config["l"][i]

        dist = [
            0,
            0,
            0,
            0,
            0.0,
            0.017,
            0.002,
            0.002,
            0.017,
            0.0,
            0.002,
            0.0,
            0.017,
            0.0,
            0.002,
            0.0,
            0.002,
        ]
        for i in range(4, 17):
            mg_next = mycalculator.copy_component(
                "mg" + str(i),
                mg1,
                AT=microGap + movable_guide_config["l"][i - 1] + dist[i],
                RELATIVE="PREVIOUS",
            )
            mg_next.set_parameters(
                l=movable_guide_config["l"][i],
                nelements=inactive_coll(i, collimation_length, movable_guide_config),
            )
            collimation_length = collimation_length - movable_guide_config["l"][i]

        # /* Gap 17 mm at 20.5 m collimation */
        gap = 0.001
        # ------------------------------ Disk 5
        disk_index = self.add_multislit(
            mycalculator,
            "disk5",
            [
                {"x": 0.045, "y": 0.082, "r": None},  # 0°
                {"x": 0.050, "y": 0.055, "r": None},  # 45°
                {"x": None, "y": None, "r": 0.030},  # 90°
                {"x": None, "y": None, "r": 0.020},  # 135°
                {"x": None, "y": None, "r": 0.010},  # 180°
                {"x": None, "y": None, "r": 0.000},  # 225°
                {"x": 0.0385, "y": 0.055, "r": None},  # 270°
                {"x": None, "y": None, "r": 0.000},  # 315°
            ],
            movable_guide_config["l"][5] + gap,
            "mg5",
            align="b",
            after="mg5",
        )
        self.add_parameter_to_master(disk_index.name, mycalculator, disk_index)
        self.master[disk_index.name] = 1

        # ------------------------------ Disk 4
        disk_index = self.add_multislit(
            mycalculator,
            "disk4",
            [
                {"x": 0.043, "y": 0.080, "r": None},  # 0°
                {"x": 0.050, "y": 0.055, "r": None},  # 45°
                {"x": None, "y": None, "r": 0.030},  # 90°
                {"x": None, "y": None, "r": 0.020},  # 135°
                {"x": None, "y": None, "r": 0.010},  # 180°
                {"x": None, "y": None, "r": 0.000},  # 225°
                {"x": 0.0395, "y": 0.055, "r": None},  # 270°
                {"x": None, "y": None, "r": 0.000},  # 315°
            ],
            movable_guide_config["l"][8] + gap,
            "mg8",
            align="b",
            after="mg8",
        )
        self.add_parameter_to_master(disk_index.name, mycalculator, disk_index)
        self.master[disk_index.name] = 1
        # ------------------------------ Disk 3
        disk_index = self.add_multislit(
            mycalculator,
            "disk3",
            [
                {"x": 0.050, "y": 0.080, "r": None},  # 0°
                {"x": 0.050, "y": 0.055, "r": None},  # 45°
                {"x": None, "y": None, "r": 0.030},  # 90°
                {"x": None, "y": None, "r": 0.020},  # 135°
                {"x": None, "y": None, "r": 0.010},  # 180°
                {"x": None, "y": None, "r": 0.005},  # 225°
                {"x": 0.038, "y": 0.055, "r": None},  # 270°
                {"x": None, "y": None, "r": 0.000},  # 315°
            ],
            movable_guide_config["l"][12] + gap,
            "mg12",
            align="b",
            after="mg12",
        )
        self.add_parameter_to_master(disk_index.name, mycalculator, disk_index)
        self.master[disk_index.name] = 6
        # ------------------------------ Disk 2
        disk_index = self.add_multislit(
            mycalculator,
            "disk2",
            [
                {"x": 0.0365, "y": 0.040, "r": None},  # 0°
                {"x": 0.050, "y": 0.055, "r": None},  # 45°
                {"x": None, "y": None, "r": 0.030},  # 90°
                {"x": None, "y": None, "r": 0.020},  # 135°
                {"x": None, "y": None, "r": 0.010},  # 180°
                {"x": None, "y": None, "r": 0.005},  # 225°
                {"x": 0.0285, "y": 0.031, "r": None},  # 270°
                {"x": None, "y": None, "r": 0.000},  # 315°
            ],
            movable_guide_config["l"][15] + 2.5 - 1.5,
            "mg15",
            align="b",
        )
        self.add_parameter_to_master(disk_index.name, mycalculator, disk_index)
        self.master[disk_index.name] = 6
        # ------------------------------ Disk 1
        disk_index = self.add_multislit(
            mycalculator,
            "disk1",
            [
                {"x": 0.035, "y": 0.035, "r": None},  # 0°
                {"x": 0.020, "y": 0.020, "r": None},  # 45°
                {"x": 0.015, "y": 0.015, "r": None},  # 90°
                {"x": None, "y": None, "r": 0.020},  # 135°
                {"x": None, "y": None, "r": 0.015},  # 180°
                {"x": None, "y": None, "r": 0.010},  # 225°
                {"x": 0.025, "y": 0.040, "r": None},  # 270°
                {"x": 0.030, "y": 0.030, "r": None},  # 315°
            ],
            movable_guide_config["l"][15] + 2.5 - 0.5,
            "mg15",
            align="b",
        )
        self.add_parameter_to_master(disk_index.name, mycalculator, disk_index)
        self.master[disk_index.name] = 2

        # ------------------------------
        sample_mcpl_arm = mycalculator.add_component(
            "sample_mcpl_arm",
            "Arm",
            AT=movable_guide_config["l"][15] + 2.5 - 0.05,
            RELATIVE="mg15",
        )

        sample_diaphragm_radius = mycalculator.add_parameter(
            "double",
            "sample_diaphragm_radius",
            unit="m",
            comment="circular aperture radius",
            value=0,
        )
        sample_diaphragm_height = mycalculator.add_parameter(
            "double",
            "sample_diaphragm_height",
            unit="m",
            comment="circular aperture height",
            value=0,
        )
        sample_diaphragm_width = mycalculator.add_parameter(
            "double",
            "sample_diaphragm_width",
            unit="m",
            comment="circular aperture width",
            value=0,
        )
        self.add_parameter_to_master(
            sample_diaphragm_width.name, mycalculator, sample_diaphragm_width
        )
        self.master[sample_diaphragm_width.name] = 0 * ureg.m
        self.add_parameter_to_master(
            sample_diaphragm_height.name, mycalculator, sample_diaphragm_height
        )
        self.master[sample_diaphragm_height.name] = 0 * ureg.m
        self.add_parameter_to_master(
            sample_diaphragm_radius.name, mycalculator, sample_diaphragm_radius
        )
        self.master[sample_diaphragm_radius.name] = 0 * ureg.m

        sample_diaphragm_x = mycalculator.add_parameter(
            "double",
            "sample_diaphragm_x",
            unit="m",
            comment="x shift w.r.t. beam center",
            value=0,
        )
        sample_diaphragm_y = mycalculator.add_parameter(
            "double",
            "sample_diaphragm_y",
            unit="m",
            comment="y shift w.r.t. beam center",
            value=0,
        )
        self.add_parameter_to_master("sample_x", mycalculator, sample_diaphragm_x)
        self.add_parameter_to_master("sample_y", mycalculator, sample_diaphragm_x)
        # sample_diaphragm_z = mycalculator.add_parameter(
        #     "double",
        #     "sample_diaphragm_z",
        #     unit="m",
        #     comment="z shift w.r.t. beam center",
        #     value=0,
        # )
        # self.add_parameter_to_master("sample_z", mycalculator, sample_diaphragm_x)

        sample_diaphragm = mycalculator.add_component(
            "sample_diaphragm",
            "Slit",
            AT=[
                sample_diaphragm_x,
                sample_diaphragm_y,
                movable_guide_config["l"][15] + 2.5 - 0.015,
            ],
            RELATIVE="mg15",
        )
        sample_diaphragm.set_parameters(
            radius=sample_diaphragm_radius,
            xwidth=sample_diaphragm_width,
            yheight=sample_diaphragm_height,
        )

        # ------------------------------------------------------------
        # this new section contains the sample and the sample environment
        mycalculator, sample_mcpl_arm = self.add_new_section(
            "SampleCalc", sample_mcpl_arm, True
        )
        # ------------------------------------------------------------
        # the sample can be moved in x and y changing the parameters
        # sample_x, sample_y, sample_z in McStas coordiate system
        # w.r.t. the sample_arm position
        self._sample_arm.set_AT(0.05, RELATIVE=sample_mcpl_arm)
        self._sample_environment_arm.set_AT(
            self._sample_arm.AT_data, RELATIVE=sample_mcpl_arm
        )
        # self._sample_environment.set_ROTATED([0, det_angle, 0])

        # self.sample_box_shape(0.02, 0.03, 0.0035, 0.00125)
        # default sample
        detpos = mycalculator.add_parameter(
            "double", "detpos", comment="Detector distance", unit="m", value=2
        )
        detpos.add_interval(1, 28, True)
        self.add_parameter_to_master(detpos.name, mycalculator, detpos)

        # self.sample_focus(8, 3, "detpos") # defined in the detector section to catch the detector size
        sample = self.set_sample_by_name("None")

        Sample_Out = mycalculator.add_component(
            "Sample_Out", "Arm", AT=0, RELATIVE=self._sample_arm
        )

        # ------------------------------------------------------------

        # things might be improved putting a slit at the minimum distance of the detector
        # such that the mcpl saves only neutrons that have a chance to be detected
        mycalculator, center_det = self.add_new_section("DetectorCalc", Sample_Out)
        attenuator_index = mycalculator.add_parameter(
            "int",
            "attenuator_index",
            comment="select the attenuation level by combining attenuator 1,2,3",
            value=6,
        )
        attenuator_index.add_interval(0, len(self.attenuation_values) - 1, True)
        self.add_parameter_to_master(
            attenuator_index.name, mycalculator, attenuator_index
        )

        mycalculator.add_declare_var(
            "double",
            "att_factor",
            array=len(self.attenuation_values),
            value=self.attenuation_values,
        )

        attenuator = mycalculator.add_component(
            "attenuator", "Attenuator", AT=0.01, RELATIVE=center_det
        )

        attenuator.set_parameters(
            scaling="1.0/att_factor[attenuator_index]",
            xwidth=1,  # very large to not miss any neutron
            yheight=1,  # very large to not miss any neutron
        )

        if not detpos.name in mycalculator.parameters:
            print(type(mycalculator.parameters))
            detpos = mycalculator.add_parameter(
                "double", "detpos", comment="Detector distance", unit="m", value=2
            )
            detpos.add_interval(1, 28, True)
            self.add_parameter_to_master(detpos.name, mycalculator, detpos)

        bs_x = mycalculator.add_parameter(
            "double", "bs_x", comment="Beamstop x position", unit="m", value=0.141
        )
        bs_y = mycalculator.add_parameter(
            "double", "bs_y", comment="Beamstop y position", unit="m", value=0.650
        )
        self.add_parameter_to_master(bs_x.name, mycalculator, bs_x)
        self.add_parameter_to_master(bs_y.name, mycalculator, bs_y)

        bs_index = mycalculator.add_parameter(
            "int",
            "bs_index",
            comment="Index to select the beamspot: 0-> beamstop width = 65mm, height=70mm; 1-> beamstop width = 75mm, height=80mm; 2-> beamstop width = 85mm, height=90mm; 3-> beamstop width = 95mm, height=100mm",
            value=1,
        )
        bs_index.add_option([-1, 0, 1, 2, 3], True)
        self.add_parameter_to_master(bs_index.name, mycalculator, bs_index)

        mycalculator.add_declare_var(
            "double", "bs_w", comment="beam stop width", unit="m", value=0
        )
        mycalculator.add_declare_var(
            "double", "bs_h", comment="beam stop width", unit="m", value=0
        )

        mycalculator.append_initialize(
            "if(bs_index == 0){\n"
            + "  bs_w = 0.065; bs_h = 0.070;\n"
            + "} else if (bs_index == 1){\n"
            + "  bs_w = 0.075; bs_h = 0.080;\n"
            + "} else if (bs_index == 2){\n"
            + "  bs_w = 0.085; bs_h = 0.090;\n"
            + "} else if (bs_index == 3){\n"
            + "  bs_w = 0.095; bs_h = 0.100;\n"
            + "} else if ((int)bs_index == -1){\n"
            + "  bs_w = 0.001; bs_h =0.001; bs_x=0; bs_y=0;\n"
            + '} else printf("ERROR: bs_index out of range [-1:-3]\\n");'
        )

        beamstop = mycalculator.add_component(
            "beamstop",
            "Beamstop",
            AT=["bs_x - 0.141", "bs_y - 0.650", "{} - 0.08".format(detpos.name)],
            RELATIVE=center_det,
        )
        beamstop.set_parameters(xwidth="bs_w", yheight="bs_h")
        det_central_ntubes = 192  # number of tubes of the central panel of the detector
        tube_width = 0.008  # width of the tubes
        tube_length = 1.024  # tube length
        det_length = 0  # length of the tubes
        det_lateral_gap = (
            -0.004
        )  # gap between the central and lateral panels: last pixels are overlapping
        det_lateral_ntubes = 32  # number of tubes of the lateral panels of the detector
        det_length_nbins = 256

        detector_central = mycalculator.add_component(
            "detector_central",
            "Monitor_nD",
            AT=[0, -0.032, detpos],
            # -0.048 found comparing beam center of mass in test 0
            RELATIVE=center_det,
        )
        detector_central.set_parameters(
            xwidth=tube_length,
            yheight=det_central_ntubes * tube_width,
            options='"parallel square x bins={} y bins={} file={}"'.format(
                det_length_nbins, det_central_ntubes, "detector_central.dat"
            ),
        )

        # right looking from the source
        detector_right = mycalculator.add_component(
            "detector_right", "Monitor_nD", AT=0, RELATIVE=center_det
        )
        detector_right.set_parameters(
            xwidth=det_lateral_ntubes * tube_width,
            yheight=tube_length,
            options='"parallel square x bins={} y bins={} file={}"'.format(
                det_lateral_ntubes, det_length_nbins, "detector_right.dat"
            ),
        )
        detector_right.set_AT(
            [
                detector_central.xwidth / 2
                + det_lateral_gap
                + detector_right.xwidth / 2,
                0,
                "{} - 0.10".format(detpos.name),
            ]
        )
        detector_left = mycalculator.copy_component(
            "detector_left",
            detector_right,
            AT=[
                -detector_right.AT_data[0],
                detector_right.AT_data[1],
                detector_right.AT_data[2],
            ],
            RELATIVE=center_det,
        )
        detector_left.options = '"parallel square x bins={} y bins={} file={}"'.format(
            det_lateral_ntubes, det_length_nbins, "detector_left.dat"
        )

        # the sample focusing will be done based on the total width of the detectors
        # and the height of the tallest detector
        # the distance is determined by the detpos parameter
        self.sample_focus(
            detector_central.xwidth + detector_left.xwidth + detector_right.xwidth,
            max(detector_central.yheight, detector_left.yheight),
            detpos,
        )

        # ------------------------------ instrument parameters

    def set_test(self, test_number: Optional[int] = None):
        myinstrument = self
        myinstrument.master["lambda"] = 6 * ureg.angstrom
        myinstrument.master["detpos"] = 2 * ureg.m
        myinstrument.master["attenuator_index"] = 0
        myinstrument.master["collimation"] = 8 * ureg.m
        myinstrument.master["bs_index"] = 0
        myinstrument.sample_holder(
            material="quartz", shape="box", w=0.02, h=0.03, d=0.0135, th=0.00125
        )
        myinstrument.sample_shape("holder")
        if test_number == 0:  # direct attenuated beam
            myinstrument.set_sample_by_name("None")
            myinstrument.sample_holder(None, None)
            myinstrument.master["attenuator_index"] = 6
            myinstrument.master["bs_index"] = -1
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
