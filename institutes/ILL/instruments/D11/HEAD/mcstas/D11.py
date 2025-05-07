""" """

# ------------------------------ For McStasscript instruments
# import mcstasscript as ms
from mcstasscript.interface import functions

# from mcstasscript.interface import instr

# this is needed to get the location of McStas executables and libraries
my_configurator = functions.Configurator()

# ------------------------------ Importing sources
from institutes.ILL.sources.HEAD.mcstas import Full as source

from institutes.ILL.sources.HEAD.mcstas import Gauss as sourcesimple

# ------------------------------ Mandatory classes to use
from libpyvinyl.Instrument import Instrument
from libpyvinyl.Parameters import Parameter
from mcstas.McStasInstrumentBase import McStasInstrumentBase

# ------------------------------ Extras
# import os  # to add the path of custom mcstas components

# list here all the common parts to be imported
from typing import List, Optional, Any

# for unit conversions
import pint
from pint import set_application_registry

ureg = pint.get_application_registry()


############## Mandatory method
def get_flavours():
    return [
        # "None",
        # "full",
        # "nosection",
        "simple",
        "simpleNS",
    ]


############## Mandatory method
def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    # if flavour in [None, "None", "", "full"]:
    #    return D11(do_section=True)
    # if flavour == "nosection":
    #    return D11(False)
    if flavour == "simple":
        return D11(do_section=True, remove_H15=True)
    if flavour == "simpleNS":
        return D11(do_section=False, remove_H15=True)
    else:
        raise RuntimeError(f"Flavour {flavour} not implement")


class D11(McStasInstrumentBase):
    """:class: Instrument class defining the D11 instrument at ILL"""

    # ------------------------------ utility methods made available for the users

    # ------------------------------ Internal methods (not available to users)
    def add_moving_guide(self, calculator, name, copy_component, AT, RELATIVE, l=None):
        newcomp = calculator.copy_component(
            name, copy_component, AT=AT, RELATIVE=RELATIVE
        )
        if l is not None:
            newcomp.l = l
        newpar = calculator.add_parameter("int", name + "_index", value=0)
        self.add_parameter_to_master(newpar.name, calculator, newpar)
        self.master[newpar.name] = 0
        newcomp.m = newpar

        return newcomp

    def add_slit(self, calc, name, AT, RELATIVE):
        slit = calc.add_component(name, "Slit", AT=AT, RELATIVE=RELATIVE)
        xwidth = calc.add_parameter(
            "double", "{}_xwidth".format(name), comment="Width of the slit", value=0
        )
        yheight = calc.add_parameter(
            "double", "{}_yheight".format(name), comment="Height of the slit", value=0
        )
        slit.set_parameters(xwidth=xwidth, yheight=yheight)
        self.add_parameter_to_master(xwidth.name, calc, xwidth)
        self.add_parameter_to_master(yheight.name, calc, yheight)
        return slit

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
    def __init__(self, do_section=True, remove_H15=True):
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
                # flux=16.0e12,
                flux=1e11,
            )
        else:
            raise RuntimeError(f"H15 not implemented")

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

        # if remove_H15 is False:
        #    mycalculator, lastcomponent = H15(mycalculator, mysource, SourceTarget)
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
        sg30.set_parameters(w1=0.03, w2=0.03, h1=0.05, h2=0.05, l=0.5 - 0.001, m=0.65)

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

        MovableGuideStart = mycalculator.add_component(
            "MovableGuideStart", "Arm", AT=0.66, RELATIVE=sg30
        )

        # ----------------------------------------
        # TODO: set the distance w.r.t. previous element
        S01 = self.add_slit(mycalculator, "S01", AT=0, RELATIVE=MovableGuideStart)

        M01 = mycalculator.add_component("M01", "Monitor_nD", AT=0.1, RELATIVE=S01)
        T01 = mycalculator.add_component("T01", "Guide_gravity", AT=0.1, RELATIVE=M01)
        T01.set_parameters(
            w1=0.03,
            w2=0.03,
            h1=0.05,
            h2=0.05,
            l=2.5,
            chamfers=0.0008,  # TODO: check with Sylvain
            wavy=8e-4,  # TODO: check with Sylvain
            # nelements=1,  # sections
            # nslit=1,  # channels,
            R0=0.995,  # TODO: check with Sylvain
            Qc=0.0218,  # TODO: check with Sylvain
            alpha=4.07,  # TODO: check with Sylvain
            m=1,  # TODO: check with Sylvain
            W=1.0 / 300.0,  # TODO: check with Sylvain
        )
        T01_index = mycalculator.add_parameter(
            "int",
            "T01_index",
            comment="0 if absorbing tube, 1 for reflective tube, 2 for no-tube",
            value=0,
        )
        self.add_parameter_to_master("T01_index", mycalculator, T01_index)
        T01.m = T01_index

        T02 = self.add_moving_guide(
            mycalculator, "T02", T01, AT=0.010 + T01.l, RELATIVE=T01
        )
        T03 = self.add_moving_guide(
            mycalculator, "T03", T01, AT=0.010 + T02.l, RELATIVE=T02
        )

        """
        D02_index, D02 = self.add_multislit(
            mycalculator,
            "D02",
            [
                {"x": 0.045, "y": 0.088, "r": None},  # t1
                {"x": None, "y": None, "r": 0.010},  # t2
                {"x": None, "y": None, "r": 0.020},  # t3
                {"x": None, "y": None, "r": 0.030},  # t4
                {"x": 0.045, "y": 0.045, "r": None},  # b1
                {"x": 0.045, "y": 0.015, "r": None},  # b2
                {"x": None, "y": None, "r": 0.005},  # b3
                {"x": 1, "y": 1, "r": None},  # b4 N/A completely open
            ],
            at=0.030,
            relative=T03,
        )
        D02_index.value = 1
        self.add_parameter_to_master(D02_index.name, mycalculator, D02_index)
        self.master[D02_index.name] = 1

        T04 = self.add_moving_guide(
            mycalculator, "T04", T01, AT=0.020, RELATIVE=D02, l=2.000
        )
        T05 = self.add_moving_guide(
            mycalculator, "T05", T04, AT=0.010 + T04.l, RELATIVE=T04
        )
        T06 = self.add_moving_guide(
            mycalculator, "T06", T05, AT=0.010 + T05.l, RELATIVE=T05, l=2.500
        )

        # ------------------------------
        D03_index, D03 = self.add_multislit(
            mycalculator,
            "D03",
            [
                {"x": 0.045, "y": 0.088, "r": None},  # t1
                {"x": None, "y": None, "r": 0.010},  # t2
                {"x": None, "y": None, "r": 0.020},  # t3
                {"x": None, "y": None, "r": 0.030},  # t4
                {"x": 0.045, "y": 0.045, "r": None},  # b1
                {"x": 0.045, "y": 0.015, "r": None},  # b2
                # {"x": None, "y": None, "r": None},  # b3 USANS
                # {"x": None, "y": None, "r": None},  # b4 USANS
            ],
            at=T06.l,
            relative=T06,
            # align="b",
            # after="mg5",
        )
        self.add_parameter_to_master(D03_index.name, mycalculator, D03_index)
        self.master[D03_index.name] = 1

        T07 = self.add_moving_guide(mycalculator, "T07", T06, AT=0.020, RELATIVE=D03)
        S04 = self.add_slit(mycalculator, "S04", AT=T07.l, RELATIVE=T07)
        T08 = self.add_moving_guide(
            mycalculator, "T08", T07, AT=0.050, RELATIVE=S04, l=1.180
        )
        # Bride épaisse 0.140m
        T09 = self.add_moving_guide(
            mycalculator, "T09", T08, AT=0.140 + T08.l, RELATIVE=T08
        )
        D05_index, D05 = self.add_multislit(
            mycalculator,
            "D05",
            [
                {"x": 0.045, "y": 0.088, "r": None},  # t1
                {"x": None, "y": None, "r": 0.010},  # t2
                {"x": None, "y": None, "r": 0.020},  # t3
                {"x": None, "y": None, "r": 0.030},  # t4
                {"x": 0.045, "y": 0.045, "r": None},  # b1
                {"x": 0.045, "y": 0.015, "r": None},  # b2
                {"x": None, "y": None, "r": 0.005},  # b3
                {"x": 1.000, "y": 1.000, "r": None},  # b4 N/A completely open
            ],
            at=T09.l,
            relative=T09,
            # align="b",
            # after="mg5",
        )
        self.add_parameter_to_master(D05_index.name, mycalculator, D05_index)
        self.master[D05_index.name] = 1

        T10 = self.add_moving_guide(
            mycalculator, "T10", T09, AT=0.020, RELATIVE=D05, l=2.500
        )
        T11 = self.add_moving_guide(
            mycalculator, "T11", T10, AT=0.010 + T10.l, RELATIVE=T10
        )
        D06_index, D06 = self.add_multislit(
            mycalculator,
            "D06",
            [
                {"x": 0.045, "y": 0.088, "r": None},  # t1
                {"x": None, "y": None, "r": 0.010},  # t2
                {"x": None, "y": None, "r": 0.020},  # t3
                {"x": None, "y": None, "r": 0.030},  # t4
                {"x": 0.045, "y": 0.045, "r": None},  # b1
                {"x": 0.045, "y": 0.015, "r": None},  # b2
                # {"x": None, "y": None, "r": None},  # b3 USANS
                # {"x": None, "y": None, "r": None},  # b4 USANS
            ],
            at=T11.l,
            relative=T11,
            # align="b",
            # after="mg5",
        )
        self.add_parameter_to_master(D06_index.name, mycalculator, D06_index)
        self.master[D06_index.name] = 1

        T12 = self.add_moving_guide(mycalculator, "T12", T11, AT=0.020, RELATIVE=D06)
        T13 = self.add_moving_guide(
            mycalculator, "T13", T12, AT=0.010 + T12.l, RELATIVE=T12, l=1.000
        )
        S07 = self.add_slit(mycalculator, "S07", AT=T13.l, RELATIVE=T13)
        T14 = self.add_moving_guide(mycalculator, "T14", T13, AT=0.050, RELATIVE=S07)
        D08_index, D08 = self.add_multislit(
            mycalculator,
            "D08",
            [
                {"x": 0.045, "y": 0.088, "r": None},  # t1
                {"x": None, "y": None, "r": 0.010},  # t2
                {"x": None, "y": None, "r": 0.020},  # t3
                {"x": None, "y": None, "r": 0.030},  # t4
                {"x": 0.045, "y": 0.045, "r": None},  # b1
                {"x": 0.045, "y": 0.015, "r": None},  # b2
                {"x": None, "y": None, "r": 0.005},  # b3
                {"x": 1.000, "y": 1.000, "r": None},  # b4 N/A completely open
            ],
            at=T14.l,
            relative=T14,
            # align="b",
            # after="mg5",
        )
        self.add_parameter_to_master(D08_index.name, mycalculator, D08_index)
        self.master[D08_index.name] = 1

        """
        T14 = mycalculator.add_component("T14", "Arm", AT=2, RELATIVE="PREVIOUS")
        # ------------------------------
        sample_mcpl_arm = mycalculator.add_component(
            "sample_mcpl_arm",
            "Arm",
            #AT=T14.l + 2.5 - 0.05,
            AT=2,
            RELATIVE=T14,
        )

        # ------------------------------------------------------------
        # this new section contains the sample and the sample environment
        mycalculator, sample_mcpl_arm = self.add_new_section(
            "SampleCalc", sample_mcpl_arm, True
        )
        # ------------------------------------------------------------
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
