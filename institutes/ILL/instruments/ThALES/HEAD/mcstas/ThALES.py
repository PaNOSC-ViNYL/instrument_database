"""

 - The monochromator do not allow transmitted neutrons

TODO: 
 - [ ] Check RV RH for monochromator and analyzer
 - [ ] Understand the A3_offset:
        ThALES.add_parameter("double", "q_x_elastic", value=1.3139)
        ThALES.add_parameter("double", "q_z_elastic", value=0.146)
        ThALES.append_initialize("  A3_offset=atan(q_z_elastic/q_x_elastic)*RAD2DEG; ")
 - [ ] Aspetto che Martin controlli i numeri dal mio grafichetto
 - [ ] Chiedere a Martin i sample environments

"""
import os

MCSTAS_PATH = os.environ["MCSTAS"]
from mcstasscript.interface import functions

my_configurator = functions.Configurator()

my_configurator.set_mcstas_path(MCSTAS_PATH)
my_configurator.set_mcrun_path(MCSTAS_PATH + "/bin/")

# list here all the common parts to be imported

from institutes.ILL.sources.HEAD.mcstas import Full as source
from libpyvinyl.Instrument import Instrument
from mcstasscript.interface import instr
import pint

# from libpyvinyl import ureg
from pint import set_application_registry

ureg = pint.get_application_registry()

# for conversion from degree to radians and viceversa
import math

# utilities


def def_instrument():
    """Return the constructed instrument description"""
    myinstr = Instrument("ThALES", instrument_base_dir=".")

    # what are these?
    gR0 = 1
    gQc = 0.0216
    gAlpha = 4.07
    gW = 1.0 / 300.0
    Al_Thickness = 0.001
    gGap = 0.001

    ThALES = instr.McStas_instr("ThALES")
    myinstr.add_calculator(ThALES)

    a2 = ThALES.add_parameter(
        "double",
        "a2",
        comment="Angle between beam reflected by monochromator and incident beam",
        unit="degree",
        value=33,
    )
    a2.add_interval(33, 128, True)

    a3 = ThALES.add_parameter(
        "double", "a3", comment="sample table rotation angle", unit="degree", value=0
    )
    a4 = ThALES.add_parameter(
        "double",
        "a4",
        comment="Angle between reflected by sample and incident beams",
        unit="degree",
        value=0,
    )
    a4.add_interval(-128, 128, True)

    a6 = ThALES.add_parameter(
        "double",
        "a6",
        comment="Angle between reflected by analyzer and incident beams",
        unit="degree",
        value=33,
    )
    monochromator_d = ThALES.add_parameter(
        "double",
        "monochromator_d",
        comment="Monochromator lattice parameter",
        value=3.355,
    )

    #   ThALES.add_parameter("double", "q_x_elastic", value=1.3139)
    #   ThALES.add_parameter("double", "q_z_elastic", value=0.146)
    ThALES.add_parameter("double", "SAMPLE", value=0.0)

    #    ThALES.add_declare_var("double", "L")
    ThALES.add_declare_var("double", "flag")
    #    ThALES.add_declare_var("double", "A5")
    #    ThALES.add_declare_var("double", "A6")
    #    ThALES.add_declare_var("double", "final_lambda")
    ThALES.add_declare_var("double", "sample_select")

    #   ThALES.add_declare_var("double", "A3_offset")
    #   ThALES.add_declare_var("int", "elastic_flag_instr_1")
    #   ThALES.add_declare_var("int", "elastic_flag_instr_2")

    ThALES_L = 2.000  # distance between monochromator and sample
    dist_sample_ana = 1.260  # distance between sample and analyzer
    dist_ana_det = 0.640  # distance between analyzer and detector

    ThALES.append_initialize("  flag         = 0; ")
    #    ThALES.append_initialize("  A3_offset=0; ")
    #    ThALES.append_initialize("  A3_offset=atan(q_z_elastic/q_x_elastic)*RAD2DEG; ")

    Origin = ThALES.add_component("Origin", "Progress_bar")
    Origin.set_AT(["0", "0", "0"], RELATIVE="ABSOLUTE")

    # imported source and associated parameters: check the source file!
    # - lambda
    # - dlambda
    # - Ei
    # - dE
    HCS = source.HCS_source(ThALES)
    # override the value of lambda by the value of the angle
    ThALES.append_initialize("lambda =2*sin(0.5*a2*DEG2RAD)*monochromator_d;")
    ThALES.append_initialize('printf("lambda: %.2f\\n", lambda);')
    ThALES.append_initialize("dlambda=0.1*lambda;")
    # the following python string is an alternative way of doing:
    ThALES.add_declare_var("double", "energy")
    ThALES.append_initialize("energy = 81.80421/(lambda*lambda);")
    ThALES.add_declare_var("double", "denergy")
    ThALES.append_initialize("denergy = 2*energy*dlambda/lambda;")
    # lambda2energy = "81.80421/(" + wavelength.name + "*" + wavelength.name + ")"

    H5 = ThALES.add_component("H5", "Arm")
    H5.set_AT([0, 0, 2.155], RELATIVE=HCS)

    H5_rect = ThALES.add_component("H5_rect", "Guide_gravity")
    H5_rect.w1 = 0.170
    H5_rect.h1 = 0.12
    H5_rect.l = 1.0
    H5_rect.R0 = gR0
    H5_rect.Qc = gQc
    H5_rect.alpha = gAlpha
    H5_rect.m = 2
    H5_rect.W = gW
    H5_rect.set_AT([0, 0, 0], RELATIVE=H5)

    H53_origin = ThALES.add_component("H53_origin", "Arm")
    # I don't understand the x
    H53_origin.set_AT([H5_rect.w1 / 2 - 0.06 / 2, 0, H5_rect.l + gGap], RELATIVE=H5)
    H53_origin.set_ROTATED([0, 1.5, 0], RELATIVE=H5)

    H53_start = ThALES.add_component("H53_start", "Arm")
    H53_start.set_AT([0, 0, 0], RELATIVE=H53_origin)

    H53_inpile = ThALES.copy_component("H53_inpile", H5_rect)
    H53_inpile.w1 = 0.06
    H53_inpile.h1 = 0.12
    H53_inpile.l = 4.930 - 3.155
    H53_inpile.m = 3
    H53_inpile.set_AT([0, 0, 0], RELATIVE=H53_start)

    H53_Obt = ThALES.copy_component("H53_Obt", H53_inpile)
    H53_Obt.l = 3
    H53_Obt.set_AT([0, 0, Al_Thickness + 0.015 + H53_inpile.l], RELATIVE=H53_inpile)

    #    H53_Obt_Out = ThALES.add_component("H53_Obt_Out", "Arm")
    #    H53_Obt_Out.set_AT([0, 0, H53_Obt.l + 0.04], RELATIVE=H53_Obt)

    H53_VSComC1 = ThALES.copy_component("H53_VSComC1", H53_inpile)
    H53_VSComC1.l = 7
    H53_VSComC1.nelements = 7
    H53_VSComC1.set_AT([0, 0, H53_Obt.l + 0.075], RELATIVE=H53_Obt)

    H53_Nose = ThALES.add_component("H53_Nose", "Guide_tapering")
    H53_Nose.option = '"parabolical"'
    H53_Nose.w1 = 0.06
    H53_Nose.h1 = 0.12
    H53_Nose.l = 2.0
    H53_Nose.linw = 0.0
    H53_Nose.loutw = 0.7
    H53_Nose.linh = 0.0
    H53_Nose.louth = 0.0
    H53_Nose.R0 = gR0
    H53_Nose.Qcx = gQc
    H53_Nose.Qcy = gQc
    H53_Nose.alphax = gAlpha
    H53_Nose.alphay = gAlpha
    H53_Nose.W = gW
    H53_Nose.mx = 3
    H53_Nose.my = 3
    H53_Nose.segno = 20
    H53_Nose.set_AT([0, 0, H53_VSComC1.l + 0.01], RELATIVE=H53_VSComC1)

    # no length of this?
    before_monochromator_slit = ThALES.add_component(
        "before_monochromator_slit", "Slit"
    )
    before_monochromator_slit.xwidth = 0.04
    before_monochromator_slit.yheight = 0.12
    before_monochromator_slit.set_AT([0, 0, H53_Nose.l + 0.3], RELATIVE=H53_Nose)

    l_monitor = ThALES.add_component("l_monitor", "L_monitor")
    l_monitor.nL = 200
    l_monitor.filename = '"lambda_slit_mono"'
    l_monitor.xwidth = 0.5
    l_monitor.yheight = 0.5
    l_monitor.Lmin = 0
    l_monitor.Lmax = 10
    l_monitor.restore_neutron = 1
    l_monitor.set_AT([0, 0, 0.001], RELATIVE="PREVIOUS")

    H53_ThALES_Monochromator_Cradle = ThALES.add_component(
        "H53_ThALES_Monochromator_Cradle", "Arm"
    )
    H53_ThALES_Monochromator_Cradle.set_AT(
        [0, 0, H53_Nose.l + 0.3 + 2], RELATIVE=H53_Nose
    )

    # double check the probability of reflection
    H53_ThALES_Monochromator = ThALES.add_component(
        "H53_ThALES_Monochromator", "Monochromator_curved"
    )
    H53_ThALES_Monochromator.reflect = '"HOPG.rfl"'
    H53_ThALES_Monochromator.transmit = '"HOPG.trm"'
    H53_ThALES_Monochromator.gap = 0.0005
    H53_ThALES_Monochromator.NH = 13
    H53_ThALES_Monochromator.NV = 13
    H53_ThALES_Monochromator.mosaich = 30
    H53_ThALES_Monochromator.mosaicv = 30
    H53_ThALES_Monochromator.r0 = 1
    H53_ThALES_Monochromator.t0 = 0  # remove transmitted neutrons
    H53_ThALES_Monochromator.RV = "2*" + str(ThALES_L) + "*sin(DEG2RAD*a2/2)"
    H53_ThALES_Monochromator.RH = "2*" + str(ThALES_L) + "/sin(DEG2RAD*a2/2)"
    H53_ThALES_Monochromator.DM = monochromator_d
    H53_ThALES_Monochromator.width = 0.25
    H53_ThALES_Monochromator.height = 0.2
    H53_ThALES_Monochromator.append_EXTEND("flag=SCATTERED;")
    H53_ThALES_Monochromator.set_AT(
        [0, 0, 0], RELATIVE="H53_ThALES_Monochromator_Cradle"
    )
    H53_ThALES_Monochromator.set_ROTATED(
        [0, "a2/2", 0], RELATIVE="H53_ThALES_Monochromator_Cradle"
    )

    H53_ThALES_Monochromator_Out = ThALES.add_component(
        "H53_ThALES_Monochromator_Out", "Arm"
    )
    H53_ThALES_Monochromator_Out.set_AT(
        [0, 0, 0], RELATIVE=H53_ThALES_Monochromator_Cradle
    )
    H53_ThALES_Monochromator_Out.set_ROTATED(
        # [0, 2 * ThALES_A1, 0], RELATIVE="H53_ThALES_Monochromator_Cradle"
        [0, "a2", 0],
        RELATIVE=H53_ThALES_Monochromator_Cradle,
    )
    PSD_cyl = ThALES.add_component("PSD_cyl", "PSDcyl_monitor")
    PSD_cyl.yheight = 0.2
    PSD_cyl.radius = 0.1
    PSD_cyl.nr = 360
    PSD_cyl.ny = 10
    PSD_cyl.filename = '"PSD_cyl.dat"'
    PSD_cyl.restore_neutron = 1
    PSD_cyl.set_AT([0, 0, 0], RELATIVE="PREVIOUS")

    before_sample_slit = ThALES.add_component("before_sample_slit", "Slit")
    before_sample_slit.xwidth = 0.03
    before_sample_slit.yheight = 0.028
    before_sample_slit.set_AT(
        [0, 0, ThALES_L - 0.250], RELATIVE="H53_ThALES_Monochromator_Out"
    )

    E_sample_mon = ThALES.add_component("E_sample_mon", "E_monitor")
    E_sample_mon.nE = 200
    E_sample_mon.filename = '"E_sample_mon"'
    E_sample_mon.xwidth = 0.05
    E_sample_mon.yheight = 0.05
    E_sample_mon.Emin = "energy-denergy"
    E_sample_mon.Emax = "energy+denergy"
    E_sample_mon.restore_neutron = 1
    E_sample_mon.set_AT(
        [0, 0, ThALES_L - 0.050], RELATIVE="H53_ThALES_Monochromator_Out"
    )

    PSD_sample_mon = ThALES.add_component("PSD_sample_mon", "PSD_monitor")
    PSD_sample_mon.nx = 200
    PSD_sample_mon.ny = 200
    PSD_sample_mon.filename = '"PSD_sample_mon.dat"'
    PSD_sample_mon.xwidth = 0.05
    PSD_sample_mon.yheight = 0.05
    PSD_sample_mon.restore_neutron = 1
    PSD_sample_mon.set_AT([0, 0, 0.001], RELATIVE="PREVIOUS")

    sample_arm = ThALES.add_component("sample_arm", "Arm")
    sample_arm.set_AT([0, 0, ThALES_L], RELATIVE="H53_ThALES_Monochromator_Out")
    sample_arm.set_ROTATED([0, "a3", 0], RELATIVE="H53_ThALES_Monochromator_Out")

    res_sample = ThALES.add_component("res_sample", "Res_sample")
    res_sample.thickness = 0.001
    res_sample.radius = 0.005
    res_sample.E0 = 5
    res_sample.dE = 0.25
    res_sample.focus_xw = 0.03
    res_sample.focus_yh = 0.04
    res_sample.yheight = 0.05
    res_sample.target_index = 4
    res_sample.set_WHEN("SAMPLE==0")
    res_sample.set_AT([0, 0, 0], RELATIVE="sample_arm")

    v_sample = ThALES.add_component("v_sample", "V_sample")
    v_sample.radius = 0.005
    v_sample.thickness = 0.001
    v_sample.focus_xw = 0.03
    v_sample.focus_yh = 0.04
    v_sample.yheight = 0.05
    v_sample.target_index = 3
    v_sample.set_WHEN("SAMPLE==1")
    v_sample.set_AT([0, 0, 0], RELATIVE="sample_arm")

    Sample_Out = ThALES.add_component("Sample_Out", "Arm")
    Sample_Out.set_AT([0, 0, 0], RELATIVE="sample_arm")
    Sample_Out.set_ROTATED(["0", "a4", " 0"], RELATIVE="sample_arm")

    after_sample_slit = ThALES.add_component("after_sample_slit", "Slit")
    after_sample_slit.xwidth = 0.03
    after_sample_slit.yheight = 0.04
    after_sample_slit.set_AT([0, 0, 0.250], RELATIVE="Sample_Out")

    Ana_Cradle = ThALES.add_component("Ana_Cradle", "Arm")
    Ana_Cradle.set_AT([0, 0, dist_sample_ana], RELATIVE="Sample_Out")

    PSD_analyzer = ThALES.add_component("PSD_analyzer", "PSD_monitor")
    PSD_analyzer.nx = 200
    PSD_analyzer.ny = 200
    PSD_analyzer.filename = '"PSD_ana.dat"'
    PSD_analyzer.xwidth = 0.25
    PSD_analyzer.yheight = 0.25
    PSD_analyzer.restore_neutron = 1
    PSD_analyzer.set_AT([0, 0, 0], RELATIVE="Ana_Cradle")

    analyzer = ThALES.add_component("analyzer", "Monochromator_curved")
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
    analyzer.set_AT([0, 0, 0], RELATIVE="Ana_Cradle")
    analyzer.set_ROTATED([0, "a6*0.5", 0], RELATIVE="Ana_Cradle")

    Ana_Out = ThALES.add_component("Ana_Out", "Arm")
    Ana_Out.set_AT([0, 0, 0], RELATIVE="Ana_Cradle")
    Ana_Out.set_ROTATED([0, "a6", 0], RELATIVE="Ana_Cradle")

    slit = ThALES.add_component("slit", "Slit")
    slit.xwidth = 0.03
    slit.yheight = 0.08
    slit.set_AT([0, 0, 0.340], RELATIVE="Ana_Out")

    PSD_det = ThALES.add_component("PSD_det", "PSD_monitor")
    PSD_det.nx = 200
    PSD_det.ny = 200
    PSD_det.filename = '"PSD_det.dat"'
    PSD_det.xwidth = 0.2
    PSD_det.yheight = 0.2
    PSD_det.restore_neutron = 1
    PSD_det.set_AT([0, 0, dist_ana_det - 0.0001], RELATIVE="Ana_Out")

    res_monitor = ThALES.add_component("res_monitor", "Res_monitor")
    res_monitor.res_sample_comp = "res_sample"
    res_monitor.filename = '"res_monitor"'
    res_monitor.xwidth = 0.05
    res_monitor.yheight = 0.12
    res_monitor.set_WHEN("SAMPLE==0")
    res_monitor.set_AT([0, 0, dist_ana_det], RELATIVE="Ana_Out")

    detector_all = ThALES.add_component("detector_all", "Monitor")
    detector_all.xwidth = 0.05
    detector_all.yheight = 0.12
    detector_all.set_AT([0, 0, dist_ana_det + 0.001], RELATIVE="Ana_Out")

    myinstr.add_master_parameter("a2", {"ThALES": "a2"}, unit="degree")
    myinstr.add_master_parameter("a3", {"ThALES": "a3"}, unit="degree")
    myinstr.add_master_parameter("a4", {"ThALES": "a4"}, unit="degree")
    myinstr.add_master_parameter("a6", {"ThALES": "a6"}, unit="degree")
    myinstr.master["a2"] = 33 * ureg.degree
    myinstr.master["a3"] = 0 * ureg.degree
    myinstr.master["a4"] = 0 * ureg.degree
    myinstr.master["a6"] = 10 * ureg.degree

    return myinstr
