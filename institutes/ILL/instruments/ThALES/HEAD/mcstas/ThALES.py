"""
This McStasScript file was generated from a
McStas instrument file. It is advised to check
the content to ensure it is as expected.


TODO: 
 - [ ] libpyvinyl parameters and master parameters
 - [ ] test run with result validation
 - [ ] 
"""
from mcstasscript.interface import instr
from mcstas.institutes.ILL.sources.HEAD import Full as source


def def_instrument():
    ThALES = instr.McStas_instr("ThALES_generated")
    ThALES.add_parameter("double", "ThALES_dE", value=0.5)
    ThALES.add_parameter("double", "A3", value=0.0)
    ThALES.add_parameter("double", "A4", value=20.0)
    ThALES.add_parameter("double", "Ei", value=5.0)
    ThALES.add_parameter("double", "Ef", value=5.0)
    ThALES.add_parameter("double", "q_x_elastic", value=1.3139)
    ThALES.add_parameter("double", "q_z_elastic", value=0.146)
    ThALES.add_parameter("double", "SAMPLE", value=0.0)
    ThALES.add_declare_var("double", "sT3")
    ThALES.add_declare_var("double", "sI3")
    ThALES.add_declare_var("double", "sT2")
    ThALES.add_declare_var("double", "sI2")
    ThALES.add_declare_var("double", "sT1")
    ThALES.add_declare_var("double", "sI1")
    ThALES.add_declare_var("double", "gR0")
    ThALES.add_declare_var("double", "gQc")
    ThALES.add_declare_var("double", "gAlpha")
    ThALES.add_declare_var("double", "gW")
    ThALES.add_declare_var("double", "Al_Thickness")
    ThALES.add_declare_var("double", "gGap")
    ThALES.add_declare_var("double", "ThALES_DM")
    ThALES.add_declare_var("double", "ThALES_A1")
    ThALES.add_declare_var("double", "ThALES_L")
    ThALES.add_declare_var("double", "L")
    ThALES.add_declare_var("double", "flag")
    ThALES.add_declare_var("double", "A5")
    ThALES.add_declare_var("double", "A6")
    ThALES.add_declare_var("double", "ThALES_lambda")
    ThALES.add_declare_var("double", "delta_E")
    ThALES.add_declare_var("double", "final_lambda")
    ThALES.add_declare_var("double", "sample_select")
    ThALES.add_declare_var("double", "ThALES_RMV")
    ThALES.add_declare_var("double", "ThALES_RMH")
    ThALES.add_declare_var("double", "A3_offset")
    ThALES.add_declare_var("double", "dist_ana_det")
    ThALES.add_declare_var("double", "dist_sample_ana")
    ThALES.add_declare_var("double", "ana_RMV")
    ThALES.add_declare_var("double", "ana_RMH")
    ThALES.add_declare_var("int", "elastic_flag_instr_1")
    ThALES.add_declare_var("int", "elastic_flag_instr_2")
    ThALES.append_initialize("  gR0 = 1; ")
    ThALES.append_initialize("  gQc = 0.0216; ")
    ThALES.append_initialize("  gAlpha = 4.07; ")
    ThALES.append_initialize("  gW = 1.0/300.0; ")
    ThALES.append_initialize("   ")
    ThALES.append_initialize("  Al_Thickness = 0.001; ")
    ThALES.append_initialize("  gGap = 0.001; ")
    ThALES.append_initialize("   ")
    ThALES.append_initialize("  ThALES_DM      = 3.355; /* PG002 */ ")
    ThALES.append_initialize("  ThALES_A1      = 0; ")
    ThALES.append_initialize("  ThALES_L       = 2.000; ")
    ThALES.append_initialize("   ")
    ThALES.append_initialize("  flag         = 0; ")
    ThALES.append_initialize("   ")
    ThALES.append_initialize("  A5=0; ")
    ThALES.append_initialize("  A6=0; ")
    ThALES.append_initialize("  delta_E=0; ")
    ThALES.append_initialize("  ThALES_RMV = -1; ")
    ThALES.append_initialize("  ThALES_RMH=-1; ")
    ThALES.append_initialize("  A3_offset=0; ")
    ThALES.append_initialize("dist_ana_det=0.640; ")
    ThALES.append_initialize("dist_sample_ana=1.260; ")
    ThALES.append_initialize("  ana_RMV=-1; ")
    ThALES.append_initialize("  ana_RMH=-1; ")
    ThALES.append_initialize("       ")
    ThALES.append_initialize("  ThALES_lambda=1/(0.1106*sqrt(Ei));  ")
    ThALES.append_initialize("  final_lambda=1/(0.1106*sqrt(Ef)); ")
    ThALES.append_initialize("  A5  = -asin(final_lambda/2/ThALES_DM)*RAD2DEG; ")
    ThALES.append_initialize("  A6  = 2*A5; ")
    ThALES.append_initialize("  ThALES_A1  = asin(ThALES_lambda/2/ThALES_DM)*RAD2DEG; ")
    ThALES.append_initialize("   ")
    ThALES.append_initialize("   ")
    ThALES.append_initialize("     ")
    ThALES.append_initialize("  if (ThALES_RMV<0) ")
    ThALES.append_initialize(
        "    ThALES_RMV = 1/(    (   1/ThALES_L + 1/7.0 )   / (2*sin(DEG2RAD*fabs(ThALES_A1)))    ); "
    )
    ThALES.append_initialize("  if (ThALES_RMH<0) ")
    ThALES.append_initialize(
        "    ThALES_RMH = 1/(    (   1/ThALES_L + 1/2.0 )  *sin(DEG2RAD*fabs(ThALES_A1)) / (2)    ); "
    )
    ThALES.append_initialize("   ")
    ThALES.append_initialize(
        '  printf("%s: ThALES: A1=%g [deg] RMV=%g [m] RMH=%g [m] lambda=%g [Angs]\\n",  '
    )
    ThALES.append_initialize(
        "    NAME_CURRENT_COMP, ThALES_A1, ThALES_RMV, ThALES_RMH, ThALES_lambda);  "
    )
    ThALES.append_initialize("     ")
    ThALES.append_initialize("  if (ana_RMV<0) ")
    ThALES.append_initialize(
        "    ana_RMV = 1/(    (   1/dist_ana_det + 1/dist_sample_ana )   / (2*sin(DEG2RAD*fabs(A5)))    ); "
    )
    ThALES.append_initialize("  if (ana_RMH<0) ")
    ThALES.append_initialize(
        "    ana_RMH = 1/(    (   1/dist_ana_det + 1/dist_sample_ana )  *sin(DEG2RAD*fabs(A5)) / (2)    ); "
    )
    ThALES.append_initialize(
        '  printf("%s: ThALES: A5=%g [deg] ana_RMV=%g [m] ana_RMH=%g [m] lambda=%g [Angs]\\n",  '
    )
    ThALES.append_initialize(
        "    NAME_CURRENT_COMP, A5, ana_RMV, ana_RMH, ThALES_lambda);  "
    )
    ThALES.append_initialize("  A3_offset=atan(q_z_elastic/q_x_elastic)*RAD2DEG; ")

    Origin = ThALES.add_component("Origin", "Progress_bar")
    Origin.set_AT(["0", "0", "0"], RELATIVE="ABSOLUTE")

    HCS = source.HCS_source(ThALES)

    H5 = ThALES.add_component("H5", "Arm")
    H5.set_AT(["0", "0", "2.155"], RELATIVE="HCS")

    H5_rect = ThALES.add_component("H5_rect", "Guide_gravity")
    H5_rect.w1 = 0.170
    H5_rect.h1 = 0.12
    H5_rect.l = 1.0
    H5_rect.R0 = "gR0"
    H5_rect.Qc = "gQc"
    H5_rect.alpha = "gAlpha"
    H5_rect.m = 2
    H5_rect.W = "gW"
    H5_rect.set_AT(["0", "0", "2.155"], RELATIVE="HCS")

    H53_origin = ThALES.add_component("H53_origin", "Arm")
    H53_origin.set_AT(["0.17/2-0.06/2", "0", "1+gGap"], RELATIVE="H5")
    H53_origin.set_ROTATED(["0", "1.5", "0"], RELATIVE="H5")

    H53_start = ThALES.add_component("H53_start", "Arm")
    H53_start.set_AT(["0", "0", "0"], RELATIVE="H53_origin")

    H53_inpile = ThALES.copy_component("H53_inpile", "H5_rect")
    H53_inpile.w1 = 0.06
    H53_inpile.h1 = 0.12
    H53_inpile.l = "4.930-3.155"
    H53_inpile.R0 = "gR0"
    H53_inpile.Qc = "gQc"
    H53_inpile.alpha = "gAlpha"
    H53_inpile.m = 3
    H53_inpile.W = "gW"
    H53_inpile.set_AT(["0", "0", "0"], RELATIVE="H53_start")

    H53_Obt = ThALES.copy_component("H53_Obt", "H5_rect")
    H53_Obt.w1 = 0.06
    H53_Obt.h1 = 0.12
    H53_Obt.l = 3
    H53_Obt.R0 = "gR0"
    H53_Obt.Qc = "gQc"
    H53_Obt.alpha = "gAlpha"
    H53_Obt.m = 3
    H53_Obt.W = "gW"
    H53_Obt.set_AT(["0", "0", "Al_Thickness+0.015+4.930-3.155"], RELATIVE="H53_inpile")

    H53_Obt_Out = ThALES.add_component("H53_Obt_Out", "Arm")
    H53_Obt_Out.set_AT(["0", "0", "3+0.04"], RELATIVE="H53_Obt")

    H53_VSComC1 = ThALES.copy_component("H53_VSComC1", "H53_inpile")
    H53_VSComC1.w1 = 0.06
    H53_VSComC1.h1 = 0.12
    H53_VSComC1.l = 7
    H53_VSComC1.R0 = "gR0"
    H53_VSComC1.Qc = "gQc"
    H53_VSComC1.alpha = "gAlpha"
    H53_VSComC1.m = 3
    H53_VSComC1.W = "gW"
    H53_VSComC1.nelements = 7
    H53_VSComC1.set_AT(["0", "0", "3+0.075"], RELATIVE="H53_Obt")

    H53_Nose = ThALES.add_component("H53_Nose", "Guide_tapering")
    H53_Nose.option = '"parabolical"'
    H53_Nose.w1 = 0.06
    H53_Nose.h1 = 0.12
    H53_Nose.l = 2.0
    H53_Nose.linw = 0.0
    H53_Nose.loutw = 0.7
    H53_Nose.linh = 0.0
    H53_Nose.louth = 0.0
    H53_Nose.R0 = "gR0"
    H53_Nose.Qcx = "gQc"
    H53_Nose.Qcy = "gQc"
    H53_Nose.alphax = "gAlpha"
    H53_Nose.alphay = "gAlpha"
    H53_Nose.W = "gW"
    H53_Nose.mx = 3
    H53_Nose.my = 3
    H53_Nose.segno = 20
    H53_Nose.set_AT(["0", "0", "7+0.01"], RELATIVE="H53_VSComC1")

    before_monochromator_slit = ThALES.add_component(
        "before_monochromator_slit", "Slit"
    )
    before_monochromator_slit.xwidth = 0.04
    before_monochromator_slit.yheight = 0.12
    before_monochromator_slit.set_AT(["0", "0", "2.3"], RELATIVE="H53_Nose")

    l_monitor = ThALES.add_component("l_monitor", "L_monitor")
    l_monitor.nL = 200
    l_monitor.filename = '"lambda_slit_mono"'
    l_monitor.xwidth = 0.5
    l_monitor.yheight = 0.5
    l_monitor.Lmin = 0
    l_monitor.Lmax = 10
    l_monitor.restore_neutron = 1
    l_monitor.set_AT(["0", " 0", " 0.001"], RELATIVE="PREVIOUS")

    H53_ThALES_Monochromator_Cradle = ThALES.add_component(
        "H53_ThALES_Monochromator_Cradle", "Arm"
    )
    H53_ThALES_Monochromator_Cradle.set_AT(["0", "0", "2+2.3"], RELATIVE="H53_Nose")

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
    H53_ThALES_Monochromator.RV = "ThALES_RMV"
    H53_ThALES_Monochromator.RH = "ThALES_RMH"
    H53_ThALES_Monochromator.DM = "ThALES_DM"
    H53_ThALES_Monochromator.width = 0.25
    H53_ThALES_Monochromator.height = 0.2
    H53_ThALES_Monochromator.append_EXTEND("flag=SCATTERED;")
    H53_ThALES_Monochromator.set_AT(
        ["0", "0", "0"], RELATIVE="H53_ThALES_Monochromator_Cradle"
    )
    H53_ThALES_Monochromator.set_ROTATED(
        ["0", "ThALES_A1", "0"], RELATIVE="H53_ThALES_Monochromator_Cradle"
    )

    H53_ThALES_Monochromator_Out = ThALES.add_component(
        "H53_ThALES_Monochromator_Out", "Arm"
    )
    H53_ThALES_Monochromator_Out.set_AT(
        ["0", "0", "0"], RELATIVE="H53_ThALES_Monochromator_Cradle"
    )
    H53_ThALES_Monochromator_Out.set_ROTATED(
        ["0", "2*ThALES_A1", "0"], RELATIVE="H53_ThALES_Monochromator_Cradle"
    )

    before_sample_slit = ThALES.add_component("before_sample_slit", "Slit")
    before_sample_slit.xwidth = 0.03
    before_sample_slit.yheight = 0.028
    before_sample_slit.set_AT(
        ["0", "0", "ThALES_L-0.250"], RELATIVE="H53_ThALES_Monochromator_Out"
    )

    E_sample_mon = ThALES.add_component("E_sample_mon", "E_monitor")
    E_sample_mon.nE = 200
    E_sample_mon.filename = '"E_sample_mon"'
    E_sample_mon.xwidth = 0.05
    E_sample_mon.yheight = 0.05
    E_sample_mon.Emin = "Ei-ThALES_dE"
    E_sample_mon.Emax = "Ei+ThALES_dE"
    E_sample_mon.restore_neutron = 1
    E_sample_mon.set_AT(
        ["0", " 0", " ThALES_L-0.050"], RELATIVE="H53_ThALES_Monochromator_Out"
    )

    PSD_sample_mon = ThALES.add_component("PSD_sample_mon", "PSD_monitor")
    PSD_sample_mon.nx = 200
    PSD_sample_mon.ny = 200
    PSD_sample_mon.filename = '"PSD_sample_mon.dat"'
    PSD_sample_mon.xwidth = 0.05
    PSD_sample_mon.yheight = 0.05
    PSD_sample_mon.restore_neutron = 1
    PSD_sample_mon.set_AT(["0", " 0", " 0.001"], RELATIVE="PREVIOUS")

    sample_arm = ThALES.add_component("sample_arm", "Arm")
    sample_arm.set_AT(["0", "0", "ThALES_L"], RELATIVE="H53_ThALES_Monochromator_Out")

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
    res_sample.set_AT(["0", " 0", " 0"], RELATIVE="sample_arm")

    v_sample = ThALES.add_component("v_sample", "V_sample")
    v_sample.radius = 0.005
    v_sample.thickness = 0.001
    v_sample.focus_xw = 0.03
    v_sample.focus_yh = 0.04
    v_sample.yheight = 0.05
    v_sample.target_index = 3
    v_sample.set_WHEN("SAMPLE==1")
    v_sample.set_AT(["0", " 0", " 0"], RELATIVE="sample_arm")

    Sample_Out = ThALES.add_component("Sample_Out", "Arm")
    Sample_Out.set_AT(["0", "0", "0"], RELATIVE="sample_arm")
    Sample_Out.set_ROTATED(["0", " -A4", " 0"], RELATIVE="sample_arm")

    after_sample_slit = ThALES.add_component("after_sample_slit", "Slit")
    after_sample_slit.xwidth = 0.03
    after_sample_slit.yheight = 0.04
    after_sample_slit.set_AT(["0", "0", "0.250"], RELATIVE="Sample_Out")

    Ana_Cradle = ThALES.add_component("Ana_Cradle", "Arm")
    Ana_Cradle.set_AT(["0", " 0", " dist_sample_ana"], RELATIVE="Sample_Out")

    PSD_analyzer = ThALES.add_component("PSD_analyzer", "PSD_monitor")
    PSD_analyzer.nx = 200
    PSD_analyzer.ny = 200
    PSD_analyzer.filename = '"PSD_ana.dat"'
    PSD_analyzer.xwidth = 0.25
    PSD_analyzer.yheight = 0.25
    PSD_analyzer.restore_neutron = 1
    PSD_analyzer.set_AT(["0", " 0", " 0"], RELATIVE="Ana_Cradle")

    analyzer = ThALES.add_component("analyzer", "Monochromator_curved")
    analyzer.gap = 0.0005
    analyzer.NH = 11
    analyzer.NV = 9
    analyzer.mosaich = 30
    analyzer.mosaicv = 30
    analyzer.r0 = 0.7
    analyzer.RV = "ana_RMV"
    analyzer.RH = "ana_RMH"
    analyzer.DM = 3.355
    analyzer.width = 0.17
    analyzer.height = 0.13
    analyzer.set_AT(["0", " 0", " 0"], RELATIVE="Ana_Cradle")
    analyzer.set_ROTATED(["0", " -A5", " 0"], RELATIVE="Ana_Cradle")

    Ana_Out = ThALES.add_component("Ana_Out", "Arm")
    Ana_Out.set_AT(["0", "0", "0"], RELATIVE="Ana_Cradle")
    Ana_Out.set_ROTATED(["0", " -A6", " 0"], RELATIVE="Ana_Cradle")

    slit = ThALES.add_component("slit", "Slit")
    slit.xwidth = 0.03
    slit.yheight = 0.08
    slit.set_AT(["0", " 0", " 0.340"], RELATIVE="Ana_Out")

    PSD_det = ThALES.add_component("PSD_det", "PSD_monitor")
    PSD_det.nx = 200
    PSD_det.ny = 200
    PSD_det.filename = '"PSD_det.dat"'
    PSD_det.xwidth = 0.2
    PSD_det.yheight = 0.2
    PSD_det.restore_neutron = 1
    PSD_det.set_AT(["0", " 0", " dist_ana_det-0.0001"], RELATIVE="Ana_Out")

    res_monitor = ThALES.add_component("res_monitor", "Res_monitor")
    res_monitor.res_sample_comp = "res_sample"
    res_monitor.filename = '"res_monitor"'
    res_monitor.xwidth = 0.05
    res_monitor.yheight = 0.12
    res_monitor.set_WHEN("SAMPLE==0")
    res_monitor.set_AT(["0", " 0", " dist_ana_det"], RELATIVE="Ana_Out")

    detector_all = ThALES.add_component("detector_all", "Monitor")
    detector_all.xwidth = 0.05
    detector_all.yheight = 0.12
    detector_all.set_AT(["0", " 0", " dist_ana_det+0.001"], RELATIVE="Ana_Out")

    return ThALES
