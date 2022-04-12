import compare_data

from mcstasscript.interface import instr, plotter, functions

def def_instrument():
    DMC = instr.McStas_instr("DMC_model")
    
    DMC.add_declare_var("double", "mono_q", value=1.8734)
    DMC.add_declare_var("double", "OMA")
    DMC.add_declare_var("double", "RV")
    DMC.add_declare_var("double", "y_mono", value=0.025)
    DMC.add_declare_var("double", "NV", value=5.0)
    DMC.add_declare_var("double", "d_phi_0")
    DMC.add_declare_var("double", "TTM")
    DMC.add_declare_var("double", "sample_radius", value=0.008/2)
    DMC.add_declare_var("double", "sample_height", value=0.03)
    DMC.add_declare_var("double", "can_radius", value=0.0083/2)
    DMC.add_declare_var("double", "can_height", value=0.0303)
    DMC.add_declare_var("double", "can_thick", value=0.00015)
    DMC.add_declare_var("double", "alpha")
    DMC.add_declare_var("double", "Qc", value=0.0217)
    DMC.add_declare_var("double", "R0", value=0.995)
    DMC.add_declare_var("double", "Mvalue", value=1.9)
    DMC.add_declare_var("double", "W", value=1.0/250.0)
    DMC.add_declare_var("double", "alpha_curve")
    DMC.add_declare_var("double", "Qc_curve", value=0.0217)
    DMC.add_declare_var("double", "R0_curve", value=0.995)
    DMC.add_declare_var("double", "Mvalue_curve", value=2.1)
    DMC.add_declare_var("double", "W_curve", value=1.0/250.0)
    DMC.add_declare_var("double", "ldiff", value=0.05)
    DMC.add_declare_var("double", "angleGuideCurved")
    
    DMC.add_parameter("double", "R", value=0.87)
    DMC.add_parameter("double", "R_curve", value=0.87)
    DMC.append_initialize("angleGuideCurved=-2.0*asin(0.4995 /2.0/3612)/PI*180;")
    DMC.append_initialize("alpha=(R0-R)/Qc/(Mvalue-1);")
    DMC.append_initialize("alpha_curve=(R0_curve-R_curve)/Qc_curve/(Mvalue_curve-1);")
    DMC.append_initialize("")

    source_arm = DMC.add_component("source_arm", "Progress_bar")
    source_arm.set_AT([0, 0, 0])

    DMC.add_parameter("double", "lambda", value=2.5666)

    source = DMC.add_component("source", "Source_Maxwell_3")
    source.yheight = 0.156
    source.xwidth = 0.126
    source.Lmin = "lambda-ldiff/2"
    source.Lmax = "lambda+ldiff/2"
    source.dist = 1.5
    source.focus_xw = 0.02
    source.focus_yh = 0.12
    source.T1 = 296.16
    source.T2 = 40.68
    source.I1 = 8.5E11
    source.I2 = 5.2E11
    source.set_AT([0, 0, 0], RELATIVE=source_arm)

    PSDbefore_guides = DMC.add_component("PSDbefore_guides", "PSD_monitor")
    PSDbefore_guides.nx = 128
    PSDbefore_guides.ny = 128
    PSDbefore_guides.filename = '"PSDbefore_guides"'
    PSDbefore_guides.xwidth = 0.02
    PSDbefore_guides.yheight = 0.12
    PSDbefore_guides.set_AT([0, 0, 1.49999], RELATIVE=source_arm)

    l_mon_source = DMC.add_component("l_mon_source", "L_monitor")
    l_mon_source.nL = 101
    l_mon_source.filename = '"lmonsource.dat"'
    l_mon_source.xwidth = 0.02
    l_mon_source.yheight = 0.12
    l_mon_source.Lmin = 0
    l_mon_source.Lmax = 20
    l_mon_source.set_AT([0, 0, 1e-9], RELATIVE="PREVIOUS")

    guide1 = DMC.add_component("guide1", "Guide")
    guide1.w1 = 0.02
    guide1.h1 = 0.12
    guide1.w2 = 0.02
    guide1.h2 = 0.12
    guide1.l = 4.66
    guide1.R0 = "R0"
    guide1.Qc = "Qc"
    guide1.alpha = "alpha"
    guide1.m = 1.8
    guide1.W = "W"
    guide1.set_AT([0, 0, 1.50], RELATIVE=source_arm)
    guide1.set_ROTATED([0, 0, 0], RELATIVE=source_arm)

    PSDbefore_curve = DMC.add_component("PSDbefore_curve", "PSD_monitor")
    PSDbefore_curve.nx = 128
    PSDbefore_curve.ny = 128
    PSDbefore_curve.filename = '"PSDbefore_curve"'
    PSDbefore_curve.xwidth = 0.02
    PSDbefore_curve.yheight = 0.12
    PSDbefore_curve.set_AT([0, 0, 4.664], RELATIVE=guide1)

    guide2 = DMC.add_component("guide2", "Bender")
    guide2.w = 0.02
    guide2.h = 0.12
    guide2.r = 3612
    guide2.l = 20
    guide2.R0a = "R0_curve"
    guide2.Qca = "Qc_curve"
    guide2.alphaa = "alpha_curve"
    guide2.ma = "Mvalue_curve"
    guide2.Wa = "W_curve"
    guide2.R0i = "R0_curve"
    guide2.Qci = "Qc_curve"
    guide2.alphai = "alpha_curve"
    guide2.mi = 1
    guide2.Wi = "W_curve"
    guide2.R0s = "R0_curve"
    guide2.Qcs = "Qc_curve"
    guide2.alphas = "alpha_curve"
    guide2.ms = "Mvalue_curve"
    guide2.Ws = "W_curve"
    guide2.set_AT([0, 0, 4.69], RELATIVE=guide1)

    PSDafter_curve = DMC.add_component("PSDafter_curve", "PSD_monitor")
    PSDafter_curve.nx = 128
    PSDafter_curve.ny = 128
    PSDafter_curve.filename = '"PSDafter_curve"'
    PSDafter_curve.xwidth = 0.02
    PSDafter_curve.yheight = 0.12
    PSDafter_curve.set_AT([0, 0, 20.0001], RELATIVE=guide2)

    bunker = DMC.add_component("bunker", "Guide")
    bunker.w1 = 0.02
    bunker.h1 = .12
    bunker.w2 = 0.02
    bunker.h2 = .12
    bunker.l = 3.43
    bunker.R0 = "R0"
    bunker.Qc = "Qc"
    bunker.alpha = "alpha"
    bunker.m = 1.6
    bunker.W = "W"
    bunker.set_AT([0, 0, 20.1502], RELATIVE=guide2)
    bunker.set_ROTATED([0, 0, 0], RELATIVE=guide2)

    guide3 = DMC.add_component("guide3", "Guide")
    guide3.w1 = 0.02
    guide3.h1 = .12
    guide3.w2 = 0.02
    guide3.h2 = .12
    guide3.l = 12.275
    guide3.R0 = "R0"
    guide3.Qc = "Qc"
    guide3.alpha = "alpha"
    guide3.m = 1.6
    guide3.W = "W"
    guide3.set_AT([0, 0, 3.56], RELATIVE=bunker)

    guide4 = DMC.add_component("guide4", "Guide")
    guide4.w1 = 0.02
    guide4.h1 = .12
    guide4.w2 = 0.02
    guide4.h2 = .12
    guide4.l = 5.66
    guide4.R0 = "R0"
    guide4.Qc = "Qc"
    guide4.alpha = "alpha"
    guide4.m = 1.6
    guide4.W = "W"
    guide4.set_AT([0, 0, 15.8555], RELATIVE=bunker)
    guide4.set_ROTATED([0, 0, 0], RELATIVE=guide3)

    window1 = DMC.add_component("window1", "Al_window")
    window1.thickness = 0.002
    window1.set_AT([0, 0, "5.66+1e-9"], RELATIVE="PREVIOUS")

    ydist_fluxpos = DMC.add_component("ydist_fluxpos", "PSDlin_monitor")
    ydist_fluxpos.nx = 11
    ydist_fluxpos.filename = '"ydist_fluxpos.dat"'
    ydist_fluxpos.xwidth = 0.120
    ydist_fluxpos.yheight = 0.02
    ydist_fluxpos.set_AT([0, 0, "5.66+1e-8+0.01"], RELATIVE=guide4)
    ydist_fluxpos.set_ROTATED([0, 0, '90'], RELATIVE="PREVIOUS")

    PSD_fluxpos = DMC.add_component("PSD_fluxpos", "PSD_monitor")
    PSD_fluxpos.nx = 100
    PSD_fluxpos.ny = 100
    PSD_fluxpos.filename = '"xdist_fluxposy.dat"'
    PSD_fluxpos.xwidth = 0.02
    PSD_fluxpos.yheight = 0.12
    PSD_fluxpos.set_AT([0, 0, "5.66+1e-7+0.01"], RELATIVE=guide4)

    xdist_flux_pos = DMC.add_component("xdist_flux_pos", "PSDlin_monitor")
    xdist_flux_pos.nx = 11
    xdist_flux_pos.filename = '"xdist_fluxpos.dat"'
    xdist_flux_pos.xwidth = 0.020
    xdist_flux_pos.yheight = 0.12
    xdist_flux_pos.set_AT([0, 0, 1e-9], RELATIVE="PREVIOUS")

    PSD_fluxposB = DMC.add_component("PSD_fluxposB", "PSD_monitor")
    PSD_fluxposB.nx = 100
    PSD_fluxposB.ny = 100
    PSD_fluxposB.filename = '"PSD_fluxposB.dat"'
    PSD_fluxposB.xwidth = 0.02
    PSD_fluxposB.yheight = 0.12
    PSD_fluxposB.set_AT([0, 0, "6.24-1e-7-0.01"], RELATIVE=guide4)

    window2 = DMC.add_component("window2", "Al_window")
    window2.thickness = 0.002
    window2.set_AT([0, 0, 1e-9], RELATIVE="PREVIOUS")

    in_slit = DMC.add_component("in_slit", "Slit")
    in_slit.xmin = -0.01
    in_slit.xmax = 0.01
    in_slit.ymin = -0.06
    in_slit.ymax = 0.06
    in_slit.set_AT([0, 0, 0.0021], RELATIVE="PREVIOUS")

    lambda_in = DMC.add_component("lambda_in", "L_monitor")
    lambda_in.nL = 128
    lambda_in.filename = '"L_in.dat"'
    lambda_in.xmin = -0.011
    lambda_in.xmax = 0.011
    lambda_in.ymin = -0.061
    lambda_in.ymax = 0.061
    lambda_in.Lmin = 0
    lambda_in.Lmax = "2*lambda"
    lambda_in.set_AT([0, 0, 0.001], RELATIVE=in_slit)

    sma = DMC.add_component("sma", "Arm")
    sma.set_AT([0, 0, 0.65], RELATIVE=in_slit)
    sma.set_ROTATED([0, 'OMA', 0], RELATIVE=in_slit)
    
    DMC.append_initialize("TTM = 2*asin(mono_q*lambda/(4*PI))*RAD2DEG;")
    DMC.append_initialize("OMA = TTM/2;")
    DMC.append_initialize("RV = fabs(2*2.82*sin(DEG2RAD*OMA));")

    foc_mono = DMC.add_component("foc_mono", "Monochromator_2foc")
    foc_mono.zwidth = 0.05
    foc_mono.yheight = 0.025
    foc_mono.gap = 0.0005
    foc_mono.NH = 1
    foc_mono.NV = 5
    foc_mono.mosaich = 38
    foc_mono.mosaicv = 38
    foc_mono.r0 = 0.7
    foc_mono.Q = "mono_q"
    foc_mono.RV = "RV"
    foc_mono.RH = 0
    foc_mono.set_SPLIT(10)
    foc_mono.set_AT([0, 0, 0], RELATIVE=sma)
    
    

    msa = DMC.add_component("msa", "Arm")
    msa.set_AT([0, 0, 0], RELATIVE=sma)
    msa.set_ROTATED([0, ' TTM', 0], RELATIVE=in_slit)

    out1_slit = DMC.add_component("out1_slit", "Slit")
    out1_slit.xmin = -0.01
    out1_slit.xmax = 0.01
    out1_slit.ymin = -0.06
    out1_slit.ymax = 0.06
    out1_slit.set_AT([0, 0, 0.2], RELATIVE=msa)
    out1_slit.set_ROTATED([0, 0, 0], RELATIVE=msa)

    Amoin_slit = DMC.add_component("Amoin_slit", "Slit")
    Amoin_slit.xmin = -0.01
    Amoin_slit.xmax = 0.01
    Amoin_slit.ymin = -0.06
    Amoin_slit.ymax = 0.06
    Amoin_slit.set_AT([0, 0, 0.325], RELATIVE=msa)
    Amoin_slit.set_ROTATED([0, 0, 0], RELATIVE=msa)

    Bmoin_slit = DMC.add_component("Bmoin_slit", "Slit")
    Bmoin_slit.xmin = -0.01
    Bmoin_slit.xmax = 0.01
    Bmoin_slit.ymin = -0.06
    Bmoin_slit.ymax = 0.06
    Bmoin_slit.set_AT([0, 0, 0.525], RELATIVE=msa)
    Bmoin_slit.set_ROTATED([0, 0, 0], RELATIVE=msa)

    out2_slit = DMC.add_component("out2_slit", "Slit")
    out2_slit.xmin = -0.01
    out2_slit.xmax = 0.01
    out2_slit.ymin = -0.06
    out2_slit.ymax = 0.06
    out2_slit.set_AT([0, 0, 0.65], RELATIVE=msa)
    out2_slit.set_ROTATED([0, 0, 0], RELATIVE=msa)

    PSD_sample = DMC.add_component("PSD_sample", "PSD_monitor")
    PSD_sample.nx = 80
    PSD_sample.ny = 80
    PSD_sample.filename = '"PSD_sample.dat"'
    PSD_sample.xmin = -0.05
    PSD_sample.xmax = 0.05
    PSD_sample.ymin = -0.07
    PSD_sample.ymax = 0.07
    PSD_sample.set_AT([0, 0, 2.77], RELATIVE=msa)

    lambda_sample = DMC.add_component("lambda_sample", "L_monitor")
    lambda_sample.nL = 128
    lambda_sample.filename = '"L_sample.dat"'
    lambda_sample.xmin = "-sample_radius"
    lambda_sample.xmax = "sample_radius"
    lambda_sample.ymin = "-sample_height/2"
    lambda_sample.ymax = "sample_height/2"
    lambda_sample.Lmin = "lambda-0.2"
    lambda_sample.Lmax = "lambda+0.2"
    lambda_sample.set_AT([0, 0, 2.81], RELATIVE=msa)

    sa_arm = DMC.add_component("sa_arm", "Arm")
    sa_arm.set_AT([0, 0, ' 2.82'], RELATIVE=msa)

    sample = DMC.add_component("sample", "PowderN")
    sample.reflections = DMC.add_parameter("string", "filename", value='"Na2Ca3Al2F14.laz"')
    sample.radius = "sample_radius"
    sample.yheight = "sample_height"
    sample.pack = DMC.add_parameter("double", "PACK", value=0.7)
    sample.p_inc = 0
    sample.p_transmit = 0
    sample.DW = DMC.add_parameter("double", "Dw", value=0.8)
    sample.d_phi = DMC.add_parameter("double", "D_PHI", value=6.0)
    sample.barns = DMC.add_parameter("double", "BARNS", value=1.0)
    sample.set_SPLIT(10)
    sample.set_AT([0, 0, 0], RELATIVE=sa_arm)

    STOP = DMC.add_component("STOP", "Beamstop")
    STOP.radius = 0.3
    STOP.set_AT([0, 0, 1.4], RELATIVE=sa_arm)
    STOP.set_ROTATED([0, 0, 0], RELATIVE=sa_arm)

    DMC.add_parameter("double", "SHIFT", value=0.0)

    Detector = DMC.add_component("Detector", "Monitor_nD")
    Detector.xwidth = 3.0
    Detector.yheight = 0.09
    Detector.bins = 400
    Detector.min = "19.9+SHIFT"
    Detector.max = "99.9+SHIFT"
    Detector.options = '"banana,theta"'
    Detector.filename = '"detector.dat"'
    Detector.set_AT([0, 0, 0], RELATIVE=sa_arm)
    Detector.set_ROTATED([0, 0, 180], RELATIVE=sa_arm)

    return DMC
