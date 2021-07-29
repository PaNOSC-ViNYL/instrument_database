"""
This McStasScript file was generated from a
McStas instrument file. It is advised to check
the content to ensure it is as expected.
"""
from mcstasscript.interface import instr, plotter, functions
import sources.QuickSource_D22

def def_instrument():
    D22_quick = instr.McStas_instr("D22_quick")
    D22_quick.add_parameter("double", "lambda", value=4.5)
    D22_quick.add_parameter("double", "dlambda", value=0.45)
    D22_quick.add_parameter("double", "D22_collimation", value=2.0)
    D22_quick.add_parameter("string", "D22_sample", value="\"H2O_liq.qSq\"")
    D22_quick.add_parameter("double", "sample_size_r", value=0.005)
    D22_quick.add_parameter("double", "sample_size_y", value=0.05)
    D22_quick.add_parameter("string", "Vin_filename", value="\"none\"")
    D22_quick.add_parameter("double", "stage", value=-1.0)
    D22_quick.add_declare_var("double", "gR0", value=1.0)
    D22_quick.add_declare_var("double", "gQc", value=0.0216)
    D22_quick.add_declare_var("double", "gAlpha", value=4.07)
    D22_quick.add_declare_var("double", "gW", value=1.0/300.0)
    D22_quick.add_declare_var("double", "Al_Thickness", value=0.001)
    D22_quick.add_declare_var("double", "gGap", value=0.001)
    D22_quick.add_declare_var("double", "D22_nu", value=0.0)
    D22_quick.add_declare_var("double", "flag", value=0.0)
    D22_quick.append_initialize("  D22_nu   = 3956*48.3*DEG2RAD/2/PI/lambda/0.25; ")
    D22_quick.append_initialize("  printf(\"%s: D22:  nu=%g [rpm] lambda=%g [Angs] sample=%s\\n\", ")
    D22_quick.append_initialize("    NAME_CURRENT_COMP, D22_nu*60, lambda, D22_sample); ")

    H512_Before_VS = D22_quick.add_component("H512_Before_VS", "Progress_bar")
    H512_Before_VS.set_AT(['0', '0', '0'], RELATIVE="ABSOLUTE")

    QuickSource = sources.QuickSource_D22.D22_quick(D22_quick)
    # following are the values changed from the default and that make sense only in this instrument description
    QuickSource.focus_xw = "sample_size_r"
    QuickSource.focus_yh = "sample_size_y"
    QuickSource.lambda0 = "lambda"
    QuickSource.dlambda = "dlambda"
    QuickSource.target_index = +5
    QuickSource.set_AT(['0', ' 0', ' 0'], RELATIVE="H512_Before_VS")
    
    D22_Sample_Pos = D22_quick.add_component("D22_Sample_Pos", "Arm")
    D22_Sample_Pos.set_AT(['0', '0', '20+0.3+0.3'], RELATIVE="H512_Before_VS")
    
    H51_D22_Sample_Div = D22_quick.add_component("H51_D22_Sample_Div", "Monitor_nD")
    H51_D22_Sample_Div.xwidth = 0.02
    H51_D22_Sample_Div.yheight = 0.05
    H51_D22_Sample_Div.bins = 100
    H51_D22_Sample_Div.restore_neutron = 1
    H51_D22_Sample_Div.options = "\"dx limits=[-2 2], dy limits=[-2 2]\""
    H51_D22_Sample_Div.set_AT(['0', '0', '0'], RELATIVE="D22_Sample_Pos")
    
    H51_D22_Sample_XY = D22_quick.add_component("H51_D22_Sample_XY", "Monitor_nD")
    H51_D22_Sample_XY.xwidth = 0.02
    H51_D22_Sample_XY.yheight = 0.05
    H51_D22_Sample_XY.bins = 50
    H51_D22_Sample_XY.restore_neutron = 1
    H51_D22_Sample_XY.options = "\"x y\""
    H51_D22_Sample_XY.set_AT(['0', '0', '0'], RELATIVE="D22_Sample_Pos")
    
    H51_D22_Sample_L = D22_quick.add_component("H51_D22_Sample_L", "Monitor_nD")
    H51_D22_Sample_L.xwidth = 0.02
    H51_D22_Sample_L.yheight = 0.05
    H51_D22_Sample_L.bins = 50
    H51_D22_Sample_L.restore_neutron = 1
    H51_D22_Sample_L.options = "\"lambda limits=[1 10]\""
    H51_D22_Sample_L.set_AT(['0', '0', '0'], RELATIVE="D22_Sample_Pos")
    
    H51_D22_Sample = D22_quick.add_component("H51_D22_Sample", "Isotropic_Sqw")
    H51_D22_Sample.Sqw_coh = "D22_sample"
    H51_D22_Sample.Sqw_inc = "NULL"
    H51_D22_Sample.radius = "sample_size_r"
    H51_D22_Sample.yheight = "sample_size_y"
    H51_D22_Sample.d_phi = "RAD2DEG*atan2(1, D22_collimation)"
    H51_D22_Sample.append_EXTEND("if (!SCATTERED) ABSORB;")
    H51_D22_Sample.set_AT(['0', '0', '0'], RELATIVE="D22_Sample_Pos")
   
    Slit = D22_quick.add_component("Slit", "Slit")
    Slit.xmin = -0.6
    Slit.xmax = 0.6
    Slit.ymin = -0.6
    Slit.ymax = 0.6
    Slit.set_AT(['0', '0', 'D22_collimation-0.011'], RELATIVE="D22_Sample_Pos")
    
    Vout = D22_quick.add_component("Vout", "MCPL_output")
    Vout.filename = "\"sDETECTOR\""
    Vout.set_AT(['0', ' 0', ' D22_collimation-0.01'], RELATIVE="PREVIOUS")
    
    Detector = D22_quick.add_component("Detector", "Monitor_nD")
    Detector.xwidth = 0.98
    Detector.yheight = 1.024
    Detector.options = "\"x bins=128 y bins=256\""
    Detector.append_EXTEND("ABSORB; ")
    Detector.set_AT(['0', '0', 'D22_collimation'], RELATIVE="D22_Sample_Pos")
    
    return D22_quick

#D22_quick = D22_quick_fun()
