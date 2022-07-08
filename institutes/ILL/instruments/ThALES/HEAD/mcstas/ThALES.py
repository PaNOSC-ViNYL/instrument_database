"""

 - The monochromator do not allow transmitted neutrons

TODO: 
 - [ ] Check RV RH for monochromator and analyzer
 - [ ] Understand the A3_offset:
        ThALES.add_parameter("double", "q_x_elastic", value=1.3139)
        ThALES.add_parameter("double", "q_z_elastic", value=0.146)
        ThALES.append_initialize("  A3_offset=atan(q_z_elastic/q_x_elastic)*RAD2DEG; ")
 - [X] Aspetto che Martin controlli i numeri dal mio grafichetto
 - [X] Chiedere a Martin i sample environments

"""
import os

# MCSTAS_PATH = os.environ["MCSTAS"]
from mcstasscript.interface import functions
import mcstasscript as ms

my_configurator = functions.Configurator()

# my_configurator.set_mcstas_path(MCSTAS_PATH)
# my_configurator.set_mcrun_path(MCSTAS_PATH + "/bin/")

# list here all the common parts to be imported

# from institutes.ILL.sources.HEAD.mcstas import Full as source
from institutes.ILL.sources.HEAD.mcstas import Gauss as source
from libpyvinyl.Instrument import Instrument
from mcstasscript.interface import instr
import pint

# from libpyvinyl import ureg
from pint import set_application_registry

ureg = pint.get_application_registry()

# for conversion from degree to radians and viceversa
import math


from libpyvinyl.Parameters import Parameter


## Global variables for DEBUG to be removed!
addMonitors = True
addCryostat = True


def addMonitor(
    instr,
    name: str,
    position=[0, 0, 0],
    width=0.10,
    height=0.15,
    bins=100,
    radius=0.30,
    energy=5,
    denergy=0.6,
):
    if addMonitors is False:
        return

    monitor = instr.add_component(name + "_DEBUG", "Monitor_nD")
    monitor.filename = '"' + name + '_DEBUG"'
    monitor.xwidth = width
    monitor.yheight = height
    monitor.bins = bins
    monitor.options = '"energy limits=[{0:.2f} {1:.2f}] x y, parallel"'.format(
        energy - denergy, energy + denergy
    )

    # monitor.options = '"energy limits=[Ei-dE, Ei+dE] x y, borders, parallel"'
    monitor.set_AT(position, RELATIVE="PREVIOUS")

    psd = instr.add_component(name + "_psd_DEBUG", "Monitor_nD")
    psd.filename = '"' + name + '_DEBUG"'
    psd.xwidth = width
    psd.yheight = height
    psd.bins = bins
    psd.options = '"x y, parallel"'
    psd.set_AT([0, 0, 0], RELATIVE="PREVIOUS")

    psd = instr.add_component(name + "_psdcyl_DEBUG", "Monitor_nD")
    psd.filename = '"' + name + '_DEBUG"'
    psd.xwidth = 2 * radius
    psd.yheight = height
    psd.options = '"theta bins=360, y bins=20, cylinder, parallel"'
    psd.set_AT([0, 0, 0], RELATIVE="PREVIOUS")

    return psd


def cryoOrange(
    instrument,
    name="cryoOrange",
    position=[0, 0, 0],
    relative="PREVIOUS",
    number_of_activations=1,
):

    if addCryostat is False:
        mycryo = instrument.add_component(name, "Arm")
        mycryo.set_AT(position, RELATIVE=relative)
        return mycryo, None

    mycryo = ms.Cryostat(name, instrument)
    mycryo.set_AT(position, RELATIVE=relative)

    mycryo.add_layer(
        inner_radius=150e-3 / 2.0,  # diameter: 150 mm
        thickness=1.0e-3,  # 1 mm
        origin_to_top=150e-3 / 2,  # guessed
        top_thickness=3e-3,
        origin_to_bottom=150e-3 / 2,  # guessed
        bottom_thickness=3e-3,
        p_interact=0.2,
    )

    #    mycryo.last_layer.add_window(
    #        inner_radius=mycryo.last_layer.inner_radius,
    #        thickness=1.0 / 1000,
    #        origin_to_top=11.0 / 1000,  # guessed
    #        origin_to_bottom=11.0 / 1000,  # guessed
    #    )

    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_layer(
        inner_radius=159e-3 / 2.0,  # diameter: 150 mm
        thickness=0.5e-3,  # 1 mm
        origin_to_top=150e-3 / 2,  # guessed
        top_thickness=3e-3,
        origin_to_bottom=150e-3 / 2,  # guessed
        bottom_thickness=3e-3,
        p_interact=0.2,
    )
    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_layer(
        inner_radius=170e-3 / 2.0,  # diameter: 150 mm
        thickness=1e-3,  # 1 mm
        origin_to_top=150e-3 / 2,  # guessed
        top_thickness=3e-3,
        origin_to_bottom=150e-3 / 2,  # guessed
        bottom_thickness=3e-3,
        p_interact=0.2,
    )
    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_layer(
        inner_radius=280e-3 / 2.0,  # diameter: 150 mm
        thickness=1e-3,  # 1 mm
        origin_to_top=150e-3 / 2,  # guessed
        top_thickness=3e-3,
        origin_to_bottom=150e-3 / 2,  # guessed
        bottom_thickness=3e-3,
        p_interact=0.2,
    )
    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_spatial_loggers()
    mycryo.build(include_master=False)

    exit = instrument.add_component("exit_volume", "Union_cylinder")
    exit.set_AT([0, 0, 0], RELATIVE=name)

    master_before_sample = instrument.add_component(
        "master_before_sample", "Union_master", RELATIVE=name
    )

    return mycryo, exit


def cryo10T(
    instrument,
    name="10T",
    position=[0, 0, 0],
    relative="PREVIOUS",
    number_of_activations=1,
):
    """Create the cryostat at the sample_arm position previously defined"""
    if addCryostat is False:
        mycryo = instrument.add_component(name, "Arm")
        mycryo.set_AT(position, RELATIVE=relative)
        return mycryo, None

    mycryo = ms.Cryostat(name, instrument)
    mycryo.set_AT(position, RELATIVE=relative)

    mycryo.add_layer(
        inner_radius=50.8 / 2.0 / 1000,  # diameter: 50.8mm
        thickness=2.0 / 1000,  # 1 mm
        origin_to_top=150.0 / 1000,  # guessed
        top_thickness=0.001,
        origin_to_bottom=150.0 / 1000,  # guessed
        bottom_thickness=0.001,
        p_interact=0.2,
    )

    mycryo.last_layer.add_window(
        inner_radius=mycryo.last_layer.inner_radius,
        thickness=1.0 / 1000,
        origin_to_top=11.0 / 1000,  # guessed
        origin_to_bottom=11.0 / 1000,  # guessed
    )

    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.add_spatial_loggers()
    mycryo.build(include_master=False)

    exit = instrument.add_component("exit_volume", "Union_cylinder")
    exit.set_AT([0, 0, 0], RELATIVE=name)

    master_before_sample = instrument.add_component(
        "master_before_sample", "Union_master", RELATIVE=name
    )

    return mycryo, exit

    ################### ignored for the moment
    mycryo.add_layer(
        inner_radius=592.0 / 2.0 / 1000,  # diameter: 592 mm
        thickness=3.0 / 1000,  # 3 mm
        origin_to_top=140.0 / 1000,  # guessed
        top_thickness=0.001,
        origin_to_bottom=140.0 / 1000,  # guessed
        bottom_thickness=0.001,
        p_interact=0.2,
        material="Al",
    )
    mycryo.last_layer.add_window(
        inner_radius=mycryo.last_layer.inner_radius,
        thickness=1.0 / 1000,
        origin_to_top=51.0 / 1000,  # guessed
        origin_to_bottom=51.0 / 1000,  # guessed
    )
    #    mycryo.add_spatial_loggers()
    for component in mycryo.last_layer.union_components:
        component.number_of_activations = number_of_activations

    mycryo.build(include_master=False)

    return mycryo


class BraggAngle(Parameter):
    """
    A parameter describing an angle and associating it with wavelength or
    energy satisfying the bragg law
    """

    def __init__(self, *args, **kwargs):
        """ """
        super().__init__(*args, **kwargs)
        __d_lattice = pint.Quantity("0 angstroms")

    @property
    def d_lattice(self):
        return self.__d_lattice

    @d_lattice.setter
    def d_lattice(self, value) -> None:
        self.__d_lattice = value

    @property
    def wavelength(self):
        return 2 * self.d_lattice * math.sin(self.pint_value.to("radians") / 2.0)

    @wavelength.setter
    def wavelength(self, value):
        self.value = math.asin(value / 2.0 / self.d_lattice) * 2 * ureg.radians

    @property
    def energy(self):
        l = self.wavelength
        return (81.80421 * ureg.meV * ureg.angstrom * ureg.angstrom) / (l * l)

    @energy.setter
    def energy(self, value):
        self.wavelength = math.sqrt(81.80421 * ureg.meV / value) * ureg.angstrom


import mcstasscript as ms


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

    ThALES_L = 2.000  # distance between monochromator and sample
    dist_sample_ana = 1.260  # distance between sample and analyzer
    dist_ana_det = 0.640  # distance between analyzer and detector

    ThALES = instr.McStas_instr("ThALES")
    myinstr.add_calculator(ThALES)

    a2 = BraggAngle(
        "a2",
        comment="Angle between beam reflected by monochromator and incident beam",
        unit="degree",
    )
    a2.type = "double"
    a2.value = 33
    ThALES.parameters.add(a2)
    #    a2 = ThALES.add_parameter(
    #        "double",
    #        "a2",
    #        comment="Angle between beam reflected by monochromator and incident beam",
    #        unit="degree",
    #        value=33,
    #    )
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

    a6 = BraggAngle(
        "a6",
        comment="Angle between reflected by analyzer and incident beams",
        unit="degree",
    )
    a6.type = "double"
    a6.value = 33
    ThALES.parameters.add(a6)

    monochromator_d = ThALES.add_parameter(
        "double",
        "monochromator_d",
        comment="Monochromator lattice parameter",
        value=3.355,
        unit="angstrom",
    )
    a2.d_lattice = monochromator_d.pint_value
    a6.d_lattice = monochromator_d.pint_value

    #   ThALES.add_parameter("double", "q_x_elastic", value=1.3139)
    #   ThALES.add_parameter("double", "q_z_elastic", value=0.146)
    ThALES.add_parameter("double", "SAMPLE", value=1)

    #    ThALES.add_declare_var("double", "L")
    ThALES.add_declare_var("double", "flag")
    #    ThALES.add_declare_var("double", "A5")
    #    ThALES.add_declare_var("double", "A6")
    #    ThALES.add_declare_var("double", "final_lambda")
    ThALES.add_declare_var("double", "sample_select")

    #   ThALES.add_declare_var("double", "A3_offset")
    #   ThALES.add_declare_var("int", "elastic_flag_instr_1")
    #   ThALES.add_declare_var("int", "elastic_flag_instr_2")

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
    HCS.E0 = "Ei"
    HCS.target_index = 1

    # override the value of lambda by the value of the angle
    ThALES.append_initialize("lambda =2*sin(0.5*a2*DEG2RAD)*monochromator_d;")
    ThALES.append_initialize('printf("lambda: %.2f\\n", lambda);')
    ThALES.append_initialize("dlambda=0.01*lambda;")
    # the following python string is an alternative way of doing:
    ThALES.add_declare_var("double", "energy")
    ThALES.append_initialize("Ei = 81.80421/(lambda*lambda);")
    ThALES.add_declare_var("double", "denergy")
    ThALES.append_initialize("dE = 2*Ei*dlambda/lambda;")
    ThALES.append_initialize('printf("lambda: %.2f +/- %.2f\\n", lambda, dlambda);')
    ThALES.append_initialize('printf("energy: %.2f +/- %.2f\\n", Ei, dE);')
    # lambda2energy = "81.80421/(" + wavelength.name + "*" + wavelength.name + ")"

    # start of the guide
    H5 = ThALES.add_component("H5", "Arm")
    H5.set_AT([0, 0, 2.155], RELATIVE=HCS)

    # adds monitors at the same position of the previous component
    addMonitor(ThALES, "H5")

    H53_1a = ThALES.add_component("H5_1", "Guide_gravity")
    # material: Al 5083
    H53_1a.w1 = 0.060
    H53_1a.h1 = 0.120
    H53_1a.l = 1.0
    #    H53_1a.R0 = gR0
    #    H53_1a.Qc = gQc
    #    H53_1a.alpha = gAlpha
    #    H53_1a.m = 3
    #    H53_1a.W = gW
    H53_1a.reflect = '"supermirror_m3.rfl"'
    H53_1a.set_AT([0, 0, 0], RELATIVE=H5)

    addMonitor(ThALES, "H5_1a", H53_1a.l)
    #    H53_origin = ThALES.add_component("H53_origin", "Arm")
    # I don't understand the x
    #    H53_origin.set_AT([H5_rect.w1 / 2 - 0.06 / 2, 0, H5_rect.l + gGap], RELATIVE=H5)

    #    H53_origin.set_ROTATED([0, 1.5, 0], RELATIVE=H5)

    #   H53_start = ThALES.add_component("H53_start", "Arm")
    #   H53_start.set_AT([0, 0, 0], RELATIVE=H53_origin)

    H53_1b = ThALES.copy_component("H5_1b", H53_1a)
    H53_1b.l = 1.775
    H53_1b.set_AT([0, 0, H53_1a.l], RELATIVE=H53_1a)

    addMonitor(ThALES, "H5_1b", [0, 0, H53_1b.l])

    # membrane
    # gap 25 mm
    H53_A = ThALES.add_component("H53_A", "Arm")
    H53_A.set_AT([0, 0, H53_1b.l + 0.025 + 0.006], RELATIVE=H53_1b)

    H53_2a = ThALES.copy_component("H53_2a", H53_1a)
    H53_2a.l = 2.317
    H53_2a.set_AT([0, 0, 0], RELATIVE=H53_A)

    H53_2b = ThALES.copy_component("H53_2b", H53_1a)
    H53_2b.l = 0.757
    H53_2b.set_AT([0, 0, H53_2a.l], RELATIVE=H53_2a)

    H53_B = ThALES.add_component("H53_B", "Arm")
    H53_B.set_AT([0, 0, H53_2b.l + 0.060 + 0.027], RELATIVE=H53_2b)

    H53_3 = ThALES.copy_component("H53_3", H53_1a)
    H53_3.l = 0.501
    H53_3.set_AT([0, 0, 0], RELATIVE=H53_B)

    H53_4 = ThALES.copy_component("H53_4", H53_1a)
    H53_4.l = 0.501
    H53_4.set_AT([0, 0, H53_3.l], RELATIVE=H53_3)

    H53_5 = ThALES.copy_component("H53_5", H53_1a)
    H53_5.l = 5.809
    H53_5.set_AT([0, 0, H53_4.l], RELATIVE=H53_4)

    H53_C = ThALES.add_component("H53_C", "Arm")
    H53_C.set_AT([0, 0, H53_5.l + 0.0078 + 0.024 + 0.110], RELATIVE=H53_5)

    H53_6 = ThALES.copy_component("H53_6", H53_1a)
    H53_6.l = 2.420
    H53_6.set_AT([0, 0, 0], RELATIVE=H53_C)

    H53_D = ThALES.add_component("H53_D", "Arm")
    H53_D.set_AT([0, 0, H53_6.l + 0.0195 + 0.002], RELATIVE=H53_6)

    # adds monitors at the same position of the previous component
    addMonitor(ThALES, "H53_D")

    H53_7 = ThALES.add_component("H53_7", "Guide_tapering")
    H53_7.option = '"elliptical"'
    H53_7.w1 = 0.06
    H53_7.h1 = 0.12
    H53_7.l = 2.4915
    H53_7_ah = 2.737445  # major semi-axis
    H53_7_bh = H53_7.w1 / 2.0  # minor semi-axis
    H53_7.linw = math.sqrt(H53_7_ah * H53_7_ah - H53_7_bh * H53_7_bh)
    H53_7.loutw = (
        H53_7.linw - H53_7.l
    )  # the length of the guide must be smaller than the coordinate of the focus
    #    H53_7_linw = 0
    #    H53_7_linh = 0
    H53_7.linh = 0  # plane mirrors (vertical)
    H53_7.louth = 0  #  plane mirrors (vertical)
    H53_7.R0 = gR0
    H53_7.Qcx = gQc
    H53_7.Qcy = gQc
    H53_7.alphax = gAlpha
    H53_7.alphay = gAlpha
    H53_7.W = gW
    H53_7.mx = 3
    H53_7.my = 3
    H53_7.segno = 5  # is this the number of segments?
    H53_7.set_AT([0, 0, 0], RELATIVE=H53_D)

    H53_7out = ThALES.add_component("H53_7out", "Arm")
    H53_7out.set_AT([0, 0, H53_7.l], RELATIVE=H53_7)

    # adds monitors at the same position of the previous component
    addMonitor(ThALES, "H53_7")

    slit_A = ThALES.add_component("slit_A", "Slit")
    slit_A.xwidth = 0.04
    slit_A.yheight = 0.12
    slit_A.set_AT([0, 0, H53_7.l + 0.3], RELATIVE=H53_7)

    # adds monitors at the same position of the previous component
    addMonitor(ThALES, "slit_A")

    Monochromator_Arm = ThALES.add_component("Monochromator_Arm", "Arm")
    Monochromator_Arm.set_AT([0, 0, 2], RELATIVE=slit_A)

    addMonitor(ThALES, "mono_in")

    # double check the probability of reflection
    Monochromator = ThALES.add_component("Monochromator", "Monochromator_curved")
    #   Monochromator.reflect = '"HOPG.rfl"'
    #   Monochromator.transmit = '"HOPG.trm"'
    Monochromator.gap = 0.0005
    Monochromator.NH = 13
    Monochromator.NV = 13
    Monochromator.mosaich = 30
    Monochromator.mosaicv = 30
    Monochromator.r0 = 1
    Monochromator.t0 = 0  # remove transmitted neutrons
    Monochromator.RV = "2*" + str(ThALES_L) + "*sin(DEG2RAD*a2/2)"
    Monochromator.RH = "2*" + str(ThALES_L) + "/sin(DEG2RAD*a2/2)"
    Monochromator.DM = monochromator_d
    Monochromator.width = 0.25
    Monochromator.height = 0.2
    Monochromator.verbose = 1
    #    Monochromator.append_EXTEND("if(flag!=SCATTERED) ABSORB;")
    Monochromator.set_AT([0, 0, 0], RELATIVE=Monochromator_Arm)
    Monochromator.set_ROTATED([0, "a2/2", 0], RELATIVE=Monochromator_Arm)

    addMonitor(ThALES, "mono_out")

    Monochromator_Out = ThALES.add_component("Monochromator_Out", "Arm")
    Monochromator_Out.set_AT([0, 0, 0], RELATIVE=Monochromator_Arm)
    Monochromator_Out.set_ROTATED([0, "a2", 0], RELATIVE=Monochromator_Arm)

    addMonitor(ThALES, "mono_out_rot", [0, 0, 0.1])

    before_sample_slit = ThALES.add_component("before_sample_slit", "Slit")
    before_sample_slit.xwidth = 0.03
    before_sample_slit.yheight = 0.028
    before_sample_slit.set_AT([0, 0, ThALES_L - 0.250], RELATIVE=Monochromator_Out)

    #    addMonitor(ThALES, "sample", 0.200)

    sample_arm = ThALES.add_component("sample_arm", "Arm")
    sample_arm.set_AT([0, 0, ThALES_L], RELATIVE=Monochromator_Out)
    sample_arm.set_ROTATED([0, "a3", 0], RELATIVE=Monochromator_Out)

    addMonitor(ThALES, "sample", [0, 0, 0])

    mycryo, exit = cryo10T(ThALES, "10T", [0, 0, 0], sample_arm, 2)

    # res_sample = ThALES.add_component("res_sample", "Res_sample")
    # res_sample.thickness = 0.001
    # res_sample.radius = 0.005
    # res_sample.E0 = 5
    # res_sample.dE = 0.25
    # res_sample.focus_xw = 0.03
    # res_sample.focus_yh = 0.04
    # res_sample.yheight = 0.05
    # res_sample.target_index = 4
    # res_sample.set_WHEN("SAMPLE==0")
    # res_sample.set_AT([0, 0, 0], RELATIVE="sample_arm")

    v_sample = ThALES.add_component("v_sample", "V_sample")
    v_sample.radius = 0.01
    v_sample.yheight = 0.05
    v_sample.thickness = 0.001
    v_sample.focus_xw = 0.04
    v_sample.focus_yh = 0.12
    v_sample.target_z = 0.25
    #    #    v_sample.target_x = 0.10
    #    #    v_sample.target_y = 0.00
    #    # v_sample.target_index = 5
    #   v_sample.set_WHEN("SAMPLE==1")
    v_sample.set_AT([0, 0, 0], RELATIVE=sample_arm)
    v_sample.set_ROTATED([0, "a4", 0], RELATIVE=sample_arm)
    #    v_sample.append_EXTEND("if(flag==SCATTERED) ABSORB;")
    # Absorption fraction           =0.0425179
    # Single   scattering intensity =1.65546e+07 (coh=1.65473e+07 inc=7331.45)
    # Multiple scattering intensity =276313

    # quartz_sample = ThALES.add_component("quartz_sample", "Isotropic_Sqw")
    # quartz_sample.radius = 0.005
    # quartz_sample.yheight = 0.05
    # # quartz_sample.classical = 0
    # quartz_sample.p_interact = 0.5
    # quartz_sample.Sqw_coh = '"SiO2_liq.qSq"'
    # quartz_sample.Sqw_inc = '"SiO2_liq.qSq"'
    # quartz_sample.sigma_coh = 10.6
    # quartz_sample.sigma_inc = 0.0056
    # quartz_sample.sigma_abs = 0.17
    # quartz_sample.density = 2.2
    # quartz_sample.weight = 60.08
    # # quartz_sample.powder_format = '"qSq"'
    # quartz_sample.verbose = 1
    # quartz_sample.T = 273.21
    # #    quartz_sample.target_index = 3
    # # quartz_sample.append_EXTEND("if(!SCATTERED) ABSORB;")
    # #    Monochromator.

    #    quartz_sample.set_WHEN("SAMPLE==2")
    #    quartz_sample.set_AT([0, "a4", 0], RELATIVE=sample_arm)
    #    quartz_sample.set_SPLIT(20)

    if addCryostat:
        exit.set_parameters(
            radius=v_sample.radius + 1e-6,
            yheight=v_sample.yheight + 1e-6,
            priority=100000,
            material_string='"Exit"',
        )

        union_master_after_sample = ThALES.add_component(
            "master_after_sample", "Union_master"
        )
        union_master_after_sample.allow_inside_start = 1
        union_master_after_sample.set_AT([0, 0, 0], RELATIVE=mycryo.name)

    Sample_Out = ThALES.add_component("Sample_Out", "Arm")
    Sample_Out.set_AT([0, 0, 0], RELATIVE=sample_arm)
    Sample_Out.set_ROTATED([0, "a4", 0], RELATIVE=sample_arm)

    addMonitor(ThALES, "sample_out", 0)

    after_sample_slit = ThALES.add_component("after_sample_slit", "Slit")
    after_sample_slit.xwidth = 0.03
    after_sample_slit.yheight = 0.04
    after_sample_slit.set_AT([0, 0, 0.250], RELATIVE=Sample_Out)

    addMonitor(ThALES, "slit_B")

    Ana_Cradle = ThALES.add_component("Ana_Cradle", "Arm")
    Ana_Cradle.set_AT([0, 0, dist_sample_ana], RELATIVE=Sample_Out)

    addMonitor(ThALES, "analyzer_IN")

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

    addMonitor(ThALES, "detector", dist_ana_det)

    #    res_monitor = ThALES.add_component("res_monitor", "Res_monitor")
    #    res_monitor.res_sample_comp = "res_sample"
    #    res_monitor.filename = '"res_monitor"'
    #    res_monitor.xwidth = 0.05
    #    res_monitor.yheight = 0.12
    #    res_monitor.set_WHEN("SAMPLE==0")
    #    res_monitor.set_AT([0, 0, dist_ana_det], RELATIVE="Ana_Out")

    detector_all = ThALES.add_component("detector_all", "Monitor")
    detector_all.xwidth = 0.05
    detector_all.yheight = 0.12
    detector_all.restore_neutron = 0
    detector_all.set_AT([0, 0, dist_ana_det + 0.001], RELATIVE="Ana_Out")

    # detector = ThALES.add_component("detector_final", "Monitor_nD")
    # detector.xwidth = 0.05
    # detector.yheight = 0.12
    # detector.bins = 1
    # detector.options = '"p n energy"'
    # detector.set_AT([0, 0, dist_ana_det + 0.002], RELATIVE="Ana_Out")

    myinstr.add_master_parameter("a2", {"ThALES": "a2"}, unit="degree")
    myinstr.add_master_parameter("a3", {"ThALES": "a3"}, unit="degree")
    myinstr.add_master_parameter("a4", {"ThALES": "a4"}, unit="degree")
    myinstr.add_master_parameter("a6", {"ThALES": "a6"}, unit="degree")
    myinstr.master["a2"] = 79.10 * ureg.degree
    myinstr.master["a3"] = 0 * ureg.degree
    myinstr.master["a4"] = 60 * ureg.degree
    myinstr.master["a6"] = 74.34 * ureg.degree

    return myinstr


def def_instrument2():
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

    a2 = BraggAngle(
        "a2",
        comment="Angle between beam reflected by monochromator and incident beam",
        unit="degree",
    )
    a2.type = "double"
    a2.value = 33
    ThALES.parameters.add(a2)
    #    a2 = ThALES.add_parameter(
    #        "double",
    #        "a2",
    #        comment="Angle between beam reflected by monochromator and incident beam",
    #        unit="degree",
    #        value=33,
    #    )
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

    a6 = BraggAngle(
        "a6",
        comment="Angle between reflected by analyzer and incident beams",
        unit="degree",
    )
    a6.type = "double"
    a6.value = 33
    ThALES.parameters.add(a6)

    monochromator_d = ThALES.add_parameter(
        "double",
        "monochromator_d",
        comment="Monochromator lattice parameter",
        value=3.355,
        unit="angstrom",
    )
    a2.d_lattice = monochromator_d.pint_value
    a6.d_lattice = monochromator_d.pint_value

    #   ThALES.add_parameter("double", "q_x_elastic", value=1.3139)
    #   ThALES.add_parameter("double", "q_z_elastic", value=0.146)
    ThALES.add_parameter("double", "SAMPLE", value=1)

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
    HCS.E0 = "Ei"
    HCS.target_index = 1

    # override the value of lambda by the value of the angle
    ThALES.append_initialize("lambda =2*sin(0.5*a2*DEG2RAD)*monochromator_d;")
    ThALES.append_initialize('printf("lambda: %.2f\\n", lambda);')
    ThALES.append_initialize("dlambda=0.01*lambda;")
    # the following python string is an alternative way of doing:
    ThALES.add_declare_var("double", "energy")
    ThALES.append_initialize("Ei = 81.80421/(lambda*lambda);")
    ThALES.add_declare_var("double", "denergy")
    ThALES.append_initialize("dE = 2*Ei*dlambda/lambda;")
    ThALES.append_initialize('printf("lambda: %.2f +/- %.2f\\n", lambda, dlambda);')
    ThALES.append_initialize('printf("energy: %.2f +/- %.2f\\n", Ei, dE);')
    # lambda2energy = "81.80421/(" + wavelength.name + "*" + wavelength.name + ")"

    H5 = ThALES.add_component("H5", "Arm")
    H5.set_AT([0, 0, 2.155], RELATIVE=HCS)

    before_sample_slit = ThALES.add_component("before_sample_slit", "Slit")
    before_sample_slit.xwidth = 0.03
    before_sample_slit.yheight = 0.028
    before_sample_slit.set_AT([0, 0, ThALES_L - 0.250], RELATIVE="PREVIOUS")

    addMonitor(ThALES, "before_slit")

    mycryo = cryo10T(ThALES)
    mycryo.set_AT([0, 0, 0.250], RELATIVE=before_sample_slit)

    v_sample = ThALES.add_component("v_sample", "V_sample")
    v_sample.radius = 0.005
    v_sample.yheight = 0.02
    v_sample.thickness = 0.005
    v_sample.focus_xw = 0.03
    v_sample.focus_yh = 0.12
    v_sample.target_z = 0.25
    #    v_sample.target_x = 0.10
    #    v_sample.target_y = 0.00
    # v_sample.target_index = 5
    #   v_sample.set_WHEN("SAMPLE==1")
    v_sample.set_AT([0, 0, 0], RELATIVE="PREVIOUS")
    v_sample.set_ROTATED([0, 0, 0], RELATIVE="PREVIOUS")
    #    v_sample.append_EXTEND("if(flag==SCATTERED) ABSORB;")

    if addCryostat:

        union_master_after_sample = ThALES.add_component(
            "master_after_sample", "Union_master"
        )
        union_master_after_sample.allow_inside_start = 1
        union_master_after_sample.set_AT([0, 0, 0], RELATIVE="PREVIOUS")

    addMonitor(ThALES, "detector", 0.25)

    # detector = ThALES.add_component("detector_final", "Monitor_nD")
    # detector.xwidth = 0.05
    # detector.yheight = 0.12
    # detector.bins = 1
    # detector.options = '"p n energy"'
    # detector.set_AT([0, 0, dist_ana_det + 0.002], RELATIVE="Ana_Out")

    myinstr.add_master_parameter("a2", {"ThALES": "a2"}, unit="degree")
    myinstr.add_master_parameter("a3", {"ThALES": "a3"}, unit="degree")
    myinstr.add_master_parameter("a4", {"ThALES": "a4"}, unit="degree")
    myinstr.add_master_parameter("a6", {"ThALES": "a6"}, unit="degree")
    myinstr.master["a2"] = 79.10 * ureg.degree
    myinstr.master["a3"] = 0 * ureg.degree
    myinstr.master["a4"] = 60 * ureg.degree
    myinstr.master["a6"] = 74.34 * ureg.degree

    return myinstr
