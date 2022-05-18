from libpyvinyl.Instrument import Instrument


def set_vanadium_sample(instrument, calculator, position=[0, 0, 0], rotation=[0, 0, 0]):

    instrument.sample_name = "vanadium"
    instrument.sample = calculator.add_component(
        instrument.sample_name,
        "V_sample",
        AT=position,
        ROTATED=rotation,  # [0, "a4", 0],
        RELATIVE=instrument._sample_arm,
        after=instrument._sample_arm,
    )
    v_sample = instrument.sample

    p_radius = calculator.add_parameter(
        "double",
        "sample_radius",
        unit="m",
        comment="Sample radius for cylindric shape",
        value=0.01,
    )
    p_height = calculator.add_parameter(
        "double",
        "sample_height",
        unit="m",
        comment="Sample height for cylindric shape",
        value=0.01,
    )
    p_thickness = calculator.add_parameter(
        "double",
        "sample_thickness",
        unit="m",
        comment="For hollow cylinders",
        value=0.001,
    )
    v_sample.radius = "sample_radius"
    v_sample.yheight = "sample_height"
    v_sample.thickness = "sample_thickness"

    # the calculator parameters are added to the master parameters
    for pname in ["sample_radius", "sample_height", "sample_thickness"]:
        p = calculator.parameters[pname]
        instrument.add_master_parameter(
            pname, {calculator.name: pname}, unit=p.unit, comment=p.comment
        )

        instrument.master[pname].value = p.value
    # the following should be set when this function is called
    # v_sample.focus_xw = 0.04
    # v_sample.focus_yh = 0.12
    # v_sample.target_z = 0.25

    return instrument.sample
