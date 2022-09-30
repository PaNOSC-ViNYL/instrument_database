from typing import List, Optional, Any
from libpyvinyl.Instrument import Instrument
from libpyvinyl.BaseData import DataCollection
from SimExLite.SampleData import SampleData, ASEFormat
from SimExLite.SourceCalculators import GaussianSourceCalculator
from SimExLite.SourceCalculators.GaussianSourceCalculator import (
    get_divergence_from_beam_diameter,
)
from SimExLite.PropagationCalculators import WPGPropagationCalculator
from SimExLite.PMICalculators import SimpleScatteringPMICalculator
from SimExLite.DiffractionCalculators import SingFELDiffractionCalculator

import pint

ureg = pint.get_application_registry()


############## Mandatory method
def get_flavours():
    return ["AGIPD_detector", "JUNGFRAU_detector", "DEMO"]


def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    if flavour in [None, "None", "", "DEMO"]:
        return SPB_SFX()
    elif flavour == "AGIPD_detector":
        return SPB_SFX_AGIPD()
    elif flavour == "JUNGFRAU_detector":
        return SPB_SFX_JUNGFRAU()
    else:
        raise RuntimeError(f"Flavour {flavour} not implement")


class SPB_SFX(Instrument):
    """:class: SPB/SFX Instrument"""

    def __init__(self):
        super().__init__("SPB_SFX")
        # self.samples = ["2NIP"]
        myinstr = self
        source = GaussianSourceCalculator("gaussian_source")
        source.parameters["photon_energy"] = 9000
        source.parameters["photon_energy"].add_interval(4000, 28000, True)

        source.parameters["photon_energy_relative_bandwidth"] = 1e-3
        source.parameters["beam_diameter_fwhm"] = 1e-4
        source.parameters["pulse_energy"] = 2e-3
        source.parameters["divergence"] = get_divergence_from_beam_diameter(
            source.parameters["photon_energy"].value,
            source.parameters["beam_diameter_fwhm"].value,
        )
        source.parameters["photon_energy_spectrum_type"] = "SASE"
        source.parameters["number_of_transverse_grid_points"] = 400
        source.parameters["number_of_time_slices"] = 12
        source.parameters["z"] = 100
        source_data = source.output

        propogation = WPGPropagationCalculator(name="WPGCalculator", input=source_data)
        # propogation.parameters["beamline_config_file"] = "./simple_beamline.py"

        # pmi_input = DataCollection(sample_data, prop_data)
        pmi = SimpleScatteringPMICalculator(
            name="PMI_calculator", input=DataCollection()
        )  # pmi_input)
        pmi_data_collection = pmi.output
        pmi_data = pmi_data_collection.to_list()[0]

        diffraction = SingFELDiffractionCalculator(
            name="Diffr_calculator", input=pmi_data
        )
        diffraction.parameters["calculate_Compton"] = False
        diffraction.parameters["slice_interval"] = 100
        diffraction.parameters["slice_index_upper"] = 1
        diffraction.parameters["pmi_start_ID"] = 1
        diffraction.parameters["pmi_stop_ID"] = 1
        diffraction.parameters["number_of_diffraction_patterns"] = 10
        diffraction.parameters["pixel_size"] = 1e-3
        diffraction.parameters["pixels_x"] = 80
        diffraction.parameters["pixels_y"] = 80
        diffraction.parameters["distance"] = 0.13
        diffraction.parameters["mpi_command"] = "mpirun -n 10"

        myinstr.add_calculator(source)
        myinstr.add_calculator(propogation)
        myinstr.add_calculator(pmi)
        myinstr.add_calculator(diffraction)

        myinstr.add_master_parameter(
            "photon_energy",
            {source.name: "photon_energy"},
            comment=source.parameters["photon_energy"].comment,
            unit="eV",
        )
        

        myinstr.add_master_parameter(
            "photon_energy_relative_bandwidth",
            {source.name: "photon_energy_relative_bandwidth"},
            comment=source.parameters["photon_energy_relative_bandwidth"].comment,
        )

        myinstr.add_master_parameter(
            "pulse_energy",
            {source.name: "pulse_energy"},
            comment=source.parameters["pulse_energy"].comment,
            unit="joule",
        )

        myinstr.add_master_parameter(
            "number_of_diffraction_patterns",
            {diffraction.name: "number_of_diffraction_patterns"},
            comment=diffraction.parameters["number_of_diffraction_patterns"].comment,
        )

        myinstr.add_master_parameter(
            "distance",
            {diffraction.name: "distance"},
            comment=diffraction.parameters["distance"].comment,
            unit="meter",
        )

        myinstr.add_master_parameter(
            "diffraction_mpi_command",
            {diffraction.name: "mpi_command"},
            comment=diffraction.parameters["mpi_command"].comment,
        )

        myinstr.master["photon_energy"] = source.parameters["photon_energy"].value_no_conversion
        myinstr.master["photon_energy_relative_bandwidth"] = source.parameters["photon_energy_relative_bandwidth"].value_no_conversion
        myinstr.master["pulse_energy"] = source.parameters["pulse_energy"].value_no_conversion
        myinstr.master["number_of_diffraction_patterns"] = diffraction.parameters["number_of_diffraction_patterns"].value_no_conversion
        myinstr.master["distance"] = diffraction.parameters["distance"].value_no_conversion
        myinstr.master["diffraction_mpi_command"] = diffraction.parameters["mpi_command"].value_no_conversion

    def set_sample_by_file(self, sample_file: str) -> None:
        propogation = self.calculators["WPGCalculator"]
        prop_data_collection = propogation.output
        prop_data = prop_data_collection.to_list()[0]

        # sample_file = "./2nip.pdb"
        sample_data = SampleData.from_file(sample_file, ASEFormat, "sample_data")

        self.calculators["PMI_calculator"].input = DataCollection(
            sample_data, prop_data
        )


class SPB_SFX_AGIPD(SPB_SFX):
    """:class: SPB/SFX Instrument with the AGIPD detector"""

    def __init__(self):
        super().__init__()
        diffraction = self.calculators["Diffr_calculator"]
        diffraction.parameters["pixel_size"] = 2e-4  # meter
        diffraction.parameters["pixels_x"] = 1024
        diffraction.parameters["pixels_y"] = 1024


class SPB_SFX_JUNGFRAU(SPB_SFX):
    """:class: SPB/SFX Instrument with the JUNGFRAU detector"""

    def __init__(self):
        super().__init__()
        diffraction = self.calculators["Diffr_calculator"]
        diffraction.parameters["pixel_size"] = 7.5e-05  # meter
        diffraction.parameters["pixels_x"] = 2048
        diffraction.parameters["pixels_y"] = 2048
