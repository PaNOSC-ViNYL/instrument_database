from libpyvinyl.Instrument import Instrument
from libpyvinyl.BaseData import DataCollection
from SimExLite.SampleData import SampleData, ASEFormat
from SimExLite.SourceCalculators import GaussianSourceCalculator, get_divergence_from_beam_diameter
from SimExLite.PropagationCalculators import WPGPropagationCalculator
from SimExLite.PMICalculators import SimpleScatteringPMICalculator
from SimExLite.DiffractionCalculators import SingFELDiffractionCalculator


import pint
from pint import set_application_registry

ureg = pint.get_application_registry()


############## Mandatory method
def get_flavours():
    return ["full", "from_sample", "merge"]

class SPB_SFX(Instrument):
    """:class: SPB/SFX Instrument"""

    def __init__(self):
        self.samples =["2NIP"]
        myinstr = self
        source = GaussianSourceCalculator("gaussian_source")
        source.parameters["photon_energy"] = 9000
        source.parameters["photon_energy_relative_bandwidth"] = 1e-3
        source.parameters["beam_diameter_fwhm"] = 1e-4
        source.parameters["pulse_energy"] = 2e-3
        source.parameters["divergence"] = get_divergence_from_beam_diameter(source.parameters["photon_energy"].value, source.parameters["beam_diameter_fwhm"].value)
        source.parameters["photon_energy_spectrum_type"] = "SASE"
        source.parameters["number_of_transverse_grid_points"] = 400
        source.parameters["number_of_time_slices"] = 12
        source.parameters["z"] = 100
        source_data = source.output

        propogation = WPGPropagationCalculator(name="WPGCalculator", input=source_data)
        propogation.parameters["beamline_config_file"] = "./simple_beamline.py"
        prop_data_collection = propogation.output
        prop_data = prop_data_collection.to_list()[0]

        sample_file = "./2nip.pdb"
        sample_data = SampleData.from_file(sample_file, ASEFormat, "sample_data")

        pmi_input = DataCollection(sample_data, prop_data)
        pmi = SimpleScatteringPMICalculator(name="PMI_calculator", input=pmi_input)
        pmi_data_collection = pmi.output
        pmi_data = pmi_data_collection.to_list()[0]

        diffraction = SingFELDiffractionCalculator(name="Diffr_calculator", input=pmi_data)
        diffraction.parameters["calculate_Compton"] = False
        diffraction.parameters["slice_inverval"] = 1
        diffraction.parameters["slice_index_upper"] = 1
        diffraction.parameters["pmi_start_ID"] = 1
        diffraction.parameters["pmi_sstop_ID"] = 1
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

