import shutil
from instrumentdatabaseapi import instrumentdatabaseapi as API
from SimExLite.WavefrontData import WPGFormat
import wpg
from wpg import wpg_uti_wf, srwlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors


repo = API.Repository(local_repo=".")

# import the units
import pint

ureg = pint.get_application_registry()

instrument_name = "SPB-SFX"

repo.ls_flavours("EuXFEL", instrument_name, "HEAD", "simex-lite")

SPB_SFX = repo.load("EuXFEL", instrument_name, "HEAD", "simex-lite")

print(SPB_SFX.master)
# print(SPB_SFX)

# shutil.rmtree("./SPB_SFX_instrument", ignore_errors=True)

SPB_SFX.set_instrument_base_dir("./SPB_SFX_instrument")

SPB_SFX.set_sample_by_file(
    "institutes/EuXFEL/instruments/SPB-SFX/HEAD/simex-lite/2nip.pdb"
)
# SPB_SFX.master["energy"] = 6000
SPB_SFX.master["energy"] = 20000

SPB_SFX.run()

# Visualization
mwf = wpg.Wavefront()
source_WPG = (
    SPB_SFX.calculators["gaussian_source"]
    .output.to_list()[0]
    .write("source.h5", WPGFormat)
)
mwf.load_hdf5(source_WPG.filename)
wpg_uti_wf.plot_intensity_map(mwf)
wpg_uti_wf.integral_intensity(mwf)
srwlib.srwl.SetRepresElecField(mwf._srwl_wf, "f")
wpg_uti_wf.integral_intensity(mwf)

mwf = wpg.Wavefront()
prop_WPG = SPB_SFX.calculators["WPGCalculator"].output.to_list()[0]
mwf.load_hdf5(prop_WPG.filename)
wpg_uti_wf.plot_intensity_map(mwf)
wpg_uti_wf.integral_intensity(mwf)
srwlib.srwl.SetRepresElecField(mwf._srwl_wf, "f")
wpg_uti_wf.integral_intensity(mwf)

diffr = SPB_SFX.calculators["Diffr_calculator"].output.to_list()[0]
diffr_data = diffr.get_data()
fig, ax = plt.subplots(2, 5, figsize=(18, 8))
for i in range(2):
    for j in range(5):
        ax[i, j].imshow(diffr_data["img_array"][j + i * 5], norm=colors.LogNorm())
        ax[i, j].get_xaxis().set_visible(False)
        ax[i, j].get_yaxis().set_visible(False)

plt.tight_layout()
plt.show()
