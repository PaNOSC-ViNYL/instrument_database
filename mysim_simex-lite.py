from instrumentdatabaseapi import instrumentdatabaseapi as API

import sys
import os


repo = API.Repository(local_repo=".")

# import the units
import pint

ureg = pint.get_application_registry()

instrument_name = "SPB-SFX"

repo.ls_flavours("EuXFEL", instrument_name, "HEAD", "simex-lite")

# I dout it's a good idea to always force the system to install the requirements. This will get the user's environment messed up easily.
SPB_SFX = repo.load("EuXFEL", instrument_name, "HEAD", "simex-lite")


print(SPB_SFX.master)

SPB_SFX.set_instrument_base_dir("./SPB_SFX_instrument")
# The required files are not in the run_time folder

SPB_SFX.set_sample_by_file(
    "institutes/EuXFEL/instruments/SPB-SFX/HEAD/simex-lite/2nip.pdb"
)
SPB_SFX.master["energy"] = 5000

# SPB_SFX.run()
# The obvious problem is about how to set the sample
