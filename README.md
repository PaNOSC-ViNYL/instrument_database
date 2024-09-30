# PANOSC Simulation Instrument Repository

The goal of this repository is to collect and store the description
of instruments at various research infrastructures following the
[libpyvinyl](https://github.com/PaNOSC-ViNYL/libpyvinyl) API.

**Familiarity with the `libpyvinyl` API is required.**

It is worth dividing the users of this repository in two categories:

 - **CONTRIBUTORS**: instrument experts writing and maintaining
   up-to-date the instrument descriptions in this repository for final
   users 
 - **SIMULATORS**: end users running X-ray and neutron simulations
   using the instrument descriptions contained in this repository

## Instruments
Instruments at various research facilities can be described either in
some simplified way or up to a very high level of details and
complexity. This repository allows to:

 - split sections of a simulation into multiple parts (e.g. source of the rays, beamline, detector),
 - maintain multiple versions of the same instrument in order to follow
   time evolution of the latter (e.g. before and after some upgrade),
 - represent different flavors of the same instrument (e.g. very simplified,
   intermediate, very detailed).
 - adopt different simulation programs (adopting the libpyvinyl API)
   to describe the same instrument
 - keep track of specific software dependencies for each instrument
   description as well as specialized/supplementary pieces of code not available in
   the upstream simulation programs

## Repository structure
   
In order to achieve this objectives, the repository has been
structured as follows:
```
institutes/<institute>/instruments/<instrument>/<version>/<simulation program>/<instrument><flavour>.py
```

Each instrument is identified by the following:

 1. name of the institute
 1. name of the instrument
 1. a version: 
	- it is `HEAD` for the up-to-date description of existing instruments
	- it is in the form of `YYYY_MM_DD`, with the date indicating the
		last validity of the description (e.g. it can be the last day before
		an upgrade)
 1. simulation program. Currently supported are 
	[McStasscript](https://github.com/PaNOSC-ViNYL/McStasScript) 
	and [Simex-lite](https://github.com/PaNOSC-ViNYL/SimEx-Lite),
	but any `libpyvinyl` compatible code can be used
 1. a flavor to identify alternative description of the same
	instrument.
   

At the same level as the `instruments` subdirectory, there can be
others to store part of a simulation that can be in common to multiple
instruments.
Those files can be reused in the single instrument descriptions,
avoiding duplicated code (easier to mantain, less errors). Some common
possibilities are `source` and `beamline`.


```
institutes/
├── ILL
│   ├── instruments
│   │   ├── D22
│   │   │   └── HEAD
│   │   │       └── mcstas
│   │   │           ├── D22_quick.py
│   │   │           ├── requirements.txt
│   │   │           └── validation
│   │   │               ├── Detector_1652433908.x_y
│   │   │               ├── Detector_1652434584.x_y
│   │   │               ├── H51_D22_Sample_Div_1652434584.hd_vd
│   │   │               ├── H51_D22_Sample_L_1652434584.L
│   │   │               ├── H51_D22_Sample_XY_1652434584.x_y
│   │   │               └── mccode.sim
│   │   └── Thales
│   │       └── HEAD
│   │           ├── mcstas
│   │           │   └── ThALES.py
│   │           └── ThALES_resolution_v2.instr
│   └── sources
│       └── HEAD
│           └── mcstas
│               ├── Full.py
│               └── QuickSource_D22.py
└── test_institute
    └── instruments
        └── test_instrument
            └── HEAD
                └── simex
                    ├── requirements.txt
                    └── test.py

```

For each instrument, there is also a `requirements.txt` file with the
specific python dependencies and a `validation/` subdirectory
containing the output of a test simulation obtained with the
instrument description in the repository. This is used for validation
of any update of the simulation software or of the instrument description.

## Instructions for SIMULATORS
### Requirements
Simulators should have git and python (at least version 3.6) installed
on their machine.
[McStas](http://www.mcstas.org/) should also be installed on the machine if willing to use an instrument implemented for mcstas.

A github account is **not** required.


### Set up 

Clone this repository directly in a local directory:
```
git clone  --recurse-submodules https://github.com/PaNOSC-ViNYL/instrument_database.git
cd instrument_database/
```
  
Now both the instrument repository as well as its associated API are
available.

Create an isolated virtual environment
```
python -m venv --system-site-packages --symlinks python_packages
```

Activate the environment (more info for different operating systems can be found in the table [here](https://docs.python.org/3/library/venv.html)
```
source python_packages/bin/activate
pip install --upgrade pip
```

Install the API and its dependencies:
```
pip install -r instrumentDataBaseAPI/requirements.txt
pip install -e instrumentDataBaseAPI/
```

#### Setting up McStas environment
##### Install McStas 3.5.1
 -# Download miniforge from https://conda-forge.org/miniforge/
 -# Load conda manually to setup automatic conda configuration at shell startup. For fish shell:
 `source /home/nourbakhsh/miniforge3/etc/fish/conf.d/conda.fish`
 `conda init fish`
 -# Create a mcstas dedicated environment
 `conda create --name mcstas`
 -# Activate environment 
 `conda activate mcstas`
 -# Install mcstas in the environment
 `mamba install mcstas compilers openmpi=4`
 -# Install McStasScript
 `mamba install pip jupyterlab ipympl`
 `mamba install --file McStasScript/requirements.txt`
 
In case the simulation uses mcstas as simulation program, the `MCSTAS` environment variable should be set, exported and pointing to the McStas root directory.

An example in fish shell:
```
set -x MCSTAS /usr/local/mcstas/2.7/
python mcstas/scripts/setup.py
```

### Quick Start: Accessing an Instrument Description

```
from instrumentdatabaseapi import instrumentdatabaseapi as API
repo = API.Repository(local_repo=".")

myinstrument = repo.load("ILL","ThALES","HEAD","mcstas","full",False)

# import the units
import pint
ureg = pint.get_application_registry()

# check the methods specifically defined for this instrument in the help
help(myinstrument)

# setting the base directory for the simulation output
basedir = "/tmp/ThALES_scan/"
myinstrument.set_instrument_base_dir(basedir)


# Getting the list of master parameters
print(myinstrument.master):

# Modify a master parameter value:
myinstrument.master["a2"] = myinstrument.energy_to_angle(4.98 * ureg.meV)
myinstrument.master["a4"] = 60 * ureg.degree
myinstrument.master["a6"] = myinstrument.master["a2"].pint_value

# check the list of implemented samples:
print(myinstrument.samples)

# choose a sample
myinstrument.set_sample_by_name("vanadium")

# set the sample shape
# - sample_box_shape(width, height, depth, thickness)
# - sample_cylinder_shape(radius, height, thickness)
# - sample_sphere_shape(radius, thickness)
myinstrument.sample_cylinder_shape(0.005, 0.01)

# set the number of MC neutrons to simulate
myinstrument.sim_neutrons(500000)

# fix the simulation seed to ensure reproducibility
myinstrument.set_seed(654321)

# running the simulation
myinstrument.run()

# retrieving the output
output = myinstrument.output()
```

In this example code, `myinstrument` is the instrument description for
the instrument ThALES at ILL, at version `HEAD`, using the simulation
software "mcstas" and the specific flavor `full`. The last boolean
is optional and allowes to install all the additional dependencies to
run the specific instrument using `pip`.


## Instructions for CONTRIBUTORS
Contributors **need** a valid github account and are supposed to
have some basic familiarity with git and python.

[McStas](http://www.mcstas.org/) should also be installed on the machine if willing to use an instrument implemented for mcstas.

### Development process
  1. first clone a fresh version of the repository
	 ```
	 git clone --recurse-submodules git@github.com:PaNOSC-ViNYL/instrument_database.git
	 ```
  1. create a new branch
	  ```
	  git branch SomeMeaningfulName
	  ```
  1. modify or create a new instrument
  1. run the test instrument script
	 ```
	 ./instrumentDataBaseAPI/scripts/test_instrument.py <institute> <instrument> <version> <simulation_program> <flavour>
	 ```
  1. copy the output files in a new subdirectory for the instrument
	 called `validation`
  1. Push the changes to a new branch on the repository
	 ```
	 git push origin HEAD:SomeMeaningfulName
	 ```
  1. Create a pull request on github 
  1. After review it will be integrated and available for everyone

### Guidelines for writing an instrument
The instrument description should be placed as explained in the
[Repository Structure](#repository-structure).

Mandatory content:

- import of Instrument class from libpyvinyl
  ```
# ------------------------------ Mandatory classes to use
from libpyvinyl.Instrument import Instrument
from libpyvinyl.Parameters import Parameter
  ```
- for mcstas simulations
  ```
# ------------------------------ For McStasscript instruments
import mcstasscript as ms
from mcstasscript.interface import functions
from mcstasscript.interface import instr

from mcstas.McStasInstrumentBase import McStasInstrumentBase

# this is needed to get the location of McStas executables and libraries
my_configurator = functions.Configurator()
  ```
- import of all other needed libraries
  We recommend to use the pint library for physical quantities:
  ```
import pint
from pint import set_application_registry
ureg = pint.get_application_registry()
  ```
- one function to return the list of implemented flavours for the instrument
```
############## Mandatory method
def get_flavours():
    return ["full", "nosection"]
```
- one function to return the instrument object given a flavour
```
############## Mandatory method
def def_instrument(flavour: Optional[str] = None):
    """Function returning the specialized instrument object based on the flavour requested"""
    if flavour not in get_flavours() and flavour != "":
        raise RuntimeError(f"Flavour {flavour} not in the flavour list")

    if flavour in [None, "None", "", "full"]:
        return ThALES()
    if flavour == "nosection":
        return ThALES(False)
    else:
        raise RuntimeError(f"Flavour {flavour} not implement")
```
- the instrument implemented as a class inheriting from libpyvinyl.Instrument
```
class ThALES(McStasInstrumentBase):
    """:class: Instrument class defining the ThALES instrument at ILL"""

    # ------------------------------ utility methods made available for the users

    # ------------------------------ The instrument definition goes in the __init__
    def __init__(self, do_section=True):
        """Here the real definition of the instrument is performed"""

        super().__init__("ThALES", do_section)
     ...
 ```
- add master parameters at the end of the *init* of the instrument to let the **SIMULATORS**
  know which are the parameters they are supposed to play with.
  ** Don't forget to set the units! **
  Example:
  ```
        # ------------------------------ instrument parameters
        myinstr.add_master_parameter("a2", {OriginCalc.name: "a2"}, unit="degree")
        myinstr.add_master_parameter(
            "a3", {SampleCalc.name: "sample_y_rotation"}, unit="degree"
        )
        myinstr.add_master_parameter("a4", {SampleCalc.name: "a4"}, unit="degree")
        myinstr.add_master_parameter("a6", {AnalyzerCalc.name: "a6"}, unit="degree")
  ```
- set default values for the master parameters with units
  ```
        myinstr.master["a2"] = 79.10 * ureg.degree
        myinstr.master["a3"] = 0 * ureg.degree
        myinstr.master["a4"] = 60 * ureg.degree
        myinstr.master["a6"] = 79.10 * ureg.degree
  ```
   
```


#### Useful commands
How to see the tree structure from unix:
```
tree  -I "__pycache__|*~|#*#|__init__.py" institutes/
```

##### Convert mcstas format into python script for McStasScripts

```
import os
MCSTAS_PATH = os.environ['MCSTAS']
mcstas_instrument_file = "/tmp/ILL_D22_quick.instr"
from mcstasscript.interface import functions
my_configurator = functions.Configurator()
#my_configurator.set_mcrun_path("/usr/bin/")
my_configurator.set_mcstas_path(MCSTAS_PATH)
#my_configurator.set_mxrun_path("/usr/bin/")
#my_configurator.set_mcxtrace_path("/usr/share/mcxtrace/1.5/")

from mcstasscript.interface import reader
Reader = reader.McStas_file(mcstas_instrument_file)
Reader.write_python_file("/tmp/myinstrument.py")
```

##### Convert McStasScripts python script into original McStas format
From the location of the python script (`D22_quick.py` in this example)
```
from instrumentdatabaseapi import instrumentdatabaseapi as API
repo = API.Repository(local_repo=".")

myinstrument = repo.load("ILL","D22","HEAD","mcstas","quick",False)
myinstrument.calculators["D22_quick"].write_full_instrument()

```

  

