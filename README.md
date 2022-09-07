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

myinstrument = repo.load("ILL","D22","HEAD","mcstas","quick",False)

# Getting the list of master parameters
for par in myinstrument.master:
    print(par)

# running the simulation
myinstrument.run()

# retrieving the output
output = myinstrument.output()
```
In this example code, `myinstrument` is the instrument description for
the instrument D22 at ILL, at version `HEAD`, using the simulation
software "mcstas" and the specific flavor `quick`. The last boolean
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
  from libpyvinyl import Instrument
  ```
- for mcstas simulations
  ```
  import os
  MCSTAS_PATH = os.environ["MCSTAS"]
  from mcstasscript.interface import functions
  my_configurator = functions.Configurator()
  my_configurator.set_mcstas_path(MCSTAS_PATH)
  my_configurator.set_mcrun_path(MCSTAS_PATH + "/bin/")
  ```
- import of all other needed libraries
  We recommend to use the pint library for physical quantities:
  ```
  import pint
  ureg = pint.UnitRegistry()
  ```
- one function `def_instrument() -> Instrument:` that returns the
  instrument object. Example:
  ```
  def def_instrument():
	   myinstr = Instrument("instrument_name", instrument_base_dir=".")
	   ...
	   ...
	   return myinstr
  ```
- add master parameters to the instrument to let the **SIMULATORS**
  know which are the parameters they are supposed to play with.
  ** Don't forget to set the units! **
  Example:
  ```
  myinstr.add_master_parameter("wavelength", {"D22_quick": "lambda"}, unit="angstrom")
  myinstr.add_master_parameter(
	  "collimation", {"D22_quick": "D22_collimation"}, unit="meter"
  )
  ```
- set default values for the master parameters with units
  ```
  myinstr.master["collimation"] = 2
  myinstr.master["wavelength"] = 4.5 * ureg.angstrom
  myinstr.master["collimation"] = 2 * ureg.meter
  ```
   
- define a class inheriting from libpyvinyl.Instrument
  this is the object that will be returned by the def_instrument() function
  
  This class should have the following private members:
```
    __sample_environment = None
    __sample = None
    __sample_environment_arm = None
    __sample_arm = None
    __calculator_name = "ThALES"
```

  This class should define the following methods:
  
  - ```def set_sample(self, name) -> None:
       """Always put a sample relative to the __sample_arm and after the __sample_arm component"""
   ```
  - ```   def set_sample_environment(self, name: str) -> None:
        """Adding a sample environment to the simulation"""
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

  

