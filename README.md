# PANOSC Simulation Instrument Repository

The goal of this repository is to collect and store the description
of instruments at various research infrastructures following the
[libpyvinyl](https://github.com/PaNOSC-ViNYL/libpyvinyl) API.


It is worth dividing the users of this repository in two categories:

 - **CONTRIBUTORS**: instrument experts writing and mantaining
   up-to-date the instrument descriptions in this repository for final
   users 
 - **SIMULATORS**: end users running X-ray and neutron simulations
   using the instrument descriptions contained in this repository

## Instruments
Instruments at various research facilities can be described either in
some simplified way or up to a very high level of details and
complexity. This repository allows to:

 - split sections of a simulation in multiple parts (e.g. source of the rays, beamline, detector),
 - mantain multiple version of the same instrument in order to follow
   time evolution of the latter (e.g. before and after some upgrade),
 - represent different flavours of the same instrument (e.g. very simplified,
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

Each instrument is identified by the following informations:

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
 1. a flavour to identify alternative description of the same
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
of any software or instrument description update.

## Instructions for SIMULATORS
### Requirements
Simulators should have git and python (at least version 3.6) installed
on their machine.
A github account is **not** required.

```
pip install -r requirements/simulators.txt
```

Simulators have two equivalent ways to access the simulation
descriptions contained in this repository:
  1. Clone this repository in a local directory:
  ```
  git clone git://github.com/PaNOSC-ViNYL/instrument_database.git
  ```
  
## Instructions for CONTRIBUTORS
Contributors should have a valid github account and are supposed to
have some basic familiarity with git and python.





Repository of instrument descriptions, components, etc.

To clone without a github account:

## Dependencies
 - [McStasScripts](https://github.com/PaNOSC-ViNYL/McStasScript)
 - Python virtual env
   - Ubuntu: python3-venv

### Testing dependencies:
 - McStas
 
  
## File content conventions
The PYTHONPATH variable should point to the facility directory. In a python script for example
sys.path.append("mcstas/ILL")


## Developer
### Examples & testing
```
cmake -S . -B /dev/shm/instrument_database/
cmake --build /dev/shm/instrument_database/
ctest --test-dir /dev/shm/instrument_database/ --rerunfailed --output-on-failure
```

## TODO
	SimEx SPB_instrument.py:
	- [ ] import from pyvinyl parameters 
	- [ ] declare the instrument class as collection of the defined parameters
	McStas:
	- [ ] define D22 as instrument class object
	- [ ] Add the parameters
	

List of institutes and instruments:
```
tree -I "__pycache__|*~|#*#|__init__.py"  institutes/
```

## Examples

### Convert mcstas format into python script for McStasScripts

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

### Convert McStasScripts python script into original McStas format
From the location of the python script (`D22_quick.py` in this example)
```
# from <filename> import <instrument_name>
from D22_quick import D22_quick
# this writes into a file called <instrument_name>.instr
D22_quick.write_full_instrument()

```
