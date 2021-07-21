# instrument_database
Repository of instrument descriptions, components, etc.

To clone without a github account:
```
git clone git://github.com/PaNOSC-ViNYL/instrument_database.git
```

## Dependencies
 - [McStasScripts](https://github.com/PaNOSC-ViNYL/McStasScript)
 - Python virtual env
   - Ubuntu: python3-venv
### Testing dependencies:
 - McStas
 
## Repository structure
  - mcstas/
    - components/
    - ILL/
      - sources/
        - HEAD/source.py
      - beamlines/
        - H51/HEAD/H51.py
      - D22/
        - HEAD/D22.py:
          - components/
          - share/
  - simex
  

## Developer
### Examples
```
cmake -S . -B /dev/shm/instrument_database/
cmake --build /dev/shm/instrument_database/
ctest --test-dir /dev/shm/instrument_database/ --rerunfailed --output-on-failure
```

## TODO
 - [ ]


## Examples

### Convert mcstas format into python script for McStasScripts

```
from mcstasscript.interface import functions
my_configurator = functions.Configurator()
my_configurator.set_mcrun_path("/usr/bin/")
my_configurator.set_mcstas_path("/usr/share/mcstas/2.5/")
my_configurator.set_mxrun_path("/usr/bin/")
my_configurator.set_mcxtrace_path("/usr/share/mcxtrace/1.5/")
my_configurator.set_mcstas_path("/usr/share/mcstas/2.7/")
my_configurator.set_mcstas_path("/usr/share/mcstas/2.6.1/")

from mcstasscript.interface import reader
Reader = reader.McStas_file("/tmp/ILL_D22_quick.instr")
Reader.write_python_file("/tmp/D22_quick.py")
```

### Convert McStasScripts python script into original McStas format
From the location of the python script (`D22_quick.py` in this example)
```
# from <filename> import <instrument_name>
from D22_quick import D22_quick
# this writes into a file called <instrument_name>.instr
D22_quick.write_full_instrument()

```
