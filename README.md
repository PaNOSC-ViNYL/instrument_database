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
