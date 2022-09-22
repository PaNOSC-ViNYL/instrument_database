from instrumentdatabaseapi import instrumentdatabaseapi as API

import sys
import os

# print(os.getenv("MCSTAS"))
repo = API.Repository(local_repo=".")
instrument_name = "ThALES"

repo.ls_flavours("ILL", instrument_name, "HEAD", "mcstas")
myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", dep=False)
# myinstrument = repo.load("ILL", instrument_name, "HEAD", "mcstas", "merge", dep=False)


# import the units
import pint

ureg = pint.get_application_registry()

# setting the base directory for the simulation output
basedir = "/tmp/ThALES_scan/"
myinstrument.set_instrument_base_dir(basedir)

# generation energy (monochromator)
myinstrument.master["a2"] = myinstrument.energy_to_angle(4.98 * ureg.meV)
myinstrument.master["a4"] = 60 * ureg.degree
myinstrument.master["a6"] = myinstrument.master["a2"].pint_value

myinstrument.calculators[myinstrument._calculator_name].settings(force_compile=False)
myinstrument.sim_neutrons(10000000)
myinstrument.set_instrument_base_dir(basedir + "generation/")

# the calculator has the same name as the instrument just for ThALES because it has been created like this. It is not true in general
# myThALES = myinstrument.calculators[instrument_name]
# myThALES.show_diagram()
# myThALES.show_components()

# myinstrument.set_sample_by_name("sqw")
print(myinstrument)
print("------------------------------")
print("Implemented samples: ", myinstrument.samples)
print("Current sample name: ", myinstrument.sample_name)
print("Current sample object: \n", myinstrument.sample)

print("------------------------------")
print("Implemented sample environments: ", myinstrument.sample_environments)
print("Current sample environment name: ", myinstrument.sample_environment_name)
print("Current sample environment object: \n", myinstrument.sample_environment)

# sys.exit(0)
s = myinstrument.sample
myinstrument.run()

# myinstrument.sim_neutrons(5000)
# myinstrument.run()
# print((vars(s).values()))
# print("hash: ", (hash(vars(s).values())))
# myinstrument.sample.thickness = 0.002
# print("hash after:", (frozenset(vars(s))))
# myinstrument.run()

sys.exit(0)

myinstrument = repo.load(
    "ILL", instrument_name, "HEAD", "mcstas", "from_sample", dep=False
)
myinstrument.set_instrument_base_dir(basedir)
# myThALES = myinstrument.calculators[instrument_name]
# myThALES.parameters["filelist"] = '"/tmp/ThALES_scan/generation/sSAMPLE.mcpl.gz"'
print(myinstrument)

sys.exit(0)
import os

for r in [0.002, 0.005, 0.01, 0.02]:
    for h in [0.005, 0.01, 0.02]:
        os.mkdir(basedir + "r{0:.3f}_h{1:.3f}/".format(r, h))
        myinstrument.master["sample_size_r"] = r
        myinstrument.master["sample_size_y"] = h
        energy = 4.48  # starting energy of the scan
        dE = 0.05  # scan steps
        myThALES.show_components()
        # sys.exit(0)

        for i in range(0, 21):
            if i == 0:
                myThALES.settings(force_compile=True)
            else:
                myThALES.settings(force_compile=False)

            a6.energy = energy * ureg.meV
            print("Scan Energy: %.2f" % energy)
            print(f"Scan Energy: {a6.energy}")
            print(f"Scan angle: {a6.value}")
            print(f"Generate energy: {a2.energy}")
            print(f"Generate angle: {a2.value}")
            myThALES.calculator_base_dir = "r{0:.3f}_h{1:.3f}/{2:.2f}/".format(
                r, h, energy
            )
            myinstrument.run()
            energy += dE

# for r in /tmp/ThALES_scan/r0.0*/;  for dir in $r/{4,5}*; set -l d (basename $dir); echo -en "$d\t";tail -1 $dir/detector_all.dat; end | sort > $r/scan.dat;end

#  for file in /tmp/ThALES_scan/r*/scan.dat; gnuplot -persistent -c scan.gpl  $file;end
# gnuplot  -c mysim.gpl /tmp/ThALES_scan/4.48_5/
# for dir in /tmp/ThALES_scan/*; set -l d (basename $dir); echo -en "$d\t";tail -1 $dir/detector_all.dat; end | sort > /tmp/scan.dat'
# gnuplot p '/tmp/scan.dat' u 1:2:3 w yerr
# A=0.1
# sigma=0.05
# mu=5
# gauss(x) = A*exp(-(x-mu)*(x-mu)/(sigma*sigma))
# fit gauss(x) '/tmp/scan.dat' u 1:2 via A,sigma,mu

# set fit errorvariables
# set samples 500

# cd "/tmp/ThALES_scan/4.48"
# set xrange [4:5]
# p  'H5.E',    'H53_D.E', 'H53_7.E','slit_mono.E',  'ThALES_mono_in.E',  'ThALES_mono_out.E'
