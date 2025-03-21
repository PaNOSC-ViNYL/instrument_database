import sys
import os
from instrumentdatabaseapi import instrumentdatabaseapi as API
from mcstasscript.interface import instr

import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--outdir", nargs="?", help="output directory", default=".")
parser.add_argument("institute", help="name of the institute")
parser.add_argument("instrument", help="name of the instrument")
parser.add_argument("version", help="version of the instrument")
parser.add_argument("software", help="simulation software used by the description")
parser.add_argument("flavour", help="flavour of the description")
parser.add_argument("sample_name", help="name of the sample")
parser.add_argument("sample_environment_name", help="name of the sample environment")
parser.add_argument("sample_holder_material", help="sample holder material")
parser.add_argument(
    "--dry_run", help="Do not write the instrument files", action="store_true"
)
args = parser.parse_args()


repo = API.Repository(local_repo=os.path.dirname(__file__) + "/../..")

try:
    myinstrument = repo.load(
        args.institute,
        args.instrument,
        args.version,
        args.software,
        args.flavour,
        dep=False,
    )
except RuntimeError as e:
    print(e)
    print(args)

    sys.exit(1)
    #    raise RuntimeError("")

myinstrument.set_sample_by_name(args.sample_name)
if args.sample_holder_material != "None":
    myinstrument.sample_holder(
        material=args.sample_holder_material,
        shape="box",
        w=0.02,
        h=0.03,
        d=0.0135,
        th=0.00125,
    )

myinstrument.set_sample_environment_by_name(args.sample_environment_name)

calc_with_sample = myinstrument._calculator_with_sample

instrfiles = []
import os

os.makedirs(args.outdir, exist_ok=True)
os.chdir(args.outdir)
for calcname in myinstrument.calculators:
    calc = myinstrument.calculators[calcname]
    # calc.show_components()
    if calcname == calc_with_sample.name:
        calc.name = (
            calcname
            + "_"
            + args.flavour
            + "_"
            + args.sample_name
            + "_"
            + args.sample_holder_material
            + "_"
            + args.sample_environment_name
        )
    else:
        calc.name = calcname + "_" + args.flavour

    if isinstance(calc, instr.McStas_instr):
        newname = args.instrument + "_" + calc.name + ".instr"
        instrfiles.append(newname)
        if args.dry_run is False:
            calc.write_full_instrument()
            os.rename(calc.name + ".instr", newname)

print(instrfiles)
