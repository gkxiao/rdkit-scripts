#!/usr/bin/env python

import os
import argparse
import subprocess

from schrodinger import structure


###############################################################################
# MacroModel fixed-column formatter
###############################################################################

def fmt(opcd,
        i1=0, i2=0, i3=0, i4=0,
        r1=0.0, r2=0.0, r3=0.0, r4=0.0):

    return (
        f" {opcd:<4}"
        f"{i1:7d}"
        f"{i2:7d}"
        f"{i3:7d}"
        f"{i4:7d}"
        f"{r1:11.4f}"
        f"{r2:11.4f}"
        f"{r3:11.4f}"
        f"{r4:11.4f}"
    )


###############################################################################
# Argument parser
###############################################################################

parser = argparse.ArgumentParser(
    description="""
Generate MacroModel LMCS conformational search jobs.

Features:
  - OPLS4 / OPLS2005
  - replicate runs
  - random seeds
  - direct bmin execution
  - MacroModel fixed-column formatting
"""
)

###############################################################################
# IO
###############################################################################

parser.add_argument(
    "-i", "--input",
    required=True,
    help="Input MAE/MAEGZ"
)

parser.add_argument(
    "-o", "--outdir",
    default="mmod_jobs",
    help="Output directory"
)

###############################################################################
# Mode
###############################################################################

parser.add_argument(
    "-m", "--mode",
    choices=["LMCS", "FLEXIBLE"],
    default="FLEXIBLE",
    help="Search mode preset"
)

###############################################################################
# Force field
###############################################################################

parser.add_argument(
    "--ff",
    choices=["opls4", "opls2005"],
    default="opls4",
    help="Force field"
)

###############################################################################
# Replicates
###############################################################################

parser.add_argument(
    "--nrep",
    type=int,
    default=3,
    help="Number of replicate searches"
)

parser.add_argument(
    "--seed",
    type=int,
    default=11111,
    help="Base random seed"
)

###############################################################################
# Optional overrides
###############################################################################

parser.add_argument("--rmsd", type=float)
parser.add_argument("--steps", type=int)
parser.add_argument("--mcnv", type=int)
parser.add_argument("--mcss", type=float)
parser.add_argument("--ewin", type=float)
parser.add_argument("--ewin2", type=float)
parser.add_argument("--maxkeep", type=int)
parser.add_argument("--conv", type=float)
parser.add_argument("--mini", type=int)
parser.add_argument("--bdco", type=float)

###############################################################################
# Other options
###############################################################################

parser.add_argument(
    "--run",
    action="store_true",
    help="Run jobs immediately using bmin"
)

parser.add_argument(
    "--first_only",
    action="store_true",
    help="Only process first structure"
)

args = parser.parse_args()

###############################################################################
# Force field mapping
###############################################################################

FF_MAP = {
    "opls2005": 14,
    "opls4": 16,
}

ffld = FF_MAP[args.ff]

###############################################################################
# Presets
###############################################################################

if args.mode == "LMCS":

    defaults = {
        "bdco": 89.4427,
        "rmsd": 0.5,
        "steps": 1000,
        "mcnv": 5,
        "mcss": 21.0,
        "ewin": 21.0,
        "ewin2": 42.0,
        "maxkeep": 833,
        "conv": 0.05,
        "mini": 2500,
    }

else:

    ###########################################################################
    # FLEXIBLE mode
    #
    # Aggressive bioactive conformer search
    ###########################################################################

    defaults = {
        "bdco": 89.4427,
        "rmsd": 0.25,
        "steps": 10000,
        "mcnv": 20,
        "mcss": 1000.0,
        "ewin": 50.0,
        "ewin2": 100.0,
        "maxkeep": 5000,
        "conv": 0.01,
        "mini": 5000,
    }

###############################################################################
# Apply defaults
###############################################################################

for key, value in defaults.items():

    if getattr(args, key) is None:
        setattr(args, key, value)

###############################################################################
# SCHRODINGER env
###############################################################################

SCHRODINGER = os.environ.get("SCHRODINGER")

if SCHRODINGER is None:
    raise RuntimeError(
        "SCHRODINGER environment variable not set."
    )

###############################################################################
# Output directory
###############################################################################

os.makedirs(args.outdir, exist_ok=True)

###############################################################################
# Read structures
###############################################################################

reader = structure.StructureReader(args.input)

mol_count = 0
job_count = 0

for idx, st in enumerate(reader, start=1):

    mol_count += 1

    if args.first_only and idx > 1:
        break

    title = st.title.strip()

    if not title:
        title = f"mol_{idx}"

    title = title.replace(" ", "_")

    ###########################################################################
    # Replicate jobs
    ###########################################################################

    for rep in range(1, args.nrep + 1):

        seed = args.seed + rep - 1

        jobname = (
            f"{title}_{args.ff}_rep{rep}"
        )

        mae_name = f"{jobname}.mae"
        com_name = f"{jobname}.com"

        mae_file = os.path.join(args.outdir, mae_name)
        com_file = os.path.join(args.outdir, com_name)

        output_maegz = f"{jobname}-out.maegz"

        print(
            f"[INFO] Generating {jobname} "
            f"(SEED={seed}, FF={args.ff})"
        )

        #######################################################################
        # Write MAE
        #######################################################################

        with structure.StructureWriter(mae_file) as writer:
            writer.append(st)

        #######################################################################
        # Keyword list
        #######################################################################

        keywords = [

            fmt("MMOD", 0, 1),

            fmt("DEBG", 55),

            fmt("FFLD", ffld, 1, 0, 0, 1.0),

            fmt("SOLV", 3, 1),

            fmt("EXNB"),

            ###################################################################
            # Long-range electrostatics cutoff
            ###################################################################

            fmt(
                "BDCO",
                0, 0, 0, 0,
                args.bdco,
                99999.0
            ),

            fmt("READ"),

            fmt("SEED", seed),

            fmt(
                "CRMS",
                0, 0, 0, 0,
                0.0,
                args.rmsd
            ),

            fmt(
                "LMCS",
                args.steps,
                0, 0, 0,
                0.0,
                0.0,
                3.0,
                6.0
            ),

            fmt("NANT"),

            fmt(
                "MCNV",
                1,
                args.mcnv
            ),

            fmt(
                "MCSS",
                2,
                0,
                0,
                0,
                args.mcss
            ),

            fmt(
                "MCOP",
                1,
                0,
                0,
                0,
                0.5
            ),

            fmt(
                "DEMX",
                0,
                args.maxkeep,
                0,
                0,
                args.ewin,
                args.ewin2
            ),

            fmt("MSYM"),

            fmt(
                "AUOP",
                0,
                0,
                0,
                0,
                100.0
            ),

            fmt(
                "AUTO",
                0,
                2,
                1,
                1,
                0.0,
                1.0,
                0.0,
                1.0
            ),

            fmt(
                "CONV",
                2,
                0,
                0,
                0,
                args.conv
            ),

            fmt(
                "MINI",
                1,
                0,
                args.mini
            ),
        ]

        #######################################################################
        # Write COM
        #######################################################################

        with open(com_file, "w") as f:

            ###################################################################
            # IMPORTANT:
            # trailing space is intentional
            ###################################################################

            f.write(f"{mae_name} \n")

            f.write(f"{output_maegz}\n")

            for line in keywords:
                f.write(line + "\n")

        #######################################################################
        # Run
        #######################################################################

        if args.run:

            print(f"[RUN] bmin {jobname}")

            cmd = [
                os.path.join(SCHRODINGER, "bmin"),
                jobname
            ]

            subprocess.run(
                cmd,
                cwd=args.outdir
            )

        job_count += 1

###############################################################################
# Summary
###############################################################################

print("")
print(f"[INFO] Structures processed : {mol_count}")
print(f"[INFO] Replicates/job       : {args.nrep}")
print(f"[INFO] Total jobs generated : {job_count}")
print("[INFO] Finished.")
