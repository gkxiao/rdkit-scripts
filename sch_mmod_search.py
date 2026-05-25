#!/usr/bin/env python

import os
import argparse
import subprocess

from schrodinger import structure


###############################################################################
# MacroModel COM formatter
###############################################################################

def fmt(opcd,
        i1=0, i2=0, i3=0, i4=0,
        r1=0.0, r2=0.0, r3=0.0, r4=0.0):
    """
    Generate fixed-width MacroModel COM line.
    """

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
Generate MacroModel LMCS conformational search inputs
with:

  - OPLS4 or OPLS2005
  - independent random seeds
  - aggressive flexible-ligand exploration
  - optional automatic BatchMin execution

Suitable for:
  - flexible ligands
  - IMHB systems
  - folded conformers
  - macrocycles
  - PROTAC-like molecules
""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "-i", "--input",
    required=True,
    help="Input MAE/MAEGZ file"
)

parser.add_argument(
    "-o", "--outdir",
    default="mmod_jobs",
    help="Output directory"
)

parser.add_argument(
    "-m", "--mode",
    choices=["LMCS", "FLEXIBLE"],
    default="FLEXIBLE",
    help="Preset search mode"
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
# Replicates / random seed
###############################################################################

parser.add_argument(
    "--nrep",
    type=int,
    default=3,
    help="Number of independent replicate runs"
)

parser.add_argument(
    "--seed",
    type=int,
    default=11111,
    help="Base random seed"
)

###############################################################################
# Search parameters
###############################################################################

parser.add_argument(
    "--dielectric",
    type=float,
    default=None,
    help="Effective dielectric constant"
)

parser.add_argument(
    "--rmsd",
    type=float,
    default=None,
    help="RMS pruning cutoff (Å)"
)

parser.add_argument(
    "--steps",
    type=int,
    default=None,
    help="LMCS search steps"
)

parser.add_argument(
    "--mcnv",
    type=int,
    default=None,
    help="Maximum torsions perturbed simultaneously"
)

parser.add_argument(
    "--mcss",
    type=float,
    default=None,
    help="Monte Carlo search scaling"
)

parser.add_argument(
    "--ewin",
    type=float,
    default=None,
    help="Primary energy window"
)

parser.add_argument(
    "--ewin2",
    type=float,
    default=None,
    help="Secondary energy window"
)

parser.add_argument(
    "--maxkeep",
    type=int,
    default=None,
    help="Maximum retained conformers"
)

parser.add_argument(
    "--conv",
    type=float,
    default=None,
    help="Minimization convergence threshold"
)

parser.add_argument(
    "--mini",
    type=int,
    default=None,
    help="Maximum minimization iterations"
)

###############################################################################
# Execution
###############################################################################

parser.add_argument(
    "--run",
    action="store_true",
    help="Automatically run BatchMin"
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
        "dielectric": 89.4427,
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

    defaults = {
        "dielectric": 4.0,
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
# Fill unset values
###############################################################################

for key, value in defaults.items():

    if getattr(args, key) is None:
        setattr(args, key, value)


###############################################################################
# Environment
###############################################################################

schrodinger = os.environ.get("SCHRODINGER")

if schrodinger is None:
    raise RuntimeError(
        "SCHRODINGER environment variable is not set."
    )

os.makedirs(args.outdir, exist_ok=True)


###############################################################################
# Read molecules
###############################################################################

reader = structure.StructureReader(args.input)

for mol_idx, st in enumerate(reader, start=1):

    title = st.title.strip()

    if not title:
        title = f"mol_{mol_idx}"

    title = title.replace(" ", "_")

    ###########################################################################
    # Replicate runs
    ###########################################################################

    for rep in range(1, args.nrep + 1):

        seed = args.seed + rep - 1

        jobname = (
            f"{title}_"
            f"{args.ff}_"
            f"rep{rep}"
        )

        mae_file = os.path.join(
            args.outdir,
            f"{jobname}.mae"
        )

        com_file = os.path.join(
            args.outdir,
            f"{jobname}.com"
        )

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
        # Build COM keywords
        #######################################################################

        keywords = [

            ###################################################################
            # Basic
            ###################################################################

            fmt("MMOD", 0, 1),

            fmt("DEBG", 55),

            ###################################################################
            # Force field
            ###################################################################

            fmt("FFLD", ffld, 1, 0, 0, 1.0),

            ###################################################################
            # Solvent
            ###################################################################

            fmt("SOLV", 3, 1),

            ###################################################################
            # Keep electrostatics enabled
            ###################################################################

            fmt("EXNB"),

            ###################################################################
            # Dielectric
            ###################################################################

            fmt(
                "BDCO",
                0, 0, 0, 0,
                args.dielectric,
                99999.0
            ),

            ###################################################################
            # Read structure
            ###################################################################

            fmt("READ"),

            ###################################################################
            # Random seed
            ###################################################################

            fmt("SEED", seed),

            ###################################################################
            # RMS pruning
            ###################################################################

            fmt(
                "CRMS",
                0, 0, 0, 0,
                0.0,
                args.rmsd
            ),

            ###################################################################
            # LMCS search
            ###################################################################

            fmt(
                "LMCS",
                args.steps,
                0, 0, 0,
                0.0, 0.0, 3.0, 6.0
            ),

            ###################################################################
            # Non-anti torsions
            ###################################################################

            fmt("NANT"),

            ###################################################################
            # Torsion perturbation
            ###################################################################

            fmt(
                "MCNV",
                1,
                args.mcnv
            ),

            ###################################################################
            # Search scaling
            ###################################################################

            fmt(
                "MCSS",
                2, 0, 0, 0,
                args.mcss
            ),

            ###################################################################
            # Monte Carlo optimization
            ###################################################################

            fmt(
                "MCOP",
                1, 0, 0, 0,
                0.5
            ),

            ###################################################################
            # Energy windows
            ###################################################################

            fmt(
                "DEMX",
                0,
                args.maxkeep,
                0,
                0,
                args.ewin,
                args.ewin2
            ),

            ###################################################################
            # Symmetry
            ###################################################################

            fmt("MSYM"),

            ###################################################################
            # Auto optimization
            ###################################################################

            fmt(
                "AUOP",
                0, 0, 0, 0,
                100.0
            ),

            ###################################################################
            # Automatic setup
            ###################################################################

            fmt(
                "AUTO",
                0, 2, 1, 1,
                0.0, 1.0, 0.0, 1.0
            ),

            ###################################################################
            # Convergence
            ###################################################################

            fmt(
                "CONV",
                2, 0, 0, 0,
                args.conv
            ),

            ###################################################################
            # Minimization
            ###################################################################

            fmt(
                "MINI",
                1, 0,
                args.mini
            ),
        ]

        #######################################################################
        # Write COM file
        #######################################################################

        with open(com_file, "w") as f:

            f.write(f"{os.path.basename(mae_file)}\n")
            f.write(f"{output_maegz}\n")

            for line in keywords:
                f.write(line + "\n")

        #######################################################################
        # Run BatchMin
        #######################################################################

        if args.run:

            cmd = [
                os.path.join(schrodinger, "bmin"),
                jobname
            ]

            subprocess.run(
                cmd,
                cwd=args.outdir
            )

print("[INFO] All jobs completed.")
