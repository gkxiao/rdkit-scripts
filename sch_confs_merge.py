#!/usr/bin/env python

import os
import argparse

from schrodinger import structure
from schrodinger.structutils import rmsd


###############################################################################
# Argument parser
###############################################################################

parser = argparse.ArgumentParser(
    description="""
Merge multiple MacroModel conformer ensembles and remove duplicates.

Features:
  - merge multiple MAE/MAEGZ files
  - symmetry-corrected RMSD
  - heavy-atom RMSD pruning
  - keep lowest-energy conformer
  - preserve ensemble diversity

Typical workflow:
  merge_confs.py \
      -i rep1-out.maegz rep2-out.maegz rep3-out.maegz \
      -o merged.maegz

Default behavior:
  - RMSD cutoff = 0.5 Å
  - heavy atoms only
  - symmetry-aware RMSD
  - keep lowest-energy conformer among duplicates
""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "-i", "--input",
    nargs="+",
    required=True,
    help="Input MAE/MAEGZ conformer files"
)

parser.add_argument(
    "-o", "--output",
    required=True,
    help="Output merged MAE/MAEGZ file"
)

parser.add_argument(
    "--rmsd",
    type=float,
    default=0.5,
    help="RMSD cutoff for duplicate removal (Å)"
)

parser.add_argument(
    "--heavy",
    action="store_true",
    help="Use heavy atoms only"
)

parser.add_argument(
    "--nosymm",
    action="store_true",
    help="Disable symmetry correction"
)

parser.add_argument(
    "--energy_prop",
    default="r_mmod_Potential_Energy-MM3*",
    help="Property name used for energy ranking"
)

args = parser.parse_args()


###############################################################################
# Read all conformers
###############################################################################

all_structures = []

for infile in args.input:

    print(f"[INFO] Reading {infile}")

    reader = structure.StructureReader(infile)

    for st in reader:
        all_structures.append(st)

print(f"[INFO] Total conformers loaded: {len(all_structures)}")


###############################################################################
# Energy extraction
###############################################################################

def get_energy(st, propname):

    if propname in st.property:
        return float(st.property[propname])

    return 999999.0


###############################################################################
# Sort by energy
###############################################################################

all_structures.sort(
    key=lambda x: get_energy(x, args.energy_prop)
)


###############################################################################
# Duplicate pruning
###############################################################################

unique_confs = []

for i, st in enumerate(all_structures):

    keep = True

    for ref in unique_confs:

        try:

            val = rmsd.calculate_in_place_rmsd(
                st,
                ref,
                use_symmetry=(not args.nosymm),
                use_heavy_atom=(args.heavy)
            )

        except Exception as e:

            print(f"[WARNING] RMSD failed: {e}")
            continue

        if val < args.rmsd:

            keep = False
            break

    if keep:
        unique_confs.append(st)

    if (i + 1) % 100 == 0:
        print(
            f"[INFO] Processed {i+1} / {len(all_structures)} "
            f" -> kept {len(unique_confs)}"
        )


###############################################################################
# Write output
###############################################################################

writer = structure.StructureWriter(args.output)

for idx, st in enumerate(unique_confs, start=1):

    st.title = f"conf_{idx}"

    writer.append(st)

writer.close()


###############################################################################
# Summary
###############################################################################

print("")
print("===================================================")
print(f"Input conformers : {len(all_structures)}")
print(f"Unique conformers: {len(unique_confs)}")
print(f"Output file      : {args.output}")
print("===================================================")
