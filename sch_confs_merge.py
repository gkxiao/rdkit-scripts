#!/usr/bin/env python

"""
merge_dedup_mae.py

Merge multiple MacroModel conformer ensembles and remove duplicates
based on heavy-atom RMSD.

Features
--------
1. Merge multiple .mae/.maegz files
2. Heavy-atom RMSD deduplication
3. Compatible with older Schrödinger APIs
4. Works well for MacroModel conformers
5. Optional energy sorting
6. Symmetry-independent robust workflow

Usage
-----
$SCHRODINGER/run merge_dedup_mae.py \
    -i rep1-out.maegz rep2-out.maegz rep3-out.maegz \
    -o merged_unique.maegz

Advanced
--------
$SCHRODINGER/run merge_dedup_mae.py \
    -i *.maegz \
    -o merged_unique.maegz \
    --rmsd 0.5 \
    --max_keep 500 \
    --sort_energy
"""

import os
import sys
import argparse

from schrodinger import structure
from schrodinger.structutils import rmsd


# ============================================================
# Argument Parser
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description="Merge and deduplicate MacroModel conformers"
    )

    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input mae/maegz conformer files"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output merged unique mae/maegz"
    )

    parser.add_argument(
        "--rmsd",
        type=float,
        default=0.5,
        help="Heavy atom RMSD cutoff (default: 0.5)"
    )

    parser.add_argument(
        "--max_keep",
        type=int,
        default=0,
        help="Maximum number of conformers to keep (0=no limit)"
    )

    parser.add_argument(
        "--sort_energy",
        action="store_true",
        help="Sort conformers by energy before deduplication"
    )

    parser.add_argument(
        "--energy_prop",
        default="r_mmod_Potential_Energy-OPLS",
        help="Energy property name "
             "(default: r_mmod_Potential_Energy-OPLS)"
    )

    return parser.parse_args()


# ============================================================
# Heavy Atom RMSD
# ============================================================

def calc_heavy_rmsd(st1, st2):

    heavy_atoms1 = [
        atom.index
        for atom in st1.atom
        if atom.atomic_number > 1
    ]

    heavy_atoms2 = [
        atom.index
        for atom in st2.atom
        if atom.atomic_number > 1
    ]

    try:

        val = rmsd.calculate_in_place_rmsd(
            st1,
            st2,
            atom_list1=heavy_atoms1,
            atom_list2=heavy_atoms2
        )

        return val

    except Exception as e:

        print(f"[WARNING] RMSD failed: {e}")

        return 999.0


# ============================================================
# Read Structures
# ============================================================

def read_structures(files):

    structs = []

    for f in files:

        print(f"[INFO] Reading {f}")

        reader = structure.StructureReader(f)

        for st in reader:
            structs.append(st)

    return structs


# ============================================================
# Energy Extraction
# ============================================================

def get_energy(st, propname):

    try:
        return float(st.property[propname])

    except Exception:
        return 999999.0


# ============================================================
# Deduplicate
# ============================================================

def deduplicate(structs, rmsd_cutoff):

    unique = []

    total = len(structs)

    for i, st in enumerate(structs):

        is_duplicate = False

        for ust in unique:

            r = calc_heavy_rmsd(st, ust)

            if r < rmsd_cutoff:
                is_duplicate = True
                break

        if not is_duplicate:

            unique.append(st)

        if (i + 1) % 10 == 0 or (i + 1) == total:

            print(
                f"[INFO] Processed "
                f"{i+1}/{total} "
                f"Unique={len(unique)}"
            )

    return unique


# ============================================================
# Write Output
# ============================================================

def write_output(structs, outfile):

    writer = structure.StructureWriter(outfile)

    for st in structs:
        writer.append(st)

    writer.close()


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    # --------------------------------------------------------
    # Read
    # --------------------------------------------------------

    structs = read_structures(args.input)

    print(f"[INFO] Total conformers read: {len(structs)}")

    if len(structs) == 0:

        print("[ERROR] No structures found")

        sys.exit(1)

    # --------------------------------------------------------
    # Sort by energy
    # --------------------------------------------------------

    if args.sort_energy:

        print("[INFO] Sorting conformers by energy")

        structs.sort(
            key=lambda x: get_energy(
                x,
                args.energy_prop
            )
        )

    # --------------------------------------------------------
    # Deduplicate
    # --------------------------------------------------------

    print(
        f"[INFO] Deduplicating "
        f"(RMSD cutoff = {args.rmsd:.2f} Å)"
    )

    unique = deduplicate(
        structs,
        args.rmsd
    )

    # --------------------------------------------------------
    # Max keep
    # --------------------------------------------------------

    if args.max_keep > 0:

        unique = unique[:args.max_keep]

    # --------------------------------------------------------
    # Write
    # --------------------------------------------------------

    write_output(unique, args.output)

    print()

    print("[INFO] Finished")
    print(f"[INFO] Final conformers: {len(unique)}")
    print(f"[INFO] Output: {args.output}")


# ============================================================

if __name__ == "__main__":
    main()
