#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate relative conformer energies from SDF properties.

Features
--------
1. Read energies from any SDF property
2. Calculate relative energies:
       Rel_E(i) = E(i) - min(E)
3. Automatic unit conversion:
       hartree <-> kcal/mol
4. Write relative energies to a user-defined property

Author: ChatGPT
"""

import argparse
from rdkit import Chem


# =========================
# Constants
# =========================

HARTREE_TO_KCAL = 627.509474
KCAL_TO_HARTREE = 1.0 / HARTREE_TO_KCAL


# =========================
# Argument Parser
# =========================

def parse_args():

    parser = argparse.ArgumentParser(
        prog="calc_rel_energy.py",

        description=(
            "Calculate relative conformer energies "
            "from an SDF energy property."
        ),

        epilog=(
            "Examples:\n\n"

            "1. Psi4 energies already in kcal/mol:\n"
            "   python calc_rel_energy.py \\\n"
            "       -i input.sdf \\\n"
            "       -o output.sdf \\\n"
            '       --prop "Psi4_Energy (kcal/mol)" \\\n'
            "       --unit kcal \\\n"
            '       --outprop "Psi4_Rel_Energy"\n\n'

            "2. xTB energies in hartree:\n"
            "   python calc_rel_energy.py \\\n"
            "       -i input.sdf \\\n"
            "       -o output.sdf \\\n"
            '       --prop "Energy_xTB" \\\n'
            "       --unit hartree \\\n"
            '       --outprop "Rel_Energy_xTB"\n\n'

            "3. Output relative energies in hartree:\n"
            "   python calc_rel_energy.py \\\n"
            "       -i input.sdf \\\n"
            "       -o output.sdf \\\n"
            '       --prop "Energy_xTB" \\\n'
            "       --unit hartree \\\n"
            "       --outunit hartree \\\n"
            '       --outprop "Rel_Energy_xTB"'
        ),

        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input SDF file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output SDF file"
    )

    parser.add_argument(
        "--prop",
        required=True,
        help=(
            "Input energy property name in SDF\n"
            'Example: "Psi4_Energy (kcal/mol)"'
        )
    )

    parser.add_argument(
        "--unit",
        choices=["kcal", "hartree"],
        required=True,
        help=(
            "Unit of input energy property\n"
            "Choices: kcal, hartree"
        )
    )

    parser.add_argument(
        "--outprop",
        required=True,
        help=(
            "Output property name for relative energies\n"
            'Example: "Psi4_Rel_Energy"'
        )
    )

    parser.add_argument(
        "--outunit",
        choices=["kcal", "hartree"],
        default="kcal",
        help=(
            "Output unit for relative energies\n"
            "Choices: kcal, hartree\n"
            "Default: kcal"
        )
    )

    parser.add_argument(
        "--digits",
        type=int,
        default=6,
        help=(
            "Number of decimal places\n"
            "Default: 6"
        )
    )

    return parser.parse_args()


# =========================
# Read Molecules
# =========================

def read_molecules(sdf_file):

    suppl = Chem.SDMolSupplier(
        sdf_file,
        removeHs=False
    )

    mols = [mol for mol in suppl if mol is not None]

    if not mols:
        raise ValueError(
            f"No valid molecules found in: {sdf_file}"
        )

    return mols


# =========================
# Unit Conversion
# =========================

def convert_to_kcal(value, unit):

    if unit == "kcal":
        return value

    elif unit == "hartree":
        return value * HARTREE_TO_KCAL

    else:
        raise ValueError(
            f"Unsupported unit: {unit}"
        )


def convert_from_kcal(value, outunit):

    if outunit == "kcal":
        return value

    elif outunit == "hartree":
        return value * KCAL_TO_HARTREE

    else:
        raise ValueError(
            f"Unsupported output unit: {outunit}"
        )


# =========================
# Read Energies
# =========================

def get_energies(
    mols,
    prop_name,
    unit
):

    energies_kcal = []

    for i, mol in enumerate(mols, start=1):

        if not mol.HasProp(prop_name):

            raise ValueError(
                f'Molecule #{i} missing property: "{prop_name}"'
            )

        try:
            raw_e = float(
                mol.GetProp(prop_name).strip()
            )

        except Exception:

            raise ValueError(
                f'Cannot parse "{prop_name}" '
                f'in molecule #{i}'
            )

        e_kcal = convert_to_kcal(
            raw_e,
            unit
        )

        energies_kcal.append(e_kcal)

    return energies_kcal


# =========================
# Calculate Relative Energies
# =========================

def add_relative_energies(
    mols,
    energies_kcal,
    outprop,
    outunit,
    digits
):

    min_energy = min(energies_kcal)

    print("\n===================================")
    print(" Relative Energy Calculation")
    print("===================================\n")

    print(
        f"Minimum energy : "
        f"{min_energy:.8f} kcal/mol\n"
    )

    fmt = f"{{:.{digits}f}}"

    for i, (mol, e_kcal) in enumerate(
        zip(mols, energies_kcal),
        start=1
    ):

        rel_kcal = e_kcal - min_energy

        rel_out = convert_from_kcal(
            rel_kcal,
            outunit
        )

        mol.SetProp(
            outprop,
            fmt.format(rel_out)
        )

        print(
            f"Mol {i:4d} | "
            f"RelE = "
            f"{fmt.format(rel_out):>12s} "
            f"{outunit}"
        )

    print()


# =========================
# Write SDF
# =========================

def write_sdf(
    mols,
    outfile
):

    writer = Chem.SDWriter(outfile)

    for mol in mols:
        writer.write(mol)

    writer.close()


# =========================
# Main
# =========================

def main():

    args = parse_args()

    print("\nReading molecules...")
    mols = read_molecules(args.input)

    print(f"Loaded {len(mols)} molecules")

    print("\nReading energies...")

    energies_kcal = get_energies(
        mols,
        args.prop,
        args.unit
    )

    add_relative_energies(
        mols=mols,
        energies_kcal=energies_kcal,
        outprop=args.outprop,
        outunit=args.outunit,
        digits=args.digits
    )

    print("Writing output SDF...")

    write_sdf(
        mols,
        args.output
    )

    print("\nDone.")
    print(f"Output written to: {args.output}\n")


if __name__ == "__main__":
    main()
