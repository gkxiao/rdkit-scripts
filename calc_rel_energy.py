#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculate relative conformer energies from SDF properties.

Supports automatic unit conversion:
    hartree <-> kcal/mol

Formula:
    Rel_E(i) = E(i) - min(E)

Author: ChatGPT
"""

import argparse
from rdkit import Chem

# Conversion factor
HARTREE_TO_KCAL = 627.509474
KCAL_TO_HARTREE = 1.0 / HARTREE_TO_KCAL


def parse_args():

    parser = argparse.ArgumentParser(
        prog="calc_rel_energy.py",
        description=(
            "Calculate relative energies for conformers "
            "stored in an SDF file."
        ),
        epilog=(
            "Examples:\n\n"

            "1. Psi4 energy already in kcal/mol:\n"
            "   python calc_rel_energy.py \\\n"
            "       -i input.sdf \\\n"
            "       -o output.sdf \\\n"
            '       --prop "Psi4_Energy (kcal/mol)" \\\n'
            "       --unit kcal\n\n"

            "2. xTB energy in hartree:\n"
            "   python calc_rel_energy.py \\\n"
            "       -i input.sdf \\\n"
            "       -o output.sdf \\\n"
            '       --prop "Energy_xTB" \\\n'
            "       --unit hartree\n\n"

            "3. Write relative energies in hartree:\n"
            "   python calc_rel_energy.py \\\n"
            "       -i input.sdf \\\n"
            "       -o output.sdf \\\n"
            '       --prop "Energy_xTB" \\\n'
            "       --unit hartree \\\n"
            "       --outunit hartree"
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
        default="Psi4_Energy (kcal/mol)",
        help=(
            "Input energy property name\n"
            'Default: "Psi4_Energy (kcal/mol)"'
        )
    )

    parser.add_argument(
        "--unit",
        choices=["kcal", "hartree"],
        default="kcal",
        help=(
            "Unit of input energy property\n"
            "Choices: kcal, hartree\n"
            "Default: kcal"
        )
    )

    parser.add_argument(
        "--outprop",
        default="Psi4_Rel_Energy",
        help=(
            "Output relative energy property name\n"
            'Default: "Psi4_Rel_Energy"'
        )
    )

    parser.add_argument(
        "--outunit",
        choices=["kcal", "hartree"],
        default="kcal",
        help=(
            "Unit of output relative energies\n"
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


def convert_to_kcal(value, unit):

    if unit == "kcal":
        return value

    elif unit == "hartree":
        return value * HARTREE_TO_KCAL

    else:
        raise ValueError(f"Unsupported unit: {unit}")


def convert_from_kcal(value, outunit):

    if outunit == "kcal":
        return value

    elif outunit == "hartree":
        return value * KCAL_TO_HARTREE

    else:
        raise ValueError(f"Unsupported output unit: {outunit}")


def get_energies(mols, prop_name, unit):

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
                f'Cannot parse "{prop_name}" in molecule #{i}'
            )

        e_kcal = convert_to_kcal(raw_e, unit)

        energies_kcal.append(e_kcal)

    return energies_kcal


def add_relative_energies(
    mols,
    energies_kcal,
    outprop,
    outunit,
    digits
):

    min_energy = min(energies_kcal)

    print("\nMinimum energy:")
    print(f"  {min_energy:.8f} kcal/mol\n")

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
            f"RelE = {fmt.format(rel_out):>12s} "
            f"{outunit}"
        )


def write_sdf(mols, outfile):

    writer = Chem.SDWriter(outfile)

    for mol in mols:
        writer.write(mol)

    writer.close()


def main():

    args = parse_args()

    mols = read_molecules(args.input)

    energies_kcal = get_energies(
        mols,
        args.prop,
        args.unit
    )

    add_relative_energies(
        mols,
        energies_kcal,
        args.outprop,
        args.outunit,
        args.digits
    )

    write_sdf(
        mols,
        args.output
    )

    print("\nDone.")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
