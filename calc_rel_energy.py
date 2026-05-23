#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculate relative Psi4 energies from multi-conformer SDF.

Formula:
    Rel_E(i) = E(i) - min(E)

Author: ChatGPT
"""

import argparse
from rdkit import Chem


def parse_args():
    parser = argparse.ArgumentParser(
        prog="calc_rel_energy.py",
        description=(
            "Calculate relative energies for conformers in an SDF file "
            "using a specified energy property."
        ),
        epilog=(
            "Example:\n"
            "  python calc_rel_energy.py "
            "-i input.sdf "
            "-o output.sdf\n\n"
            "Custom property name:\n"
            "  python calc_rel_energy.py "
            "-i input.sdf "
            "-o output.sdf "
            '--prop "Psi4_Energy (kcal/mol)" '
            '--outprop "Psi4_Rel_Energy"'
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input multi-conformer SDF file"
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
            "Energy property name in SDF\n"
            'Default: "Psi4_Energy (kcal/mol)"'
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

    return parser.parse_args()


def read_molecules(sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mols = [mol for mol in suppl if mol is not None]

    if not mols:
        raise ValueError(f"No valid molecules found in: {sdf_file}")

    return mols


def get_energies(mols, prop_name):
    energies = []

    for i, mol in enumerate(mols, start=1):

        if not mol.HasProp(prop_name):
            raise ValueError(
                f'Molecule #{i} missing property: "{prop_name}"'
            )

        try:
            e = float(mol.GetProp(prop_name).strip())
        except Exception:
            raise ValueError(
                f'Cannot parse "{prop_name}" in molecule #{i}'
            )

        energies.append(e)

    return energies


def add_relative_energies(
    mols,
    energies,
    outprop
):
    min_energy = min(energies)

    print(f"\nMinimum energy:")
    print(f"  {min_energy:.6f} kcal/mol\n")

    for i, (mol, e) in enumerate(zip(mols, energies), start=1):

        rel_e = e - min_energy

        mol.SetProp(outprop, f"{rel_e:.6f}")

        print(
            f"Mol {i:4d} | "
            f"E = {e:15.6f} | "
            f"RelE = {rel_e:12.6f}"
        )


def write_sdf(mols, outfile):
    writer = Chem.SDWriter(outfile)

    for mol in mols:
        writer.write(mol)

    writer.close()


def main():

    args = parse_args()

    mols = read_molecules(args.input)

    energies = get_energies(
        mols,
        args.prop
    )

    add_relative_energies(
        mols,
        energies,
        args.outprop
    )

    write_sdf(
        mols,
        args.output
    )

    print("\nDone.")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
