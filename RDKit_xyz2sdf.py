#!/usr/bin/env python

import argparse
from rdkit import Chem
from rdkit.Geometry import Point3D


def read_xyz(xyz_file):

    elements = []
    coords = []

    with open(xyz_file, "r") as f:
        lines = f.readlines()

    natoms = int(lines[0].strip())

    for line in lines[2:2 + natoms]:

        fields = line.split()

        elem = fields[0]
        x = float(fields[1])
        y = float(fields[2])
        z = float(fields[3])

        elements.append(elem)
        coords.append((x, y, z))

    return elements, coords


def main():

    parser = argparse.ArgumentParser(
        description="Transfer XYZ coordinates onto topology SDF"
    )

    parser.add_argument(
        "-i", "--input_sdf",
        required=True,
        help="Input topology SDF"
    )

    parser.add_argument(
        "-x", "--xyz",
        required=True,
        help="Optimized XYZ"
    )

    parser.add_argument(
        "-o", "--output_sdf",
        required=True,
        help="Output optimized SDF"
    )

    args = parser.parse_args()

    # =========================
    # Read SDF
    # =========================
    suppl = Chem.SDMolSupplier(
        args.input_sdf,
        removeHs=False
    )

    mol = suppl[0]

    if mol is None:
        raise ValueError("Failed to read SDF.")

    # =========================
    # Read XYZ
    # =========================
    xyz_elements, xyz_coords = read_xyz(args.xyz)

    # =========================
    # Atom count check
    # =========================
    if mol.GetNumAtoms() != len(xyz_coords):

        raise ValueError(
            f"Atom count mismatch:\n"
            f"SDF atoms: {mol.GetNumAtoms()}\n"
            f"XYZ atoms: {len(xyz_coords)}"
        )

    # =========================
    # Element order check
    # =========================
    sdf_elements = [
        atom.GetSymbol()
        for atom in mol.GetAtoms()
    ]

    if sdf_elements != xyz_elements:

        raise ValueError(
            "Element order mismatch between SDF and XYZ."
        )

    # =========================
    # Update coordinates
    # =========================
    conf = mol.GetConformer()

    for i, (x, y, z) in enumerate(xyz_coords):

        conf.SetAtomPosition(
            i,
            Point3D(x, y, z)
        )

    # =========================
    # Write output
    # =========================
    writer = Chem.SDWriter(args.output_sdf)
    writer.write(mol)
    writer.close()

    print(f"Done: {args.output_sdf}")


if __name__ == "__main__":
    main()
