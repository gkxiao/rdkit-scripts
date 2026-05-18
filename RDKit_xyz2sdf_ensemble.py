#!/usr/bin/env python

import argparse
from rdkit import Chem
from rdkit.Geometry import Point3D


HARTREE_TO_KCAL = 627.509474


def read_multi_xyz(xyz_file):
    """
    Read multi-structure XYZ file.

    Returns:
        conformers = [
            {
                "elements": [...],
                "coords": [...],
                "energy_hartree": float,
                "comment": str
            },
            ...
        ]
    """

    conformers = []

    with open(xyz_file, "r") as f:
        lines = f.readlines()

    i = 0

    while i < len(lines):

        line = lines[i].strip()

        if not line:
            i += 1
            continue

        natoms = int(line)

        comment = lines[i + 1].strip()

        # =========================
        # Try to parse energy
        # =========================
        energy = None

        for token in comment.replace("=", " ").split():

            try:
                energy = float(token)
                break
            except ValueError:
                continue

        elements = []
        coords = []

        for j in range(i + 2, i + 2 + natoms):

            fields = lines[j].split()

            elem = fields[0]
            x = float(fields[1])
            y = float(fields[2])
            z = float(fields[3])

            elements.append(elem)
            coords.append((x, y, z))

        conformers.append(
            {
                "elements": elements,
                "coords": coords,
                "energy_hartree": energy,
                "comment": comment
            }
        )

        i += natoms + 2

    return conformers


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Transfer multi-XYZ conformer coordinates "
            "onto topology SDF"
        )
    )

    parser.add_argument(
        "-i", "--input_sdf",
        required=True,
        help="Input topology SDF"
    )

    parser.add_argument(
        "-x", "--xyz",
        required=True,
        help="Multi-conformer XYZ file"
    )

    parser.add_argument(
        "-o", "--output_sdf",
        required=True,
        help="Output ensemble SDF"
    )

    parser.add_argument(
        "--energy_kcal",
        action="store_true",
        help="Also save Energy_xTB_kcal_mol property"
    )

    args = parser.parse_args()

    # =========================
    # Read topology molecule
    # =========================
    suppl = Chem.SDMolSupplier(
        args.input_sdf,
        removeHs=False
    )

    template_mol = suppl[0]

    if template_mol is None:
        raise ValueError("Failed to read input SDF.")

    template_elements = [
        atom.GetSymbol()
        for atom in template_mol.GetAtoms()
    ]

    # =========================
    # Read multi-XYZ
    # =========================
    conformers = read_multi_xyz(args.xyz)

    if len(conformers) == 0:
        raise ValueError("No conformers found in XYZ.")

    # =========================
    # Write output SDF
    # =========================
    writer = Chem.SDWriter(args.output_sdf)

    for idx, conf_data in enumerate(conformers, start=1):

        xyz_elements = conf_data["elements"]
        xyz_coords = conf_data["coords"]
        energy_hartree = conf_data["energy_hartree"]

        # =========================
        # Atom count check
        # =========================
        if len(template_elements) != len(xyz_elements):

            raise ValueError(
                f"Atom count mismatch in conformer {idx}"
            )

        # =========================
        # Element order check
        # =========================
        if template_elements != xyz_elements:

            raise ValueError(
                f"Element order mismatch in conformer {idx}"
            )

        # =========================
        # Copy molecule
        # =========================
        mol = Chem.Mol(template_mol)

        conf = mol.GetConformer()

        # =========================
        # Update coordinates
        # =========================
        for atom_idx, (x, y, z) in enumerate(xyz_coords):

            conf.SetAtomPosition(
                atom_idx,
                Point3D(x, y, z)
            )

        # =========================
        # Set title
        # =========================
        mol.SetProp("_Name", f"CONF_{idx}")

        # =========================
        # Save xTB energy
        # =========================
        if energy_hartree is not None:

            mol.SetProp(
                "Energy_xTB",
                f"{energy_hartree:.10f}"
            )

            if args.energy_kcal:

                energy_kcal = (
                    energy_hartree * HARTREE_TO_KCAL
                )

                mol.SetProp(
                    "Energy_xTB_kcal_mol",
                    f"{energy_kcal:.6f}"
                )

        # =========================
        # Save original comment
        # =========================
        mol.SetProp(
            "XYZ_Comment",
            conf_data["comment"]
        )

        writer.write(mol)

    writer.close()

    print(
        f"Done: wrote {len(conformers)} conformers "
        f"to {args.output_sdf}"
    )


if __name__ == "__main__":
    main()
