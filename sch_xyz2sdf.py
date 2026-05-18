#!/usr/bin/env python

import argparse
from schrodinger import structure


def read_xyz_coords(xyz_file):
    """
    Read coordinates and element symbols from XYZ file.
    """
    coords = []
    elements = []

    with open(xyz_file, "r") as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise ValueError(f"Invalid XYZ file: {xyz_file}")

    for line in lines[2:]:
        fields = line.split()

        if len(fields) < 4:
            continue

        elem = fields[0]
        x = float(fields[1])
        y = float(fields[2])
        z = float(fields[3])

        elements.append(elem)
        coords.append((x, y, z))

    return elements, coords


def main():

    parser = argparse.ArgumentParser(
        description="Transfer optimized XYZ coordinates onto topology SDF."
    )

    parser.add_argument(
        "-i", "--input_sdf",
        required=True,
        help="Input topology SDF file"
    )

    parser.add_argument(
        "-x", "--xyz",
        required=True,
        help="Optimized XYZ file"
    )

    parser.add_argument(
        "-o", "--output_sdf",
        required=True,
        help="Output optimized SDF file"
    )

    args = parser.parse_args()

    # =========================
    # Read topology structure
    # =========================
    st = next(structure.StructureReader(args.input_sdf))

    # =========================
    # Read XYZ coordinates
    # =========================
    xyz_elements, xyz_coords = read_xyz_coords(args.xyz)

    # =========================
    # Atom count check
    # =========================
    if st.atom_total != len(xyz_coords):
        raise ValueError(
            f"Atom count mismatch:\n"
            f"SDF atoms : {st.atom_total}\n"
            f"XYZ atoms : {len(xyz_coords)}"
        )

    # =========================
    # Element order check
    # =========================
    sdf_elements = [atom.element for atom in st.atom]

    if sdf_elements != xyz_elements:
        raise ValueError(
            "Atom element order mismatch between SDF and XYZ.\n"
            "Coordinate transfer aborted."
        )

    # =========================
    # Update coordinates
    # =========================
    for atom, (x, y, z) in zip(st.atom, xyz_coords):
        atom.x = x
        atom.y = y
        atom.z = z

    # =========================
    # Write output
    # =========================
    with structure.StructureWriter(args.output_sdf) as writer:
        writer.append(st)

    print(f"Done: {args.output_sdf}")


if __name__ == "__main__":
    main()
