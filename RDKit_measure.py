#!/usr/bin/env python3

import sys
import csv
import gzip
import argparse

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

try:
    from rdkit.Chem.rdmolfiles import MaeMolSupplier
    HAVE_MAE = True
except Exception:
    HAVE_MAE = False


##############################################################################
# IO
##############################################################################

def load_molecules(infile):

    if infile.lower().endswith(".sdf"):

        suppl = Chem.SDMolSupplier(
            infile,
            removeHs=False
        )

        return [m for m in suppl if m is not None]

    elif infile.lower().endswith(".sdf.gz"):

        fh = gzip.open(infile, "rb")

        suppl = Chem.ForwardSDMolSupplier(
            fh,
            removeHs=False
        )

        return [m for m in suppl if m is not None]

    elif infile.lower().endswith(".mae") or \
         infile.lower().endswith(".maegz"):

        if not HAVE_MAE:
            sys.exit(
                "ERROR: This RDKit build does not support MaeMolSupplier."
            )

        suppl = MaeMolSupplier(
            infile,
            removeHs=False
        )

        return [m for m in suppl if m is not None]

    else:

        sys.exit(
            f"Unsupported file format: {infile}"
        )


##############################################################################
# Property Handling
##############################################################################

def resolve_property_name(requested, available_props):

    # exact match
    if requested in available_props:
        return requested

    requested_lower = requested.lower()

    matches = []

    for p in available_props:

        if p.lower() == requested_lower:
            matches.append(p)

    if len(matches) == 1:
        return matches[0]

    matches = []

    for p in available_props:

        if p.lower().startswith(requested_lower):
            matches.append(p)

    if len(matches) == 1:
        return matches[0]

    if len(matches) > 1:

        sys.exit(
            "\nERROR: ambiguous property name:\n"
            f"  {requested}\n\n"
            "Possible matches:\n  "
            + "\n  ".join(matches)
        )

    sys.exit(
        f"\nERROR: property not found: {requested}"
    )


##############################################################################
# Geometry
##############################################################################

def calc_distance(conf, atom1, atom2):

    return rdMolTransforms.GetBondLength(
        conf,
        atom1 - 1,
        atom2 - 1
    )


def calc_dihedral(conf, a1, a2, a3, a4):

    return rdMolTransforms.GetDihedralDeg(
        conf,
        a1 - 1,
        a2 - 1,
        a3 - 1,
        a4 - 1
    )


##############################################################################
# Main
##############################################################################

def main():

    parser = argparse.ArgumentParser(
        description="Measure distances and dihedrals from conformer ensembles"
    )

    parser.add_argument(
        "input",
        help="Input molecule file"
    )

    parser.add_argument(
        "-d",
        "--distance",
        nargs=2,
        metavar=("A1", "A2"),
        action="append",
        default=[],
        type=int,
        help="Distance measurement"
    )

    parser.add_argument(
        "-t",
        "--dihedral",
        nargs=4,
        metavar=("A1", "A2", "A3", "A4"),
        action="append",
        default=[],
        type=int,
        help="Dihedral measurement"
    )

    parser.add_argument(
        "--property",
        dest="properties",
        action="append",
        default=[],
        help="Property to output"
    )

    parser.add_argument(
        "--list-properties",
        action="store_true"
    )

    parser.add_argument(
        "--sort-by"
    )

    parser.add_argument(
        "-o",
        "--output"
    )

    args = parser.parse_args()

    mols = load_molecules(args.input)

    if len(mols) == 0:
        sys.exit("No molecules found.")

    available_props = list(
        mols[0].GetPropNames()
    )

    ##########################################################################
    # list properties
    ##########################################################################

    if args.list_properties:

        print("\nAvailable properties:\n")

        for p in sorted(available_props):
            print(p)

        return

    ##########################################################################
    # check measurements
    ##########################################################################

    if len(args.distance) == 0 and len(args.dihedral) == 0:

        sys.exit(
            "ERROR: At least one -d or -t measurement is required."
        )

    ##########################################################################
    # resolve properties
    ##########################################################################

    resolved_properties = []

    for p in args.properties:

        resolved_properties.append(
            resolve_property_name(
                p,
                available_props
            )
        )

    ##########################################################################
    # column names
    ##########################################################################

    headers = ["ConfID"]

    headers.extend(resolved_properties)

    distance_columns = []

    for atoms in args.distance:

        col = f"Dist_{atoms[0]}_{atoms[1]}"
        distance_columns.append(col)

    dihedral_columns = []

    for atoms in args.dihedral:

        col = (
            f"Dih_{atoms[0]}_"
            f"{atoms[1]}_"
            f"{atoms[2]}_"
            f"{atoms[3]}"
        )

        dihedral_columns.append(col)

    headers.extend(distance_columns)
    headers.extend(dihedral_columns)

    ##########################################################################
    # generate rows
    ##########################################################################

    rows = []

    for idx, mol in enumerate(mols, start=1):

        row = {}

        row["ConfID"] = idx

        for p in resolved_properties:

            if mol.HasProp(p):
                row[p] = mol.GetProp(p)
            else:
                row[p] = ""

        conf = mol.GetConformer()

        for atoms, col in zip(
                args.distance,
                distance_columns):

            val = calc_distance(
                conf,
                atoms[0],
                atoms[1]
            )

            row[col] = f"{val:.2f}"

        for atoms, col in zip(
                args.dihedral,
                dihedral_columns):

            val = calc_dihedral(
                conf,
                atoms[0],
                atoms[1],
                atoms[2],
                atoms[3]
            )

            row[col] = f"{val:.2f}"

        rows.append(row)

    ##########################################################################
    # sort
    ##########################################################################

    if args.sort_by:

        sort_col = args.sort_by

        if sort_col not in headers:

            resolved_sort = None

            try:
                resolved_sort = resolve_property_name(
                    sort_col,
                    headers
                )
            except Exception:
                pass

            if resolved_sort:
                sort_col = resolved_sort

        if sort_col not in headers:

            sys.exit(
                f"ERROR: sort field not found: {args.sort_by}"
            )

        def sort_key(row):

            try:
                return float(row[sort_col])
            except Exception:
                return float("inf")

        rows.sort(key=sort_key)

    ##########################################################################
    # print
    ##########################################################################

    widths = {}

    for h in headers:

        widths[h] = max(
            len(h),
            max(
                len(str(r.get(h, "")))
                for r in rows
            )
        )

    print()

    print(
        "  ".join(
            f"{h:<{widths[h]}}"
            for h in headers
        )
    )

    for row in rows:

        print(
            "  ".join(
                f"{str(row.get(h,'')):<{widths[h]}}"
                for h in headers
            )
        )

    ##########################################################################
    # csv output
    ##########################################################################

    if args.output:

        with open(
            args.output,
            "w",
            newline=""
        ) as f:

            writer = csv.DictWriter(
                f,
                fieldnames=headers
            )

            writer.writeheader()

            writer.writerows(rows)

        print(
            f"\nCSV written: {args.output}"
        )


if __name__ == "__main__":
    main()
