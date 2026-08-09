#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
sch_glidescore.py

Schrödinger Glide workflow for protein-ligand complexes.

Modes
-----

1. build

   Input:
       protein-ligand complex

   Steps:
       complex
          |
          +-- extract ligand
          |
          +-- identify contacting chains
          |
          +-- extract complete receptor chains
          |
          +-- generate Glide grid with -WAIT
          |
          +-- score original ligand pose

   Output:
       *_receptor.mae
       *_ligand.mae
       *_grid.in
       *_grid.zip
       *_score.in
       Glide score output


2. score

   Input:
       protein-ligand complex
       existing grid.zip

   Steps:
       complex
          |
          +-- extract ligand
          |
          +-- score ligand in-place against existing grid


3. dock

   Input:
       ligand file
       existing grid.zip

   Steps:
       ligand + grid
          |
          +-- Glide docking


Receptor definition
-------------------

The receptor is defined at CHAIN level.

A chain is selected if at least one non-ligand
atom of that chain is within the specified
distance from any ligand atom.

Once a chain is selected, the ENTIRE chain
is retained.

The receptor is NOT spatially cropped.

Water is NOT removed.

Other components belonging to a selected chain
are retained.

Only the ligand is excluded from the receptor.

Grid generation
---------------

Glide is executed as:

    glide -WAIT grid.in

Therefore the Python process waits until the
grid generation job is actually finished.

This avoids the situation where Glide submits
the grid-generation job to the background and
Python immediately checks for grid.zip before
the file has been generated.
"""

import argparse
import os
import sys
import subprocess

import numpy as np

from schrodinger import structure
from schrodinger.structutils import analyze


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Glide grid generation, scoring and docking "
            "for protein-ligand complexes."
        )
    )

    parser.add_argument(
        "--mode",
        required=True,
        choices=[
            "build",
            "score",
            "dock"
        ],
        help=(
            "build: generate grid and score complex ligand; "
            "score: score complex ligand using existing grid; "
            "dock: dock ligand using existing grid."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=(
            "Input structure. "
            "For build/score: protein-ligand complex. "
            "For dock: ligand file."
        )
    )

    parser.add_argument(
        "-l",
        "--ligand-asl",
        help=(
            "ASL expression selecting ligand. "
            "Required for build/score."
        )
    )

    parser.add_argument(
        "-g",
        "--grid",
        help=(
            "Existing Glide grid ZIP file. "
            "Required for score/dock."
        )
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output prefix."
    )

    parser.add_argument(
        "-d",
        "--distance",
        type=float,
        default=6.0,
        help=(
            "Distance cutoff for receptor-chain selection "
            "in Angstrom. Default: 6.0"
        )
    )

    parser.add_argument(
        "-p",
        "--padding",
        type=float,
        default=8.0,
        help=(
            "Padding around ligand bounding box "
            "in Angstrom. Default: 8.0"
        )
    )

    parser.add_argument(
        "--precision",
        choices=[
            "HTVS",
            "SP",
            "XP"
        ],
        default="SP",
        help="Glide docking precision. Default: SP."
    )

    return parser.parse_args()


# ============================================================
# Schrödinger paths
# ============================================================

def get_schrodinger_path():

    schrodinger = os.environ.get(
        "SCHRODINGER"
    )

    if not schrodinger:

        raise RuntimeError(
            "SCHRODINGER environment variable is not set."
        )

    return schrodinger


def get_glide_executable():

    return os.path.join(
        get_schrodinger_path(),
        "glide"
    )


# ============================================================
# Read structure
# ============================================================

def read_structure_file(filename):

    print()
    print(
        f"Reading structure: {filename}"
    )

    st = structure.Structure.read(
        filename
    )

    print(
        f"Number of atoms: {st.atom_total}"
    )

    return st


# ============================================================
# Extract ligand
# ============================================================

def extract_ligand(
    st,
    ligand_asl
):

    ligand_indices = analyze.evaluate_asl(
        st,
        ligand_asl
    )

    if not ligand_indices:

        raise ValueError(
            "No atoms found matching ligand ASL:\n"
            f"    {ligand_asl}"
        )

    ligand_st = st.extract(
        ligand_indices
    )

    print()
    print(
        f"Selected ligand atoms: "
        f"{len(ligand_indices)}"
    )

    return (
        ligand_st,
        ligand_indices
    )


# ============================================================
# Find contacting chains
# ============================================================

def find_receptor_chains(
    st,
    ligand_indices,
    cutoff
):
    """
    Determine receptor chains by ligand proximity.

    Chain-level definition:

        ligand
           |
           | <= cutoff
           v
        any atom
           |
           v
        chain

    Once a chain is selected, the complete chain
    is retained.
    """

    ligand_set = set(
        ligand_indices
    )

    # --------------------------------------------------------
    # Group non-ligand atoms by chain
    # --------------------------------------------------------

    chain_atoms = {}

    for atom in st.atom:

        if atom.index in ligand_set:
            continue

        chain = getattr(
            atom,
            "chain",
            ""
        )

        if not chain:
            continue

        chain_atoms.setdefault(
            chain,
            []
        ).append(
            atom.index
        )

    # --------------------------------------------------------
    # Find contacting chains
    # --------------------------------------------------------

    receptor_chains = set()

    for chain, atom_indices in chain_atoms.items():

        chain_found = False

        for ligand_idx in ligand_indices:

            for atom_idx in atom_indices:

                distance = st.measure(
                    ligand_idx,
                    atom_idx
                )

                if distance <= cutoff:

                    receptor_chains.add(
                        chain
                    )

                    chain_found = True
                    break

            if chain_found:
                break

    return receptor_chains


# ============================================================
# Extract receptor
# ============================================================

def extract_receptor(
    st,
    ligand_indices,
    ligand_asl,
    cutoff
):

    receptor_chains = find_receptor_chains(
        st,
        ligand_indices,
        cutoff
    )

    if not receptor_chains:

        raise RuntimeError(
            f"No receptor chain found within "
            f"{cutoff:.2f} A of ligand."
        )

    print()
    print(
        f"Receptor chains within {cutoff:.2f} A:"
    )

    for chain in sorted(
        receptor_chains
    ):

        print(
            f"    {chain}"
        )

    # --------------------------------------------------------
    # Construct receptor ASL
    # --------------------------------------------------------

    chain_asls = [
        f"chain = '{chain}'"
        for chain in sorted(
            receptor_chains
        )
    ]

    chain_asl = " or ".join(
        chain_asls
    )

    receptor_asl = (
        f"({chain_asl}) "
        f"and not ({ligand_asl})"
    )

    print()
    print(
        "Receptor ASL:"
    )

    print(
        f"    {receptor_asl}"
    )

    # --------------------------------------------------------
    # Evaluate ASL
    # --------------------------------------------------------

    receptor_indices = analyze.evaluate_asl(
        st,
        receptor_asl
    )

    if not receptor_indices:

        raise RuntimeError(
            "Receptor ASL selected zero atoms."
        )

    # --------------------------------------------------------
    # Extract complete chains
    # --------------------------------------------------------

    receptor_st = st.extract(
        receptor_indices
    )

    receptor_st.title = (
        f"Receptor from {st.title}"
    )

    print()
    print(
        f"Receptor atoms: "
        f"{receptor_st.atom_total}"
    )

    return (
        receptor_st,
        receptor_chains
    )


# ============================================================
# Grid parameters
# ============================================================

def calculate_grid_parameters(
    ligand_st,
    padding
):

    coords = np.array(
        [
            [
                atom.x,
                atom.y,
                atom.z
            ]
            for atom in ligand_st.atom
        ],
        dtype=float
    )

    if coords.size == 0:

        raise ValueError(
            "Ligand contains zero atoms."
        )

    min_xyz = np.min(
        coords,
        axis=0
    )

    max_xyz = np.max(
        coords,
        axis=0
    )

    ligand_size = (
        max_xyz - min_xyz
    )

    # --------------------------------------------------------
    # Grid center
    # --------------------------------------------------------

    center = analyze.center_of_mass(
        ligand_st
    )

    # --------------------------------------------------------
    # Outer box
    # --------------------------------------------------------

    outer = (
        ligand_size
        + 2.0 * padding
    )

    # --------------------------------------------------------
    # Inner box
    #
    # Glide 2024 requires INNERBOX to be an int_list.
    # --------------------------------------------------------

    inner = np.maximum(
        4.0,
        outer - 6.0
    )

    inner = np.ceil(
        inner
    ).astype(int)

    # --------------------------------------------------------
    # Print
    # --------------------------------------------------------

    print()
    print(
        "Ligand bounding box:"
    )

    print(
        f"    X = {ligand_size[0]:.3f} A"
    )

    print(
        f"    Y = {ligand_size[1]:.3f} A"
    )

    print(
        f"    Z = {ligand_size[2]:.3f} A"
    )

    print()
    print(
        "Grid center:"
    )

    print(
        f"    X = {center[0]:.3f}"
    )

    print(
        f"    Y = {center[1]:.3f}"
    )

    print(
        f"    Z = {center[2]:.3f}"
    )

    print()
    print(
        "Outer box:"
    )

    print(
        f"    X = {outer[0]:.3f} A"
    )

    print(
        f"    Y = {outer[1]:.3f} A"
    )

    print(
        f"    Z = {outer[2]:.3f} A"
    )

    print()
    print(
        "Inner box:"
    )

    print(
        f"    X = {inner[0]} A"
    )

    print(
        f"    Y = {inner[1]} A"
    )

    print(
        f"    Z = {inner[2]} A"
    )

    return (
        center,
        outer,
        inner
    )


# ============================================================
# Generate Glide grid
# ============================================================

def generate_grid(
    receptor_file,
    ligand_st,
    output,
    padding
):

    (
        center,
        outer,
        inner
    ) = calculate_grid_parameters(
        ligand_st,
        padding
    )

    grid_input = (
        f"{output}_grid.in"
    )

    grid_file = (
        f"{output}_grid.zip"
    )

    # --------------------------------------------------------
    # Remove previous grid
    # --------------------------------------------------------

    if os.path.exists(
        grid_file
    ):

        print()
        print(
            f"Removing existing grid: {grid_file}"
        )

        os.remove(
            grid_file
        )

    # --------------------------------------------------------
    # Write Glide grid input
    # --------------------------------------------------------

    with open(
        grid_input,
        "w"
    ) as f:

        f.write(
            "USECOMPMAE YES\n"
        )

        # IMPORTANT:
        # INNERBOX is int_list in Glide 2024.

        f.write(
            f"INNERBOX "
            f"{inner[0]}, "
            f"{inner[1]}, "
            f"{inner[2]}\n"
        )

        # IMPORTANT:
        # These three lines are required by the
        # user's Schrödinger 2024 workflow.

        f.write(
            f"ACTXRANGE "
            f"{outer[0]:.6f}\n"
        )

        f.write(
            f"ACTYRANGE "
            f"{outer[1]:.6f}\n"
        )

        f.write(
            f"ACTZRANGE "
            f"{outer[2]:.6f}\n"
        )

        f.write(
            f"GRID_CENTER "
            f"{center[0]:.6f},"
            f"{center[1]:.6f},"
            f"{center[2]:.6f}\n"
        )

        f.write(
            f"OUTERBOX "
            f"{outer[0]:.6f}, "
            f"{outer[1]:.6f}, "
            f"{outer[2]:.6f}\n"
        )

        f.write(
            f"GRIDFILE "
            f"{grid_file}\n"
        )

        f.write(
            f"RECEP_FILE "
            f"{receptor_file}\n"
        )

    # --------------------------------------------------------
    # Show input
    # --------------------------------------------------------

    print()
    print(
        f"Grid input: {grid_input}"
    )

    print()
    print(
        "Grid input contents:"
    )

    with open(
        grid_input,
        "r"
    ) as f:

        print(
            f.read()
        )

    # --------------------------------------------------------
    # Run Glide in WAIT mode
    #
    # IMPORTANT:
    #
    # Do NOT submit the grid generation job and then
    # immediately check for grid.zip.
    #
    # -WAIT makes Glide wait until the grid-generation
    # calculation has actually completed.
    # --------------------------------------------------------

    glide = get_glide_executable()

    command = [
        glide,
        "-WAIT",
        grid_input
    ]

    print()
    print(
        "Generating Glide grid in WAIT mode..."
    )

    print(
        "    " + " ".join(command)
    )

    result = subprocess.run(
        command,
        text=True,
        capture_output=True
    )

    if result.stdout:

        print()
        print(
            result.stdout
        )

    if result.stderr:

        print()
        print(
            result.stderr
        )

    # --------------------------------------------------------
    # Check Glide return code
    # --------------------------------------------------------

    if result.returncode != 0:

        raise RuntimeError(
            "Glide grid generation failed.\n"
            f"Return code: {result.returncode}"
        )

    # --------------------------------------------------------
    # At this point Glide has completed.
    # --------------------------------------------------------

    print()
    print(
        "Glide grid-generation process finished."
    )

    # --------------------------------------------------------
    # Verify grid
    # --------------------------------------------------------

    if not os.path.exists(
        grid_file
    ):

        raise RuntimeError(
            "Glide completed without producing "
            "the expected grid file:\n"
            f"    {grid_file}"
        )

    print()
    print(
        "=" * 70
    )

    print(
        "Glide grid generated successfully"
    )

    print(
        "=" * 70
    )

    print(
        f"Grid file: {grid_file}"
    )

    return grid_file


# ============================================================
# Write ligand
# ============================================================

def write_ligand(
    ligand_st,
    output
):

    ligand_file = (
        f"{output}_ligand.mae"
    )

    ligand_st.write(
        ligand_file
    )

    print()
    print(
        f"Ligand written: {ligand_file}"
    )

    return ligand_file


# ============================================================
# Score ligand in-place
# ============================================================

def score_ligand(
    ligand_file,
    grid_file,
    output
):

    if not os.path.exists(
        grid_file
    ):

        raise FileNotFoundError(
            "Grid file does not exist:\n"
            f"    {grid_file}"
        )

    score_input = (
        f"{output}_score.in"
    )

    jobname = (
        f"{output}_score"
    )

    # --------------------------------------------------------
    # Glide score-in-place input
    # --------------------------------------------------------

    with open(
        score_input,
        "w"
    ) as f:

        f.write(
            f"JOBNAME {jobname}\n"
        )

        f.write(
            f"GRIDFILE {grid_file}\n"
        )

        f.write(
            f"LIGANDFILE {ligand_file}\n"
        )

        f.write(
            "DOCKING_METHOD inplace\n"
        )

    print()
    print(
        f"Score input: {score_input}"
    )

    print()
    print(
        "Score input contents:"
    )

    with open(
        score_input,
        "r"
    ) as f:

        print(
            f.read()
        )

    glide = get_glide_executable()

    # --------------------------------------------------------
    # Use WAIT here as well.
    #
    # This guarantees that the score calculation itself
    # has completed before the Python script continues.
    # --------------------------------------------------------

    command = [
        glide,
        "-WAIT",
        score_input
    ]

    print()
    print(
        "Running Glide score in WAIT mode..."
    )

    print(
        "    " + " ".join(command)
    )

    result = subprocess.run(
        command,
        text=True,
        capture_output=True
    )

    if result.stdout:

        print()
        print(
            result.stdout
        )

    if result.stderr:

        print()
        print(
            result.stderr
        )

    if result.returncode != 0:

        raise RuntimeError(
            "Glide score calculation failed.\n"
            f"Return code: {result.returncode}"
        )

    print()
    print(
        "Glide score calculation finished."
    )

    return jobname


# ============================================================
# Dock
# ============================================================

def dock_ligand(
    ligand_file,
    grid_file,
    output,
    precision
):

    if not os.path.exists(
        grid_file
    ):

        raise FileNotFoundError(
            "Grid file does not exist:\n"
            f"    {grid_file}"
        )

    dock_input = (
        f"{output}_dock.in"
    )

    with open(
        dock_input,
        "w"
    ) as f:

        f.write(
            f"JOBNAME {output}\n"
        )

        f.write(
            f"GRIDFILE {grid_file}\n"
        )

        f.write(
            f"LIGANDFILE {ligand_file}\n"
        )

        f.write(
            f"PRECISION {precision}\n"
        )

    print()
    print(
        f"Docking input: {dock_input}"
    )

    print()
    print(
        "Docking input contents:"
    )

    with open(
        dock_input,
        "r"
    ) as f:

        print(
            f.read()
        )

    glide = get_glide_executable()

    # --------------------------------------------------------
    # WAIT mode
    #
    # This guarantees that docking has actually completed
    # before the Python script exits.
    # --------------------------------------------------------

    command = [
        glide,
        "-WAIT",
        dock_input
    ]

    print()
    print(
        "Running Glide docking in WAIT mode..."
    )

    print(
        "    " + " ".join(command)
    )

    result = subprocess.run(
        command,
        text=True,
        capture_output=True
    )

    if result.stdout:

        print()
        print(
            result.stdout
        )

    if result.stderr:

        print()
        print(
            result.stderr
        )

    if result.returncode != 0:

        raise RuntimeError(
            "Glide docking failed.\n"
            f"Return code: {result.returncode}"
        )

    print()
    print(
        "Glide docking finished."
    )


# ============================================================
# BUILD mode
# ============================================================

def mode_build(args):

    if not args.ligand_asl:

        raise ValueError(
            "--ligand-asl is required in build mode."
        )

    # --------------------------------------------------------
    # Read complex
    # --------------------------------------------------------

    st = read_structure_file(
        args.input
    )

    # --------------------------------------------------------
    # Extract ligand
    # --------------------------------------------------------

    (
        ligand_st,
        ligand_indices
    ) = extract_ligand(
        st,
        args.ligand_asl
    )

    ligand_file = write_ligand(
        ligand_st,
        args.output
    )

    # --------------------------------------------------------
    # Extract receptor
    # --------------------------------------------------------

    (
        receptor_st,
        receptor_chains
    ) = extract_receptor(
        st,
        ligand_indices,
        args.ligand_asl,
        args.distance
    )

    receptor_file = (
        f"{args.output}_receptor.mae"
    )

    receptor_st.write(
        receptor_file
    )

    print()
    print(
        f"Receptor written: {receptor_file}"
    )

    # --------------------------------------------------------
    # Generate grid
    #
    # generate_grid() does NOT return until Glide has
    # actually finished because it uses -WAIT.
    # --------------------------------------------------------

    grid_file = generate_grid(
        receptor_file,
        ligand_st,
        args.output,
        args.padding
    )

    # --------------------------------------------------------
    # Score original ligand
    #
    # At this point grid_file is guaranteed to exist.
    # --------------------------------------------------------

    score_ligand(
        ligand_file,
        grid_file,
        args.output
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print(
        "=" * 70
    )

    print(
        "BUILD completed successfully"
    )

    print(
        "=" * 70
    )

    print(
        f"Receptor : {receptor_file}"
    )

    print(
        f"Ligand   : {ligand_file}"
    )

    print(
        f"Grid     : {grid_file}"
    )

    print(
        "Chains   : "
        + ", ".join(
            sorted(
                receptor_chains
            )
        )
    )


# ============================================================
# SCORE mode
# ============================================================

def mode_score(args):

    if not args.grid:

        raise ValueError(
            "--grid is required in score mode."
        )

    if not args.ligand_asl:

        raise ValueError(
            "--ligand-asl is required in score mode."
        )

    # --------------------------------------------------------
    # IMPORTANT:
    # Check grid BEFORE starting score.
    # --------------------------------------------------------

    if not os.path.exists(
        args.grid
    ):

        raise FileNotFoundError(
            "Specified Glide grid does not exist:\n"
            f"    {args.grid}"
        )

    # --------------------------------------------------------
    # Read complex
    # --------------------------------------------------------

    st = read_structure_file(
        args.input
    )

    # --------------------------------------------------------
    # Extract ligand
    # --------------------------------------------------------

    (
        ligand_st,
        ligand_indices
    ) = extract_ligand(
        st,
        args.ligand_asl
    )

    ligand_file = write_ligand(
        ligand_st,
        args.output
    )

    # --------------------------------------------------------
    # Score
    # --------------------------------------------------------

    score_ligand(
        ligand_file,
        args.grid,
        args.output
    )

    print()
    print(
        "=" * 70
    )

    print(
        "SCORE completed successfully"
    )

    print(
        "=" * 70
    )

    print(
        f"Grid   : {args.grid}"
    )

    print(
        f"Ligand : {ligand_file}"
    )


# ============================================================
# DOCK mode
# ============================================================

def mode_dock(args):

    if not args.grid:

        raise ValueError(
            "--grid is required in dock mode."
        )

    # --------------------------------------------------------
    # IMPORTANT:
    # Check grid BEFORE starting docking.
    # --------------------------------------------------------

    if not os.path.exists(
        args.grid
    ):

        raise FileNotFoundError(
            "Specified Glide grid does not exist:\n"
            f"    {args.grid}"
        )

    # --------------------------------------------------------
    # The input is the ligand file in dock mode.
    # --------------------------------------------------------

    ligand_file = args.input

    if not os.path.exists(
        ligand_file
    ):

        raise FileNotFoundError(
            f"Ligand file not found:\n"
            f"    {ligand_file}"
        )

    # --------------------------------------------------------
    # Dock
    # --------------------------------------------------------

    dock_ligand(
        ligand_file,
        args.grid,
        args.output,
        args.precision
    )

    print()
    print(
        "=" * 70
    )

    print(
        "DOCK completed successfully"
    )

    print(
        "=" * 70
    )

    print(
        f"Grid     : {args.grid}"
    )

    print(
        f"Ligand   : {ligand_file}"
    )

    print(
        f"Precision: {args.precision}"
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    try:

        if args.mode == "build":

            mode_build(
                args
            )

        elif args.mode == "score":

            mode_score(
                args
            )

        elif args.mode == "dock":

            mode_dock(
                args
            )

    except KeyboardInterrupt:

        print()
        print(
            "Interrupted by user."
        )

        sys.exit(130)

    except Exception as exc:

        print()
        print(
            "ERROR:"
        )

        print(
            str(exc)
        )

        sys.exit(1)


if __name__ == "__main__":

    main()
