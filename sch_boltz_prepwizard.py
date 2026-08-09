#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
boltz_prepwizard.py

Prepare Boltz co-folding results using Schrödinger PrepWizard.

Example
-------
$SCHRODINGER/run boltz_prepwizard.py \
    -i incyte.yaml \
    -ref 7S8N_X2.maegz \
    -oprefox incyte \
    -renumber_shift 27


Expected Boltz directory
------------------------
boltz_results_incyte/
└── predictions/
    └── incyte/
        ├── incyte_model_0.cif
        ├── incyte_model_1.cif
        ├── ...
        ├── confidence_incyte_model_0.json
        ├── confidence_incyte_model_1.json
        └── ...


Output
------
incyte_model_0.maegz
incyte_model_1.maegz
...

Each MAEGZ contains one complex Structure.

Structure title:
    model_0
    model_1
    ...

Boltz confidence properties:
    r_user_confidence_score
    r_user_ptm
    r_user_iptm
    r_user_ligand_iptm
    r_user_protein_iptm
    r_user_complex_plddt
    r_user_complex_iplddt
    r_user_complex_pde
    r_user_complex_ipde

Nested JSON objects are stored as string properties:
    s_user_chains_ptm
    s_user_pair_chains_iptm


Residue renumbering
-------------------
If:

    -renumber_shift 27

then, AFTER PrepWizard:

    residue 1 -> 28
    residue 2 -> 29
    residue 3 -> 30
    ...

Only protein residues are renumbered.


Serial JobControl
-----------------
Only one PrepWizard job is submitted at a time:

    model_0 -> wait -> model_1 -> wait -> model_2 -> ...

This is intentional to avoid consuming multiple Schrödinger
license tokens simultaneously.
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

from schrodinger import structure
from schrodinger.job import jobcontrol


# ============================================================================
# Defaults
# ============================================================================

DEFAULT_HOST = "localhost:4"


# ============================================================================
# Command line
# ============================================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Prepare Boltz co-folding CIF results "
            "with Schrödinger PrepWizard."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Boltz input YAML file, e.g. incyte.yaml",
    )

    parser.add_argument(
        "-ref",
        "--reference",
        default=None,
        help=(
            "Optional reference structure. "
            "Passed to -reference_st_file."
        ),
    )

    # Keep the original user's spelling -oprefox.
    # Also accept the corrected -oprefix.
    parser.add_argument(
        "-oprefox",
        "-oprefix",
        "--oprefox",
        "--oprefix",
        dest="output_prefix",
        default=None,
        help=(
            "Output filename prefix. "
            "Default: input YAML stem."
        ),
    )

    parser.add_argument(
        "-renumber_shift",
        "--renumber_shift",
        type=int,
        default=0,
        help=(
            "Residue-number offset applied AFTER PrepWizard. "
            "Example: -renumber_shift 27 "
            "changes residue 1 to 28."
        ),
    )

    parser.add_argument(
        "-HOST",
        "--host",
        default=DEFAULT_HOST,
        help=(
            "PrepWizard JobControl host. "
            "Default: localhost:4"
        ),
    )

    parser.add_argument(
        "--no-samplewater",
        action="store_true",
        help="Do not use PrepWizard -samplewater.",
    )

    return parser.parse_args()


# ============================================================================
# Basic utilities
# ============================================================================

def get_input_stem(input_yaml):
    """
    incyte.yaml -> incyte
    """
    return Path(input_yaml).stem


def find_prediction_directory(input_yaml):
    """
    Locate:

        boltz_results_${input}/predictions/${input}
    """

    input_stem = get_input_stem(input_yaml)

    result_dir = (
        Path.cwd()
        / f"boltz_results_{input_stem}"
    )

    prediction_dir = (
        result_dir
        / "predictions"
        / input_stem
    )

    if not result_dir.is_dir():

        raise FileNotFoundError(
            "Boltz result directory not found:\n"
            f"    {result_dir}"
        )

    if not prediction_dir.is_dir():

        raise FileNotFoundError(
            "Boltz prediction directory not found:\n"
            f"    {prediction_dir}"
        )

    return prediction_dir


def extract_model_number(filename):
    """
    Extract model number from:

        incyte_model_0.cif
        incyte_model_1.cif
        incyte_model_10.cif
    """

    basename = Path(filename).name

    match = re.search(
        r"_model_(\d+)\.cif$",
        basename,
    )

    if match is None:
        return None

    return int(match.group(1))


def find_cif_files(prediction_dir):
    """
    Find all *_model_*.cif files and sort numerically.
    """

    files = []

    for cif_file in prediction_dir.glob(
        "*_model_*.cif"
    ):

        model_number = extract_model_number(
            cif_file
        )

        if model_number is not None:

            files.append(
                (
                    model_number,
                    cif_file,
                )
            )

    files.sort(
        key=lambda x: x[0]
    )

    return [
        cif_file
        for model_number, cif_file
        in files
    ]


def check_file_exists(filename, description):
    """
    Validate an input file.
    """

    path = (
        Path(filename)
        .expanduser()
        .resolve()
    )

    if not path.is_file():

        raise FileNotFoundError(
            f"{description} not found:\n"
            f"    {path}"
        )

    return path


# ============================================================================
# Confidence JSON
# ============================================================================

def find_confidence_json(cif_file):
    """
    For:

        incyte_model_0.cif

    find:

        confidence_incyte_model_0.json
    """

    cif_file = Path(cif_file)

    json_name = (
        "confidence_"
        + cif_file.stem
        + ".json"
    )

    json_file = (
        cif_file.parent
        / json_name
    )

    if not json_file.is_file():

        raise FileNotFoundError(
            "Confidence JSON not found for:\n"
            f"    {cif_file.name}\n"
            "Expected:\n"
            f"    {json_file.name}"
        )

    return json_file


def read_confidence_json(json_file):
    """
    Read and validate Boltz confidence JSON.
    """

    with open(
        json_file,
        "r",
        encoding="utf-8",
    ) as f:

        data = json.load(f)

    if not isinstance(data, dict):

        raise ValueError(
            "Confidence JSON must contain "
            "a JSON object:\n"
            f"    {json_file}"
        )

    return data


# ============================================================================
# Structure property handling
# ============================================================================

def make_user_property_name(
    value,
    key,
):
    """
    Convert a Python value and JSON key into
    a Maestro user-property name.

    float:
        confidence_score
        ->
        r_user_confidence_score

    int:
        some_integer
        ->
        i_user_some_integer

    bool:
        some_flag
        ->
        b_user_some_flag

    string:
        some_text
        ->
        s_user_some_text

    dict/list:
        chains_ptm
        ->
        s_user_chains_ptm
    """

    # bool must be checked before int because bool
    # is a subclass of int in Python.
    if isinstance(value, bool):

        property_type = "b"

    elif isinstance(value, int):

        property_type = "i"

    elif isinstance(value, float):

        property_type = "r"

    elif isinstance(value, str):

        property_type = "s"

    elif isinstance(
        value,
        (dict, list, tuple),
    ):

        property_type = "s"

    else:

        raise TypeError(
            f"Unsupported JSON value type "
            f"for '{key}': "
            f"{type(value).__name__}"
        )

    safe_key = re.sub(
        r"[^A-Za-z0-9_]+",
        "_",
        str(key),
    )

    return (
        f"{property_type}_user_{safe_key}"
    )


def convert_property_value(value):
    """
    Convert a JSON value into a value suitable
    for a Maestro Structure property.

    Nested objects/lists are preserved as
    compact JSON strings.
    """

    if isinstance(
        value,
        (dict, list, tuple),
    ):

        return json.dumps(
            value,
            ensure_ascii=False,
            separators=(",", ":"),
        )

    return value


def add_confidence_properties(
    st,
    confidence_data,
):
    """
    Write every top-level Boltz confidence JSON item
    as a Structure-level user property.
    """

    for key, value in confidence_data.items():

        property_name = (
            make_user_property_name(
                value,
                key,
            )
        )

        property_value = (
            convert_property_value(
                value
            )
        )

        st.property[
            property_name
        ] = property_value

        print(
            "    Property: "
            f"{property_name} = "
            f"{property_value}"
        )


# ============================================================================
# Post-PrepWizard structure processing
# ============================================================================

def read_single_structure(maegz_file):
    """
    Read exactly one Structure from a MAEGZ file.

    A Boltz model represents one complex,
    therefore one Structure is expected.
    """

    with structure.StructureReader(
        str(maegz_file)
    ) as reader:

        structures = list(reader)

    if len(structures) != 1:

        raise RuntimeError(
            "Expected exactly one Structure in:\n"
            f"    {maegz_file}\n"
            f"Found: {len(structures)}"
        )

    return structures[0]


def renumber_protein_residues(
    st,
    shift,
):
    """
    Add `shift` to protein residue numbers.

    This operation is intentionally performed
    AFTER PrepWizard.

    Example:

        shift = 27

        1 -> 28
        2 -> 29
        3 -> 30

    Non-protein residues are not changed.
    """

    if shift == 0:

        return 0

    changed = 0

    for residue in st.residue:

        if not residue.is_protein:

            continue

        old_resnum = residue.resnum

        residue.resnum = (
            old_resnum + shift
        )

        changed += 1

    return changed


def finalize_structure(
    maegz_file,
    model_number,
    confidence_data,
    renumber_shift,
):
    """
    Finalize a PrepWizard output:

    1. Read the single complex Structure.
    2. Set title to model_i.
    3. Optionally renumber protein residues.
    4. Add all Boltz confidence JSON values.
    5. Write the Structure back to the same MAEGZ.
    """

    maegz_file = Path(
        maegz_file
    )

    st = read_single_structure(
        maegz_file
    )

    # --------------------------------------------------------------
    # Structure title
    # --------------------------------------------------------------

    title = (
        f"model_{model_number}"
    )

    st.title = title

    print(
        f"    Structure title: {title}"
    )

    # --------------------------------------------------------------
    # Optional residue renumbering
    # --------------------------------------------------------------

    if renumber_shift != 0:

        changed = (
            renumber_protein_residues(
                st,
                renumber_shift,
            )
        )

        print(
            f"    Renumber shift: "
            f"{renumber_shift}"
        )

        print(
            f"    Protein residues changed: "
            f"{changed}"
        )

    else:

        print(
            "    Renumber shift: 0 "
            "(no renumbering)"
        )

    # --------------------------------------------------------------
    # Boltz confidence properties
    # --------------------------------------------------------------

    print(
        "    Adding Boltz confidence properties:"
    )

    add_confidence_properties(
        st,
        confidence_data,
    )

    # --------------------------------------------------------------
    # IMPORTANT:
    #
    # Schrödinger 2024 structure API does NOT provide:
    #
    #     structure.write(...)
    #
    # Instead, write the Structure itself:
    #
    #     st.write(...)
    #
    # --------------------------------------------------------------

    st.write(
        str(maegz_file)
    )

    print(
        f"    Finalized: {maegz_file}"
    )


# ============================================================================
# Output filename
# ============================================================================

def make_output_filename(
    cif_file,
    output_prefix,
):
    """
    Convert:

        incyte_model_0.cif

    to:

        incyte_model_0.maegz
    """

    model_number = (
        extract_model_number(
            cif_file
        )
    )

    if model_number is None:

        raise RuntimeError(
            "Cannot determine model number from:\n"
            f"    {cif_file}"
        )

    return (
        Path.cwd()
        / (
            f"{output_prefix}"
            f"_model_{model_number}.maegz"
        )
    )


# ============================================================================
# PrepWizard command
# ============================================================================

def build_prepwizard_command(
    cif_file,
    output_file,
    reference,
    host,
    samplewater,
):
    """
    Build the PrepWizard command.
    """

    schrodinger = os.environ.get(
        "SCHRODINGER"
    )

    if not schrodinger:

        raise EnvironmentError(
            "SCHRODINGER environment variable "
            "is not set.\n"
            "Run this script using:\n"
            "    $SCHRODINGER/run "
            "boltz_prepwizard.py ..."
        )

    prepwizard = (
        Path(schrodinger)
        / "utilities"
        / "prepwizard"
    )

    if not prepwizard.is_file():

        raise FileNotFoundError(
            "PrepWizard executable not found:\n"
            f"    {prepwizard}"
        )

    jobname = (
        Path(output_file).stem
    )

    cmd = [
        str(prepwizard),

        str(cif_file),
        str(output_file),

        "-fillsidechains",
        "-disulfides",
        "-assign_all_residues",
        "-rehtreat",

        "-max_states",
        "1",

        "-epik_pH",
        "7.4",

        "-epik_pHt",
        "2.0",

        "-antibody_cdr_scheme",
        "Kabat",
    ]

    # --------------------------------------------------------------
    # Optional reference structure
    # --------------------------------------------------------------

    if reference is not None:

        cmd.extend([
            "-reference_st_file",
            str(reference),
        ])

    # --------------------------------------------------------------
    # Optional sample water
    # --------------------------------------------------------------

    if samplewater:

        cmd.append(
            "-samplewater"
        )

    cmd.extend([
        "-propka_pH",
        "7.4",

        "-f",
        "S-OPLS",

        "-rmsd",
        "0.3",

        "-watdist",
        "5.0",

        "-JOBNAME",
        jobname,

        "-HOST",
        host,
    ])

    return cmd


# ============================================================================
# Run one PrepWizard job
# ============================================================================

def run_prepwizard(
    cif_file,
    output_file,
    reference,
    host,
    samplewater,
):
    """
    Launch one PrepWizard job and WAIT for completion.

    Only after this function returns successfully
    is the next model submitted.

    This intentionally limits concurrent
    Schrödinger license-token usage.
    """

    cmd = build_prepwizard_command(
        cif_file=cif_file,
        output_file=output_file,
        reference=reference,
        host=host,
        samplewater=samplewater,
    )

    print()
    print("=" * 78)

    print(
        f"PrepWizard: "
        f"{Path(cif_file).name}"
    )

    print("=" * 78)

    print(
        f"Input : {cif_file}"
    )

    print(
        f"Output: {output_file}"
    )

    if reference is not None:

        print(
            f"Reference: {reference}"
        )

    print(
        f"Host: {host}"
    )

    print()
    print(
        "Command:"
    )

    print(
        " ".join(cmd)
    )

    print()
    print(
        "Launching PrepWizard ..."
    )

    # --------------------------------------------------------------
    # Launch
    # --------------------------------------------------------------

    job = jobcontrol.launch_job(
        cmd,
        print_output=True,
        launch_dir=str(
            Path.cwd()
        ),
        show_failure_dialog=False,
    )

    print(
        f"JobId: {job.JobId}"
    )

    print()
    print(
        "Waiting for completion "
        "(serial license-token control) ..."
    )

    # --------------------------------------------------------------
    # WAIT before the next model.
    # --------------------------------------------------------------

    job.wait(
        throw_on_failure=True
    )

    print(
        "PrepWizard completed successfully."
    )

    # --------------------------------------------------------------
    # Verify output
    # --------------------------------------------------------------

    if not Path(
        output_file
    ).is_file():

        raise RuntimeError(
            "PrepWizard reported success, "
            "but output MAEGZ was not found:\n"
            f"    {output_file}"
        )

    return job


# ============================================================================
# Main
# ============================================================================

def main():

    args = parse_args()

    # --------------------------------------------------------------
    # Input YAML
    # --------------------------------------------------------------

    input_yaml = (
        check_file_exists(
            args.input,
            "Boltz input YAML",
        )
    )

    input_stem = (
        get_input_stem(
            input_yaml
        )
    )

    # --------------------------------------------------------------
    # Output prefix
    # --------------------------------------------------------------

    if args.output_prefix is None:

        output_prefix = input_stem

    else:

        output_prefix = (
            args.output_prefix
        )

    # --------------------------------------------------------------
    # Reference
    # --------------------------------------------------------------

    reference = None

    if args.reference is not None:

        reference = (
            check_file_exists(
                args.reference,
                "Reference structure",
            )
        )

    # --------------------------------------------------------------
    # Boltz prediction directory
    # --------------------------------------------------------------

    prediction_dir = (
        find_prediction_directory(
            input_yaml
        )
    )

    # --------------------------------------------------------------
    # Find CIF files
    # --------------------------------------------------------------

    cif_files = (
        find_cif_files(
            prediction_dir
        )
    )

    if not cif_files:

        raise RuntimeError(
            "No Boltz model CIF files found in:\n"
            f"    {prediction_dir}\n\n"
            "Expected filenames such as:\n"
            f"    {input_stem}_model_0.cif\n"
            f"    {input_stem}_model_1.cif"
        )

    # --------------------------------------------------------------
    # Header
    # --------------------------------------------------------------

    print()
    print("#" * 78)
    print(
        "# boltz_prepwizard.py"
    )
    print("#" * 78)

    print(
        f"Input YAML       : {input_yaml}"
    )

    print(
        f"Input stem       : {input_stem}"
    )

    print(
        f"Prediction dir   : {prediction_dir}"
    )

    print(
        f"Output prefix    : {output_prefix}"
    )

    print(
        "Reference        : "
        f"{reference if reference else 'None'}"
    )

    print(
        "Renumber shift   : "
        f"{args.renumber_shift}"
    )

    print(
        f"PrepWizard host  : {args.host}"
    )

    print(
        "Sample water     : "
        f"{not args.no_samplewater}"
    )

    print(
        f"Models found     : {len(cif_files)}"
    )

    print()
    print(
        "Models:"
    )

    for cif_file in cif_files:

        model_number = (
            extract_model_number(
                cif_file
            )
        )

        confidence_file = (
            find_confidence_json(
                cif_file
            )
        )

        print(
            f"model_{model_number}: "
            f"{cif_file.name}"
        )

        print(
            "    confidence: "
            f"{confidence_file.name}"
        )

    print("#" * 78)

    # --------------------------------------------------------------
    # Serial processing
    # --------------------------------------------------------------

    successful = []
    failed = []

    for index, cif_file in enumerate(
        cif_files,
        start=1,
    ):

        model_number = (
            extract_model_number(
                cif_file
            )
        )

        output_file = (
            make_output_filename(
                cif_file,
                output_prefix,
            )
        )

        print()
        print()

        print(
            f"[{index}/{len(cif_files)}] "
            f"Starting model_{model_number}"
        )

        try:

            # ------------------------------------------------------
            # Read confidence JSON before consuming license token.
            # ------------------------------------------------------

            confidence_file = (
                find_confidence_json(
                    cif_file
                )
            )

            confidence_data = (
                read_confidence_json(
                    confidence_file
                )
            )

            # ------------------------------------------------------
            # PrepWizard
            # ------------------------------------------------------

            run_prepwizard(
                cif_file=cif_file,
                output_file=output_file,
                reference=reference,
                host=args.host,
                samplewater=(
                    not args.no_samplewater
                ),
            )

            # ------------------------------------------------------
            # Post-processing
            # ------------------------------------------------------

            print()
            print(
                f"Post-processing "
                f"model_{model_number} ..."
            )

            finalize_structure(
                maegz_file=output_file,
                model_number=model_number,
                confidence_data=confidence_data,
                renumber_shift=(
                    args.renumber_shift
                ),
            )

            successful.append(
                (
                    model_number,
                    output_file,
                )
            )

            print()
            print(
                f"model_{model_number} "
                "completed successfully."
            )

        except Exception as exc:

            failed.append(
                (
                    model_number,
                    cif_file,
                    str(exc),
                )
            )

            print()
            print("#" * 78)

            print(
                f"ERROR processing "
                f"model_{model_number}"
            )

            print("#" * 78)

            print(
                f"CIF   : {cif_file}"
            )

            print(
                f"Error : {exc}"
            )

            print()
            print(
                "Workflow stopped because "
                "this model failed."
            )

            break

    # --------------------------------------------------------------
    # Final summary
    # --------------------------------------------------------------

    print()
    print()
    print("#" * 78)
    print(
        "# FINAL SUMMARY"
    )
    print("#" * 78)

    print(
        f"Total models : {len(cif_files)}"
    )

    print(
        f"Successful   : {len(successful)}"
    )

    print(
        f"Failed       : {len(failed)}"
    )

    if successful:

        print()
        print(
            "Successful outputs:"
        )

        for (
            model_number,
            output_file,
        ) in successful:

            print(
                f"    model_{model_number}: "
                f"{output_file.name}"
            )

    if failed:

        print()
        print(
            "Failed models:"
        )

        for (
            model_number,
            cif_file,
            error,
        ) in failed:

            print(
                f"    model_{model_number}: "
                f"{cif_file.name}"
            )

            print(
                f"        {error}"
            )

        print()
        print(
            "Workflow finished with errors."
        )

        sys.exit(1)

    print()
    print(
        "All Boltz models were "
        "prepared successfully."
    )

    print("#" * 78)


if __name__ == "__main__":
    main()
