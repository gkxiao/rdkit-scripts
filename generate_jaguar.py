#!/usr/bin/env python
import os
import sys
import subprocess

def print_help():
    help_text = """
================================================================================
===  Jaguar Input Generator for r2SCAN-3c / COSMO Solvent  ===
================================================================================
Function:
  1. Generate Jaguar .in file from SDF via Schrodinger obabel
  2. Predefined settings: r2SCAN-3c, COSMO solvent(isolv=7), nogas=2
  3. Set molecular charge by command line; multip=1 fixed as default
  4. Generate executable bash script for Jaguar calculation

Usage:
  python generate_jaguar.py molecule.sdf molecular_charge
  python generate_jaguar.py -h | --help

Example:
  python generate_jaguar.py test.sdf 1     # charge +1
  python generate_jaguar.py test.sdf 0     # neutral
  python generate_jaguar.py test.sdf -1    # charge -1

Output files:
  1. jag_${name}_spe_r2SCAN-3c.in    Jaguar input file
  2. jag_${name}_spe_r2SCAN-3c.sh    Job submission script
================================================================================
"""
    print(help_text)
    sys.exit(0)

def main():
    # Show help if -h / --help
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print_help()

    # Check command line arguments
    if len(sys.argv) != 3:
        print("Error: Invalid arguments!")
        print("Usage: python generate_jaguar.py molecule.sdf molecular_charge")
        print("Example: python generate_jaguar.py test.sdf 1")
        sys.exit(1)

    # Parse input arguments
    sdf_filename = sys.argv[1]
    mol_charge = sys.argv[2]
    multip_value = 1      # fixed spin multiplicity: multip=1

    # Get base name without .sdf suffix
    base_name = os.path.splitext(sdf_filename)[0]
    in_file = f"jag_{base_name}_spe_r2SCAN-3c.in"
    sh_file = f"jag_{base_name}_spe_r2SCAN-3c.sh"

    # Get SCHRODINGER environment variable
    schrodinger_env = os.getenv("SCHRODINGER")
    if not schrodinger_env:
        print("Error: Environment variable $SCHRODINGER not found!")
        sys.exit(1)

    # Path to Schrodinger obabel
    obabel_path = os.path.join(schrodinger_env, "utilities", "obabel")

    # Convert SDF to Jaguar .in file using obabel
    print(f"Generating Jaguar input file: {in_file}")
    subprocess.run(
        [obabel_path, "-isdf", sdf_filename, "-ojin", "-O", in_file],
        check=True
    )

    # Modify Jaguar input file header & gen block
    with open(in_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    new_content = []
    skip_old_header = True

    for line in lines:
        # Skip original header until &zmat appears
        if skip_old_header:
            if "&zmat" in line.strip():
                # Write customized &gen block
                new_content.append("&gen\n")
                new_content.append("isolv=7\n")
                new_content.append("pcm_model=cosmo\n")
                new_content.append("dftname=r2SCAN-3c\n")
                new_content.append("nogas=2\n")
                new_content.append(f"molchg={mol_charge}\n")
                new_content.append(f"multip={multip_value}\n")
                new_content.append("&\n")
                new_content.append(f"entry_name: {base_name}\n")
                new_content.append("&zmat\n")
                skip_old_header = False
            continue
        # Keep all coordinates and remaining lines
        new_content.append(line)

    # Write modified content back to .in file
    with open(in_file, "w", encoding="utf-8") as f:
        f.writelines(new_content)

    # Generate submission script, keep ${SCHRODINGER} raw without escape backslash
    sh_content = f"${{SCHRODINGER}}/jaguar run -jobname=jag_{base_name}_spe_r2SCAN-3c {in_file} -HOST localhost -PARALLEL 1 -max_threads 32 -TMPLAUNCHDIR\n"
    with open(sh_file, "w", encoding="utf-8") as f:
        f.write(sh_content)

    # Add execute permission to shell script
    os.chmod(sh_file, 0o755)

    print("Task completed successfully!")
    print(f"Settings applied: molchg={mol_charge}, multip={multip_value}")
    print(f"Run calculation with: bash {sh_file}")

if __name__ == "__main__":
    main()
