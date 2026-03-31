# Binaries

This directory contains the primary runtime executable and helper utilities.

## Contents

- `u-dales`: Main u-DALES executable.
- `ud_nml2json`: Convert namelist input to JSON (no mapping applied).
- `ud_nml2v3`: Apply mapping rules from `docs/schemas/nml_mapping.txt` directly to a namelist file (in-place by default).
- `ud_testenvironment`: Convenience script to set up the test environment.

## Usage Notes

- Most utilities assume `UD_TOPDIR` is set or are run from within the repo.
- `ud_nam2v3 <input.nml> [output.nml]` overwrites the input file if no output is provided.
