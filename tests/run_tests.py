#!/usr/bin/env python

# uDALES (https://github.com/uDALES/u-dales).
# Copyright (C) 2019 D. Meyer.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""Test uDALES

This program is used to compare the results obtained from executables
produced from two different branches.
"""

import os
import sys
from pathlib import Path
import platform
import subprocess
import shutil
import argparse
import warnings

import numpy as np
import matplotlib.pyplot as plt


from scripts import compare_outputs, build_model


PROJ_DIR = Path(__file__).resolve().parents[1]


def main(branch_a: str, branch_b: str, build_type: str):
    if platform.system() not in ['Linux', 'Darwin']:
        raise RuntimeError(
            f'The operating system {platform.system()} is not currently suppoorted.')

    if branch_a == branch_b:
        warnings.warn(
            'branch_a and branch_b are the same. Skipping regression tests')
        _ = build_model.build_from_branch(branch_a, PROJ_DIR, build_type)
        return

    # Build executables
    path_to_exes = []
    for branch in [branch_a, branch_b]:
        path_to_exe = build_model.build_from_branch(
            branch, PROJ_DIR, build_type)
        # We always compare between two branches -- i.e. two executables.
        path_to_exes.append(path_to_exe)

    # Run model and store outputs
    for test_case_dir in (PROJ_DIR / 'tests' / 'cases').iterdir():
        outputs_case_dir = PROJ_DIR / 'tests' / 'outputs' / test_case_dir.name
        # Always start afresh.
        shutil.rmtree(outputs_case_dir, ignore_errors=True)

        model_output_dirs = []
        for path_to_exe in path_to_exes:
            # Create path to out folder
            model_output_dir = outputs_case_dir / path_to_exe.name
            shutil.copytree(test_case_dir, model_output_dir)
            namelist = "namoptions." + test_case_dir.name

            run_udales(path_to_exe, namelist, model_output_dir, model_output_dirs)

            if test_case_dir.name == '101':
                # Run again with outputs from precursor simualton
                namelist = "namoptions.driven"
                run_udales(path_to_exe, namelist, model_output_dir, model_output_dirs)

        if test_case_dir.name == '001':
            # Scalar outputs are only outputted for this case
            quantities = ['u', 'v', 'w', 'sca1']
        else:
            quantities = ['u', 'v', 'w']
        # TODO: concatenate filedumps?
        compare_outputs.compare(model_output_dirs[0] / f'fielddump.000.{test_case_dir.name}.nc',
                                model_output_dirs[1] / f'fielddump.000.{test_case_dir.name}.nc',
                                model_output_dirs[0].parent, quantities)

def run_udales(path_to_exe: Path, namelist: str, model_output_dir: str, 
               model_output_dirs: list, cpu_count=None) -> None:
    if cpu_count is None:
        cpu_count = str(os.cpu_count())
    print(f'Running uDALES in: {path_to_exe}')
    try:
        subprocess.run(['mpiexec', '-np', cpu_count, path_to_exe / 'u-dales',
                        namelist], cwd=model_output_dir, check=True)
    except:
        print(f'Could not run case uDALES in {path_to_exe} for namelist {namelist}')
        sys.exit()
    model_output_dirs.append(model_output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='TODO')
    parser.add_argument('branch_a', help='TODO')
    parser.add_argument('branch_b', help='TODO')
    parser.add_argument('build_type', help='TODO')
    args = parser.parse_args()
    main(args.branch_a, args.branch_b, args.build_type)
