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

from pathlib import Path
import platform
import subprocess
import shutil
import argparse

import numpy as np
import matplotlib.pyplot as plt


from scripts import compare_outputs, build_model


PROJ_DIR = Path(__file__).resolve().parents[1]


def main(branch_a: str, branch_b: str, build_type: str):
    if platform.system() not in ['Linux', 'Darwin']:
        raise RuntimeError(
            f'The operating system {platform.system()} is not currently suppoorted.')

    # Build executables
    path_to_exes = []
    for branch in [branch_a, branch_b]:
        path_to_exe = build_model.build_from_branch(
            branch, PROJ_DIR, build_type)
        # We always compare between two branches -- i.e. two executables.
        path_to_exes.append(path_to_exe)

    # Run model and store outputs
    for case_dir in (PROJ_DIR / 'tests' / 'cases').iterdir():
        output_dirs = []
        for path_to_exe in path_to_exes:
            # Create path to out folder
            output_dir = PROJ_DIR / 'tests' / 'outputs' / case_dir.name / path_to_exe.name
            # Always start afresh. Remove all files and folders under the case dir.
            shutil.rmtree(output_dir.parent(), ignore_errors=True)
            shutil.copytree(case_dir, output_dir)
            namelist = "namoptions." + case_dir.name

            # Run model
            # FIXME: make num proc func of sys used.
            subprocess.run(['mpiexec', '-np', '2', path_to_exe / 'u-dales',
                            namelist], cwd=output_dir)
            output_dirs.append(output_dir)

        # FIXME: concatenate filedumps.
        compare_outputs.compare(output_dirs[0] / f'fielddump.001.{case_dir.name}.nc',
                                output_dirs[1] / f'fielddump.001.{case_dir.name}.nc',
                                output_dirs[0].parent)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='TODO')
    parser.add_argument('branch_a', help='TODO')
    parser.add_argument('branch_b', help='TODO')
    parser.add_argument('build_type', help='TODO')
    args = parser.parse_args()
    main(args.branch_a, args.branch_b, args.build_type)
