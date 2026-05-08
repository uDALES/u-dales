#!/usr/bin/env python

# uDALES (https://github.com/uDALES/u-dales).

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

# Copyright (C) 2019 the uDALES Team.


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


import f90nml

from scripts import compare_outputs, build_model

PROJ_DIR = Path(__file__).resolve().parents[3]


def _current_head(path_to_repo: Path) -> str:
    result = subprocess.run(
        ['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
        cwd=path_to_repo,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    head = result.stdout.strip()
    if head != 'HEAD':
        return head

    detached = subprocess.run(
        ['git', 'rev-parse', 'HEAD'],
        cwd=path_to_repo,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return detached.stdout.strip()


def _repo_has_uncommitted_changes(path_to_repo: Path) -> bool:
    result = subprocess.run(
        ['git', 'status', '--porcelain'],
        cwd=path_to_repo,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return bool(result.stdout.strip())


def _require_clean_worktree_for_branch_switch(path_to_repo: Path, branch_a: str, branch_b: str) -> None:
    if branch_a == branch_b:
        return
    if not _repo_has_uncommitted_changes(path_to_repo):
        return
    raise RuntimeError(
        "Regression tests require switching the repository checkout between "
        f"'{branch_a}' and '{branch_b}'. The current worktree has uncommitted "
        "changes, so Git checkout would overwrite local edits. This is okay "
        "for normal development; commit or stash your current changes before "
        "running the supported regression path locally."
    )


def main(branch_a: str, branch_b: str, build_type: str):
    if platform.system() not in ['Linux', 'Darwin']:
        raise RuntimeError(
            f'The operating system {platform.system()} is not currently suppoorted.')

    if branch_a != branch_b and _repo_has_uncommitted_changes(PROJ_DIR):
        warnings.warn(
            "Skipping supported regression branch-switch harness because the "
            "worktree has uncommitted changes. This harness checks out both "
            "refs in-place, so local development edits would be overwritten. "
            "CI still runs this path from a clean checkout."
        )
        return
    original_head = _current_head(PROJ_DIR)

    try:
        if branch_a == branch_b:
            warnings.warn(
                'branch_a and branch_b are the same. Skipping regression tests')
            _ = build_model.build_from_branch(branch_a, PROJ_DIR, build_type)
            return

        # Build executables
        path_to_exes = []
        for branch in [branch_a, branch_b]:
            path_to_exe = build_model.build_from_branch(
                branch, PROJ_DIR, build_type, skip_build=False)
            # We always compare between two branches -- i.e. two executables.
            path_to_exes.append(path_to_exe)

        # Run model and store outputs - currently impossible for uDALES 2 because it requires different input files, and tests are run in master repo.
        #test_cases = (PROJ_DIR / 'tests' / 'regression' / 'cases').iterdir()
        #patched_example_cases = (PROJ_DIR / 'tests' / 'regression' / 'patches').iterdir()
        #run_and_compare(test_cases, path_to_exes, is_patch=False)
        #run_and_compare(patched_example_cases, path_to_exes, is_patch=True)
    finally:
        subprocess.run(['git', 'checkout', original_head], cwd=PROJ_DIR, check=True)


def run_and_compare(cases_dir, path_to_exes, is_patch=False):
    excluded_cases = ['501', '502']
    excluded_platforms = ['Darwin']
    precursor_sims = ['501']
    driver_sims = ['502']

    for case_path in sorted(cases_dir):
        case_id = case_path.stem

        if case_id in excluded_cases:
            if platform.system() in excluded_platforms:
                print(f'Skipping tests for case {case_id} on {excluded_platforms}')
                continue

        print(f'Running tests for example {case_id}')

        if is_patch:
            test_case_dir = PROJ_DIR / 'examples'/ case_id
        else:
            test_case_dir = case_path

        outputs_case_dir = PROJ_DIR / 'tests' / 'outputs' / test_case_dir.name
        # Always start afresh.
        shutil.rmtree(outputs_case_dir, ignore_errors=True)

        model_output_dirs = []
        for path_to_exe in path_to_exes:
            # Create path to out folder
            model_output_dir = outputs_case_dir / path_to_exe.name
            shutil.copytree(test_case_dir, model_output_dir)

            if is_patch:
                # Apply test namelist patches to examples to reduce runtime
                nml = model_output_dir / f'namoptions.{case_id}'
                nml_patch = f90nml.read(case_path)
                nml_patched = model_output_dir / f'namoptions.{case_id}.patch'
                f90nml.patch(nml, nml_patch, nml_patched)
                namelist = nml_patched.name
            else:
                 namelist = f'namoptions.{case_id}'

            cpu_count = get_mpi_process_count(model_output_dir / namelist)

            # For driver sims we need to copy all files in first from the precursor simulation.
            if case_id in driver_sims:
                for f_name in (model_output_dir.parents[1] / '501' / model_output_dir.name).glob('*driver*'):
                    shutil.copy(f_name, model_output_dir.parents[1] / '502' / model_output_dir.name)

            run_udales(path_to_exe, namelist, model_output_dir, model_output_dirs,
                       cpu_count=cpu_count)

        # We do not compare precursor sims
        if not case_id in precursor_sims:
            # TODO: concatenate filedumps?
            compare_outputs.compare(model_output_dirs[0] / f'fielddump.000.{test_case_dir.name}.nc',
                                    model_output_dirs[1] / f'fielddump.000.{test_case_dir.name}.nc',
                                    model_output_dirs[0].parent)

def get_mpi_process_count(namelist_path: Path) -> int:
    run_nml = f90nml.read(namelist_path).get('run', {})
    nprocx = int(run_nml.get('nprocx', 1))
    nprocy = int(run_nml.get('nprocy', 1))
    cpu_count = nprocx * nprocy

    if cpu_count > 4:
        raise RuntimeError(
            f'Refusing to run regression case with {cpu_count} MPI ranks '
            f'from {namelist_path.name}; limit is 4.')

    return cpu_count


def run_udales(path_to_exe: Path, namelist: str, model_output_dir: str,
               model_output_dirs: list, cpu_count: int) -> None:
    print(f'Running uDALES in: {path_to_exe}')
    try:
        subprocess.run(['mpiexec', '-np', str(cpu_count), path_to_exe / 'u-dales',
                        namelist], cwd=model_output_dir, check=True,
                        stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print(f'Could not run case uDALES in {path_to_exe} for namelist {namelist}')
        sys.exit(1)
    model_output_dirs.append(model_output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='TODO')
    parser.add_argument('branch_a', help='TODO')
    parser.add_argument('branch_b', help='TODO')
    parser.add_argument('build_type', help='TODO')
    args = parser.parse_args()
    main(args.branch_a, args.branch_b, args.build_type)
