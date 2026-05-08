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

"""Build uDALES with CMake.

"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def build_from_branch(branch_name: str, path_to_proj_dir: Path, build_type: str, clean_build_dir=False, skip_build=False) -> str:
    subprocess.run(['git', 'checkout', branch_name], cwd=path_to_proj_dir, check=True)
    # Common branch names use / as user separator.
    path_to_build_dir = path_to_proj_dir / 'build' / branch_name.replace('/', '_')
    if not skip_build:
        build(path_to_proj_dir, path_to_build_dir, build_type, clean_build_dir=False)
    return path_to_build_dir


def build(path_to_proj_dir: Path, path_to_build_dir: Path, build_type: str, clean_build_dir=False) -> None:
    if clean_build_dir:
        shutil.rmtree(path_to_build_dir, ignore_errors=True)
    if not path_to_build_dir.is_dir():
        path_to_build_dir.mkdir(parents=True)

    cmake_configure_cmd = ['cmake', f'-DCMAKE_BUILD_TYPE={build_type}']
    cmake_args = os.environ.get('UDALES_CMAKE_ARGS', '').strip()
    if cmake_args:
        cmake_configure_cmd.extend(cmake_args.split())
    cmake_configure_cmd.extend([str(path_to_proj_dir), '-LA'])

    try:
        subprocess.run(
            cmake_configure_cmd,
            cwd=path_to_build_dir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        print(
            f"CMake configure failed in {path_to_build_dir} for {path_to_proj_dir} "
            f"({build_type}).",
            file=sys.stderr,
        )
        if exc.stdout:
            print(exc.stdout, file=sys.stderr, end="" if exc.stdout.endswith("\n") else "\n")
        raise

    cpu_count = str(os.cpu_count())
    try:
        subprocess.run(
            ['cmake', '--build', '.', '--', '-j', cpu_count],
            cwd=path_to_build_dir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        print(
            f"CMake build failed in {path_to_build_dir} for {path_to_proj_dir} "
            f"({build_type}).",
            file=sys.stderr,
        )
        if exc.stdout:
            print(exc.stdout, file=sys.stderr, end="" if exc.stdout.endswith("\n") else "\n")
        raise
    return None


if __name__ == "__main__":
    pass
