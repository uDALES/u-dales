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
from pathlib import Path


def build_from_branch(branch_name: str, path_to_proj_dir: Path, build_type: str, clean_build_dir=False, skip_build=False) -> str:
    subprocess.run(['git', 'checkout', branch_name])
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

    subprocess.run(
        ['cmake', f'-DCMAKE_BUILD_TYPE={build_type}', path_to_proj_dir, '-LA'], cwd=path_to_build_dir)

    cpu_count = str(os.cpu_count())
    subprocess.run(['cmake', '--build', '.', '--', '-j', cpu_count], cwd=path_to_build_dir)
    return None


if __name__ == "__main__":
    pass
