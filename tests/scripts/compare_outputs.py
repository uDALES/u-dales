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


"""Compare two netCDF datasets using approximate error.

This module compares any two dataset of arbitrary dimension and
calculates their approximate error (i.e. A - B).
Running main will save a boxplot of approximate errors in the dataset
folder and print values of maximum approximate error and standard deviation
to standard output.

"""

from pathlib import Path

from typing import List

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


def compare(path_to_ds_a: Path, path_to_ds_b: Path, path_to_fig_dir: Path, quantities: list) -> None:
    diffs = calc_diff(path_to_ds_a, path_to_ds_b, quantities)
    plot(diffs, path_to_fig_dir, quantities)


def calc_diff(path_to_ds_a: Path, path_to_ds_b: Path, quantities: List[str]) -> List[nc.Dataset]:
    ds_a = nc.Dataset(path_to_ds_a)
    ds_b = nc.Dataset(path_to_ds_b)
    diffs = []
    for quantity in quantities:
        diff = (ds_a[quantity][:] - ds_b[quantity][:]).ravel()
        print(f'Max approximate error for {quantity}: {diff.max()}')
        print(
            f'Standard deviation of approximate error for {quantity}: {diff.std()}')
        diffs.append(diff)
    return diffs


def plot(diffs: List[nc.Dataset], path_to_fig_dir: Path, quantities:list) -> None:
    plt.boxplot(diffs)
    plt.xticks(list(range(1, len(quantities) +1 )), quantities)
    plt.ylabel('Approximation error')
    plt.xlabel('Quantity')
    plt.savefig(path_to_fig_dir / 'compare_outputs.png')


if __name__ == "__main__":
    pass
