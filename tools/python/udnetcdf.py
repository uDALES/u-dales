"""NetCDF loading helpers for uDALES output files.

Reads a NetCDF file, reverses dimension order to the MATLAB (column-major)
convention, and returns either a requested variable as a numpy array
(load-and-close) or the whole dataset in memory (browse). Extracted from
UDBase so file IO is separable from case state; UDBase keeps a thin wrapper.
"""
from __future__ import annotations

from pathlib import Path
import operator
from typing import Optional, Sequence, Union

import numpy as np
import xarray as xr


def load_ncdata(
    filename: Path, var: Optional[str], *, time_index: Optional[int] = None,
    vertical_indices: Optional[Sequence[int]] = None,
) -> Union[xr.Dataset, np.ndarray]:
    """
    Helper method to load NetCDF data using xarray.
    
    Automatically reverses dimension order to match MATLAB conventions.
    When a specific variable is requested, returns it as a numpy array for
    memory efficiency with large datasets.
    
    NetCDF files use C-style (row-major) dimension ordering, while MATLAB
    uses Fortran-style (column-major) ordering. Due to this fundamental
    difference, MATLAB automatically reverses dimension order when reading
    NetCDF files. This method reverses dimensions to match MATLAB's behavior.
    
    Examples of transformations:
    - (time, zt) -> (zt, time)
    - (time, zt, yt, xm) -> (xm, yt, zt, time)
    - (time, fct) -> (fct, time)
    - (time, lyr, fct) -> (fct, lyr, time)
    
    Parameters
    ----------
    filename : Path
        Path to NetCDF file
    var : str, optional
        Variable to extract. If None, displays available variables.
    time_index : int, optional
        Select a single record along ``time`` before loading the variable.
    vertical_indices : sequence of int, optional
        Select zero-based ``zt`` levels before loading; avoids reading a full
        3-D statistics record when only a few planes are needed.
    
    Returns
    -------
    xarray.Dataset or numpy.ndarray
        Full dataset if var is None (for browsing), numpy array if var is specified.
    """
    if not filename.exists():
        raise FileNotFoundError(f"File not found: {filename}")

    # Open inside a context manager so the file handle is always released,
    # rather than kept alive implicitly by the returned object.
    if (time_index is not None or vertical_indices is not None) and var is None:
        raise ValueError("time_index and vertical_indices require a variable name")
    with xr.open_dataset(filename) as ds:
        if time_index is not None or vertical_indices is not None:
            if var not in ds:
                raise KeyError(f"Variable '{var}' not found in {filename.name}")
            data_var = ds[var]
            if time_index is not None:
                if "time" not in data_var.dims:
                    raise ValueError(f"Variable '{var}' has no time dimension")
                index = operator.index(time_index)
                if not 0 <= index < ds.sizes["time"]:
                    raise IndexError(f"time_index {index} outside 0..{ds.sizes['time'] - 1}")
                data_var = data_var.isel(time=index)
            if vertical_indices is not None:
                if "zt" not in data_var.dims:
                    raise ValueError(f"Variable '{var}' has no zt dimension")
                levels = [operator.index(value) for value in vertical_indices]
                if not levels or any(index < 0 or index >= ds.sizes["zt"] for index in levels):
                    raise IndexError("vertical_indices lie outside the zt dimension")
                data_var = data_var.isel(zt=levels)
            return data_var.transpose(*reversed(data_var.dims)).values
        # Transpose all data variables to match MATLAB's column-major
        # convention (reverse dimension order for variables with 2+ dims).
        transposed_vars = {}
        for var_name in ds.data_vars:
            data_var = ds[var_name]
            if len(data_var.dims) >= 2:
                transposed_vars[var_name] = data_var.transpose(*reversed(data_var.dims))
            else:
                transposed_vars[var_name] = data_var

        ds_transposed = xr.Dataset(transposed_vars, coords=ds.coords, attrs=ds.attrs)

        if var is None:
            # Browse mode: read into memory so the returned dataset owns no
            # open file handle (the whole file is loaded — intended for
            # interactive inspection, not bulk field reads).
            display_ncinfo(ds_transposed, filename.name)
            return ds_transposed.load()

        if var not in ds_transposed:
            raise KeyError(f"Variable '{var}' not found in {filename.name}")
        # Load-and-close: materialise the one requested array before the file
        # handle is closed on exit from the with-block.
        return ds_transposed[var].values

def display_ncinfo(ds: xr.Dataset, filename: str):
    """Display information about NetCDF dataset."""
    print(f"\nContents of {filename}:")
    print(f"{'Name':<20} {'Dimensions':<30} {'Shape':<20}")
    print("-" * 70)
    
    for var_name in sorted(ds.data_vars):
        var = ds[var_name]
        dims = ', '.join(var.dims)
        shape = ' x '.join(map(str, var.shape))
        print(f"{var_name:<20} {dims:<30} {shape:<20}")
