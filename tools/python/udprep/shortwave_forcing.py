"""Archive atmospheric shortwave inputs before they are mapped to facets."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np


def write_shortwave_forcing(
    path: Path,
    times: np.ndarray,
    dni: np.ndarray,
    dsky: np.ndarray,
    zenith: np.ndarray,
    azimuth_local: np.ndarray,
    *,
    source: str,
    start_time: str,
    ghi: np.ndarray | None = None,
    overwrite: bool = False,
) -> None:
    """Save the atmospheric series used for facet shortwave calculations.

    Values are written as float64 without clipping or recomputing the direct
    and diffuse components. Time is in seconds since the simulation start.
    """
    from netCDF4 import Dataset

    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Shortwave forcing file already exists: {path}")

    series = {
        "time": np.asarray(times, dtype=np.float64),
        "dni": np.asarray(dni, dtype=np.float64),
        "dsky": np.asarray(dsky, dtype=np.float64),
        "solar_zenith": np.asarray(zenith, dtype=np.float64),
        "solar_azimuth_local": np.asarray(azimuth_local, dtype=np.float64),
    }
    if ghi is not None:
        series["ghi"] = np.asarray(ghi, dtype=np.float64)
    count = series["time"].size
    if count == 0 or any(values.ndim != 1 or values.size != count or not np.isfinite(values).all()
                         for values in series.values()):
        raise ValueError("Shortwave forcing arrays must be nonempty, finite 1D series of equal length")
    if np.any(np.diff(series["time"]) <= 0):
        raise ValueError("Shortwave forcing times must be strictly increasing")

    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        with Dataset(temp_path, "w", format="NETCDF4") as dataset:
            dataset.createDimension("time", count)
            dataset.setncattr("source", source)
            dataset.setncattr("simulation_start", start_time)
            descriptions = {
                "time": ("Seconds since simulation start", "s"),
                "dni": ("Direct normal irradiance", "W m-2"),
                "dsky": ("Diffuse horizontal sky irradiance", "W m-2"),
                "solar_zenith": ("Solar zenith angle", "degree"),
                "solar_azimuth_local": ("Solar azimuth after subtracting xazimuth", "degree"),
                "ghi": ("Global horizontal irradiance", "W m-2"),
            }
            for name, values in series.items():
                variable = dataset.createVariable(name, "f8", ("time",))
                variable.long_name, variable.units = descriptions[name]
                variable[:] = values
        if path.exists() and not overwrite:
            raise FileExistsError(f"Shortwave forcing file already exists: {path}")
        os.replace(temp_path, path)
    finally:
        temp_path.unlink(missing_ok=True)
