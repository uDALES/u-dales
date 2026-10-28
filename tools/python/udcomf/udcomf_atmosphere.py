"""Extract height-specific atmospheric planes from gathered uDALES statistics."""

from __future__ import annotations

from pathlib import Path
import os
import tempfile
from typing import TYPE_CHECKING

import f90nml
from netCDF4 import Dataset
import numpy as np

from .heights import configured_receptor_heights, height_tag, validate_heights

if TYPE_CHECKING:
    from udbase import UDBase


_FIELD_UNITS = {"ta": "K", "pabs": "Pa", "qv": "kg kg-1", "ws_local": "m s-1", "ws_10": "m s-1"}


class UDComfAtmosphere:
    """Use UDBase's gathered statistics loaders; never reconstruct mean wind."""

    def __init__(self, sim: UDBase):
        self.sim = sim

    def write_hourly_heights(
        self, target_times: np.ndarray, *, source: str = "stats_kslice",
        heights: tuple[float, ...] | None = None,
        output_dir: Path | None = None, overwrite: bool = False,
        time_tolerance: float = 1.0,
    ) -> dict[float, Path]:
        """Write ta, pabs, qv and saved mean speeds for every selected height.

        `target_times` are seconds since simulation start. A restart inside an
        averaging window cannot be detected from stats files alone.
        """
        if source not in ("stats_t", "stats_kslice"):
            raise ValueError("source must be stats_t or stats_kslice")
        sim = self.sim
        configured = configured_receptor_heights(sim.path)
        selected = configured if heights is None else validate_heights(heights)
        output_options = f90nml.read(
            Path(sim.path) / f"namoptions.{sim.expnr}"
        ).get("output", {})
        if not np.isclose(float(output_options.get("tstatsdump", float("nan"))),
                          900.0, rtol=0, atol=1e-6):
            raise ValueError("tstatsdump must be 900 s for preceding-15-minute atmospheric means")
        load = sim.load_stat_t if source == "stats_t" else sim.load_stat_kslice
        source_path = Path(sim.path) / f"{source}.{sim.expnr}.nc"
        if not source_path.is_file():
            raise FileNotFoundError(f"Gathered statistics file is missing: {source_path}")

        saved_heights = np.atleast_1d(np.asarray(load("receptor_height"), dtype=float))
        height_tolerance = 100 * np.finfo(np.float32).eps * np.maximum(1.0, np.abs(configured))
        if saved_heights.shape != (len(configured),) or np.any(
            abs(saved_heights - configured) > height_tolerance
        ):
            raise ValueError("Saved ws_local receptor heights disagree with namoptions")
        height_indices = {}
        for height in selected:
            tolerance = 100 * np.finfo(np.float32).eps * max(1.0, abs(height))
            matches = np.flatnonzero(np.isclose(saved_heights, height, rtol=0, atol=tolerance))
            if len(matches) != 1:
                raise ValueError(f"No unique saved ws_local plane for height {height:g} m")
            height_indices[height] = int(matches[0])

        z_model = np.asarray(sim.zt, dtype=float)
        z_saved = np.asarray(load("zt"), dtype=float)
        if z_model.ndim != 1 or z_saved.ndim != 1 or np.any(np.diff(z_model) <= 0):
            raise ValueError("Statistics or model scalar-height coordinates are invalid")
        source_levels = {}
        for source_index, z in enumerate(z_saved):
            matches = np.flatnonzero(np.isclose(z_model, z, rtol=1e-6, atol=1e-4))
            if len(matches) != 1 or int(matches[0]) in source_levels:
                raise ValueError("Statistics zt coordinates do not uniquely match the model grid")
            source_levels[int(matches[0])] = source_index
        brackets = {}
        for height in selected:
            upper = int(np.searchsorted(z_model, height))
            if upper == len(z_model) or height < z_model[0]:
                raise ValueError(f"Height {height:g} m is outside scalar-centre levels")
            if np.isclose(z_model[upper], height, rtol=0, atol=1e-5):
                lower = upper
                weight = 0.0
            else:
                lower = upper - 1
                weight = (height - z_model[lower]) / (z_model[upper] - z_model[lower])
            if lower not in source_levels or upper not in source_levels:
                raise ValueError(f"Statistics lack bracketing scalar levels for {height:g} m")
            brackets[height] = (lower, upper, weight)
        needed = sorted({source_levels[k] for lo, hi, _ in brackets.values() for k in (lo, hi)})
        positions = {source_index: offset for offset, source_index in enumerate(needed)}

        x = np.asarray(sim.xt, dtype=float)
        y = np.asarray(sim.yt, dtype=float)
        saved_x = np.asarray(load("xt"), dtype=float)
        saved_y = np.asarray(load("yt"), dtype=float)
        if (saved_x.shape != x.shape or saved_y.shape != y.shape
                or not np.allclose(saved_x, x, rtol=0, atol=1e-4)
                or not np.allclose(saved_y, y, rtol=0, atol=1e-4)):
            raise ValueError("Statistics x/y coordinates disagree with the model grid")
        solid = getattr(sim, "Sc", None)
        if solid is None or solid.shape != (len(x), len(y), len(z_model)):
            raise ValueError("The scalar solid mask is required for atmospheric planes")
        masks = {}
        for height, (lo, hi, _) in brackets.items():
            receptor = sim.comf.radiation.receptor_grid(height)
            masks[height] = receptor.valid & ~solid[:, :, lo] & ~solid[:, :, hi]

        times = np.asarray(load("time"), dtype=float)
        targets = np.asarray(target_times, dtype=float)
        if (times.ndim != 1 or targets.ndim != 1 or not len(times) or not len(targets)
                or not np.isfinite(times).all() or not np.isfinite(targets).all()
                or np.any(np.diff(times) <= 0) or np.any(np.diff(targets) <= 0)
                or not np.isfinite(time_tolerance) or time_tolerance < 0):
            raise ValueError("Statistics and requested times must be finite and strictly increasing")
        records = np.searchsorted(times, targets)
        records = np.minimum(records, len(times) - 1)
        previous = np.maximum(records - 1, 0)
        records = np.where(abs(times[previous] - targets) < abs(times[records] - targets), previous, records)
        if np.any(abs(times[records] - targets) > time_tolerance) or len(set(records)) != len(records):
            raise ValueError("Statistics do not contain a unique record at each requested hour")

        directory = Path(sim.path) if output_dir is None else Path(output_dir)
        paths = {
            height: directory / f"pedestrian_atmosphere.{height_tag(height)}.{sim.expnr}.nc"
            for height in selected
        }
        if not overwrite:
            for path in paths.values():
                if path.exists():
                    raise FileExistsError(f"Atmospheric output exists: {path}")
        directory.mkdir(parents=True, exist_ok=True)

        temporary: dict[float, Path] = {}
        datasets: dict[float, Dataset] = {}
        try:
            for height, path in paths.items():
                fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=directory)
                os.close(fd)
                temporary[height] = Path(name)
                ds = Dataset(name, "w", format="NETCDF4")
                datasets[height] = ds
                self._create_output(ds, x, y, targets, times[records], height, masks[height], source)

            def interpolate(values: np.ndarray, height: float) -> np.ndarray:
                lo, hi, weight = brackets[height]
                lower = values[:, :, positions[source_levels[lo]]]
                upper = values[:, :, positions[source_levels[hi]]]
                return (1.0 - weight) * lower + weight * upper

            for output_index, record in enumerate(records):
                ws_local = np.asarray(load("ws_local", time_index=int(record)), dtype=float)
                ws_10 = np.asarray(load("ws_10", time_index=int(record)), dtype=float)
                for height, ds in datasets.items():
                    wind = ws_local if ws_local.ndim == 2 else ws_local[:, :, height_indices[height]]
                    self._write_field(ds, "ws_local", wind, masks[height], output_index, nonnegative=True)
                    self._write_field(ds, "ws_10", ws_10, masks[height], output_index, nonnegative=True)

                for source_name, output_name in (("tha", "ta"), ("pabs", "pabs")):
                    values = np.asarray(
                        load(source_name, time_index=int(record), vertical_indices=needed), dtype=float
                    )
                    for height, ds in datasets.items():
                        self._write_field(ds, output_name, interpolate(values, height),
                                          masks[height], output_index, positive=True)
                qt = np.asarray(load("qt", time_index=int(record), vertical_indices=needed), dtype=float)
                ql = np.asarray(load("ql", time_index=int(record), vertical_indices=needed), dtype=float)
                for height, ds in datasets.items():
                    qv = interpolate(qt, height) - interpolate(ql, height)
                    self._write_field(ds, "qv", qv, masks[height],
                                      output_index, nonnegative=True)
            for ds in datasets.values():
                ds.close()
            datasets.clear()
            for height, path in paths.items():
                temporary[height].replace(path)
        finally:
            for ds in datasets.values():
                ds.close()
            for path in temporary.values():
                path.unlink(missing_ok=True)
        return paths

    @staticmethod
    def _create_output(ds: Dataset, x: np.ndarray, y: np.ndarray,
                       targets: np.ndarray, source_times: np.ndarray,
                       height: float, mask: np.ndarray, source: str) -> None:
        ds.createDimension("x", len(x))
        ds.createDimension("y", len(y))
        ds.createDimension("time", len(targets))
        ds.createDimension("bounds", 2)
        ds.createVariable("x", "f8", ("x",))[:] = x
        ds.createVariable("y", "f8", ("y",))[:] = y
        ds["x"].units = ds["y"].units = "m"
        ds.createVariable("time", "f8", ("time",))[:] = targets
        ds["time"].units = "s since simulation start"
        ds.createVariable("source_time", "f8", ("time",))[:] = source_times
        ds["source_time"].units = "s since simulation start"
        ds.createVariable("time_bounds", "f8", ("time", "bounds"))[:] = np.column_stack(
            (targets - 900.0, targets)
        )
        ds["time"].bounds = "time_bounds"
        ds["time_bounds"].units = "s since simulation start"
        ds.createVariable("receptor_height", "f8").assignValue(height)
        ds["receptor_height"].units = "m above model ground"
        ds.createVariable("valid_mask", "i1", ("x", "y"))[:] = mask.astype(np.int8)
        for name, units in _FIELD_UNITS.items():
            variable = ds.createVariable(
                name, "f4", ("x", "y", "time"), fill_value=np.float32(np.nan),
                zlib=True, complevel=3, chunksizes=(min(len(x), 64), min(len(y), 64), 1),
            )
            variable.units = units
        ds.source = source
        ds.setncattr("window_note", "Requested preceding 900 s mean; verify warm-start boundaries separately")

    @staticmethod
    def _write_field(ds: Dataset, name: str, values: np.ndarray, mask: np.ndarray,
                     index: int, *, positive: bool = False, nonnegative: bool = False) -> None:
        values = np.asarray(values, dtype=float)
        if values.shape != mask.shape or not np.isfinite(values[mask]).all():
            raise ValueError(f"{name} contains missing or nonfinite values in fluid receptors")
        if positive and np.any(values[mask] <= 0):
            raise ValueError(f"{name} must be positive")
        if nonnegative and np.any(values[mask] < 0):
            raise ValueError(f"{name} must be nonnegative")
        ds[name][:, :, index] = np.where(mask, values, np.nan).astype(np.float32)
