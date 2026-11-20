"""Extract height-specific atmospheric planes from gathered uDALES statistics."""

from __future__ import annotations

from pathlib import Path
import os
import tempfile
from typing import TYPE_CHECKING

from netCDF4 import Dataset
import numpy as np

from .heights import bracketing_levels, configured_receptor_heights, height_tag, validate_heights
from .schedule import read_stats_schedule, validate_statistics_times
from .sources import STATISTICS_SOURCES, load_statistics

if TYPE_CHECKING:
    from udbase import UDBase


_FIELD_UNITS = {"ta": "K", "pabs": "Pa", "qv": "kg kg-1", "ws_local": "m s-1", "ws_10": "m s-1"}


class UDComfAtmosphere:
    """Use UDBase's gathered statistics loaders; never reconstruct mean wind."""

    def __init__(self, sim: UDBase):
        self.sim = sim

    def write_heights(
        self, *, source: str = "stats_t",
        heights: tuple[float, ...] | None = None,
        output_dir: Path | None = None, overwrite: bool = False,
    ) -> dict[float, Path]:
        """Write ta, pabs, qv and saved mean speeds for every selected height.

        ``source`` is ``stats_t``, ``stats_tkslice`` or the warm-start
        ``merged_stats_tkslice`` written by :func:`merge_statistics_outputs`.
        Every saved statistics record becomes one output record: the mean
        over the ``tstatsdump`` seconds preceding its timestamp.
        """
        if source not in STATISTICS_SOURCES:
            raise ValueError(f"source must be one of {list(STATISTICS_SOURCES)}")
        sim = self.sim
        configured = configured_receptor_heights(sim.path)
        selected = configured if heights is None else validate_heights(heights)
        schedule = read_stats_schedule(sim.path)

        def load(var, **kwargs):
            return load_statistics(sim, source, var, **kwargs)

        source_path = Path(sim.path) / f"{source}.{sim.expnr}.nc"
        if not source_path.is_file():
            raise FileNotFoundError(f"Statistics file is missing: {source_path}")

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
        brackets = bracketing_levels(z_model, selected)
        for height, (lower, upper, _) in brackets.items():
            if lower not in source_levels or upper not in source_levels:
                raise ValueError(f"Statistics lack bracketing scalar levels for {height:g} m")
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

        saved_times = validate_statistics_times(load("time"), schedule)
        records = np.arange(saved_times.size)
        record_times = saved_times[records]

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
                self._create_output(ds, x, y, record_times, schedule.window_bounds(record_times),
                                    height, masks[height], source, schedule.tstatsdump)

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
                       record_times: np.ndarray, bounds: np.ndarray,
                       height: float, mask: np.ndarray, source: str, duration: float) -> None:
        ds.createDimension("x", len(x))
        ds.createDimension("y", len(y))
        ds.createDimension("time", len(record_times))
        ds.createDimension("bounds", 2)
        ds.createVariable("x", "f8", ("x",))[:] = x
        ds.createVariable("y", "f8", ("y",))[:] = y
        ds["x"].units = ds["y"].units = "m"
        ds.createVariable("time", "f8", ("time",))[:] = record_times
        ds["time"].units = "s since simulation start"
        ds["time"].long_name = "end of the statistics averaging window"
        ds.createVariable("time_bounds", "f8", ("time", "bounds"))[:] = bounds
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
            variable.cell_methods = "time: mean"
        ds.source = source
        ds.setncattr("averaging_window_s", float(duration))
        ds.setncattr("window_note", "Each record is the solver's mean over the tstatsdump seconds "
                                    "preceding its timestamp")

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
