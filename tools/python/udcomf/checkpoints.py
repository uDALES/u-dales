"""Gap-aware integration and atomic pedestrian-radiation checkpoints."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import tempfile

from netCDF4 import Dataset
import numpy as np


@dataclass(frozen=True)
class Window:
    end: float
    indices: np.ndarray
    weights: np.ndarray

    @property
    def complete(self) -> bool:
        return bool(self.indices.size)


def hourly_windows(source_times: np.ndarray, target_times: np.ndarray,
                   *, duration: float = 900.0, max_gap: float | None = None) -> tuple[list[Window], float]:
    """Integrate piecewise-linear source samples over preceding windows.

    A window is incomplete if its endpoints are not bracketed or any source
    interval inside it exceeds max_gap. No extrapolation or gap filling.
    """
    source = np.asarray(source_times, dtype=float)
    targets = np.asarray(target_times, dtype=float)
    if (source.ndim != 1 or source.size < 2 or not np.isfinite(source).all()
            or np.any(np.diff(source) <= 0)):
        raise ValueError("Source times must be finite, increasing, and contain two records")
    if (targets.ndim != 1 or targets.size == 0 or not np.isfinite(targets).all()
            or np.any(np.diff(targets) <= 0)
            or not np.allclose(targets / 3600.0, np.round(targets / 3600.0), rtol=0, atol=1e-9)):
        raise ValueError("Target times must be increasing exact-hour endpoints")
    if not np.isfinite(duration) or duration <= 0 or duration > 3600:
        raise ValueError("Window duration must be in (0, 3600] seconds")
    intervals = np.diff(source)
    if max_gap is None:
        max_gap = 1.1 * float(np.median(intervals))
    if not np.isfinite(max_gap) or max_gap <= 0:
        raise ValueError("max_gap must be finite and positive")

    windows = []
    for end in targets:
        start = float(end - duration)
        if start < source[0] or end > source[-1]:
            windows.append(Window(float(end), np.empty(0, dtype=int), np.empty(0)))
            continue
        segments = np.flatnonzero((source[:-1] < end) & (source[1:] > start))
        if segments.size == 0 or np.any(intervals[segments] > max_gap):
            windows.append(Window(float(end), np.empty(0, dtype=int), np.empty(0)))
            continue
        weights = np.zeros(source.size, dtype=float)
        for i in segments:
            a = max(start, source[i])
            b = min(float(end), source[i + 1])
            fraction = (0.5 * (a + b) - source[i]) / intervals[i]
            weights[i] += (b - a) * (1.0 - fraction) / duration
            weights[i + 1] += (b - a) * fraction / duration
        if not np.isclose(weights.sum(), 1.0, rtol=0, atol=1e-10):
            raise ValueError("Source times do not fully cover an averaging window")
        indices = np.flatnonzero(weights)
        windows.append(Window(float(end), indices, weights[indices]))
    return windows, float(max_gap)


def file_signature(path: Path) -> dict:
    path = Path(path).resolve()
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def selection_hash(indices: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(indices, dtype="<i8").tobytes()).hexdigest()


class RadiationCheckpoints:
    """One atomic checkpoint per output hour and receptor tile."""

    def __init__(self, directory: Path, manifest: dict, *, resume: bool):
        self.directory = Path(directory)
        self.manifest = manifest
        path = self.directory / "manifest.json"
        if path.exists():
            if not resume:
                raise FileExistsError(f"Radiation checkpoints exist: {self.directory}")
            with path.open("r", encoding="ascii") as handle:
                if json.load(handle) != manifest:
                    raise ValueError("Checkpoint inputs or settings changed; use a new directory")
        else:
            if self.directory.exists() and any(self.directory.iterdir()):
                raise ValueError(f"Checkpoint directory has no manifest: {self.directory}")
            self.directory.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(prefix=".manifest.", suffix=".tmp", dir=self.directory)
            try:
                with os.fdopen(fd, "w", encoding="ascii") as handle:
                    json.dump(manifest, handle, sort_keys=True)
                os.replace(name, path)
            finally:
                Path(name).unlink(missing_ok=True)

    def path(self, time_index: int, tile_index: int) -> Path:
        return self.directory / f"time_{time_index:04d}_tile_{tile_index:05d}.npz"

    def read(self, time_index: int, tile_index: int, shape: tuple[int, int]) -> np.ndarray | None:
        path = self.path(time_index, tile_index)
        if not path.exists():
            return None
        try:
            with np.load(path, allow_pickle=False) as archive:
                values = archive["values"]
            if values.shape != shape or not np.all(np.isfinite(values) | np.isnan(values)):
                raise ValueError(f"Invalid checkpoint data: {path}")
            return values
        except (OSError, KeyError, ValueError) as exc:
            raise ValueError(f"Corrupt checkpoint: {path}") from exc

    def write(self, time_index: int, tile_index: int, values: np.ndarray) -> None:
        path = self.path(time_index, tile_index)
        fd, name = tempfile.mkstemp(prefix=".radiation.", suffix=".npz", dir=self.directory)
        os.close(fd)
        try:
            np.savez_compressed(name, values=np.asarray(values, dtype=np.float32))
            os.replace(name, path)
        finally:
            Path(name).unlink(missing_ok=True)


def consolidate_radiation(path: Path, checkpoints: RadiationCheckpoints,
                          grid, indices: np.ndarray, names: tuple[str, ...],
                          windows: list[Window], tile_size: int, *, overwrite: bool) -> Path:
    """Write native-grid (x,y,time) fields after all tile checkpoints exist."""
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Radiation output exists: {path}")
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temp_path = Path(name)
    try:
        with Dataset(temp_path, "w", format="NETCDF4") as ds:
            nx, ny = grid.valid.shape
            ds.createDimension("x", nx)
            ds.createDimension("y", ny)
            ds.createDimension("time", len(windows))
            ds.createDimension("bounds", 2)
            x = ds.createVariable("x", "f8", ("x",))
            y = ds.createVariable("y", "f8", ("y",))
            x[:] = grid.x
            y[:] = grid.y
            x.units = y.units = "m"
            time = ds.createVariable("time", "f8", ("time",))
            time[:] = [window.end for window in windows]
            time.units = "s since simulation start"
            time.bounds = "time_bounds"
            bounds = ds.createVariable("time_bounds", "f8", ("time", "bounds"))
            bounds[:] = [[window.end - 900.0, window.end] for window in windows]
            complete = ds.createVariable("time_complete", "i1", ("time",))
            complete[:] = [int(window.complete) for window in windows]
            ds.createVariable("valid_mask", "i1", ("x", "y"))[:] = grid.valid.astype(np.int8)
            processed = np.zeros(grid.valid.shape, dtype=np.int8)
            processed.ravel()[indices] = 1
            ds.createVariable("processed_mask", "i1", ("x", "y"))[:] = processed
            height = ds.createVariable("receptor_height", "f8")
            height.assignValue(float(grid.z.flat[0]))
            height.units = "m above model ground"
            variables = {}
            for field_name in names:
                var = ds.createVariable(
                    field_name, "f4", ("x", "y", "time"), fill_value=np.float32(np.nan),
                    zlib=True, complevel=3, chunksizes=(min(nx, 64), min(ny, 64), 1),
                )
                var.units = "W m-2"
                variables[field_name] = var
            ds.setncattr("integration", "Piecewise-linear source interpolation; exact 900 s integral")
            ds.setncattr("missing_windows", "NaN where source coverage is absent or a cadence gap is detected")
            ds.setncattr("provenance", json.dumps(checkpoints.manifest, sort_keys=True))
            for time_index, _window in enumerate(windows):
                planes = np.full((len(names), nx * ny), np.nan, dtype=np.float32)
                for tile_index, offset in enumerate(range(0, len(indices), tile_size)):
                    tile = indices[offset:offset + tile_size]
                    values = checkpoints.read(time_index, tile_index, (len(tile), len(names)))
                    if values is None:
                        raise ValueError(f"Missing checkpoint: time {time_index}, tile {tile_index}")
                    planes[:, tile] = values.T
                for field_name, plane in zip(names, planes):
                    variables[field_name][:, :, time_index] = plane.reshape(nx, ny)
        if path.exists() and not overwrite:
            raise FileExistsError(f"Radiation output appeared during consolidation: {path}")
        os.replace(temp_path, path)
    finally:
        temp_path.unlink(missing_ok=True)
    return path
