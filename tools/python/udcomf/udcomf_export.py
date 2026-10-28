"""Preparation checks for model-neutral thermal-comfort export.

The checks use the namelist parser shared with UDBase and only load the
one-dimensional profile. They do not initialise the Paris geometry or alter
simulation inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from contextlib import ExitStack
import os
import tempfile

import numpy as np
import f90nml
from netCDF4 import Dataset

from udconfig import parse_namoptions
from udprep.solar import solar_position_python
from .heights import configured_receptor_heights, height_tag, validate_heights


@dataclass(frozen=True)
class ComfortPreflight:
    case_dir: Path
    receptor_height_m: float
    receptor_heights_m: tuple[float, ...]
    first_scalar_height_m: float
    analysis_times_s: tuple[float, ...]
    estimated_stats_bytes: int
    blockers: tuple[str, ...]
    cautions: tuple[str, ...]

    @property
    def ready(self) -> bool:
        return not self.blockers


def check_comfort_preflight(case_dir: Path, analysis_day: date) -> ComfortPreflight:
    """Audit full-3D or gathered k-slice statistics before a simulation.

    The byte count is a lower-bound raw payload estimate for 28 3-D float32
    variables with temperature and moisture enabled. It excludes metadata,
    NetCDF overhead, slices, probes, restart files, and radiation output.
    """
    case = Path(case_dir).expanduser().resolve()
    options_path = case / f"namoptions.{case.name}"
    options = {key.lower(): value for key, value in parse_namoptions(options_path).items()}
    profile = np.loadtxt(case / f"prof.inp.{case.name}", comments="#", ndmin=2)
    if profile.ndim != 2 or profile.shape[1] < 1 or not np.isfinite(profile[:, 0]).all():
        raise ValueError("Scalar-height profile is missing or invalid")
    heights = np.asarray(profile[:, 0], dtype=float)
    if np.any(np.diff(heights) <= 0):
        raise ValueError("Scalar-height profile must increase strictly")

    blockers: list[str] = []
    cautions: list[str] = []
    try:
        start = datetime(
            int(options["year"]), int(options["month"]), int(options["day"]),
            int(options["hour"]), int(options.get("minute", 0)),
            int(options.get("second", 0)), tzinfo=timezone.utc,
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("A valid UTC simulation start date is required in namoptions") from exc
    if float(options.get("timezone", 0.0)) != 0.0:
        blockers.append("Simulation timezone must be UTC for the agreed output timeline")
    midnight = datetime.combine(analysis_day, datetime.min.time(), tzinfo=timezone.utc)
    targets = tuple(
        (midnight + timedelta(hours=i) - start).total_seconds()
        for i in range(1, 25)
    )
    if targets[0] <= 0:
        blockers.append("Analysis day does not follow the simulation start")

    try:
        receptor_heights = configured_receptor_heights(case)
    except (TypeError, ValueError) as exc:
        blockers.append(str(exc))
        receptor_heights = ()
    for height in receptor_heights:
        if not heights[0] <= height <= heights[-1]:
            blockers.append(
                f"receptor_height={height:g} m is outside scalar-centre levels "
                f"[{heights[0]:g}, {heights[-1]:g}] m; ws_local cannot be interpolated"
            )
    if not heights[0] <= 10.0 <= heights[-1]:
        blockers.append("10 m is outside scalar-centre levels; ws_10 cannot be interpolated")
    full_stats = bool(options.get("ltdump", False))
    slice_stats = bool(options.get("ltkslicedump", False))
    if not (full_stats or slice_stats):
        blockers.append("Enable ltdump or ltkslicedump to save atmospheric fields")
    if slice_stats and not full_stats:
        output = f90nml.read(options_path).get("output", {})
        count = int(output.get("nkslice", 0))
        raw = output.get("kslice", [])
        indices = raw if isinstance(raw, list) else [raw]
        indices = set(int(index) for index in indices[:count] if index is not None)
        if count <= 0 or len(indices) != count or any(index < 1 or index > len(heights) for index in indices):
            blockers.append("ltkslicedump needs valid, distinct kslice indices")
        else:
            for height in receptor_heights:
                if not heights[0] <= height <= heights[-1]:
                    continue
                upper = int(np.searchsorted(heights, height, side="left"))
                needed = {upper + 1}
                if not np.isclose(heights[upper], height, rtol=0, atol=1e-5):
                    needed.add(upper)
                if not needed.issubset(indices):
                    blockers.append(
                        f"kslice is missing scalar levels {sorted(needed)} for receptor_height={height:g} m"
                    )
    if not options.get("ltempeq", False):
        blockers.append("ltempeq is disabled: tha cannot be saved")
    if not options.get("lmoist", False):
        blockers.append("lmoist is disabled: qt and ql cannot be saved")
    if not options.get("leb", False) or not options.get("lwriteebfiles", False):
        blockers.append("Facet energy balance output is required for pedestrian longwave")
    if not options.get("ltimedepsw", False) or not options.get("ltimedeplw", False):
        blockers.append("Time-dependent shortwave and longwave forcing must be enabled")

    dump = float(options.get("tstatsdump", 0))
    gap = float(options.get("tstatsgap", 0))
    stat_start = float(options.get("tstatstart", 0))
    period = dump + gap
    if not np.isclose(dump, 900.0, rtol=0, atol=1e-6):
        blockers.append("tstatsdump must be 900 s for preceding-15-minute means")
    if not np.isclose(period, 3600.0, rtol=0, atol=1e-6):
        blockers.append("tstatsdump + tstatsgap must be 3600 s for hourly means")
    if not np.isclose(stat_start + period, targets[0], rtol=0, atol=1e-6):
        blockers.append(
            f"First stats dump would be at {stat_start + period:g} s; "
            f"the first required hour is {targets[0]:g} s"
        )
    if float(options.get("runtime", 0)) < targets[-1]:
        blockers.append("Configured runtime ends before the last required hourly plane")
    sample = float(options.get("tsample", 0))
    if not np.isfinite(sample) or sample <= 0 or sample > 900:
        blockers.append("tsample must be finite, positive, and no longer than 900 s")
    if sample > 5:
        cautions.append("tsample exceeds 5 s; 15-minute means may under-sample rapid changes")

    itot = int(options["itot"])
    jtot = int(options["jtot"])
    ktot = int(options["ktot"])
    if len(heights) != ktot:
        blockers.append(f"prof has {len(heights)} scalar levels but ktot={ktot}")
    vertical_count = ktot if full_stats else int(options.get("nkslice", 0))
    estimated = itot * jtot * vertical_count * 28 * 4 * 24
    if full_stats:
        cautions.append(
            "Full stats_t writes at least 28 three-dimensional fields per dump; "
            f"24 dumps are roughly {estimated / 1e12:.2f} TB of raw float32 payload"
        )
    elif slice_stats:
        cautions.append("K-slice-only output saves disk space but still allocates 3D statistics in the solver")
    cautions.append(
        "Warm starts inside a 15-minute window can yield partial means even if dump "
        "timestamps match; check restart times against all requested windows"
    )
    forcing_archive = case / f"shortwave_forcing.{case.name}.nc"
    if not forcing_archive.is_file():
        cautions.append(
            f"{forcing_archive.name} is absent; preserve or generate the original "
            "DNI and diffuse-sky archive before pedestrian shortwave processing"
        )
    if "receptor_height" not in f90nml.read(options_path).get("output", {}):
        cautions.append("receptor_height is implicit (UDBase/Fortran default 1.1 m); set it explicitly")
    return ComfortPreflight(
        case_dir=case,
        receptor_height_m=receptor_heights[0] if receptor_heights else float("nan"),
        receptor_heights_m=receptor_heights,
        first_scalar_height_m=float(heights[0]),
        analysis_times_s=targets,
        estimated_stats_bytes=estimated,
        blockers=tuple(blockers),
        cautions=tuple(cautions),
    )


_ATMOSPHERE = ("ta", "pabs", "qv", "ws_local", "ws_10")
_SHORTWAVE = ("sw_direct_normal",) + tuple(
    f"sw_nondirect_{side}face" for side in
    ("up", "down", "north", "south", "east", "west")
)
_LONGWAVE = tuple(
    f"lw_{side}face" for side in ("up", "down", "north", "south", "east", "west")
)
_UNITS = {"ta": "K", "pabs": "Pa", "qv": "kg kg-1",
          "ws_local": "m s-1", "ws_10": "m s-1"}
_UNITS.update({name: "W m-2" for name in _SHORTWAVE + _LONGWAVE})


@dataclass(frozen=True)
class ExchangeMetadata:
    """Verified georeferencing and run identity not recoverable from uDALES outputs.

    x/y and their bounds are projected *model-axis* coordinates. For a rotated
    LES grid, the supplied CRS must describe those axes, not merely EPSG:2154.
    Longitude/latitude must be transformed from the same grid definition.
    """

    x: np.ndarray
    y: np.ndarray
    x_bounds: np.ndarray
    y_bounds: np.ndarray
    longitude: np.ndarray
    latitude: np.ndarray
    z_ground: np.ndarray
    crs_attributes: dict[str, str | float]
    vertical_datum: str
    model_version: str
    institution: str
    building_representation: str
    terrain_convention: str
    spinup_start_utc: str
    spinup_end_utc: str


def _check_exchange_source(ds: Dataset, path: Path, x: np.ndarray, y: np.ndarray,
                           times: np.ndarray, height: float, fields: tuple[str, ...],
                           *, radiation: bool) -> np.ndarray:
    """Validate a previously generated native-grid plane file without copying it."""
    for name in ("x", "y", "time", "time_bounds", "valid_mask", "receptor_height", *fields):
        if name not in ds.variables:
            raise ValueError(f"{path.name} lacks {name}")
    for name, expected in (("x", x), ("y", y), ("time", times)):
        value = np.asarray(ds[name][:], dtype=float)
        if value.shape != expected.shape or not np.allclose(value, expected, rtol=0, atol=1e-5):
            raise ValueError(f"{path.name} has a different {name} grid")
    bounds = np.column_stack((times - 900.0, times))
    if (ds["time_bounds"].dimensions != ("time", "bounds")
            or not np.allclose(ds["time_bounds"][:], bounds, rtol=0, atol=1e-5)):
        raise ValueError(f"{path.name} has incorrect 15-minute time bounds")
    if not np.isclose(float(ds["receptor_height"][...]), height, rtol=0, atol=1e-5):
        raise ValueError(f"{path.name} has a different receptor height")
    mask = np.asarray(ds["valid_mask"][:])
    if (ds["valid_mask"].dimensions != ("x", "y")
            or mask.shape != (len(x), len(y)) or not np.isin(mask, (0, 1)).all()):
        raise ValueError(f"{path.name} has an invalid receptor mask")
    if not radiation:
        if "source_time" not in ds.variables or ds["source_time"].dimensions != ("time",):
            raise ValueError(f"{path.name} lacks actual statistics record times")
        actual = np.asarray(ds["source_time"][:], dtype=float)
        if actual.shape != times.shape or not np.isfinite(actual).all() or np.any(abs(actual - times) > 1.0):
            raise ValueError(f"{path.name} statistics records miss the requested hours")
    for name in fields:
        var = ds[name]
        if var.dimensions != ("x", "y", "time") or var.shape != (len(x), len(y), len(times)):
            raise ValueError(f"{path.name}:{name} must have dimensions (x,y,time)")
        if getattr(var, "units", None) != _UNITS[name]:
            raise ValueError(f"{path.name}:{name} has unexpected units")
    if radiation:
        if "time_complete" not in ds.variables or "processed_mask" not in ds.variables:
            raise ValueError(f"{path.name} lacks radiation completion masks")
        if (ds["time_complete"].dimensions != ("time",)
                or not np.array_equal(np.asarray(ds["time_complete"][:]), np.ones(len(times), dtype=int))):
            raise ValueError(f"{path.name} contains incomplete radiation windows")
        processed = np.asarray(ds["processed_mask"][:])
        if (ds["processed_mask"].dimensions != ("x", "y")
                or processed.shape != mask.shape or not np.isin(processed, (0, 1)).all()):
            raise ValueError(f"{path.name} has an invalid processed mask")
    return mask.astype(bool)


def _check_georeference(meta: ExchangeMetadata, sim, shape: tuple[int, int]) -> None:
    nx, ny = shape
    for name, expected_shape in (("x", (nx,)), ("y", (ny,)),
                                 ("x_bounds", (nx, 2)), ("y_bounds", (ny, 2)),
                                 ("longitude", shape), ("latitude", shape),
                                 ("z_ground", shape)):
        value = np.asarray(getattr(meta, name), dtype=float)
        if value.shape != expected_shape or not np.isfinite(value).all():
            raise ValueError(f"Georeference {name} must be finite with shape {expected_shape}")
    for axis, spacing in (("x", sim.dx), ("y", sim.dy)):
        centers = np.asarray(getattr(meta, axis), dtype=float)
        bounds = np.asarray(getattr(meta, f"{axis}_bounds"), dtype=float)
        if (not np.allclose(bounds[:, 1] - bounds[:, 0], spacing, rtol=0, atol=1e-5)
                or not np.allclose(bounds.mean(axis=1), centers, rtol=0, atol=1e-5)
                or not np.allclose(np.diff(centers), spacing, rtol=0, atol=1e-5)
                or not np.allclose(bounds[1:, 0], bounds[:-1, 1], rtol=0, atol=1e-5)):
            raise ValueError(f"Georeference {axis} coordinates/bounds disagree with LES spacing")
    if (np.any(np.abs(meta.longitude) > 180) or np.any(np.abs(meta.latitude) > 90)
            or not meta.crs_attributes.get("grid_mapping_name")
            or not (meta.crs_attributes.get("crs_wkt") or meta.crs_attributes.get("spatial_ref"))):
        raise ValueError("Georeference needs valid WGS84 coordinates and a complete CF CRS/WKT")
    if not np.allclose(meta.z_ground, np.asarray(meta.z_ground)[0, 0], rtol=0, atol=1e-5):
        raise ValueError("uDALES has flat model ground; z_ground must be a constant datum height")
    for name in ("vertical_datum", "model_version", "institution", "building_representation",
                 "terrain_convention", "spinup_start_utc", "spinup_end_utc"):
        if not isinstance(getattr(meta, name), str) or not getattr(meta, name).strip():
            raise ValueError(f"Exchange metadata {name} must be explicitly supplied")


def write_exchange_heights(sim, analysis_day: date, metadata: ExchangeMetadata, *,
                           heights: tuple[float, ...] | None = None,
                           input_dir: Path | None = None, output_dir: Path | None = None,
                           overwrite: bool = False) -> dict[float, Path]:
    """Assemble strict model-neutral 24-hour input files from saved planes.

    This is an I/O operation; atmospheric and radiative quantities are not
    recalculated. Solar angles alone are evaluated at exact window midpoints
    with the existing SPA routine because those instants are not archived.
    """
    case = Path(sim.path)
    source_dir = case if input_dir is None else Path(input_dir)
    target_dir = case if output_dir is None else Path(output_dir)
    selected = configured_receptor_heights(case) if heights is None else validate_heights(heights)
    configured = configured_receptor_heights(case)
    if any(not any(np.isclose(height, saved, rtol=0, atol=1e-5) for saved in configured)
           for height in selected):
        raise ValueError("Requested heights must be saved receptor heights from namoptions")
    options = {key.lower(): value for key, value in
               parse_namoptions(case / f"namoptions.{sim.expnr}").items()}
    if "timezone" not in options or float(options["timezone"]) != 0:
        raise ValueError("Exchange output requires a UTC simulation clock")
    start = datetime(int(options["year"]), int(options["month"]), int(options["day"]),
                     int(options["hour"]), int(options.get("minute", 0)),
                     int(options.get("second", 0)), tzinfo=timezone.utc)
    midnight = datetime.combine(analysis_day, datetime.min.time(), tzinfo=timezone.utc)
    expected = np.arange(1, 25, dtype=float) * 3600.0
    model_times = expected + (midnight - start).total_seconds()
    _check_georeference(metadata, sim, (len(sim.xt), len(sim.yt)))
    if (datetime.fromisoformat(metadata.spinup_start_utc.replace("Z", "+00:00")) != start
            or datetime.fromisoformat(metadata.spinup_end_utc.replace("Z", "+00:00")) != midnight):
        raise ValueError("UTC spin-up dates must match simulation start and analysis-day midnight")
    for key in ("longitude", "latitude", "elevation"):
        if key not in options or not np.isfinite(float(options[key])):
            raise ValueError(f"Solar {key} must be explicitly configured")
    solar = [solar_position_python(midnight + timedelta(seconds=float(t - 450)),
                                   float(options["longitude"]), float(options["latitude"]),
                                   0.0, float(options["elevation"]))
             for t in expected]
    zenith = np.deg2rad([float(state["zenith"]) for state in solar])
    azimuth = np.mod(np.deg2rad([float(state["azimuth"]) for state in solar]), 2 * np.pi)
    if not np.isfinite(zenith).all() or not np.isfinite(azimuth).all():
        raise ValueError("Solar geometry is not finite")
    paths = {h: target_dir / f"thermal_comfort_inputs_udales.{height_tag(h)}.{sim.expnr}.nc"
             for h in selected}
    if not overwrite and any(path.exists() for path in paths.values()):
        raise FileExistsError("One or more exchange files already exist")
    target_dir.mkdir(parents=True, exist_ok=True)
    for height, path in paths.items():
        sources = {
            kind: source_dir / f"pedestrian_{kind}.{height_tag(height)}.{sim.expnr}.nc"
            for kind in ("atmosphere", "shortwave", "longwave")
        }
        with ExitStack() as stack:
            datasets = {}
            for kind, source in sources.items():
                datasets[kind] = stack.enter_context(Dataset(source, "r"))
                datasets[kind].set_auto_mask(False)
            x = np.asarray(sim.xt, dtype=float)
            y = np.asarray(sim.yt, dtype=float)
            atm = _check_exchange_source(datasets["atmosphere"], sources["atmosphere"],
                                         x, y, model_times, height, _ATMOSPHERE, radiation=False)
            if not atm.any():
                raise ValueError("No valid pedestrian cells")
            for kind, fields in (("shortwave", _SHORTWAVE), ("longwave", _LONGWAVE)):
                ds = datasets[kind]
                mask = _check_exchange_source(ds, sources[kind], x, y, model_times,
                                              height, fields, radiation=True)
                if np.any(atm & ~mask) or np.any(atm & (np.asarray(ds["processed_mask"][:]) != 1)):
                    raise ValueError(f"{kind} does not cover every atmospheric-valid receptor")
            fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=target_dir)
            os.close(fd)
            temporary = Path(name)
            try:
                with Dataset(temporary, "w", format="NETCDF4") as out:
                    nx, ny = len(x), len(y)
                    for dim, length in (("x", nx), ("y", ny), ("time", 24), ("bounds", 2)):
                        out.createDimension(dim, length)
                    for axis in ("x", "y"):
                        out.createVariable(axis, "f8", (axis,))[:] = getattr(metadata, axis)
                        out[axis].units = "m"
                        out[axis].standard_name = f"projection_{axis}_coordinate"
                        out[axis].axis = axis.upper()
                        out.createVariable(f"{axis}_bounds", "f8", (axis, "bounds"))[:] = getattr(metadata, f"{axis}_bounds")
                        out[f"{axis}_bounds"].units = "m"
                        out[axis].bounds = f"{axis}_bounds"
                    for coord, units in (("longitude", "degrees_east"), ("latitude", "degrees_north"),
                                         ("z_ground", "m")):
                        out.createVariable(coord, "f8", ("x", "y"))[:] = getattr(metadata, coord)
                        out[coord].units = units
                    out["longitude"].standard_name = "longitude"
                    out["latitude"].standard_name = "latitude"
                    out["z_ground"].standard_name = "surface_altitude"
                    out["z_ground"].vertical_datum = metadata.vertical_datum
                    out.createVariable("receptor_height", "f8").assignValue(height)
                    out["receptor_height"].units = "m above local terrain"
                    out.createVariable("pedestrian_mask", "i1", ("x", "y"))[:] = atm.astype("i1")
                    out["pedestrian_mask"].long_name = "valid pedestrian receptor"
                    out["pedestrian_mask"].flag_values = np.array([0, 1], dtype="i1")
                    out["pedestrian_mask"].flag_meanings = "invalid valid"
                    crs = out.createVariable("crs", "i4")
                    crs.assignValue(0)
                    crs.setncatts(metadata.crs_attributes)
                    time = out.createVariable("time", "f8", ("time",))
                    time[:] = expected
                    time.units = f"seconds since {analysis_day.isoformat()} 00:00:00 UTC"
                    time.calendar = "gregorian"
                    time.bounds = "time_bounds"
                    time.standard_name = "time"
                    time.axis = "T"
                    out.createVariable("time_bounds", "f8", ("time", "bounds"))[:] = np.column_stack((expected - 900, expected))
                    out["time_bounds"].units = time.units
                    out.createVariable("solar_zenith", "f8", ("time",))[:] = zenith
                    out.createVariable("solar_azimuth", "f8", ("time",))[:] = azimuth
                    out["solar_zenith"].units = out["solar_azimuth"].units = "rad"
                    out["solar_azimuth"].long_name = "solar azimuth clockwise from true north"
                    out["solar_zenith"].long_name = "solar zenith at midpoint of averaging window"
                    out.model_name = "uDALES"
                    out.Conventions = "CF-1.10"
                    out.model_version = metadata.model_version
                    out.institution = metadata.institution
                    out.simulation_identifier = str(sim.expnr)
                    out.horizontal_grid_spacing_m = f"x={sim.dx:g}, y={sim.dy:g}"
                    out.vertical_grid_spacing_m = ",".join(f"{dz:.9g}" for dz in np.asarray(sim.dzt))
                    out.terrain_convention = metadata.terrain_convention
                    out.building_representation = metadata.building_representation
                    out.model_axis_rotation_degrees_from_true_north = float(sim.xazimuth)
                    out.vertical_interpolation_method = "linear in physical height between scalar centres"
                    out.spinup_start_utc = metadata.spinup_start_utc
                    out.spinup_end_utc = metadata.spinup_end_utc
                    out.solar_geometry = "udprep.solar SPA, exact 15-minute window midpoint"
                    out.source_files = ", ".join(path.name for path in sources.values())
                    for kind, names in (("atmosphere", _ATMOSPHERE),
                                        ("shortwave", _SHORTWAVE), ("longwave", _LONGWAVE)):
                        for field in names:
                            variable = out.createVariable(
                                field, "f4", ("x", "y", "time"),
                                fill_value=np.float32(-9.96921e36), zlib=True, complevel=3,
                                chunksizes=(min(nx, 64), min(ny, 64), 1),
                            )
                            variable.units = _UNITS[field]
                            variable.cell_methods = "time: mean"
                            variable.grid_mapping = "crs"
                            variable.coordinates = "longitude latitude"
                    for index in range(24):
                        for kind, names in (("atmosphere", _ATMOSPHERE),
                                            ("shortwave", _SHORTWAVE), ("longwave", _LONGWAVE)):
                            for field in names:
                                values = np.asarray(datasets[kind][field][:, :, index], dtype=float)
                                if not np.isfinite(values[atm]).all():
                                    raise ValueError(f"{sources[kind].name}:{field} has missing values in valid cells at hour {index + 1}")
                                if field in ("ta", "pabs") and np.any(values[atm] <= 0):
                                    raise ValueError(f"{field} must be positive in valid cells")
                                if field not in ("ta", "pabs") and np.any(values[atm] < 0):
                                    raise ValueError(f"{field} must be nonnegative in valid cells")
                                out[field][:, :, index] = np.where(atm, values, out[field]._FillValue)
                if path.exists() and not overwrite:
                    raise FileExistsError(f"Exchange output appeared during export: {path}")
                os.replace(temporary, path)
            finally:
                temporary.unlink(missing_ok=True)
    return paths
