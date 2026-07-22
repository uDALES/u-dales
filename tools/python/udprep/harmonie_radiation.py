"""HARMONIE radiation helpers for uDALES preprocessing.

This module bridges HARMONIE accumulated surface radiation fields to uDALES
time-dependent radiation inputs. The shortwave workflow uses HARMONIE ``ssrd``
for the domain-mean global horizontal irradiance, then maps that atmospheric
forcing onto uDALES facets with the existing direct-shortwave, view-factor, and
shortwave-reflection machinery.
"""

from __future__ import annotations

import math
import runpy
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from . import _radiation_compute
from .solar import nsun_from_angles, solar_position_python


DEFAULT_NWP_ROOT = Path("/data/simulateParis/NWP_demo_data")
DEFAULT_VERSION = "8.0"
DEFAULT_URL = "http://exporter.nsc.liu.se/aebc1d3690d441cf82818d9893fa9e57"
DEFAULT_CYCLE_TIME = "2023-08-20T00:00:00Z"
DEFAULT_OUTPUT_FREQUENCY = "PT15M"
DEFAULT_GRID_SHAPE = (989, 989)  # ny, nx
DEFAULT_BBOX_LAMBERT93 = (649375.0, 652600.0, 6861175.0, 6864020.0)

SOLAR_CONSTANT_W_M2 = 1367.0
MIN_DIRECT_COS_ZENITH = 1.0e-2
ACCUMULATED_FLUX_TOLERANCE_W_M2 = 1.0e-2


@dataclass(frozen=True)
class VersionInfo:
    """Metadata needed to locate and decode a HARMONIE data version."""

    url: str
    cycle_time: datetime
    output_frequency_seconds: int
    grid_shape: tuple[int, int]


@dataclass(frozen=True)
class ShortwaveAtmosphere:
    """Atmospheric shortwave forcing prepared for uDALES facet mapping."""

    times: np.ndarray
    ghi: np.ndarray
    dni: np.ndarray
    dsky: np.ndarray
    zenith: np.ndarray
    azimuth_local: np.ndarray


@dataclass(frozen=True)
class NamelistRadiationConfig:
    """Minimal radiation-like object built without loading UDPrep geometry."""

    runtime: float
    dtSP: float
    ntimedepsw: int
    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: int
    longitude: float
    latitude: float
    timezone: float
    elevation: float
    xazimuth: float


@dataclass(frozen=True)
class NamelistLongwaveConfig:
    """Minimal namelist state needed for HARMONIE STRD longwave forcing."""

    expnr: str
    runtime: float
    ntimedeplw: int
    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: int


@dataclass(frozen=True)
class TimedepShortwaveResult:
    """Generated uDALES shortwave fields."""

    times: np.ndarray
    atmosphere: ShortwaveAtmosphere
    sdir: np.ndarray
    knet: np.ndarray
    sveg: np.ndarray | None
    timedepsw_path: Path
    sdir_nc_path: Path | None
    timedepsveg_path: Path | None


@dataclass(frozen=True)
class LongwaveSeries:
    """Scalar sky longwave forcing prepared for uDALES."""

    times: np.ndarray
    lwsky: np.ndarray
    mask_points: int


@dataclass(frozen=True)
class TimedepLongwaveResult:
    """Generated uDALES longwave input."""

    times: np.ndarray
    lwsky: np.ndarray
    timedeplw_path: Path
    mask_points: int


def import_grib_dependencies() -> tuple[Any, Any, Any]:
    """Import GRIB dependencies lazily so pure helpers remain lightweight."""
    try:
        import cfgrib  # noqa: F401
        import xarray as xr
        from pyproj import Transformer
    except ImportError as exc:
        raise RuntimeError(
            "Reading HARMONIE GRIB files requires xarray, cfgrib, and pyproj. "
            "Use the NWP demo environment, for example "
            "/data/simulateParis/NWP_demo_data/.venv/bin/python."
        ) from exc
    return np, xr, Transformer


def parse_iso_time(text: str) -> datetime:
    clean = text.strip()
    if clean.endswith("Z"):
        clean = clean[:-1] + "+00:00"
    value = datetime.fromisoformat(clean)
    if value.tzinfo is None:
        value = value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def parse_iso_duration_seconds(text: str | None, default: int) -> int:
    if not text:
        return default
    import re

    match = re.fullmatch(
        r"P(?:\d+D)?T(?:(\d+)H)?(?:(\d+)M)?(?:(\d+)S)?", text.strip()
    )
    if not match:
        return default
    hours = int(match.group(1) or 0)
    minutes = int(match.group(2) or 0)
    seconds = int(match.group(3) or 0)
    total = hours * 3600 + minutes * 60 + seconds
    return total or default


def find_namoptions(case_dir: Path) -> Path:
    candidate = Path(case_dir) / f"namoptions.{Path(case_dir).name}"
    if candidate.exists():
        return candidate

    matches = sorted(Path(case_dir).glob("namoptions.*"))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(f"No namoptions.* file found in {case_dir}")
    raise ValueError(
        f"Multiple namoptions.* files found in {case_dir}; use a canonical case directory."
    )


def _required(values: dict[str, Any], key: str) -> Any:
    try:
        return values[key]
    except KeyError as exc:
        raise KeyError(f"Required namoptions key {key!r} is missing") from exc


def radiation_config_from_case_dir(case_dir: Path) -> NamelistRadiationConfig:
    """Build the minimal radiation config needed for atmospheric conversion."""
    from udconfig import parse_namoptions

    values = parse_namoptions(find_namoptions(Path(case_dir)))
    return NamelistRadiationConfig(
        runtime=float(_required(values, "runtime")),
        dtSP=float(_required(values, "dtSP")),
        ntimedepsw=int(_required(values, "ntimedepsw")),
        year=int(_required(values, "year")),
        month=int(_required(values, "month")),
        day=int(_required(values, "day")),
        hour=int(_required(values, "hour")),
        minute=int(_required(values, "minute")),
        second=int(_required(values, "second")),
        longitude=float(_required(values, "longitude")),
        latitude=float(_required(values, "latitude")),
        timezone=float(_required(values, "timezone")),
        elevation=float(_required(values, "elevation")),
        xazimuth=float(_required(values, "xazimuth")),
    )


def longwave_config_from_case_dir(case_dir: Path) -> NamelistLongwaveConfig:
    """Build the minimal longwave config from a uDALES case directory."""
    from udconfig import parse_namoptions

    case_dir = Path(case_dir)
    values = parse_namoptions(find_namoptions(case_dir))
    return NamelistLongwaveConfig(
        expnr=case_dir.name,
        runtime=float(_required(values, "runtime")),
        ntimedeplw=int(_required(values, "ntimedeplw")),
        year=int(_required(values, "year")),
        month=int(_required(values, "month")),
        day=int(_required(values, "day")),
        hour=int(_required(values, "hour")),
        minute=int(_required(values, "minute")),
        second=int(_required(values, "second")),
    )


def load_version_info(nwp_root: Path, version: str = DEFAULT_VERSION) -> VersionInfo:
    """Read HARMONIE version metadata from the NWP demo ``versions.py`` file."""
    url = DEFAULT_URL
    metadata: dict[str, Any] = {
        "date": DEFAULT_CYCLE_TIME,
        "output_frequency": DEFAULT_OUTPUT_FREQUENCY,
        "nx": DEFAULT_GRID_SHAPE[1],
        "ny": DEFAULT_GRID_SHAPE[0],
    }

    versions_py = nwp_root / "versions.py"
    if versions_py.exists():
        namespace = runpy.run_path(str(versions_py))
        urban_air_data = namespace.get("UrbanAirData")
        urls = getattr(urban_air_data, "urls", {})
        try:
            version_entry = urls[str(version)]
        except KeyError as exc:
            raise KeyError(
                f"Version {version!r} is not listed in {versions_py}"
            ) from exc
        url = version_entry["url"]
        metadata = version_entry.get("metadata", metadata)

    nx = int(metadata.get("nx", DEFAULT_GRID_SHAPE[1]))
    ny = int(metadata.get("ny", DEFAULT_GRID_SHAPE[0]))
    return VersionInfo(
        url=url,
        cycle_time=parse_iso_time(str(metadata.get("date", DEFAULT_CYCLE_TIME))),
        output_frequency_seconds=parse_iso_duration_seconds(
            metadata.get("output_frequency"), default=900
        ),
        grid_shape=(ny, nx),
    )


def data_dir_from_url(nwp_root: Path, url: str) -> Path:
    """Return the local path convention used by NWP_demo_data/download.py."""
    return nwp_root / "data" / url


def day_directory(data_dir: Path, cycle_time: datetime) -> Path:
    if any(data_dir.glob("GRIBPFDEOD+*")):
        return data_dir
    return data_dir / f"{cycle_time.year:04d}" / f"{cycle_time.month:02d}" / (
        f"{cycle_time.day:02d}"
    )


def ensure_data_available(
    day_dir: Path, nwp_root: Path, version: str, *, download: bool
) -> None:
    if day_dir.exists():
        return
    if not download:
        raise FileNotFoundError(
            f"NWP day directory not found: {day_dir}\n"
            "Pass --data-dir if the files live elsewhere, or pass --download "
            "to run the NWP demo downloader."
        )

    download_py = nwp_root / "download.py"
    if not download_py.exists():
        raise FileNotFoundError(f"Cannot find downloader: {download_py}")
    subprocess.run(
        [sys.executable, str(download_py), "-v", str(version)],
        cwd=nwp_root,
        check=True,
    )


def format_forecast_offset(seconds: int) -> str:
    if seconds < 0:
        raise ValueError(f"Negative forecast offset: {seconds} seconds")
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:04d}h{minutes:02d}m{seconds:02d}s"


def make_model_times(runtime: float, dt: float, nt: int | None = None) -> np.ndarray:
    """Return uDALES radiation times, matching the existing dtSP convention."""
    if nt is not None:
        return np.arange(int(nt), dtype=np.float64) * float(dt)
    return np.arange(0.0, float(runtime) + 0.5 * float(dt), float(dt), dtype=float)


def nearest_common_interval(seconds: float) -> int:
    common = (300, 600, 900, 1800, 3600, 7200, 10800)
    nearest = min(common, key=lambda item: abs(item - seconds))
    if abs(nearest - seconds) <= max(2.0, 0.02 * nearest):
        return nearest
    return int(round(seconds))


def make_longwave_times(
    runtime: float, nt: int, output_interval_seconds: int | None = None
) -> np.ndarray:
    """Return timedeplw timestamps, preserving the hourly case-300 convention."""
    if nt <= 0:
        raise ValueError("ntimedeplw must be positive")
    if output_interval_seconds is None:
        if nt == 1:
            output_interval_seconds = 0
        else:
            output_interval_seconds = nearest_common_interval(runtime / float(nt - 1))
    return np.arange(int(nt), dtype=np.float64) * float(output_interval_seconds)


class HarmonieSurfaceFieldReader:
    """Read and spatially average accumulated HARMONIE surface fields."""

    def __init__(
        self,
        *,
        day_dir: Path,
        field: str,
        type_of_level: str = "surface",
        grid_shape: tuple[int, int] = DEFAULT_GRID_SHAPE,
        bbox_lambert93: tuple[float, float, float, float] = DEFAULT_BBOX_LAMBERT93,
        grib_prefix: str = "GRIBPFDEOD+",
        verbose: bool = True,
    ) -> None:
        _, self.xr, transformer_cls = import_grib_dependencies()
        self.day_dir = Path(day_dir)
        self.field = field
        self.type_of_level = type_of_level
        self.grid_shape = grid_shape
        self.bbox_lambert93 = bbox_lambert93
        self.grib_prefix = grib_prefix
        self.verbose = verbose
        self.transformer = transformer_cls.from_crs(
            "EPSG:4326", "EPSG:2154", always_xy=True
        )
        self.mask: Any | None = None
        self.mask_points = 0
        self.cache: dict[int, float] = {}

    def path_for_offset(self, offset_seconds: int) -> Path:
        return self.day_dir / f"{self.grib_prefix}{format_forecast_offset(offset_seconds)}"

    def reshape_grid(self, values: Any, name: str) -> np.ndarray:
        arr = np.asarray(values).squeeze()
        if arr.ndim == 2:
            return arr
        expected = self.grid_shape[0] * self.grid_shape[1]
        if arr.size != expected:
            raise ValueError(
                f"Cannot reshape {name}: got {arr.size} values, expected {expected} "
                f"for grid shape {self.grid_shape}."
            )
        return arr.reshape(self.grid_shape)

    def ensure_mask(self, dataset: Any) -> None:
        if self.mask is not None:
            return

        lon = self.reshape_grid(dataset.longitude.values, "longitude")
        lat = self.reshape_grid(dataset.latitude.values, "latitude")
        east, north = self.transformer.transform(lon, lat)
        xmin, xmax, ymin, ymax = self.bbox_lambert93
        self.mask = (
            (east >= xmin) & (east <= xmax) & (north >= ymin) & (north <= ymax)
        )
        self.mask_points = int(self.mask.sum())
        if self.mask_points == 0:
            raise ValueError(
                "The Lambert-93 bounding box did not select any HARMONIE grid cells."
            )

    def mean_accumulation(self, offset_seconds: int) -> float:
        if offset_seconds in self.cache:
            return self.cache[offset_seconds]

        path = self.path_for_offset(offset_seconds)
        if not path.exists():
            raise FileNotFoundError(f"Missing GRIB file for offset {path}")
        if self.verbose:
            print(f"Reading {path.name}", file=sys.stderr, flush=True)

        dataset = self.xr.open_dataset(
            path,
            engine="cfgrib",
            backend_kwargs={
                "filter_by_keys": {
                    "shortName": self.field,
                    "typeOfLevel": self.type_of_level,
                },
                "indexpath": "",
            },
        )
        try:
            var_name = self.field if self.field in dataset.data_vars else None
            if var_name is None:
                data_vars = list(dataset.data_vars)
                if len(data_vars) != 1:
                    raise ValueError(
                        f"Could not identify {self.field!r} in {path}; "
                        f"data variables are {data_vars}."
                    )
                var_name = data_vars[0]
            self.ensure_mask(dataset)
            field_values = self.reshape_grid(dataset[var_name].values, var_name)
            mean_value = float(
                np.nanmean(field_values[self.mask], dtype=np.float64)
            )
        finally:
            dataset.close()

        self.cache[offset_seconds] = mean_value
        return mean_value


def accumulated_flux_series(
    reader: HarmonieSurfaceFieldReader,
    *,
    start_offset_seconds: int,
    end_offset_seconds: int,
    difference_interval_seconds: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert accumulated J/m2 values to a native-time flux series in W/m2."""
    if difference_interval_seconds <= 0:
        raise ValueError("difference_interval_seconds must be positive")
    if end_offset_seconds < start_offset_seconds:
        raise ValueError("end_offset_seconds must be >= start_offset_seconds")

    offsets = np.arange(
        int(start_offset_seconds),
        int(end_offset_seconds) + difference_interval_seconds,
        int(difference_interval_seconds),
        dtype=int,
    )
    model_times = offsets.astype(np.float64) - float(start_offset_seconds)
    fluxes = np.empty(offsets.size, dtype=np.float64)

    for idx, offset in enumerate(offsets):
        previous = int(offset) - int(difference_interval_seconds)
        if previous < 0:
            raise ValueError(
                "Cannot difference accumulated radiation before forecast start. "
                f"Needed {previous} seconds for offset {offset}."
            )
        current_accum = reader.mean_accumulation(int(offset))
        previous_accum = reader.mean_accumulation(previous)
        flux = (current_accum - previous_accum) / float(difference_interval_seconds)
        if flux < -ACCUMULATED_FLUX_TOLERANCE_W_M2:
            raise ValueError(
                f"Negative {reader.field} flux at forecast offset "
                f"{offset} seconds: {flux:.6f} W/m2."
            )
        fluxes[idx] = max(0.0, flux)

    return model_times, fluxes


def interpolate_flux_to_times(
    source_times: np.ndarray, source_fluxes: np.ndarray, target_times: np.ndarray
) -> np.ndarray:
    if source_times.size == 0:
        raise ValueError("source_times is empty")
    return np.interp(target_times, source_times, source_fluxes)


def extraterrestrial_horizontal_irradiance(
    when: datetime, cos_zenith: float
) -> float:
    """Extraterrestrial horizontal irradiance used for clearness-index splits."""
    if cos_zenith <= 0.0:
        return 0.0
    day_of_year = int(when.timetuple().tm_yday)
    eccentricity = 1.0 + 0.033 * math.cos(2.0 * math.pi * day_of_year / 365.0)
    return SOLAR_CONSTANT_W_M2 * eccentricity * cos_zenith


def erbs_diffuse_fraction(clearness_index: float) -> float:
    """Diffuse fraction for global horizontal irradiance using Erbs et al."""
    kt = max(0.0, float(clearness_index))
    if kt <= 0.22:
        return 1.0 - 0.09 * kt
    if kt <= 0.80:
        return (
            0.9511
            - 0.1604 * kt
            + 4.388 * kt**2
            - 16.638 * kt**3
            + 12.336 * kt**4
        )
    return 0.165


def split_global_horizontal_erbs(
    ghi: float,
    zenith_deg: float,
    when: datetime,
    *,
    min_direct_cos_zenith: float = MIN_DIRECT_COS_ZENITH,
) -> tuple[float, float]:
    """Split GHI into direct normal irradiance and diffuse horizontal irradiance.

    Returns ``(dni, dsky)`` in W/m2. For very low sun angles, all positive GHI is
    treated as diffuse because the direct-shortwave ray tracer intentionally
    skips near-horizontal sun vectors.
    """
    ghi = max(0.0, float(ghi))
    if ghi <= 0.0:
        return 0.0, 0.0

    cos_zenith = math.cos(math.radians(float(zenith_deg)))
    if zenith_deg >= 90.0 or cos_zenith <= 0.0:
        return 0.0, 0.0
    if cos_zenith < min_direct_cos_zenith:
        return 0.0, ghi

    i0h = extraterrestrial_horizontal_irradiance(when, cos_zenith)
    if i0h <= 0.0:
        return 0.0, ghi

    diffuse_fraction = np.clip(erbs_diffuse_fraction(ghi / i0h), 0.0, 1.0)
    dsky = float(np.clip(diffuse_fraction * ghi, 0.0, ghi))
    dni = max(0.0, (ghi - dsky) / cos_zenith)
    return dni, dsky


def prepare_shortwave_atmosphere(
    *,
    radiation: Any,
    times: np.ndarray,
    ghi: np.ndarray,
) -> ShortwaveAtmosphere:
    """Compute solar geometry and DNI/Dsky for target uDALES times."""
    start = datetime(
        int(radiation.year),
        int(radiation.month),
        int(radiation.day),
        int(radiation.hour),
        int(radiation.minute),
        int(radiation.second),
    )
    dni = np.zeros_like(ghi, dtype=np.float64)
    dsky = np.zeros_like(ghi, dtype=np.float64)
    zenith = np.zeros_like(ghi, dtype=np.float64)
    azimuth_local = np.zeros_like(ghi, dtype=np.float64)

    for idx, time_value in enumerate(times):
        when = start + timedelta(seconds=float(time_value))
        sp = solar_position_python(
            when,
            float(radiation.longitude),
            float(radiation.latitude),
            float(radiation.timezone),
            float(radiation.elevation),
        )
        zen = float(sp["zenith"])
        az_local = float(sp["azimuth"]) - float(radiation.xazimuth)
        this_dni, this_dsky = split_global_horizontal_erbs(float(ghi[idx]), zen, when)
        zenith[idx] = zen
        azimuth_local[idx] = az_local
        dni[idx] = this_dni
        dsky[idx] = this_dsky

    return ShortwaveAtmosphere(
        times=times,
        ghi=ghi,
        dni=dni,
        dsky=dsky,
        zenith=zenith,
        azimuth_local=azimuth_local,
    )


def prepare_harmonie_ssrd_atmosphere(
    *,
    radiation: Any,
    nwp_root: Path = DEFAULT_NWP_ROOT,
    version: str = DEFAULT_VERSION,
    data_dir: Path | None = None,
    download: bool = False,
    cycle_time: datetime | None = None,
    field: str = "ssrd",
    type_of_level: str = "surface",
    bbox_lambert93: tuple[float, float, float, float] = DEFAULT_BBOX_LAMBERT93,
    grib_prefix: str = "GRIBPFDEOD+",
    verbose: bool = True,
) -> tuple[ShortwaveAtmosphere, HarmonieSurfaceFieldReader]:
    """Read HARMONIE ``ssrd`` and prepare a uDALES shortwave atmosphere series."""
    nwp_root = Path(nwp_root).expanduser().resolve()
    version_info = load_version_info(nwp_root, version)
    cycle = cycle_time or version_info.cycle_time
    data_root = (
        Path(data_dir).expanduser().resolve()
        if data_dir is not None
        else data_dir_from_url(nwp_root, version_info.url)
    )
    day_dir = day_directory(data_root, cycle)
    ensure_data_available(day_dir, nwp_root, str(version), download=download)

    start = datetime(
        int(radiation.year),
        int(radiation.month),
        int(radiation.day),
        int(radiation.hour),
        int(radiation.minute),
        int(radiation.second),
        tzinfo=timezone.utc,
    )
    start_offset_seconds = int(round((start - cycle).total_seconds()))
    if start_offset_seconds < 0:
        raise ValueError(f"Case starts before HARMONIE cycle: {start} < {cycle}")

    ntimedep = getattr(radiation, "ntimedepsw", None)
    times = make_model_times(float(radiation.runtime), float(radiation.dtSP), ntimedep)
    end_offset_seconds = start_offset_seconds + int(round(float(times[-1])))

    reader = HarmonieSurfaceFieldReader(
        day_dir=day_dir,
        field=field,
        type_of_level=type_of_level,
        grid_shape=version_info.grid_shape,
        bbox_lambert93=bbox_lambert93,
        grib_prefix=grib_prefix,
        verbose=verbose,
    )
    native_times, native_ghi = accumulated_flux_series(
        reader,
        start_offset_seconds=start_offset_seconds,
        end_offset_seconds=end_offset_seconds,
        difference_interval_seconds=version_info.output_frequency_seconds,
    )
    ghi = interpolate_flux_to_times(native_times, native_ghi, times)
    atmosphere = prepare_shortwave_atmosphere(
        radiation=radiation,
        times=times,
        ghi=ghi,
    )
    return atmosphere, reader


def prepare_harmonie_strd_longwave(
    *,
    config: NamelistLongwaveConfig,
    nwp_root: Path = DEFAULT_NWP_ROOT,
    version: str = DEFAULT_VERSION,
    data_dir: Path | None = None,
    download: bool = False,
    cycle_time: datetime | None = None,
    output_interval_seconds: int | None = None,
    difference_interval_seconds: int | None = None,
    field: str = "strd",
    type_of_level: str = "surface",
    bbox_lambert93: tuple[float, float, float, float] = DEFAULT_BBOX_LAMBERT93,
    grib_prefix: str = "GRIBPFDEOD+",
    verbose: bool = True,
) -> LongwaveSeries:
    """Read HARMONIE ``strd`` and prepare a scalar uDALES LWsky series."""
    nwp_root = Path(nwp_root).expanduser().resolve()
    version_info = load_version_info(nwp_root, version)
    cycle = cycle_time or version_info.cycle_time
    data_root = (
        Path(data_dir).expanduser().resolve()
        if data_dir is not None
        else data_dir_from_url(nwp_root, version_info.url)
    )
    day_dir = day_directory(data_root, cycle)
    ensure_data_available(day_dir, nwp_root, str(version), download=download)

    start = datetime(
        int(config.year),
        int(config.month),
        int(config.day),
        int(config.hour),
        int(config.minute),
        int(config.second),
        tzinfo=timezone.utc,
    )
    start_offset_seconds = int(round((start - cycle).total_seconds()))
    if start_offset_seconds < 0:
        raise ValueError(f"Case starts before HARMONIE cycle: {start} < {cycle}")

    times = make_longwave_times(
        float(config.runtime),
        int(config.ntimedeplw),
        output_interval_seconds=output_interval_seconds,
    )
    end_offset_seconds = start_offset_seconds + int(round(float(times[-1])))
    diff_seconds = (
        int(difference_interval_seconds)
        if difference_interval_seconds is not None
        else int(version_info.output_frequency_seconds)
    )

    reader = HarmonieSurfaceFieldReader(
        day_dir=day_dir,
        field=field,
        type_of_level=type_of_level,
        grid_shape=version_info.grid_shape,
        bbox_lambert93=bbox_lambert93,
        grib_prefix=grib_prefix,
        verbose=verbose,
    )
    native_times, native_lwsky = accumulated_flux_series(
        reader,
        start_offset_seconds=start_offset_seconds,
        end_offset_seconds=end_offset_seconds,
        difference_interval_seconds=diff_seconds,
    )
    lwsky = interpolate_flux_to_times(native_times, native_lwsky, times)
    return LongwaveSeries(times=times, lwsky=lwsky, mask_points=reader.mask_points)


def write_timedepsw(path: Path, times: np.ndarray, knet: np.ndarray, *, overwrite: bool) -> None:
    """Write the uDALES timedepsw format."""
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Pass --overwrite to replace it.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(
            "# time-dependent net shortwave on facets [W/m2]. "
            "First line: times (1 x nt), then netsw (nfcts x nt)\n"
        )
        np.savetxt(handle, times[None, :], fmt="%9.2f")
        np.savetxt(handle, knet, fmt="%9.4f")


def write_timedepsveg(
    path: Path, times: np.ndarray, sveg: np.ndarray, *, overwrite: bool
) -> None:
    """Write optional time-dependent vegetation shortwave absorption."""
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Pass --overwrite to replace it.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(
            "# time-dependent vegetation shortwave absorption [W/m3]. "
            "First line: times (1 x nt), then sveg (nveg x nt)\n"
        )
        np.savetxt(handle, times[None, :], fmt="%9.2f")
        np.savetxt(handle, sveg, fmt="%9.4f")


def write_timedeplw(path: Path, times: np.ndarray, lwsky: np.ndarray, *, overwrite: bool) -> None:
    """Write the uDALES timedeplw format."""
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Pass --overwrite to replace it.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("Time-varying sky longwave flux from HARMONIE strd\n")
        handle.write("time             LWsky\n")
        for time_value, flux_value in zip(times, lwsky):
            handle.write(f"{time_value:13.6f} {flux_value:12.6f}\n")


def _default_case_output_paths(sim: Any) -> tuple[Path, Path, Path]:
    out_dir = Path(sim.path)
    expnr = str(sim.expnr)
    return (
        out_dir / f"timedepsw.inp.{expnr}",
        out_dir / "Sdir.nc",
        out_dir / f"timedepsveg.inp.{expnr}",
    )


def _check_output_paths(paths: Iterable[Path], *, overwrite: bool) -> None:
    if overwrite:
        return
    existing = [str(path) for path in paths if Path(path).exists()]
    if existing:
        raise FileExistsError(
            "Output file already exists. Pass --overwrite to replace: "
            + ", ".join(existing)
        )


def map_atmosphere_to_facets(
    prep: Any,
    atmosphere: ShortwaveAtmosphere,
    *,
    method: str | None = None,
    verbose: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Map atmospheric shortwave forcing to uDALES facets."""
    sim = prep.sim
    radiation = prep.radiation
    lscatter = bool(radiation.lEB)
    albedo = sim.assign_prop_to_fac("al")
    face_normals = sim.geom.stl.face_normals
    fss = (1.0 + face_normals[:, 2]) * 0.5 if not lscatter else None

    default_method, resolution = radiation._shortwave_method()
    method_name = method or default_method

    vf = None
    svf = None
    if lscatter:
        vf, svf, _ = radiation.calc_view_factors(maxD=float(radiation.maxD))

    nt = atmosphere.times.size
    nfacets = albedo.size
    sdir_all = np.zeros((nfacets, nt), dtype=np.float32)
    knet_all = np.zeros((nfacets, nt), dtype=np.float32)
    sveg_all: np.ndarray | None = None

    for idx, time_value in enumerate(atmosphere.times):
        dni = float(atmosphere.dni[idx])
        dsky = float(atmosphere.dsky[idx])
        zenith = float(atmosphere.zenith[idx])
        if verbose:
            print(
                f"[{idx + 1:3d}/{nt}] t={time_value:8.1f}s "
                f"GHI={atmosphere.ghi[idx]:8.3f} DNI={dni:8.3f} "
                f"Dsky={dsky:8.3f} zenith={zenith:7.3f}",
                flush=True,
            )

        if dni > 0.0 and zenith < 90.0:
            nsun = nsun_from_angles(zenith, float(atmosphere.azimuth_local[idx]))
            if abs(float(nsun[2])) >= MIN_DIRECT_COS_ZENITH:
                sdir, knet, sveg = radiation._compute_knet(
                    nsun,
                    dni,
                    dsky,
                    method_name,
                    resolution,
                    lscatter,
                    albedo,
                    vf,
                    svf,
                    fss,
                )
                sdir_all[:, idx] = np.asarray(sdir, dtype=np.float32)
                knet_all[:, idx] = np.asarray(knet, dtype=np.float32)
                if sveg is not None and np.asarray(sveg).size > 0:
                    if sveg_all is None:
                        sveg_all = np.zeros((np.asarray(sveg).size, nt), dtype=np.float32)
                    sveg_all[:, idx] = np.asarray(sveg, dtype=np.float32)
                continue

        if dsky > 0.0:
            zero_sdir = np.zeros(nfacets, dtype=np.float64)
            if lscatter:
                if vf is None or svf is None:
                    raise ValueError("View factors are required for diffuse shortwave")
                knet = radiation.calc_reflections_sw(zero_sdir, dsky, vf, svf, albedo)
            else:
                if fss is None:
                    raise ValueError("Fss is required for non-scattering shortwave")
                knet = _radiation_compute.net_shortwave_nonscattering(
                    zero_sdir, dsky, fss, albedo
                )
            knet_all[:, idx] = np.asarray(knet, dtype=np.float32)

    return sdir_all, knet_all, sveg_all


def generate_timedepsw_from_harmonie(
    *,
    case_dir: Path,
    nwp_root: Path = DEFAULT_NWP_ROOT,
    version: str = DEFAULT_VERSION,
    data_dir: Path | None = None,
    download: bool = False,
    output: Path | None = None,
    sdir_nc: Path | None = None,
    write_sdir_nc: bool = True,
    overwrite: bool = False,
    atmos_only: bool = False,
    method: str | None = None,
    verbose: bool = True,
) -> TimedepShortwaveResult | ShortwaveAtmosphere:
    """Generate ``timedepsw.inp.<expnr>`` from HARMONIE ``ssrd``.

    Set ``atmos_only=True`` to stop after the cheap atmospheric conversion
    without loading geometry or running the facet ray tracer.
    """
    case_dir = Path(case_dir).expanduser().resolve()
    expnr = case_dir.name

    prep = None
    if atmos_only:
        radiation = radiation_config_from_case_dir(case_dir)
    else:
        try:
            from .udprep import UDPrep
        except ImportError as exc:
            raise RuntimeError(
                "Facet mapping requires the uDALES preprocessing environment "
                "(geometry/radiation dependencies such as trimesh, triangle, "
                "netCDF4, and directshortwave_f2py) in addition to the HARMONIE "
                "GRIB dependencies (xarray, cfgrib, pyproj)."
            ) from exc

        prep = UDPrep(expnr, case_dir, load_geometry=True)
        radiation = prep.radiation

    atmosphere, reader = prepare_harmonie_ssrd_atmosphere(
        radiation=radiation,
        nwp_root=nwp_root,
        version=version,
        data_dir=data_dir,
        download=download,
        verbose=verbose,
    )
    if verbose:
        print(f"Spatial mask points: {reader.mask_points}")
        print(
            f"GHI range: {float(atmosphere.ghi.min()):.3f} "
            f"to {float(atmosphere.ghi.max()):.3f} W/m2"
        )
        print(
            f"DNI range: {float(atmosphere.dni.min()):.3f} "
            f"to {float(atmosphere.dni.max()):.3f} W/m2"
        )
        print(
            f"Dsky range: {float(atmosphere.dsky.min()):.3f} "
            f"to {float(atmosphere.dsky.max()):.3f} W/m2"
        )

    if atmos_only:
        return atmosphere

    if prep is None:
        raise RuntimeError("Internal error: UDPrep was not loaded for facet mapping")

    timedepsw_path, default_sdir_nc_path, timedepsveg_path = _default_case_output_paths(
        prep.sim
    )
    if output is not None:
        timedepsw_path = Path(output).expanduser().resolve()
    if sdir_nc is not None:
        default_sdir_nc_path = Path(sdir_nc).expanduser().resolve()
    sdir_nc_path = default_sdir_nc_path if write_sdir_nc else None
    paths_to_check = [timedepsw_path]
    if sdir_nc_path is not None:
        paths_to_check.append(sdir_nc_path)
    _check_output_paths(paths_to_check, overwrite=overwrite)

    sdir, knet, sveg = map_atmosphere_to_facets(
        prep, atmosphere, method=method, verbose=verbose
    )

    write_timedepsw(timedepsw_path, atmosphere.times, knet, overwrite=overwrite)
    if sdir_nc_path is not None:
        prep.radiation._write_sdir_nc(sdir_nc_path, atmosphere.times, sdir)
    written_sveg_path = None
    if sveg is not None:
        _check_output_paths([timedepsveg_path], overwrite=overwrite)
        write_timedepsveg(timedepsveg_path, atmosphere.times, sveg, overwrite=overwrite)
        written_sveg_path = timedepsveg_path

    return TimedepShortwaveResult(
        times=atmosphere.times,
        atmosphere=atmosphere,
        sdir=sdir,
        knet=knet,
        sveg=sveg,
        timedepsw_path=timedepsw_path,
        sdir_nc_path=sdir_nc_path,
        timedepsveg_path=written_sveg_path,
    )


def generate_timedeplw_from_harmonie(
    *,
    case_dir: Path,
    nwp_root: Path = DEFAULT_NWP_ROOT,
    version: str = DEFAULT_VERSION,
    data_dir: Path | None = None,
    download: bool = False,
    output: Path | None = None,
    overwrite: bool = False,
    output_interval_seconds: int | None = None,
    difference_interval_seconds: int | None = None,
    verbose: bool = True,
) -> TimedepLongwaveResult:
    """Generate ``timedeplw.inp.<expnr>`` from HARMONIE ``strd``."""
    case_dir = Path(case_dir).expanduser().resolve()
    config = longwave_config_from_case_dir(case_dir)
    output_path = (
        Path(output).expanduser().resolve()
        if output is not None
        else case_dir / f"timedeplw.inp.{config.expnr}"
    )
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"{output_path} already exists. Pass --overwrite to replace it.")

    series = prepare_harmonie_strd_longwave(
        config=config,
        nwp_root=nwp_root,
        version=version,
        data_dir=data_dir,
        download=download,
        output_interval_seconds=output_interval_seconds,
        difference_interval_seconds=difference_interval_seconds,
        verbose=verbose,
    )

    if verbose:
        print(f"Spatial mask points: {series.mask_points}")
        print(
            f"LWsky range: {float(series.lwsky.min()):.6f} "
            f"to {float(series.lwsky.max()):.6f} W/m2"
        )

    write_timedeplw(output_path, series.times, series.lwsky, overwrite=overwrite)
    return TimedepLongwaveResult(
        times=series.times,
        lwsky=series.lwsky,
        timedeplw_path=output_path,
        mask_points=series.mask_points,
    )
