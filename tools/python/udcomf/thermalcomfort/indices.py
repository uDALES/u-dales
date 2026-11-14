"""Model-neutral MRT, PET, UTCI and WBGT from exchange NetCDF files."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib.metadata import version
import json
import operator
import os
from pathlib import Path
import re
import tempfile

from netCDF4 import Dataset
import numpy as np

from .physics import (
    EPSILON, mean_radiant_temperature, saturation_pressure_pa, sensor_temperatures,
    vapour_pressure_pa, wet_bulb_globe_temperature,
)


_ATMOSPHERE = {"ta": "K", "pabs": "Pa", "qv": "kg kg-1",
               "ws_local": "m s-1", "ws_10": "m s-1"}
_RADIATION = {"sw_direct_normal": "W m-2"}
for _side in ("upface", "downface", "northface", "southface", "eastface", "westface"):
    _RADIATION[f"sw_nondirect_{_side}"] = "W m-2"
    _RADIATION[f"lw_{_side}"] = "W m-2"
_FIELDS = {**_ATMOSPHERE, **_RADIATION}
_RESULTS = {"mrt": "mean radiant temperature", "pet": "physiological equivalent temperature",
            "utci": "universal thermal climate index", "wbgt": "wet-bulb globe temperature",
            "globe_temperature": "black-globe sensor temperature",
            "natural_wet_bulb_temperature": "natural wet-bulb sensor temperature"}


@dataclass(frozen=True)
class ComfortParameters:
    """Explicit common person and sensor choices applied to every model.

    PET metabolic activity is 1.37 met = 79.7 W m-2, not 80 W total-body
    power. This distinguishes the library's units from the informal PET
    reference-person wording in the project note.
    """

    human_shortwave_absorptivity: float = 0.70
    human_longwave_emissivity: float = 0.97
    human_horizontal_factor: float = 0.06
    human_vertical_factor: float = 0.22
    human_direct_side_factor: float = 0.28
    pet_met: float = 1.37
    pet_clo: float = 0.9
    pet_age_years: int = 35
    pet_sex: str = "male"
    pet_weight_kg: float = 75.0
    pet_height_m: float = 1.75
    supersaturation_tolerance_percent: float = 5.0

    def __post_init__(self) -> None:
        values = (self.human_shortwave_absorptivity, self.human_longwave_emissivity,
                  self.human_horizontal_factor, self.human_vertical_factor,
                  self.human_direct_side_factor, self.pet_met, self.pet_clo,
                  self.pet_weight_kg, self.pet_height_m,
                  self.supersaturation_tolerance_percent)
        if (not np.isfinite(values).all()
                or not 0 < self.human_shortwave_absorptivity <= 1
                or not 0 < self.human_longwave_emissivity <= 1
                or not np.isclose(2 * self.human_horizontal_factor +
                                  4 * self.human_vertical_factor, 1.0, atol=1e-12)
                or self.human_horizontal_factor <= 0 or self.human_vertical_factor <= 0
                or self.human_direct_side_factor <= 0
                or self.pet_met <= 0 or self.pet_clo < 0 or self.pet_age_years <= 0
                or self.pet_weight_kg <= 0 or self.pet_height_m <= 0
                or self.supersaturation_tolerance_percent < 0
                or self.pet_sex not in ("male", "female")):
            raise ValueError("Comfort person and radiation parameters are invalid")


def _validate_exchange(ds: Dataset) -> np.ndarray | None:
    """Check what the index physics needs; dynamic values are checked by tile.

    Required: the native (x, y) grid, at least one time record on a strictly
    increasing axis, per-record solar angles, and the physical fields with
    the documented units. Everything else - time bounds, grid bounds,
    longitude/latitude, ground elevation, CRS, receptor height, pedestrian
    mask - is optional and validated only when present. Returns the
    pedestrian mask, or None when the file carries none (the caller then
    derives one from data finiteness).
    """
    for dim in ("x", "y", "time"):
        if dim not in ds.dimensions:
            raise ValueError(f"Exchange file lacks dimension {dim}")
    nx, ny, nt = (len(ds.dimensions[name]) for name in ("x", "y", "time"))
    if nt < 1:
        raise ValueError("Exchange file needs at least one time record")
    if "bounds" in ds.dimensions and len(ds.dimensions["bounds"]) != 2:
        raise ValueError("The bounds dimension must have size two")
    schema = {
        "time": ("time",),
        "time_bounds": ("time", "bounds"),
        "x": ("x",), "y": ("y",),
        "x_bounds": ("x", "bounds"),
        "y_bounds": ("y", "bounds"),
        "longitude": ("x", "y"),
        "latitude": ("x", "y"),
        "z_ground": ("x", "y"),
        "receptor_height": (),
        "solar_zenith": ("time",),
        "solar_azimuth": ("time",),
    }
    required = {"time", "x", "y", "solar_zenith", "solar_azimuth"}
    for name, dims in schema.items():
        if name not in ds.variables:
            if name in required:
                raise ValueError(f"Exchange file has missing or invalid {name} dimensions")
            continue
        if ds[name].dimensions != dims:
            raise ValueError(f"Exchange file has missing or invalid {name} dimensions")
    for name in _FIELDS:
        if (name not in ds.variables or ds[name].dimensions != ("x", "y", "time")
                or getattr(ds[name], "units", None) != _FIELDS[name]):
            raise ValueError(f"Exchange file lacks a valid {name} field")
    mask = None
    if "pedestrian_mask" in ds.variables:
        if ds["pedestrian_mask"].dimensions != ("x", "y"):
            raise ValueError("The pedestrian mask must live on the native grid")
        mask = np.asarray(ds["pedestrian_mask"][:])
        if not np.isin(mask, (0, 1)).all() or not np.any(mask):
            raise ValueError("Exchange file has an empty or invalid pedestrian mask")
        mask = mask.astype(bool)
    if "receptor_height" in ds.variables:
        if (not np.isfinite(float(ds["receptor_height"][...]))
                or float(ds["receptor_height"][...]) <= 0):
            raise ValueError("Receptor height must be finite and positive")
    time = np.asarray(ds["time"][:], dtype=float)
    if not np.isfinite(time).all() or (nt > 1 and not np.all(np.diff(time) > 0)):
        raise ValueError("Exchange time axis must be finite and strictly increasing")
    if "time_bounds" in ds.variables:
        bounds = np.asarray(ds["time_bounds"][:], dtype=float)
        if not np.isfinite(bounds).all() or np.any(bounds[:, 0] >= bounds[:, 1]):
            raise ValueError("Exchange time bounds must be finite windows of positive length")
    for axis, n in (("x", nx), ("y", ny)):
        centers = np.asarray(ds[axis][:], dtype=float)
        if centers.shape != (n,) or not np.isfinite(centers).all() or not np.all(np.diff(centers) > 0):
            raise ValueError(f"Exchange {axis} grid and bounds disagree")
        if f"{axis}_bounds" in ds.variables:
            bounds = np.asarray(ds[f"{axis}_bounds"][:], dtype=float)
            if (bounds.shape != (n, 2) or not np.isfinite(bounds).all()
                    or not np.allclose(bounds.mean(axis=1), centers, rtol=0, atol=1e-5)
                    or not np.allclose(bounds[1:, 0], bounds[:-1, 1], rtol=0, atol=1e-5)):
                raise ValueError(f"Exchange {axis} grid and bounds disagree")
    for name in ("longitude", "latitude", "z_ground"):
        if name not in ds.variables:
            continue
        values = np.asarray(ds[name][:], dtype=float)
        if values.shape != (nx, ny) or not np.isfinite(values).all():
            raise ValueError(f"Exchange {name} must be finite on the native grid")
    for name, maximum in (("longitude", 180), ("latitude", 90)):
        if name in ds.variables and np.any(abs(ds[name][:]) > maximum):
            raise ValueError(f"Exchange {name} is outside WGS84 bounds")
    for name, maximum in (("solar_zenith", np.pi), ("solar_azimuth", 2 * np.pi)):
        values = np.asarray(ds[name][:], dtype=float)
        if not np.isfinite(values).all() or np.any((values < 0) | (values > maximum)):
            raise ValueError(f"Exchange {name} is outside the physical range")
    return mask


def _load_physics():
    try:
        from pythermalcomfort.models import pet_steady, utci
    except ImportError as exc:
        raise ImportError(
            "Step 11 requires pythermalcomfort; install tools/python[comfort]"
        ) from exc
    return pet_steady, utci


def _library_values(function, *, result: str, size: int, **kwargs) -> np.ndarray:
    values = np.asarray(getattr(function(**kwargs), result), dtype=float)
    if values.size != size:
        raise ValueError(
            f"Scientific library returned {values.size} {result} values for {size} cells"
        )
    return values.reshape(size)


def _pet_values(pet_steady, tdb: np.ndarray, mrt: np.ndarray, wind: np.ndarray,
                rh: np.ndarray, pabs: np.ndarray, config: ComfortParameters) -> np.ndarray:
    kwargs = dict(tdb=tdb, tr=mrt, v=wind, rh=rh, p_atm=pabs / 100.0,
                  met=config.pet_met, clo=config.pet_clo, age=config.pet_age_years,
                  sex=config.pet_sex, weight=config.pet_weight_kg,
                  height=config.pet_height_m, position="standing")
    try:
        return _library_values(pet_steady, result="pet", size=len(tdb), **kwargs)
    except (ArithmeticError, RuntimeError, ValueError):
        # The same solver is retried per cell only to isolate nonconvergent cells.
        values = np.full(len(tdb), np.nan)
        for i in range(len(tdb)):
            try:
                one = {key: (float(value[i]) if isinstance(value, np.ndarray) else value)
                       for key, value in kwargs.items()}
                values[i] = _library_values(pet_steady, result="pet", size=1, **one)[0]
            except (ArithmeticError, RuntimeError, ValueError):
                continue
        return values


def calculate_indices(input_path: Path, *, output_path: Path | None = None,
                      parameters: ComfortParameters | None = None,
                      overwrite: bool = False, tile_x: int = 64) -> Path:
    """Calculate common-model hourly maps, preserving the input's native grid.

    No horizontal interpolation or model-specific adjustment occurs. The one
    documented exception to "no clipping" is slight supersaturation: relative
    humidity above 100 % but within ``supersaturation_tolerance_percent`` is
    treated as exactly saturated and flagged (bit 32) rather than rejected,
    because averaging the convex saturation curve and interpolating between
    saturated levels legitimately produce it. Larger overshoots remain invalid
    (bit 1). Non-applicable indices are missing and flagged, while valid MRT
    remains available independently of PET/UTCI/WBGT.

    Any number of time records is accepted, from a single instant upward.
    When the file carries no pedestrian mask, the computation mask is
    derived from data finiteness at the first record and written to the
    output. Per-cell nonphysical inputs are flagged (validity bit 64)
    rather than failing the whole file.
    """
    try:
        tile_x = operator.index(tile_x)
    except TypeError as exc:
        raise TypeError("tile_x must be an integer") from exc
    if tile_x <= 0:
        raise ValueError("tile_x must be positive")
    config = ComfortParameters() if parameters is None else parameters
    if not isinstance(config, ComfortParameters):
        raise TypeError("parameters must be a ComfortParameters instance")
    pet_steady, utci = _load_physics()
    source_path = Path(input_path)
    with Dataset(source_path) as source:
        source.set_auto_mask(False)
        mask = _validate_exchange(source)
        raw_name = str(getattr(source, "model_name", "") or "")
        if not raw_name:
            raw_name = source_path.stem
        model = re.sub(r"[^a-z0-9_-]+", "_", raw_name.lower()).strip("_")
        if not model:
            raise ValueError("Model name cannot be used in an output filename")
        nt = len(source.dimensions["time"])
        if "receptor_height" in source.variables:
            height = float(source["receptor_height"][...])
            tag = "." + f"h{height:g}".replace(".", "p")
        else:
            tag = ""
        if mask is None:
            nx0, ny0 = (len(source.dimensions[d]) for d in ("x", "y"))
            derived = np.ones((nx0, ny0), dtype=bool)
            for field in _FIELDS:
                var = source[field]
                values = np.asarray(var[:, :, 0], dtype=float)
                fill = getattr(var, "_FillValue", None)
                if fill is not None and np.isfinite(fill):
                    values = np.where(values == fill, np.nan, values)
                derived &= np.isfinite(values)
            if not derived.any():
                raise ValueError("No computable cells in the exchange file")
            mask = derived
        output = (source_path.with_name(f"thermal_comfort_indices_{model}{tag}.nc")
                  if output_path is None else Path(output_path))
        if output.resolve() == source_path.resolve():
            raise ValueError("Comfort output must not replace its exchange input")
        if output.exists() and not overwrite:
            raise FileExistsError(f"Comfort output exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent)
        os.close(fd)
        temporary = Path(name)
        counts = {name: 0 for name in ("invalid_humidity", "humidity_capped",
                                          "utci_outside_domain", "wbgt_outside_domain",
                                          "wbgt_nonconverged", "pet_failed",
                                          "invalid_input")}
        try:
            with Dataset(temporary, "w", format="NETCDF4") as target:
                nx, ny = mask.shape
                for dim in ("x", "y", "time", "bounds"):
                    if dim in source.dimensions:
                        target.createDimension(dim, len(source.dimensions[dim]))
                for name in ("x", "y", "x_bounds", "y_bounds", "longitude", "latitude",
                             "z_ground", "receptor_height", "pedestrian_mask", "crs",
                             "time", "time_bounds", "solar_zenith", "solar_azimuth"):
                    if name not in source.variables:
                        continue
                    original = source[name]
                    attrs = {
                        key: original.getncattr(key)
                        for key in original.ncattrs() if key != "_FillValue"
                    }
                    copy = target.createVariable(name, original.datatype, original.dimensions,
                                                 fill_value=getattr(original, "_FillValue", None))
                    copy[...] = original[...]
                    copy.setncatts(attrs)
                target.setncatts({name: source.getncattr(name) for name in source.ncattrs()})
                target.input_file = source_path.name
                if "pedestrian_mask" not in source.variables:
                    derived_var = target.createVariable("pedestrian_mask", "i1", ("x", "y"))
                    derived_var[:] = mask.astype(np.int8)
                    derived_var.long_name = ("computation mask derived from finite "
                                             "inputs at the first time record")
                target.comfort_method = (
                    "standing six-direction MRT; MEMI PET; operational UTCI; "
                    "directional Liljegren-type WBGT"
                )
                target.comfort_parameters = json.dumps(asdict(config), sort_keys=True)
                target.pet_utci_library = f"pythermalcomfort {version('pythermalcomfort')}"
                target.wbgt_sensor_convention = (
                    "0.0508 m globe; 0.007 x 0.0254 m wetted cylinder; "
                    "direct-exposed ISO weighting"
                )
                target.wbgt_radiation_convention = (
                    "sphere: six-plane mean plus DNI/4; wick: lateral/end area weighting"
                )
                target.wbgt_model_limit = "forced-convection wind >= 0.13 m/s; no wind clipping"
                target.validity_flag_meanings = (
                    "1 invalid_relative_humidity; 2 utci_outside_domain; "
                    "4 wbgt_outside_domain; 8 wbgt_sensor_nonconverged; 16 pet_failed; "
                    "32 humidity_capped_to_saturation; 64 invalid_input"
                )
                variables = {}
                for field, long_name in _RESULTS.items():
                    var = target.createVariable(field, "f4", ("x", "y", "time"),
                                                fill_value=np.float32(-9.96921e36),
                                                zlib=True, complevel=3,
                                                chunksizes=(min(nx, tile_x), min(ny, 64), 1))
                    var.long_name = long_name
                    var.units = "degree_Celsius"
                    var.averaging_convention = (
                        "index calculated from preceding 15-minute mean inputs"
                    )
                    var.coordinates = "longitude latitude"
                    var.grid_mapping = "crs"
                    variables[field] = var
                flags = target.createVariable("validity_flags", "i2", ("x", "y", "time"),
                                              fill_value=np.int16(-1), zlib=True, complevel=3,
                                              chunksizes=(min(nx, tile_x), min(ny, 64), 1))
                flags.long_name = "bitwise index validity reasons"
                flags.flag_masks = np.array([1, 2, 4, 8, 16, 32, 64], dtype="i2")
                flags.flag_meanings = (
                    "invalid_relative_humidity utci_outside_domain wbgt_outside_domain "
                    "wbgt_sensor_nonconverged pet_failed humidity_capped_to_saturation "
                    "invalid_input"
                )
                for t in range(nt):
                    zenith = float(source["solar_zenith"][t])
                    for lo in range(0, nx, tile_x):
                        hi = min(nx, lo + tile_x)
                        valid = mask[lo:hi]
                        if not valid.any():
                            for var in variables.values():
                                var[lo:hi, :, t] = var._FillValue
                            flags[lo:hi, :, t] = flags._FillValue
                            continue
                        tile_mask = valid
                        usable = valid.copy()
                        fields = {}
                        for field in _FIELDS:
                            var = source[field]
                            values = np.asarray(var[lo:hi, :, t], dtype=float)
                            fill = getattr(var, "_FillValue", None)
                            if fill is not None and np.isfinite(fill):
                                values = np.where(values == fill, np.nan, values)
                            with np.errstate(invalid="ignore"):
                                good = np.isfinite(values)
                                if field in ("ta", "pabs"):
                                    good &= values > 0
                                else:
                                    good &= values >= 0
                            usable &= good
                            fields[field] = np.where(valid, values, np.nan)
                        invalid_input = tile_mask & ~usable
                        valid = usable
                        if zenith >= np.pi / 2:
                            # direct beam with the sun at or below the horizon is
                            # inconsistent input: flag the cells, keep the plane
                            with np.errstate(invalid="ignore"):
                                sunlit = fields["sw_direct_normal"] > 1.0e-6
                            invalid_input = invalid_input | (valid & sunlit)
                            valid = valid & ~sunlit
                        ta = fields["ta"]
                        pressure = fields["pabs"]
                        wind = fields["ws_local"]
                        es = saturation_pressure_pa(ta)
                        e_raw = vapour_pressure_pa(fields["qv"], pressure)
                        rh_raw = 100.0 * e_raw / es
                        # Slight supersaturation is expected from averaging the
                        # convex saturation curve and from vertical interpolation
                        # between saturated levels; within the configured
                        # tolerance the cell is treated as exactly saturated in
                        # every index input, not rejected. Bit 32 marks it, but
                        # only when the overshoot exceeds 0.01 % - float32
                        # storage round-off alone must not raise flags.
                        capped = (np.isfinite(rh_raw) & (rh_raw > 100.0)
                                  & (rh_raw <= 100.0 + config.supersaturation_tolerance_percent)
                                  & (e_raw < pressure) & (es < pressure))
                        humidity_ok = ((np.isfinite(rh_raw) & (rh_raw >= 0)
                                        & (rh_raw <= 100.0) & (e_raw < pressure))
                                       | capped)
                        e = np.where(capped, es, e_raw)
                        rh = np.where(capped, 100.0, rh_raw)
                        qv_saturated = EPSILON * es / (pressure - (1.0 - EPSILON) * es)
                        qv_used = np.where(capped, qv_saturated, fields["qv"])
                        mrt_k = mean_radiant_temperature(
                            fields, zenith, alpha_sw=config.human_shortwave_absorptivity,
                            emissivity=config.human_longwave_emissivity,
                            factor_horizontal=config.human_horizontal_factor,
                            factor_vertical=config.human_vertical_factor,
                            direct_side_factor=config.human_direct_side_factor,
                        )
                        mrt = mrt_k - 273.15
                        state = np.zeros(tile_mask.shape, dtype=np.int16)
                        state[invalid_input] |= 64
                        state[valid & ~humidity_ok] |= 1
                        state[valid & capped & (rh_raw > 100.01)] |= 32
                        tdb = ta - 273.15
                        utci_ok = (valid & humidity_ok & (tdb >= -50) & (tdb <= 50)
                                   & ((mrt - tdb) >= -30) & ((mrt - tdb) <= 70)
                                   & (fields["ws_10"] >= 0.5) & (fields["ws_10"] <= 17)
                                   & (e <= 5000) & np.isfinite(mrt))
                        state[valid & ~utci_ok] |= 2
                        utci_values = np.full(valid.shape, np.nan)
                        if utci_ok.any():
                            utci_values[utci_ok] = _library_values(
                                utci, result="utci", size=int(utci_ok.sum()),
                                tdb=tdb[utci_ok], tr=mrt[utci_ok], v=fields["ws_10"][utci_ok],
                                rh=rh[utci_ok], limit_inputs=True, round_output=False,
                            )
                        state[utci_ok & ~np.isfinite(utci_values)] |= 2
                        pet_ok = valid & humidity_ok & np.isfinite(mrt)
                        pet_values = np.full(valid.shape, np.nan)
                        if pet_ok.any():
                            pet_values[pet_ok] = _pet_values(
                                pet_steady, tdb[pet_ok], mrt[pet_ok], wind[pet_ok],
                                rh[pet_ok], pressure[pet_ok], config,
                            )
                        state[pet_ok & ~np.isfinite(pet_values)] |= 16
                        wbgt_domain = (valid & humidity_ok & (wind >= 0.13)
                                       & (ta >= 283.15) & (ta <= 313.15))
                        state[valid & ~wbgt_domain] |= 4
                        tg = np.full(valid.shape, np.nan)
                        tnwb = np.full(valid.shape, np.nan)
                        wbgt_values = np.full(valid.shape, np.nan)
                        if wbgt_domain.any():
                            sensor_fields = {key: np.where(wbgt_domain, value, np.nan)
                                             for key, value in fields.items() if key in _RADIATION}
                            tg, tnwb, converged = sensor_temperatures(
                                np.where(wbgt_domain, ta, np.nan),
                                np.where(wbgt_domain, pressure, np.nan),
                                np.where(wbgt_domain, qv_used, np.nan),
                                np.where(wbgt_domain, wind, np.nan), sensor_fields, zenith,
                            )
                            state[wbgt_domain & ~converged] |= 8
                            wbgt_values = wet_bulb_globe_temperature(
                                ta, tg, tnwb, fields["sw_direct_normal"]
                            )
                        outputs = {"mrt": mrt, "pet": pet_values, "utci": utci_values,
                                   "wbgt": wbgt_values, "globe_temperature": tg - 273.15,
                                   "natural_wet_bulb_temperature": tnwb - 273.15}
                        for field, values in outputs.items():
                            variables[field][lo:hi, :, t] = np.where(
                                valid & np.isfinite(values), values, variables[field]._FillValue
                            ).astype("f4")
                        flags[lo:hi, :, t] = np.where(tile_mask, state, flags._FillValue)
                        for label, bit in (("invalid_humidity", 1), ("humidity_capped", 32),
                                           ("utci_outside_domain", 2),
                                           ("wbgt_outside_domain", 4), ("wbgt_nonconverged", 8),
                                           ("pet_failed", 16), ("invalid_input", 64)):
                            counts[label] += int(np.count_nonzero(state & bit))
                target.validity_counts = json.dumps(counts, sort_keys=True)
            if output.exists() and not overwrite:
                raise FileExistsError(f"Comfort output appeared during calculation: {output}")
            os.replace(temporary, output)
        finally:
            temporary.unlink(missing_ok=True)
    return output
