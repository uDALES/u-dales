"""Model-neutral MRT, PET, UTCI and WBGT from exchange NetCDF files."""

from __future__ import annotations

from importlib.metadata import version
import json
import os
from pathlib import Path
import tempfile

from netCDF4 import Dataset
import numpy as np

from .parameters import ComfortParameters
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
_HUMIDITY = ("ta", "pabs", "qv")

INDEX_NAMES = ("mrt", "pet", "utci", "wbgt")

# Exchange fields each index reads. MRT itself is an input of PET and UTCI,
# so selecting either of them also evaluates MRT from the radiation planes.
_INDEX_FIELDS = {
    "mrt": tuple(_RADIATION),
    "utci": tuple(_RADIATION) + _HUMIDITY + ("ws_10",),
    "pet": tuple(_RADIATION) + _HUMIDITY + ("ws_local",),
    "wbgt": tuple(_RADIATION) + _HUMIDITY + ("ws_local",),
}
_RESULTS = {"mrt": "mean radiant temperature", "pet": "physiological equivalent temperature",
            "utci": "universal thermal climate index", "wbgt": "wet-bulb globe temperature",
            "globe_temperature": "black-globe sensor temperature",
            "natural_wet_bulb_temperature": "natural wet-bulb sensor temperature"}
_FLAGS = (("invalid_humidity", 1), ("utci_outside_domain", 2), ("wbgt_outside_domain", 4),
          ("wbgt_nonconverged", 8), ("pet_failed", 16), ("humidity_capped", 32),
          ("invalid_input", 64))
_FLAG_MEANINGS = ("invalid_relative_humidity utci_outside_domain wbgt_outside_domain "
                  "wbgt_sensor_nonconverged pet_failed humidity_capped_to_saturation "
                  "invalid_input")
_INPUT_PREFIX = "thermal_comfort_inputs"
_OUTPUT_PREFIX = "thermal_comfort_indices"


def selected_indices(indices) -> tuple[str, ...]:
    """Normalise a selection of index names into canonical order."""
    if indices is None:
        return INDEX_NAMES
    names = tuple(str(name).lower() for name in indices)
    unknown = sorted(set(names) - set(INDEX_NAMES))
    if unknown:
        raise ValueError(f"Unknown comfort indices {unknown}; choose from {list(INDEX_NAMES)}")
    if not names:
        raise ValueError("At least one comfort index must be selected")
    return tuple(name for name in INDEX_NAMES if name in names)


def required_fields(indices: tuple[str, ...]) -> tuple[str, ...]:
    """Exchange fields the selected indices need, in schema order."""
    needed = set()
    for name in indices:
        needed.update(_INDEX_FIELDS[name])
    return tuple(field for field in _FIELDS if field in needed)


def default_output_path(input_path: Path) -> Path:
    """``thermal_comfort_inputs<rest>`` beside the input becomes ``thermal_comfort_indices<rest>``."""
    input_path = Path(input_path)
    if not input_path.name.startswith(_INPUT_PREFIX):
        raise ValueError(
            f"Exchange file name must start with {_INPUT_PREFIX}; the result is written "
            f"beside it with that prefix replaced by {_OUTPUT_PREFIX}: {input_path.name}"
        )
    return input_path.with_name(_OUTPUT_PREFIX + input_path.name[len(_INPUT_PREFIX):])


def _validate_exchange(ds: Dataset, fields: tuple[str, ...]) -> np.ndarray:
    """Check the exchange schema; dynamic values are checked per record.

    Required: the native (x, y) grid with cell bounds, at least one time
    record on a strictly increasing axis with window bounds, the receptor
    height, the pedestrian mask, per-record solar angles, and the physical
    ``fields`` with the documented units. Returns the pedestrian mask.
    """
    for dim in ("x", "y", "time", "bounds"):
        if dim not in ds.dimensions:
            raise ValueError(f"Exchange file lacks dimension {dim}")
    nx, ny, nt = (len(ds.dimensions[name]) for name in ("x", "y", "time"))
    if nt < 1:
        raise ValueError("Exchange file needs at least one time record")
    if len(ds.dimensions["bounds"]) != 2:
        raise ValueError("The bounds dimension must have size two")
    schema = {
        "time": ("time",),
        "time_bounds": ("time", "bounds"),
        "x": ("x",), "y": ("y",),
        "x_bounds": ("x", "bounds"),
        "y_bounds": ("y", "bounds"),
        "receptor_height": (),
        "pedestrian_mask": ("x", "y"),
        "solar_zenith": ("time",),
        "solar_azimuth": ("time",),
    }
    for name, dims in schema.items():
        if name not in ds.variables or ds[name].dimensions != dims:
            raise ValueError(f"Exchange file has missing or invalid {name} dimensions")
    for name in fields:
        if (name not in ds.variables or ds[name].dimensions != ("x", "y", "time")
                or getattr(ds[name], "units", None) != _FIELDS[name]):
            raise ValueError(f"Exchange file lacks a valid {name} field")
    mask = np.asarray(ds["pedestrian_mask"][:])
    if not np.isin(mask, (0, 1)).all() or not np.any(mask):
        raise ValueError("Exchange file has an empty or invalid pedestrian mask")
    mask = mask.astype(bool)
    height = float(ds["receptor_height"][...])
    if not np.isfinite(height) or height <= 0:
        raise ValueError("Receptor height must be finite and positive")
    time = np.asarray(ds["time"][:], dtype=float)
    if not np.isfinite(time).all() or (nt > 1 and not np.all(np.diff(time) > 0)):
        raise ValueError("Exchange time axis must be finite and strictly increasing")
    bounds = np.asarray(ds["time_bounds"][:], dtype=float)
    if not np.isfinite(bounds).all() or np.any(bounds[:, 0] >= bounds[:, 1]):
        raise ValueError("Exchange time bounds must be finite windows of positive length")
    for axis, n in (("x", nx), ("y", ny)):
        centers = np.asarray(ds[axis][:], dtype=float)
        bounds = np.asarray(ds[f"{axis}_bounds"][:], dtype=float)
        if (centers.shape != (n,) or not np.isfinite(centers).all() or not np.all(np.diff(centers) > 0)
                or bounds.shape != (n, 2) or not np.isfinite(bounds).all()
                or not np.allclose(bounds.mean(axis=1), centers, rtol=0, atol=1e-5)
                or not np.allclose(bounds[1:, 0], bounds[:-1, 1], rtol=0, atol=1e-5)):
            raise ValueError(f"Exchange {axis} grid and bounds disagree")
    for name, maximum in (("solar_zenith", np.pi), ("solar_azimuth", 2 * np.pi)):
        values = np.asarray(ds[name][:], dtype=float)
        if not np.isfinite(values).all() or np.any((values < 0) | (values > maximum)):
            raise ValueError(f"Exchange {name} is outside the physical range")
    return mask


def _load_physics(indices: tuple[str, ...]):
    """Import the PET/UTCI library only when one of them is selected."""
    if not ({"pet", "utci"} & set(indices)):
        return None, None, None
    try:
        from pythermalcomfort.models import pet_steady, utci
    except ImportError as exc:
        raise ImportError(
            "PET and UTCI require pythermalcomfort; install tools/python[comfort]"
        ) from exc
    return pet_steady, utci, version("pythermalcomfort")


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
                  height=config.pet_height_m, position=config.pet_position)
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


def _method_description(indices: tuple[str, ...]) -> str:
    parts = {"mrt": "standing six-direction MRT", "pet": "MEMI PET",
             "utci": "operational UTCI", "wbgt": "directional Liljegren-type WBGT"}
    return "; ".join(parts[name] for name in indices)


def calculate_indices(input_path: Path, *, output_path: Path | None = None,
                      parameters: ComfortParameters | None = None,
                      indices=None, overwrite: bool = False) -> Path:
    """Calculate index maps for every record of an exchange file on its native grid.

    ``indices`` selects any subset of ``mrt``, ``pet``, ``utci`` and ``wbgt``
    (default: all four). Only the exchange fields those indices need are
    required and read: MRT needs the radiation planes and the solar zenith;
    UTCI adds air temperature, pressure, humidity and ``ws_10``; PET adds the
    same with ``ws_local``; WBGT needs humidity, pressure, ``ws_local`` and
    the radiation planes but not MRT. Globe and natural wet-bulb temperatures
    are written with WBGT.

    ``parameters`` defaults to the values of ``comfort_defaults.json``.

    No horizontal interpolation or model-specific adjustment occurs. The one
    documented exception to "no clipping" is slight supersaturation: relative
    humidity above 100 % but within ``supersaturation_tolerance_percent`` is
    treated as exactly saturated and flagged (bit 32) rather than rejected,
    because averaging the convex saturation curve and interpolating between
    saturated levels legitimately produce it. Larger overshoots remain invalid
    (bit 1). Non-applicable indices are missing and flagged, while valid MRT
    remains available independently of PET/UTCI/WBGT.

    Any number of time records is accepted, from a single instant upward.
    One record (the full map) is processed at a time; every cell is
    independent, so the result does not depend on that choice. Per-cell
    nonphysical inputs are flagged (validity bit 64) rather than failing the
    whole file.
    """
    config = ComfortParameters.defaults() if parameters is None else parameters
    if not isinstance(config, ComfortParameters):
        raise TypeError("parameters must be a ComfortParameters instance")
    selected = selected_indices(indices)
    fields_needed = required_fields(selected)
    need_mrt = bool({"mrt", "pet", "utci"} & set(selected))
    need_humidity = bool({"pet", "utci", "wbgt"} & set(selected))
    results = [name for name in INDEX_NAMES if name in selected]
    if "wbgt" in selected:
        results += ["globe_temperature", "natural_wet_bulb_temperature"]
    pet_steady, utci, library_version = _load_physics(selected)
    source_path = Path(input_path)
    with Dataset(source_path) as source:
        source.set_auto_mask(False)
        mask = _validate_exchange(source, fields_needed)
        nt = len(source.dimensions["time"])
        output = default_output_path(source_path) if output_path is None else Path(output_path)
        if output.resolve() == source_path.resolve():
            raise ValueError("Comfort output must not replace its exchange input")
        if output.exists() and not overwrite:
            raise FileExistsError(f"Comfort output exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent)
        os.close(fd)
        temporary = Path(name)
        counts = {label: 0 for label, _bit in _FLAGS}
        try:
            with Dataset(temporary, "w", format="NETCDF4") as target:
                nx, ny = mask.shape
                for dim in ("x", "y", "time", "bounds"):
                    target.createDimension(dim, len(source.dimensions[dim]))
                for name in ("x", "y", "x_bounds", "y_bounds", "receptor_height",
                             "pedestrian_mask", "time", "time_bounds",
                             "solar_zenith", "solar_azimuth"):
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
                target.indices = " ".join(selected)
                target.comfort_method = _method_description(selected)
                target.comfort_parameters = config.to_json()
                if library_version is not None:
                    target.pet_utci_library = f"pythermalcomfort {library_version}"
                if "wbgt" in selected:
                    target.wbgt_sensor_convention = (
                        f"{config.wbgt_globe_diameter_m:g} m globe; "
                        f"{config.wbgt_wick_diameter_m:g} x {config.wbgt_wick_length_m:g} m "
                        "wetted cylinder; direct-exposed ISO weighting"
                    )
                    target.wbgt_radiation_convention = (
                        "sphere: six-plane mean plus DNI/4; wick: lateral/end area weighting"
                    )
                    target.wbgt_model_limit = (
                        f"forced-convection wind >= {config.wbgt_wind_min_m_s:g} m/s; "
                        "no wind clipping"
                    )
                target.validity_flag_meanings = (
                    "1 invalid_relative_humidity; 2 utci_outside_domain; "
                    "4 wbgt_outside_domain; 8 wbgt_sensor_nonconverged; 16 pet_failed; "
                    "32 humidity_capped_to_saturation; 64 invalid_input"
                )
                variables = {}
                for field in results:
                    var = target.createVariable(field, "f4", ("x", "y", "time"),
                                                fill_value=np.float32(-9.96921e36),
                                                zlib=True, complevel=3,
                                                chunksizes=(nx, ny, 1))
                    var.long_name = _RESULTS[field]
                    var.units = "degree_Celsius"
                    var.averaging_convention = (
                        "index calculated from the window-mean inputs of each record"
                    )
                    variables[field] = var
                flags = target.createVariable("validity_flags", "i2", ("x", "y", "time"),
                                              fill_value=np.int16(-1), zlib=True, complevel=3,
                                              chunksizes=(nx, ny, 1))
                flags.long_name = "bitwise index validity reasons"
                flags.flag_masks = np.array([bit for _label, bit in _FLAGS], dtype="i2")
                flags.flag_meanings = _FLAG_MEANINGS
                for t in range(nt):
                    zenith = float(source["solar_zenith"][t])
                    valid = mask
                    usable = mask.copy()
                    fields = {}
                    for field in fields_needed:
                        var = source[field]
                        values = np.asarray(var[:, :, t], dtype=float)
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
                    invalid_input = mask & ~usable
                    valid = usable
                    if zenith >= np.pi / 2:
                        # direct beam with the sun at or below the horizon is
                        # inconsistent input: flag the cells, keep the plane
                        with np.errstate(invalid="ignore"):
                            sunlit = fields["sw_direct_normal"] > 1.0e-6
                        invalid_input = invalid_input | (valid & sunlit)
                        valid = valid & ~sunlit
                    state = np.zeros(mask.shape, dtype=np.int16)
                    state[invalid_input] |= 64

                    if need_humidity:
                        ta = fields["ta"]
                        pressure = fields["pabs"]
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
                        tdb = ta - 273.15
                        state[valid & ~humidity_ok] |= 1
                        state[valid & capped & (rh_raw > 100.01)] |= 32

                    outputs = {}
                    if need_mrt:
                        mrt_k = mean_radiant_temperature(
                            fields, zenith, alpha_sw=config.human_shortwave_absorptivity,
                            emissivity=config.human_longwave_emissivity,
                            factor_horizontal=config.human_horizontal_factor,
                            factor_vertical=config.human_vertical_factor,
                            direct_side_factor=config.human_direct_side_factor,
                        )
                        mrt = mrt_k - 273.15
                        if "mrt" in selected:
                            outputs["mrt"] = mrt

                    if "utci" in selected:
                        ws_10 = fields["ws_10"]
                        utci_ok = (valid & humidity_ok
                                   & (tdb >= config.utci_ta_min_c)
                                   & (tdb <= config.utci_ta_max_c)
                                   & ((mrt - tdb) >= config.utci_mrt_minus_ta_min_k)
                                   & ((mrt - tdb) <= config.utci_mrt_minus_ta_max_k)
                                   & (ws_10 >= config.utci_ws_10_min_m_s)
                                   & (ws_10 <= config.utci_ws_10_max_m_s)
                                   & (e <= config.utci_vapour_pressure_max_pa)
                                   & np.isfinite(mrt))
                        state[valid & ~utci_ok] |= 2
                        utci_values = np.full(valid.shape, np.nan)
                        if utci_ok.any():
                            utci_values[utci_ok] = _library_values(
                                utci, result="utci", size=int(utci_ok.sum()),
                                tdb=tdb[utci_ok], tr=mrt[utci_ok], v=ws_10[utci_ok],
                                rh=rh[utci_ok], limit_inputs=True, round_output=False,
                            )
                        state[utci_ok & ~np.isfinite(utci_values)] |= 2
                        outputs["utci"] = utci_values

                    if "pet" in selected:
                        pet_ok = valid & humidity_ok & np.isfinite(mrt)
                        pet_values = np.full(valid.shape, np.nan)
                        if pet_ok.any():
                            pet_values[pet_ok] = _pet_values(
                                pet_steady, tdb[pet_ok], mrt[pet_ok],
                                fields["ws_local"][pet_ok], rh[pet_ok],
                                pressure[pet_ok], config,
                            )
                        state[pet_ok & ~np.isfinite(pet_values)] |= 16
                        outputs["pet"] = pet_values

                    if "wbgt" in selected:
                        wind = fields["ws_local"]
                        wbgt_domain = (valid & humidity_ok
                                       & (wind >= config.wbgt_wind_min_m_s)
                                       & (ta >= config.wbgt_ta_min_k)
                                       & (ta <= config.wbgt_ta_max_k))
                        state[valid & ~wbgt_domain] |= 4
                        tg = np.full(valid.shape, np.nan)
                        tnwb = np.full(valid.shape, np.nan)
                        wbgt_values = np.full(valid.shape, np.nan)
                        if wbgt_domain.any():
                            sensor_fields = {key: np.where(wbgt_domain, value, np.nan)
                                             for key, value in fields.items()
                                             if key in _RADIATION}
                            tg, tnwb, converged = sensor_temperatures(
                                np.where(wbgt_domain, ta, np.nan),
                                np.where(wbgt_domain, pressure, np.nan),
                                np.where(wbgt_domain, qv_used, np.nan),
                                np.where(wbgt_domain, wind, np.nan), sensor_fields, zenith,
                                globe_diameter=config.wbgt_globe_diameter_m,
                                globe_emissivity=config.wbgt_globe_emissivity,
                                wick_diameter=config.wbgt_wick_diameter_m,
                                wick_length=config.wbgt_wick_length_m,
                                wick_emissivity=config.wbgt_wick_emissivity,
                                wick_shortwave_absorptivity=config.wbgt_wick_shortwave_absorptivity,
                                wind_min=config.wbgt_wind_min_m_s,
                                iterations=config.wbgt_max_iterations,
                                tolerance_k=config.wbgt_tolerance_k,
                            )
                            state[wbgt_domain & ~converged] |= 8
                            wbgt_values = wet_bulb_globe_temperature(
                                ta, tg, tnwb, fields["sw_direct_normal"],
                                sun_weights=config.wbgt_sun_weights,
                                shade_weights=config.wbgt_shade_weights,
                            )
                        outputs["wbgt"] = wbgt_values
                        outputs["globe_temperature"] = tg - 273.15
                        outputs["natural_wet_bulb_temperature"] = tnwb - 273.15

                    for field, values in outputs.items():
                        variables[field][:, :, t] = np.where(
                            valid & np.isfinite(values), values, variables[field]._FillValue
                        ).astype("f4")
                    flags[:, :, t] = np.where(mask, state, flags._FillValue)
                    for label, bit in _FLAGS:
                        counts[label] += int(np.count_nonzero(state & bit))
                target.validity_counts = json.dumps(counts, sort_keys=True)
            if output.exists() and not overwrite:
                raise FileExistsError(f"Comfort output appeared during calculation: {output}")
            os.replace(temporary, output)
        finally:
            temporary.unlink(missing_ok=True)
    return output
