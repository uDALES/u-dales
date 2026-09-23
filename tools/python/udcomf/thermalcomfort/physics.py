"""Radiation and directional-sensor physics for model-neutral comfort maps.

The globe/wick heat and mass balances follow Liljegren et al. (2008). Their
flat-ground incoming-radiation model is replaced by saved six-plane irradiance
and a stated, orientation-averaged sensor geometry. This is an urban adaptation,
not a claim of bitwise equivalence to Liljegren's reference WBGT program.
"""

from __future__ import annotations

import numpy as np


SIGMA = 5.670374419e-8
EPSILON = 18.015 / 28.97
_CP = 1003.5
_R_AIR = 8314.34 / 28.97
_PR = _CP / (_CP + 1.25 * _R_AIR)
_RATIO = _CP / EPSILON
_SIDES = ("upface", "downface", "northface", "southface", "eastface", "westface")


def vapour_pressure_pa(qv: np.ndarray, pabs: np.ndarray) -> np.ndarray:
    """Partial pressure from vapour mass fraction of total moist air."""
    return qv * pabs / (EPSILON + (1.0 - EPSILON) * qv)


def saturation_pressure_pa(temperature_k: np.ndarray) -> np.ndarray:
    """Buck liquid-water saturation formula used by Liljegren's wick model."""
    return 100.0 * 1.004 * 6.1121 * np.exp(
        17.502 * (temperature_k - 273.15) / (temperature_k - 32.18)
    )


def direct_body_factor(zenith: float, *, side_factor: float = 0.28,
                       top_factor: float = 0.06) -> float:
    """Orientation-averaged standing-cylinder projected area per body area.

    The lateral coefficient is the standing-person SOLWEIG/VDI convention;
    the overhead coefficient is the six-direction up-face factor. The result
    varies continuously with solar elevation, but not human azimuth.
    """
    if not np.isfinite(zenith) or not 0 <= zenith <= np.pi:
        raise ValueError("Solar zenith must be within [0, pi]")
    if zenith >= np.pi / 2:
        return 0.0
    return side_factor * np.sin(zenith) + top_factor * np.cos(zenith)


def mean_radiant_temperature(fields: dict[str, np.ndarray], zenith: float,
                             *, alpha_sw: float = 0.70, emissivity: float = 0.97,
                             factor_horizontal: float = 0.06,
                             factor_vertical: float = 0.22,
                             direct_side_factor: float = 0.28) -> np.ndarray:
    """Six-plane standing-person MRT in kelvin, with separate direct beam."""
    if not np.isclose(2 * factor_horizontal + 4 * factor_vertical, 1.0, atol=1e-12):
        raise ValueError("Six-direction body factors must sum to one")
    sw = factor_horizontal * (fields["sw_nondirect_upface"] + fields["sw_nondirect_downface"])
    lw = factor_horizontal * (fields["lw_upface"] + fields["lw_downface"])
    for side in _SIDES[2:]:
        sw = sw + factor_vertical * fields[f"sw_nondirect_{side}"]
        lw = lw + factor_vertical * fields[f"lw_{side}"]
    sw = sw + direct_body_factor(zenith, side_factor=direct_side_factor,
                                 top_factor=factor_horizontal) * fields["sw_direct_normal"]
    absorbed = alpha_sw * sw + emissivity * lw
    with np.errstate(invalid="ignore"):
        return np.where(absorbed > 0, (absorbed / (emissivity * SIGMA)) ** 0.25, np.nan)


def _air_viscosity(temperature_k: np.ndarray) -> np.ndarray:
    reduced = temperature_k / 97.0
    omega = 1.048 - 0.034 * (reduced - 2.9) / 0.4
    return 2.6693e-6 * np.sqrt(28.97 * temperature_k) / (3.617**2 * omega)


def _air_conductivity(temperature_k: np.ndarray) -> np.ndarray:
    return (_CP + 1.25 * _R_AIR) * _air_viscosity(temperature_k)


def _vapour_diffusivity(temperature_k: np.ndarray, pressure_hpa: np.ndarray) -> np.ndarray:
    pcrit = (36.4 * 218.0) ** (1.0 / 3.0)
    tcrit = (132.0 * 647.3) ** (5.0 / 12.0)
    troot = np.sqrt(132.0 * 647.3)
    mmix = np.sqrt(1.0 / 28.97 + 1.0 / 18.015)
    return 3.640e-4 * (temperature_k / troot) ** 2.334 * pcrit * tcrit * mmix / (
        pressure_hpa / 1013.25
    ) * 1e-4


def _convective_coefficient(temperature_k: np.ndarray, pressure_pa: np.ndarray,
                            wind: np.ndarray, diameter: float, *, sphere: bool) -> np.ndarray:
    density = pressure_pa / (_R_AIR * temperature_k)
    reynolds = wind * density * diameter / _air_viscosity(temperature_k)
    nusselt = (2.0 + 0.6 * np.sqrt(reynolds) * _PR ** (1.0 / 3.0) if sphere else
               0.281 * reynolds**0.6 * _PR**0.44)
    return nusselt * _air_conductivity(temperature_k) / diameter


def _sensor_irradiance(fields: dict[str, np.ndarray], zenith: float,
                       *, sphere: bool) -> tuple[np.ndarray, np.ndarray]:
    if sphere:
        sw = sum(fields[f"sw_nondirect_{side}"] for side in _SIDES) / 6.0
        lw = sum(fields[f"lw_{side}"] for side in _SIDES) / 6.0
        return sw + fields["sw_direct_normal"] / 4.0, lw
    end_to_side = 0.007 / (4.0 * 0.0254)
    norm = 1.0 + 2.0 * end_to_side
    sw = sum(fields[f"sw_nondirect_{side}"] for side in _SIDES[2:]) / 4.0
    lw = sum(fields[f"lw_{side}"] for side in _SIDES[2:]) / 4.0
    sw = (
        sw + end_to_side * (
            fields["sw_nondirect_upface"] + fields["sw_nondirect_downface"]
        )
    ) / norm
    lw = (lw + end_to_side * (fields["lw_upface"] + fields["lw_downface"])) / norm
    if zenith < np.pi / 2:
        direct = fields["sw_direct_normal"] * (
            np.sin(zenith) / np.pi + end_to_side * np.cos(zenith)
        ) / norm
        sw = sw + direct
    return sw, lw


def sensor_temperatures(ta_k: np.ndarray, pressure_pa: np.ndarray, qv: np.ndarray,
                        wind: np.ndarray, fields: dict[str, np.ndarray], zenith: float,
                        *, iterations: int = 80, tolerance_k: float = 0.02,
                        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Iterate globe and natural-wet-bulb balances on a bounded grid tile.

    No meteorological input is clipped. Cells outside Liljegren's forced-flow
    range (wind < 0.13 m/s) or with a nonconvergent sensor balance return NaN.
    """
    shape = ta_k.shape
    if any(array.shape != shape for array in (pressure_pa, qv, wind)):
        raise ValueError("Sensor input shapes disagree")
    globe_sw, globe_lw = _sensor_irradiance(fields, zenith, sphere=True)
    wick_sw, wick_lw = _sensor_irradiance(fields, zenith, sphere=False)
    e_air = vapour_pressure_pa(qv, pressure_pa)
    valid = (np.isfinite(ta_k) & np.isfinite(pressure_pa) & np.isfinite(wind)
             & np.isfinite(e_air) & (wind >= 0.13) & (e_air >= 0)
             & (e_air < pressure_pa) & (ta_k > 0))
    globe = np.where(valid, ta_k, np.nan).astype(float)
    wick = np.where(valid, ta_k, np.nan).astype(float)
    globe_ok = np.zeros(shape, dtype=bool)
    wick_ok = np.zeros(shape, dtype=bool)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for _ in range(iterations):
            ref = 0.5 * (globe + ta_k)
            h = _convective_coefficient(ref, pressure_pa, wind, 0.0508, sphere=True)
            rhs = (0.95 * globe_lw + 0.95 * globe_sw + h * (ta_k - globe)) / (0.95 * SIGMA)
            proposed = np.where(rhs > 0, rhs**0.25, np.nan)
            globe_ok |= valid & np.isfinite(proposed) & (np.abs(proposed - globe) < tolerance_k)
            globe = np.where(globe_ok, proposed, 0.9 * globe + 0.1 * proposed)
            if np.all(globe_ok | ~valid):
                break
        for _ in range(iterations):
            ref = 0.5 * (wick + ta_k)
            h = _convective_coefficient(ref, pressure_pa, wind, 0.007, sphere=False)
            density = pressure_pa / (_R_AIR * ref)
            schmidt = _air_viscosity(ref) / (
                density * _vapour_diffusivity(ref, pressure_pa / 100.0)
            )
            e_wick = saturation_pressure_pa(wick)
            latent = (2.4073e6 - 71100.0 * (313.15 - ref) / 30.0) / _RATIO
            evaporation = (
                latent * (e_wick - e_air) / (pressure_pa - e_wick)
                * (_PR / schmidt) ** 0.56
            )
            radiative = 0.95 * (wick_lw - SIGMA * wick**4) + 0.60 * wick_sw
            proposed = ta_k - evaporation + radiative / h
            wick_ok |= valid & np.isfinite(proposed) & (e_wick < pressure_pa) & (
                np.abs(proposed - wick) < tolerance_k
            )
            wick = np.where(wick_ok, proposed, 0.9 * wick + 0.1 * proposed)
            if np.all(wick_ok | ~valid):
                break
    return np.where(globe_ok, globe, np.nan), np.where(wick_ok, wick, np.nan), globe_ok & wick_ok


def wet_bulb_globe_temperature(ta_k: np.ndarray, globe_k: np.ndarray,
                               wick_k: np.ndarray, direct_normal: np.ndarray) -> np.ndarray:
    """ISO 7243 solar/no-direct-solar expressions, selected per grid cell."""
    ta_c = ta_k - 273.15
    tg_c = globe_k - 273.15
    tw_c = wick_k - 273.15
    return np.where(direct_normal > 0, 0.7 * tw_c + 0.2 * tg_c + 0.1 * ta_c,
                    0.7 * tw_c + 0.3 * tg_c)
