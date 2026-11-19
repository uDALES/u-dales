"""Person, sensor and domain parameters of the comfort indices.

Every numerical choice of the index calculation that is not a physical
constant lives in ``comfort_defaults.json`` beside this module: the
six-direction body factors and optical properties for MRT, the PET
reference person, the UTCI operational domain, the WBGT globe and wick
sensors with their forced-flow and temperature domain, the sensor solver
settings and the ISO 7243 weights. A user who wants other values copies
that file, edits it and passes the copy to the index calculation. A
parameter file must list every key exactly once; nothing is defaulted
behind the user's back.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
import json
from pathlib import Path

import numpy as np


DEFAULT_PARAMETERS_PATH = Path(__file__).with_name("comfort_defaults.json")


@dataclass(frozen=True)
class ComfortParameters:
    """Explicit person, sensor and domain choices applied to every file.

    PET metabolic activity is in met: 1.37 met = 79.7 W m-2, not 80 W
    total-body power. WBGT weights are ordered (natural wet-bulb, globe, air)
    for the solar-load expression and (natural wet-bulb, globe) for the
    no-direct-solar expression.
    """

    human_shortwave_absorptivity: float
    human_longwave_emissivity: float
    human_horizontal_factor: float
    human_vertical_factor: float
    human_direct_side_factor: float
    supersaturation_tolerance_percent: float
    pet_met: float
    pet_clo: float
    pet_age_years: int
    pet_sex: str
    pet_weight_kg: float
    pet_height_m: float
    pet_position: str
    utci_ta_min_c: float
    utci_ta_max_c: float
    utci_mrt_minus_ta_min_k: float
    utci_mrt_minus_ta_max_k: float
    utci_ws_10_min_m_s: float
    utci_ws_10_max_m_s: float
    utci_vapour_pressure_max_pa: float
    wbgt_globe_diameter_m: float
    wbgt_globe_emissivity: float
    wbgt_wick_diameter_m: float
    wbgt_wick_length_m: float
    wbgt_wick_emissivity: float
    wbgt_wick_shortwave_absorptivity: float
    wbgt_wind_min_m_s: float
    wbgt_ta_min_k: float
    wbgt_ta_max_k: float
    wbgt_max_iterations: int
    wbgt_tolerance_k: float
    wbgt_sun_weights: tuple[float, float, float]
    wbgt_shade_weights: tuple[float, float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "wbgt_sun_weights", tuple(float(w) for w in self.wbgt_sun_weights))
        object.__setattr__(self, "wbgt_shade_weights", tuple(float(w) for w in self.wbgt_shade_weights))
        scalars = tuple(
            getattr(self, field.name) for field in fields(self)
            if field.name not in ("pet_sex", "pet_position", "wbgt_sun_weights", "wbgt_shade_weights")
        )
        if not all(isinstance(value, (int, float)) and not isinstance(value, bool) for value in scalars):
            raise ValueError("Comfort parameters must be numbers")
        if not np.isfinite(np.asarray(scalars, dtype=float)).all():
            raise ValueError("Comfort parameters must be finite")
        if (not 0 < self.human_shortwave_absorptivity <= 1
                or not 0 < self.human_longwave_emissivity <= 1
                or self.human_horizontal_factor <= 0 or self.human_vertical_factor <= 0
                or not np.isclose(2 * self.human_horizontal_factor
                                  + 4 * self.human_vertical_factor, 1.0, atol=1e-12)
                or self.human_direct_side_factor <= 0):
            raise ValueError("MRT body factors or optical properties are invalid")
        if self.supersaturation_tolerance_percent < 0:
            raise ValueError("The supersaturation tolerance must not be negative")
        if (self.pet_met <= 0 or self.pet_clo < 0 or self.pet_age_years <= 0
                or self.pet_weight_kg <= 0 or self.pet_height_m <= 0
                or self.pet_sex not in ("male", "female")
                or self.pet_position not in ("standing", "sitting")):
            raise ValueError("PET person parameters are invalid")
        if (self.utci_ta_min_c >= self.utci_ta_max_c
                or self.utci_mrt_minus_ta_min_k >= self.utci_mrt_minus_ta_max_k
                or not 0 <= self.utci_ws_10_min_m_s < self.utci_ws_10_max_m_s
                or self.utci_vapour_pressure_max_pa <= 0):
            raise ValueError("UTCI domain limits are invalid")
        if (self.wbgt_globe_diameter_m <= 0 or self.wbgt_wick_diameter_m <= 0
                or self.wbgt_wick_length_m <= 0
                or not 0 < self.wbgt_globe_emissivity <= 1
                or not 0 < self.wbgt_wick_emissivity <= 1
                or not 0 <= self.wbgt_wick_shortwave_absorptivity <= 1
                or self.wbgt_wind_min_m_s <= 0
                or self.wbgt_ta_min_k >= self.wbgt_ta_max_k or self.wbgt_ta_min_k <= 0
                or self.wbgt_max_iterations <= 0 or self.wbgt_tolerance_k <= 0):
            raise ValueError("WBGT sensor or domain parameters are invalid")
        for name, weights, count in (("wbgt_sun_weights", self.wbgt_sun_weights, 3),
                                     ("wbgt_shade_weights", self.wbgt_shade_weights, 2)):
            if (len(weights) != count or not np.isfinite(weights).all()
                    or np.any(np.asarray(weights) < 0)
                    or not np.isclose(sum(weights), 1.0, atol=1e-12)):
                raise ValueError(f"{name} must be {count} nonnegative weights summing to one")

    @classmethod
    def from_file(cls, path: Path) -> "ComfortParameters":
        """Read a complete parameter file; every key is required, none extra."""
        path = Path(path)
        with path.open(encoding="utf-8") as stream:
            values = json.load(stream)
        if not isinstance(values, dict):
            raise ValueError(f"Comfort parameter file must contain one object: {path}")
        expected = {field.name for field in fields(cls)}
        unknown = sorted(values.keys() - expected)
        missing = sorted(expected - values.keys())
        if unknown:
            raise ValueError(f"Unknown comfort parameters in {path}: {unknown}")
        if missing:
            raise ValueError(f"Missing comfort parameters in {path}: {missing}")
        return cls(**values)

    @classmethod
    def defaults(cls) -> "ComfortParameters":
        """The values shipped in ``comfort_defaults.json``."""
        return cls.from_file(DEFAULT_PARAMETERS_PATH)

    def to_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True)
