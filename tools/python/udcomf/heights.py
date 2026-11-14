"""Validated receptor heights from the solver's OUTPUT namelist."""

from __future__ import annotations

from pathlib import Path

import f90nml
import numpy as np


def configured_receptor_heights(case_dir: Path) -> tuple[float, ...]:
    """Return the solver's effective wind-height list without editing the case."""
    case = Path(case_dir)
    output = f90nml.read(case / f"namoptions.{case.name}").get("output", {})
    primary = float(output.get("receptor_height", 1.1))
    count = int(output.get("nreceptor_heights", 0))
    if count < 0 or count > 1000:
        raise ValueError("nreceptor_heights must be between 0 and 1000")
    if count:
        raw = output.get("receptor_heights", [])
        values = raw if isinstance(raw, list) else [raw]
        if len(values) != count or any(value is None for value in values):
            raise ValueError("receptor_heights must contain exactly nreceptor_heights values")
        heights = tuple(float(value) for value in values)
        tolerance = 100 * np.finfo(np.float32).eps * max(1.0, abs(primary))
        if not np.isclose(heights[0], primary, rtol=0, atol=tolerance):
            raise ValueError("First receptor_heights entry must equal receptor_height")
    else:
        heights = (primary,)
    validate_heights(heights)
    return heights


def validate_heights(heights) -> tuple[float, ...]:
    """Validate a nonempty increasing sequence of physical heights in metres."""
    values = tuple(float(height) for height in heights)
    if not values or not np.isfinite(values).all() or any(height <= 0 for height in values):
        raise ValueError("Receptor heights must be finite and positive")
    if any(right <= left for left, right in zip(values, values[1:])):
        raise ValueError("Receptor heights must be strictly increasing")
    return values


def height_tag(height: float) -> str:
    """A round-trip-unique, filesystem-safe tag for one positive float height."""
    value = validate_heights((height,))[0]
    return "h" + repr(value).replace("-", "m").replace("+", "").replace(".", "p")
