"""Namelist (namoptions) parsing helpers for uDALES.

Pure file-in/dict-out parsing, extracted from UDBase so config reading is
separable from the case object. UDBase applies the returned mapping to itself.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Union

# Scalars that must exist even when the namelist omits them.
NAMOPTIONS_SCALAR_DEFAULTS: Dict[str, Any] = {
    "nsv": 0,
    "lscasrc": False,
    "lscasrcl": False,
    "lscasrcr": False,
    "nscasrc": 0,
    "nscasrcl": 0,
}


def parse_value(text: str) -> Union[bool, int, float, str]:
    """Convert one namelist value string to bool / int / float / str."""
    low = text.lower()
    if low == ".true.":
        return True
    if low == ".false.":
        return False
    try:
        return float(text) if ("." in text or "e" in low) else int(text)
    except ValueError:
        return text.strip("'\"")


def parse_namoptions(filepath: Union[str, Path]) -> Dict[str, Any]:
    """Parse a namoptions file into ``{key: typed value}``.

    Namelist-group headers (``&...``), ``!`` comments, and lines without ``=``
    are skipped; inline ``!`` comments are stripped. The scalar defaults in
    :data:`NAMOPTIONS_SCALAR_DEFAULTS` are included unless overridden by the file.
    """
    values: Dict[str, Any] = dict(NAMOPTIONS_SCALAR_DEFAULTS)
    # Force UTF-8 so parsing is deterministic across platforms (Windows would
    # otherwise default to cp1252 and could mis-read non-ASCII bytes).
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("&") or stripped.startswith("!") or "=" not in line:
                continue
            key, _, rest = line.partition("=")
            values[key.strip()] = parse_value(rest.split("!")[0].strip())
    return values
