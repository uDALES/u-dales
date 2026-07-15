"""Boundary condition preprocessing for uDALES."""
from __future__ import annotations

from typing import Any, Dict, List

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("bcs", {})
FIELDS: List[str] = list(DEFAULTS.keys())

SPEC = SectionSpec(name="bcs", fields=FIELDS, defaults=DEFAULTS, section_cls=Section)
