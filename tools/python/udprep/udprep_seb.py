from __future__ import annotations

from typing import Any, Dict, List

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("seb", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class SEBSection(Section):
    def run_all(self) -> None:
        """Run SEB preprocessing steps."""
        steps = [
            ("write_Tfacinit", self.write_Tfacinit),
            ("write_Tfacinit_layers", self.write_Tfacinit_layers),
        ]
        self.run_steps("seb", steps)

    def write_Tfacinit(self) -> None:
        """Write initial facet temperatures (write_Tfacinit in preprocessing.m)."""

    def write_Tfacinit_layers(self) -> None:
        """Write initial facet temperatures per layer (write_Tfacinit_layers)."""


SPEC = SectionSpec(name="seb", fields=FIELDS, defaults=DEFAULTS, section_cls=SEBSection)
