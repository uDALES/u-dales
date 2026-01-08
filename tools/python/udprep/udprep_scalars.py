from __future__ import annotations

from typing import Any, Dict, List

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("scalars", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class ScalarsSection(Section):
    def run_all(self) -> None:
        """Run scalar preprocessing steps."""
        steps = [
            ("generate_scalar", self.generate_scalar),
            ("write_scalar", self.write_scalar),
        ]
        if self.lscasrc or self.lscasrcl:
            steps.extend(
                [
                    ("generate_scalarsources", self.generate_scalarsources),
                    ("write_scalarsources", self.write_scalarsources),
                ]
            )
        self.run_steps("scalars", steps)

    def generate_scalar(self) -> None:
        """Generate scalar initial conditions (generate_scalar in MATLAB)."""

    def write_scalar(self) -> None:
        """Write scalar.inp file."""

    def generate_scalarsources(self) -> None:
        """Generate scalar sources (generate_scalarsources in MATLAB)."""

    def write_scalarsources(self) -> None:
        """Write scalarsources.inp file."""

    def plot_scalarsources(self) -> None:
        """Plot scalar sources for inspection (plot_scalarsources in MATLAB)."""


SPEC = SectionSpec(name="scalars", fields=FIELDS, defaults=DEFAULTS, section_cls=ScalarsSection)
