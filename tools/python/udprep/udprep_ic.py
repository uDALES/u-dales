from __future__ import annotations

from typing import Any, Dict, List

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("ic", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class ICSection(Section):
    def run_all(self) -> None:
        """Run initial condition preprocessing steps."""
        steps = [
            ("generate_prof", self.generate_prof),
            ("write_prof", self.write_prof),
        ]
        self.run_steps("ic", steps)

    def generate_prof(self) -> None:
        """Generate initial vertical profiles (generate_prof in MATLAB)."""

    def write_prof(self) -> None:
        """Write prof.inp file."""

    def plot_profiles(self) -> None:
        """Plot initial profiles for inspection (plot_profiles in MATLAB)."""


SPEC = SectionSpec(name="ic", fields=FIELDS, defaults=DEFAULTS, section_cls=ICSection)
