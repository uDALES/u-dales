from __future__ import annotations

from typing import Any, Callable, Dict, List

from .udprep_common import SKIP, Section, SectionSpec

FIELDS: List[str] = [
    "u0",
    "v0",
    "tke",
    "thl0",
    "qt0",
    "dpdx",
    "dpdy",
    "lzstretch",
    "stretchconst",
    "lstretchexp",
    "lstretchexpcheck",
    "lstretchtanh",
    "lstretch2tanh",
    "hlin",
    "dzlin",
    "dz",
]

DEFAULTS: Dict[str, Any | Callable[[Any], Any]] = {
    "u0": 0,
    "v0": 0,
    "tke": 0,
    "thl0": 288,
    "qt0": 0,
    "dpdx": 0,
    "dpdy": 0,
    "lzstretch": 0,
    "stretchconst": lambda self: 0.01 if self.lzstretch else SKIP,
    "lstretchexp": lambda self: 0 if self.lzstretch else SKIP,
    "lstretchexpcheck": lambda self: 0 if self.lzstretch else SKIP,
    "lstretchtanh": lambda self: 0 if self.lzstretch else SKIP,
    "lstretch2tanh": lambda self: 0 if self.lzstretch else SKIP,
    "hlin": lambda self: 0 if self.lzstretch else SKIP,
    "dzlin": lambda self: 0 if self.lzstretch else SKIP,
    "dz": lambda self: self.dzlin if self.lzstretch else (self.zsize / self.ktot),
}

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
