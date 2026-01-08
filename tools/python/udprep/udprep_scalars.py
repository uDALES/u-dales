from __future__ import annotations

from typing import Any, Callable, Dict, List

from .udprep import Section, SectionSpec

FIELDS: List[str] = [
    "nsv",
    "sv10",
    "sv20",
    "sv30",
    "sv40",
    "sv50",
    "lscasrc",
    "lscasrcl",
    "lscasrcr",
    "xS",
    "yS",
    "zS",
    "SSp",
    "sigSp",
    "nscasrc",
    "xSb",
    "ySb",
    "zSb",
    "xSe",
    "ySe",
    "zSe",
    "SSl",
    "sigSl",
    "nscasrcl",
]

DEFAULTS: Dict[str, Any | Callable[[Any], Any]] = {
    "nsv": 0,
    "sv10": 0,
    "sv20": 0,
    "sv30": 0,
    "sv40": 0,
    "sv50": 0,
    "lscasrc": 0,
    "lscasrcl": 0,
    "lscasrcr": 0,
    "xS": -1,
    "yS": -1,
    "zS": -1,
    "SSp": -1,
    "sigSp": -1,
    "nscasrc": 0,
    "xSb": -1,
    "ySb": -1,
    "zSb": -1,
    "xSe": -1,
    "ySe": -1,
    "zSe": -1,
    "SSl": -1,
    "sigSl": -1,
    "nscasrcl": 0,
}

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
