from __future__ import annotations

from typing import Any, Callable, Dict, List

from .udprep_common import SKIP, Section, SectionSpec

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
    "sv10": lambda self: 0 if self.nsv > 0 else SKIP,
    "sv20": lambda self: 0 if self.nsv > 0 else SKIP,
    "sv30": lambda self: 0 if self.nsv > 0 else SKIP,
    "sv40": lambda self: 0 if self.nsv > 0 else SKIP,
    "sv50": lambda self: 0 if self.nsv > 0 else SKIP,
    "lscasrc": lambda self: 0 if self.nsv > 0 else SKIP,
    "lscasrcl": lambda self: 0 if self.nsv > 0 else SKIP,
    "lscasrcr": lambda self: 0 if self.nsv > 0 else SKIP,
    "xS": lambda self: -1 if self.nsv > 0 else SKIP,
    "yS": lambda self: -1 if self.nsv > 0 else SKIP,
    "zS": lambda self: -1 if self.nsv > 0 else SKIP,
    "SSp": lambda self: -1 if self.nsv > 0 else SKIP,
    "sigSp": lambda self: -1 if self.nsv > 0 else SKIP,
    "nscasrc": lambda self: 0 if self.nsv > 0 else SKIP,
    "xSb": lambda self: -1 if self.nsv > 0 else SKIP,
    "ySb": lambda self: -1 if self.nsv > 0 else SKIP,
    "zSb": lambda self: -1 if self.nsv > 0 else SKIP,
    "xSe": lambda self: -1 if self.nsv > 0 else SKIP,
    "ySe": lambda self: -1 if self.nsv > 0 else SKIP,
    "zSe": lambda self: -1 if self.nsv > 0 else SKIP,
    "SSl": lambda self: -1 if self.nsv > 0 else SKIP,
    "sigSl": lambda self: -1 if self.nsv > 0 else SKIP,
    "nscasrcl": lambda self: 0 if self.nsv > 0 else SKIP,
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
