from __future__ import annotations

from typing import Any, Callable, Dict, List

from .udprep_common import SKIP, Section, SectionSpec

FIELDS: List[str] = [
    "lEB",
    "lfacTlyrs",
    "facT",
    "nfaclyrs",
    "facT_file",
    "dtEB",
]

DEFAULTS: Dict[str, Any | Callable[[Any], Any]] = {
    "lEB": 0,
    "lfacTlyrs": 0,
    "facT": 288.0,
    "nfaclyrs": 3,
    "facT_file": "",
    "dtEB": lambda self: 10.0 if self.lEB else SKIP,
}

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
