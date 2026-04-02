from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np

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
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        sc = np.zeros((len(sim.zf), int(self.nsv) + 1), dtype=float)
        sc[:, 0] = np.asarray(sim.zf, dtype=float)
        for idx, name in enumerate(("sv10", "sv20", "sv30", "sv40", "sv50"), start=1):
            if idx <= int(self.nsv):
                sc[:, idx] = float(getattr(self, name))
        self.sc = sc

    def write_scalar(self) -> None:
        """Write scalar.inp file."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")
        if not hasattr(self, "sc"):
            self.generate_scalar()

        path = Path(sim.path) / f"scalar.inp.{sim.expnr}"
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# SDBL flow\n")
            f.write("# z scaN,  N=1,2...nsv\n")
            for row in self.sc:
                values = [f"{row[0]:-20.15f}"]
                values.extend(f"{val:-14.10f}" for val in row[1:])
                f.write(" ".join(values) + "\n")

    def generate_scalarsources(self) -> None:
        """Generate scalar sources (generate_scalarsources in MATLAB)."""
        if bool(self.lscasrcr):
            raise ValueError("Network of point sources not currently implemented")

        if bool(self.lscasrc):
            n = int(self.nscasrc)
            src = np.zeros((n, 5), dtype=float)
            if n == 1:
                src[0, :] = [
                    float(self.xS),
                    float(self.yS),
                    float(self.zS),
                    float(self.SSp),
                    float(self.sigSp),
                ]
            self.scasrcp = src

        if bool(self.lscasrcl):
            n = int(self.nscasrcl)
            src = np.zeros((n, 8), dtype=float)
            if n == 1:
                src[0, :] = [
                    float(self.xSb),
                    float(self.ySb),
                    float(self.zSb),
                    float(self.xSe),
                    float(self.ySe),
                    float(self.zSe),
                    float(self.SSl),
                    float(self.sigSl),
                ]
            self.scasrcl = src

    def write_scalarsources(self) -> None:
        """Write scalarsources.inp file."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")
        if bool(self.lscasrc) or bool(self.lscasrcl):
            if not hasattr(self, "scasrcp") and not hasattr(self, "scasrcl"):
                self.generate_scalarsources()

        for ii in range(1, int(self.nsv) + 1):
            if bool(self.lscasrc):
                path = Path(sim.path) / f"scalarsourcep.inp.{ii}.{sim.expnr}"
                with path.open("w", encoding="ascii", newline="\n") as f:
                    f.write("# Scalar point source data\n")
                    f.write("#xS yS zS SS sigS\n")
                    for row in self.scasrcp:
                        f.write(" ".join(f"{val:-12.6f}" for val in row) + "\n")
            if bool(self.lscasrcl):
                path = Path(sim.path) / f"scalarsourcel.inp.{ii}.{sim.expnr}"
                with path.open("w", encoding="ascii", newline="\n") as f:
                    f.write("# Scalar line source data\n")
                    f.write("#xSb ySb zSb xSe ySe zSe SS sigS\n")
                    for row in self.scasrcl:
                        f.write(" ".join(f"{val:-12.6f}" for val in row) + "\n")

    def plot_scalarsources(self) -> None:
        """Plot scalar sources for inspection (plot_scalarsources in MATLAB)."""


SPEC = SectionSpec(name="scalars", fields=FIELDS, defaults=DEFAULTS, section_cls=ScalarsSection)
