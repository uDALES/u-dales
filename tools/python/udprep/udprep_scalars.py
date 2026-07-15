"""Scalar field preprocessing for uDALES.

Generates scalar.inp initial profiles and optional point/line scalar
source files (scalarsources.inp) for passive or reactive tracer
simulations.
"""
from __future__ import annotations

import warnings
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
        """Generate scalar initial conditions"""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")

        if self.nsv > 0:
            if self.nsv > 5:
                raise ValueError("nsv (number of scalar variables) cannot exceed 5 in current implementation")
            
            sc = np.zeros((self.ktot, self.nsv + 1), dtype=float)
            sc[:, 0] = self.zt
            for idx, name in enumerate(("sv10", "sv20", "sv30", "sv40", "sv50"), start=1):
                if idx <= self.nsv:
                    sc[:, idx] = float(getattr(self, name))
            self.sim.sc = sc

    def write_scalar(self) -> None:
        """Write scalar.inp file."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        
        if self.nsv > 0:
            if not hasattr(self.sim, "sc"):
                self.generate_scalar()

            path = Path(self.path) / f"scalar.inp.{self.expnr}"
            with path.open("w", encoding="ascii", newline="\n") as f:
                f.write("# SDBL flow\n")
                f.write("# z scaN,  N=1,2...nsv\n")
                for row in self.sc:
                    values = [f"{row[0]:-20.15f}"]
                    values.extend(f"{val:-14.10f}" for val in row[1:])
                    f.write(" ".join(values) + "\n")

    def generate_scalarsources(self) -> None:
        """Generate scalar sources"""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        
        if self.lscasrcr:
            raise ValueError("Network of point sources not currently implemented... "
                             "only support lscasrc (for point sources) or lscasrcl (for line sources)")

        if (self.lscasrc or self.lscasrcl) and self.nsv < 1:
            raise ValueError( "Must set non-zero positive nsv under &SCALARS, since lscasrc or lscasrcl or both are enabled." )

        if self.lscasrc:
            if (self.nscasrc < 1):
                raise ValueError( "Must set non-zero positive nscasrc under &SCALARS, since lscasrc is set to true." )
            if (self.nscasrc < 2) and (
                self.xS == -1.0
                or self.yS == -1.0
                or self.zS == -1.0
                or self.SSp == -1.0
                or self.sigSp == -1.0 ):
                raise ValueError("Must set appropriate xS, yS, zS, SSp and sigSp under &INPS for a scalar point source.")
            src = np.zeros((self.nscasrc, 5), dtype=float)
            if self.nscasrc == 1:
                src[0, :] = [
                    float(self.xS),
                    float(self.yS),
                    float(self.zS),
                    float(self.SSp),
                    float(self.sigSp),
                ]
            self.sim.scasrcp = src
            if self.nscasrc > 1 or self.nsv > 1:
                warnings.warn( "Manually set appropriate xS, yS, zS, SSp and sigSp "
                               "for scalar source points in scalarsourcep.inp. files",
                               stacklevel=2 )

        if self.lscasrcl:
            if (self.nscasrcl < 1):
                raise ValueError( "Must set non-zero positive nscasrcl under &SCALARS, since lscasrcl is set to true." )
            if (self.nscasrcl < 2) and (
                self.xSb == -1.0
                or self.ySb == -1.0
                or self.zSb == -1.0
                or self.xSe == -1.0
                or self.ySe == -1.0
                or self.zSe == -1.0
                or self.SSl == -1.0
                or self.sigSl == -1.0 ):
                raise ValueError("Must set appropriate xSb, ySb, zSb, xSe, ySe, zSe, SSl and sigSl under &INPS for a scalar line source.")
            src = np.zeros((self.nscasrcl, 8), dtype=float)
            if self.nscasrcl == 1:
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
            self.sim.scasrcl = src
            if self.nscasrcl > 1 or self.nsv > 1:
                warnings.warn( "Manually set appropriate xSb, ySb, zSb, xSe, ySe, zSe, SSl and sigSl "
                               "for scalar source lines in scalarsourcel.inp. files",
                               stacklevel=2 )

    def write_scalarsources(self) -> None:
        """Write scalarsources.inp file."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        
        if self.lscasrc or self.lscasrcl:
            if not hasattr(self.sim, "scasrcp") and not hasattr(self.sim, "scasrcl"):
                self.generate_scalarsources()
            for ii in range(1, self.nsv + 1):
                if self.lscasrc:
                    path = Path(self.path) / f"scalarsourcep.inp.{ii}.{self.expnr}"
                    if path.exists():
                        warnings.warn(f"{path} already exists; NOT overwriting.", stacklevel=2)
                        return
                    with path.open("w", encoding="ascii", newline="\n") as f:
                        f.write("# Scalar point source data\n")
                        f.write("#xS yS zS SS sigS\n")
                        for row in self.scasrcp:
                            f.write(" ".join(f"{val:-12.6f}" for val in row) + "\n")
                if self.lscasrcl:
                    path = Path(self.path) / f"scalarsourcel.inp.{ii}.{self.expnr}"
                    if path.exists():
                        warnings.warn(f"{path} already exists; NOT overwriting.", stacklevel=2)
                        return
                    with path.open("w", encoding="ascii", newline="\n") as f:
                        f.write("# Scalar line source data\n")
                        f.write("#xSb ySb zSb xSe ySe zSe SS sigS\n")
                        for row in self.scasrcl:
                            f.write(" ".join(f"{val:-12.6f}" for val in row) + "\n")
            warnings.warn( "Ensure scalar sources do not intersect any building !! "
                           "Check by plotting the sources along with the geometry using prep.sim.vis.plot_scalar_source(). ", stacklevel=2 )

SPEC = SectionSpec(name="scalars", fields=FIELDS, defaults=DEFAULTS, section_cls=ScalarsSection)
