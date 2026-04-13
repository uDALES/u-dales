"""Large-scale forcing and initial condition preprocessing for uDALES.

Generates prof.inp and lscale.inp containing initial profiles of
velocity, temperature, moisture, and TKE together with vertical profiles
of large-scale forcings (pressure-gradient, volume-flow-rate, profile, or
Coriolis forcing, subsidence, and radiation-cooling).
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("forcing", {})
FIELDS: List[str] = list(DEFAULTS.keys())


class ForcingSection(Section):
    def run_all(self) -> None:
        """Run initial condition and large-scale forcing preprocessing steps."""
        steps = [
            ("generate_prof", self.generate_prof),
            ("write_prof", self.write_prof),
            ("generate_lscale", self.generate_lscale),
            ("write_lscale", self.write_lscale),
        ]
        self.run_steps("forcing", steps)

    def generate_prof(self) -> None:
        """Generate initial vertical profiles (generate_prof in MATLAB)."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")

        zf = np.asarray(self.zf, dtype=float)
        pr = np.zeros((len(zf), 6), dtype=float)
        pr[:, 0] = zf

        lapse = float(self.lapse)
        if lapse:
            thl = np.zeros(int(self.ktot), dtype=float)
            thl[0] = float(self.thl0)
            dz = float(self.zsize) / float(self.ktot)
            for k in range(int(self.ktot) - 1):
                thl[k + 1] = thl[k] + lapse * dz
            pr[:, 1] = thl
        else:
            pr[:, 1] = float(self.thl0)

        pr[:, 2] = float(self.qt0)
        pr[:, 3] = float(self.u0)
        pr[:, 4] = float(self.v0)
        pr[:, 5] = float(self.tke)
        self.pr = pr

    def write_prof(self) -> None:
        """Write prof.inp file."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        if not hasattr(self, "pr"):
            self.generate_prof()

        path = Path(self.path) / f"prof.inp.{self.expnr}"
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# SDBL flow \n")
            f.write("# z thl qt u v tke                                          \n")
            for row in self.pr:
                f.write(
                    f"{row[0]:-20.15f} {row[1]:-12.6f} {row[2]:-12.6f} "
                    f"{row[3]:-12.6f} {row[4]:-12.6f} {row[5]:-12.6f}\n"
                )

    def plot_profiles(self) -> None:
        """Plot initial profiles for inspection (plot_profiles in MATLAB)."""

    def generate_lscale(self) -> None:
        """Compute the large-scale forcing array (generate_lscale in MATLAB)."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")

        zf = np.asarray(self.zf, dtype=float)
        ls = np.zeros((len(zf), 10), dtype=float)
        ls[:, 0] = zf
        ls[:, 5] = float(self.w_s)
        ls[:, 6] = float(self.dqtdxls)
        ls[:, 7] = float(self.dqtdyls)
        ls[:, 8] = float(self.dqtdtls)
        ls[:, 9] = float(self.R)

        idriver = int(self.idriver)
        no_forcing = not any([self.luoutflowr, self.lvoutflowr, self.luvolflowr,
                               self.lvvolflowr, self.lprofforc, self.lcoriol])
        ldp = no_forcing and (idriver != 2)
        if ldp:
            warnings.warn(
                "No forcing switch config. setup and not a driven simulation so "
                "initial velocities and/or pressure gradients applied.",
                stacklevel=2,
            )
        forcing_flags = (
            bool(self.luoutflowr or self.lvoutflowr)
            + bool(self.luvolflowr or self.lvvolflowr)
            + bool(self.lprofforc)
            + bool(self.lcoriol)
        )
        if forcing_flags + int(ldp) > 1:
            raise ValueError("More than one forcing type specified")
        if self.lprofforc or self.lcoriol:
            ls[:, 1] = float(self.u0)
            ls[:, 2] = float(self.v0)
        elif ldp:
            ls[:, 3] = float(self.dpdx)
            ls[:, 4] = float(self.dpdy)
        self.ls = ls

    def write_lscale(self) -> None:
        """Write lscale.inp file."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        if not hasattr(self, "ls"):
            self.generate_lscale()

        path = Path(self.path) / f"lscale.inp.{self.expnr}"
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# SDBL flow \n")
            f.write("# z uq vq pqx pqy wfls dqtdxls dqtdyls dqtdtls dthlrad      \n")
            for row in self.ls:
                f.write(
                    f"{row[0]:-20.15f} {row[1]:-12.6f} {row[2]:-12.6f} {row[3]:-12.9f} "
                    f"{row[4]:-12.9f} {row[5]:-15.9f} {row[6]:-12.6f} {row[7]:-12.6f} "
                    f"{row[8]:-12.6f} {row[9]:-17.12f}\n"
                )


SPEC = SectionSpec(name="forcing", fields=FIELDS, defaults=DEFAULTS, section_cls=ForcingSection)
