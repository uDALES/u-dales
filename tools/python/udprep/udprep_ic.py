from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np

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
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        zf = np.asarray(sim.zf, dtype=float)
        pr = np.zeros((len(zf), 6), dtype=float)
        pr[:, 0] = zf

        lapse = float(getattr(self, "lapse", 0.0))
        if lapse:
            thl = np.zeros(int(sim.ktot), dtype=float)
            thl[0] = float(self.thl0)
            dz = float(sim.zsize) / float(sim.ktot)
            for k in range(int(sim.ktot) - 1):
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
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")
        if not hasattr(self, "pr"):
            self.generate_prof()

        path = Path(sim.path) / f"prof.inp.{sim.expnr}"
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


SPEC = SectionSpec(name="ic", fields=FIELDS, defaults=DEFAULTS, section_cls=ICSection)
