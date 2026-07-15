"""Surface energy balance (SEB) preprocessing for uDALES.

Writes initial facet-temperature files (Tfacinit, Tfacinit_layers)
used by the energy balance solver at startup.
"""
from __future__ import annotations

from typing import Any, Dict, List

from pathlib import Path
import numpy as np

from exceptions import DependencyError
from ._section import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("seb", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class SEBSection(Section):
    def run_all(self) -> None:
        """Run SEB preprocessing steps."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        needs_facet_temperature = (sim.lEB or sim.iwallmom == 2 or sim.iwalltemp == 2 or sim.iwallmoist == 2)
        if not needs_facet_temperature:
            return

        if self.lfacTlyrs:
            steps = [("write_Tfacinit_layers", self.write_Tfacinit_layers)]
        else:
            steps = [("write_Tfacinit", self.write_Tfacinit)]
        self.run_steps("seb", steps)

    def write_Tfacinit(self) -> None:
        """Write initial facet temperatures (write_Tfacinit in preprocessing.m)."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        if sim.nfcts <= 0:
            raise ValueError("nfcts must be positive to write Tfacinit")

        Tfacinit = np.full(sim.nfcts, self.facT, dtype=float)

        path = Path(sim.path) / f"Tfacinit.inp.{sim.expnr}"
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# Initial facet temperatures in radiative equilibrium\n")
        with path.open("a", encoding="ascii", newline="\n") as f:
            np.savetxt(f, Tfacinit, fmt="%.4f")

    def write_Tfacinit_layers(self) -> None:
        """Write initial facet temperatures per layer (write_Tfacinit_layers)."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        facT_file = self.facT_file
        if not facT_file:
            raise ValueError("facT_file is required for layered facet temperatures")

        facT_path = Path(facT_file)
        if not facT_path.is_absolute():
            facT_path = Path(sim.path) / facT_path
        if not facT_path.exists():
            raise FileNotFoundError(f"facT_file not found: {facT_path}")

        try:
            from netCDF4 import Dataset
        except ImportError as exc:
            raise DependencyError("netCDF4 is required to read facT_file") from exc

        with Dataset(facT_path, "r") as ds:
            if "T" not in ds.variables:
                raise ValueError("facT_file missing variable 'T'")
            Tfac = ds.variables["T"][:]

        # Tfacinit_layers takes the last layer via Tfac[:, :, -1], which assumes
        # the (nfcts, ntimes, nlayers) axis order. Require exactly 3D so that a
        # differently-shaped array is rejected rather than silently mis-sliced
        # (P27).
        if Tfac.ndim != 3:
            raise ValueError(
                "facT_file variable 'T' must be 3D with axis order "
                f"(nfcts, ntimes, nlayers); got {Tfac.ndim}D shape {tuple(Tfac.shape)}"
            )

        Tfacinit_layers = np.asarray(Tfac[:, :, -1], dtype=float)

        path = Path(sim.path) / f"Tfacinit_layers.inp.{sim.expnr}"
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# Initial facet temperatures in radiative equilibrium\n")
        with path.open("a", encoding="ascii", newline="\n") as f:
            np.savetxt(f, Tfacinit_layers, fmt="%.4f")


SPEC = SectionSpec(name="seb", fields=FIELDS, defaults=DEFAULTS, section_cls=SEBSection)
