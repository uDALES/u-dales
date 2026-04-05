from __future__ import annotations

from typing import Any, Dict, List

from pathlib import Path
import numpy as np

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("seb", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class SEBSection(Section):
    def run_all(self) -> None:
        """Run SEB preprocessing steps."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        if not (bool(getattr(self, "lEB", False)) or getattr(sim, "iwallmom", 0) == 2 or getattr(sim, "iwalltemp", 0) == 2):
            return

        if bool(getattr(self, "lfacTlyrs", False)):
            steps = [("write_Tfacinit_layers", self.write_Tfacinit_layers)]
        else:
            steps = [("write_Tfacinit", self.write_Tfacinit)]
        self.run_steps("seb", steps)

    def write_Tfacinit(self) -> None:
        """Write initial facet temperatures (write_Tfacinit in preprocessing.m)."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        nfcts = int(getattr(sim, "nfcts", 0))
        if nfcts <= 0:
            raise ValueError("nfcts must be positive to write Tfacinit")

        facT = float(getattr(self, "facT", 288.0))
        Tfacinit = np.full(nfcts, facT, dtype=float)

        expnr = getattr(sim, "expnr", "")
        out_dir = Path(sim.path) if getattr(sim, "path", None) is not None else Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"Tfacinit.inp.{expnr}"
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# Initial facet tempereatures in radiative equilibrium\n")
        with path.open("a", encoding="ascii", newline="\n") as f:
            np.savetxt(f, Tfacinit, fmt="%4f")

    def write_Tfacinit_layers(self) -> None:
        """Write initial facet temperatures per layer (write_Tfacinit_layers)."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        facT_file = str(getattr(self, "facT_file", ""))
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
            raise ImportError("netCDF4 is required to read facT_file") from exc

        with Dataset(facT_path, "r") as ds:
            if "T" not in ds.variables:
                raise ValueError("facT_file missing variable 'T'")
            Tfac = ds.variables["T"][:]

        if Tfac.ndim < 3:
            raise ValueError("facT_file variable 'T' must be at least 3D")

        Tfacinit_layers = np.asarray(Tfac[:, :, -1], dtype=float)

        expnr = getattr(sim, "expnr", "")
        out_dir = Path(sim.path) if getattr(sim, "path", None) is not None else Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"Tfacinit_layers.inp.{expnr}"
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# Initial facet tempereatures in radiative equilibrium\n")
        with path.open("a", encoding="ascii", newline="\n") as f:
            np.savetxt(f, Tfacinit_layers, fmt="%4f")


SPEC = SectionSpec(name="seb", fields=FIELDS, defaults=DEFAULTS, section_cls=SEBSection)
