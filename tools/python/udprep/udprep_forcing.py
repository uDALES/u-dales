"""Large-scale forcing and initial condition preprocessing for uDALES.

Generates prof.inp and lscale.inp containing initial profiles of
velocity, temperature, moisture, and TKE together with vertical profiles
of large-scale forcings (pressure-gradient, volume-flow-rate, profile, or
Coriolis forcing, subsidence, and radiation-cooling).
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from scipy.interpolate import CubicSpline
from ._section import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("forcing", {})
FIELDS: List[str] = list(DEFAULTS.keys())


class ForcingSection(Section):
    def run_all(self) -> None:
        """Run initial condition and large-scale forcing preprocessing steps."""
        if self.idriver == 2:
            steps = [
                ("generate_prof", self.generate_prof),
                ("write_prof", self.write_prof),
                ("update_prof", self._update_prof_from_driver),
                ("generate_lscale", self.generate_lscale),
                ("write_lscale", self.write_lscale),
            ]
        elif self.ltimedepnudge:
            steps = [
                ("generate_prof", self.generate_prof),
                ("write_prof", self.write_prof),
                ("update_prof", self._update_prof_from_nudge_data),
                ("generate_lscale", self.generate_lscale),
                ("write_lscale", self.write_lscale),
            ]
        else:
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

        pr = np.zeros((self.ktot, 6), dtype=float)
        
        pr[:, 0] = self.zt

        lapse = float(self.lapse)
        if lapse:
            thl = np.zeros(self.ktot, dtype=float)
            thl[0] = float(self.thl0)
            for k in range(self.ktot - 1):
                thl[k + 1] = thl[k] + lapse * (self.dzt[k] + self.dzt[k+1]) * 0.5
            pr[:, 1] = thl
        else:
            pr[:, 1] = float(self.thl0)

        pr[:, 2] = float(self.qt0)
        pr[:, 3] = float(self.u0)
        pr[:, 4] = float(self.v0)
        pr[:, 5] = float(self.tke)
        
        self.sim.pr = pr

    def _update_prof_from_driver(self) -> None:
        """Update self.sim.pr with driver simulation output if idriver=2."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        self.update_prof_from_driver(self.driverjobnr, self.driveroutpath, self.drivertimeidx)

    def _update_prof_from_nudge_data(self) -> None:
        """Update self.sim.pr with nudging data."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        self.update_prof_from_nudge_data(self.profsourcefile)

    def update_prof_from_nudge_data(
            self,
            profsourcefile: str) -> None:
        """Update prof.inp with nudging data."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")

        def _warn_keep_original(reason: str) -> None:
            warnings.warn(
                f"Profile sourcefile {profsourcefile!r} is {reason}; "
                "original prof.inp is kept without updating although time dependent nudging is true.",
                stacklevel=1,
            )
        
        path = Path(self.path) / f"prof.inp.{self.expnr}"
        if path.exists():
            pr = self.sim.load_prof()
        else:
            raise FileNotFoundError(f"prof.inp file {path} not found for updating from nudge data.")
        
        if profsourcefile is None or not str(profsourcefile).strip():
            _warn_keep_original("empty or not set")
            return

        path = Path(profsourcefile)
        if not path.is_absolute():
            path = Path(self.path) / path
        if not path.is_file():
            _warn_keep_original(f"not a valid file at {path}")
            return

        try:
            prdata = self.sim.read_matrix(path, 1)
        except Exception as exc:
            _warn_keep_original(f"invalid ({exc})")
            return

        prdata = np.asarray(prdata, dtype=float)
        if prdata.size == 0:
            _warn_keep_original("empty")
            return
        prdata = np.atleast_2d(prdata)
        if prdata.shape[1] < 5:
            _warn_keep_original(f"invalid; expected at least 5 columns, found {prdata.shape[1]}")
            return
        prdata = prdata[:, :5]

        if prdata[0, 0] > 0:
            surface_thl = float(self.thl0)
            prdata = np.vstack(([0.0, 0.0, 0.0, surface_thl, 0.0], prdata))
        if prdata.shape[0] < 2 or np.any(np.diff(prdata[:, 0]) <= 0):
            _warn_keep_original("invalid; height column must contain at least two strictly increasing values")
            return

        try:
            pr[:, 1] = CubicSpline(prdata[:, 0], prdata[:, 3])(pr[:, 0])
            pr[:, 2] = CubicSpline(prdata[:, 0], prdata[:, 4])(pr[:, 0])
            pr[:, 3] = CubicSpline(prdata[:, 0], prdata[:, 1])(pr[:, 0])
            pr[:, 4] = CubicSpline(prdata[:, 0], prdata[:, 2])(pr[:, 0])
        except Exception as exc:
            _warn_keep_original(f"invalid ({exc})")
            return

        warnings.warn(f"Using profile sourcefile {path} for prof.inp generation with time dependent nudging.", stacklevel=1)
        self.sim.pr = pr
        self.write_prof(force=True)  # overwrite prof.inp with updated profiles

    def update_prof_from_driver(
        self,
        driverjobnr: str,
        driveroutpath: str,
        drivertimeidx: Optional[int]) -> None:
        """Overwrite prof.inp columns 1-5 with time-averaged driver simulation output.

        Loads xytdump.<driverjobnr>.nc from *driveroutpath* and replaces the thl,
        qt, u, v, and tke columns of *pr* with the slice at *drivertimeidx*.
        Issues a warning and leaves *pr* unchanged if the file is missing or the
        index is out of range.
        """
        from udbase import UDBase

        if self.sim is None:
            raise ValueError("UDBase instance must be provided")

        path = Path(self.path) / f"prof.inp.{self.expnr}"
        if path.exists():
            pr = self.sim.load_prof()
        else:
            raise FileNotFoundError(f"prof.inp file {path} not found for updating from driver output.")
        
        path = Path(driveroutpath) / f"xytdump.{driverjobnr}.nc"
        if not path.exists():
            warnings.warn(
                f"Driver output file {path} not found; original prof.inp is kept without updating.",
                stacklevel=1,
            )
            return

        simdriver = UDBase(driverjobnr, driveroutpath, load_geometry=False, suppress_load_warnings=True)
        u   = simdriver.load_stat_xyt('uxyt')
        v   = simdriver.load_stat_xyt('vxyt')
        thl = simdriver.load_stat_xyt('thlxyt')
        qt  = simdriver.load_stat_xyt('qtxyt')
        tke = simdriver.load_stat_xyt('tketxyc')

        if drivertimeidx is not None and 0 < drivertimeidx <= u.shape[1]:
            warnings.warn(
                f"Using driver simulation output xytdump.{driverjobnr}.nc data "
                "for prof.inp generation.",
                stacklevel=1,
            )
            pr[:, 1] = thl[:, drivertimeidx - 1]
            pr[:, 2] = qt[:,  drivertimeidx - 1]
            pr[:, 3] = u[:,   drivertimeidx - 1]
            pr[:, 4] = v[:,   drivertimeidx - 1]
            pr[:, 5] = tke[:, drivertimeidx - 1]
        else:
            warnings.warn(
                "drivertimeidx is not set or out of bounds for driver output; "
                "original prof.inp is kept without updating.",
                stacklevel=1,
            )
            return

        self.sim.pr = pr
        self.write_prof(force=True)  # overwrite prof.inp with updated profiles

    def write_prof(self, force: bool = False) -> None:
        """Write prof.inp file."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")

        path = Path(self.path) / f"prof.inp.{self.expnr}"
        if path.exists() and not force:
            warnings.warn(f"{path} already exists; NOT overwriting.", stacklevel=1)
            return
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# SDBL flow \n")
            f.write("# z thl qt u v tke                                          \n")
            for row in self.pr:
                f.write(
                    f"{row[0]:-20.15f} {row[1]:-12.6f} {row[2]:-12.6f} "
                    f"{row[3]:-12.6f} {row[4]:-12.6f} {row[5]:-12.6f}\n"
                )

    def generate_lscale(self) -> None:
        """Compute the large-scale forcing array (generate_lscale in MATLAB)."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")

        ls = np.zeros((self.ktot, 10), dtype=float)

        ls[:, 0] = self.zt
        ls[:, 5] = float(self.w_s)
        ls[:, 6] = float(self.dqtdxls)
        ls[:, 7] = float(self.dqtdyls)
        ls[:, 8] = float(self.dqtdtls)
        ls[:, 9] = float(self.R)

        idriver = int(self.idriver)
        no_forcing = not any([self.luoutflowr, self.lvoutflowr, self.luvolflowr,
                               self.lvvolflowr, self.lprofforc, self.lcoriol, self.lnudge])
        
        ldp = no_forcing and (idriver != 2)
        if ldp:
            warnings.warn(
                "No forcing switch config. setup and not a driven simulation so "
                "initial velocities and/or pressure gradients applied.",
                stacklevel=1,
            )

        forcing_flags = (
            bool(self.luoutflowr or self.lvoutflowr)
            + bool(self.luvolflowr or self.lvvolflowr)
            + bool(self.lprofforc)
            + bool(self.lcoriol)
        )

        if forcing_flags + int(ldp) > 1:
            raise ValueError("More than one forcing type specified in namoptions which is not allowed")
        
        if self.lprofforc or self.lcoriol:
            ls[:, 1] = float(self.u0)
            ls[:, 2] = float(self.v0)
        elif ldp:
            ls[:, 3] = float(self.dpdx)
            ls[:, 4] = float(self.dpdy)
        
        self.sim.ls = ls

    def write_lscale(self) -> None:
        """Write lscale.inp file."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        
        if not hasattr(self.sim, "ls"):
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
