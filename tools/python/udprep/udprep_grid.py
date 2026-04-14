"""Grid generation for uDALES preprocessing.

Builds uniform x/y grids and stretched or uniform z grids from
namoptions parameters.  Supports exponential and tanh stretching
schemes with automatic stretch-constant fitting.
"""
from __future__ import annotations

from typing import Any, Dict, List
import numpy as np
from scipy.optimize import fsolve
import warnings

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("grid", {})
FIELDS: List[str] = list(DEFAULTS.keys())


class GridSection(Section):
    def run_all(self) -> None:
        """Run grid generation preprocessing steps in the standard order."""
        self._refresh_derived_grid_params()
        steps = [
            ("generate_xygrid", self.generate_xygrid),
            ("generate_zgrid", self.generate_zgrid),
        ]
        self.run_steps("grid", steps)

    def _refresh_derived_grid_params(self) -> None:
        """Recompute grid-derived quantities from the active case settings."""
        self.sim.dx = float(self.xlen / self.itot)
        self.sim.dy = float(self.ylen / self.jtot)
        self.sim.dz = float(self.zsize / self.ktot)
        if self.lzstretch:
            if self.dzlin is None:
                self.dzlin = self.dz
                warnings.warn(f"dzlin has not been set in namoptions while creating stretched z-grid; using default value (zsize/ktot) = {self.dz}")
            if self.hlin is None:
                self.hlin = 0.1 * self.zsize
                warnings.warn(f"hlin has not been set in namoptions while creating stretched z-grid; using default value (0.1*zsize) = {self.hlin}")

    def generate_xygrid(self) -> None:
        """Create staggered x/y grids for preprocessing.
        
        Generates cell-centered (xt, yt) and cell-face (xm, ym) grids
        based on domain size and grid spacing.
        """
        self.sim.xt = np.arange(0.5 * self.dx, self.xlen, self.dx)
        self.sim.yt = np.arange(0.5 * self.dy, self.ylen, self.dy)
        self.sim.xm = np.arange(0, self.xlen, self.dx)
        self.sim.ym = np.arange(0, self.ylen, self.dy)

    def generate_zgrid(self) -> None:
        """Create vertical grid.
        
        Generates cell-centered (zt) and cell-faces (zm) vertical grids
        along with cell spacing (dzt). Uses uniform or stretched grids
        depending on configuration flags.
        """
        if not self.lzstretch:
            self.sim.zt = np.arange(0.5 * self.dz, self.zsize, self.dz)
            self.sim.zm = np.arange(0, self.zsize + self.dz, self.dz)
            self.sim.dzt = self.sim.zm[1:] - self.sim.zm[:-1]
        else:
            # Validate exactly one stretch method is selected
            stretch_methods = {
                "lstretchexp": bool(self.lstretchexp),
                "lstretchexpcheck": bool(self.lstretchexpcheck),
                "lstretchtanh": bool(self.lstretchtanh),
                "lstretch2tanh": bool(self.lstretch2tanh),
            }
            active = [name for name, flag in stretch_methods.items() if flag]
            if len(active) == 0:
                raise ValueError(
                    "lzstretch is true but no stretch method is selected. "
                    "Set exactly one of: lstretchexp, lstretchexpcheck, lstretchtanh, lstretch2tanh."
                )
            if len(active) > 1:
                raise ValueError(
                    f"lzstretch is true but multiple stretch methods are active: {active}. "
                    "Set exactly one of: lstretchexp, lstretchexpcheck, lstretchtanh, lstretch2tanh."
                )
            # Stretched grid - call appropriate stretching method
            if self.lstretchexp:
                self._stretch_exp()
            elif self.lstretchexpcheck:
                self._stretch_exp_check()
            elif self.lstretchtanh:
                self._stretch_tanh()
            elif self.lstretch2tanh:
                self._stretch_2tanh()
        # Add derived z quantities
        self.sim.zf = self.sim.zt.copy()
        self.sim.zh = self.sim.zm.copy()
        self.sim.dzm = np.concatenate([[2 * self.sim.zt[0]], np.diff(self.sim.zt)])

    def _set_vertical_grid_from_faces(self, zm: np.ndarray) -> None:
        """Store vertical face/center coordinates from a face grid in real space."""
        self.sim.zm = np.asarray(zm, dtype=float)
        self.sim.zt = 0.5 * (self.sim.zm[:-1] + self.sim.zm[1:])
        self.sim.dzt = np.diff(self.sim.zm)

    def _linear_prefix_faces(self) -> tuple[int, int, np.ndarray]:
        """Return the linear near-wall prefix of the vertical face grid."""
        il = int(round(self.hlin / self.dzlin))
        ir = int(self.ktot) - il
        zm = np.zeros(int(self.ktot) + 1, dtype=float)
        zm[: il + 1] = np.arange(0.0, self.hlin + self.dzlin, self.dzlin)
        return il, ir, zm

    def _warn_large_top_spacing(self, dz: np.ndarray) -> None:
        if dz[-1] > 3 * self.dzlin:
            warnings.warn("Final grid spacing large; consider reducing domain height", RuntimeWarning)

    def _stretch_from_computational_coordinate(self, transform) -> None:
        """Map a uniform computational coordinate to a stretched real-space z-grid."""
        il, ir, zm = self._linear_prefix_faces()
        if ir <= 0:
            self._set_vertical_grid_from_faces(zm)
            return

        gf = float(self.stretchconst)
        xi = np.arange(0, ir + 1, dtype=float) / ir

        while True:
            stretched = zm[il] + (self.zsize - zm[il]) * transform(gf, xi)
            zm[il:] = stretched
            if (zm[il + 1] - zm[il]) < self.dzlin:
                gf -= 0.01
                continue
            self._warn_large_top_spacing(np.diff(zm))
            break

        self._set_vertical_grid_from_faces(zm)

    def _stretch_exp(self) -> None:
        """Generate exponential stretch grid.
        
        Creates a vertical grid with linear spacing in the lower part and 
        exponentially stretched spacing in the upper part.
        """
        self._stretch_from_computational_coordinate(
            lambda gf, xi: (np.exp(gf * xi) - 1.0) / (np.exp(gf) - 1.0)
        )

    def _stretch_exp_check(self) -> None:
        """Generate exponential stretch grid after validating exponential stretch.
        
        Creates a vertical grid with exponential stretching using a root-finding
        procedure to ensure smooth grid transition and quality checks.
        """
        il = round(self.hlin / self.dzlin)
        ir = self.ktot - il
        z0 = il * self.dzlin  # hlin will be modified as z0
        
        L = self.zsize - z0
        dxi = 1 / ir
        xi = np.arange(0, 1 + dxi, dxi)
        
        # Determine alpha using root finding
        # alpha / (exp(alpha)-1) = (dzlin*ir)/L
        def func(alpha):
            return alpha - (self.dzlin * ir) / L * (np.exp(alpha) - 1)
        
        alpha = fsolve(func, 1.0)[0]
        A = 1 / (np.exp(alpha) - 1)
        
        # Define grid functions
        zhat = A * (np.exp(alpha * xi) - 1)
        z = z0 + zhat * L
        
        zm = np.zeros(self.ktot + 1, dtype=float)
        zm[:il + 1] = np.arange(0, z0 + self.dzlin, self.dzlin)
        zm[il:] = z
        
        # Perform grid quality checks
        dz = np.diff(zm)
        stretch = dz[1:] / dz[:-1]
        
        if (np.min(stretch) < 0.95 or np.max(stretch) > 1.05):
            warnings.warn(
                "Generated grid is of bad quality; stretching factor dz(n+1)/dz(n) "
                f"should stay between 0.95 and 1.05 (min={np.min(stretch):.3f}, "
                f"max={np.max(stretch):.3f})",
                RuntimeWarning,
            )
        
        # Warning if grid is refined near the top
        if alpha < 0:
            warnings.warn(
                "Calculated alpha is negative, implying refinement towards the domain top.",
                RuntimeWarning,
            )

        self._set_vertical_grid_from_faces(zm)

    def _stretch_tanh(self) -> None:
        """Generate tanh-based stretch grid.
        
        Creates a vertical grid with linear spacing in the lower part and
        tanh-based stretched spacing in the upper part.
        """
        self._stretch_from_computational_coordinate(
            lambda gf, xi: 1.0 - np.tanh(gf * (1.0 - xi)) / np.tanh(gf)
        )

    def _stretch_2tanh(self) -> None:
        """Generate double-tanh stretch grid.
        
        Creates a vertical grid with linear spacing in the lower part and
        double-tanh stretched spacing in the upper part.
        """
        self._stretch_from_computational_coordinate(
            lambda gf, xi: 0.5 * (1.0 - np.tanh(gf * (1.0 - 2.0 * xi)) / np.tanh(gf))
        )


SPEC = SectionSpec(name="grid", fields=FIELDS, defaults=DEFAULTS, section_cls=GridSection)
