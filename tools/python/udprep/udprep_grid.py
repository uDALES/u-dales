"""Grid generation for uDALES preprocessing.

Builds uniform x/y grids and stretched or uniform z grids from
namoptions parameters.  Supports exponential and tanh stretching
schemes with automatic stretch-constant fitting.
"""
from __future__ import annotations

from typing import Any, Dict, List
import numpy as np
from scipy.optimize import brentq
import warnings

from ._section import Section, SectionSpec

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
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        
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
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")

        self.sim.xt = np.arange(0.5 * self.dx, self.xlen, self.dx)
        self.sim.yt = np.arange(0.5 * self.dy, self.ylen, self.dy)
        self.sim.xm = np.arange(0.0, self.xlen, self.dx)
        self.sim.ym = np.arange(0.0, self.ylen, self.dy)

    def generate_zgrid(self) -> None:
        """Create vertical grid.
        
        Generates cell-centered (zt) and cell-faces (zm) vertical grids
        along with cell spacing (dzt). Uses uniform or stretched grids
        depending on configuration flags.
        """
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        
        if not self.lzstretch:
            self.sim.zt = np.arange(0.5 * self.dz, self.zsize, self.dz)
            self.sim.zm = np.arange(0.0, self.zsize, self.dz)
        else:
            # Validate exactly one stretch method is selected
            stretch_methods = {
                "lstretchexp": self.lstretchexp,
                "lstretchexpcheck": self.lstretchexpcheck,
                "lstretchtanh": self.lstretchtanh,
                "lstretch2tanh": self.lstretch2tanh,
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
        self.sim.dzt = np.diff(np.append(self.zm, self.zsize))

    def _set_vertical_grid_from_faces(self, zm: np.ndarray) -> None:
        """Store vertical face/center coordinates from a face grid in real space.
        
        zm must have ktot+1 elements (cell faces including top boundary).
        Produces ktot-length zm (bottom face of each cell), zt, and dzt.
        """
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        
        zm = np.asarray(zm, dtype=float)
        self.sim.zm = zm[:-1]
        self.sim.zt = 0.5 * (zm[:-1] + zm[1:])

    @staticmethod
    def _matlab_round_int(value: float) -> int:
        """Round halves away from zero, matching MATLAB's round for grid counts."""
        return int(np.sign(value) * np.floor(abs(value) + 0.5))

    def _linear_prefix_faces(self) -> tuple[int, int, np.ndarray]:
        """Return the linear near-wall prefix of the vertical face grid (ktot+1 elements)."""
        il = self._matlab_round_int(self.hlin / self.dzlin)
        ir = self.ktot - il
        zm = np.zeros(self.ktot + 1, dtype=float)
        zm[: il + 1] = np.linspace(0.0, self.hlin, il + 1)
        return il, ir, zm

    def _warn_large_top_spacing(self, dz: np.ndarray, reference_dz: float | None = None) -> None:
        if reference_dz is None:
            reference_dz = self.dzlin
        if dz[-1] > 3 * reference_dz:
            warnings.warn("Final grid spacing large; consider reducing domain height", RuntimeWarning)

    @staticmethod
    def _solve_exponential_alpha(ratio: float) -> float:
        """Solve alpha/(exp(alpha)-1)=ratio without converging to the zero root."""
        if not np.isfinite(ratio) or ratio <= 0.0:
            raise ValueError(f"Invalid exponential stretch ratio: {ratio}")

        if np.isclose(ratio, 1.0, rtol=1e-12, atol=1e-12):
            return 0.0

        def func(alpha: float) -> float:
            return alpha - ratio * np.expm1(alpha)

        eps = 1e-12
        if ratio < 1.0:
            lower = eps
            upper = 1.0
            while func(upper) > 0.0:
                upper *= 2.0
                if upper > 700.0:
                    raise ValueError(f"Unable to bracket exponential stretch alpha for ratio={ratio}")
        else:
            lower = -1.0
            upper = -eps
            while func(lower) > 0.0:
                lower *= 2.0
                if lower < -1.0e6:
                    raise ValueError(f"Unable to bracket exponential stretch alpha for ratio={ratio}")

        return float(brentq(func, lower, upper))

    def _stretch_from_computational_coordinate(self, transform) -> None:
        """Map a uniform computational coordinate to a stretched real-space z-grid."""
        il, ir, zm = self._linear_prefix_faces()
        if ir <= 0:
            self._set_vertical_grid_from_faces(zm)
            return

        gf = float(self.stretchconst)
        linear_dz = self.hlin / il if il > 0 else self.dzlin
        xi = np.arange(0, ir + 1, dtype=float) / ir

        while True:
            stretched = zm[il] + (self.zsize - zm[il]) * transform(gf, xi)
            zm[il:] = stretched
            if (zm[il + 1] - zm[il]) < linear_dz:
                gf -= 0.01
                continue
            self._warn_large_top_spacing(np.diff(zm), linear_dz)
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
        il = self._matlab_round_int(self.hlin / self.dzlin)
        ir = self.ktot - il
        z0 = il * self.dzlin  # hlin will be modified as z0
        
        L = self.zsize - z0
        xi = np.linspace(0.0, 1.0, ir + 1)
        
        # Determine alpha using root finding
        # alpha / (exp(alpha)-1) = (dzlin*ir)/L
        ratio = (self.dzlin * ir) / L
        alpha = self._solve_exponential_alpha(ratio)
        
        # Define grid functions
        if alpha == 0.0:
            zhat = xi
        else:
            zhat = np.expm1(alpha * xi) / np.expm1(alpha)
        z = z0 + zhat * L
        
        zm = np.zeros(self.ktot + 1, dtype=float)
        zm[:il + 1] = np.linspace(0.0, z0, il + 1)
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
