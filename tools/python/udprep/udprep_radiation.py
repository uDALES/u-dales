from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np
from .udprep import Section, SectionSpec
from .directshortwave import DirectShortwaveSolver

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("radiation", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class RadiationSection(Section):
    def __init__(
        self,
        name: str,
        values: Dict[str, Any],
        sim: Any | None = None,
        defaults: Dict[str, Any] | None = None,
    ) -> None:
        super().__init__(name, values, sim=sim, defaults=defaults)
        self._direct_sw_solver: DirectShortwaveSolver | None = None
        self._direct_sw_solver_cfg: Dict[str, Any] | None = None

    def run_all(self) -> None:
        """Run radiation preprocessing steps."""
        steps = [
            ("compute_direct_shortwave", self.compute_direct_shortwave),
            ("compute_view_factors", self.compute_view_factors),
            ("write_vf", self.write_vf),
            ("write_svf", self.write_svf),
            ("write_netsw", self.write_netsw),
            ("write_timedepsw", self.write_timedepsw),
        ]
        self.run_steps("radiation", steps)

    def calc_direct_sw(
        self,
        nsun: np.ndarray,
        irradiance: float,
        method: str = "facsec",
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """
        Compute direct shortwave using the requested method.

        Parameters
        ----------
        nsun : array-like, shape (3,)
            Vector pointing toward the sun.
        irradiance : float
            Direct normal irradiance [W/m^2].
        method : {"moller", "facsec", "scanline"}
            Dispatches to the corresponding implementation.
            - "moller": ray casting with Moller-Trumbore triangle hits (most accurate, most expensive)
            - "facsec": ray casting with solid mask + facet-section reconstruction (accurate, cheaper)
            - "scanline": f2py scanline rasterization (no vegetation, fastest)
        kwargs : dict
            Optional arguments that define the solver configuration:
            - ray_density (float): ray grid density (default 4)
            - periodic_xy (bool): periodic boundaries in x/y (default False)
            - ray_jitter (float): jitter factor for ray positions (default 1)
            - veg_data (dict): vegetation data (default: loaded from sim)
            - resolution (float): scanline resolution override (scanline only)
            The solver is cached; changing any of these options recreates it.
        """
        method_key = method.strip().lower()
        if method_key in ("moller", "raycast", "mt"):
            method_key = "moller"
        elif method_key in ("facsec", "section"):
            method_key = "facsec"
        elif method_key in ("scanline", "f2py"):
            method_key = "scanline"
        else:
            raise ValueError(f"Unknown direct shortwave method: {method}")

        veg_data = kwargs.pop("veg_data", self._get_veg_data())
        ray_density = kwargs.pop("ray_density", 4.0)
        periodic_xy = kwargs.pop("periodic_xy", False)
        ray_jitter = kwargs.pop("ray_jitter", 1.0)
        resolution = kwargs.pop("resolution", None)
        if kwargs:
            raise ValueError(f"Unknown direct shortwave options: {', '.join(sorted(kwargs.keys()))}")

        cfg = {
            "method": method_key,
            "ray_density": float(ray_density),
            "ray_jitter": float(ray_jitter),
            "veg_key": id(veg_data) if veg_data is not None else None,
        }
        if self._direct_sw_solver is None or self._direct_sw_solver_cfg != cfg:
            sim = self._require_sim()
            self._direct_sw_solver = DirectShortwaveSolver(
                sim,
                method_key,
                ray_density=ray_density,
                ray_jitter=ray_jitter,
                veg_data=veg_data,
            )
            self._direct_sw_solver_cfg = cfg

        return self._direct_sw_solver.compute(
            nsun,
            irradiance,
            periodic_xy=periodic_xy,
            resolution=resolution,
        )

    def calc_direct_sw_moller(
        self,
        nsun: np.ndarray,
        irradiance: float,
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """Most accurate and most expensive: ray casting with Moller-Trumbore hits."""
        return self.calc_direct_sw(nsun, irradiance, method="moller", **kwargs)

    def calc_direct_sw_facsec(
        self,
        nsun: np.ndarray,
        irradiance: float,
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """Accurate and cheaper: ray casting on cell solids + facsec reconstruction."""
        return self.calc_direct_sw(nsun, irradiance, method="facsec", **kwargs)

    def calc_direct_sw_scanline(
        self,
        nsun: np.ndarray,
        irradiance: float,
        ray_density: float = 4.0,
        resolution: float | None = None,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """
        Scanline polygon rasterization via f2py-wrapped Fortran (no vegetation).
        """
        return self.calc_direct_sw(
            nsun,
            irradiance,
            method="scanline",
            ray_density=ray_density,
            resolution=resolution,
        )

    def solar_position(self) -> None:
        """Compute solar position (SPA/solarPosition.m)."""

    def ashrae_coeffs(self) -> None:
        """Load ASHRAE coefficients (ASHRAE.m)."""

    def compute_direct_shortwave(self) -> None:
        """Compute direct shortwave on facets (directShortwave.m)."""

    def compute_net_shortwave(self) -> None:
        """Compute net shortwave radiation (netShortwave.m)."""

    def compute_facet_areas(self) -> None:
        """Compute facet areas (facetAreas.m)."""

    def export_view3d_geometry(self) -> None:
        """Export geometry for View3D (STLtoView3D.m)."""

    def compute_view_factors(self) -> None:
        """Compute view factors (View3D external)."""

    def write_vf(self) -> None:
        """Write full view factors (write_vf in preprocessing.m)."""

    def write_vfsparse(self) -> None:
        """Write sparse view factors (write_vfsparse in preprocessing.m)."""

    def write_svf(self) -> None:
        """Write sky view factors (write_svf in preprocessing.m)."""

    def write_netsw(self) -> None:
        """Write net shortwave (write_netsw in preprocessing.m)."""

    def write_timedepsw(self) -> None:
        """Write time-dependent shortwave (write_timedepsw in preprocessing.m)."""

    def poly2mask_ids(self) -> None:
        """Rasterize polygons for IDs (poly2maskIDs.m)."""

    def _require_sim(self):
        if self.sim is None:
            raise ValueError("RadiationSection requires a UDBase instance (sim).")
        return self.sim

    def _get_veg_data(self) -> Dict[str, Any] | None:
        sim = self._require_sim()
        ltree = bool(getattr(sim, "ltrees", False))
        if not ltree:
            return None
        if hasattr(sim, "veg") and sim.veg is not None:
            return sim.veg
        if not hasattr(sim, "load_veg"):
            raise ValueError("Vegetation data not available; sim.load_veg() is required")
        return sim.load_veg(zero_based=True, cache=True)


SPEC = SectionSpec(
    name="radiation",
    fields=FIELDS,
    defaults=DEFAULTS,
    section_cls=RadiationSection,
)
