from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np
import warnings

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("radiation", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class RadiationSection(Section):
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
            Method-specific optional arguments forwarded to the implementation.
        """
        method_key = method.strip().lower()
        if method_key in ("moller", "raycast", "mt"):
            return self.calc_direct_sw_moller(nsun, irradiance, **kwargs)
        if method_key in ("facsec", "section"):
            return self.calc_direct_sw_facsec(nsun, irradiance, **kwargs)
        if method_key in ("scanline", "f2py"):
            return self.calc_direct_sw_scanline(nsun, irradiance, **kwargs)
        raise ValueError(f"Unknown direct shortwave method: {method}")

    def calc_direct_sw_moller(
        self,
        nsun: np.ndarray,
        irradiance: float,
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """Most accurate and most expensive: ray casting with Moller-Trumbore hits."""
        sim = self._require_sim()
        from .directshortwave_moller import directshortwave

        veg_data = kwargs.pop("veg_data", self._get_veg_data())
        return directshortwave(sim, nsun=nsun, irradiance=irradiance, veg_data=veg_data, **kwargs)

    def calc_direct_sw_facsec(
        self,
        nsun: np.ndarray,
        irradiance: float,
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """Accurate and cheaper: ray casting on cell solids + facsec reconstruction."""
        sim = self._require_sim()
        from .directshortwave_facsec import directshortwave

        veg_data = kwargs.pop("veg_data", self._get_veg_data())
        return directshortwave(sim, nsun=nsun, irradiance=irradiance, veg_data=veg_data, **kwargs)

    def calc_direct_sw_scanline(
        self,
        nsun: np.ndarray,
        irradiance: float,
        ray_density: float = 1.0,
        resolution: float | None = None,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """
        Scanline polygon rasterization via f2py-wrapped Fortran (no vegetation).
        """
        sim = self._require_sim()
        if bool(getattr(sim, "ltree", 0)):
            warnings.warn(
                "Scanline (f2py) direct shortwave does not include vegetation; "
                "results ignore tree attenuation.",
                RuntimeWarning,
            )
        try:
            from udprep.directshortwave_f2py import directshortwave_f2py_mod as _dsmod
        except ImportError as exc:
            raise RuntimeError(
                "directshortwave_f2py module not available; "
                "build it with tools/python/fortran/build_f2py.ps1"
            ) from exc

        if sim.geom is None or sim.geom.stl is None:
            raise ValueError("Geometry not loaded; UDBase must be created with load_geometry=True")
        if ray_density <= 0.0:
            raise ValueError("ray_density must be > 0")

        mesh = sim.geom.stl
        vertices = np.asfortranarray(mesh.vertices, dtype=float)
        faces = np.asfortranarray(mesh.faces, dtype=np.int32) + 1  # 1-based
        incenter = np.asfortranarray(mesh.triangles_center, dtype=float)
        face_normal = np.asfortranarray(mesh.face_normals, dtype=float)
        nsun_f = np.asarray(nsun, dtype=float)
        if resolution is None:
            cell_min = min(sim.dx, sim.dy, float(np.min(sim.dzt)))
            resolution = 0.25 * cell_min / ray_density
        if resolution <= 0.0:
            raise ValueError("resolution must be > 0")

        sdir = _dsmod.calculate_direct_shortwave(
            faces,
            incenter,
            face_normal,
            vertices,
            nsun_f,
            float(irradiance),
            float(resolution),
        )
        sdir = np.asarray(sdir, dtype=float)
        areas = mesh.area_faces
        if hasattr(sim, "facs") and "area" in sim.facs and sim.facs["area"].shape == areas.shape:
            areas = sim.facs["area"]
        bud = {
            "fac": float(np.sum(sdir * areas)),
            "veg": 0.0,
        }
        return sdir, np.zeros(0, dtype=float), bud

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
