from __future__ import annotations

from typing import Any, Callable, Dict, List

from .udprep import Section, SectionSpec

FIELDS: List[str] = [
    "iwallmom",
    "iwalltemp",
    "lwritefac",
    "xazimuth",
    "ltimedepsw",
    "ishortwave",
    "isolar",
    "runtime",
    "dtSP",
    "solarazimuth",
    "solarzenith",
    "I",
    "Dsky",
    "longitude",
    "latitude",
    "timezone",
    "elevation",
    "hour",
    "minute",
    "second",
    "year",
    "month",
    "day",
    "weatherfname",
    "psc_res",
    "lvfsparse",
    "calc_vf",
    "maxD",
    "vf_path",
    "view3d_out",
]

DEFAULTS: Dict[str, Any | Callable[[Any], Any]] = {
    "iwallmom": 3,
    "iwalltemp": 1,
    "lwritefac": 0,
    "xazimuth": 90,
    "ltimedepsw": 0,
    "ishortwave": 1,
    "isolar": 1,
    "runtime": 0,
    "dtSP": 10.0,
    "solarazimuth": 135,
    "solarzenith": 28.4066,
    "I": 800,
    "Dsky": 418.8041,
    "longitude": -0.13,
    "latitude": 51.5,
    "timezone": 0,
    "elevation": 0,
    "hour": 6,
    "minute": 0,
    "second": 0,
    "year": 2011,
    "month": 9,
    "day": 30,
    "weatherfname": "",
    "psc_res": 0.1,
    "lvfsparse": False,
    "calc_vf": True,
    "maxD": float("inf"),
    "vf_path": "",
    "view3d_out": 0,
}

class RadiationSection(Section):
    def run_all(self) -> None:
        """Run radiation preprocessing steps."""
        steps = [
            ("compute_facet_areas", self.compute_facet_areas),
            ("compute_direct_shortwave", self.compute_direct_shortwave),
            ("compute_view_factors", self.compute_view_factors),
            ("write_vf", self.write_vf),
            ("write_svf", self.write_svf),
            ("write_netsw", self.write_netsw),
            ("write_timedepsw", self.write_timedepsw),
        ]
        self.run_steps("radiation", steps)

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


SPEC = SectionSpec(
    name="radiation",
    fields=FIELDS,
    defaults=DEFAULTS,
    section_cls=RadiationSection,
)
