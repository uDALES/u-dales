from __future__ import annotations

from typing import Any, Callable, Dict, List

from .udprep_common import SKIP, Section, SectionSpec

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
    "xazimuth": lambda self: 90 if self.lEB else SKIP,
    "ltimedepsw": lambda self: 0 if self.lEB else SKIP,
    "ishortwave": lambda self: 1 if self.lEB else SKIP,
    "isolar": lambda self: 1 if self.lEB else SKIP,
    "runtime": lambda self: 0 if self.lEB else SKIP,
    "dtSP": lambda self: (self.dtEB if self.lEB else SKIP),
    "solarazimuth": lambda self: 135 if (self.lEB and self.isolar == 1) else SKIP,
    "solarzenith": lambda self: 28.4066 if (self.lEB and self.isolar == 1) else SKIP,
    "I": lambda self: 800 if (self.lEB and self.isolar == 1) else SKIP,
    "Dsky": lambda self: 418.8041 if (self.lEB and self.isolar == 1) else SKIP,
    "longitude": lambda self: -0.13 if (self.lEB and self.isolar == 2) else SKIP,
    "latitude": lambda self: 51.5 if (self.lEB and self.isolar == 2) else SKIP,
    "timezone": lambda self: 0 if (self.lEB and self.isolar == 2) else SKIP,
    "elevation": lambda self: 0 if (self.lEB and self.isolar == 2) else SKIP,
    "hour": lambda self: 6 if (self.lEB and self.isolar == 2) else (0 if (self.lEB and self.isolar == 3) else SKIP),
    "minute": lambda self: 0 if (self.lEB and self.isolar in (2, 3)) else SKIP,
    "second": lambda self: 0 if (self.lEB and self.isolar in (2, 3)) else SKIP,
    "year": lambda self: 2011 if (self.lEB and self.isolar == 2) else (0 if (self.lEB and self.isolar == 3) else SKIP),
    "month": lambda self: 9 if (self.lEB and self.isolar == 2) else (6 if (self.lEB and self.isolar == 3) else SKIP),
    "day": lambda self: 30 if (self.lEB and self.isolar == 2) else (1 if (self.lEB and self.isolar == 3) else SKIP),
    "weatherfname": lambda self: "" if (self.lEB and self.isolar == 3) else SKIP,
    "psc_res": lambda self: 0.1 if self.lEB else SKIP,
    "lvfsparse": lambda self: False if self.lEB else SKIP,
    "calc_vf": lambda self: True if self.lEB else SKIP,
    "maxD": lambda self: float("inf") if self.lEB else SKIP,
    "vf_path": lambda self: "" if (self.lEB and not self.calc_vf) else SKIP,
    "view3d_out": lambda self: 0 if self.lEB else SKIP,
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
