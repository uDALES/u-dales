from __future__ import annotations

from typing import Any, Callable, Dict, List

import numpy as np

from .udprep_common import SKIP, Section, SectionSpec

FIELDS: List[str] = [
    "itot",
    "xlen",
    "jtot",
    "ylen",
    "ktot",
    "zsize",
    "dx",
    "dy",
    "stl_file",
    "gen_geom",
    "geom_path",
    "diag_neighbs",
    "stl_ground",
    "libm",
    "isolid_bound",
    "ifacsec",
    "read_types",
    "types_path",
    "nfcts",
    "factypes",
    "maxlen",
]

DEFAULTS: Dict[str, Any | Callable[[Any], Any]] = {
    "itot": 64,
    "xlen": 64,
    "jtot": 64,
    "ylen": 64,
    "ktot": 96,
    "zsize": 96,
    "dx": lambda self: self.xlen / self.itot,
    "dy": lambda self: self.ylen / self.jtot,
    "stl_file": "",
    "gen_geom": True,
    "geom_path": "",
    "diag_neighbs": True,
    "stl_ground": True,
    "libm": 1,
    "isolid_bound": 1,
    "ifacsec": 1,
    "read_types": 0,
    "types_path": lambda self: 0 if self.read_types else SKIP,
    "nfcts": 0,
    "factypes": lambda self: SKIP,
    "maxlen": lambda self: 10 if self.lEB else np.inf,
}


class IBMSection(Section):
    def run_all(self) -> None:
        """Run IBM preprocessing steps in the standard order."""
        steps = [
            ("generate_xygrid", self.generate_xygrid),
            ("generate_zgrid", self.generate_zgrid),
            ("write_xgrid", self.write_xgrid),
            ("write_zgrid", self.write_zgrid),
            ("generate_lscale", self.generate_lscale),
            ("write_lscale", self.write_lscale),
            ("generate_factypes", self.generate_factypes),
            ("write_factypes", self.write_factypes),
            ("write_facets", self.write_facets),
            ("write_facetarea", self.write_facetarea),
        ]
        self.run_steps("ibm", steps)

    def generate_xygrid(self) -> None:
        """Create staggered x/y grids for IBM preprocessing."""

    def write_xgrid(self) -> None:
        """Write x-grid definitions to disk."""

    def generate_zgrid(self) -> None:
        """Create vertical grid and stretching based on z parameters."""

    def stretch_exp(self) -> None:
        """Generate exponential stretch profile (stretch_exp in MATLAB)."""

    def stretch_exp_check(self) -> None:
        """Validate exponential stretch (stretch_exp_check in MATLAB)."""

    def stretch_tanh(self) -> None:
        """Generate tanh-based stretch profile (stretch_tanh in MATLAB)."""

    def stretch_2tanh(self) -> None:
        """Generate double-tanh stretch profile (stretch_2tanh in MATLAB)."""

    def write_zgrid(self) -> None:
        """Write z-grid definitions to disk."""

    def generate_lscale(self) -> None:
        """Compute the length scale input (generate_lscale in MATLAB)."""

    def write_lscale(self) -> None:
        """Write the length scale input to disk."""

    def generate_factypes(self) -> None:
        """Construct default factypes array (generate_factypes in MATLAB)."""

    def write_factypes(self) -> None:
        """Write factypes.inp file."""

    def write_facets(self) -> None:
        """Write facets.inp file from STL and facet types."""

    def write_facetarea(self) -> None:
        """Write facetarea.inp file (facet areas)."""

    def compute_boundary_cells(self) -> None:
        """Compute fluid/solid boundary points (IBM getBoundaryCells)."""

    def match_facets_to_cells(self) -> None:
        """Match facets to fluid cells (IBM matchFacetsToCells)."""

    def classify_solid_points(self) -> None:
        """Classify solid points using inpolyhedron/in_mypoly methods."""

    def compute_facet_sections(self) -> None:
        """Compute facet sections for u/v/w/c grids."""

    def compute_fluid_boundary(self) -> None:
        """Compute fluid boundary locations for IBM grids."""


SPEC = SectionSpec(name="ibm", fields=FIELDS, defaults=DEFAULTS, section_cls=IBMSection)
