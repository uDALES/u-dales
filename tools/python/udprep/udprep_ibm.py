from __future__ import annotations

from typing import Any, Dict, List

from .udprep import Section, SectionSpec

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("ibm", {})
FIELDS: List[str] = list(DEFAULTS.keys())


class IBMSection(Section):
    def run_all(self) -> None:
        """Run IBM preprocessing steps in the standard order."""
        steps = [
            ("generate_lscale", self.generate_lscale),
            ("write_lscale", self.write_lscale),
            ("generate_factypes", self.generate_factypes),
            ("write_factypes", self.write_factypes),
            ("write_facets", self.write_facets),
            ("write_facetarea", self.write_facetarea),
        ]
        self.run_steps("ibm", steps)

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
