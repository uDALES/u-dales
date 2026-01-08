from __future__ import annotations

from typing import Any, Callable, Dict, List

from .udprep import Section, SectionSpec

FIELDS: List[str] = [
    "ltrees",
    "ltreesfile",
    "treesfile",
    "tree_dz",
    "tree_dx",
    "tree_dy",
    "tree_h",
    "tree_w",
    "tree_b",
    "nrows",
]

DEFAULTS: Dict[str, Any | Callable[[Any], Any]] = {
    "ltrees": 0,
    "ltreesfile": 0,
    "treesfile": "",
    "tree_dz": 0,
    "tree_dx": 0,
    "tree_dy": 0,
    "tree_h": 0,
    "tree_w": 0,
    "tree_b": 0,
    "nrows": 0,
}

class VegetationSection(Section):
    def run_all(self) -> None:
        """Run vegetation preprocessing steps."""
        steps = [
            ("generate_trees_from_namoptions", self.generate_trees_from_namoptions),
            ("write_trees", self.write_trees),
            ("convert_sparse_veg_inputs", self.convert_sparse_veg_inputs),
        ]
        self.run_steps("vegetation", steps)

    def generate_trees_from_namoptions(self) -> None:
        """Generate trees array from namoptions (generate_trees_from_namoptions)."""

    def write_trees(self) -> None:
        """Write trees.inp file."""

    def plot_trees(self) -> None:
        """Plot tree layout for inspection (plot_trees in MATLAB)."""

    def convert_sparse_veg_inputs(self) -> None:
        """Convert block vegetation to sparse veg inputs (convert_block_to_sparse)."""


SPEC = SectionSpec(
    name="vegetation",
    fields=FIELDS,
    defaults=DEFAULTS,
    section_cls=VegetationSection,
)
