from __future__ import annotations

from typing import Any, Callable, Dict, List

from .udprep_common import SKIP, Section, SectionSpec

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
    "treesfile": lambda self: "" if (self.ltrees and self.ltreesfile) else SKIP,
    "tree_dz": lambda self: 0 if (self.ltrees and not self.ltreesfile) else SKIP,
    "tree_dx": lambda self: 0 if (self.ltrees and not self.ltreesfile) else SKIP,
    "tree_dy": lambda self: 0 if (self.ltrees and not self.ltreesfile) else SKIP,
    "tree_h": lambda self: 0 if (self.ltrees and not self.ltreesfile) else SKIP,
    "tree_w": lambda self: 0 if (self.ltrees and not self.ltreesfile) else SKIP,
    "tree_b": lambda self: 0 if (self.ltrees and not self.ltreesfile) else SKIP,
    "nrows": lambda self: 0 if (self.ltrees and not self.ltreesfile) else SKIP,
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
