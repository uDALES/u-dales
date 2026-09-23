"""Facade for thermal-comfort postprocessing of a UDBase case."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from udbase import UDBase


class UDComf:
    """Hold the existing case state for future comfort calculations."""

    def __init__(self, sim: UDBase):
        self.sim = sim
        from .udcomf_radiation import UDComfRadiation
        from .udcomf_atmosphere import UDComfAtmosphere

        self.radiation = UDComfRadiation(sim)
        self.atmosphere = UDComfAtmosphere(sim)

    def prepare_facet_outputs(self, case_paths, *, overwrite: bool = False):
        """Merge facEB/facT warm starts into this case for comfort postprocessing."""
        from .udcomf_io import merge_facet_outputs

        return merge_facet_outputs(Path(self.sim.path), case_paths, overwrite=overwrite)

    def export_exchange(self, analysis_day, metadata, **kwargs):
        """Write one strict model-neutral input file per requested height."""
        from .udcomf_export import write_exchange_heights

        return write_exchange_heights(self.sim, analysis_day, metadata, **kwargs)
