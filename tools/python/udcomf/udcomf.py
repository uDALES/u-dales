"""Facade for thermal-comfort postprocessing of a UDBase case."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from udbase import UDBase


class UDComf:
    """Hold the existing case state for comfort postprocessing."""

    def __init__(self, sim: UDBase):
        self.sim = sim
        from .udcomf_radiation import UDComfRadiation
        from .udcomf_atmosphere import UDComfAtmosphere

        self.radiation = UDComfRadiation(sim)
        self.atmosphere = UDComfAtmosphere(sim)

    @property
    def schedule(self):
        """Statistics averaging schedule (tstatstart, tstatsdump, tstatsgap, tsample)."""
        from .schedule import read_stats_schedule

        return read_stats_schedule(Path(self.sim.path))

    def statistics_times(self, source: str = "stats_tkslice") -> np.ndarray:
        """Saved record times of the gathered statistics file, validated for use."""
        from .schedule import validate_statistics_times

        if source not in ("stats_t", "stats_tkslice"):
            raise ValueError("source must be stats_t or stats_tkslice")
        load = self.sim.load_stat_t if source == "stats_t" else self.sim.load_stat_tkslice
        return validate_statistics_times(load("time"), self.schedule)

    def prepare_facet_outputs(self, case_paths, *, overwrite: bool = False):
        """Merge facEB/facT warm starts into this case for comfort postprocessing."""
        from .udcomf_io import merge_facet_outputs

        return merge_facet_outputs(Path(self.sim.path), case_paths, overwrite=overwrite)

    def prepare_statistics_outputs(self, case_paths, *, overwrite: bool = False):
        """Merge the comfort subset of warm-start statistics into stats_tkslice of this case."""
        from .udcomf_io import merge_statistics_outputs

        return merge_statistics_outputs(Path(self.sim.path), case_paths, overwrite=overwrite)

    def export_exchange(self, metadata, **kwargs):
        """Write one exchange input file per requested height."""
        from .udcomf_export import write_exchange_heights

        return write_exchange_heights(self.sim, metadata, **kwargs)
