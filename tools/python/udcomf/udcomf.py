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

    def statistics_times(self, source: str = "stats_t") -> np.ndarray:
        """Saved record times of the chosen statistics file, validated for use.

        ``source`` is ``stats_t``, ``stats_tkslice`` or ``merged_stats_tkslice``.
        """
        from .schedule import validate_statistics_times
        from .sources import load_statistics

        return validate_statistics_times(load_statistics(self.sim, source, "time"), self.schedule)

    def prepare_facet_outputs(self, case_paths):
        """Merge facEB t and LWout of warm-start segments into merged_facEB of this case."""
        from .udcomf_io import merge_facet_outputs

        return merge_facet_outputs(Path(self.sim.path), case_paths)

    def prepare_statistics_outputs(self, case_paths, *, source: str = "stats_t"):
        """Merge the comfort subset of each segment's ``source`` file into merged_stats_tkslice."""
        from .udcomf_io import merge_statistics_outputs

        return merge_statistics_outputs(Path(self.sim.path), case_paths, source=source)

    def export_exchange(self, **kwargs):
        """Write one exchange input file per requested height."""
        from .udcomf_export import write_exchange_heights

        return write_exchange_heights(self.sim, **kwargs)
