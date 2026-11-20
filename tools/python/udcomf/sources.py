"""Names of the case files the comfort pipeline can read, and one loader for them.

Standard uDALES outputs (``stats_t``, ``stats_tkslice``, ``facEB``) are read
through the UDBase loaders. The warm-start merges of :mod:`udcomf_io` write
their own files beside them, ``merged_stats_tkslice.<expnr>.nc`` and
``merged_facEB.<expnr>.nc``, which hold only what the comfort indices need;
those are read directly here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from udnetcdf import load_ncdata

MERGED_STATISTICS = "merged_stats_tkslice"
MERGED_FACETS = "merged_facEB"
STATISTICS_SOURCES = ("stats_t", "stats_tkslice", MERGED_STATISTICS)
FACET_SOURCES = ("facEB", MERGED_FACETS)


def source_path(sim, prefix: str) -> Path:
    return Path(sim.path) / f"{prefix}.{sim.expnr}.nc"


def load_statistics(sim, source: str, var: str, **kwargs) -> np.ndarray:
    """Read ``var`` from the chosen statistics file of the case."""
    if source not in STATISTICS_SOURCES:
        raise ValueError(f"source must be one of {list(STATISTICS_SOURCES)}")
    if source == "stats_t":
        return sim.load_stat_t(var, **kwargs)
    if source == "stats_tkslice":
        return sim.load_stat_tkslice(var, **kwargs)
    return load_ncdata(source_path(sim, source), var, **kwargs)


def load_facets(sim, source: str, var: str, **kwargs) -> np.ndarray:
    """Read ``var`` from the chosen facet energy-balance file of the case."""
    if source not in FACET_SOURCES:
        raise ValueError(f"source must be one of {list(FACET_SOURCES)}")
    if source == "facEB":
        return sim.load_fac_eb(var, **kwargs)
    return load_ncdata(source_path(sim, source), var, **kwargs)
