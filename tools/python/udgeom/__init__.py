"""
udgeom - Geometry utilities package for uDALES

This package provides geometry manipulation and analysis tools for uDALES simulations,
including mesh operations, building separation, and outline calculation.

Copyright (C) 2016-2025 the uDALES Team.
"""

from .udgeom import UDGeom
from .calculate_outline import calculate_outline
from .split_buildings import split_buildings

__all__ = [
    'UDGeom',
    'calculate_outline',
    'split_buildings'
]
