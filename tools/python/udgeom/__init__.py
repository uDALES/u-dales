"""
udgeom - Geometry utilities package for uDALES

This package provides geometry manipulation and analysis tools for uDALES simulations,
including mesh operations, building separation, and outline calculation.

Copyright (C) 2016-2025 the uDALES Team.
"""

from .udgeom import UDGeom
from .calculate_outline import calculate_outline
from .check_mesh import (
    calculate_independent_surfaces,
    check,
    find_internal_touching_wall_regions,
    find_nonmanifold_regions,
    find_unstitched_touching_regions,
    identify_ground_faces,
)
from .delete_ground import delete_ground
from .fix_mesh import fix, repair_adjacent_buildings, resolve_vertical_coplanar_overlaps, weld_touching_boundaries
from .split_buildings import split_buildings
from .geometry_generation import (
    add_ground,
    create_flat_surface,
    create_canyons,
    create_cubes,
    plot_ground_generation_step1,
)

__all__ = [
    'UDGeom',
    'calculate_outline',
    'calculate_independent_surfaces',
    'check',
    'identify_ground_faces',
    'find_nonmanifold_regions',
    'find_internal_touching_wall_regions',
    'find_unstitched_touching_regions',
    'delete_ground',
    'extrude_to_ground',
    'fix',
    'repair_adjacent_buildings',
    'resolve_vertical_coplanar_overlaps',
    'split_buildings',
    'truncate_below_ground',
    'weld_touching_boundaries',
    'add_ground',
    'create_flat_surface',
    'create_canyons',
    'create_cubes',
    'plot_ground_generation_step1',
]


def truncate_below_ground(*args, **kwargs):
    from .truncate_below_ground import truncate_below_ground as _truncate_below_ground

    return _truncate_below_ground(*args, **kwargs)


def extrude_to_ground(*args, **kwargs):
    from .extrude_to_ground import extrude_to_ground as _extrude_to_ground

    return _extrude_to_ground(*args, **kwargs)
