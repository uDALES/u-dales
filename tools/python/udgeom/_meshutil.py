"""
_meshutil - Shared low-level mesh helpers for udgeom.

Small, byte-identical helpers that were previously duplicated across
``check_mesh``, ``fix_mesh`` and ``truncate_below_ground``. This module is kept
dependency-light (numpy / shapely / trimesh only, no intra-package imports) so
any udgeom module can import from it without introducing an import cycle.
"""

from __future__ import annotations

import numpy as np

from shapely.geometry import GeometryCollection, MultiPolygon, Polygon

try:
    import trimesh
except ImportError:  # pragma: no cover - trimesh availability is handled by callers
    trimesh = None


def _copy_mesh(mesh: "trimesh.Trimesh") -> "trimesh.Trimesh":
    return trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices, dtype=float).copy(),
        faces=np.asarray(mesh.faces, dtype=int).copy(),
        process=False,
    )


def _iter_polygons(geom):
    if geom.is_empty:
        return
    if isinstance(geom, Polygon):
        yield geom
    elif isinstance(geom, MultiPolygon):
        for poly in geom.geoms:
            yield from _iter_polygons(poly)
    elif isinstance(geom, GeometryCollection):
        for item in geom.geoms:
            yield from _iter_polygons(item)


def _project_vertical_face(vertices: np.ndarray, axis: int) -> Polygon | None:
    if axis == 0:
        coords = vertices[:, [1, 2]]
    else:
        coords = vertices[:, [0, 2]]
    poly = Polygon(coords)
    if not poly.is_valid:
        poly = poly.buffer(0)
    if poly.is_empty or poly.area <= 1.0e-12:
        return None
    return poly
