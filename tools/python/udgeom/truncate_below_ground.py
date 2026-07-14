"""
truncate_below_ground - Trim geometry below a planar identified ground surface.

This targets the common CFD import case where buildings extend below a flat
ground plane. The helper removes the old ground, clips all non-ground faces to
the identified ground elevation, clips the existing ground locally against the
building footprint, and welds the result back into a stitched surface mesh.
"""

from __future__ import annotations

from typing import Dict, Tuple, Union

import numpy as np
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon
from shapely.ops import triangulate

from .check_mesh import _as_trimesh, identify_ground_faces
from .delete_ground import delete_ground
try:
    import trimesh
    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False


def _copy_mesh(mesh: "trimesh.Trimesh") -> "trimesh.Trimesh":
    return trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices, dtype=float).copy(),
        faces=np.asarray(mesh.faces, dtype=int).copy(),
        process=False,
    )


def _estimate_planar_ground_level(mesh: "trimesh.Trimesh", ground_mask: np.ndarray, tolerance: float) -> float:
    ground_vertices = np.unique(np.asarray(mesh.faces, dtype=int)[ground_mask].reshape(-1))
    z = np.asarray(mesh.vertices, dtype=float)[ground_vertices, 2]
    if len(z) == 0:
        raise ValueError("No ground vertices were identified.")
    z0 = float(np.median(z))
    if float(np.max(np.abs(z - z0))) > tolerance:
        raise NotImplementedError(
            "truncate_below_ground currently supports only near-planar identified ground surfaces."
        )
    return z0


def _estimate_ground_spacing(mesh: "trimesh.Trimesh", ground_mask: np.ndarray) -> float:
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)[ground_mask]
    if len(faces) == 0:
        return 1.0

    lengths = []
    for face in faces:
        tri = vertices[face][:, :2]
        for a, b in ((0, 1), (1, 2), (2, 0)):
            length = float(np.linalg.norm(tri[a] - tri[b]))
            if length > 1.0e-8:
                lengths.append(length)
    if not lengths:
        return 1.0
    return float(np.median(lengths))


def _clip_triangle_above_z(tri: np.ndarray, z0: float, tol: float) -> list[np.ndarray]:
    distances = tri[:, 2] - z0
    if np.all(distances >= -tol):
        clipped = tri.copy()
        clipped[np.abs(clipped[:, 2] - z0) <= tol, 2] = z0
        return [clipped]
    if np.all(distances < -tol):
        return []

    polygon = [tri[0], tri[1], tri[2]]
    output = []
    for i in range(3):
        current = polygon[i]
        prev = polygon[i - 1]
        current_inside = current[2] >= z0 - tol
        prev_inside = prev[2] >= z0 - tol

        if current_inside != prev_inside:
            dz = current[2] - prev[2]
            if abs(dz) > tol:
                t = (z0 - prev[2]) / dz
                point = prev + t * (current - prev)
                point[2] = z0
                output.append(point)
        if current_inside:
            point = current.copy()
            if abs(point[2] - z0) <= tol:
                point[2] = z0
            output.append(point)

    if len(output) < 3:
        return []

    clipped = np.asarray(output, dtype=float)
    result = []
    for i in range(1, len(clipped) - 1):
        tri_out = np.vstack([clipped[0], clipped[i], clipped[i + 1]])
        if np.linalg.norm(np.cross(tri_out[1] - tri_out[0], tri_out[2] - tri_out[0])) > 1.0e-12:
            result.append(tri_out)
    return result


def _build_mesh_from_triangles(triangles: list[np.ndarray]) -> "trimesh.Trimesh":
    if not triangles:
        return trimesh.Trimesh(
            vertices=np.empty((0, 3), dtype=float),
            faces=np.empty((0, 3), dtype=int),
            process=False,
        )

    vertex_map: dict[tuple[float, float, float], int] = {}
    vertices: list[list[float]] = []
    faces: list[list[int]] = []

    def add_vertex(pt: np.ndarray) -> int:
        key = tuple(float(x) for x in pt)
        idx = vertex_map.get(key)
        if idx is None:
            idx = len(vertices)
            vertex_map[key] = idx
            vertices.append([float(pt[0]), float(pt[1]), float(pt[2])])
        return idx

    for tri in triangles:
        faces.append([add_vertex(tri[0]), add_vertex(tri[1]), add_vertex(tri[2])])

    mesh = trimesh.Trimesh(vertices=np.asarray(vertices, dtype=float), faces=np.asarray(faces, dtype=int), process=False)
    mesh.remove_unreferenced_vertices()
    return mesh


def _iter_polygon_pieces(geometry, area_tol: float = 1.0e-12):
    if geometry is None or geometry.is_empty:
        return
    if isinstance(geometry, Polygon):
        if geometry.area > area_tol:
            yield geometry
        return
    if isinstance(geometry, (MultiPolygon, GeometryCollection)):
        for geom in geometry.geoms:
            yield from _iter_polygon_pieces(geom, area_tol=area_tol)


def _clip_ground_against_footprint(
    ground_mesh: "trimesh.Trimesh",
    footprint_union,
    *,
    z0: float = 0.0,
    tol: float = 1.0e-12,
) -> Tuple[list[np.ndarray], int, int]:
    if len(ground_mesh.faces) == 0:
        return [], 0, 0

    triangles_out: list[np.ndarray] = []
    preserved_faces = 0
    retriangulated_faces = 0
    vertices = np.asarray(ground_mesh.vertices, dtype=float)

    for face in np.asarray(ground_mesh.faces, dtype=int):
        tri2d = vertices[face][:, :2]
        tri_poly = Polygon(tri2d)
        if tri_poly.area <= tol:
            continue

        if footprint_union is None or footprint_union.is_empty:
            diff = tri_poly
        else:
            diff = tri_poly.difference(footprint_union)

        pieces = list(_iter_polygon_pieces(diff, area_tol=tol))
        if not pieces:
            continue

        unchanged = (
            len(pieces) == 1
            and abs(pieces[0].area - tri_poly.area) <= 1.0e-12
            and pieces[0].symmetric_difference(tri_poly).area <= 1.0e-12
        )
        if unchanged:
            preserved_faces += 1
        else:
            retriangulated_faces += 1

        for piece in pieces:
            for tri in triangulate(piece):
                if tri.area <= tol:
                    continue
                if not piece.covers(tri.representative_point()):
                    continue
                coords = np.asarray(tri.exterior.coords[:-1], dtype=float)
                if len(coords) != 3:
                    continue
                tri3d = np.column_stack([coords, np.full(3, float(z0), dtype=float)])
                if np.linalg.norm(np.cross(tri3d[1] - tri3d[0], tri3d[2] - tri3d[0])) <= 1.0e-12:
                    continue
                triangles_out.append(tri3d)

    return triangles_out, int(preserved_faces), int(retriangulated_faces)


def truncate_below_ground(
    mesh_or_geom,
    *,
    ground_planarity_tolerance: float = 1.0e-6,
    clip_tolerance: float = 1.0e-9,
    edgelength: float | None = None,
    return_trimesh: bool = False,
) -> Tuple[Union["trimesh.Trimesh", "UDGeom"], Dict[str, Union[int, float, bool]]]:
    """
    Remove geometry below an identified planar ground and restitch to the existing ground mesh.

    This routine is intended for the common CFD import problem where a building
    extends below an otherwise valid planar ground surface. The repair strategy
    is:
    1. identify the planar ground surface
    2. clip all non-ground faces at that elevation
    3. preserve as much of the original ground mesh as possible
    4. retriangulate only the local ground patch cut by the building footprint
    5. weld the resulting ground-building seam

    Parameters
    ----------
    mesh_or_geom : trimesh.Trimesh or UDGeom
        Mesh to clean.
    ground_planarity_tolerance : float, default=1e-6
        Maximum allowed deviation from a single ground plane.
    clip_tolerance : float, default=1e-9
        Numerical tolerance for triangle-plane clipping.
    edgelength : float, optional
        Optional reference spacing stored in the report. The current
        implementation preserves and locally restitches the existing ground
        mesh rather than remeshing the whole domain.
    return_trimesh : bool, default=False
        If True, return the repaired ``trimesh.Trimesh`` directly. Otherwise
        return a ``UDGeom`` object.

    Returns
    -------
    cleaned_geom : UDGeom or trimesh.Trimesh
        Geometry with below-ground geometry removed and welded back to the
        clipped ground. By default this is returned as ``UDGeom``.
    report : dict
        Summary of the truncation and local ground restitching operation.
    """
    if not TRIMESH_AVAILABLE:
        raise ImportError("trimesh is required. Install with: pip install trimesh")

    mesh = _copy_mesh(_as_trimesh(mesh_or_geom))
    if len(mesh.faces) == 0:
        report = {
            "ground_face_count": 0,
            "ground_z": 0.0,
            "removed_below_ground_faces": 0,
            "clipped_faces": 0,
            "edgelength": float(edgelength) if edgelength is not None else 0.0,
        }
        if return_trimesh:
            return mesh, report
        from .udgeom import UDGeom, DEFAULT_BACKEND
        return UDGeom(stl=mesh, backend=getattr(mesh_or_geom, "backend", DEFAULT_BACKEND)), report

    ground_mask = identify_ground_faces(mesh)
    ground_face_count = int(np.count_nonzero(ground_mask))
    z0 = _estimate_planar_ground_level(mesh, ground_mask, ground_planarity_tolerance)

    vertices = np.asarray(mesh.vertices, dtype=float).copy()
    vertices[:, 2] -= z0
    shifted = trimesh.Trimesh(vertices=vertices, faces=np.asarray(mesh.faces, dtype=int).copy(), process=False)

    building_mesh, _ = delete_ground(shifted)
    ground_faces = np.asarray(shifted.faces, dtype=int)[ground_mask]
    ground_mesh = trimesh.Trimesh(vertices=np.asarray(shifted.vertices, dtype=float).copy(), faces=ground_faces.copy(), process=False)
    spacing = _estimate_ground_spacing(shifted, ground_mask)

    triangles_out: list[np.ndarray] = []
    removed_below_ground_faces = 0
    clipped_faces = 0
    for face in np.asarray(building_mesh.faces, dtype=int):
        tri = np.asarray(building_mesh.vertices, dtype=float)[face]
        clipped = _clip_triangle_above_z(tri, 0.0, clip_tolerance)
        if not clipped:
            removed_below_ground_faces += 1
            continue
        if len(clipped) != 1 or not np.allclose(clipped[0], tri):
            clipped_faces += 1
        triangles_out.extend(clipped)

    clipped_buildings = _build_mesh_from_triangles(triangles_out)

    if len(clipped_buildings.faces) > 0:
        from .geometry_generation import _ground_footprint_union

        footprint_union = _ground_footprint_union(clipped_buildings)
    else:
        footprint_union = None

    clipped_ground_triangles, preserved_ground_faces, retriangulated_ground_faces = _clip_ground_against_footprint(
        ground_mesh,
        footprint_union,
        z0=0.0,
    )
    combined_triangles = clipped_ground_triangles + triangles_out
    cleaned = _build_mesh_from_triangles(combined_triangles)
    cleaned_vertices = np.asarray(cleaned.vertices, dtype=float)
    cleaned_vertices[:, 2] += z0
    cleaned = trimesh.Trimesh(vertices=cleaned_vertices, faces=np.asarray(cleaned.faces, dtype=int).copy(), process=False)
    cleaned.merge_vertices(digits_vertex=8)
    cleaned.remove_unreferenced_vertices()
    from .fix_mesh import weld_touching_boundaries
    cleaned, weld_report = weld_touching_boundaries(cleaned)

    report = {
        "ground_face_count": ground_face_count,
        "ground_z": float(z0),
        "removed_below_ground_faces": int(removed_below_ground_faces),
        "clipped_faces": int(clipped_faces),
        "edgelength": float(spacing),
        "preserved_ground_faces": int(preserved_ground_faces),
        "retriangulated_ground_faces": int(retriangulated_ground_faces),
        "welded_touching_boundary_faces": int(weld_report["welded_touching_boundary_faces"]),
    }
    if return_trimesh:
        return cleaned, report

    from .udgeom import UDGeom, DEFAULT_BACKEND

    return UDGeom(stl=cleaned, backend=getattr(mesh_or_geom, "backend", DEFAULT_BACKEND)), report


__all__ = ["truncate_below_ground"]
