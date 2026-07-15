"""
fix_mesh - Conservative mesh repair helpers for udgeom.

This module intentionally limits itself to unambiguous, auditable repairs:
- remove duplicate faces
- remove degenerate / zero-area faces
- remove unused vertices
- optionally remove tiny disconnected components
- optionally merge nearby vertices
"""

from __future__ import annotations

from typing import Dict, Tuple, Union

import numpy as np

try:
    import trimesh
    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False

from exceptions import DependencyError
from .check_mesh import (
    _as_trimesh,
    _boundary_edge_face_map,
    _find_unstitched_touching_edge_groups,
    check,
    identify_ground_faces,
)
from .delete_ground import delete_ground
from ._meshgraph import build_face_adjacency, connected_components

from shapely.geometry import GeometryCollection, MultiPolygon, Polygon
from shapely.ops import triangulate, unary_union
import triangle as triangle_lib


def _copy_mesh(mesh: "trimesh.Trimesh") -> "trimesh.Trimesh":
    return trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices, dtype=float).copy(),
        faces=np.asarray(mesh.faces, dtype=int).copy(),
        process=False,
    )


def _face_components(mesh: "trimesh.Trimesh") -> list[np.ndarray]:
    if len(mesh.faces) == 0:
        return []
    return connected_components(build_face_adjacency(mesh))


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


def _triangle_to_3d(coords2d: np.ndarray, axis: int, plane_value: float, sign: int) -> np.ndarray:
    if axis == 0:
        tri = np.column_stack([np.full(3, plane_value, dtype=float), coords2d[:, 0], coords2d[:, 1]])
    else:
        tri = np.column_stack([coords2d[:, 0], np.full(3, plane_value, dtype=float), coords2d[:, 1]])

    normal = np.cross(tri[1] - tri[0], tri[2] - tri[0])
    comp = normal[axis]
    if comp * sign < 0.0:
        tri = tri[[0, 2, 1]]
    return tri


def _triangulate_polygon_region(poly: "Polygon") -> list[np.ndarray]:
    if poly.is_empty or poly.area <= 1.0e-12:
        return []

    vertices = []
    segments = []
    holes = []

    def add_ring(coords):
        start = len(vertices)
        ring = np.asarray(coords[:-1], dtype=float)
        vertices.extend(ring.tolist())
        n = len(ring)
        for i in range(n):
            segments.append([start + i, start + ((i + 1) % n)])

    add_ring(poly.exterior.coords)
    for interior in poly.interiors:
        ring_poly = Polygon(interior)
        if ring_poly.area <= 1.0e-12:
            continue
        add_ring(interior.coords)
        rep = ring_poly.representative_point()
        holes.append([float(rep.x), float(rep.y)])

    pslg = {
        "vertices": np.asarray(vertices, dtype=float),
        "segments": np.asarray(segments, dtype=int),
    }
    if holes:
        pslg["holes"] = np.asarray(holes, dtype=float)

    result = triangle_lib.triangulate(pslg, "pQ")
    tri_vertices = np.asarray(result.get("vertices", []), dtype=float)
    tri_faces = np.asarray(result.get("triangles", []), dtype=int)
    triangles_out = []
    for simplex in tri_faces:
        tri2d = tri_vertices[np.asarray(simplex, dtype=int)]
        tri_poly = Polygon(tri2d)
        if tri_poly.area > 1.0e-12 and poly.covers(tri_poly.representative_point()):
            triangles_out.append(tri2d)
    if triangles_out:
        return triangles_out

    triangles_out = []
    for tri2d in triangulate(poly):
        rep = tri2d.representative_point()
        if not poly.covers(rep):
            continue
        coords2d = np.asarray(tri2d.exterior.coords[:3], dtype=float)
        triangles_out.append(coords2d)
    return triangles_out


def _resolve_axis_aligned_vertical_overlaps(
    mesh: "trimesh.Trimesh",
    plane_tolerance: float = 1.0e-8,
) -> tuple["trimesh.Trimesh", int]:
    """
    Remove overlapping opposite-facing coplanar vertical wall patches.

    This targets the common urban case of adjacent buildings with different
    heights, where one building wall partially overlaps the neighboring wall.
    """
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)
    normals = np.asarray(mesh.face_normals, dtype=float)

    plane_groups: dict[tuple[int, int], dict[int, list[int]]] = {}
    for face_id, face in enumerate(faces):
        tri = vertices[face]
        span_x = float(np.ptp(tri[:, 0]))
        span_y = float(np.ptp(tri[:, 1]))
        axis = None
        plane_value = None
        if span_x <= plane_tolerance and abs(normals[face_id, 0]) > 0.5:
            axis = 0
            plane_value = float(np.mean(tri[:, 0]))
        elif span_y <= plane_tolerance and abs(normals[face_id, 1]) > 0.5:
            axis = 1
            plane_value = float(np.mean(tri[:, 1]))
        if axis is None:
            continue
        sign = 1 if normals[face_id, axis] >= 0.0 else -1
        key = (axis, int(np.round(plane_value / plane_tolerance)))
        plane_groups.setdefault(key, {}).setdefault(sign, []).append(face_id)

    remove_face_ids: set[int] = set()
    new_triangles: list[np.ndarray] = []
    resolved_pairs = 0

    for (axis, plane_key), by_sign in plane_groups.items():
        if 1 not in by_sign or -1 not in by_sign:
            continue
        plane_value = float(plane_key) * plane_tolerance

        def build_union(face_ids: list[int]):
            polys = []
            for fid in face_ids:
                poly = _project_vertical_face(vertices[faces[fid]], axis)
                if poly is not None:
                    polys.append(poly)
            if not polys:
                return GeometryCollection()
            return unary_union(polys)

        pos_union = build_union(by_sign[1])
        neg_union = build_union(by_sign[-1])
        overlap = pos_union.intersection(neg_union)
        if overlap.is_empty or overlap.area <= 1.0e-8:
            continue

        resolved_pairs += 1
        pos_remaining = pos_union.difference(overlap)
        neg_remaining = neg_union.difference(overlap)
        remove_face_ids.update(by_sign[1])
        remove_face_ids.update(by_sign[-1])

        for sign, geom in ((1, pos_remaining), (-1, neg_remaining)):
            for poly in _iter_polygons(geom):
                for coords2d in _triangulate_polygon_region(poly):
                    tri3d = _triangle_to_3d(coords2d, axis, plane_value, sign)
                    new_triangles.append(tri3d)

    if not remove_face_ids:
        return mesh, 0

    keep_mask = np.ones(len(faces), dtype=bool)
    keep_mask[np.fromiter(sorted(remove_face_ids), dtype=int)] = False
    kept_faces = faces[keep_mask]

    vertex_map: dict[tuple[float, float, float], int] = {}
    new_vertices: list[list[float]] = []

    def add_vertex(pt: np.ndarray) -> int:
        key = tuple(float(x) for x in pt)
        idx = vertex_map.get(key)
        if idx is None:
            idx = len(new_vertices)
            vertex_map[key] = idx
            new_vertices.append([float(pt[0]), float(pt[1]), float(pt[2])])
        return idx

    rebuilt_faces: list[list[int]] = []
    for face in kept_faces:
        rebuilt_faces.append([add_vertex(vertices[int(face[0])]), add_vertex(vertices[int(face[1])]), add_vertex(vertices[int(face[2])])])
    for tri in new_triangles:
        rebuilt_faces.append([add_vertex(tri[0]), add_vertex(tri[1]), add_vertex(tri[2])])

    repaired = trimesh.Trimesh(
        vertices=np.asarray(new_vertices, dtype=float),
        faces=np.asarray(rebuilt_faces, dtype=int),
        process=False,
    )
    repaired.remove_unreferenced_vertices()
    return repaired, resolved_pairs


def _point_strictly_inside_segment(point: np.ndarray, a: np.ndarray, b: np.ndarray, tol: float) -> bool:
    ab = b - a
    ap = point - a
    denom = float(np.dot(ab, ab))
    if denom <= tol:
        return False
    cross = np.linalg.norm(np.cross(ab, ap))
    if cross > tol * max(1.0, np.linalg.norm(ab)):
        return False
    t = float(np.dot(ap, ab) / denom)
    return tol < t < 1.0 - tol


def _weld_touching_boundary_t_junctions(
    mesh: "trimesh.Trimesh",
    point_tolerance: float = 1.0e-8,
) -> tuple["trimesh.Trimesh", int]:
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)
    edge_to_faces = _boundary_edge_face_map(mesh)
    groups = _find_unstitched_touching_edge_groups(mesh, point_tolerance=point_tolerance)
    if not groups:
        return mesh, 0

    split_points: dict[tuple[int, int], list[np.ndarray]] = {}
    for group in groups:
        edges = [tuple(map(int, edge)) for edge in group["edges"]]
        for edge in edges:
            a = vertices[edge[0]]
            b = vertices[edge[1]]
            for other in edges:
                if other == edge:
                    continue
                for vid in other:
                    point = vertices[int(vid)]
                    if _point_strictly_inside_segment(point, a, b, point_tolerance):
                        split_points.setdefault(edge, []).append(point.copy())

    if not split_points:
        return mesh, 0

    vertex_map: dict[tuple[float, float, float], int] = {}
    new_vertices: list[list[float]] = []

    def add_vertex(pt: np.ndarray) -> int:
        key = tuple(round(float(x), 12) for x in pt)
        idx = vertex_map.get(key)
        if idx is None:
            idx = len(new_vertices)
            vertex_map[key] = idx
            new_vertices.append([float(pt[0]), float(pt[1]), float(pt[2])])
        return idx

    rebuilt_faces: list[list[int]] = []
    split_face_count = 0
    face_edges = (
        lambda face: [
            (tuple(sorted((int(face[0]), int(face[1])))), (int(face[0]), int(face[1])), int(face[2])),
            (tuple(sorted((int(face[1]), int(face[2])))), (int(face[1]), int(face[2])), int(face[0])),
            (tuple(sorted((int(face[2]), int(face[0])))), (int(face[2]), int(face[0])), int(face[1])),
        ]
    )

    for face in faces:
        applicable = [(undirected, oriented, third) for undirected, oriented, third in face_edges(face) if undirected in split_points]
        if len(applicable) != 1:
            rebuilt_faces.append([add_vertex(vertices[int(face[0])]), add_vertex(vertices[int(face[1])]), add_vertex(vertices[int(face[2])])])
            continue

        undirected, oriented, third = applicable[0]
        start = vertices[oriented[0]]
        end = vertices[oriented[1]]
        edge_vec = end - start
        denom = float(np.dot(edge_vec, edge_vec))
        pts = []
        for pt in split_points[undirected]:
            t = float(np.dot(pt - start, edge_vec) / denom) if denom > 0.0 else 0.0
            if point_tolerance < t < 1.0 - point_tolerance:
                pts.append((t, pt))
        if not pts:
            rebuilt_faces.append([add_vertex(vertices[int(face[0])]), add_vertex(vertices[int(face[1])]), add_vertex(vertices[int(face[2])])])
            continue

        split_face_count += 1
        pts.sort(key=lambda item: item[0])
        unique_pts = []
        for t, pt in pts:
            if unique_pts and abs(t - unique_pts[-1][0]) <= point_tolerance:
                continue
            unique_pts.append((t, pt))
        ordered = [vertices[oriented[0]]] + [pt for _, pt in unique_pts] + [vertices[oriented[1]]]
        third_pt = vertices[third]
        original_normal = np.cross(vertices[int(face[1])] - vertices[int(face[0])], vertices[int(face[2])] - vertices[int(face[0])])
        for p0, p1 in zip(ordered[:-1], ordered[1:]):
            tri = np.vstack([p0, p1, third_pt])
            normal = np.cross(tri[1] - tri[0], tri[2] - tri[0])
            if float(np.dot(normal, original_normal)) < 0.0:
                tri = tri[[0, 2, 1]]
            rebuilt_faces.append([add_vertex(tri[0]), add_vertex(tri[1]), add_vertex(tri[2])])

    repaired = trimesh.Trimesh(
        vertices=np.asarray(new_vertices, dtype=float),
        faces=np.asarray(rebuilt_faces, dtype=int),
        process=False,
    )
    repaired.merge_vertices(digits_vertex=8)
    repaired.remove_unreferenced_vertices()
    return repaired, split_face_count


def fix(
    mesh_or_geom,
    *,
    merge_tolerance: float | None = None,
    remove_small_components: bool = False,
    min_component_faces: int = 0,
    min_component_area_fraction: float = 0.0,
    resolve_vertical_coplanar_overlaps: bool = False,
    weld_touching_boundaries: bool = False,
) -> Tuple["trimesh.Trimesh", Dict[str, Union[int, float, bool, list, dict]]]:
    """
    Apply conservative small-scale mesh cleanup and return the repaired mesh plus a report.

    This routine is intended for low-risk housekeeping operations that do not
    change the geometric meaning of the mesh. Typical uses are:
    - removing duplicate faces
    - removing degenerate or zero-area triangles
    - removing unused vertices
    - optionally removing tiny disconnected fragments
    - optionally merging nearly coincident vertices

    For larger geometric repairs, such as resolving internal shared walls or
    welding unstitched touching seams, use the dedicated helpers
    ``resolve_vertical_coplanar_overlaps(...)`` and
    ``weld_touching_boundaries(...)``.

    Parameters
    ----------
    mesh_or_geom : trimesh.Trimesh or UDGeom
        Mesh to repair.
    merge_tolerance : float, optional
        If provided and > 0, merge vertices after rounding coordinates to this tolerance.
    remove_small_components : bool, default=False
        If True, remove disconnected components below the requested thresholds.
    min_component_faces : int, default=0
        Minimum face count for components to keep when `remove_small_components=True`.
    min_component_area_fraction : float, default=0.0
        Minimum component area as a fraction of total mesh area to keep.
    resolve_vertical_coplanar_overlaps : bool, default=False
        If True, also resolve opposite-facing coplanar vertical overlaps.
        This option is supported for convenience, but the dedicated
        ``resolve_vertical_coplanar_overlaps(...)`` routine is usually clearer.
    weld_touching_boundaries : bool, default=False
        If True, also weld touching but unstitched boundary seams. This option
        is supported for convenience, but the dedicated
        ``weld_touching_boundaries(...)`` routine is usually clearer.

    Returns
    -------
    fixed_mesh : trimesh.Trimesh
        Conservatively cleaned mesh.
    report : dict
        Summary of what was changed, including ``before`` and ``after`` mesh
        checks.
    """
    if not TRIMESH_AVAILABLE:
        raise DependencyError("trimesh is required. Install with: pip install trimesh")

    original = _as_trimesh(mesh_or_geom)
    mesh = _copy_mesh(original)

    report: Dict[str, Union[int, float, bool, list, dict]] = {
        "removed_duplicate_faces": 0,
        "removed_degenerate_faces": 0,
        "removed_zero_area_faces": 0,
        "removed_small_component_faces": 0,
        "removed_unused_vertices": 0,
        "merged_vertices": 0,
        "resolved_vertical_overlap_regions": 0,
        "welded_touching_boundary_faces": 0,
    }

    before = check(mesh, require_single_component=False)
    report["before"] = before

    # Remove duplicate faces, keeping first occurrence.
    faces = np.asarray(mesh.faces, dtype=int)
    canonical = np.sort(faces, axis=1)
    _, keep_idx = np.unique(canonical, axis=0, return_index=True)
    keep_mask = np.zeros(len(faces), dtype=bool)
    keep_mask[np.sort(keep_idx)] = True
    report["removed_duplicate_faces"] = int(np.count_nonzero(~keep_mask))
    if not np.all(keep_mask):
        mesh.update_faces(keep_mask)
        mesh.remove_unreferenced_vertices()

    # Remove degenerate / zero-area faces.
    faces = np.asarray(mesh.faces, dtype=int)
    degenerate_mask = np.array([len(np.unique(face)) < 3 for face in faces], dtype=bool)
    zero_area_mask = np.asarray(mesh.area_faces, dtype=float) <= 1.0e-12
    drop_mask = degenerate_mask | zero_area_mask
    report["removed_degenerate_faces"] = int(np.count_nonzero(degenerate_mask))
    report["removed_zero_area_faces"] = int(np.count_nonzero(zero_area_mask & ~degenerate_mask))
    if np.any(drop_mask):
        mesh.update_faces(~drop_mask)
        mesh.remove_unreferenced_vertices()

    # Optional tiny disconnected component removal.
    if remove_small_components and len(mesh.faces) > 0:
        components = _face_components(mesh)
        if len(components) > 1:
            face_areas = np.asarray(mesh.area_faces, dtype=float)
            total_area = float(np.sum(face_areas))
            keep_mask = np.ones(len(mesh.faces), dtype=bool)
            removed_faces = 0
            for comp in components:
                comp_area = float(np.sum(face_areas[comp]))
                keep = True
                if min_component_faces > 0 and len(comp) < min_component_faces:
                    keep = False
                if keep and min_component_area_fraction > 0.0 and total_area > 0.0:
                    if comp_area / total_area < min_component_area_fraction:
                        keep = False
                if not keep:
                    keep_mask[comp] = False
                    removed_faces += int(len(comp))
            if removed_faces > 0 and np.any(keep_mask):
                mesh.update_faces(keep_mask)
                mesh.remove_unreferenced_vertices()
                report["removed_small_component_faces"] = int(removed_faces)

    if resolve_vertical_coplanar_overlaps and len(mesh.faces) > 0:
        mesh, n_regions = _resolve_axis_aligned_vertical_overlaps(mesh)
        report["resolved_vertical_overlap_regions"] = int(n_regions)

    if (weld_touching_boundaries or resolve_vertical_coplanar_overlaps) and len(mesh.faces) > 0:
        mesh, n_faces = _weld_touching_boundary_t_junctions(mesh)
        report["welded_touching_boundary_faces"] = int(n_faces)

    # Optional vertex merge by tolerance.
    if merge_tolerance is not None and merge_tolerance > 0.0 and len(mesh.vertices) > 0:
        before_vertices = int(len(mesh.vertices))
        digits = max(0, int(np.ceil(-np.log10(float(merge_tolerance)))))
        mesh.merge_vertices(digits_vertex=digits)
        mesh.remove_unreferenced_vertices()
        report["merged_vertices"] = int(before_vertices - len(mesh.vertices))

    # Final cleanup
    before_vertices = int(len(mesh.vertices))
    mesh.remove_unreferenced_vertices()
    report["removed_unused_vertices"] = int(before_vertices - len(mesh.vertices))

    after = check(mesh, require_single_component=False)
    report["after"] = after
    return mesh, report


def resolve_vertical_coplanar_overlaps(
    mesh_or_geom,
    *,
    weld_touching_boundaries: bool = True,
    plane_tolerance: float = 1.0e-8,
    point_tolerance: float = 1.0e-8,
) -> Tuple["trimesh.Trimesh", Dict[str, Union[int, dict]]]:
    """
    Resolve opposite-facing coplanar vertical wall overlaps.

    This targets the common adjacent-building case where two shells both keep
    faces on the same shared vertical plane, producing an internal touching
    wall rather than one clean exterior surface.

    The routine removes the overlapping internal wall patch, retriangulates the
    remaining exposed wall region, and can optionally weld any touching
    boundary T-junctions created by that repair.
    """
    if not TRIMESH_AVAILABLE:
        raise DependencyError("trimesh is required. Install with: pip install trimesh")

    mesh = _copy_mesh(_as_trimesh(mesh_or_geom))
    before = check(mesh, require_single_component=False)
    mesh, n_regions = _resolve_axis_aligned_vertical_overlaps(mesh, plane_tolerance=plane_tolerance)
    welded_faces = 0
    if weld_touching_boundaries and len(mesh.faces) > 0:
        mesh, welded_faces = _weld_touching_boundary_t_junctions(mesh, point_tolerance=point_tolerance)
    after = check(mesh, require_single_component=False)
    report = {
        "resolved_vertical_overlap_regions": int(n_regions),
        "welded_touching_boundary_faces": int(welded_faces),
        "before": before,
        "after": after,
    }
    return mesh, report


def weld_touching_boundaries(
    mesh_or_geom,
    *,
    point_tolerance: float = 1.0e-8,
) -> Tuple["trimesh.Trimesh", Dict[str, Union[int, dict]]]:
    """
    Weld touching boundary seams by splitting boundary T-junctions.

    This routine targets the case where surfaces touch geometrically along a
    line segment but do not yet share the same edge discretization. It repairs
    those seams by splitting longer boundary edges at touching intermediate
    vertices and rebuilding the affected triangles.
    """
    if not TRIMESH_AVAILABLE:
        raise DependencyError("trimesh is required. Install with: pip install trimesh")

    mesh = _copy_mesh(_as_trimesh(mesh_or_geom))
    before = check(mesh, require_single_component=False)
    mesh, n_faces = _weld_touching_boundary_t_junctions(mesh, point_tolerance=point_tolerance)
    after = check(mesh, require_single_component=False)
    report = {
        "welded_touching_boundary_faces": int(n_faces),
        "before": before,
        "after": after,
    }
    return mesh, report


def repair_adjacent_buildings(
    mesh_or_geom,
    *,
    merge_digits_vertex: int = 8,
    plane_tolerance: float = 1.0e-8,
    point_tolerance: float = 1.0e-8,
    return_trimesh: bool = False,
) -> Tuple[Union["trimesh.Trimesh", "UDGeom"], Dict[str, Union[int, dict]]]:
    """
    Repair adjacent-building overlap defects while preserving the original ground.

    This user-facing helper is intended for the common case where an imported
    geometry contains:
    - a valid ground surface
    - adjacent building shells
    - an internal shared vertical wall where neighboring buildings both keep
      faces on the same plane

    The routine keeps the original ground, repairs the non-ground shell with
    ``resolve_vertical_coplanar_overlaps(...)``, rejoins shell and ground, and
    welds the stitched result.

    Parameters
    ----------
    mesh_or_geom : trimesh.Trimesh or UDGeom
        Geometry to repair.
    merge_digits_vertex : int, default=8
        Decimal rounding precision used when welding coincident vertices after
        shell and ground are rejoined.
    plane_tolerance : float, default=1e-8
        Tolerance used to detect axis-aligned coplanar overlap regions.
    point_tolerance : float, default=1e-8
        Tolerance used to weld touching boundary seams after rejoining.
    return_trimesh : bool, default=False
        If True, return the repaired ``trimesh.Trimesh`` directly. Otherwise
        return a ``UDGeom`` object.

    Returns
    -------
    cleaned_geom, report
        Repaired full geometry as ``UDGeom`` by default, plus a compact repair
        summary.
    """
    if not TRIMESH_AVAILABLE:
        raise DependencyError("trimesh is required. Install with: pip install trimesh")

    mesh = _copy_mesh(_as_trimesh(mesh_or_geom))
    ground_mask = identify_ground_faces(mesh)
    ground_face_count = int(np.count_nonzero(ground_mask))
    if ground_face_count == 0:
        raise ValueError("No ground surface could be identified in the input mesh")

    ground = trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices, dtype=float).copy(),
        faces=np.asarray(mesh.faces, dtype=int)[ground_mask].copy(),
        process=False,
    )
    shell, _ = delete_ground(mesh)
    shell, shell_report = resolve_vertical_coplanar_overlaps(
        shell,
        weld_touching_boundaries=False,
        plane_tolerance=plane_tolerance,
        point_tolerance=point_tolerance,
    )

    cleaned = trimesh.util.concatenate([ground, shell])
    cleaned.merge_vertices(digits_vertex=int(merge_digits_vertex))
    cleaned.remove_unreferenced_vertices()
    cleaned, welded_face_count = _weld_touching_boundary_t_junctions(
        cleaned,
        point_tolerance=point_tolerance,
    )
    report = {
        "ground_face_count": ground_face_count,
        "resolved_vertical_overlap_regions": int(shell_report["resolved_vertical_overlap_regions"]),
        "welded_touching_boundary_faces": int(welded_face_count),
    }
    if return_trimesh:
        return cleaned, report

    from .udgeom import UDGeom, DEFAULT_BACKEND

    return UDGeom(stl=cleaned, backend=getattr(mesh_or_geom, "backend", DEFAULT_BACKEND)), report


__all__ = ["fix", "resolve_vertical_coplanar_overlaps", "weld_touching_boundaries", "repair_adjacent_buildings"]
