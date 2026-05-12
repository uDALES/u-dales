"""
Geometry generation helpers ported from MATLAB udgeom.

These functions mirror the MATLAB implementations in tools/matlab/+udgeom:
- createFlatSurface.m
- createCanyons.m
- createCubes.m
- createRealistic.m
"""

from pathlib import Path
from functools import reduce
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union
import warnings

import numpy as np

try:
    import trimesh
except ImportError as exc:
    raise ImportError("trimesh is required for geometry generation; install with `pip install trimesh`.") from exc

try:
    from scipy.spatial import Delaunay
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    from shapely.geometry import LineString, Polygon
    from shapely.ops import polygonize, triangulate, unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False

from .udgeom import UDGeom


# -----------------------------------------------------------------------------#
# Utility helpers
# -----------------------------------------------------------------------------#

def _remove_unused(vertices: np.ndarray, faces: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Remove vertices that are not referenced by any face."""
    used = np.unique(faces.flatten())
    new_index = -np.ones(len(vertices), dtype=int)
    new_index[used] = np.arange(len(used))
    new_vertices = vertices[used]
    new_faces = new_index[faces]
    return new_vertices, new_faces


def _divide_faces(mesh: "trimesh.Trimesh", times: int = 1) -> "trimesh.Trimesh":
    """
    Subdivide each triangular face into four by adding midpoints.
    Repeats `times` times (MATLAB divideFaces equivalent).
    """
    if times <= 0:
        return mesh.copy()

    current = mesh
    for _ in range(times):
        verts = current.vertices
        faces = current.faces

        edge_mid_cache = {}

        def midpoint(a: int, b: int) -> int:
            key = tuple(sorted((a, b)))
            if key in edge_mid_cache:
                return edge_mid_cache[key]
            mid = (verts[a] + verts[b]) / 2.0
            edge_mid_cache[key] = len(verts_list)
            verts_list.append(mid)
            return edge_mid_cache[key]

        verts_list = verts.tolist()
        new_faces = []
        for f in faces:
            a, b, c = f
            ab = midpoint(a, b)
            bc = midpoint(b, c)
            ca = midpoint(c, a)
            new_faces.extend(
                [
                    (a, ab, ca),
                    (ab, b, bc),
                    (ca, bc, c),
                    (ab, bc, ca),
                ]
            )

        verts_arr = np.asarray(verts_list, dtype=float)
        faces_arr = np.asarray(new_faces, dtype=int)
        verts_arr, faces_arr = _remove_unused(verts_arr, faces_arr)
        current = trimesh.Trimesh(vertices=verts_arr, faces=faces_arr, process=False)

    return current


def _unit_cube() -> "trimesh.Trimesh":
    """
    Unit cube centered at origin with side length 1 (open bottom to match MATLAB).
    Vertices are in [-0.5, 0.5] for each axis.
    """
    verts = np.array([
        [-0.5, -0.5, -0.5],
        [-0.5,  0.5, -0.5],
        [ 0.5,  0.5, -0.5],
        [ 0.5, -0.5, -0.5],
        [ 0.5, -0.5,  0.5],
        [-0.5, -0.5,  0.5],
        [-0.5,  0.5,  0.5],
        [ 0.5,  0.5,  0.5],
    ], dtype=float)

    faces = np.array([
        [0, 1, 2], [2, 3, 0],  # bottom (open in MATLAB; keep for watertightness)
        [0, 3, 4], [4, 5, 0],  # side
        [0, 5, 6], [6, 1, 0],  # side
        [7, 4, 3], [3, 2, 7],  # side
        [7, 2, 1], [1, 6, 7],  # side
        [7, 6, 5], [5, 4, 7],  # top
    ], dtype=int)

    return trimesh.Trimesh(vertices=verts, faces=faces, process=False)


def _operate_unit_cube(scale: Sequence[float], shifts: np.ndarray, divisions: int) -> "trimesh.Trimesh":
    """
    Scale and shift subdivided unit cube copies, then concatenate.
    """
    base = _divide_faces(_unit_cube(), times=divisions)
    scaled_vertices = base.vertices * np.asarray(scale, dtype=float)
    base_faces = base.faces

    meshes = []
    for shift in shifts:
        verts = scaled_vertices + np.asarray(shift, dtype=float)
        meshes.append(trimesh.Trimesh(vertices=verts, faces=base_faces, process=False))

    return trimesh.util.concatenate(meshes)


def _structured_ground_points(xsize: float, ysize: float, edgelength: float) -> np.ndarray:
    """Create regularly spaced ground points including domain edges."""
    nx = int(round(xsize / edgelength))
    ny = int(round(ysize / edgelength))
    xs = np.linspace(0.0, xsize, nx + 1)
    ys = np.linspace(0.0, ysize, ny + 1)
    pts = np.array([[x, y] for x in xs for y in ys], dtype=float)
    return pts


def _grid_axis(length: float, edgelength: float) -> np.ndarray:
    """Return grid coordinates including both domain endpoints."""
    n = int(round(length / edgelength))
    return np.linspace(0.0, length, n + 1)


def _build_ground_edge_lines(existing: "trimesh.Trimesh") -> list:
    """Extract z=0 mesh edges as 2D line constraints."""
    verts = np.asarray(existing.vertices, dtype=float)
    edges = np.asarray(existing.edges_unique, dtype=int)
    if len(edges) == 0:
        return []

    mask = np.isclose(verts[edges][:, :, 2], 0.0).all(axis=1)
    ground_edges = edges[mask]

    lines = []
    for edge in ground_edges:
        p0 = tuple(verts[edge[0], :2])
        p1 = tuple(verts[edge[1], :2])
        if np.linalg.norm(np.asarray(p0) - np.asarray(p1)) > 1e-12:
            lines.append(LineString([p0, p1]))
    return lines


def _ground_footprint_union(existing: "trimesh.Trimesh"):
    """Construct 2D footprint polygons from z=0 building edges."""
    if not SHAPELY_AVAILABLE:
        return None

    edge_lines = _build_ground_edge_lines(existing)
    if not edge_lines:
        return None

    polygons = list(polygonize(unary_union(edge_lines)))
    if not polygons:
        return None

    filtered = []
    for i, poly in enumerate(polygons):
        is_duplicate_hole = False
        for j, other in enumerate(polygons):
            if i == j or other.area <= poly.area:
                continue
            other_shell = Polygon(other.exterior.coords)
            if other_shell.covers(poly) and other.intersection(poly).area <= 1e-12:
                is_duplicate_hole = True
                break
        if not is_duplicate_hole:
            filtered.append(poly)

    return unary_union(filtered)


def _triangles_to_mesh(triangles: list) -> "trimesh.Trimesh":
    """Convert a list of shapely triangles into a trimesh surface on z=0."""
    vertices = []
    faces = []
    vertex_map = {}

    for tri in triangles:
        if hasattr(tri, "exterior"):
            coords = list(tri.exterior.coords)[:-1]
        else:
            coords = [tuple(row) for row in np.asarray(tri, dtype=float)]
        if len(coords) != 3:
            continue
        face = []
        for x, y in coords:
            key = (round(float(x), 12), round(float(y), 12), 0.0)
            if key not in vertex_map:
                vertex_map[key] = len(vertices)
                vertices.append(key)
            face.append(vertex_map[key])
        faces.append(face)

    if not faces:
        return trimesh.Trimesh(vertices=np.empty((0, 3)), faces=np.empty((0, 3), dtype=int), process=False)

    return trimesh.Trimesh(vertices=np.asarray(vertices, dtype=float), faces=np.asarray(faces, dtype=int), process=False)


def _signed_area_2d(points: np.ndarray) -> float:
    """Signed polygon area; positive for counter-clockwise rings."""
    x = points[:, 0]
    y = points[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _point_in_triangle(pt: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray, eps: float = 1e-12) -> bool:
    """Return True when pt lies inside or on the boundary of triangle abc."""
    v0 = c - a
    v1 = b - a
    v2 = pt - a

    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)
    denom = dot00 * dot11 - dot01 * dot01
    if abs(denom) <= eps:
        return False
    inv = 1.0 / denom
    u = (dot11 * dot02 - dot01 * dot12) * inv
    v = (dot00 * dot12 - dot01 * dot02) * inv
    return u >= -eps and v >= -eps and (u + v) <= 1.0 + eps


def _triangulate_simple_polygon(piece: "Polygon") -> List[np.ndarray]:
    """
    Triangulate a simple polygon while preserving every boundary edge exactly.

    The polygons produced by the constrained ground cell construction are simple
    rings without holes. Ear clipping is sufficient here and ensures the ground
    shares the same boundary edges as the building base.
    """
    if len(piece.interiors) > 0:
        return []

    coords = np.asarray(piece.exterior.coords[:-1], dtype=float)
    if len(coords) < 3:
        return []

    # Remove duplicate consecutive vertices and collinear slivers.
    cleaned = []
    for pt in coords:
        if not cleaned or np.linalg.norm(pt - cleaned[-1]) > 1e-12:
            cleaned.append(pt)
    coords = np.asarray(cleaned, dtype=float)
    if len(coords) < 3:
        return []

    if _signed_area_2d(coords) < 0.0:
        coords = coords[::-1]

    indices = list(range(len(coords)))
    triangles: List[np.ndarray] = []
    eps = 1e-12

    while len(indices) > 3:
        ear_found = False
        n = len(indices)
        for pos in range(n):
            ia = indices[(pos - 1) % n]
            ib = indices[pos]
            ic = indices[(pos + 1) % n]
            a, b, c = coords[ia], coords[ib], coords[ic]

            cross = np.cross(b - a, c - b)
            if cross <= eps:
                continue

            tri = np.array([a, b, c], dtype=float)
            contains_other = False
            for idx in indices:
                if idx in (ia, ib, ic):
                    continue
                if _point_in_triangle(coords[idx], a, b, c, eps=eps):
                    contains_other = True
                    break
            if contains_other:
                continue

            triangles.append(tri)
            indices.pop(pos)
            ear_found = True
            break

        if not ear_found:
            # Fallback for unexpected degeneracy; better to return a valid but
            # less exact triangulation than fail outright.
            return [np.asarray(list(t.exterior.coords)[:-1], dtype=float) for t in triangulate(piece) if piece.covers(t)]

    if len(indices) == 3:
        triangles.append(coords[np.asarray(indices, dtype=int)])
    return triangles


def _sample_linestring_points(coords: np.ndarray, spacing: float) -> np.ndarray:
    """Sample points along a polyline including both endpoints."""
    if len(coords) == 0:
        return np.empty((0, 2), dtype=float)

    sampled = [coords[0]]
    spacing = max(float(spacing), 1e-9)
    for start, end in zip(coords[:-1], coords[1:]):
        seg = end - start
        length = float(np.linalg.norm(seg))
        if length <= 1e-12:
            continue
        nseg = max(1, int(np.ceil(length / spacing)))
        for step in range(1, nseg + 1):
            sampled.append(start + (step / nseg) * seg)
    ordered = []
    for pt in np.round(np.asarray(sampled, dtype=float), 12):
        if not ordered or np.linalg.norm(pt - ordered[-1]) > 1e-12:
            ordered.append(pt)
    if len(ordered) > 1 and np.linalg.norm(ordered[0] - ordered[-1]) <= 1e-12:
        ordered = ordered[:-1]
    return np.asarray(ordered, dtype=float)


def _triangle_quality(tri: np.ndarray) -> float:
    """Return a simple triangle quality metric in [0, 1]."""
    a = float(np.linalg.norm(tri[1] - tri[0]))
    b = float(np.linalg.norm(tri[2] - tri[1]))
    c = float(np.linalg.norm(tri[0] - tri[2]))
    denom = a * a + b * b + c * c
    if denom <= 1e-12:
        return 0.0
    area = 0.5 * abs(float(np.cross(tri[1] - tri[0], tri[2] - tri[0])))
    return float(4.0 * np.sqrt(3.0) * area / denom)


def _refine_skewed_triangles(
    triangles: List[np.ndarray],
    target_spacing: float,
    min_quality: float = 0.18,
    max_passes: int = 2,
) -> List[np.ndarray]:
    """
    Improve very skewed triangles by midpoint insertion on the longest edge.

    This is intentionally conservative: it keeps the original exact cell
    triangulation and only adds points where a triangle is clearly poor and
    still large enough relative to the target spacing.
    """
    refined = [np.asarray(tri, dtype=float) for tri in triangles]
    target_spacing = max(float(target_spacing), 1e-6)

    for _ in range(max_passes):
        changed = False
        next_tris: List[np.ndarray] = []
        for tri in refined:
            edges = [
                (0, 1, float(np.linalg.norm(tri[1] - tri[0]))),
                (1, 2, float(np.linalg.norm(tri[2] - tri[1]))),
                (2, 0, float(np.linalg.norm(tri[0] - tri[2]))),
            ]
            longest = max(edges, key=lambda item: item[2])
            quality = _triangle_quality(tri)
            if quality >= min_quality or longest[2] <= 1.35 * target_spacing:
                next_tris.append(tri)
                continue

            i, j, _ = longest
            k = next(idx for idx in (0, 1, 2) if idx not in (i, j))
            midpoint = 0.5 * (tri[i] + tri[j])
            next_tris.append(np.asarray([tri[i], midpoint, tri[k]], dtype=float))
            next_tris.append(np.asarray([midpoint, tri[j], tri[k]], dtype=float))
            changed = True

        refined = next_tris
        if not changed:
            break

    return refined


def _triangulate_polygon_piece(piece: "Polygon", spacing: float) -> List[np.ndarray]:
    """
    Triangulate a polygon cell exactly, then locally refine poor triangles.

    This keeps the topology of the original working approach while adding
    midpoint insertion only where individual triangles are obviously skewed.
    """
    if piece.is_empty or piece.area <= 1e-12:
        return []
    if len(piece.interiors) > 0:
        base = _triangulate_simple_polygon(piece)
        return _refine_skewed_triangles(base, spacing) if base else []

    base = _triangulate_simple_polygon(piece)
    if not base:
        return []
    return _refine_skewed_triangles(base, spacing)


def _point_on_segment_2d(point: np.ndarray, a: np.ndarray, b: np.ndarray, tol: float = 1e-9) -> bool:
    """Return True when a 2D point lies on the closed segment a-b."""
    ab = b - a
    ap = point - a
    cross = ab[0] * ap[1] - ab[1] * ap[0]
    if abs(cross) > tol:
        return False
    dot = float(np.dot(ap, ab))
    if dot < -tol:
        return False
    if dot > float(np.dot(ab, ab)) + tol:
        return False
    return True


def _split_existing_ground_boundary(existing: "trimesh.Trimesh", ground_mesh: "trimesh.Trimesh") -> "trimesh.Trimesh":
    """
    Split wall faces along new ground-boundary vertices so ground and buildings share edges.

    When the ground generator introduces extra points on the z=0 building perimeter,
    the adjacent wall triangles must be subdivided as well; otherwise the meshes only
    touch geometrically and remain topologically disconnected.
    """
    verts = np.asarray(existing.vertices, dtype=float)
    faces = np.asarray(existing.faces, dtype=int)
    gverts = np.asarray(ground_mesh.vertices, dtype=float)

    if len(gverts) == 0:
        return existing.copy()

    new_vertices = verts.tolist()
    vertex_map = {tuple(np.round(v, 12)): i for i, v in enumerate(verts)}
    new_faces: List[Tuple[int, int, int]] = []

    def get_vid(pt: np.ndarray) -> int:
        key = tuple(np.round(np.asarray(pt, dtype=float), 12))
        if key in vertex_map:
            return vertex_map[key]
        idx = len(new_vertices)
        new_vertices.append([float(key[0]), float(key[1]), float(key[2])])
        vertex_map[key] = idx
        return idx

    for face in faces:
        tri = verts[face]
        ground_mask = np.isclose(tri[:, 2], 0.0)
        if ground_mask.sum() != 2:
            new_faces.append(tuple(face.tolist()))
            continue

        bottom_local = np.where(ground_mask)[0]
        top_local = int(np.where(~ground_mask)[0][0])
        a = tri[bottom_local[0], :2]
        b = tri[bottom_local[1], :2]
        edge_points = [tri[bottom_local[0]], tri[bottom_local[1]]]

        for gv in gverts[np.isclose(gverts[:, 2], 0.0)]:
            if _point_on_segment_2d(gv[:2], a, b):
                edge_points.append(gv)

        if len(edge_points) <= 2:
            new_faces.append(tuple(face.tolist()))
            continue

        ab = b - a
        denom = float(np.dot(ab, ab))
        if denom <= 1e-12:
            new_faces.append(tuple(face.tolist()))
            continue

        unique_points = np.unique(np.round(np.asarray(edge_points, dtype=float), 12), axis=0)
        tvals = [float(np.dot(pt[:2] - a, ab) / denom) for pt in unique_points]
        order = np.argsort(tvals)
        ordered = unique_points[order]

        top_vid = int(face[top_local])
        original_normal = np.cross(tri[1] - tri[0], tri[2] - tri[0])
        for p0, p1 in zip(ordered[:-1], ordered[1:]):
            v0 = get_vid(p0)
            v1 = get_vid(p1)
            candidate = np.array([new_vertices[v0], new_vertices[v1], new_vertices[top_vid]], dtype=float)
            normal = np.cross(candidate[1] - candidate[0], candidate[2] - candidate[0])
            if np.dot(normal, original_normal) < 0.0:
                new_faces.append((v1, v0, top_vid))
            else:
                new_faces.append((v0, v1, top_vid))

    verts_arr = np.asarray(new_vertices, dtype=float)
    faces_arr = np.asarray(new_faces, dtype=int)
    verts_arr, faces_arr = _remove_unused(verts_arr, faces_arr)
    return trimesh.Trimesh(vertices=verts_arr, faces=faces_arr, process=False)


def _generate_ground_with_constraints(
    existing: "trimesh.Trimesh", xsize: float, ysize: float, edgelength: float
) -> tuple[Optional["trimesh.Trimesh"], Dict[str, Any]]:
    """
    Build a constrained ground mesh using a polygonized grid plus building footprint edges.

    This follows the MATLAB `generateGround` intent more closely than an unconstrained
    Delaunay triangulation: building footprint edges remain hard constraints, so ground
    triangles cannot cut across building bases.
    """
    if not SHAPELY_AVAILABLE:
        return None, {}

    domain = Polygon([(0.0, 0.0), (xsize, 0.0), (xsize, ysize), (0.0, ysize)])
    xs = _grid_axis(xsize, edgelength)
    ys = _grid_axis(ysize, edgelength)

    footprint_union = _ground_footprint_union(existing)
    outside_region = domain.difference(footprint_union) if footprint_union is not None else domain

    grid_lines = [LineString([(float(x), 0.0), (float(x), ysize)]) for x in xs]
    grid_lines.extend(LineString([(0.0, float(y)), (xsize, float(y))]) for y in ys)
    clipped_grid_lines = []
    for line in grid_lines:
        clipped = line.intersection(outside_region)
        if clipped.is_empty:
            continue
        if clipped.geom_type == "LineString":
            clipped_grid_lines.append(clipped)
        else:
            clipped_grid_lines.extend([geom for geom in getattr(clipped, "geoms", []) if geom.geom_type == "LineString"])

    boundary_lines = [LineString(domain.exterior.coords)]
    if footprint_union is not None:
        if footprint_union.geom_type == "Polygon":
            footprint_geoms = [footprint_union]
        else:
            footprint_geoms = list(footprint_union.geoms)
        for geom in footprint_geoms:
            boundary_lines.append(LineString(geom.exterior.coords))
            for ring in geom.interiors:
                boundary_lines.append(LineString(ring.coords))

    lines = boundary_lines + clipped_grid_lines

    planar_cells = list(polygonize(unary_union(lines)))
    if not planar_cells:
        return None, {
            "domain": domain,
            "footprint_union": footprint_union,
            "outside_region": outside_region,
            "grid_lines": clipped_grid_lines,
            "planar_cells": [],
        }

    triangles = []
    for poly in planar_cells:
        if poly.is_empty or poly.area <= 1e-12:
            continue
        cell = poly.intersection(domain)
        if cell.is_empty:
            continue
        if footprint_union is not None:
            cell = cell.difference(footprint_union)
        if cell.is_empty:
            continue

        pieces = [cell] if cell.geom_type == "Polygon" else list(getattr(cell, "geoms", []))
        for piece in pieces:
            if piece.is_empty or piece.area <= 1e-12:
                continue
            tris_piece = _triangulate_polygon_piece(piece, edgelength)
            if tris_piece:
                triangles.extend(tris_piece)
            else:
                for tri in triangulate(piece):
                    if tri.area > 1e-12 and piece.covers(tri):
                        triangles.append(tri)

    ground_mesh = _triangles_to_mesh(triangles)
    if len(ground_mesh.faces) == 0:
        return None, {
            "domain": domain,
            "footprint_union": footprint_union,
            "outside_region": outside_region,
            "grid_lines": clipped_grid_lines,
            "planar_cells": planar_cells,
        }

    existing_split = _split_existing_ground_boundary(existing, ground_mesh)
    combined = trimesh.util.concatenate([existing_split, ground_mesh])
    combined.merge_vertices(digits_vertex=12)
    verts, faces = _remove_unused(np.asarray(combined.vertices), np.asarray(combined.faces))
    return trimesh.Trimesh(vertices=verts, faces=faces, process=False), {
        "domain": domain,
        "footprint_union": footprint_union,
        "outside_region": outside_region,
        "grid_lines": clipped_grid_lines,
        "planar_cells": planar_cells,
    }


def _triangulate_ground(points_xy: np.ndarray, xsize: float, ysize: float) -> "trimesh.Trimesh":
    """
    Triangulate 2D points on the ground plane.
    Uses SciPy Delaunay if available, otherwise falls back to a structured grid.
    """
    if SCIPY_AVAILABLE:
        tri = Delaunay(points_xy)
        simplices = tri.simplices
        centroids = points_xy[simplices].mean(axis=1)
        keep = (
            (centroids[:, 0] >= -1e-9)
            & (centroids[:, 0] <= xsize + 1e-9)
            & (centroids[:, 1] >= -1e-9)
            & (centroids[:, 1] <= ysize + 1e-9)
        )
        simplices = simplices[keep]
        vertices_3d = np.column_stack([points_xy, np.zeros(len(points_xy))])
        return trimesh.Trimesh(vertices=vertices_3d, faces=simplices, process=False)

    # Fallback: build structured grid
    nx = int(round(xsize / (points_xy[:, 0].ptp() / max(1, len(np.unique(points_xy[:, 0])) - 1))))
    ny = int(round(ysize / (points_xy[:, 1].ptp() / max(1, len(np.unique(points_xy[:, 1])) - 1))))
    xs = np.linspace(0.0, xsize, nx + 1)
    ys = np.linspace(0.0, ysize, ny + 1)
    vertices = np.array([(x, y, 0.0) for x in xs for y in ys], dtype=float)
    faces: List[Tuple[int, int, int]] = []

    def vid(i: int, j: int) -> int:
        return i * (ny + 1) + j

    for i in range(nx):
        for j in range(ny):
            v0 = vid(i, j)
            v1 = vid(i + 1, j)
            v2 = vid(i + 1, j + 1)
            v3 = vid(i, j + 1)
            faces.append((v0, v1, v2))
            faces.append((v0, v2, v3))

    return trimesh.Trimesh(vertices=vertices, faces=np.array(faces), process=False)


def _generate_ground(
    existing: Optional["trimesh.Trimesh"],
    xsize: float,
    ysize: float,
    edgelength: float,
    preserve_existing_edges: bool = False,
) -> "trimesh.Trimesh":
    """
    Port of MATLAB generateGround: add ground facets and merge with existing mesh.
    """
    if preserve_existing_edges and existing is not None and len(existing.faces) > 0:
        constrained, _ = _generate_ground_with_constraints(existing, xsize, ysize, edgelength)
        if constrained is not None:
            return constrained

    ground_points = np.empty((0, 2))
    if existing is not None and len(existing.faces) > 0:
        ground_vert_ids = np.where(np.isclose(existing.vertices[:, 2], 0.0))[0]
        ground_points = existing.vertices[ground_vert_ids][:, :2]

    grid_pts = _structured_ground_points(xsize, ysize, edgelength)
    all_points = np.vstack([ground_points, grid_pts])
    ground_mesh = _triangulate_ground(all_points, xsize, ysize)

    if existing is None or len(existing.faces) == 0:
        return ground_mesh

    combined = trimesh.util.concatenate([existing, ground_mesh])

    # Remove ground faces under buildings (requires watertightness)
    verts = combined.vertices
    faces = combined.faces
    ground_mask = np.all(np.isclose(verts[faces][:, :, 2], 0.0), axis=1)

    if existing.is_watertight:
        ground_centers = verts[faces[ground_mask]].mean(axis=1)
        try:
            inside = existing.contains(ground_centers)
            ground_mask_indices = np.where(ground_mask)[0]
            faces = np.delete(faces, ground_mask_indices[inside], axis=0)
        except Exception:
            pass

    verts, faces = _remove_unused(verts, faces)
    return trimesh.Trimesh(vertices=verts, faces=faces, process=False)


# -----------------------------------------------------------------------------#
# Public API
# -----------------------------------------------------------------------------#

def create_flat_surface(xsize: float, ysize: float, edgelength: float) -> UDGeom:
    """
    Create a flat surface consisting of triangular facets (MATLAB createFlatSurface).
    """
    ground = _generate_ground(None, xsize, ysize, edgelength)
    return UDGeom(stl=ground)


def create_canyons(
    xsize: float,
    ysize: float,
    B: float,
    W: float,
    H: float,
    shift: float,
    edgelength: float,
    rotate90: bool = False,
) -> UDGeom:
    """
    Create one-dimensional street canyons (MATLAB createCanyons).
    """
    L = B + W
    divisions = int(round(B / edgelength))
    if divisions > 2:
        raise ValueError("divisions must be <= 2 (B/edgelength <= 2)")
    if divisions == 2:
        L = L / 2.0

    Nx = xsize / (B + W)
    Ny = ysize / L
    if abs(Nx - round(Nx)) > 1e-9 or abs(Ny - round(Ny)) > 1e-9:
        raise ValueError("The domain size should be a multiple of canyon width/length")
    Nx = int(round(Nx))
    Ny = int(round(Ny))

    # Base canyon vertices (floors, walls, roof)
    points_floor1 = np.array([[0, 0, 0], [W / 2, 0, 0], [W / 2, L, 0], [0, L, 0]], dtype=float)
    points_floor2 = points_floor1 + np.array([W / 2 + B, 0, 0])
    points_wall1 = np.array([[0, 0, 0], [0, 0, H], [0, L, H], [0, L, 0]], dtype=float) + np.array([W / 2, 0, 0])
    points_wall2 = np.flip(points_wall1, axis=0) + np.array([B, 0, 0])
    points_roof = np.array([[0, 0, 0], [B, 0, 0], [B, L, 0], [0, L, 0]], dtype=float) + np.array([0, 0, H]) + np.array([W / 2, 0, 0])
    points_unit = np.vstack([points_floor1, points_wall1, points_roof, points_wall2, points_floor2])

    conn: List[Tuple[int, int, int]] = []
    points: List[Tuple[float, float, float]] = []

    def add_point(pt: np.ndarray) -> int:
        p = tuple(pt.tolist())
        if p in point_index:
            return point_index[p]
        idx = len(points)
        points.append(p)
        point_index[p] = idx
        return idx

    point_index = {}

    for i in range(Nx):
        for j in range(Ny):
            for n in range(5):
                locs = []
                for m in range(4):
                    pt = points_unit[4 * n + m] + np.array([(i) * (B + W), (j) * L, 0])
                    locs.append(add_point(pt))
                conn.append((locs[0], locs[1], locs[2]))
                conn.append((locs[0], locs[2], locs[3]))

    vertices = np.array(points, dtype=float)
    faces = np.array(conn, dtype=int)

    # shift interior points in x (matching MATLAB behavior)
    x_max = vertices[:, 0].max() if len(vertices) else 0.0
    mask_shift = (vertices[:, 0] > 0.0) & (vertices[:, 0] < x_max)
    vertices[mask_shift, 0] += shift

    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    if divisions >= 1:
        mesh = _divide_faces(mesh, times=1)

    if rotate90:
        angle_rad = np.pi / 2.0
        R = np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0.0],
                      [np.sin(angle_rad), np.cos(angle_rad), 0.0],
                      [0.0, 0.0, 1.0]])
        rot = np.block([[R, np.zeros((3, 1))], [np.zeros((1, 3)), np.ones((1, 1))]])
        # Match MATLAB: rotate first, then shift by xsize in +x
        mesh.apply_transform(rot)
        mesh.apply_translation([xsize, 0.0, 0.0])

    # Remove existing ground facets (z == 0) before adding merged ground
    #verts = mesh.vertices
    #faces = mesh.faces
    #ground_mask = np.all(np.isclose(verts[faces][:, :, 2], 0.0), axis=1)
    #if np.any(ground_mask):
    #    faces = faces[~ground_mask]
    #    verts, faces = _remove_unused(verts, faces)
    #    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)

    # Ground generation disabled (request): return canyon mesh directly
    return UDGeom(stl=mesh)


def create_cubes(
    xsize: float,
    ysize: float,
    Hx: float,
    Hy: float,
    Hz: float,
    Cx: float,
    Cy: float,
    geom_option: str,
    edgelength: float,
) -> UDGeom:
    """
    Create cubes: single, aligned, or staggered (MATLAB createCubes).

    Parameters
    ----------
    xsize, ysize : float
        Domain size in meters (x, y).
    Hx, Hy, Hz : float
        Cube lengths in x, y and height z.
    Cx, Cy : float
        Spacing between cubes in x and y directions (ignored for single cube).
    geom_option : {'S','AC','SC'}
        'S'  - single cube centred in the domain.
        'AC' - aligned array on a regular grid.
        'SC' - staggered array (alternate rows shifted by half spacing).
    edgelength : float
        Target facet size; controls subdivision level.

    Raises
    ------
    ValueError
        If domain is not an integer multiple of cube+spacing for AC/SC,
        or if geom_option is invalid.
    """
    divisions = int(round(Hx / edgelength))
    tol = 0.0

    opt = geom_option.upper()
    if opt not in {"S", "AC", "SC"}:
        raise ValueError("geom_option must be 'S', 'AC', or 'SC'")

    # Build cube arrays
    if opt == "SC":
        Nx = xsize / (Hx + Cx)
        Ny = ysize / (Hy + Cy)
        if abs(Nx - round(Nx)) > 1e-9 or abs(Ny - round(Ny)) > 1e-9:
            raise ValueError("The domain size should be a multiple of cube width + canyon width")
        Nx = int(round(Nx))
        Ny = int(round(Ny))

        shifts = []
        for i in range(1, Nx + 1):
            for j in range(1, Ny + 1):
                if i % 2 == 1:
                    shifts.append([i * (Cx + Hx) - Hx / 2 - Cx / 2, j * (Cy + Hy) - Hy / 2 - Cy / 2, Hz / 2])
                else:
                    shifts.append([i * (Cx + Hx) - Hx / 2 - Cx / 2, j * (Cy + Hy), Hz / 2])
        array = _operate_unit_cube([Hx, Hy, Hz], np.array(shifts), divisions)
        mesh = array

        # Trim outside domain in y
        verts = mesh.vertices.copy()
        faces = mesh.faces.copy()
        if divisions < 2:
            verts[:, 1] = np.clip(verts[:, 1], 0.0, ysize)

        points_remove = np.where((verts[:, 1] <= 0.0) | (verts[:, 1] >= ysize))[0]
        if len(points_remove) > 0:
            mask = np.all(np.isin(faces, points_remove), axis=1)
            faces = faces[~mask]
        verts, faces = _remove_unused(verts, faces)
        mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)

    elif opt == "AC":
        Nx = xsize / (Hx + Cx)
        Ny = ysize / (Hy + Cy)
        if abs(Nx - round(Nx)) > 1e-9 or abs(Ny - round(Ny)) > 1e-9:
            raise ValueError("The domain size should be a multiple of cube width + canyon width")
        Nx = int(round(Nx))
        Ny = int(round(Ny))

        shifts = []
        for i in range(1, Nx + 1):
            for j in range(1, Ny + 1):
                shifts.append([i * (Cx + Hx) - Hx / 2 - Cx / 2, j * (Cy + Hy) - Hy / 2 - Cy / 2, Hz / 2])
        mesh = _operate_unit_cube([Hx, Hy, Hz], np.array(shifts), divisions)

    else:  # 'S'
        shift = np.array([[xsize / 2.0, ysize / 2.0, 0.0]]) + np.array([0.0, 0.0, Hz / 2.0])
        mesh = _operate_unit_cube([Hx, Hy, Hz], shift, divisions)

    # Apply tiny normal shifts (tol) except -z faces (tol is zero here, kept for parity)
    if tol != 0.0:
        normals = mesh.face_normals
        verts = mesh.vertices.copy()
        faces = mesh.faces
        for idx, nrm in enumerate(normals):
            ids = faces[idx]
            if np.allclose(nrm, [1, 0, 0]):
                verts[ids] += tol * np.array([1, 0, 0])
            elif np.allclose(nrm, [-1, 0, 0]):
                verts[ids] += tol * np.array([-1, 0, 0])
            elif np.allclose(nrm, [0, 1, 0]):
                verts[ids] += tol * np.array([0, 1, 0])
            elif np.allclose(nrm, [0, -1, 0]):
                verts[ids] += tol * np.array([0, -1, 0])
            elif np.allclose(nrm, [0, 0, 1]):
                verts[ids] += tol * np.array([0, 0, 1])
        mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)

    ground = _generate_ground(mesh, xsize, ysize, edgelength)
    return UDGeom(stl=ground)


def create_realistic(
    stlfile: Path,
    xsize: float,
    ysize: float,
    shift: Tuple[float, float, float],
    edgelength: float,
    return_debug: bool = False,
) -> Union[UDGeom, tuple[UDGeom, Dict[str, Any]]]:
    """
    Create a realistic urban surface by loading buildings from STL and adding ground.
    Mirrors MATLAB createRealistic.
    """
    building_mesh = trimesh.load_mesh(str(stlfile))
    building_mesh.apply_translation(np.asarray(shift, dtype=float))
    if return_debug:
        ground, debug = _generate_ground_with_constraints(building_mesh, xsize, ysize, edgelength)
        if ground is None:
            ground = _generate_ground(building_mesh, xsize, ysize, edgelength, preserve_existing_edges=False)
        return UDGeom(stl=ground), debug
    ground = _generate_ground(building_mesh, xsize, ysize, edgelength, preserve_existing_edges=True)
    return UDGeom(stl=ground)


def plot_ground_generation_step1(debug: Dict[str, Any], show: bool = True):
    """Visualize the constrained-ground step-1 result (outside regions + clipped grid)."""
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise ImportError("plotly is required for this visualization. Install with: pip install plotly") from exc

    def _iter_polygons(geom):
        if geom is None or getattr(geom, "is_empty", True):
            return []
        if geom.geom_type == "Polygon":
            return [geom]
        return [g for g in getattr(geom, "geoms", []) if g.geom_type == "Polygon"]

    fig = go.Figure()

    outside_region = debug.get("outside_region")
    for idx, poly in enumerate(_iter_polygons(outside_region)):
        x, y = poly.exterior.xy
        fig.add_trace(
            go.Scatter(
                x=list(x),
                y=list(y),
                mode="lines",
                fill="toself",
                line=dict(color="rgba(120,120,120,0.8)", width=1),
                fillcolor="rgba(210,210,210,0.55)",
                name="outside region" if idx == 0 else None,
                showlegend=idx == 0,
            )
        )
        for ring in poly.interiors:
            rx, ry = ring.xy
            fig.add_trace(
                go.Scatter(
                    x=list(rx),
                    y=list(ry),
                    mode="lines",
                    line=dict(color="rgba(0,0,0,0.9)", width=1.5),
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    for idx, line in enumerate(debug.get("grid_lines", [])):
        x, y = line.xy
        fig.add_trace(
            go.Scatter(
                x=list(x),
                y=list(y),
                mode="lines",
                line=dict(color="rgba(60,90,160,0.7)", width=1),
                name="clipped grid" if idx == 0 else None,
                showlegend=idx == 0,
                hoverinfo="skip",
            )
        )

    fig.update_layout(
        title="Ground Generation Step 1",
        xaxis_title="x (m)",
        yaxis_title="y (m)",
        yaxis=dict(scaleanchor="x", scaleratio=1),
        margin=dict(l=20, r=20, b=20, t=50, pad=0),
        template="plotly_white",
    )
    if show:
        fig.show()
    return fig


__all__ = [
    "create_flat_surface",
    "create_canyons",
    "create_cubes",
    "create_realistic",
    "plot_ground_generation_step1",
]
