"""
Geometry generation helpers ported from MATLAB udgeom.

These functions mirror the MATLAB implementations in tools/matlab/+udgeom:
- createFlatSurface.m
- createCanyons.m
- createCubes.m
- createRealistic.m
"""

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union
import warnings

import numpy as np

from exceptions import DependencyError
try:
    import trimesh
except ImportError as exc:
    raise DependencyError("trimesh is required for geometry generation; install with `pip install trimesh`.") from exc

try:
    from scipy.spatial import Delaunay
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from shapely.ops import polygonize, triangulate, unary_union
import triangle as triangle_lib

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


def _linework_to_triangle_pslg(
    lines: List["LineString"],
    hole_points: Optional[List[Tuple[float, float]]] = None,
    free_points: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """Convert split linework to Triangle PSLG input."""
    vertices: List[List[float]] = []
    segments: List[List[int]] = []
    segment_set: set[Tuple[int, int]] = set()
    vertex_map: Dict[Tuple[float, float], int] = {}

    def add_vertex(pt: Sequence[float]) -> int:
        key = (round(float(pt[0]), 12), round(float(pt[1]), 12))
        if key in vertex_map:
            return vertex_map[key]
        idx = len(vertices)
        vertices.append([key[0], key[1]])
        vertex_map[key] = idx
        return idx

    for geom in lines:
        if getattr(geom, "geom_type", None) != "LineString":
            continue
        coords = np.asarray(geom.coords, dtype=float)
        if len(coords) < 2:
            continue
        ids = [add_vertex(pt) for pt in coords]
        for i in range(len(ids) - 1):
            if ids[i] == ids[i + 1]:
                continue
            seg = tuple(sorted((ids[i], ids[i + 1])))
            if seg in segment_set:
                continue
            segment_set.add(seg)
            segments.append([ids[i], ids[i + 1]])

    if free_points is not None:
        for pt in np.asarray(free_points, dtype=float):
            add_vertex(pt)

    data: Dict[str, np.ndarray] = {
        "vertices": np.asarray(vertices, dtype=float),
        "segments": np.asarray(segments, dtype=int),
    }
    if hole_points:
        data["holes"] = np.asarray(hole_points, dtype=float)
    return data


def _points_constraints_to_triangle_pslg(
    constrained_points: np.ndarray,
    constraints: np.ndarray,
    free_points: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """Build Triangle PSLG input from explicit point and segment arrays."""
    vertices: List[List[float]] = []
    segments: List[List[int]] = []
    vertex_map: Dict[Tuple[float, float], int] = {}
    remap: Dict[int, int] = {}

    def add_vertex(pt: Sequence[float]) -> int:
        key = (round(float(pt[0]), 12), round(float(pt[1]), 12))
        if key in vertex_map:
            return vertex_map[key]
        idx = len(vertices)
        vertices.append([key[0], key[1]])
        vertex_map[key] = idx
        return idx

    constrained_points = np.asarray(constrained_points, dtype=float)
    constraints = np.asarray(constraints, dtype=int)

    for idx, pt in enumerate(constrained_points):
        remap[idx] = add_vertex(pt)

    for seg in constraints:
        a = remap[int(seg[0])]
        b = remap[int(seg[1])]
        if a != b:
            segments.append([a, b])

    if free_points is not None:
        for pt in np.asarray(free_points, dtype=float):
            add_vertex(pt)

    data: Dict[str, np.ndarray] = {
        "vertices": np.asarray(vertices, dtype=float),
        "segments": np.asarray(segments, dtype=int),
    }
    return data


def _mesh_edges_from_simplices(simplices: np.ndarray) -> set[Tuple[int, int]]:
    """Return the unique undirected edges of a triangulation."""
    edges: set[Tuple[int, int]] = set()
    for tri in np.asarray(simplices, dtype=int):
        edges.add(tuple(sorted((int(tri[0]), int(tri[1])))))
        edges.add(tuple(sorted((int(tri[1]), int(tri[2])))))
        edges.add(tuple(sorted((int(tri[2]), int(tri[0])))))
    return edges


def _segment_chain_exists(
    vertices: np.ndarray,
    edges: set[Tuple[int, int]],
    seg: np.ndarray,
    tol: float = 1e-9,
) -> bool:
    """Check that a PSLG segment survives as a chain of collinear mesh edges."""
    i0, i1 = int(seg[0]), int(seg[1])
    p0 = np.asarray(vertices[i0], dtype=float)
    p1 = np.asarray(vertices[i1], dtype=float)
    d = p1 - p0
    denom = float(np.dot(d, d))
    if denom <= tol:
        return False

    on_segment: List[Tuple[float, int]] = []
    for idx, pt in enumerate(np.asarray(vertices, dtype=float)):
        v = pt - p0
        cross = d[0] * v[1] - d[1] * v[0]
        if abs(cross) > tol:
            continue
        t = float(np.dot(v, d) / denom)
        if -tol <= t <= 1.0 + tol:
            on_segment.append((t, idx))

    if len(on_segment) < 2:
        return False

    on_segment.sort(key=lambda item: item[0])
    chain = [idx for _, idx in on_segment]
    if chain[0] != i0 or chain[-1] != i1:
        return False

    for a, b in zip(chain[:-1], chain[1:]):
        if tuple(sorted((a, b))) not in edges:
            return False
    return True


def _validate_triangle_boundary(result: Dict[str, np.ndarray], pslg: Dict[str, np.ndarray]) -> None:
    """Verify that every PSLG segment is represented in the output mesh."""
    vertices = np.asarray(result.get("vertices", []), dtype=float)
    simplices = np.asarray(result.get("triangles", []), dtype=int)
    segments = np.asarray(pslg.get("segments", []), dtype=int)
    if len(vertices) == 0 or len(simplices) == 0 or len(segments) == 0:
        raise RuntimeError("Triangle output is incomplete; cannot validate boundary constraints.")

    edges = _mesh_edges_from_simplices(simplices)
    missing = [seg.tolist() for seg in segments if not _segment_chain_exists(vertices, edges, seg)]
    if missing:
        raise RuntimeError(f"Triangle output does not preserve {len(missing)} boundary segments.")


def _triangle_incenters(triangles: np.ndarray) -> np.ndarray:
    """Return incenters for an array of triangles with shape (n, 3, d)."""
    edge_a = np.linalg.norm(triangles[:, 1] - triangles[:, 2], axis=1)
    edge_b = np.linalg.norm(triangles[:, 0] - triangles[:, 2], axis=1)
    edge_c = np.linalg.norm(triangles[:, 0] - triangles[:, 1], axis=1)
    perimeter = edge_a + edge_b + edge_c
    valid = perimeter > 0.0
    incenters = np.zeros((len(triangles), triangles.shape[2]), dtype=float)
    if np.any(valid):
        incenters[valid] = (
            edge_a[valid, None] * triangles[valid, 0]
            + edge_b[valid, None] * triangles[valid, 1]
            + edge_c[valid, None] * triangles[valid, 2]
        ) / perimeter[valid, None]
    return incenters


def _validate_triangle_domain(result: Dict[str, np.ndarray], piece: "Polygon", tol: float = 1e-8) -> None:
    """Verify that the triangle union matches the intended polygonal domain."""
    vertices = np.asarray(result.get("vertices", []), dtype=float)
    simplices = np.asarray(result.get("triangles", []), dtype=int)
    if len(vertices) == 0 or len(simplices) == 0:
        raise RuntimeError("Triangle output is incomplete; cannot validate domain coverage.")

    polys = []
    for simplex in simplices:
        tri = vertices[np.asarray(simplex, dtype=int)]
        poly = Polygon(tri)
        if poly.area > tol:
            polys.append(poly)
    union = unary_union(polys)
    mismatch = union.symmetric_difference(piece).area
    if mismatch > max(tol, 1e-8 * max(1.0, piece.area)):
        raise RuntimeError(f"Triangle output does not match target domain; mismatch area={mismatch:.6g}.")


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


def _matlab_ground_constraints(existing: Optional["trimesh.Trimesh"]) -> Tuple[np.ndarray, np.ndarray]:
    """Return exact z=0 points and segment constraints from the existing mesh."""
    if existing is None or len(existing.faces) == 0:
        return np.empty((0, 2), dtype=float), np.empty((0, 2), dtype=int)

    verts = np.asarray(existing.vertices, dtype=float)
    ground_vertex_ids = np.where(np.isclose(verts[:, 2], 0.0))[0]
    if len(ground_vertex_ids) == 0:
        return np.empty((0, 2), dtype=float), np.empty((0, 2), dtype=int)

    ground_points = verts[ground_vertex_ids][:, :2]
    remap = {int(old): new for new, old in enumerate(ground_vertex_ids)}

    constraints: List[Tuple[int, int]] = []
    for edge in np.asarray(existing.edges_unique, dtype=int):
        a, b = int(edge[0]), int(edge[1])
        if a in remap and b in remap:
            ra = remap[a]
            rb = remap[b]
            if ra != rb:
                constraints.append((ra, rb))

    if not constraints:
        return ground_points, np.empty((0, 2), dtype=int)
    return ground_points, np.asarray(constraints, dtype=int)


def _matlab_domain_boundary(xsize: float, ysize: float, edgelength: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return MATLAB-style domain perimeter points and constraints."""
    x0 = 0.0
    y0 = 0.0
    nedges_x = int(round(xsize / edgelength))
    nedges_y = int(round(ysize / edgelength))

    points: List[Tuple[float, float]] = []
    constraints: List[Tuple[int, int]] = []

    def add_point(x: float, y: float) -> None:
        points.append((float(x), float(y)))

    for i in range(nedges_x + 1):
        add_point(x0 + i * edgelength, y0)
        n = len(points)
        constraints.append((n - 1, n))
    for i in range(1, nedges_y + 1):
        add_point(x0 + xsize, y0 + i * edgelength)
        n = len(points)
        constraints.append((n - 1, n))
    for i in range(1, nedges_x + 1):
        add_point(x0 + xsize - i * edgelength, y0 + ysize)
        n = len(points)
        constraints.append((n - 1, n))
    for i in range(1, nedges_y):
        add_point(x0, y0 + ysize - i * edgelength)
        n = len(points)
        if i < nedges_y - 1:
            constraints.append((n - 1, n))

    if points:
        constraints.append((len(points) - 1, 0))

    return np.asarray(points, dtype=float), np.asarray(constraints, dtype=int)


def _matlab_interior_points(xsize: float, ysize: float, edgelength: float) -> np.ndarray:
    """Return MATLAB-style regular interior points."""
    nedges_x = int(round(xsize / edgelength))
    nedges_y = int(round(ysize / edgelength))
    points = []
    for i in range(1, nedges_x + 1):
        for j in range(1, nedges_y + 1):
            points.append((float(i * edgelength), float(j * edgelength)))
    return np.asarray(points, dtype=float) if points else np.empty((0, 2), dtype=float)


def _structured_quad_mesh(
    p00: np.ndarray,
    p10: np.ndarray,
    p11: np.ndarray,
    p01: np.ndarray,
    nu: int,
    nv: int,
) -> "trimesh.Trimesh":
    """Create a structured triangular mesh on a quadrilateral patch."""
    nu = max(int(nu), 1)
    nv = max(int(nv), 1)
    us = np.linspace(0.0, 1.0, nu + 1)
    vs = np.linspace(0.0, 1.0, nv + 1)

    vertices = []
    for i, u in enumerate(us):
        for j, v in enumerate(vs):
            vertices.append(
                (1.0 - u) * (1.0 - v) * p00
                + u * (1.0 - v) * p10
                + u * v * p11
                + (1.0 - u) * v * p01
            )

    def vid(i: int, j: int) -> int:
        return i * (nv + 1) + j

    faces: List[Tuple[int, int, int]] = []
    for i in range(nu):
        for j in range(nv):
            v00 = vid(i, j)
            v10 = vid(i + 1, j)
            v11 = vid(i + 1, j + 1)
            v01 = vid(i, j + 1)
            faces.append((v00, v10, v11))
            faces.append((v00, v11, v01))

    return trimesh.Trimesh(
        vertices=np.asarray(vertices, dtype=float),
        faces=np.asarray(faces, dtype=int),
        process=False,
    )


def _remove_under_building_ground_faces(
    combined: "trimesh.Trimesh",
    existing: "trimesh.Trimesh",
) -> "trimesh.Trimesh":
    """Mirror MATLAB ground removal: remove z=0 facets whose incenters lie under buildings."""
    footprint_union = _ground_footprint_union(existing)
    if footprint_union is None or footprint_union.is_empty:
        return combined

    verts = np.asarray(combined.vertices, dtype=float)
    faces = np.asarray(combined.faces, dtype=int)
    if len(faces) == 0:
        return combined

    ground_mask = np.all(np.isclose(verts[faces][:, :, 2], 0.0), axis=1)
    ground_indices = np.where(ground_mask)[0]
    if len(ground_indices) == 0:
        return combined

    ground_tris = verts[faces[ground_indices]][:, :, :2]
    incenters = _triangle_incenters(ground_tris)
    remove_local = np.array(
        [footprint_union.contains(Point(float(pt[0]), float(pt[1]))) for pt in incenters],
        dtype=bool,
    )
    if not np.any(remove_local):
        return combined

    keep_faces = np.ones(len(faces), dtype=bool)
    keep_faces[ground_indices[remove_local]] = False
    filtered_faces = faces[keep_faces]

    if len(filtered_faces) == 0:
        return trimesh.Trimesh(
            vertices=np.empty((0, 3), dtype=float),
            faces=np.empty((0, 3), dtype=int),
            process=False,
        )

    degenerate = np.array([len(np.unique(face)) < 3 for face in filtered_faces], dtype=bool)
    if np.any(degenerate):
        filtered_faces = filtered_faces[~degenerate]

    verts_arr, faces_arr = _remove_unused(verts, filtered_faces)
    return trimesh.Trimesh(vertices=verts_arr, faces=faces_arr, process=False)


def _generate_ground_matlab_style(
    existing: Optional["trimesh.Trimesh"],
    xsize: float,
    ysize: float,
    edgelength: float,
    validate_boundary: bool = True,
    split_existing_boundary: bool = True,
) -> Tuple["trimesh.Trimesh", Dict[str, Any]]:
    """Generate ground using the same high-level workflow as MATLAB generateGround."""
    existing_points, existing_constraints = _matlab_ground_constraints(existing)
    edge_points, edge_constraints = _matlab_domain_boundary(xsize, ysize, edgelength)
    interior_points = _matlab_interior_points(xsize, ysize, edgelength)

    if len(existing_points) > 0:
        constrained_points = np.vstack([existing_points, edge_points])
        constrained_constraints = np.vstack(
            [
                existing_constraints,
                edge_constraints + len(existing_points),
            ]
        ) if len(edge_constraints) else existing_constraints.copy()
    else:
        constrained_points = edge_points.copy()
        constrained_constraints = edge_constraints.copy()

    debug: Dict[str, Any] = {
        "ground_points": existing_points,
        "ground_constraints": existing_constraints,
        "ground_points_edges": edge_points,
        "ground_constraints_edges": edge_constraints,
        "ground_points_inside": interior_points,
    }

    domain = Polygon([(0.0, 0.0), (xsize, 0.0), (xsize, ysize), (0.0, ysize)])
    footprint_union = _ground_footprint_union(existing) if existing is not None else None
    outside_region = domain.difference(footprint_union) if footprint_union is not None else domain
    xs = _grid_axis(xsize, edgelength)
    ys = _grid_axis(ysize, edgelength)
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
    debug.update(
        {
            "domain": domain,
            "footprint_union": footprint_union,
            "outside_region": outside_region,
            "grid_lines": clipped_grid_lines,
            "planar_cells": [],
        }
    )

    if len(constrained_points) == 0:
        empty = trimesh.Trimesh(
            vertices=np.empty((0, 3), dtype=float),
            faces=np.empty((0, 3), dtype=int),
            process=False,
        )
        return empty, debug

    pslg = _points_constraints_to_triangle_pslg(
        constrained_points,
        constrained_constraints,
        free_points=interior_points,
    )
    result = triangle_lib.triangulate(pslg, "pQe")
    vertices = np.asarray(result.get("vertices", []), dtype=float)
    simplices = np.asarray(result.get("triangles", []), dtype=int)
    if validate_boundary:
        _validate_triangle_boundary(result, pslg)
    if len(vertices) == 0 or len(simplices) == 0:
        ground_mesh = trimesh.Trimesh(
            vertices=np.empty((0, 3), dtype=float),
            faces=np.empty((0, 3), dtype=int),
            process=False,
        )
    else:
        ground_mesh = trimesh.Trimesh(
            vertices=np.column_stack([vertices, np.zeros(len(vertices))]),
            faces=simplices,
            process=False,
        )
    debug["triangle_pslg"] = pslg

    if existing is None or len(existing.faces) == 0:
        return ground_mesh, debug

    existing_for_merge = _split_existing_ground_boundary(existing, ground_mesh) if split_existing_boundary else existing
    combined = trimesh.util.concatenate([existing_for_merge, ground_mesh])
    combined.merge_vertices(digits_vertex=12)
    combined = _remove_under_building_ground_faces(combined, existing)
    verts_arr, faces_arr = _remove_unused(np.asarray(combined.vertices), np.asarray(combined.faces))
    return trimesh.Trimesh(vertices=verts_arr, faces=faces_arr, process=False), debug


def _generate_ground_with_constraints(
    existing: "trimesh.Trimesh",
    xsize: float,
    ysize: float,
    edgelength: float,
    validate_boundary: bool = True,
    footprint_union_override=None,
    split_existing_boundary: bool = True,
    triangle_options: Optional[str] = None,
) -> tuple[Optional["trimesh.Trimesh"], Dict[str, Any]]:
    """
    Build a constrained ground mesh using a polygonized grid plus building footprint edges.

    This follows the MATLAB `generateGround` intent more closely than an unconstrained
    Delaunay triangulation: building footprint edges remain hard constraints, so ground
    triangles cannot cut across building bases.
    """
    domain = Polygon([(0.0, 0.0), (xsize, 0.0), (xsize, ysize), (0.0, ysize)])
    xs = _grid_axis(xsize, edgelength)
    ys = _grid_axis(ysize, edgelength)

    footprint_union = footprint_union_override if footprint_union_override is not None else _ground_footprint_union(existing)
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

    lines = boundary_lines

    hole_points: List[Tuple[float, float]] = []
    if footprint_union is not None:
        footprint_geoms = [footprint_union] if footprint_union.geom_type == "Polygon" else list(footprint_union.geoms)
        for geom in footprint_geoms:
            rp = geom.representative_point()
            hole_points.append((float(rp.x), float(rp.y)))

    grid_pts = _structured_ground_points(xsize, ysize, edgelength)
    free_points = []
    for pt in grid_pts:
        if outside_region.contains(Point(float(pt[0]), float(pt[1]))):
            free_points.append(pt)
    free_points_arr = np.asarray(free_points, dtype=float) if free_points else np.empty((0, 2), dtype=float)

    pslg = _linework_to_triangle_pslg(lines, hole_points=hole_points, free_points=free_points_arr)
    max_area = np.sqrt(3.0) * max(float(edgelength), 1e-6) ** 2 / 4.0
    tri_opts = triangle_options or f"pq18a{max_area:.12g}e"
    result = triangle_lib.triangulate(pslg, tri_opts)
    if validate_boundary:
        _validate_triangle_boundary(result, pslg)
    _validate_triangle_domain(result, outside_region)

    vertices = np.asarray(result.get("vertices", []), dtype=float)
    simplices = np.asarray(result.get("triangles", []), dtype=int)
    if len(vertices) == 0 or len(simplices) == 0:
        return None, {
            "domain": domain,
            "footprint_union": footprint_union,
            "outside_region": outside_region,
            "grid_lines": clipped_grid_lines,
            "planar_cells": [],
        }

    ground_mesh = trimesh.Trimesh(
        vertices=np.column_stack([vertices, np.zeros(len(vertices))]),
        faces=simplices,
        process=False,
    )
    if len(ground_mesh.faces) == 0:
        return None, {
            "domain": domain,
            "footprint_union": footprint_union,
            "outside_region": outside_region,
            "grid_lines": clipped_grid_lines,
            "planar_cells": [],
        }

    existing_for_merge = _split_existing_ground_boundary(existing, ground_mesh) if split_existing_boundary else existing
    combined = trimesh.util.concatenate([existing_for_merge, ground_mesh])
    combined.merge_vertices(digits_vertex=12)
    verts, faces = _remove_unused(np.asarray(combined.vertices), np.asarray(combined.faces))
    return trimesh.Trimesh(vertices=verts, faces=faces, process=False), {
        "domain": domain,
        "footprint_union": footprint_union,
        "outside_region": outside_region,
        "grid_lines": clipped_grid_lines,
        "planar_cells": [],
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
    nx = int(round(xsize / (np.ptp(points_xy[:, 0]) / max(1, len(np.unique(points_xy[:, 0])) - 1))))
    ny = int(round(ysize / (np.ptp(points_xy[:, 1]) / max(1, len(np.unique(points_xy[:, 1])) - 1))))
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
        constrained, _ = _generate_ground_matlab_style(existing, xsize, ysize, edgelength)
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
    Nx = xsize / (B + W)
    Ny = ysize / L
    if abs(Nx - round(Nx)) > 1e-9 or abs(Ny - round(Ny)) > 1e-9:
        raise ValueError("The domain size should be a multiple of canyon width/length")
    Nx = int(round(Nx))
    Ny = int(round(Ny))

    ny_edges = max(1, int(np.ceil(ysize / edgelength)))
    nz_edges = max(1, int(np.ceil(H / edgelength)))
    nx_edges = max(1, int(np.ceil(B / edgelength)))

    x_lefts = np.array([i * (B + W) + W / 2.0 for i in range(Nx)], dtype=float)
    x_rights = x_lefts + B
    mask_left = (x_lefts > 0.0) & (x_lefts < xsize)
    mask_right = (x_rights > 0.0) & (x_rights < xsize)
    x_lefts = x_lefts.copy()
    x_rights = x_rights.copy()
    x_lefts[mask_left] += shift
    x_rights[mask_right] += shift

    strip_meshes = []
    for x0, x1 in zip(x_lefts, x_rights):
        left_wall = _structured_quad_mesh(
            np.array([x0, 0.0, 0.0], dtype=float),
            np.array([x0, ysize, 0.0], dtype=float),
            np.array([x0, ysize, H], dtype=float),
            np.array([x0, 0.0, H], dtype=float),
            ny_edges,
            nz_edges,
        )
        right_wall = _structured_quad_mesh(
            np.array([x1, ysize, 0.0], dtype=float),
            np.array([x1, 0.0, 0.0], dtype=float),
            np.array([x1, 0.0, H], dtype=float),
            np.array([x1, ysize, H], dtype=float),
            ny_edges,
            nz_edges,
        )
        roof = _structured_quad_mesh(
            np.array([x0, 0.0, H], dtype=float),
            np.array([x1, 0.0, H], dtype=float),
            np.array([x1, ysize, H], dtype=float),
            np.array([x0, ysize, H], dtype=float),
            nx_edges,
            ny_edges,
        )
        strip_meshes.extend([left_wall, right_wall, roof])

    mesh = trimesh.util.concatenate(strip_meshes) if strip_meshes else trimesh.Trimesh(
        vertices=np.empty((0, 3), dtype=float),
        faces=np.empty((0, 3), dtype=int),
        process=False,
    )
    mesh.merge_vertices(digits_vertex=12)

    shifted_extent = max(float(shift), 0.0)
    ground_xsize = max(float(np.max(x_rights) + W / 2.0), float(xsize + shifted_extent))
    ground_ysize = float(ysize)

    x_breaks = [0.0]
    for x0, x1 in zip(x_lefts, x_rights):
        x_breaks.extend([float(x0), float(x1)])
    x_breaks.append(float(ground_xsize))
    x_breaks = sorted(dict.fromkeys(x_breaks))

    ground_meshes = []
    for a, b in zip(x_breaks[:-1], x_breaks[1:]):
        if b - a <= 1e-9:
            continue
        inside_building = any((a >= x0 - 1e-9) and (b <= x1 + 1e-9) for x0, x1 in zip(x_lefts, x_rights))
        if inside_building:
            continue
        nu = max(1, int(np.ceil((b - a) / edgelength)))
        nv = max(1, int(np.ceil(ground_ysize / edgelength)))
        ground_meshes.append(
            _structured_quad_mesh(
                np.array([a, 0.0, 0.0], dtype=float),
                np.array([b, 0.0, 0.0], dtype=float),
                np.array([b, ground_ysize, 0.0], dtype=float),
                np.array([a, ground_ysize, 0.0], dtype=float),
                nu,
                nv,
            )
        )

    ground_mesh = trimesh.util.concatenate(ground_meshes) if ground_meshes else trimesh.Trimesh(
        vertices=np.empty((0, 3), dtype=float),
        faces=np.empty((0, 3), dtype=int),
        process=False,
    )

    combined = trimesh.util.concatenate([mesh, ground_mesh])
    combined.merge_vertices(digits_vertex=12)
    verts, faces = _remove_unused(np.asarray(combined.vertices), np.asarray(combined.faces))
    combined_mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)

    if rotate90:
        angle_rad = np.pi / 2.0
        R = np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0.0],
                      [np.sin(angle_rad), np.cos(angle_rad), 0.0],
                      [0.0, 0.0, 1.0]])
        rot = np.block([[R, np.zeros((3, 1))], [np.zeros((1, 3)), np.ones((1, 1))]])
        combined_mesh.apply_transform(rot)
        combined_mesh.apply_translation([xsize, 0.0, 0.0])

    return UDGeom(stl=combined_mesh)


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
            if i % 2 == 0:
                shifts.append([i * (Cx + Hx) - Hx / 2 - Cx / 2, 0.0, Hz / 2])
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

    ground, _ = _generate_ground_matlab_style(mesh, xsize, ysize, edgelength)
    return UDGeom(stl=ground)


def add_ground(
    geom_or_stlfile: Union[Path, str, "trimesh.Trimesh", UDGeom],
    xsize: float,
    ysize: float,
    shift: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    edgelength: float = 1.0,
    return_debug: bool = False,
    preserve_existing_edges: bool = True,
) -> Union[UDGeom, tuple[UDGeom, Dict[str, Any]]]:
    """
    Add a flat ground surface around an existing building mesh or STL geometry.

    Parameters
    ----------
    geom_or_stlfile : Path, str, trimesh.Trimesh, or UDGeom
        Building geometry to surround with ground. This may be:
        - a path to an STL file
        - an existing ``trimesh.Trimesh``
        - an existing ``UDGeom``
    xsize, ysize : float
        Ground-domain size in the x and y directions.
    shift : tuple of float, default=(0, 0, 0)
        Translation applied to the input geometry before ground is added.
        This is mainly useful when loading geometry from an STL file.
    edgelength : float, default=1.0
        Target ground facet size.
    return_debug : bool, default=False
        If True, also return the clipped-grid debug information used in the
        step-1 ground-generation visualisation.
    preserve_existing_edges : bool, default=True
        If True, preserve the building-footprint boundary discretisation when
        stitching the flat ground back into the mesh.

    Returns
    -------
    geom : UDGeom
        Geometry with stitched flat ground added.
    debug : dict, optional
        Returned only when ``return_debug=True``.
    """
    if isinstance(geom_or_stlfile, (str, Path)):
        building_mesh = trimesh.load_mesh(str(geom_or_stlfile))
    elif isinstance(geom_or_stlfile, UDGeom):
        building_mesh = geom_or_stlfile.stl.copy()
    else:
        building_mesh = geom_or_stlfile.copy()

    building_mesh.apply_translation(np.asarray(shift, dtype=float))
    if return_debug:
        ground = _generate_ground(
            building_mesh,
            xsize,
            ysize,
            edgelength,
            preserve_existing_edges=preserve_existing_edges,
        )
        debug = {}
        domain = Polygon([(0.0, 0.0), (xsize, 0.0), (xsize, ysize), (0.0, ysize)])
        footprint_union = _ground_footprint_union(building_mesh)
        outside_region = domain.difference(footprint_union) if footprint_union is not None else domain
        xs = _grid_axis(xsize, edgelength)
        ys = _grid_axis(ysize, edgelength)
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
        debug = {
            "domain": domain,
            "footprint_union": footprint_union,
            "outside_region": outside_region,
            "grid_lines": clipped_grid_lines,
            "planar_cells": [],
        }
        return UDGeom(stl=ground), debug
    ground = _generate_ground(
        building_mesh,
        xsize,
        ysize,
        edgelength,
        preserve_existing_edges=preserve_existing_edges,
    )
    return UDGeom(stl=ground)


def plot_ground_generation_step1(debug: Dict[str, Any], show: bool = True):
    """Visualize the constrained-ground step-1 result (outside regions + clipped grid)."""
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise DependencyError("plotly is required for this visualization. Install with: pip install plotly") from exc

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
    "add_ground",
    "plot_ground_generation_step1",
]
