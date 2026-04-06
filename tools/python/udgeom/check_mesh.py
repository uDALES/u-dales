"""
check_mesh - Validate udgeom/trimesh meshes for common topology and geometry issues.

The checks here are intentionally pragmatic rather than exhaustive. They target
the failure modes that matter most for uDALES geometry preparation:
- unused / disconnected vertices
- multiple disconnected face components
- non-manifold edges
- duplicated faces (a practical proxy for internal back-to-back walls)
- non-finite coordinates
- degenerate / zero-area faces
- invalid face normals
"""

from __future__ import annotations

from typing import Dict, Union

import numpy as np

try:
    import trimesh
    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False

try:
    from shapely.geometry import GeometryCollection, MultiPolygon, Polygon
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:  # pragma: no cover - optional runtime dependency
    SHAPELY_AVAILABLE = False


def _as_trimesh(mesh_or_geom) -> "trimesh.Trimesh":
    if not TRIMESH_AVAILABLE:
        raise ImportError("trimesh is required. Install with: pip install trimesh")

    if isinstance(mesh_or_geom, trimesh.Trimesh):
        return mesh_or_geom

    stl = getattr(mesh_or_geom, "stl", None)
    if isinstance(stl, trimesh.Trimesh):
        return stl

    raise TypeError("Expected a trimesh.Trimesh or a UDGeom-like object with a .stl mesh")


def identify_ground_faces(
    mesh_or_geom,
    max_slope_deg: float = 20.0,
    min_component_area_fraction: float = 0.05,
    elevation_tolerance: float = 1.0e-6,
) -> np.ndarray:
    """
    Identify likely ground faces from geometry rather than absolute height.

    The heuristic is intentionally simple:
    - take near-horizontal faces (upward or downward) using ``abs(n_z)``
    - split them into face-connected components
    - keep only components with non-trivial area
    - choose the lowest remaining component by median face-center z

    This allows topography or elevated domains without assuming ``z == 0``.

    Parameters
    ----------
    mesh_or_geom : trimesh.Trimesh or UDGeom
        Mesh to analyze.
    max_slope_deg : float, default=20.0
        Maximum slope from horizontal for ground candidates.
    min_component_area_fraction : float, default=0.05
        Minimum fraction of total near-horizontal candidate area for a component
        to be considered as possible ground.
    elevation_tolerance : float, default=1e-6
        Components whose median elevation is within this tolerance of the
        lowest candidate are all treated as ground.

    Returns
    -------
    mask : ndarray of bool, shape (n_faces,)
        Boolean mask of likely ground faces.
    """
    mesh = _as_trimesh(mesh_or_geom)
    faces = np.asarray(mesh.faces, dtype=int)
    n_faces = len(faces)
    if n_faces == 0:
        return np.zeros(0, dtype=bool)

    face_normals = np.asarray(mesh.face_normals, dtype=float)
    face_areas = np.asarray(mesh.area_faces, dtype=float)
    face_centers = np.asarray(mesh.triangles_center, dtype=float)

    cos_thresh = float(np.cos(np.deg2rad(max_slope_deg)))
    horizontal_mask = np.isfinite(face_normals).all(axis=1) & np.isfinite(face_areas) & (face_areas > 1.0e-12)
    horizontal_mask &= np.abs(face_normals[:, 2]) >= cos_thresh
    if not np.any(horizontal_mask):
        return np.zeros(n_faces, dtype=bool)

    adjacency = {i: set() for i in range(n_faces)}
    for face_a, face_b in np.asarray(mesh.face_adjacency, dtype=int):
        adjacency[int(face_a)].add(int(face_b))
        adjacency[int(face_b)].add(int(face_a))

    visited = set()
    components = []
    candidate_ids = np.flatnonzero(horizontal_mask)
    for face_id in candidate_ids:
        face_id = int(face_id)
        if face_id in visited:
            continue
        stack = [face_id]
        visited.add(face_id)
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            for nb in adjacency[current]:
                if nb not in visited and horizontal_mask[nb]:
                    visited.add(nb)
                    stack.append(nb)
        components.append(np.asarray(component, dtype=int))

    if not components:
        return np.zeros(n_faces, dtype=bool)

    total_candidate_area = float(np.sum(face_areas[horizontal_mask]))
    area_cutoff = max(0.0, float(min_component_area_fraction)) * total_candidate_area

    viable_components = []
    for comp in components:
        comp_area = float(np.sum(face_areas[comp]))
        if comp_area + 1.0e-12 < area_cutoff:
            continue
        comp_z = float(np.median(face_centers[comp, 2]))
        comp_mean_nz = float(np.mean(face_normals[comp, 2]))
        viable_components.append((comp_z, comp_area, comp_mean_nz, comp))

    if not viable_components:
        viable_components = [
            (
                float(np.median(face_centers[comp, 2])),
                float(np.sum(face_areas[comp])),
                float(np.mean(face_normals[comp, 2])),
                comp,
            )
            for comp in components
        ]

    upward_components = [item for item in viable_components if item[2] > 1.0e-6]
    candidates = upward_components if upward_components else viable_components
    candidates.sort(key=lambda item: (item[0], -item[1]))
    lowest_z = candidates[0][0]
    ground_components = [item[3] for item in candidates if abs(item[0] - lowest_z) <= elevation_tolerance]

    mask = np.zeros(n_faces, dtype=bool)
    for comp in ground_components:
        mask[comp] = True
    return mask


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


def find_internal_touching_wall_regions(
    mesh_or_geom,
    plane_tolerance: float = 1.0e-8,
) -> list[dict]:
    """
    Detect opposite-facing coplanar vertical wall overlaps.

    These are common in adjacent-building exports where both buildings keep the
    shared wall, producing an internal touching surface rather than a clean
    exterior mesh.
    """
    if not SHAPELY_AVAILABLE:
        return []

    mesh = _as_trimesh(mesh_or_geom)
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
        plane_groups.setdefault(key, {}).setdefault(sign, []).append(int(face_id))

    regions = []
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
        if overlap.is_empty:
            continue

        overlap_polys = [poly for poly in _iter_polygons(overlap) if poly.area > 1.0e-8]
        if not overlap_polys:
            continue

        face_ids = sorted(set(by_sign[1]) | set(by_sign[-1]))
        bbox = _bbox_from_face_ids(vertices, faces, np.asarray(face_ids, dtype=int))
        overlap_area = float(sum(poly.area for poly in overlap_polys))
        regions.append(
            {
                "axis": int(axis),
                "plane_value": float(plane_value),
                "n_positive_faces": int(len(by_sign[1])),
                "n_negative_faces": int(len(by_sign[-1])),
                "n_faces": int(len(face_ids)),
                "face_ids": face_ids,
                "bbox": bbox,
                "overlap_area": overlap_area,
            }
        )

    regions.sort(key=lambda item: (-item["overlap_area"], -item["n_faces"]))
    return regions


def _boundary_edge_face_map(mesh: "trimesh.Trimesh") -> dict[tuple[int, int], list[int]]:
    faces = np.asarray(mesh.faces, dtype=int)
    edge_to_faces: dict[tuple[int, int], list[int]] = {}
    for face_id, face in enumerate(faces):
        a, b, c = map(int, face)
        for edge in (tuple(sorted((a, b))), tuple(sorted((b, c))), tuple(sorted((c, a)))):
            edge_to_faces.setdefault(edge, []).append(int(face_id))
    return edge_to_faces


def _canonical_line_direction(vec: np.ndarray, tol: float = 1.0e-12) -> np.ndarray:
    length = float(np.linalg.norm(vec))
    if length <= tol:
        raise ValueError("Zero-length edge cannot define a line direction.")
    direction = np.asarray(vec, dtype=float) / length
    for comp in direction:
        if abs(comp) > tol:
            if comp < 0.0:
                direction = -direction
            break
    return direction


def _line_group_key(p0: np.ndarray, p1: np.ndarray, tol: float) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    direction = _canonical_line_direction(p1 - p0)
    moment = np.cross(p0, direction)
    d_key = tuple(int(x) for x in np.round(direction / tol).astype(np.int64).tolist())
    m_key = tuple(int(x) for x in np.round(moment / tol).astype(np.int64).tolist())
    return d_key, m_key


def _find_unstitched_touching_edge_groups(mesh: "trimesh.Trimesh", point_tolerance: float = 1.0e-8) -> list[dict]:
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)
    edge_to_faces = _boundary_edge_face_map(mesh)
    boundary_edges = [edge for edge, face_ids in edge_to_faces.items() if len(face_ids) == 1]
    if not boundary_edges:
        return []

    line_groups: dict[tuple[tuple[int, int, int], tuple[int, int, int]], list[dict]] = {}
    for edge in boundary_edges:
        p0 = vertices[edge[0]]
        p1 = vertices[edge[1]]
        if not np.isfinite(np.vstack([p0, p1])).all():
            continue
        try:
            key = _line_group_key(p0, p1, point_tolerance)
        except ValueError:
            continue
        direction = _canonical_line_direction(p1 - p0)
        t0 = float(np.dot(p0, direction))
        t1 = float(np.dot(p1, direction))
        line_groups.setdefault(key, []).append(
            {
                "edge": tuple(map(int, edge)),
                "direction": direction,
                "t0": min(t0, t1),
                "t1": max(t0, t1),
            }
        )

    groups: list[dict] = []
    tol = float(point_tolerance)
    for entries in line_groups.values():
        if len(entries) <= 1:
            continue

        adjacency = {i: set() for i in range(len(entries))}
        for i in range(len(entries)):
            for j in range(i + 1, len(entries)):
                overlap = min(entries[i]["t1"], entries[j]["t1"]) - max(entries[i]["t0"], entries[j]["t0"])
                if overlap > tol:
                    adjacency[i].add(j)
                    adjacency[j].add(i)

        visited: set[int] = set()
        for start in range(len(entries)):
            if start in visited or not adjacency[start]:
                continue
            stack = [start]
            visited.add(start)
            idxs: list[int] = []
            while stack:
                current = stack.pop()
                idxs.append(current)
                for nb in adjacency[current]:
                    if nb not in visited:
                        visited.add(nb)
                        stack.append(nb)

            edges = sorted({entries[idx]["edge"] for idx in idxs})
            if len(edges) <= 1:
                continue
            face_ids = sorted({int(edge_to_faces[edge][0]) for edge in edges})
            bbox = _bbox_from_face_ids(vertices, faces, np.asarray(face_ids, dtype=int))
            groups.append(
                {
                    "edges": edges,
                    "face_ids": face_ids,
                    "bbox": bbox,
                }
            )

    groups.sort(key=lambda item: (-len(item["edges"]), -len(item["face_ids"])))
    return groups


def find_unstitched_touching_regions(
    mesh_or_geom,
    point_tolerance: float = 1.0e-8,
) -> list[dict]:
    """
    Detect boundary seams that touch geometrically but are not stitched topologically.

    This targets the case where two surfaces meet along the same segment in space
    but do not share the same mesh edge/vertices.
    """
    mesh = _as_trimesh(mesh_or_geom)
    vertices = np.asarray(mesh.vertices, dtype=float)
    regions = []
    for group in _find_unstitched_touching_edge_groups(mesh, point_tolerance=point_tolerance):
        edges = group["edges"]
        regions.append(
            {
                "n_edges": int(len(edges)),
                "edge_vertex_ids": [tuple(map(int, edge)) for edge in edges],
                "edge_points": [
                    (vertices[edge[0]].astype(float).tolist(), vertices[edge[1]].astype(float).tolist()) for edge in edges
                ],
                "face_ids": group["face_ids"],
                "bbox": group["bbox"],
            }
        )

    regions.sort(key=lambda item: (-item["n_edges"], -len(item["face_ids"])))
    return regions


def find_nonmanifold_regions(mesh_or_geom) -> list[dict]:
    """
    Locate and cluster non-manifold edges for inspection.

    Parameters
    ----------
    mesh_or_geom : trimesh.Trimesh or UDGeom
        Mesh to analyze.

    Returns
    -------
    regions : list of dict
        One entry per clustered non-manifold defect region. Each region contains:
        - ``edge_vertex_ids``: list of vertex-id pairs
        - ``edge_points``: corresponding endpoint coordinates
        - ``face_ids``: sorted incident face ids touching the region
        - ``n_edges``: number of non-manifold edges in the region
        - ``n_faces``: number of incident faces in the region
        - ``bbox``: region bounding box ``[[xmin,ymin,zmin],[xmax,ymax,zmax]]``
    """
    mesh = _as_trimesh(mesh_or_geom)
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)

    edge_to_faces: dict[tuple[int, int], list[int]] = {}
    for face_id, face in enumerate(faces):
        a, b, c = map(int, face)
        for edge in (tuple(sorted((a, b))), tuple(sorted((b, c))), tuple(sorted((c, a)))):
            edge_to_faces.setdefault(edge, []).append(int(face_id))

    nonmanifold = {edge: face_ids for edge, face_ids in edge_to_faces.items() if len(face_ids) > 2}
    if not nonmanifold:
        return []

    edges = list(nonmanifold.keys())
    edge_adjacency = {edge: set() for edge in edges}
    vertex_to_edges: dict[int, list[tuple[int, int]]] = {}
    for edge in edges:
        for vid in edge:
            vertex_to_edges.setdefault(int(vid), []).append(edge)
    for edge_list in vertex_to_edges.values():
        for edge in edge_list:
            edge_adjacency[edge].update(other for other in edge_list if other != edge)

    visited: set[tuple[int, int]] = set()
    regions: list[dict] = []
    for edge in edges:
        if edge in visited:
            continue
        stack = [edge]
        visited.add(edge)
        region_edges: list[tuple[int, int]] = []
        region_faces: set[int] = set()
        region_vertices: set[int] = set()
        while stack:
            current = stack.pop()
            region_edges.append(current)
            region_faces.update(nonmanifold[current])
            region_vertices.update(current)
            for nb in edge_adjacency[current]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)

        region_face_vertex_ids = np.unique(faces[np.asarray(sorted(region_faces), dtype=int)].reshape(-1))
        region_points = vertices[region_face_vertex_ids]
        bbox = np.array([region_points.min(axis=0), region_points.max(axis=0)], dtype=float)
        regions.append(
            {
                "edge_vertex_ids": [tuple(map(int, e)) for e in sorted(region_edges)],
                "edge_points": [
                    (vertices[e[0]].astype(float).tolist(), vertices[e[1]].astype(float).tolist())
                    for e in sorted(region_edges)
                ],
                "face_ids": sorted(int(fid) for fid in region_faces),
                "n_edges": int(len(region_edges)),
                "n_faces": int(len(region_faces)),
                "bbox": bbox,
            }
        )

    regions.sort(key=lambda item: (-item["n_edges"], -item["n_faces"]))
    return regions


def _bbox_from_face_ids(vertices: np.ndarray, faces: np.ndarray, face_ids: np.ndarray) -> np.ndarray:
    face_vertex_ids = np.unique(faces[np.asarray(face_ids, dtype=int)].reshape(-1))
    points = vertices[face_vertex_ids]
    return np.array([points.min(axis=0), points.max(axis=0)], dtype=float)


def calculate_independent_surfaces(mesh_or_geom) -> Dict[str, Union[int, list, np.ndarray]]:
    """
    Partition the mesh into face-connected independent surfaces.

    In this context an "independent surface" means one face-connected piece of
    the STL. For many valid urban meshes this count is 1, but floating
    buildings, detached ground patches, or small disconnected fragments will
    increase the count.

    Parameters
    ----------
    mesh_or_geom : trimesh.Trimesh or UDGeom
        Mesh to analyze.

    Returns
    -------
    result : dict
        Dictionary containing:
        - ``n_surfaces``: number of independent face-connected surfaces
        - ``face_surface_ids``: integer surface id for each face (1-based)
        - ``surfaces``: list of per-surface dictionaries with ``surface_id``,
          ``face_ids``, ``n_faces``, and ``bbox``
    """
    mesh = _as_trimesh(mesh_or_geom)
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)
    n_faces = len(faces)

    if n_faces == 0:
        return {
            "n_surfaces": 0,
            "face_surface_ids": np.zeros(0, dtype=int),
            "surfaces": [],
        }

    adjacency = {i: set() for i in range(n_faces)}
    for face_a, face_b in np.asarray(mesh.face_adjacency, dtype=int):
        adjacency[int(face_a)].add(int(face_b))
        adjacency[int(face_b)].add(int(face_a))

    visited = set()
    face_surface_ids = np.zeros(n_faces, dtype=int)
    surfaces = []
    surface_id = 0
    for face_id in range(n_faces):
        if face_id in visited:
            continue
        surface_id += 1
        stack = [face_id]
        visited.add(face_id)
        component_face_ids = []
        while stack:
            current = stack.pop()
            component_face_ids.append(current)
            for nb in adjacency[current]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)

        component_face_ids = np.asarray(component_face_ids, dtype=int)
        face_surface_ids[component_face_ids] = int(surface_id)
        surfaces.append(
            {
                "surface_id": int(surface_id),
                "face_ids": component_face_ids.astype(int).tolist(),
                "n_faces": int(len(component_face_ids)),
                "bbox": _bbox_from_face_ids(vertices, faces, component_face_ids),
            }
        )

    surfaces.sort(key=lambda item: (-item["n_faces"], item["surface_id"]))
    return {
        "n_surfaces": int(surface_id),
        "face_surface_ids": face_surface_ids,
        "surfaces": surfaces,
    }


def _format_check_summary(report: Dict[str, Union[bool, int, float, list, dict]]) -> str:
    lines = ["Mesh check"]
    lines.append(f"- valid: {bool(report['valid'])}")
    lines.append(f"- faces: {int(report['n_faces'])}")
    lines.append(f"- vertices: {int(report['n_vertices'])}")
    lines.append(f"- independent surfaces: {int(report.get('n_independent_surfaces', 0))}")
    issues = list(report.get("issues", []))
    n_independent_surfaces = int(report.get("n_independent_surfaces", 0))
    if n_independent_surfaces > 1 and not any("independent surfaces" in str(issue) for issue in issues):
        issues.append(f"mesh has {n_independent_surfaces} independent surfaces")
    if issues:
        lines.append("- issues:")
        lines.extend(f"  - {issue}" for issue in issues)
    else:
        lines.append("- issues: none")
    return "\n".join(lines)


def _estimate_planar_ground_level_for_check(
    vertices: np.ndarray,
    faces: np.ndarray,
    ground_mask: np.ndarray,
    tolerance: float,
) -> float | None:
    ground_face_ids = np.flatnonzero(np.asarray(ground_mask, dtype=bool))
    if len(ground_face_ids) == 0:
        return None
    ground_vertex_ids = np.unique(faces[ground_face_ids].reshape(-1))
    if len(ground_vertex_ids) == 0:
        return None
    z = np.asarray(vertices, dtype=float)[ground_vertex_ids, 2]
    z0 = float(np.median(z))
    if float(np.max(np.abs(z - z0))) > tolerance:
        return None
    return z0


def check(mesh_or_geom, require_single_component: bool = True) -> Dict[str, Union[bool, int, float, list]]:
    """
    Validate a surface mesh for common udgeom topology problems.

    The routine combines lightweight geometric checks with spatial diagnostics.
    It is intended as the first step in a repair workflow:
    - call ``check(...)``
    - inspect ``issues`` and ``details``
    - apply one of the repair helpers if appropriate
    - call ``check(...)`` again to confirm the mesh is now clean

    Parameters
    ----------
    mesh_or_geom : trimesh.Trimesh or UDGeom
        Mesh to validate.
    require_single_component : bool, default=True
        If True, face connectivity must form a single component.

    Returns
    -------
    report : dict
        Validation report with:
        - ``valid``: overall boolean pass/fail
        - ``issues``: short human-readable issue strings
        - ``summary``: formatted multi-line summary for user-facing output
        - summary counts such as ``n_nonmanifold_edges``
        - ``details``: face ids, regions, and bounding boxes for diagnosis
    """
    mesh = _as_trimesh(mesh_or_geom)
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)

    report: Dict[str, Union[bool, int, list, dict, str]] = {
        "valid": True,
        "issues": [],
        "summary": "",
        "details": {
            "nonfinite_vertex_ids": [],
            "duplicate_face_groups": [],
            "degenerate_face_ids": [],
            "zero_area_face_ids": [],
            "zero_normal_face_ids": [],
            "bad_normal_length_face_ids": [],
            "below_ground_vertex_ids": [],
            "below_ground_bbox": None,
            "downward_ground_face_ids": [],
            "downward_ground_bbox": None,
            "nonmanifold_regions": [],
            "internal_touching_wall_regions": [],
            "unstitched_touching_regions": [],
            "connected_components": [],
            "independent_surfaces": [],
        },
        "n_vertices": int(len(vertices)),
        "n_faces": int(len(faces)),
        "n_nonfinite_vertices": 0,
        "n_nonfinite_face_centers": 0,
        "n_nonfinite_face_areas": 0,
        "n_unused_vertices": 0,
        "n_degenerate_faces": 0,
        "n_zero_area_faces": 0,
        "n_zero_normals": 0,
        "n_bad_normal_lengths": 0,
        "n_below_ground_vertices": 0,
        "n_connected_components": 0,
        "n_nonmanifold_edges": 0,
        "n_internal_touching_wall_regions": 0,
        "n_unstitched_touching_regions": 0,
        "n_duplicate_faces": 0,
        "n_downward_ground_faces": 0,
        "n_independent_surfaces": 0,
    }

    if len(faces) == 0:
        report["valid"] = False
        report["issues"].append("mesh has no faces")
        report["summary"] = _format_check_summary(report)
        return report

    # Finite values
    nonfinite_vertices = ~np.isfinite(vertices).all(axis=1)
    n_nonfinite_vertices = int(np.count_nonzero(nonfinite_vertices))
    report["n_nonfinite_vertices"] = n_nonfinite_vertices
    if n_nonfinite_vertices > 0:
        report["details"]["nonfinite_vertex_ids"] = np.flatnonzero(nonfinite_vertices).astype(int).tolist()
        report["valid"] = False
        report["issues"].append(f"mesh has {n_nonfinite_vertices} non-finite vertices")

    face_centers = np.asarray(mesh.triangles_center, dtype=float)
    nonfinite_face_centers = ~np.isfinite(face_centers).all(axis=1)
    n_nonfinite_face_centers = int(np.count_nonzero(nonfinite_face_centers))
    report["n_nonfinite_face_centers"] = n_nonfinite_face_centers
    if n_nonfinite_face_centers > 0:
        report["valid"] = False
        report["issues"].append(f"mesh has {n_nonfinite_face_centers} non-finite face centers")

    face_areas = np.asarray(mesh.area_faces, dtype=float)
    n_nonfinite_face_areas = int(np.count_nonzero(~np.isfinite(face_areas)))
    report["n_nonfinite_face_areas"] = n_nonfinite_face_areas
    if n_nonfinite_face_areas > 0:
        report["valid"] = False
        report["issues"].append(f"mesh has {n_nonfinite_face_areas} non-finite face areas")

    # Unused vertices
    used_vertices = np.unique(faces.reshape(-1))
    n_unused_vertices = int(len(vertices) - len(used_vertices))
    report["n_unused_vertices"] = n_unused_vertices
    if n_unused_vertices > 0:
        report["valid"] = False
        report["issues"].append(f"mesh has {n_unused_vertices} unused vertices")

    # Duplicate faces: same geometric triangle regardless of winding.
    canonical_faces = np.sort(faces, axis=1)
    _, inverse, counts = np.unique(canonical_faces, axis=0, return_inverse=True, return_counts=True)
    n_duplicate_faces = int(np.sum(np.clip(counts - 1, 0, None)))
    report["n_duplicate_faces"] = n_duplicate_faces
    if n_duplicate_faces > 0:
        duplicate_groups = []
        for group_id, count in enumerate(counts):
            if count > 1:
                duplicate_groups.append(np.flatnonzero(inverse == group_id).astype(int).tolist())
        report["details"]["duplicate_face_groups"] = duplicate_groups
        report["valid"] = False
        report["issues"].append(f"mesh has {n_duplicate_faces} duplicate faces")

    # Degenerate faces / zero-area faces
    n_degenerate_faces = int(np.count_nonzero(np.array([len(np.unique(face)) < 3 for face in faces], dtype=bool)))
    report["n_degenerate_faces"] = n_degenerate_faces
    if n_degenerate_faces > 0:
        report["details"]["degenerate_face_ids"] = np.flatnonzero(
            np.array([len(np.unique(face)) < 3 for face in faces], dtype=bool)
        ).astype(int).tolist()
        report["valid"] = False
        report["issues"].append(f"mesh has {n_degenerate_faces} degenerate faces")

    finite_area_mask = np.isfinite(face_areas)
    zero_area_mask = finite_area_mask & (face_areas <= 1.0e-12)
    n_zero_area_faces = int(np.count_nonzero(zero_area_mask))
    report["n_zero_area_faces"] = n_zero_area_faces
    if n_zero_area_faces > 0:
        report["details"]["zero_area_face_ids"] = np.flatnonzero(zero_area_mask).astype(int).tolist()
        report["valid"] = False
        report["issues"].append(f"mesh has {n_zero_area_faces} zero-area faces")

    # Normals
    face_normals = np.asarray(mesh.face_normals, dtype=float)
    normal_lengths = np.linalg.norm(face_normals, axis=1)
    zero_normal_mask = np.isfinite(normal_lengths) & (normal_lengths <= 1.0e-12)
    n_zero_normals = int(np.count_nonzero(zero_normal_mask))
    report["n_zero_normals"] = n_zero_normals
    if n_zero_normals > 0:
        report["details"]["zero_normal_face_ids"] = np.flatnonzero(zero_normal_mask).astype(int).tolist()
        report["valid"] = False
        report["issues"].append(f"mesh has {n_zero_normals} zero normals")

    bad_normal_length_mask = np.isfinite(normal_lengths) & ~zero_normal_mask & (np.abs(normal_lengths - 1.0) > 1.0e-6)
    n_bad_normal_lengths = int(np.count_nonzero(bad_normal_length_mask))
    report["n_bad_normal_lengths"] = n_bad_normal_lengths
    if n_bad_normal_lengths > 0:
        report["details"]["bad_normal_length_face_ids"] = np.flatnonzero(bad_normal_length_mask).astype(int).tolist()
        report["valid"] = False
        report["issues"].append(f"mesh has {n_bad_normal_lengths} non-unit face normals")

    # Orientation sanity: downward-pointing ground faces are ignored by IBM and are almost always accidental.
    ground_mask = identify_ground_faces(mesh)
    planar_ground_z = _estimate_planar_ground_level_for_check(
        vertices,
        faces,
        ground_mask,
        tolerance=1.0e-6,
    )
    if planar_ground_z is not None:
        ground_vertex_ids = np.unique(faces[np.flatnonzero(ground_mask)].reshape(-1))
        non_ground_vertex_mask = np.ones(len(vertices), dtype=bool)
        non_ground_vertex_mask[ground_vertex_ids] = False
        non_ground_vertex_ids = np.flatnonzero(non_ground_vertex_mask)
        if len(non_ground_vertex_ids) > 0:
            ground_xy = vertices[ground_vertex_ids][:, :2]
            non_ground_xy = vertices[non_ground_vertex_ids][:, :2]
            ground_min = np.min(ground_xy, axis=0)
            ground_max = np.max(ground_xy, axis=0)
            non_ground_min = np.min(non_ground_xy, axis=0)
            non_ground_max = np.max(non_ground_xy, axis=0)
            extends_beyond_non_ground = bool(
                np.any(ground_min < non_ground_min - 1.0e-9)
                or np.any(ground_max > non_ground_max + 1.0e-9)
            )
            if extends_beyond_non_ground:
                below_ground_mask = non_ground_vertex_mask & np.isfinite(vertices).all(axis=1) & (vertices[:, 2] < planar_ground_z - 1.0e-9)
                n_below_ground_vertices = int(np.count_nonzero(below_ground_mask))
                report["n_below_ground_vertices"] = n_below_ground_vertices
                if n_below_ground_vertices > 0:
                    below_vertex_ids = np.flatnonzero(below_ground_mask).astype(int)
                    report["details"]["below_ground_vertex_ids"] = below_vertex_ids.tolist()
                    points = vertices[below_vertex_ids]
                    report["details"]["below_ground_bbox"] = np.array([points.min(axis=0), points.max(axis=0)], dtype=float)
                    report["valid"] = False
                    report["issues"].append(f"mesh has {n_below_ground_vertices} non-ground vertices below planar ground")

    downward_ground_mask = ground_mask & np.isfinite(face_normals).all(axis=1) & (face_normals[:, 2] < -1.0e-6)
    n_downward_ground_faces = int(np.count_nonzero(downward_ground_mask))
    report["n_downward_ground_faces"] = n_downward_ground_faces
    if n_downward_ground_faces > 0:
        downward_face_ids = np.flatnonzero(downward_ground_mask).astype(int)
        report["details"]["downward_ground_face_ids"] = downward_face_ids.tolist()
        report["details"]["downward_ground_bbox"] = _bbox_from_face_ids(vertices, faces, downward_face_ids)
        report["valid"] = False
        report["issues"].append(f"mesh has {n_downward_ground_faces} downward-facing ground faces")

    # Non-manifold edges
    edge_counts = {}
    for face in faces:
        a, b, c = map(int, face)
        for edge in (tuple(sorted((a, b))), tuple(sorted((b, c))), tuple(sorted((c, a)))):
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    n_nonmanifold_edges = int(sum(1 for count in edge_counts.values() if count > 2))
    report["n_nonmanifold_edges"] = n_nonmanifold_edges
    if n_nonmanifold_edges > 0:
        report["details"]["nonmanifold_regions"] = find_nonmanifold_regions(mesh)
        report["valid"] = False
        report["issues"].append(f"mesh has {n_nonmanifold_edges} non-manifold edges")

    internal_touching_wall_regions = find_internal_touching_wall_regions(mesh)
    n_internal_touching_wall_regions = int(len(internal_touching_wall_regions))
    report["n_internal_touching_wall_regions"] = n_internal_touching_wall_regions
    if n_internal_touching_wall_regions > 0:
        report["details"]["internal_touching_wall_regions"] = internal_touching_wall_regions
        report["valid"] = False
        report["issues"].append(
            f"mesh has {n_internal_touching_wall_regions} internal touching wall regions"
        )

    unstitched_touching_regions = find_unstitched_touching_regions(mesh)
    n_unstitched_touching_regions = int(len(unstitched_touching_regions))
    report["n_unstitched_touching_regions"] = n_unstitched_touching_regions
    if n_unstitched_touching_regions > 0:
        report["details"]["unstitched_touching_regions"] = unstitched_touching_regions
        report["valid"] = False
        report["issues"].append(
            f"mesh has {n_unstitched_touching_regions} unstitched touching regions"
        )

    # Face connectivity / independent surfaces
    surface_result = calculate_independent_surfaces(mesh)
    components = int(surface_result["n_surfaces"])
    component_summaries = []
    for item in surface_result["surfaces"]:
        component_summaries.append(
            {
                "surface_id": int(item["surface_id"]),
                "face_ids_preview": item["face_ids"][:20],
                "n_faces": int(item["n_faces"]),
                "bbox": item["bbox"],
            }
        )

    report["n_connected_components"] = int(components)
    report["n_independent_surfaces"] = int(components)
    report["details"]["connected_components"] = component_summaries
    report["details"]["independent_surfaces"] = surface_result["surfaces"]
    if require_single_component and components != 1:
        report["valid"] = False
        report["issues"].append(f"mesh has {components} disconnected face components")

    report["summary"] = _format_check_summary(report)
    return report


__all__ = [
    "check",
    "calculate_independent_surfaces",
    "identify_ground_faces",
    "find_nonmanifold_regions",
    "find_internal_touching_wall_regions",
    "find_unstitched_touching_regions",
]
