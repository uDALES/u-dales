"""
uniform_facets - make the facet (triangle) sizes of an STL geometry uniform.

uDALES treats every STL triangle as one facet with a single surface temperature
and a single set of fluxes. Real-world STL files mix triangles of very different
sizes: a few huge ground or roof triangles next to clusters of sliver triangles
left behind by CAD boolean operations. Both ends are bad for the surface energy
balance: a huge facet spreads one temperature over a whole wall, while a sliver
gets all-or-nothing sun and hardly any view-factor neighbours.

This module rebuilds the mesh so that facet areas fall within a band around a
target area while preserving the geometry, the facet types and the conformity
of the mesh (no T-junctions, no unused vertices):

* faces are grouped into coplanar, edge-connected patches of equal facet type;
* a patch that contains faces below the band is retriangulated as a whole,
  which merges the tiny faces into their coplanar neighbours;
* faces above the band are split, as part of the retriangulation of their patch;
* faces inside the band are copied unchanged, so a mesh that is already uniform
  passes through untouched (the operation is idempotent);
* boundary vertices that exist only because of the faces being merged (the
  corners of slivers along a wall base, say) are dropped where the straightened
  boundary stays within the coplanarity tolerance, and boundary edges between
  two retriangulated patches are subdivided to the target edge length, both in
  one shared edge table so neighbouring patches stay conforming; an edge shared
  with an untouched face is never changed.

Patch interiors are triangulated with Jonathan Shewchuk's Triangle (constrained
Delaunay refinement with maximum-area and minimum-angle constraints, and no new
vertices on the constraint segments).

uDALES (https://github.com/uDALES/u-dales)

Copyright (C) 2016-2026 the uDALES Team.
Licensed under GNU General Public License v3.0
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import logging
import math
import multiprocessing
import time
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from exceptions import DependencyError

try:
    import trimesh
except ImportError as exc:  # pragma: no cover - import guard
    raise DependencyError("trimesh is required for uniform_facets; install with `pip install trimesh`.") from exc

import triangle as triangle_lib
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely import (
    STRtree,
    contains_xy as shapely_contains_xy,
    is_simple as shapely_is_simple,
    multilinestrings as shapely_multilinestrings,
    points as shapely_points,
    polygons as shapely_polygons,
    unary_union as shapely_unary_union,
)

logger = logging.getLogger(__name__)

# Area of the equilateral triangle with unit edge length.
_EQUILATERAL_AREA = math.sqrt(3.0) / 4.0

# Triangle refines a region until every triangle is below the maximum-area
# constraint; the mean area of the result is roughly this fraction of the
# constraint. It converts the target (mean) facet area into the constraint
# handed to Triangle and is used to predict the facet count before remeshing.
_MEAN_AREA_FRACTION = 0.6


class FacetCountError(RuntimeError):
    """Raised when uniformising would produce more facets than ``max_faces``."""


def target_edge_from_area(area: float) -> float:
    """Edge length of the equilateral triangle with the given area."""
    return math.sqrt(float(area) / _EQUILATERAL_AREA)


def target_area_from_edge(edge: float) -> float:
    """Area of the equilateral triangle with the given edge length."""
    return _EQUILATERAL_AREA * float(edge) ** 2


@dataclass
class UniformFacetsResult:
    """Outcome of :func:`uniformize_facets`.

    Attributes
    ----------
    mesh : trimesh.Trimesh
        The uniformised mesh (welded, no unused vertices).
    face_types : ndarray of int or None
        Facet type of every output face (inherited from the face it replaces);
        ``None`` when no types were given.
    face_map : ndarray of int
        For every output face, the index of the input face it came from
        (the input face containing its centroid for retriangulated patches).
    report : dict
        Counts, areas and settings; see :meth:`summary`.
    """

    mesh: "trimesh.Trimesh"
    face_types: Optional[np.ndarray]
    face_map: np.ndarray
    report: Dict[str, Any]

    def summary(self) -> str:
        r = self.report
        lines = [
            "Facet uniformisation",
            f"  target edge / area      : {r['target_edge']:.3f} m / {r['target_area']:.3f} m^2"
            + ("  (automatic: total area / input faces)" if r["target_auto"] else ""),
            f"  accepted area band      : [{r['area_low']:.4g}, {r['area_high']:.4g}] m^2  (band {r['size_band']:g})",
            f"  split / merge allowed   : {r['allow_split']} / {r['allow_merge']}",
            f"  faces   in -> out       : {r['n_faces_in']} -> {r['n_faces_out']}  ({r['n_faces_out'] / max(r['n_faces_in'], 1):.3f}x)",
            f"  vertices in -> out      : {r['n_vertices_in']} -> {r['n_vertices_out']}",
            f"  total area in -> out    : {r['area_in']:.6g} -> {r['area_out']:.6g} m^2",
            f"  coplanar patches        : {r['n_patches']} ({r['n_patches_remeshed']} retriangulated, "
            f"{r['n_patches_noncoplanar']} dissolved as not coplanar, {r['n_patches_kept']} kept unchanged)",
            f"  faces merged / split    : {r['n_faces_merged']} / {r['n_faces_split']}",
            f"  boundary chains split   : {r['n_edges_split']}  (collinear boundary vertices removed: {r['n_boundary_vertices_removed']})",
            f"  small isolated faces    : {r['n_small_isolated']}  (below the band, no coplanar neighbour; kept as is)",
            f"  degenerate faces removed: {r['n_degenerate_removed']}",
            f"  out of band (out)       : {r['n_out_below_band']} below, {r['n_out_above_band']} above",
        ]
        pin = r["area_percentiles_in"]
        pout = r["area_percentiles_out"]
        lines.append(f"  boundary points of retriangulated patches: {r['n_boundary_points_remeshed']}"
                     "  (segment length percentiles [1 10 50 90 99] m: "
                     + " ".join(f"{v:.3g}" for v in r["boundary_segment_percentiles"]) + ")")
        lines.append("  area percentiles [1 10 50 90 99] m^2")
        lines.append("    in : " + " ".join(f"{v:.4g}" for v in pin))
        lines.append("    out: " + " ".join(f"{v:.4g}" for v in pout))
        lines.append(f"  elapsed                 : {r['elapsed_s']:.2f} s")
        return "\n".join(lines)


# -----------------------------------------------------------------------------#
# Helpers
# -----------------------------------------------------------------------------#

def _remove_unused(vertices: np.ndarray, faces: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    used = np.unique(faces.reshape(-1))
    new_index = -np.ones(len(vertices), dtype=int)
    new_index[used] = np.arange(len(used))
    return vertices[used], new_index[faces]


def _face_geometry(vertices: np.ndarray, faces: np.ndarray):
    tri = vertices[faces]
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    double_area = np.linalg.norm(cross, axis=1)
    areas = 0.5 * double_area
    normals = np.zeros_like(cross)
    ok = double_area > 0.0
    normals[ok] = cross[ok] / double_area[ok, None]
    centroids = tri.mean(axis=1)
    return areas, normals, centroids


def _edge_table(faces: np.ndarray):
    """Unique undirected edges of a face array.

    Returns ``(edges, face_edge_ids, edge_face_ptr, edge_face_idx)`` where
    ``edges`` has shape (E, 2) with sorted vertex ids, ``face_edge_ids`` (n, 3)
    gives the edge id of each face edge, and ``edge_face_idx[edge_face_ptr[e]:
    edge_face_ptr[e + 1]]`` lists the faces incident to edge ``e`` (CSR).
    """
    raw = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=0)
    raw.sort(axis=1)
    edges, inverse = np.unique(raw, axis=0, return_inverse=True)
    inverse = np.asarray(inverse).reshape(-1)
    n = len(faces)
    face_edge_ids = np.column_stack([inverse[:n], inverse[n:2 * n], inverse[2 * n:]])
    face_of_entry = np.concatenate([np.arange(n), np.arange(n), np.arange(n)])
    order = np.argsort(inverse, kind="stable")
    edge_face_idx = face_of_entry[order]
    counts = np.bincount(inverse, minlength=len(edges))
    edge_face_ptr = np.concatenate([[0], np.cumsum(counts)])
    return edges, face_edge_ids, edge_face_ptr, edge_face_idx


def _adjacent_face_pairs(edge_face_ptr: np.ndarray, edge_face_idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """All (i, j) face pairs sharing an edge, including non-manifold edges."""
    counts = np.diff(edge_face_ptr)
    pairs_i: List[np.ndarray] = []
    pairs_j: List[np.ndarray] = []
    two = np.flatnonzero(counts == 2)
    if two.size:
        start = edge_face_ptr[two]
        pairs_i.append(edge_face_idx[start])
        pairs_j.append(edge_face_idx[start + 1])
    for e in np.flatnonzero(counts > 2):
        members = edge_face_idx[edge_face_ptr[e]:edge_face_ptr[e + 1]]
        a, b = np.triu_indices(len(members), k=1)
        pairs_i.append(members[a])
        pairs_j.append(members[b])
    if not pairs_i:
        return np.empty(0, dtype=int), np.empty(0, dtype=int)
    return np.concatenate(pairs_i), np.concatenate(pairs_j)


def _plane_bases(normals: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Orthonormal in-plane bases (u, v) with u x v = normal, one per row of ``normals``."""
    normals = np.atleast_2d(np.asarray(normals, dtype=float))
    helper = np.zeros_like(normals)
    helper[np.arange(len(normals)), np.argmin(np.abs(normals), axis=1)] = 1.0
    u = np.cross(normals, helper)
    u /= np.sqrt(np.einsum("ij,ij->i", u, u))[:, None]
    v = np.cross(normals, u)
    return u, v


def _prescreen_patch(
    kept_points: np.ndarray,
    points: np.ndarray,
    segments: np.ndarray,
    patch_area: float,
    plane_tol: float,
) -> Optional[str]:
    """Reason why a patch must not be handed to Triangle, or None when it is fine.

    ``kept_points`` are the boundary points that survive the removal of
    collinear vertices; ``points``/``segments`` the raw boundary. Triangle
    segfaults on (numerically) collinear input and would have to add vertices
    on the constraint segments when they cross or overlap, which conformity
    forbids; such patches are kept as they are.
    """
    if len(kept_points) < 3:
        return "fewer than three boundary points"
    if not np.isfinite(points).all():
        return "non-finite boundary coordinates"
    centred = kept_points - kept_points.mean(axis=0)
    sv = np.linalg.svd(centred, compute_uv=False)
    extent = float(np.ptp(kept_points, axis=0).max())
    # Triangle's region flood fill is unreliable for aspect ratios in the
    # thousands, hence the margin on the width (ten coplanarity tolerances).
    if sv[0] <= 0.0 or sv[1] <= 1.0e-9 * sv[0] or extent <= 0.0 or 2.0 * patch_area / extent <= 10.0 * plane_tol:
        return "degenerate (needle) patch"
    if segments.size == 0:
        return "no constraint segments"
    if np.any(segments[:, 0] == segments[:, 1]):
        return "zero-length constraint segment"
    if np.any(np.linalg.norm(points[segments[:, 0]] - points[segments[:, 1]], axis=1) <= 0.0):
        return "constraint segment between coincident points"
    _, counts = np.unique(np.round(kept_points, 9), axis=0, return_counts=True)
    if np.any(counts > 1):
        return "duplicated boundary points"
    if not shapely_is_simple(shapely_multilinestrings(points[segments])):
        return "self-intersecting boundary"
    return None


def _triangulate_region(pslg: Dict[str, np.ndarray], opts: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    result = triangle_lib.triangulate(pslg, opts)
    tris = np.asarray(result.get("triangles", np.zeros((0, 3), dtype=int)), dtype=int)
    attr = np.asarray(result.get("triangle_attributes", np.zeros((len(tris), 1)))).reshape(-1)
    tris = tris[attr > 0.5]
    rverts = np.asarray(result["vertices"], dtype=float)
    p2 = rverts[tris]
    signed = (p2[:, 1, 0] - p2[:, 0, 0]) * (p2[:, 2, 1] - p2[:, 0, 1]) - (p2[:, 2, 0] - p2[:, 0, 0]) * (p2[:, 1, 1] - p2[:, 0, 1])
    # Triangle returns counter-clockwise triangles in (u, v), i.e. the patch
    # normal for a right-handed (u, v, n) basis; make sure of it.
    flip = signed < 0.0
    if flip.any():
        tris = tris.copy()
        tris[flip] = tris[flip][:, [0, 2, 1]]
    return tris, rverts, np.abs(signed)


def _remesh_patch(
    patch: int,
    points: np.ndarray,
    segments: np.ndarray,
    regions: np.ndarray,
    orig2: np.ndarray,
    patch_area: float,
    plane_tol: float,
    min_angle: float,
) -> Tuple[np.ndarray, np.ndarray, bool]:
    """Constrained Delaunay refinement of one planar patch.

    Returns ``(triangles, vertices2d, overlap)``; the first ``len(points)``
    vertices are the input points in order. ``overlap`` is True when the input
    faces of the patch did not cover its polygon exactly once (overlapping or
    tilted slivers), in which case the output covers the polygon once.
    """
    pslg = {"vertices": points, "segments": segments, "regions": regions}
    # First pass without refinement to find the holes: every triangle outside
    # the marked regions (in a hole or outside the boundary) yields a hole point.
    first = triangle_lib.triangulate(pslg, "pAQ")
    tris0 = np.asarray(first.get("triangles", np.zeros((0, 3), dtype=int)), dtype=int)
    attr0 = np.asarray(first.get("triangle_attributes", np.zeros((len(tris0), 1)))).reshape(-1)
    if tris0.size and np.any(attr0 < 0.5):
        v0 = np.asarray(first["vertices"], dtype=float)
        pslg["holes"] = v0[tris0[attr0 < 0.5]].mean(axis=1)

    quality = f"q{float(min_angle):g}" if min_angle > 0.0 else ""
    tris, rverts, areas2 = _triangulate_region(pslg, f"p{quality}aAYYQ")
    if tris.size == 0:
        # Region marking failed (very thin region): fall back to the plain
        # constrained triangulation of the boundary points and keep the
        # triangles whose centroids lie in the patch polygon. No vertex is
        # added, so neighbours stay conforming.
        result = triangle_lib.triangulate({"vertices": points, "segments": segments}, "pQ")
        tris_all = np.asarray(result.get("triangles", np.zeros((0, 3), dtype=int)), dtype=int)
        rverts = np.asarray(result["vertices"], dtype=float)
        polygon = shapely_unary_union(shapely_polygons(orig2))
        inside = shapely_contains_xy(polygon, rverts[tris_all].mean(axis=1))
        tris = tris_all[inside]
        p2 = rverts[tris]
        signed = (p2[:, 1, 0] - p2[:, 0, 0]) * (p2[:, 2, 1] - p2[:, 0, 1]) - (p2[:, 2, 0] - p2[:, 0, 0]) * (p2[:, 1, 1] - p2[:, 0, 1])
        flip = signed < 0.0
        if flip.any():
            tris = tris.copy()
            tris[flip] = tris[flip][:, [0, 2, 1]]
        areas2 = np.abs(signed)
    new_area = 0.5 * float(areas2.sum())
    # Boundary vertices carry the coordinate noise of the STL (single precision
    # for most files); the area of the polygon is uncertain by about that
    # noise times its perimeter.
    perimeter = float(np.linalg.norm(points[segments[:, 0]] - points[segments[:, 1]], axis=1).sum())
    tol = max(1.0e-6 * max(1.0, patch_area), plane_tol * perimeter)
    if tris.size and abs(new_area - patch_area) <= tol:
        return tris, rverts, False
    # The sum of the input face areas is not the polygon area: the faces overlap
    # or are tilted within the coplanarity tolerance. Compare against the
    # polygon itself.
    polygon_area = float(shapely_unary_union(shapely_polygons(orig2)).area)
    if tris.size and abs(new_area - polygon_area) <= tol:
        return tris, rverts, True
    raise RuntimeError(
        f"retriangulation of patch {patch} failed: the new triangles cover {new_area:.6g} m^2, "
        f"the input faces {patch_area:.6g} m^2 and their union {polygon_area:.6g} m^2 "
        f"(patch bounding box in plane coordinates {points.min(axis=0)} .. {points.max(axis=0)})"
    )


def _percentiles(values: np.ndarray) -> List[float]:
    if values.size == 0:
        return [float("nan")] * 5
    return [float(x) for x in np.percentile(values, [1, 10, 50, 90, 99])]


# -----------------------------------------------------------------------------#
# Per-patch tasks
# -----------------------------------------------------------------------------#
# The pre-screen and the retriangulation treat every patch on its own, given
# the read-only tables that uniformize_facets builds once (vertices, faces,
# edges, the per-patch CSR tables, the plane bases). The tables are published
# to the tasks through this module-level dictionary: set directly on the serial
# path, by the pool initializer in worker processes (inherited for free with the
# fork start method). A task takes a patch id and returns small arrays, and the
# results are consumed in patch order, so the output does not depend on the
# number of workers.
_STATE: Dict[str, Any] = {}

# Each worker process must have at least this many patches to do: starting a
# process costs about as much as a few dozen patches, so beyond this point more
# workers only add start-up time (on a 380k-face city, 32 workers run the
# per-patch phases as fast as 128).
_MIN_PATCHES_PER_WORKER = 512


def _set_state(state: Dict[str, Any]) -> None:
    _STATE.clear()
    _STATE.update(state)


def _map_patches(task: Callable[[int], Any], ids, state: Dict[str, Any], workers: int) -> List[Any]:
    """``[task(p) for p in ids]`` with ``state`` published to the task, run in
    the calling process or in up to ``workers`` processes; the results come
    back in the order of ``ids`` either way."""
    ids = [int(p) for p in ids]
    workers = max(1, min(int(workers), len(ids) // _MIN_PATCHES_PER_WORKER))
    if workers <= 1:
        _set_state(state)
        try:
            return [task(p) for p in ids]
        finally:
            _STATE.clear()
    methods = multiprocessing.get_all_start_methods()
    ctx = multiprocessing.get_context("fork") if "fork" in methods else multiprocessing.get_context()
    chunksize = max(1, len(ids) // (16 * workers))
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx, initializer=_set_state, initargs=(state,)) as pool:
        return list(pool.map(task, ids, chunksize=chunksize))


def _prescreen_task(p: int) -> Optional[str]:
    """Pre-screen one scheduled patch; see :func:`_prescreen_patch`."""
    s = _STATE
    edges = s["edges"]
    pe = s["pce_idx"][s["pce_ptr"][p]:s["pce_ptr"][p + 1]]
    pv = s["pv_idx"][s["pv_ptr"][p]:s["pv_ptr"][p + 1]]
    point_ids = np.unique(np.concatenate([edges[pe].reshape(-1), pv[s["vertex_shared"][pv]]]))
    rel = s["vertices"][point_ids] - s["patch_centre"][p]
    pts2 = np.column_stack([rel @ s["patch_u"][p], rel @ s["patch_v"][p]])
    seg2 = np.searchsorted(point_ids, edges[pe]) if pe.size else np.empty((0, 2), dtype=int)
    return _prescreen_patch(pts2[~s["removable"][point_ids]], pts2, seg2, float(s["patch_area"][p]), s["plane_tol"])


def _remesh_task(p: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray], bool]:
    """Retriangulate one scheduled patch.

    Returns ``(point_ids, triangles, steiner, face_map, segment_lengths,
    overlap)``: the global ids of the boundary points; the new triangles as
    indices into the boundary points followed by the new (Steiner) vertices;
    those vertices in 3-D; the index within the patch of the input face every
    triangle comes from; the lengths of the constraint segments (None when
    there are none); and the overlap flag of :func:`_remesh_patch`.
    """
    s = _STATE
    pf = s["face_order"][s["face_ptr"][p]:s["face_ptr"][p + 1]]
    u, v, centre = s["patch_u"][p], s["patch_v"][p], s["patch_centre"][p]

    # Constraint chains of this patch, with their split vertices.
    pe = s["pce_idx"][s["pce_ptr"][p]:s["pce_ptr"][p + 1]]
    chain_ends = s["chain_ends"]
    chain_split_vertices = s["chain_split_vertices"]
    seg_list: List[np.ndarray] = []
    for c in np.unique(s["chain_of_edge"][pe]):
        a, b = chain_ends[c]
        mids = chain_split_vertices.get(int(c))
        if mids is None:
            seg_list.append(np.array([[a, b]]))
        else:
            chain = np.concatenate([[a], mids, [b]])
            seg_list.append(np.column_stack([chain[:-1], chain[1:]]))
    segments_global = np.concatenate(seg_list, axis=0) if seg_list else np.empty((0, 2), dtype=int)
    pv = s["pv_idx"][s["pv_ptr"][p]:s["pv_ptr"][p + 1]]
    fixed = pv[s["vertex_shared"][pv]]
    point_ids = np.unique(np.concatenate([segments_global.reshape(-1), fixed]))
    segments_local = np.searchsorted(point_ids, segments_global) if segments_global.size else np.empty((0, 2), dtype=int)
    rel = s["base_vertices"][point_ids] - centre
    pts2 = np.column_stack([rel @ u, rel @ v])

    # Original faces in plane coordinates; region markers are the 2-D
    # centroids of the largest face of every sub-region.
    faces = s["faces"]
    areas = s["areas"]
    region_of_face = s["region_of_face"]
    orig_rel = s["vertices"][faces[pf]] - centre
    orig2 = np.stack([orig_rel @ u, orig_rel @ v], axis=-1)       # (m, 3, 2)
    regions = []
    for r in np.unique(region_of_face[pf]):
        local = np.flatnonzero(region_of_face[pf] == r)
        f0 = local[np.argmax(areas[pf[local]])]
        cx, cy = orig2[f0].mean(axis=0)
        regions.append([float(cx), float(cy), 1.0, float(s["max_area"])])
    regions_arr = np.asarray(regions, dtype=float)

    lengths = None
    if segments_local.size:
        lengths = np.linalg.norm(pts2[segments_local[:, 0]] - pts2[segments_local[:, 1]], axis=1)
    tris, rverts, overlap = _remesh_patch(int(p), pts2, segments_local, regions_arr, orig2,
                                          float(s["patch_area"][p]), s["plane_tol"], s["min_angle"])

    # Lift new (Steiner) vertices back to 3-D on the patch plane.
    n_pts = len(point_ids)
    steiner = centre + rverts[n_pts:, 0:1] * u + rverts[n_pts:, 1:2] * v

    # Face map: input face containing each new centroid (2-D search).
    tree = STRtree(shapely_polygons(orig2))
    cents = shapely_points(rverts[tris].mean(axis=1))
    q_in, q_tree = tree.query(cents, predicate="within")
    fmap_local = np.full(len(tris), -1, dtype=int)
    fmap_local[q_in] = q_tree
    missing = fmap_local < 0
    if missing.any():
        fmap_local[missing] = tree.nearest(cents[missing])
    return point_ids, tris, steiner, fmap_local, lengths, overlap


# -----------------------------------------------------------------------------#
# Main routine
# -----------------------------------------------------------------------------#

def uniformize_facets(
    mesh: "trimesh.Trimesh",
    *,
    face_types: Optional[np.ndarray] = None,
    target_edge: Optional[float] = None,
    size_band: float = 4.0,
    allow_split: bool = True,
    allow_merge: bool = True,
    max_faces: Optional[int] = None,
    plane_tol: Optional[float] = None,
    min_angle: float = 0.0,
    weld_digits: int = 6,
    workers: Optional[int] = None,
    _debug: Optional[Dict[str, Any]] = None,
) -> UniformFacetsResult:
    """Split large facets and merge tiny coplanar facets of a triangle mesh.

    Parameters
    ----------
    mesh : trimesh.Trimesh
        Input geometry. It is not modified.
    face_types : array of int, optional
        Facet type per input face. Faces of different types are never merged,
        and every output face inherits the type of the face it replaces.
    target_edge : float, optional
        Target facet edge length [m]; the target facet area is that of the
        equilateral triangle with this edge. ``None`` (or 0) selects the
        fixed-count default: the mean area of the input faces, so that the
        number of facets stays about the same while their sizes become
        uniform.
    size_band : float, default 4.0
        Faces with an area between ``target_area / size_band`` and
        ``target_area * size_band`` are accepted as they are.
    allow_split : bool, default True
        Allow faces above the band to be split.
    allow_merge : bool, default True
        Allow faces below the band to be merged with coplanar neighbours of
        the same type. When False, faces are never grouped into patches and
        only splitting of individual faces takes place.
    max_faces : int, optional
        Refuse (raise :class:`FacetCountError`) when the output would have more
        faces than this. Defaults to four times the input face count.
    plane_tol : float, optional
        Distance [m] within which vertices count as lying in a neighbour's
        plane. Default ``4e-6 * max(1, largest |coordinate|)``, i.e. a few
        single-precision ulps of the STL coordinates.
    min_angle : float, default 0.0
        Minimum angle [deg] requested from Triangle for new triangles; 0 means
        no angle constraint (refinement by area alone), which gives the fewest
        facets. A value such as 20 gives better-shaped triangles next to short
        boundary edges at the price of more of them.
    weld_digits : int, default 6
        Decimal digits used to weld coincident input vertices.
    workers : int, optional
        Maximum number of worker processes for the per-patch work
        (pre-screening and retriangulating the patches); fewer are started
        when there are not enough patches to keep them busy. ``None`` or 1
        does everything in the calling process. The result does not depend
        on this setting.

    Returns
    -------
    UniformFacetsResult

    Raises
    ------
    FacetCountError
        When the predicted or actual output face count exceeds ``max_faces``.
    ValueError
        For invalid settings or an empty mesh.
    """
    t0 = time.perf_counter()
    if size_band < 1.0:
        raise ValueError("size_band must be >= 1")
    if not (allow_split or allow_merge):
        raise ValueError("at least one of allow_split and allow_merge must be true")
    workers = 1 if workers is None else int(workers)
    if workers < 1:
        raise ValueError("workers must be >= 1")

    # --- 0. weld, drop degenerate faces -------------------------------------
    work = trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices, dtype=float).copy(),
        faces=np.asarray(mesh.faces, dtype=int).copy(),
        process=False,
    )
    n_in = len(work.faces)
    if n_in == 0:
        raise ValueError("mesh has no faces")
    if face_types is not None:
        face_types = np.asarray(face_types, dtype=int).reshape(-1)
        if len(face_types) != n_in:
            raise ValueError(f"face_types has {len(face_types)} entries but the mesh has {n_in} faces")
    work.merge_vertices(digits_vertex=int(weld_digits))
    vertices = np.asarray(work.vertices, dtype=float)
    faces_all = np.asarray(work.faces, dtype=int)
    areas_all, _, _ = _face_geometry(vertices, faces_all)
    nondegenerate = (areas_all > 1.0e-12) & (faces_all[:, 0] != faces_all[:, 1]) \
        & (faces_all[:, 1] != faces_all[:, 2]) & (faces_all[:, 0] != faces_all[:, 2])
    if len(faces_all) != n_in:  # merge_vertices never drops faces, but be safe
        raise RuntimeError("unexpected change of face count while welding")
    keep_index = np.flatnonzero(nondegenerate)          # -> input face index
    n_degenerate = int(n_in - keep_index.size)
    faces = faces_all[keep_index]
    types = face_types[keep_index] if face_types is not None else np.ones(len(faces), dtype=int)
    areas, normals, centroids = _face_geometry(vertices, faces)
    n = len(faces)
    if n == 0:
        raise ValueError("mesh has no non-degenerate faces")
    total_area = float(areas.sum())

    if plane_tol is None:
        plane_tol = 4.0e-6 * max(1.0, float(np.abs(vertices).max()))

    # --- 1. target -----------------------------------------------------------
    target_auto = target_edge is None or float(target_edge) <= 0.0
    if target_auto:
        target_area = total_area / n
        target_edge = target_edge_from_area(target_area)
    else:
        target_edge = float(target_edge)
        target_area = target_area_from_edge(target_edge)
    area_low = target_area / size_band
    area_high = target_area * size_band
    max_area = target_area / _MEAN_AREA_FRACTION
    if max_faces is None:
        max_faces = 4 * n_in
    max_faces = int(max_faces)

    too_small = areas < area_low
    too_large = areas > area_high

    # --- 2. edges and coplanar patches --------------------------------------
    edges, face_edge_ids, edge_ptr, edge_idx = _edge_table(faces)
    n_edges = len(edges)

    if allow_merge:
        pi, pj = _adjacent_face_pairs(edge_ptr, edge_idx)
        if pi.size:
            same_type = types[pi] == types[pj]
            same_side = np.einsum("ij,ij->i", normals[pi], normals[pj]) > 0.0
            larger_is_i = areas[pi] >= areas[pj]
            ref = np.where(larger_is_i, pi, pj)
            oth = np.where(larger_is_i, pj, pi)
            dist = np.abs(np.einsum("ijk,ik->ij", vertices[faces[oth]] - centroids[ref][:, None, :], normals[ref]))
            coplanar = dist.max(axis=1) <= plane_tol
            ok = same_type & same_side & coplanar
            pi, pj = pi[ok], pj[ok]
        graph = coo_matrix((np.ones(pi.size), (pi, pj)), shape=(n, n))
        n_patches, patch_of_face = connected_components(graph, directed=False)
    else:
        n_patches = n
        patch_of_face = np.arange(n)

    # Patch plane (area-weighted) and verification against drift.
    patch_area = np.bincount(patch_of_face, weights=areas, minlength=n_patches)
    patch_normal = np.zeros((n_patches, 3))
    patch_centre = np.zeros((n_patches, 3))
    for k in range(3):
        patch_normal[:, k] = np.bincount(patch_of_face, weights=areas * normals[:, k], minlength=n_patches)
        patch_centre[:, k] = np.bincount(patch_of_face, weights=areas * centroids[:, k], minlength=n_patches)
    patch_normal /= np.maximum(np.linalg.norm(patch_normal, axis=1), 1e-300)[:, None]
    patch_centre /= patch_area[:, None]
    patch_count = np.bincount(patch_of_face, minlength=n_patches)

    n_noncoplanar = 0
    if allow_merge and n_patches < n:
        dev = np.abs(np.einsum("ijk,ik->ij", vertices[faces] - patch_centre[patch_of_face][:, None, :],
                               patch_normal[patch_of_face])).max(axis=1)
        bad_patch = np.zeros(n_patches, dtype=bool)
        np.maximum.at(bad_patch, patch_of_face, dev > 10.0 * plane_tol)
        bad_patch &= patch_count > 1
        if bad_patch.any():
            # Dissolve drifting patches into single-face patches.
            n_noncoplanar = int(bad_patch.sum())
            dissolve = bad_patch[patch_of_face]
            new_ids = np.arange(int(dissolve.sum())) + n_patches
            patch_of_face = patch_of_face.copy()
            patch_of_face[dissolve] = new_ids
            _, patch_of_face = np.unique(patch_of_face, return_inverse=True)
            patch_of_face = np.asarray(patch_of_face).reshape(-1)
            n_patches = int(patch_of_face.max()) + 1
            patch_area = np.bincount(patch_of_face, weights=areas, minlength=n_patches)
            patch_normal = np.zeros((n_patches, 3))
            patch_centre = np.zeros((n_patches, 3))
            for k in range(3):
                patch_normal[:, k] = np.bincount(patch_of_face, weights=areas * normals[:, k], minlength=n_patches)
                patch_centre[:, k] = np.bincount(patch_of_face, weights=areas * centroids[:, k], minlength=n_patches)
            patch_normal /= np.maximum(np.linalg.norm(patch_normal, axis=1), 1e-300)[:, None]
            patch_centre /= patch_area[:, None]
            patch_count = np.bincount(patch_of_face, minlength=n_patches)

    # --- 3. schedule patches -------------------------------------------------
    patch_has_small = np.zeros(n_patches, dtype=bool)
    patch_has_large = np.zeros(n_patches, dtype=bool)
    np.logical_or.at(patch_has_small, patch_of_face, too_small)
    np.logical_or.at(patch_has_large, patch_of_face, too_large)
    scheduled_patch = np.zeros(n_patches, dtype=bool)
    if allow_merge:
        scheduled_patch |= patch_has_small & (patch_count > 1)
    if allow_split:
        scheduled_patch |= patch_has_large

    # Edge classification. An edge is interior to a patch when exactly two
    # faces of that same patch share it; every other edge is a constraint.
    edge_counts = np.diff(edge_ptr)
    first_face = edge_idx[edge_ptr[:-1]]
    second_face = np.where(edge_counts >= 2, edge_idx[np.minimum(edge_ptr[:-1] + 1, len(edge_idx) - 1)], first_face)
    edge_interior = (edge_counts == 2) & (patch_of_face[first_face] == patch_of_face[second_face])

    # Vertices used by more than one patch are fixed points of every patch
    # using them (vertex-only contacts).
    vp = np.unique(np.column_stack([faces.reshape(-1), np.repeat(patch_of_face, 3)]), axis=0)
    vertex_shared = np.bincount(vp[:, 0], minlength=len(vertices)) > 1

    # Constraint-edge incidence per vertex (CSR), used below.
    constraint_edges = np.flatnonzero(~edge_interior)
    cons_deg = np.bincount(edges[constraint_edges].reshape(-1), minlength=len(vertices))
    ce_vertex = edges[constraint_edges].reshape(-1)
    ce_edge = np.repeat(constraint_edges, 2)
    order = np.argsort(ce_vertex, kind="stable")
    ce_edge_sorted = ce_edge[order]
    ce_ptr = np.concatenate([[0], np.cumsum(np.bincount(ce_vertex[order], minlength=len(vertices)))])

    def _removable_vertices(scheduled_face_now: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Boundary vertices that may be dropped: all their faces are scheduled,
        they have exactly two constraint edges, and the straight chain that
        replaces them stays within ``plane_tol`` of every dropped vertex
        (Douglas-Peucker on each run of such vertices, so a gently curving or
        stepped boundary is never flattened beyond the tolerance).
        Returns (vertex_owned, candidate, removable)."""
        vertex_unscheduled = np.zeros(len(vertices), dtype=bool)
        if np.any(~scheduled_face_now):
            vertex_unscheduled[faces[~scheduled_face_now].reshape(-1)] = True
        owned = ~vertex_unscheduled
        cand = owned & (cons_deg == 2)
        rem = np.zeros(len(vertices), dtype=bool)

        def other_end(edge_id: int, vertex_id: int) -> int:
            return int(edges[edge_id, 1]) if int(edges[edge_id, 0]) == vertex_id else int(edges[edge_id, 0])

        def walk(start: int, first_edge: int) -> Tuple[List[int], bool]:
            """Follow constraint edges from ``start`` through candidate vertices.
            Returns the vertices passed (excluding start) ending at an anchor,
            and whether the walk closed on ``start``."""
            chain: List[int] = []
            v_cur, e_cur = start, first_edge
            while True:
                v_next = other_end(e_cur, v_cur)
                chain.append(v_next)
                if v_next == start:
                    return chain, True
                if not cand[v_next]:
                    return chain, False
                inc = ce_edge_sorted[ce_ptr[v_next]:ce_ptr[v_next + 1]]
                e_cur = int(inc[0]) if int(inc[0]) != e_cur else int(inc[1])
                v_cur = v_next

        def douglas_peucker(run: List[int]) -> None:
            """Mark the interior vertices of ``run`` removable where the straight
            segment between its ends stays within plane_tol of them."""
            stack = [(0, len(run) - 1)]
            while stack:
                i, j = stack.pop()
                if j - i < 2:
                    continue
                a = vertices[run[i]]
                b = vertices[run[j]]
                ab = b - a
                ab_len2 = float(ab @ ab)
                inner = np.asarray(run[i + 1:j])
                pv = vertices[inner]
                if ab_len2 <= 0.0:
                    continue  # ends coincide: keep everything in between
                t = ((pv - a) @ ab) / ab_len2
                dist = np.linalg.norm(pv - (a + t[:, None] * ab), axis=1)
                ok = (dist <= plane_tol) & (t > 0.0) & (t < 1.0)
                if ok.all():
                    rem[inner] = True
                    continue
                k = int(np.argmax(np.where(ok, -1.0, dist)))
                stack.append((i, i + 1 + k))
                stack.append((i + 1 + k, j))

        visited = np.zeros(len(vertices), dtype=bool)
        for v0 in np.flatnonzero(cand):
            if visited[v0]:
                continue
            inc = ce_edge_sorted[ce_ptr[v0]:ce_ptr[v0 + 1]]
            back, closed = walk(int(v0), int(inc[0]))
            if closed:
                # A loop of candidates with no anchor: anchor it at the vertex
                # farthest from v0, then treat the two halves as runs.
                loop = [int(v0)] + back[:-1]
                visited[loop] = True
                far = int(np.argmax(np.linalg.norm(vertices[loop] - vertices[int(v0)], axis=1)))
                if far == 0:
                    continue
                douglas_peucker(loop[: far + 1])
                douglas_peucker(loop[far:] + [int(v0)])
                continue
            fwd, _ = walk(int(v0), int(inc[1]))
            run = back[::-1] + [int(v0)] + fwd       # anchor ... anchor
            visited[[x for x in run if cand[x]]] = True
            douglas_peucker(run)
        return owned, cand, rem

    # Per-patch tables, built once and shared by the per-patch tasks: the
    # faces, vertices and constraint edges of every patch (CSR, sorted within
    # a patch), the patches of every vertex, and the in-plane bases.
    face_order = np.argsort(patch_of_face, kind="stable")
    face_ptr = np.concatenate([[0], np.cumsum(patch_count)])
    vpat_ptr = np.concatenate([[0], np.cumsum(np.bincount(vp[:, 0], minlength=len(vertices)))])
    vpat_idx = vp[:, 1]
    order = np.lexsort((vp[:, 0], vp[:, 1]))
    pv_idx = vp[order, 0]
    pv_ptr = np.concatenate([[0], np.cumsum(np.bincount(vp[order, 1], minlength=n_patches))])
    pe_key = np.unique(np.repeat(patch_of_face, 3).astype(np.int64) * n_edges + face_edge_ids.reshape(-1))
    pe_patch, pe_edge = np.divmod(pe_key, n_edges)
    keep_pe = ~edge_interior[pe_edge]
    pce_idx = pe_edge[keep_pe]
    pce_ptr = np.concatenate([[0], np.cumsum(np.bincount(pe_patch[keep_pe], minlength=n_patches))])
    with np.errstate(divide="ignore", invalid="ignore"):
        patch_u, patch_v = _plane_bases(patch_normal)
    tables = dict(
        vertices=vertices, faces=faces, edges=edges, areas=areas, plane_tol=plane_tol,
        face_order=face_order, face_ptr=face_ptr, pv_ptr=pv_ptr, pv_idx=pv_idx,
        pce_ptr=pce_ptr, pce_idx=pce_idx, patch_centre=patch_centre, patch_area=patch_area,
        patch_u=patch_u, patch_v=patch_v,
    )

    # Pre-screen: a scheduled patch whose boundary Triangle cannot be trusted
    # with (a needle of slivers, a self-overlapping boundary, too few points
    # once collinear boundary vertices are dropped, ...) is taken off the
    # schedule here, before any edge is touched, so that its untouched faces
    # stay conforming with their neighbours. Taking a patch off the schedule
    # changes which vertices are owned, so iterate until nothing changes.
    skipped: Dict[str, int] = {}
    to_check = set(int(p) for p in np.flatnonzero(scheduled_patch))
    while to_check:
        scheduled_face = scheduled_patch[patch_of_face]
        vertex_owned, candidate, removable = _removable_vertices(scheduled_face)
        ids = sorted(to_check)
        reasons = _map_patches(_prescreen_task, ids, dict(tables, vertex_shared=vertex_shared, removable=removable), workers)
        unscheduled_now: List[int] = []
        for p, reason in zip(ids, reasons):
            if reason is not None:
                scheduled_patch[p] = False
                unscheduled_now.append(p)
                skipped[reason] = skipped.get(reason, 0) + 1
                logger.debug("[uniform] patch %d kept as is: %s", int(p), reason)
        if not unscheduled_now:
            break
        # Neighbours of the patches just taken off the schedule lose owned
        # vertices; re-check them (and only them).
        affected = np.zeros(n_patches, dtype=bool)
        for p in unscheduled_now:
            for vtx in pv_idx[pv_ptr[p]:pv_ptr[p + 1]]:
                affected[vpat_idx[vpat_ptr[vtx]:vpat_ptr[vtx + 1]]] = True
        to_check = set(int(p) for p in np.flatnonzero(affected & scheduled_patch))

    scheduled_face = scheduled_patch[patch_of_face]
    vertex_owned, candidate, removable = _removable_vertices(scheduled_face)
    n_vertices_removed = int(removable.sum())
    n_unscheduled = int(np.count_nonzero(~scheduled_face))
    n_faces_merged = int(np.count_nonzero(too_small & scheduled_face))
    n_faces_split = int(np.count_nonzero(too_large & scheduled_face))
    n_small_isolated = int(np.count_nonzero(too_small & ~scheduled_face))
    # An edge is owned when all incident faces are scheduled: only then may it be split.
    sched_per_entry = scheduled_face[edge_idx].astype(int)
    edge_owned = np.add.reduceat(sched_per_entry, edge_ptr[:-1]) == edge_counts

    # --- 4. constraint chains: simplify, then split ----------------------------
    # Removable vertices are dropped and their two edges become one straight
    # chain; the chain is then re-split at the target edge length. Both happen
    # in one global table, so every patch along the chain sees the same
    # vertices and the mesh stays conforming.
    # Chains: walk from non-removable endpoints through removable vertices.
    # chain_of_edge maps every constraint edge to its chain; chain_ends holds
    # the end vertices of each chain. An edge without a removable endpoint is a
    # chain of its own.
    chain_of_edge = -np.ones(n_edges, dtype=int)
    chain_ends: List[Tuple[int, int]] = []
    chain_owned: List[bool] = []
    for e in constraint_edges:
        if chain_of_edge[e] >= 0:
            continue
        a, b = int(edges[e, 0]), int(edges[e, 1])
        if not (removable[a] or removable[b]):
            chain_of_edge[e] = len(chain_ends)
            chain_ends.append((a, b))
            chain_owned.append(bool(edge_owned[e]))
            continue
        cid = len(chain_ends)
        chain_of_edge[e] = cid
        ends = []
        for start_v in (a, b):
            v_cur, e_cur = start_v, int(e)
            while removable[v_cur]:
                inc = ce_edge_sorted[ce_ptr[v_cur]:ce_ptr[v_cur + 1]]
                e_next = int(inc[0]) if int(inc[0]) != e_cur else int(inc[1])
                if chain_of_edge[e_next] >= 0:  # closed-loop safeguard
                    break
                chain_of_edge[e_next] = cid
                v_cur = int(edges[e_next, 0]) if int(edges[e_next, 1]) == v_cur else int(edges[e_next, 1])
                e_cur = e_next
            ends.append(v_cur)
        chain_ends.append((min(ends), max(ends)))
        chain_owned.append(True)  # a removable vertex only has scheduled faces
    chain_ends_arr = np.asarray(chain_ends, dtype=int).reshape(-1, 2)
    chain_owned_arr = np.asarray(chain_owned, dtype=bool)
    # Chains with the same end vertices are the same straight segment (a chain
    # of collinear sliver edges next to one long edge of the neighbouring
    # patch): give them one id so both sides receive identical split vertices.
    if len(chain_ends_arr):
        _, canonical, inverse = np.unique(chain_ends_arr, axis=0, return_index=True, return_inverse=True)
        inverse = np.asarray(inverse).reshape(-1)
        if len(canonical) < len(chain_ends_arr):
            remap = canonical[inverse]
            # A merged chain is owned only if every member was owned.
            owned_all = np.ones(len(chain_ends_arr), dtype=bool)
            np.logical_and.at(owned_all, remap, chain_owned_arr)
            chain_owned_arr = owned_all
            valid_edges = chain_of_edge >= 0
            chain_of_edge[valid_edges] = remap[chain_of_edge[valid_edges]]
    n_chains = len(chain_ends_arr)

    # Split owned chains longer than the target edge length.
    if n_chains:
        chain_len = np.linalg.norm(vertices[chain_ends_arr[:, 0]] - vertices[chain_ends_arr[:, 1]], axis=1)
        chain_segments = np.maximum(1, np.rint(chain_len / target_edge).astype(int))
    else:
        chain_segments = np.zeros(0, dtype=int)
    split_chain = chain_owned_arr & (chain_segments >= 2)
    n_edges_split = int(split_chain.sum())

    new_vertices: List[np.ndarray] = [vertices]
    n_vertices_now = len(vertices)
    chain_split_vertices: Dict[int, np.ndarray] = {}
    for c in np.flatnonzero(split_chain):
        k = int(chain_segments[c])
        a, b = chain_ends_arr[c]
        t = np.arange(1, k)[:, None] / k
        pts = vertices[a] * (1.0 - t) + vertices[b] * t
        ids = np.arange(n_vertices_now, n_vertices_now + k - 1)
        n_vertices_now += k - 1
        new_vertices.append(pts)
        chain_split_vertices[int(c)] = ids

    # Removed (collinear) boundary vertices are no longer fixed points.
    vertex_shared = vertex_shared & ~removable

    if _debug is not None:  # internal arrays for diagnostics and tests
        _debug.update(dict(vertices=vertices, faces=faces, keep_index=keep_index, areas=areas,
                           patch_of_face=patch_of_face, scheduled_face=scheduled_face,
                           edges=edges, edge_interior=edge_interior, edge_owned=edge_owned,
                           cons_deg=cons_deg, vertex_owned=vertex_owned, removable=removable,
                           candidate=candidate, chain_of_edge=chain_of_edge, chain_ends=chain_ends_arr,
                           area_low=area_low, area_high=area_high))

    # --- 5. count guard (prediction) ---------------------------------------
    sched_ids = np.flatnonzero(scheduled_patch)
    predicted = n_unscheduled
    if sched_ids.size:
        predicted += int(np.sum(np.maximum(1.0, patch_area[sched_ids] / (max_area * _MEAN_AREA_FRACTION))))
    if predicted > max_faces:
        raise FacetCountError(
            f"uniformising is predicted to produce about {predicted} facets from {n_in}, "
            f"more than the allowed {max_faces}. Increase the target edge length "
            f"({target_edge:.3g} m now; the fixed-count choice is "
            f"{target_edge_from_area(total_area / n):.3g} m) or raise the facet count limit."
        )

    # Vertex array after chain splitting; Steiner vertices of the patches are
    # appended separately (no segment or fixed point ever refers to them).
    base_vertices = np.concatenate(new_vertices, axis=0) if len(new_vertices) > 1 else vertices
    steiner_vertices: List[np.ndarray] = []

    # --- 6. retriangulate scheduled patches ---------------------------------
    out_faces: List[np.ndarray] = [faces[~scheduled_face]]
    out_types: List[np.ndarray] = [types[~scheduled_face]]
    out_map: List[np.ndarray] = [keep_index[~scheduled_face]]
    out_from_patch_count = 0
    n_boundary_points = 0
    n_overlap_patches = 0
    segment_lengths: List[np.ndarray] = []

    if sched_ids.size:
        # Sub-regions inside a patch connected through interior edges (a
        # non-manifold constraint edge inside a patch separates two regions,
        # each of which needs its own region marker).
        ie = np.flatnonzero(edge_interior)
        region_graph = coo_matrix((np.ones(ie.size), (first_face[ie], second_face[ie])), shape=(n, n))
        _, region_of_face = connected_components(region_graph, directed=False)

        results = _map_patches(
            _remesh_task, sched_ids,
            dict(tables, vertex_shared=vertex_shared, base_vertices=base_vertices,
                 chain_of_edge=chain_of_edge, chain_ends=chain_ends_arr,
                 chain_split_vertices=chain_split_vertices, region_of_face=region_of_face,
                 max_area=max_area, min_angle=min_angle),
            workers,
        )
        for p, (point_ids, tris, steiner, fmap_local, lengths, overlap) in zip(sched_ids, results):
            pf = face_order[face_ptr[p]:face_ptr[p + 1]]
            n_boundary_points += len(point_ids)
            if lengths is not None:
                segment_lengths.append(lengths)
            n_overlap_patches += int(overlap)
            # Steiner vertices are numbered in patch order, after the chain-split vertices.
            n_new = len(steiner)
            global_ids = np.concatenate([point_ids, np.arange(n_vertices_now, n_vertices_now + n_new)])
            if n_new > 0:
                steiner_vertices.append(steiner)
                n_vertices_now += n_new
            out_faces.append(global_ids[tris])
            out_types.append(np.full(len(tris), types[pf[0]], dtype=int))
            out_map.append(keep_index[pf[fmap_local]])
            out_from_patch_count += len(tris)

    n_patches_kept = int(sum(skipped.values()))
    if n_patches_kept:
        logger.warning(
            "[uniform] %d patch(es) due for retriangulation kept unchanged: %s",
            n_patches_kept, ", ".join(f"{k}: {v}" for k, v in sorted(skipped.items())),
        )
    if n_overlap_patches:
        logger.warning(
            "[uniform] %d patch(es) had overlapping or tilted faces; their area is now covered once",
            n_overlap_patches,
        )

    # --- 7. assemble ---------------------------------------------------------
    all_vertices = np.concatenate([base_vertices] + steiner_vertices, axis=0) if steiner_vertices else base_vertices
    faces_out = np.concatenate(out_faces, axis=0)
    types_out = np.concatenate(out_types, axis=0)
    map_out = np.concatenate(out_map, axis=0)
    if len(faces_out) > max_faces:
        raise FacetCountError(
            f"uniformising produced {len(faces_out)} facets from {n_in}, more than the allowed {max_faces}. "
            f"Increase the target edge length ({target_edge:.3g} m now) or raise the facet count limit."
        )
    mesh_out = trimesh.Trimesh(vertices=all_vertices, faces=faces_out, process=False)
    mesh_out.merge_vertices(digits_vertex=int(weld_digits))
    verts_out, faces_out = _remove_unused(np.asarray(mesh_out.vertices, dtype=float), np.asarray(mesh_out.faces, dtype=int))
    mesh_out = trimesh.Trimesh(vertices=verts_out, faces=faces_out, process=False)
    areas_out = np.asarray(mesh_out.area_faces, dtype=float)

    report: Dict[str, Any] = {
        "target_edge": float(target_edge),
        "target_area": float(target_area),
        "target_auto": bool(target_auto),
        "size_band": float(size_band),
        "area_low": float(area_low),
        "area_high": float(area_high),
        "allow_split": bool(allow_split),
        "allow_merge": bool(allow_merge),
        "max_faces": int(max_faces),
        "plane_tol": float(plane_tol),
        "min_angle": float(min_angle),
        "n_faces_in": int(n_in),
        "n_faces_out": int(len(faces_out)),
        "n_vertices_in": int(len(np.asarray(mesh.vertices))),
        "n_vertices_out": int(len(verts_out)),
        "n_degenerate_removed": n_degenerate,
        "area_in": float(np.asarray(mesh.area_faces, dtype=float).sum()),
        "area_out": float(areas_out.sum()),
        "n_patches": int(n_patches),
        "n_patches_remeshed": int(sched_ids.size),
        "n_patches_noncoplanar": int(n_noncoplanar),
        "n_patches_kept": n_patches_kept,
        "n_patches_overlapping": int(n_overlap_patches),
        "patches_kept_reasons": dict(sorted(skipped.items())),
        "n_faces_merged": n_faces_merged,
        "n_faces_split": n_faces_split,
        "n_faces_from_patches": int(out_from_patch_count),
        "n_small_isolated": n_small_isolated,
        "n_edges_split": n_edges_split,
        "n_boundary_points_remeshed": int(n_boundary_points),
        "boundary_segment_percentiles": _percentiles(np.concatenate(segment_lengths) if segment_lengths else np.zeros(0)),
        "n_boundary_vertices_removed": n_vertices_removed,
        "n_out_below_band": int(np.count_nonzero(areas_out < area_low)),
        "n_out_above_band": int(np.count_nonzero(areas_out > area_high)),
        "area_percentiles_in": _percentiles(areas),
        "area_percentiles_out": _percentiles(areas_out),
        "elapsed_s": time.perf_counter() - t0,
    }
    return UniformFacetsResult(
        mesh=mesh_out,
        face_types=types_out if face_types is not None else None,
        face_map=map_out,
        report=report,
    )
