"""
extrude_to_ground - Drop a selected independent surface to planar ground and restitch it.

This targets the common import case where one disconnected surface, typically a
floating building, should be extended vertically down to an identified planar
ground surface and then welded back into the existing ground mesh.
"""

from __future__ import annotations

from typing import Dict, Tuple, Union

import numpy as np

from .check_mesh import _as_trimesh, calculate_independent_surfaces, identify_ground_faces
from .truncate_below_ground import (
    _build_mesh_from_triangles,
    _clip_ground_against_footprint,
    _copy_mesh,
    _estimate_planar_ground_level,
)

try:
    import trimesh

    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False

try:
    from exceptions import DependencyError
except ImportError:
    from ..exceptions import DependencyError


def extrude_to_ground(
    mesh_or_geom,
    surface_id: int,
    *,
    ground_planarity_tolerance: float = 1.0e-6,
    bottom_vertex_tolerance: float = 1.0e-8,
    return_trimesh: bool = False,
) -> Tuple[Union["trimesh.Trimesh", "UDGeom"], Dict[str, Union[int, float, bool, list]]]:
    """
    Move the lowest vertices of one independent surface down to the ground plane.

    This routine is intended for simple floating-surface cases, for example a
    building shell that is disconnected from the ground. The selected surface
    is identified by its independent-surface id from
    ``calculate_independent_surfaces(...)``.

    The repair strategy is:
    1. identify the planar ground level
    2. select one non-ground independent surface
    3. move that surface's lowest vertices vertically down to the ground
    4. clip the existing ground locally against the updated footprint
    5. weld the resulting seam back into a single stitched surface mesh

    Parameters
    ----------
    mesh_or_geom : trimesh.Trimesh or UDGeom
        Mesh to modify.
    surface_id : int
        Independent-surface id from ``calculate_independent_surfaces``.
    ground_planarity_tolerance : float, default=1e-6
        Maximum allowed deviation from a single ground plane.
    bottom_vertex_tolerance : float, default=1e-8
        Tolerance used to identify the surface's bottom vertices.
    return_trimesh : bool, default=False
        If True, return the repaired ``trimesh.Trimesh`` directly. Otherwise
        return a ``UDGeom`` object.

    Returns
    -------
    cleaned_geom, report
        Repaired geometry as ``UDGeom`` by default, plus a summary report
        describing the moved vertices and local ground restitching.
    """
    if not TRIMESH_AVAILABLE:
        raise DependencyError("trimesh is required. Install with: pip install trimesh")

    mesh = _copy_mesh(_as_trimesh(mesh_or_geom))
    if len(mesh.faces) == 0:
        raise ValueError("Cannot extrude an empty mesh.")

    surface_result = calculate_independent_surfaces(mesh)
    surfaces = {int(item["surface_id"]): item for item in surface_result["surfaces"]}
    if int(surface_id) not in surfaces:
        raise ValueError(f"Unknown surface_id {surface_id}. Valid ids are {sorted(surfaces)}")

    ground_mask = identify_ground_faces(mesh)
    ground_face_count = int(np.count_nonzero(ground_mask))
    z0 = _estimate_planar_ground_level(mesh, ground_mask, ground_planarity_tolerance)

    selected_face_ids = np.asarray(surfaces[int(surface_id)]["face_ids"], dtype=int)
    if np.any(ground_mask[selected_face_ids]):
        raise ValueError("Selected surface overlaps the identified ground. Choose a non-ground surface id.")

    vertices = np.asarray(mesh.vertices, dtype=float).copy()
    vertices[:, 2] -= z0
    shifted = trimesh.Trimesh(vertices=vertices, faces=np.asarray(mesh.faces, dtype=int).copy(), process=False)

    shifted_faces = np.asarray(shifted.faces, dtype=int)
    selected_vertex_ids = np.unique(shifted_faces[selected_face_ids].reshape(-1))
    selected_vertices = np.asarray(shifted.vertices, dtype=float)[selected_vertex_ids]
    bottom_z = float(np.min(selected_vertices[:, 2]))
    bottom_vertex_ids = selected_vertex_ids[
        np.abs(selected_vertices[:, 2] - bottom_z) <= float(bottom_vertex_tolerance)
    ]

    shifted_vertices = np.asarray(shifted.vertices, dtype=float).copy()
    shifted_vertices[bottom_vertex_ids, 2] = 0.0
    shifted = trimesh.Trimesh(vertices=shifted_vertices, faces=shifted_faces.copy(), process=False)

    selected_faces_shifted = shifted_faces[selected_face_ids].copy()
    selected_triangles = np.asarray(shifted.vertices, dtype=float)[selected_faces_shifted]
    selected_normals = trimesh.Trimesh(
        vertices=np.asarray(shifted.vertices, dtype=float).copy(),
        faces=selected_faces_shifted.copy(),
        process=False,
    ).face_normals
    bottom_face_mask = np.all(np.isclose(selected_triangles[:, :, 2], 0.0, atol=bottom_vertex_tolerance), axis=1)
    bottom_face_mask &= selected_normals[:, 2] < -1.0e-6
    removed_bottom_faces = int(np.count_nonzero(bottom_face_mask))
    if removed_bottom_faces > 0:
        selected_faces_shifted = selected_faces_shifted[~bottom_face_mask]

    selected_surface_mesh = trimesh.Trimesh(
        vertices=np.asarray(shifted.vertices, dtype=float).copy(),
        faces=selected_faces_shifted,
        process=False,
    )
    other_face_mask = np.ones(len(shifted_faces), dtype=bool)
    other_face_mask[selected_face_ids] = False
    ground_faces = shifted_faces[ground_mask].copy()
    other_non_ground_faces = shifted_faces[other_face_mask & ~ground_mask].copy()

    from .geometry_generation import _ground_footprint_union

    footprint_union = None
    if len(selected_surface_mesh.faces) > 0:
        footprint_union = _ground_footprint_union(selected_surface_mesh)

    ground_mesh = trimesh.Trimesh(
        vertices=np.asarray(shifted.vertices, dtype=float).copy(),
        faces=ground_faces,
        process=False,
    )
    clipped_ground_triangles, preserved_ground_faces, retriangulated_ground_faces = _clip_ground_against_footprint(
        ground_mesh,
        footprint_union,
        z0=0.0,
    )

    non_ground_triangles = [
        np.asarray(shifted.vertices, dtype=float)[face].copy()
        for face in np.vstack([selected_surface_mesh.faces, other_non_ground_faces]) if len(face) == 3
    ]
    combined_triangles = clipped_ground_triangles + non_ground_triangles
    cleaned = _build_mesh_from_triangles(combined_triangles)
    cleaned_vertices = np.asarray(cleaned.vertices, dtype=float).copy()
    cleaned_vertices[:, 2] += z0
    cleaned = trimesh.Trimesh(
        vertices=cleaned_vertices,
        faces=np.asarray(cleaned.faces, dtype=int).copy(),
        process=False,
    )
    cleaned.merge_vertices(digits_vertex=8)
    cleaned.remove_unreferenced_vertices()

    from .fix_mesh import weld_touching_boundaries

    cleaned, weld_report = weld_touching_boundaries(cleaned)
    report = {
        "surface_id": int(surface_id),
        "ground_face_count": int(ground_face_count),
        "ground_z": float(z0),
        "bottom_z_before": float(bottom_z + z0),
        "moved_vertices": int(len(bottom_vertex_ids)),
        "bottom_vertex_ids": bottom_vertex_ids.astype(int).tolist(),
        "removed_bottom_faces": int(removed_bottom_faces),
        "preserved_ground_faces": int(preserved_ground_faces),
        "retriangulated_ground_faces": int(retriangulated_ground_faces),
        "welded_touching_boundary_faces": int(weld_report["welded_touching_boundary_faces"]),
    }
    if return_trimesh:
        return cleaned, report

    from .udgeom import UDGeom, DEFAULT_BACKEND

    return UDGeom(stl=cleaned, backend=getattr(mesh_or_geom, "backend", DEFAULT_BACKEND)), report


__all__ = ["extrude_to_ground"]
