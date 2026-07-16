"""
delete_ground - Remove likely ground faces from a triangulated mesh.

Historically this mirrored MATLAB `deleteGround` by removing faces on ``z=0``.
For topography-aware Python workflows, it now relies on ``identify_ground_faces``.
"""

from typing import Tuple

import numpy as np

from .check_mesh import identify_ground_faces

try:
    import trimesh
    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False

from exceptions import DependencyError
def delete_ground(mesh: "trimesh.Trimesh") -> Tuple["trimesh.Trimesh", np.ndarray]:
    """
    Remove likely ground faces and compact the mesh.

    Parameters
    ----------
    mesh : trimesh.Trimesh
        Input surface mesh.

    Returns
    -------
    filtered_mesh : trimesh.Trimesh
        Mesh with all-ground faces removed and vertices compacted.
    face_mapping : ndarray
        Mapping from filtered face indices back to the original face indices.
    """
    if not TRIMESH_AVAILABLE:
        raise DependencyError("trimesh is required. Install with: pip install trimesh")

    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)

    if len(faces) == 0:
        empty = trimesh.Trimesh(
            vertices=np.empty((0, 3), dtype=float),
            faces=np.empty((0, 3), dtype=int),
            process=False,
        )
        return empty, np.empty((0,), dtype=int)

    ground_mask = identify_ground_faces(mesh)
    kept_face_indices = np.where(~ground_mask)[0]
    kept_faces = faces[kept_face_indices]

    if len(kept_faces) == 0:
        empty = trimesh.Trimesh(
            vertices=np.empty((0, 3), dtype=float),
            faces=np.empty((0, 3), dtype=int),
            process=False,
        )
        return empty, kept_face_indices

    used_vertex_indices = np.unique(kept_faces.reshape(-1))
    new_vertices = vertices[used_vertex_indices]

    new_index = -np.ones(len(vertices), dtype=int)
    new_index[used_vertex_indices] = np.arange(len(used_vertex_indices))
    new_faces = new_index[kept_faces]

    filtered = trimesh.Trimesh(vertices=new_vertices, faces=new_faces, process=False)
    return filtered, kept_face_indices


__all__ = ["delete_ground"]
