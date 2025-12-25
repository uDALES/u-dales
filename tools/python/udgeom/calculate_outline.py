"""
calculate_outline - Calculate boundary and sharp edges for mesh outlining

This module provides functionality to analyze triangulated meshes and identify
edges that should be outlined, including both true boundary edges and sharp
geometric features based on face normal angles.

Copyright (C) 2025 the uDALES Team.
"""

import numpy as np
from typing import Tuple, Dict, Optional
import warnings

try:
    import trimesh
    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False


def calculate_outline(mesh: 'trimesh.Trimesh', angle_threshold: float = 45.0) -> Tuple[np.ndarray, Dict]:
    """
    Calculate boundary and sharp edges for mesh outlining.
    
    Analyzes a triangulated mesh to identify edges that should be outlined,
    including both true boundary edges and sharp geometric features based
    on face normal angles.
    
    Parameters
    ----------
    mesh : trimesh.Trimesh
        Triangulated mesh containing vertices and faces
    angle_threshold : float, default=45.0
        Angle threshold in degrees for detecting sharp edges.
        Edges where adjacent faces meet at angles greater than this
        threshold are considered outline edges.
    
    Returns
    -------
    boundary_edges : ndarray, shape (K, 2)
        Array of vertex index pairs defining edges to outline.
        Each row contains [vertex_idx1, vertex_idx2] for one edge.
    mesh_data : dict
        Dictionary containing mesh information:
        - 'faces': Triangle connectivity matrix
        - 'vertices': Vertex coordinates
        - 'normals': Face normal vectors
        - 'edge_map': Dict mapping edges to adjacent face indices
        - 'angle_threshold': The threshold used
    
    Algorithm
    ---------
    1. Computes face normal vectors for each triangle
    2. Maps each edge to its adjacent faces
    3. Identifies boundary edges (1 adjacent face)
    4. Identifies sharp edges (angle between normals > threshold)
    
    Examples
    --------
    >>> import trimesh
    >>> mesh = trimesh.load('building.stl')
    >>> edges, data = calculate_outline(mesh, angle_threshold=45.0)
    >>> print(f"Found {len(edges)} outline edges")
    
    >>> # Use with custom threshold for more sensitive detection
    >>> edges, data = calculate_outline(mesh, angle_threshold=30.0)
    
    See Also
    --------
    UDGeom.show_outline : Visualize outline edges
    split_buildings : Separate individual buildings from mesh
    """
    if not TRIMESH_AVAILABLE:
        raise ImportError("trimesh is required. Install with: pip install trimesh")
    
    # Extract mesh data
    vertices = mesh.vertices
    faces = mesh.faces
    
    # Step 1: Compute face normal vectors
    # trimesh already provides normalized face normals
    normals = mesh.face_normals
    
    # Step 2: Build edge-to-face mapping
    # Create a dictionary mapping each edge to the faces that contain it
    edge_map = {}
    
    for face_idx, face in enumerate(faces):
        # Extract the three edges of the current triangle
        # Sort vertex indices to ensure consistent edge representation
        edges = [
            tuple(sorted([face[0], face[1]])),  # Edge 0-1
            tuple(sorted([face[1], face[2]])),  # Edge 1-2
            tuple(sorted([face[2], face[0]]))   # Edge 2-0
        ]
        
        # Add current face to each edge's face list
        for edge in edges:
            if edge not in edge_map:
                edge_map[edge] = []
            edge_map[edge].append(face_idx)
    
    # Step 3: Identify boundary and sharp edges
    # Determine which edges should be outlined based on two criteria:
    # 1. Boundary edges: belong to only one face (mesh boundary)
    # 2. Sharp edges: angle between adjacent face normals > threshold
    
    boundary_edges = []
    
    for edge, face_ids in edge_map.items():
        if len(face_ids) == 1:
            # TRUE BOUNDARY EDGE: belongs to only one face
            boundary_edges.append(edge)
        else:
            # INTERNAL EDGE: check if it's a sharp geometric feature
            # Calculate maximum angle between any pair of adjacent faces
            max_angle = 0.0
            
            for i in range(len(face_ids) - 1):
                for j in range(i + 1, len(face_ids)):
                    n1 = normals[face_ids[i]]
                    n2 = normals[face_ids[j]]
                    
                    # Calculate angle between normal vectors
                    cos_theta = np.dot(n1, n2)
                    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # Clamp for numerical stability
                    angle = np.degrees(np.arccos(cos_theta))
                    max_angle = max(max_angle, angle)
            
            # Mark as boundary if angle exceeds threshold (sharp feature)
            if max_angle > angle_threshold:
                boundary_edges.append(edge)
    
    # Convert to numpy array
    boundary_edges = np.array(boundary_edges, dtype=int)
    
    # Package mesh data for output
    mesh_data = {
        'faces': faces,
        'vertices': vertices,
        'normals': normals,
        'edge_map': edge_map,
        'angle_threshold': angle_threshold
    }
    
    return boundary_edges, mesh_data
