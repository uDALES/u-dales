"""
split_buildings - Separate connected building components from triangulated mesh

This module provides functionality to identify and separate individual
connected building components from complex urban geometries.

Copyright (C) 2025 the uDALES Team.
"""

import numpy as np
from typing import List, Tuple, Optional
import warnings

try:
    import trimesh
    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False


def split_buildings(mesh: 'trimesh.Trimesh', remove_ground: bool = True) -> Tuple[List['trimesh.Trimesh'], np.ndarray]:
    """
    Separate connected building components from a triangulated mesh.
    
    Analyzes a triangulated mesh to identify and separate individual
    connected building components. The function automatically removes
    ground-level faces (Z=0) by default, then uses graph connectivity
    analysis to group triangular faces that share edges.
    
    Parameters
    ----------
    mesh : trimesh.Trimesh
        Triangulated mesh containing the geometry to analyze.
        Can include ground faces (they will be automatically removed if remove_ground=True)
    remove_ground : bool, default=True
        If True, removes faces with all vertices at Z=0 before splitting
    
    Returns
    -------
    building_components : list of trimesh.Trimesh
        List where each element is a trimesh object representing one
        connected building component. Buildings are sorted by face count
        (largest first).
    face_to_building_map : ndarray, shape (n_faces,)
        Array mapping each face index to its building ID (1-indexed).
        Ground faces (if retained) are assigned building ID 0.
        Building IDs start from 1.
    
    Algorithm
    ---------
    1. Optional ground removal (faces with Z≤0 vertices)
    2. Edge-based graph construction (shared edges = connections)
    3. Connected component analysis to group faces
    4. Triangulation object creation for each component
    5. Component sorting by face count (largest first)
    
    Connectivity Analysis
    ---------------------
    - Edge-based connectivity: Faces sharing edges are considered connected
    - Graph components: Groups of transitively connected faces
    - Component separation: Each connected group becomes a separate building
    - Topology preservation: Maintains mesh integrity within each component
    
    Examples
    --------
    >>> import trimesh
    >>> urban_mesh = trimesh.load('urban_scene.stl')
    >>> buildings, face_map = split_buildings(urban_mesh)
    >>> for i, building in enumerate(buildings):
    ...     print(f"Building {i+1}: {len(building.faces)} faces")
    
    >>> # Keep ground faces if needed
    >>> buildings, face_map = split_buildings(urban_mesh, remove_ground=False)
    
    Applications
    ------------
    - Individual building analysis in urban environments
    - Separation of merged CAD models
    - Building-specific flow simulation setup
    - Urban morphology studies
    - Facet-based energy balance calculations per building
    
    See Also
    --------
    calculate_outline : Calculate outline edges for visualization
    UDGeom.get_buildings : Cached building access method
    """
    if not TRIMESH_AVAILABLE:
        raise ImportError("trimesh is required. Install with: pip install trimesh")
    
    vertices = mesh.vertices
    faces = mesh.faces
    
    # Step 1: Optional ground removal
    if remove_ground:
        # Identify faces with all vertices at ground level (Z == 0)
        # This matches MATLAB's deleteGround implementation
        face_z_values = vertices[faces, 2]  # Shape: (n_faces, 3)
        ground_mask = np.all(face_z_values == 0, axis=1)
        building_faces = faces[~ground_mask]
        building_face_indices = np.where(~ground_mask)[0]
    else:
        building_faces = faces
        building_face_indices = np.arange(len(faces))
    
    if len(building_faces) == 0:
        warnings.warn("No building faces found after ground removal")
        return [], np.zeros(len(faces), dtype=int)
    
    # Step 2: Build edge-to-face mapping for connectivity
    edge_to_faces = {}
    
    for idx, (orig_idx, face) in enumerate(zip(building_face_indices, building_faces)):
        # Create edges (sorted pairs of vertices)
        edges = [
            tuple(sorted([face[0], face[1]])),
            tuple(sorted([face[1], face[2]])),
            tuple(sorted([face[2], face[0]]))
        ]
        
        for edge in edges:
            if edge not in edge_to_faces:
                edge_to_faces[edge] = []
            edge_to_faces[edge].append(idx)  # Use local index
    
    # Step 3: Build adjacency graph
    # Each face is a node, edges connect faces that share an edge
    n_building_faces = len(building_faces)
    adjacency = {i: set() for i in range(n_building_faces)}
    
    for edge, face_list in edge_to_faces.items():
        if len(face_list) > 1:
            # Faces sharing this edge are connected
            for i in range(len(face_list)):
                for j in range(i + 1, len(face_list)):
                    adjacency[face_list[i]].add(face_list[j])
                    adjacency[face_list[j]].add(face_list[i])
    
    # Step 4: Find connected components using depth-first search
    visited = set()
    components = []
    
    def dfs(node, component):
        """Depth-first search to find connected component"""
        visited.add(node)
        component.append(node)
        for neighbor in adjacency[node]:
            if neighbor not in visited:
                dfs(neighbor, component)
    
    for i in range(n_building_faces):
        if i not in visited:
            component = []
            dfs(i, component)
            components.append(component)
    
    # Step 5: Create triangulation objects for each component
    building_components = []
    face_to_building_map = np.zeros(len(faces), dtype=int)
    
    # Sort components by size (largest first)
    components.sort(key=lambda x: len(x), reverse=True)
    
    for building_id, component in enumerate(components, start=1):
        # Get original face indices for this component
        component_face_indices = building_face_indices[component]
        component_faces = faces[component_face_indices]
        
        # Extract submesh for this building
        try:
            # Use trimesh's submesh functionality
            building_mesh = mesh.submesh([component_face_indices], append=True)
            
            # If submesh returns a list (when append=True), take the first element
            if isinstance(building_mesh, list):
                if len(building_mesh) > 0:
                    building_mesh = building_mesh[0]
                else:
                    continue
            
            building_components.append(building_mesh)
            
            # Map original faces to building ID
            face_to_building_map[component_face_indices] = building_id
            
        except Exception as e:
            warnings.warn(f"Could not create submesh for building {building_id}: {e}")
            continue
    
    return building_components, face_to_building_map
