"""
udgeom - Geometry class for uDALES

This module provides the UDGeom class for handling triangulated surface geometries
used in uDALES simulations. It wraps trimesh functionality for STL file I/O and
visualization.

uDALES (https://github.com/uDALES/u-dales)

Copyright (C) 2016-2024 the uDALES Team.
Licensed under GNU General Public License v3.0
"""

from pathlib import Path
from typing import List, Optional, Union
import numpy as np
import warnings

try:
    import trimesh
    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False
    warnings.warn("trimesh not installed. Geometry functionality will be limited.")

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    warnings.warn("matplotlib not installed. Visualization functionality will be limited.")

try:
    from scipy.spatial import ConvexHull
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

# Import package functions
from .calculate_outline import calculate_outline
from .split_buildings import split_buildings


class UDGeom:
    """
    Geometry class for uDALES.
    
    This class handles triangulated surface geometries (STL files) used in
    uDALES simulations. It provides methods for loading, saving, and visualizing
    3D geometries.
    
    Parameters
    ----------
    path : str or Path, optional
        Path to the directory containing geometry files. Defaults to current directory.
    stl : trimesh.Trimesh, optional
        Pre-existing trimesh object to initialize with.
        
    Attributes
    ----------
    path : Path
        Directory path for geometry files
    stl : trimesh.Trimesh or None
        Triangulated surface geometry
    
    Examples
    --------
    Load a geometry from STL file:
    >>> geom = UDGeom('path/to/simulation')
    >>> geom.load('geometry.stl')
    >>> geom.show()
    
    Create from existing trimesh:
    >>> import trimesh
    >>> mesh = trimesh.creation.box()
    >>> geom = UDGeom(stl=mesh)
    >>> geom.save('box.stl')
    """
    
    def __init__(self, path: Optional[Union[str, Path]] = None, 
                 stl: Optional['trimesh.Trimesh'] = None):
        """
        Initialize UDGeom object.
        
        Parameters
        ----------
        path : str or Path, optional
            Path to geometry directory. If None, uses current directory.
        stl : trimesh.Trimesh, optional
            Pre-existing trimesh object.
        """
        if not TRIMESH_AVAILABLE:
            raise ImportError("trimesh is required for UDGeom. Install with: pip install trimesh")
        
        # Set path
        if path is not None:
            self.path = Path(path).resolve()
        else:
            self.path = Path.cwd()
        
        # Initialize geometry
        if stl is not None:
            if not isinstance(stl, trimesh.Trimesh):
                raise TypeError("stl must be a trimesh.Trimesh object")
            self.stl = stl
        else:
            self.stl = None
        
        # Cached properties for lazy loading
        self._outline2d = None
        self._outline3d = None
        self._buildings = None
        self._face_to_building_map = None
    
    def load(self, filename: str):
        """
        Load geometry from an STL file.
        
        Parameters
        ----------
        filename : str
            Name of the STL file to load (relative to self.path)
            
        Raises
        ------
        FileNotFoundError
            If the STL file does not exist
        ValueError
            If the file cannot be loaded as a valid mesh
            
        Examples
        --------
        >>> geom = UDGeom('experiments/001')
        >>> geom.load('geometry.001')
        """
        filepath = self.path / filename
        
        if not filepath.exists():
            raise FileNotFoundError(f"STL file not found: {filepath}")
        
        try:
            self.stl = trimesh.load_mesh(str(filepath))
            
            # Ensure it's a Trimesh object (not Scene)
            if isinstance(self.stl, trimesh.Scene):
                # Convert scene to single mesh
                self.stl = trimesh.util.concatenate(
                    tuple(trimesh.Trimesh(vertices=g.vertices, faces=g.faces)
                          for g in self.stl.geometry.values())
                )
            
            print(f"Loaded geometry: {len(self.stl.faces)} faces, {len(self.stl.vertices)} vertices")
            
            # Invalidate cached properties
            self._outline2d = None
            self._outline3d = None
            self._buildings = None
            self._face_to_building_map = None
            
        except Exception as e:
            raise ValueError(f"Error loading STL file {filepath}: {e}")
    
    def save(self, filename: str):
        """
        Save geometry to an STL file.
        
        Parameters
        ----------
        filename : str
            Name of the STL file to save (relative to self.path)
            
        Raises
        ------
        ValueError
            If no geometry is loaded
            
        Examples
        --------
        >>> geom.save('modified_geometry.stl')
        """
        if self.stl is None:
            raise ValueError("No geometry loaded. Cannot save.")
        
        filepath = self.path / filename
        
        try:
            self.stl.export(str(filepath))
            print(f"Saved geometry to: {filepath}")
        except Exception as e:
            raise ValueError(f"Error saving STL file {filepath}: {e}")
    
    def show(self, color_buildings: bool = True, plot_quiver: bool = True, 
             normal_scale: float = 0.2, figsize: tuple = (10, 8)):
        """
        Visualize the geometry.
        
        Creates a 3D plot of the triangulated surface. Ground facets (z <= 0) are
        shown in light gray, while building facets (z > 0) can be colored blue.
        Normal vectors can be displayed as arrows.
        
        Parameters
        ----------
        color_buildings : bool, default=True
            If True, color building facets (z > 0) blue. Ground remains gray.
            Set to False for very large geometries to improve performance.
        plot_quiver : bool, default=True
            If True, display normal vectors as arrows from face centers.
            (Renamed from show_normals to match MATLAB interface)
        normal_scale : float, default=0.2
            Scaling factor for normal vector arrow lengths.
        figsize : tuple, default=(10, 8)
            Figure size in inches (width, height).
            
        Raises
        ------
        ValueError
            If no geometry is loaded
        ImportError
            If matplotlib is not installed
            
        Examples
        --------
        >>> geom.show()  # Default: color buildings, show normals
        >>> geom.show(color_buildings=False)  # Faster for large meshes
        >>> geom.show(True, False)  # Color buildings but don't show normals
        """
        if self.stl is None:
            raise ValueError("No geometry loaded. Cannot visualize.")
        
        if not MATPLOTLIB_AVAILABLE:
            raise ImportError("matplotlib is required for visualization. Install with: pip install matplotlib")
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        
        # Get face centers and normals
        face_centers = self.stl.triangles_center
        face_normals = self.stl.face_normals
        
        # Separate ground and building facets
        is_building = face_centers[:, 2] > 0
        
        if color_buildings:
            # Plot ground facets (light gray)
            if np.any(~is_building):
                ground_faces = self.stl.faces[~is_building]
                ground_verts = self.stl.vertices
                ground_triangles = ground_verts[ground_faces]
                
                ground_collection = Poly3DCollection(
                    ground_triangles,
                    facecolors=[0.85, 0.85, 0.85],
                    edgecolors='none',
                    alpha=1.0
                )
                ax.add_collection3d(ground_collection)
            
            # Plot building facets (blue)
            if np.any(is_building):
                building_faces = self.stl.faces[is_building]
                building_verts = self.stl.vertices
                building_triangles = building_verts[building_faces]
                
                building_collection = Poly3DCollection(
                    building_triangles,
                    facecolors=[0.73, 0.83, 0.96],
                    edgecolors='none',
                    alpha=1.0
                )
                ax.add_collection3d(building_collection)
        else:
            # Plot all facets in single color (faster for large geometries)
            all_triangles = self.stl.vertices[self.stl.faces]
            
            collection = Poly3DCollection(
                all_triangles,
                facecolors=[0.85, 0.85, 0.85],
                edgecolors='none',
                alpha=1.0
            )
            ax.add_collection3d(collection)
        
        # Plot normal vectors
        if plot_quiver:
            ax.quiver(
                face_centers[:, 0], face_centers[:, 1], face_centers[:, 2],
                face_normals[:, 0], face_normals[:, 1], face_normals[:, 2],
                length=normal_scale,
                normalize=True,
                color='red',
                arrow_length_ratio=0.3,
                linewidth=0.5
            )
        
        # Set axis properties
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title(f'Geometry: {len(self.stl.faces)} facets')
        
        # Equal aspect ratio
        vertices = self.stl.vertices
        max_range = np.array([
            vertices[:, 0].max() - vertices[:, 0].min(),
            vertices[:, 1].max() - vertices[:, 1].min(),
            vertices[:, 2].max() - vertices[:, 2].min()
        ]).max() / 2.0
        
        mid_x = (vertices[:, 0].max() + vertices[:, 0].min()) * 0.5
        mid_y = (vertices[:, 1].max() + vertices[:, 1].min()) * 0.5
        mid_z = (vertices[:, 2].max() + vertices[:, 2].min()) * 0.5
        
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        # Set viewing angle
        ax.view_init(elev=30, azim=45)
        
        plt.tight_layout()
        plt.show()
    
    @property
    def n_faces(self) -> int:
        """Number of triangular faces in the geometry."""
        return len(self.stl.faces) if self.stl is not None else 0
    
    @property
    def n_vertices(self) -> int:
        """Number of vertices in the geometry."""
        return len(self.stl.vertices) if self.stl is not None else 0
    
    @property
    def bounds(self) -> np.ndarray:
        """
        Bounding box of the geometry.
        
        Returns
        -------
        bounds : ndarray, shape (2, 3)
            [[xmin, ymin, zmin], [xmax, ymax, zmax]]
        """
        if self.stl is None:
            return np.array([[0, 0, 0], [0, 0, 0]])
        return self.stl.bounds
    
    @property
    def face_centers(self) -> np.ndarray:
        """
        Centers of all triangular faces.
        
        Returns
        -------
        centers : ndarray, shape (n_faces, 3)
            (x, y, z) coordinates of face centers
        """
        if self.stl is None:
            return np.array([])
        return self.stl.triangles_center

    def get_buildings(self) -> List[trimesh.Trimesh]:
        """
        Get individual building components with lazy loading.
        
        Separates the geometry into individual disconnected building components.
        Buildings are identified as faces with z > 0 and separated using
        connected component analysis.
        
        Returns
        -------
        buildings : list of trimesh.Trimesh
            List of individual building meshes
            
        Examples
        --------
        >>> buildings = geom.get_buildings()
        >>> num_buildings = len(buildings)
        >>> print(f"Found {num_buildings} buildings")
        """
        if self.stl is None:
            warnings.warn('No STL geometry loaded. Load geometry first.')
            return []
        
        # Lazy load buildings if not already computed
        if self._buildings is None:
            self._split_buildings()
        
        return self._buildings
    
    def get_face_to_building_map(self) -> np.ndarray:
        """
        Get mapping from face indices to building IDs.
        
        Returns
        -------
        face_map : ndarray
            Array where face_map[i] is the building ID for face i.
            Ground faces (z <= 0) are assigned building ID 0.
            Building IDs start from 1.
            
        Examples
        --------
        >>> face_map = geom.get_face_to_building_map()
        >>> building_1_faces = np.where(face_map == 1)[0]
        """
        if self.stl is None:
            warnings.warn('No STL geometry loaded. Load geometry first.')
            return np.array([])
        
        # Ensure buildings and mapping are computed
        if self._face_to_building_map is None:
            self.get_buildings()  # This will compute both buildings and mapping
        
        return self._face_to_building_map
    
    def _split_buildings(self):
        """
        Internal method to split geometry into individual buildings.
        Updates _buildings and _face_to_building_map.
        """
        centers = self.stl.triangles_center
        building_faces = np.where(centers[:, 2] > 0)[0]
        
        if building_faces.size == 0:
            self._buildings = []
            self._face_to_building_map = np.zeros(len(self.stl.faces), dtype=int)
            return
        
        # Build adjacency graph restricted to building faces
        face_adj = self.stl.face_adjacency
        # Keep only adjacency where both faces are buildings
        mask = np.isin(face_adj, building_faces).all(axis=1)
        face_adj = face_adj[mask]
        
        # Connected components over building faces
        component_labels = -np.ones(len(self.stl.faces), dtype=int)
        current_label = 0
        
        for f in building_faces:
            if component_labels[f] != -1:
                continue
            stack = [f]
            component_labels[f] = current_label
            while stack:
                face_idx = stack.pop()
                # neighbors where face_idx participates
                neighbors = face_adj[(face_adj[:, 0] == face_idx) | (face_adj[:, 1] == face_idx)].ravel()
                for nb in neighbors:
                    if component_labels[nb] == -1 and nb in building_faces:
                        component_labels[nb] = current_label
                        stack.append(nb)
            current_label += 1
        
        # Extract individual building meshes
        buildings = []
        # Building IDs start from 1 (0 is reserved for ground)
        face_to_building_map = np.zeros(len(self.stl.faces), dtype=int)
        
        for comp in range(current_label):
            face_idxs = np.where(component_labels == comp)[0]
            if face_idxs.size == 0:
                continue
            
            # Extract submesh for this building
            try:
                building_mesh = self.stl.submesh([face_idxs], append=True)
                buildings.append(building_mesh)
                # Assign building ID (1-indexed)
                face_to_building_map[face_idxs] = len(buildings)
            except Exception as e:
                warnings.warn(f"Could not extract building {comp}: {e}")
        
        self._buildings = buildings
        self._face_to_building_map = face_to_building_map
        
        # Invalidate cached outlines when buildings are recomputed
        self._outline2d = None
        self._outline3d = None
    
    def calculate_outline2d(self) -> List[dict]:
        """
        Calculate 2D building outlines and centroids for all buildings.
        
        Uses lazy loading - returns cached results if available.

        Returns
        -------
        outlines : list of dict
            Each item contains:
            - 'polygon': ndarray, shape (n_vertices, 3) ordered convex hull in x-y
            - 'centroid': ndarray, shape (3,) centroid of the polygon (z=0)
            
        Examples
        --------
        >>> outlines = geom.calculate_outline2d()
        >>> for i, outline in enumerate(outlines):
        ...     print(f"Building {i+1}: centroid at {outline['centroid'][:2]}")
        """
        # Return cached results if available
        if self._outline2d is not None:
            return self._outline2d
        
        if self.stl is None:
            return []
        if not SCIPY_AVAILABLE:
            warnings.warn("scipy is required for calculate_outline2d; install scipy to enable building outlines.")
            return []
        
        # Get individual buildings
        buildings = self.get_buildings()
        if not buildings:
            return []
        
        outlines: List[dict] = []
        for building in buildings:
            verts = building.vertices
            if len(verts) < 3:
                outlines.append({'polygon': np.array([]), 'centroid': np.array([np.nan, np.nan, 0.0])})
                continue
            
            verts_xy = np.unique(verts[:, :2], axis=0)
            if len(verts_xy) < 3:
                centroid = np.array([verts[:, 0].mean(), verts[:, 1].mean(), 0.0])
                outlines.append({'polygon': np.array([]), 'centroid': centroid})
                continue

            try:
                hull = ConvexHull(verts_xy)
                polygon_xy = verts_xy[hull.vertices]
            except Exception:
                # Degenerate (nearly colinear); fallback to all points
                polygon_xy = verts_xy

            # Ensure closed polygon ordering
            polygon_xy = np.vstack([polygon_xy, polygon_xy[0]])
            polygon = np.column_stack([polygon_xy[:, 0], polygon_xy[:, 1], np.zeros(len(polygon_xy))])

            # Polygon centroid in 2D using shoelace; fallback to mean if degenerate
            x = polygon_xy[:, 0]
            y = polygon_xy[:, 1]
            cross = x[:-1] * y[1:] - x[1:] * y[:-1]
            area = 0.5 * np.sum(cross)
            if np.isclose(area, 0):
                centroid_xy = polygon_xy.mean(axis=0)
            else:
                cx = np.sum((x[:-1] + x[1:]) * cross) / (6 * area)
                cy = np.sum((y[:-1] + y[1:]) * cross) / (6 * area)
                centroid_xy = np.array([cx, cy])

            centroid = np.array([centroid_xy[0], centroid_xy[1], 0.0])
            outlines.append({'polygon': polygon, 'centroid': centroid})

        # Sort outlines by y then x (bottom-left to top-right) for consistent IDs
        outlines.sort(key=lambda o: (o['centroid'][1], o['centroid'][0]))
        
        # Cache the result
        self._outline2d = outlines
        
        return outlines
    
    @property
    def face_normals(self) -> np.ndarray:
        """
        Normal vectors for all faces.
        
        Returns
        -------
        normals : ndarray, shape (n_faces, 3)
            Unit normal vectors for each face
        """
        if self.stl is None:
            return np.array([])
        return self.stl.face_normals
    
    @property
    def face_areas(self) -> np.ndarray:
        """
        Areas of all triangular faces.
        
        Returns
        -------
        areas : ndarray, shape (n_faces,)
            Area of each face in m²
        """
        if self.stl is None:
            return np.array([])
        return self.stl.area_faces
    
    @property
    def total_area(self) -> float:
        """Total surface area of the geometry in m²."""
        if self.stl is None:
            return 0.0
        return self.stl.area
    
    @property
    def volume(self) -> float:
        """
        Volume enclosed by the geometry in m³.
        
        Note: Only valid for watertight meshes.
        """
        if self.stl is None:
            return 0.0
        if not self.stl.is_watertight:
            warnings.warn("Mesh is not watertight. Volume may be incorrect.")
        return self.stl.volume
    
    @property
    def is_watertight(self) -> bool:
        """True if the mesh is watertight (closed, no holes)."""
        if self.stl is None:
            return False
        return self.stl.is_watertight
    
    def show_outline(self, angle_threshold: float = 45.0, figsize: tuple = (10, 8)):
        """
        Plot the geometry with outline edges highlighted.
        
        Displays the mesh in gray with black outline edges overlaid.
        Outline edges are detected based on the angle between adjacent faces.
        
        Parameters
        ----------
        angle_threshold : float, default=45.0
            Angle threshold in degrees for edge detection.
            Edges where adjacent faces meet at angles greater than this
            threshold are considered outline edges.
        figsize : tuple, default=(10, 8)
            Figure size in inches (width, height).
            
        Raises
        ------
        ValueError
            If no geometry is loaded
        ImportError
            If matplotlib is not installed
            
        Examples
        --------
        >>> geom.show_outline()  # Default 45° threshold
        >>> geom.show_outline(angle_threshold=30)  # More sensitive
        """
        if self.stl is None:
            raise ValueError("No geometry loaded. Cannot visualize.")
        
        if not MATPLOTLIB_AVAILABLE:
            raise ImportError("matplotlib is required for visualization. Install with: pip install matplotlib")
        
        # Get outline edges
        outline_edges = self._calculate_outline_edges(angle_threshold)
        
        if len(outline_edges) == 0:
            warnings.warn('No outline edges found.')
            return
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot mesh without edges
        all_triangles = self.stl.vertices[self.stl.faces]
        collection = Poly3DCollection(
            all_triangles,
            facecolors=[0.85, 0.85, 0.85],
            edgecolors='none',
            alpha=0.8
        )
        ax.add_collection3d(collection)
        
        # Plot outline edges (including ground facet edges)
        for edge in outline_edges:
            v1 = self.stl.vertices[edge[0]]
            v2 = self.stl.vertices[edge[1]]
            ax.plot([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]], 'k-', linewidth=1)
        
        # Set axis properties
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title(f'Geometry Outline ({len(outline_edges)} edges)')
        
        # Equal aspect ratio
        vertices = self.stl.vertices
        max_range = np.array([
            vertices[:, 0].max() - vertices[:, 0].min(),
            vertices[:, 1].max() - vertices[:, 1].min(),
            vertices[:, 2].max() - vertices[:, 2].min()
        ]).max() / 2.0
        
        mid_x = (vertices[:, 0].max() + vertices[:, 0].min()) * 0.5
        mid_y = (vertices[:, 1].max() + vertices[:, 1].min()) * 0.5
        mid_z = (vertices[:, 2].max() + vertices[:, 2].min()) * 0.5
        
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        ax.set_box_aspect([1, 1, 1])
        ax.grid(True)
        ax.view_init(elev=30, azim=45)
        
        plt.tight_layout()
        plt.show()
    
    def _calculate_outline_edges(self, angle_threshold: float = 45.0) -> List[tuple]:
        """
        Internal method to calculate outline edges based on face angle threshold.
        
        Matches MATLAB get_outline() behavior: processes the entire geometry including
        ground faces. This ensures boundary edges where buildings meet the ground are
        properly detected as outline edges.
        
        Parameters
        ----------
        angle_threshold : float
            Angle threshold in degrees
            
        Returns
        -------
        edges : list of tuple
            List of (vertex_idx1, vertex_idx2) tuples defining outline edges
        """
        if self.stl is None:
            return []
        
        # Process the entire geometry (including ground faces) to match MATLAB get_outline()
        # MATLAB only filters ground when splitting buildings, not when computing overall outline
        vertices = self.stl.vertices
        faces = self.stl.faces
        
        # Get face adjacency info
        face_adjacency = self.stl.face_adjacency
        face_adjacency_angles = self.stl.face_adjacency_angles
        
        # Convert angle threshold to radians
        threshold_rad = np.deg2rad(angle_threshold)
        
        # Find edges where angle exceeds threshold (sharp edges)
        sharp_edges_mask = face_adjacency_angles > threshold_rad
        sharp_edge_pairs = face_adjacency[sharp_edges_mask]
        
        outline_edges = []
        outline_edge_set = set()  # Track unique edges
        
        # Add sharp-angle edges (edges between faces with large angle difference)
        for face_pair in sharp_edge_pairs:
            face1, face2 = face_pair
            
            # Find shared vertices between the two faces
            verts1 = set(self.stl.faces[face1])
            verts2 = set(self.stl.faces[face2])
            shared = verts1.intersection(verts2)
            
            if len(shared) == 2:
                # Found the shared edge
                edge = tuple(sorted(shared))
                if edge not in outline_edge_set:
                    outline_edge_set.add(edge)
                    outline_edges.append(edge)
        
        # Add boundary edges (edges that belong to only one face)
        # Build edge-to-face mapping
        edge_counts = {}
        for face_idx in range(len(faces)):
            face_verts = faces[face_idx]
            face_edges = [
                tuple(sorted([face_verts[0], face_verts[1]])),
                tuple(sorted([face_verts[1], face_verts[2]])),
                tuple(sorted([face_verts[2], face_verts[0]]))
            ]
            for edge in face_edges:
                edge_counts[edge] = edge_counts.get(edge, 0) + 1
        
        # Boundary edges appear only once (not shared between faces)
        boundary_edges_added = 0
        ground_boundary_edges = 0
        for edge, count in edge_counts.items():
            if count == 1:  # Boundary edge
                if edge not in outline_edge_set:
                    outline_edge_set.add(edge)
                    outline_edges.append(edge)
                    boundary_edges_added += 1
                    # Check if this is a ground-level edge
                    if vertices[edge[0]][2] == 0 and vertices[edge[1]][2] == 0:
                        ground_boundary_edges += 1
        
        print(f"Outline calculation: {len(outline_edges)} total edges ({len(sharp_edge_pairs)} sharp, {boundary_edges_added} boundary, {ground_boundary_edges} at ground level)")
        
        # Debug: Show sample of ground-level edges
        if ground_boundary_edges > 0:
            ground_edges_sample = [(edge, vertices[edge[0]], vertices[edge[1]]) 
                                   for edge in outline_edges 
                                   if vertices[edge[0]][2] == 0 and vertices[edge[1]][2] == 0][:3]
            print(f"Sample ground edges: {[(e[0], f'({e[1][0]:.1f},{e[1][1]:.1f},{e[1][2]:.1f})', f'({e[2][0]:.1f},{e[2][1]:.1f},{e[2][2]:.1f})') for e in ground_edges_sample]}")
        
        return outline_edges
    
    def __repr__(self) -> str:
        """String representation of the UDGeom object."""
        if self.stl is None:
            return f"UDGeom(path={self.path}, stl=None)"
        return (f"UDGeom(path={self.path}, "
                f"n_faces={self.n_faces}, n_vertices={self.n_vertices})")
