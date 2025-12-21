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
from typing import Optional, Union
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
    
    def show(self, color_buildings: bool = True, show_normals: bool = True, 
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
        show_normals : bool, default=True
            If True, display normal vectors as arrows from face centers.
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
        >>> geom.show(show_normals=False, figsize=(12, 10))
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
        if show_normals:
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
    
    def __repr__(self) -> str:
        """String representation of the UDGeom object."""
        if self.stl is None:
            return f"UDGeom(path={self.path}, stl=None)"
        return (f"UDGeom(path={self.path}, "
                f"n_faces={self.n_faces}, n_vertices={self.n_vertices})")
