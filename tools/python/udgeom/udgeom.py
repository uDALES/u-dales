"""
udgeom - Geometry class for uDALES

This module provides the UDGeom class for handling triangulated surface geometries
used in uDALES simulations. It wraps trimesh functionality for STL file I/O and
visualization.

uDALES (https://github.com/uDALES/u-dales)

Copyright (C) 2016-2024 the uDALES Team.
Licensed under GNU General Public License v3.0
"""

from __future__ import annotations

from pathlib import Path
import struct
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
    from scipy.spatial import ConvexHull
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

# Import package functions
from .calculate_outline import calculate_outline
from .check_mesh import (
    calculate_independent_surfaces as calculate_independent_surfaces_impl,
    check as check_mesh,
    find_internal_touching_wall_regions as find_internal_touching_wall_regions_impl,
    find_nonmanifold_regions as find_nonmanifold_regions_impl,
    find_unstitched_touching_regions as find_unstitched_touching_regions_impl,
    identify_ground_faces as identify_ground_faces_impl,
)
from .fix_mesh import (
    fix as fix_mesh,
    repair_adjacent_buildings as repair_adjacent_buildings_impl,
    resolve_vertical_coplanar_overlaps as resolve_vertical_coplanar_overlaps_impl,
    weld_touching_boundaries as weld_touching_boundaries_impl,
)
from .split_buildings import split_buildings

from udvis import UDVis, DEFAULT_BACKEND


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
                 stl: Optional['trimesh.Trimesh'] = None,
                 backend: str = DEFAULT_BACKEND):
        """
        Initialize UDGeom object.

        Parameters
        ----------
        path : str or Path, optional
            Path to geometry directory. If None, uses current directory.
        stl : trimesh.Trimesh, optional
            Pre-existing trimesh object.
        backend : {"plotly", "pyvista"}, default ``udvis.DEFAULT_BACKEND``
            Default rendering backend for the 3-D geometry plots. Overridable
            per instance via ``geom.backend`` or per call via ``backend=``.
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
        
        # Cached lazy-loaded properties (see _invalidate_cache).
        self._invalidate_cache()
        self.vis = UDVis(self, backend=backend)

    def _invalidate_cache(self) -> None:
        """Reset all lazily-computed geometry caches (outlines, buildings,
        face->building map). Call after any mutation of ``self.stl``."""
        self._outline2d = None
        self._outline3d = None
        self._buildings = None
        self._face_to_building_map = None

    @property
    def backend(self) -> str:
        """Default rendering backend for 3-D plots (``"plotly"`` or ``"pyvista"``).

        Settable on an existing geometry, so a geometry produced by a factory
        (``create_cubes`` ...) or a repair step can be switched without passing
        ``backend=`` to every ``show``/``plot`` call::

            geom = create_cubes(...)
            geom.backend = "pyvista"
            geom.show()
        """
        return self.vis.backend

    @backend.setter
    def backend(self, value: str) -> None:
        self.vis.backend = value

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

            self._align_loaded_stl_to_file_normals(filepath)
            
            print(f"Loaded geometry: {len(self.stl.faces)} faces, {len(self.stl.vertices)} vertices")
            
            # Invalidate cached properties
            self._invalidate_cache()
            
        except Exception as e:
            raise ValueError(f"Error loading STL file {filepath}: {e}")

    @staticmethod
    def _read_stl_file_normals(filepath: Path) -> Optional[np.ndarray]:
        """Read facet normals stored in an STL file, if available."""
        if filepath.suffix.lower() != ".stl":
            return None

        try:
            file_size = filepath.stat().st_size
            with filepath.open("rb") as f:
                header = f.read(84)
                if len(header) == 84:
                    n_faces = struct.unpack("<I", header[80:84])[0]
                    if file_size == 84 + 50 * n_faces:
                        record_dtype = np.dtype(
                            [
                                ("normal", "<f4", (3,)),
                                ("vertices", "<f4", (3, 3)),
                                ("attribute", "<u2"),
                            ]
                        )
                        records = np.fromfile(f, dtype=record_dtype, count=n_faces)
                        if len(records) == n_faces:
                            return np.asarray(records["normal"], dtype=float)

            normals = []
            with filepath.open("r", encoding="ascii", errors="ignore") as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 5 and parts[0].lower() == "facet" and parts[1].lower() == "normal":
                        normals.append([float(parts[2]), float(parts[3]), float(parts[4])])
            if normals:
                return np.asarray(normals, dtype=float)
        except (OSError, UnicodeError, ValueError, struct.error):
            return None
        return None

    def _align_loaded_stl_to_file_normals(self, filepath: Path) -> None:
        """Match in-memory STL face winding to normals stored in the STL file."""
        if self.stl is None:
            return

        file_normals = self._read_stl_file_normals(filepath)
        if file_normals is None or len(file_normals) != len(self.stl.faces):
            return

        lengths = np.linalg.norm(file_normals, axis=1)
        valid = np.isfinite(file_normals).all(axis=1) & (lengths > 0.0)
        if not np.any(valid):
            return

        normalized_file_normals = np.zeros_like(file_normals, dtype=float)
        normalized_file_normals[valid] = file_normals[valid] / lengths[valid, None]
        loaded_normals = np.asarray(self.stl.face_normals, dtype=float)
        dots = np.einsum("ij,ij->i", loaded_normals, normalized_file_normals)
        flip = valid & np.isfinite(dots) & (dots < 0.0)
        if not np.any(flip):
            return

        faces = np.asarray(self.stl.faces, dtype=int).copy()
        faces[flip, 1], faces[flip, 2] = faces[flip, 2].copy(), faces[flip, 1].copy()
        self.stl.faces = faces
    
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
    
    def show(self, color_buildings: bool = True, plot_quiver: bool = False,
             normal_scale: float = 0.2,
             show_edges: bool = True, show_ground: bool = True, show: bool = True,
             backend: Optional[str] = None):
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
        plot_quiver : bool, default=False
            If True, display normal vectors as arrows from face centers.
            (Renamed from show_normals to match the MATLAB interface; note the
            MATLAB toolbox defaults this to true, whereas Python defaults to
            False to avoid cluttering large meshes.)
        normal_scale : float, default=0.2
            Scaling factor for normal vector arrow lengths.
        show_edges : bool, default=True
            If True, draw facet edges; set False for cleaner surfaces.
        show_ground : bool, default=True
            If True, include z=0 ground facets in the geometry view.
        show : bool, default=True
            If True, display the figure immediately. If False, only return the
            figure object.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to the backend chosen when this UDGeom
            was constructed (``UDGeom(..., backend=...)``).

        Returns
        -------
        plotly.graph_objects.Figure or pyvista.Plotter or None
            When ``show=False``, the Plotly ``Figure`` or the PyVista
            ``Plotter`` depending on ``backend``; ``None`` when ``show=True``
            (the figure/window is displayed instead).

        Raises
        ------
        ValueError
            If no geometry is loaded
        ImportError
            If trimesh (plotly backend) or pyvista (pyvista backend) is not
            installed

        Examples
        --------
        >>> geom.show()  # uses the geometry's default backend
        >>> geom.show(backend="pyvista")  # override for this call
        >>> geom.show(color_buildings=False)  # Faster for large meshes
        >>> fig = geom.show(show=False)  # Build the figure without displaying it
        """
        return self.vis.show_geometry(
            color_buildings=color_buildings,
            plot_quiver=plot_quiver,
            normal_scale=normal_scale,
            show_edges=show_edges,
            show_ground=show_ground,
            show=show,
            backend=backend,
        )
    
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
    def face_incenters(self) -> np.ndarray:
        """
        Incenters of all triangular faces.

        Returns
        -------
        incenters : ndarray, shape (n_faces, 3)
            Incenter coordinates for each triangular face.
        """
        if self.stl is None:
            return np.array([])

        vertices = np.asarray(self.stl.vertices, dtype=float)
        faces = np.asarray(self.stl.faces, dtype=int)
        if faces.size == 0:
            return np.empty((0, 3), dtype=float)

        triangles = vertices[faces]
        edge_a = np.linalg.norm(triangles[:, 1] - triangles[:, 2], axis=1)
        edge_b = np.linalg.norm(triangles[:, 0] - triangles[:, 2], axis=1)
        edge_c = np.linalg.norm(triangles[:, 0] - triangles[:, 1], axis=1)
        perimeter = edge_a + edge_b + edge_c

        incenters = (
            edge_a[:, None] * triangles[:, 0]
            + edge_b[:, None] * triangles[:, 1]
            + edge_c[:, None] * triangles[:, 2]
        )
        valid = perimeter > 0.0
        if np.any(valid):
            incenters[valid] = incenters[valid] / perimeter[valid, None]
        if np.any(~valid):
            count = int(np.count_nonzero(~valid))
            warnings.warn(
                f"{count} fully collapsed triangle(s) encountered while computing face incenters; "
                "using triangle centroids for those faces.",
                stacklevel=2,
            )
            incenters[~valid] = np.mean(triangles[~valid], axis=1)
        return incenters

    def get_buildings(self) -> List[trimesh.Trimesh]:
        """
        Get individual building components with lazy loading.
        
        Separates the geometry into individual disconnected building components.
        Buildings are identified as faces with z > 0 and separated using
        connected component analysis. Buildings are ordered by their centroid
        from southwest to northeast (ascending x+y projection).
        
        Returns
        -------
        buildings : list of trimesh.Trimesh
            List of individual building meshes. Each mesh has a `.points` 
            property (alias for `.vertices`) for MATLAB compatibility.
            
        Examples
        --------
        >>> buildings = geom.get_buildings()
        >>> num_buildings = len(buildings)
        >>> for i, bld in enumerate(buildings):
        ...     max_height = np.max(bld.points[:, 2])
        ...     print(f"Building {i+1}: max height = {max_height:.1f} m")
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
        Buildings are sorted by centroid from southwest to northeast.
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
        
        # Extract individual building meshes with component labels
        buildings_with_labels = []
        
        for comp in range(current_label):
            face_idxs = np.where(component_labels == comp)[0]
            if face_idxs.size == 0:
                continue
            
            # Extract submesh for this building
            try:
                building_mesh = self.stl.submesh([face_idxs], append=True)
                # Add .points property as alias for .vertices (MATLAB compatibility)
                building_mesh.points = building_mesh.vertices
                buildings_with_labels.append((building_mesh, face_idxs, comp))
            except Exception as e:
                warnings.warn(f"Could not extract building {comp}: {e}")
        
        # Sort buildings by centroid (southwest to northeast: ascending x+y)
        if buildings_with_labels:
            centroids = np.array([bld[0].vertices[:, :2].mean(axis=0) for bld in buildings_with_labels])
            proj = centroids[:, 0] + centroids[:, 1]  # x + y projection
            # Match MATLAB's ascending sort on x+y while preserving the
            # original connected-component order for ties.
            order = np.argsort(proj, kind="stable")
            
            # Reorder buildings
            buildings_with_labels = [buildings_with_labels[i] for i in order]
            
            # Create new component ID mapping (buildings now 1-indexed)
            buildings = []
            face_to_building_map = np.zeros(len(self.stl.faces), dtype=int)
            
            for new_id, (building_mesh, face_idxs, old_comp) in enumerate(buildings_with_labels, start=1):
                buildings.append(building_mesh)
                face_to_building_map[face_idxs] = new_id
            
            self._buildings = buildings
            self._face_to_building_map = face_to_building_map
        else:
            self._buildings = []
            self._face_to_building_map = np.zeros(len(self.stl.faces), dtype=int)
        
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
        
        # Get individual buildings
        buildings = self.get_buildings()
        if not buildings:
            return []
        
        outlines: List[dict] = []
        for building in buildings:
            verts = np.asarray(building.vertices, dtype=float)
            faces = np.asarray(building.faces, dtype=int)
            if len(verts) < 3 or len(faces) == 0:
                outlines.append({'polygon': np.array([]), 'centroid': np.array([np.nan, np.nan, 0.0])})
                continue

            pts2 = verts[:, :2]
            if len(pts2) < 3:
                centroid = np.array([verts[:, 0].mean(), verts[:, 1].mean(), 0.0])
                outlines.append({'polygon': np.array([]), 'centroid': centroid})
                continue

            try:
                edge_counts = {}
                directed_boundary_edges = []
                for face in faces:
                    directed_edges = (
                        (int(face[0]), int(face[1])),
                        (int(face[1]), int(face[2])),
                        (int(face[2]), int(face[0])),
                    )
                    for start, end in directed_edges:
                        edge = tuple(sorted((start, end)))
                        edge_counts[edge] = edge_counts.get(edge, 0) + 1
                    directed_boundary_edges.extend(directed_edges)

                boundary_directed = [
                    edge for edge in directed_boundary_edges
                    if edge_counts[tuple(sorted(edge))] == 1
                ]

                if boundary_directed:
                    next_vertices = {}
                    for start, end in boundary_directed:
                        next_vertices.setdefault(start, []).append(end)

                    polygon_indices = []
                    remaining = {start: list(ends) for start, ends in next_vertices.items()}

                    while remaining:
                        start = min(remaining)
                        current = start

                        while True:
                            polygon_indices.append(current)
                            ends = remaining.get(current)
                            if not ends:
                                break

                            next_vertex = ends.pop(0)
                            if not ends:
                                remaining.pop(current, None)

                            current = next_vertex
                            if current == start:
                                break

                    # Match MATLAB's `tri2.Points(boundary_edges(:,1), :)` by
                    # preserving one ordered start vertex per boundary edge over
                    # every loop in the projected free boundary.
                    polygon = pts2[np.asarray(polygon_indices, dtype=int)]
                    centroid_xy = polygon.mean(axis=0)
                else:
                    polygon = np.array([])
                    centroid_xy = pts2.mean(axis=0)
            except Exception:
                polygon = np.array([])
                centroid_xy = pts2.mean(axis=0)

            if polygon.size == 0:
                polygon3d = np.array([])
            else:
                polygon3d = np.column_stack(
                    [polygon[:, 0], polygon[:, 1], np.zeros(len(polygon))]
                )

            centroid = np.array([centroid_xy[0], centroid_xy[1], 0.0])
            outlines.append({'polygon': polygon3d, 'centroid': centroid})

        # Cache the result
        self._outline2d = outlines
        
        return outlines

    def get_building_outlines(self, angle_threshold: float = 45.0) -> List[np.ndarray]:
        """
        Get 3D outline edges for each disconnected building component.

        Parameters
        ----------
        angle_threshold : float, default=45.0
            Angle threshold in degrees for edge detection.

        Returns
        -------
        outlines : list of ndarray
            Per-building edge arrays with shape (n_edges, 2). Edge indices refer
            to each building submesh, matching the MATLAB method contract.
        """
        if self.stl is None:
            warnings.warn('No STL geometry loaded. Load geometry first.')
            return []

        if self._outline3d is None:
            self._outline3d = {}

        cache_key = float(angle_threshold)
        if cache_key not in self._outline3d:
            outlines: List[np.ndarray] = []
            for building in self.get_buildings():
                edges, _ = calculate_outline(building, angle_threshold)
                outlines.append(edges)
            if self._outline3d is None:
                self._outline3d = {}
            self._outline3d[cache_key] = outlines
        return self._outline3d[cache_key]

    def get_outline(self, angle_threshold: float = 45.0) -> np.ndarray:
        """
        Get 3D outline edges for the entire geometry.

        Parameters
        ----------
        angle_threshold : float, default=45.0
            Angle threshold in degrees for edge detection.

        Returns
        -------
        outline_edges : ndarray, shape (n_edges, 2)
            Edge vertex pairs referring to ``self.stl.vertices``.
        """
        if self.stl is None:
            warnings.warn('No STL geometry loaded. Load geometry first.')
            return np.empty((0, 2), dtype=int)

        edges, _ = calculate_outline(self.stl, angle_threshold)
        return edges
    
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
    
    def show_outline(
        self,
        angle_threshold: float = 45.0,
        show_ground: bool = True,
        color_buildings: bool = False,
        show: bool = True,
        backend: Optional[str] = None,
    ):
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
        show_ground : bool, default=True
            If True, include ground faces in the visualization.
        color_buildings : bool, default=False
            If True, colour buildings blue and ground grey. If False (default),
            use two-tone grey shading that makes the outline edges stand out.
            Honoured by both backends.
        show : bool, default=True
            If True, display the figure immediately. If False, only return the
            figure object.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to the backend chosen when this UDGeom
            was constructed (``UDGeom(..., backend=...)``).

        Raises
        ------
        ValueError
            If no geometry is loaded
        ImportError
            If trimesh (plotly backend) or pyvista (pyvista backend) is not
            installed

        Examples
        --------
        >>> geom.show_outline()  # uses the geometry's default backend
        >>> geom.show_outline(backend="pyvista")  # override for this call
        >>> geom.show_outline(angle_threshold=30)  # More sensitive
        """
        return self.vis.show_geometry_outline(
            angle_threshold=angle_threshold,
            show_ground=show_ground,
            color_buildings=color_buildings,
            show=show,
            backend=backend,
        )

    def check(self, require_single_component: bool = True):
        """
        Validate the current mesh for common udgeom mesh-quality issues.

        Use this as the first step in a repair workflow. The returned report
        contains both summary counts and detailed defect locations.

        Parameters
        ----------
        require_single_component : bool, default=True
            If True, the mesh must form a single face-connected component.

        Returns
        -------
        report : dict
            Validation report from `udgeom.check`.
        """
        return check_mesh(self, require_single_component=require_single_component)

    def add_ground(
        self,
        xsize: float,
        ysize: float,
        edgelength: float,
        *,
        shift: tuple[float, float, float] = (0.0, 0.0, 0.0),
        preserve_existing_edges: bool = True,
        return_debug: bool = False,
        inplace: bool = False,
    ):
        """
        Add a flat ground surface around the current geometry.

        This is the public wrapper for stitching a flat ground mesh around an
        existing building shell or imported STL geometry.
        """
        from .geometry_generation import add_ground as add_ground_impl

        result = add_ground_impl(
            self,
            xsize,
            ysize,
            shift=shift,
            edgelength=edgelength,
            return_debug=return_debug,
            preserve_existing_edges=preserve_existing_edges,
        )
        if return_debug:
            geom, debug = result
            if inplace:
                self.stl = geom.stl
                self._outline2d = None
                self._outline3d = None
                self._buildings = None
                self._face_to_building_map = None
            return geom, debug

        geom = result
        if inplace:
            self.stl = geom.stl
            self._invalidate_cache()
        return geom

    def calculate_independent_surfaces(self):
        """
        Partition the mesh into face-connected independent surfaces.

        This is useful for diagnosing disconnected geometry and for selecting
        one surface id for a targeted repair such as `extrude_to_ground(...)`.

        Returns
        -------
        result : dict
            Surface partition with ``n_surfaces``, ``face_surface_ids``, and
            per-surface face-id lists.
        """
        return calculate_independent_surfaces_impl(self)

    def identify_ground_faces(
        self,
        max_slope_deg: float = 20.0,
        min_component_area_fraction: float = 0.05,
    ) -> np.ndarray:
        """
        Identify likely ground faces without assuming a fixed absolute z level.

        Parameters
        ----------
        max_slope_deg : float, default=20.0
            Maximum slope from horizontal for ground candidates.
        min_component_area_fraction : float, default=0.05
            Minimum area fraction for a horizontal component to be considered.

        Returns
        -------
        mask : ndarray of bool
            Boolean mask over faces identifying likely ground facets.
        """
        return identify_ground_faces_impl(
            self,
            max_slope_deg=max_slope_deg,
            min_component_area_fraction=min_component_area_fraction,
        )

    def find_nonmanifold_regions(self):
        """
        Locate clustered non-manifold defect regions.

        Returns
        -------
        regions : list of dict
            Clustered non-manifold regions with edge ids, face ids, and bounding boxes.
        """
        return find_nonmanifold_regions_impl(self)

    def find_internal_touching_wall_regions(self):
        """
        Locate coplanar opposite-facing internal wall overlaps.

        Returns
        -------
        regions : list of dict
            Overlap regions with face ids, plane location, overlap area, and bounding boxes.
        """
        return find_internal_touching_wall_regions_impl(self)

    def find_unstitched_touching_regions(self):
        """
        Locate touching seams that are geometrically coincident but not stitched.

        Returns
        -------
        regions : list of dict
            Seam regions with edge ids, face ids, and bounding boxes.
        """
        return find_unstitched_touching_regions_impl(self)

    def fix(
        self,
        *,
        merge_tolerance: float | None = None,
        remove_small_components: bool = False,
        min_component_faces: int = 0,
        min_component_area_fraction: float = 0.0,
        resolve_vertical_coplanar_overlaps: bool = False,
        inplace: bool = False,
    ):
        """
        Apply conservative small-scale mesh cleanup.

        This method is intended for low-risk housekeeping operations such as
        removing duplicate faces, zero-area triangles, unused vertices, and
        optionally tiny disconnected fragments. For larger geometric edits, use
        the dedicated repair helpers.

        Parameters
        ----------
        merge_tolerance : float, optional
            Merge nearby vertices after rounding to this tolerance.
        remove_small_components : bool, default=False
            Remove tiny disconnected components.
        min_component_faces : int, default=0
            Minimum face count to keep when removing small components.
        min_component_area_fraction : float, default=0.0
            Minimum component area fraction to keep when removing small components.
        resolve_vertical_coplanar_overlaps : bool, default=False
            Also resolve opposite-facing coplanar vertical wall overlaps.
            Supported here for convenience; the dedicated
            ``resolve_vertical_coplanar_overlaps()`` method is usually clearer.
        inplace : bool, default=False
            If True, replace ``self.stl`` with the fixed mesh.

        Returns
        -------
        fixed_mesh, report
            The cleaned mesh plus a change report.
        """
        fixed_mesh, report = fix_mesh(
            self,
            merge_tolerance=merge_tolerance,
            remove_small_components=remove_small_components,
            min_component_faces=min_component_faces,
            min_component_area_fraction=min_component_area_fraction,
            resolve_vertical_coplanar_overlaps=resolve_vertical_coplanar_overlaps,
        )
        if inplace:
            self.stl = fixed_mesh
            self._invalidate_cache()
        return fixed_mesh, report

    def resolve_vertical_coplanar_overlaps(
        self,
        *,
        weld_touching_boundaries: bool = True,
        inplace: bool = False,
    ):
        """
        Resolve adjacent-building overlap defects on vertical shared planes.

        This helper removes internal shared wall patches where two shells both
        contain faces on the same vertical plane, then optionally welds any new
        touching boundary seams created by the repair.
        """
        fixed_mesh, report = resolve_vertical_coplanar_overlaps_impl(
            self,
            weld_touching_boundaries=weld_touching_boundaries,
        )
        if inplace:
            self.stl = fixed_mesh
            self._invalidate_cache()
        return fixed_mesh, report

    def weld_touching_boundaries(
        self,
        *,
        inplace: bool = False,
    ):
        """
        Weld touching but unstitched mesh boundaries.

        This helper targets adjacency defects where two surfaces touch
        geometrically but do not yet share the same boundary discretization.
        """
        fixed_mesh, report = weld_touching_boundaries_impl(self)
        if inplace:
            self.stl = fixed_mesh
            self._invalidate_cache()
        return fixed_mesh, report

    def repair_adjacent_buildings(
        self,
        *,
        merge_digits_vertex: int = 8,
        plane_tolerance: float = 1.0e-8,
        point_tolerance: float = 1.0e-8,
        return_trimesh: bool = False,
        inplace: bool = False,
    ):
        """
        Repair adjacent-building shared-wall defects while preserving original ground.

        This helper keeps the current ground mesh, repairs the non-ground shell
        for vertical coplanar overlap defects, then rejoins and welds the full
        geometry.
        """
        cleaned_geom, report = repair_adjacent_buildings_impl(
            self,
            merge_digits_vertex=merge_digits_vertex,
            plane_tolerance=plane_tolerance,
            point_tolerance=point_tolerance,
            return_trimesh=False,
        )
        if inplace:
            self.stl = cleaned_geom.stl
            self._invalidate_cache()
            return (self.stl, report) if return_trimesh else (self, report)
        return (cleaned_geom.stl, report) if return_trimesh else (cleaned_geom, report)

    def truncate_below_ground(
        self,
        *,
        ground_planarity_tolerance: float = 1.0e-6,
        clip_tolerance: float = 1.0e-9,
        edgelength: float | None = None,
        return_trimesh: bool = False,
        inplace: bool = False,
    ):
        """
        Remove below-ground geometry for planar-ground cases and restitch the ground locally.

        Parameters
        ----------
        ground_planarity_tolerance : float, default=1e-6
            Maximum allowed deviation from a single ground plane.
        clip_tolerance : float, default=1e-9
            Numerical tolerance for triangle-plane clipping.
        edgelength : float, optional
            Optional reference spacing stored in the report. The routine
            preserves the current ground mesh and restitches it locally.
        return_trimesh : bool, default=False
            If True, return the repaired ``trimesh.Trimesh``. Otherwise return
            ``UDGeom``.
        inplace : bool, default=False
            If True, replace ``self.stl`` with the cleaned mesh.

        Returns
        -------
        cleaned_geom, report
            The cleaned geometry plus a truncation report.
        """
        from .truncate_below_ground import truncate_below_ground as truncate_below_ground_impl

        cleaned_geom, report = truncate_below_ground_impl(
            self,
            ground_planarity_tolerance=ground_planarity_tolerance,
            clip_tolerance=clip_tolerance,
            edgelength=edgelength,
            return_trimesh=False,
        )
        if inplace:
            self.stl = cleaned_geom.stl
            self._invalidate_cache()
            return (self.stl, report) if return_trimesh else (self, report)
        return (cleaned_geom.stl, report) if return_trimesh else (cleaned_geom, report)

    def extrude_to_ground(
        self,
        surface_id: int,
        *,
        ground_planarity_tolerance: float = 1.0e-6,
        bottom_vertex_tolerance: float = 1.0e-8,
        return_trimesh: bool = False,
        inplace: bool = False,
    ):
        """
        Drop one independent surface vertically to the identified planar ground.

        This is intended for simple floating-surface cases, such as a building
        shell that is disconnected from the ground and should be extended
        downward to meet it.

        Parameters
        ----------
        surface_id : int
            Independent-surface id from ``calculate_independent_surfaces``.
        ground_planarity_tolerance : float, default=1e-6
            Maximum allowed deviation from a single ground plane.
        bottom_vertex_tolerance : float, default=1e-8
            Tolerance used to identify the selected surface's bottom vertices.
        return_trimesh : bool, default=False
            If True, return the repaired ``trimesh.Trimesh``. Otherwise return
            ``UDGeom``.
        inplace : bool, default=False
            If True, replace ``self.stl`` with the cleaned mesh.

        Returns
        -------
        cleaned_geom, report
            Extruded geometry plus a repair report.
        """
        from .extrude_to_ground import extrude_to_ground as extrude_to_ground_impl

        cleaned_geom, report = extrude_to_ground_impl(
            self,
            surface_id=surface_id,
            ground_planarity_tolerance=ground_planarity_tolerance,
            bottom_vertex_tolerance=bottom_vertex_tolerance,
            return_trimesh=False,
        )
        if inplace:
            self.stl = cleaned_geom.stl
            self._invalidate_cache()
            return (self.stl, report) if return_trimesh else (self, report)
        return (cleaned_geom.stl, report) if return_trimesh else (cleaned_geom, report)

    def plot_independent_surfaces(self, show: bool = True, return_result: bool = False,
                                  backend: Optional[str] = None):
        """
        Visualize face-connected independent surfaces using face ids as colors.

        Parameters
        ----------
        show : bool, default=True
            If True, display the figure immediately.
        return_result : bool, default=False
            If True, also return the independent-surface partition result.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to the backend chosen when this UDGeom
            was constructed (``UDGeom(..., backend=...)``).

        Returns
        -------
        fig or (fig, result)
            The figure/plotter when ``show=False``; ``None`` when ``show=True``
            (displayed instead). When ``return_result=True``, the surface
            partition result is returned alongside, i.e. ``(fig, result)`` or
            ``(None, result)``.
        """
        return self.vis.plot_independent_surfaces(show=show, return_result=return_result, backend=backend)
    
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
        
        return outline_edges
    
    def __repr__(self) -> str:
        """String representation of the UDGeom object."""
        if self.stl is None:
            return f"UDGeom(path={self.path}, stl=None)"
        return (f"UDGeom(path={self.path}, "
                f"n_faces={self.n_faces}, n_vertices={self.n_vertices})")
