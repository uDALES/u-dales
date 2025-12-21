"""
Geometry Generation Module for uDALES

This module provides functions to create synthetic geometries for urban
canopy simulations, including flat surfaces, cubes, canyons, and realistic
building configurations.

Copyright (C) 2024 the uDALES Team.
"""

import numpy as np
from typing import Optional, Tuple, Literal
try:
    import trimesh
    TRIMESH_AVAILABLE = True
except ImportError:
    TRIMESH_AVAILABLE = False
    raise ImportError("trimesh is required for geometry generation. "
                     "Install with: pip install trimesh")

from udgeom import UDGeom


def create_flat_surface(xsize: float, ysize: float, 
                        edgelength: Optional[float] = None) -> UDGeom:
    """
    Create a flat ground surface with triangular facets.
    
    Parameters
    ----------
    xsize : float
        Domain length in x-direction [m]
    ysize : float
        Domain length in y-direction [m]
    edgelength : float, optional
        Approximate edge length of triangular facets [m].
        If None, uses xsize / 10.
        
    Returns
    -------
    geom : UDGeom
        Geometry object with flat surface
        
    Examples
    --------
    >>> geom = create_flat_surface(100, 100, edgelength=10)
    >>> geom.save('flat_surface.stl')
    """
    if edgelength is None:
        edgelength = xsize / 10
    
    # Calculate number of divisions
    nx = int(np.ceil(xsize / edgelength))
    ny = int(np.ceil(ysize / edgelength))
    
    # Create grid of vertices
    x = np.linspace(0, xsize, nx + 1)
    y = np.linspace(0, ysize, ny + 1)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    
    # Flatten to vertex list
    vertices = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    
    # Create triangular faces (two triangles per grid cell)
    faces = []
    for j in range(ny):
        for i in range(nx):
            # Vertex indices for this cell
            v0 = j * (nx + 1) + i
            v1 = v0 + 1
            v2 = v0 + (nx + 1)
            v3 = v2 + 1
            
            # Two triangles per cell
            faces.append([v0, v1, v2])
            faces.append([v1, v3, v2])
    
    faces = np.array(faces)
    
    # Create trimesh object
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
    
    # Create UDGeom from mesh
    geom = UDGeom(stl=mesh)
    
    return geom


def create_cubes(xsize: float, ysize: float,
                Hx: float, Hy: float, Hz: float,
                Cx: float = 0, Cy: float = 0,
                geom_option: Literal['S', 'AC', 'SC'] = 'S',
                edgelength: Optional[float] = None,
                add_ground: bool = True) -> UDGeom:
    """
    Create cube geometries (single, aligned array, or staggered array).
    
    Parameters
    ----------
    xsize : float
        Domain length in x-direction [m]
    ysize : float
        Domain length in y-direction [m]
    Hx : float
        Cube length in x-direction [m]
    Hy : float
        Cube length in y-direction [m]
    Hz : float
        Cube height [m]
    Cx : float, default=0
        Spacing between cubes in x-direction [m] (for arrays)
    Cy : float, default=0
        Spacing between cubes in y-direction [m] (for arrays)
    geom_option : {'S', 'AC', 'SC'}, default='S'
        Type of geometry:
        - 'S': Single cube (centered in domain)
        - 'AC': Aligned cube array
        - 'SC': Staggered cube array
    edgelength : float, optional
        Ground facet edge length. If None, uses xsize / 20.
    add_ground : bool, default=True
        Whether to add ground surface
        
    Returns
    -------
    geom : UDGeom
        Geometry object with cubes
        
    Examples
    --------
    Single cube:
    >>> geom = create_cubes(100, 100, 10, 10, 15, geom_option='S')
    
    Aligned array:
    >>> geom = create_cubes(100, 100, 10, 10, 15, Cx=10, Cy=10, geom_option='AC')
    
    Staggered array:
    >>> geom = create_cubes(100, 100, 10, 10, 15, Cx=10, Cy=10, geom_option='SC')
    """
    if edgelength is None:
        edgelength = xsize / 20
    
    # Create cubes based on option
    if geom_option == 'S':
        # Single cube centered in domain
        center = [xsize / 2, ysize / 2, Hz / 2]
        cube = trimesh.creation.box(extents=[Hx, Hy, Hz])
        cube.apply_translation(center)
        meshes = [cube]
        
    elif geom_option == 'AC':
        # Aligned cube array
        Nx = int(xsize / (Hx + Cx))
        Ny = int(ysize / (Hy + Cy))
        
        if Nx * (Hx + Cx) > xsize or Ny * (Hy + Cy) > ysize:
            raise ValueError("Domain size must accommodate cube array")
        
        meshes = []
        for i in range(Nx):
            for j in range(Ny):
                center_x = (i + 0.5) * (Hx + Cx)
                center_y = (j + 0.5) * (Hy + Cy)
                center_z = Hz / 2
                
                cube = trimesh.creation.box(extents=[Hx, Hy, Hz])
                cube.apply_translation([center_x, center_y, center_z])
                meshes.append(cube)
                
    elif geom_option == 'SC':
        # Staggered cube array
        Nx = int(xsize / (Hx + Cx))
        Ny = int(ysize / (Hy + Cy))
        
        if Nx * (Hx + Cx) > xsize or Ny * (Hy + Cy) > ysize:
            raise ValueError("Domain size must accommodate cube array")
        
        meshes = []
        for i in range(Nx):
            for j in range(Ny):
                # Stagger every other row
                offset_y = (Hy + Cy) / 2 if i % 2 == 1 else 0
                
                center_x = (i + 0.5) * (Hx + Cx)
                center_y = (j + 0.5) * (Hy + Cy) + offset_y
                center_z = Hz / 2
                
                # Check if cube is within domain
                if center_y - Hy/2 >= 0 and center_y + Hy/2 <= ysize:
                    cube = trimesh.creation.box(extents=[Hx, Hy, Hz])
                    cube.apply_translation([center_x, center_y, center_z])
                    meshes.append(cube)
    else:
        raise ValueError(f"Invalid geom_option: {geom_option}. Must be 'S', 'AC', or 'SC'")
    
    # Combine all cube meshes
    if len(meshes) > 1:
        combined_mesh = trimesh.util.concatenate(meshes)
    else:
        combined_mesh = meshes[0]
    
    # Add ground if requested
    if add_ground:
        ground_geom = create_flat_surface(xsize, ysize, edgelength)
        # Combine with cubes
        final_mesh = trimesh.util.concatenate([combined_mesh, ground_geom.stl])
    else:
        final_mesh = combined_mesh
    
    # Create UDGeom
    geom = UDGeom(stl=final_mesh)
    
    return geom


def create_canyons(xsize: float, ysize: float,
                  Hx: float, Hy: float, Hz: float,
                  Cx: float, Cy: float,
                  orientation: Literal['x', 'y'] = 'y',
                  edgelength: Optional[float] = None,
                  add_ground: bool = True) -> UDGeom:
    """
    Create street canyon geometry (parallel building rows).
    
    Parameters
    ----------
    xsize : float
        Domain length in x-direction [m]
    ysize : float
        Domain length in y-direction [m]
    Hx : float
        Building length in x-direction [m]
    Hy : float
        Building length in y-direction [m]
    Hz : float
        Building height [m]
    Cx : float
        Street width in x-direction [m]
    Cy : float
        Street width in y-direction [m]
    orientation : {'x', 'y'}, default='y'
        Canyon orientation:
        - 'x': Canyons aligned along x-axis (buildings in y-direction)
        - 'y': Canyons aligned along y-axis (buildings in x-direction)
    edgelength : float, optional
        Ground facet edge length. If None, uses min(xsize, ysize) / 20.
    add_ground : bool, default=True
        Whether to add ground surface
        
    Returns
    -------
    geom : UDGeom
        Geometry object with canyon configuration
        
    Examples
    --------
    >>> geom = create_canyons(200, 100, 50, 10, 20, Cx=50, Cy=10, orientation='x')
    """
    if edgelength is None:
        edgelength = min(xsize, ysize) / 20
    
    meshes = []
    
    if orientation == 'y':
        # Canyons along y-axis, buildings repeated in x-direction
        N_buildings = int(xsize / (Hx + Cx))
        
        for i in range(N_buildings):
            center_x = (i + 0.5) * (Hx + Cx)
            center_y = ysize / 2
            center_z = Hz / 2
            
            # Create building extending full length in y
            building = trimesh.creation.box(extents=[Hx, ysize, Hz])
            building.apply_translation([center_x, center_y, center_z])
            meshes.append(building)
            
    elif orientation == 'x':
        # Canyons along x-axis, buildings repeated in y-direction
        N_buildings = int(ysize / (Hy + Cy))
        
        for j in range(N_buildings):
            center_x = xsize / 2
            center_y = (j + 0.5) * (Hy + Cy)
            center_z = Hz / 2
            
            # Create building extending full length in x
            building = trimesh.creation.box(extents=[xsize, Hy, Hz])
            building.apply_translation([center_x, center_y, center_z])
            meshes.append(building)
    else:
        raise ValueError(f"Invalid orientation: {orientation}. Must be 'x' or 'y'")
    
    # Combine building meshes
    if len(meshes) > 1:
        combined_mesh = trimesh.util.concatenate(meshes)
    else:
        combined_mesh = meshes[0]
    
    # Add ground if requested
    if add_ground:
        ground_geom = create_flat_surface(xsize, ysize, edgelength)
        final_mesh = trimesh.util.concatenate([combined_mesh, ground_geom.stl])
    else:
        final_mesh = combined_mesh
    
    # Create UDGeom
    geom = UDGeom(stl=final_mesh)
    
    return geom


def create_realistic(xsize: float, ysize: float,
                    building_configs: list,
                    edgelength: Optional[float] = None,
                    add_ground: bool = True) -> UDGeom:
    """
    Create realistic building configuration with varying heights and positions.
    
    Parameters
    ----------
    xsize : float
        Domain length in x-direction [m]
    ysize : float
        Domain length in y-direction [m]
    building_configs : list of dict
        List of building configurations. Each dict should contain:
        - 'position': (x, y) tuple for building center [m]
        - 'size': (length_x, length_y, height) tuple [m]
        Optional:
        - 'rotation': rotation angle in degrees (default: 0)
    edgelength : float, optional
        Ground facet edge length. If None, uses min(xsize, ysize) / 20.
    add_ground : bool, default=True
        Whether to add ground surface
        
    Returns
    -------
    geom : UDGeom
        Geometry object with realistic building layout
        
    Examples
    --------
    >>> buildings = [
    ...     {'position': (25, 25), 'size': (20, 20, 30)},
    ...     {'position': (75, 25), 'size': (15, 25, 25)},
    ...     {'position': (25, 75), 'size': (18, 18, 35), 'rotation': 15},
    ...     {'position': (75, 75), 'size': (22, 16, 28)},
    ... ]
    >>> geom = create_realistic(100, 100, buildings)
    """
    if edgelength is None:
        edgelength = min(xsize, ysize) / 20
    
    meshes = []
    
    for config in building_configs:
        pos = config['position']
        size = config['size']
        rotation = config.get('rotation', 0)
        
        Hx, Hy, Hz = size
        center_x, center_y = pos
        center_z = Hz / 2
        
        # Create building
        building = trimesh.creation.box(extents=[Hx, Hy, Hz])
        
        # Apply rotation if specified
        if rotation != 0:
            rotation_matrix = trimesh.transformations.rotation_matrix(
                np.radians(rotation), [0, 0, 1], [0, 0, center_z]
            )
            building.apply_transform(rotation_matrix)
        
        # Translate to position
        building.apply_translation([center_x, center_y, center_z])
        
        # Check bounds
        if (center_x - Hx/2 < 0 or center_x + Hx/2 > xsize or
            center_y - Hy/2 < 0 or center_y + Hy/2 > ysize):
            import warnings
            warnings.warn(f"Building at ({center_x}, {center_y}) extends beyond domain bounds")
        
        meshes.append(building)
    
    # Combine all buildings
    if len(meshes) > 1:
        combined_mesh = trimesh.util.concatenate(meshes)
    elif len(meshes) == 1:
        combined_mesh = meshes[0]
    else:
        raise ValueError("No buildings provided")
    
    # Add ground if requested
    if add_ground:
        ground_geom = create_flat_surface(xsize, ysize, edgelength)
        final_mesh = trimesh.util.concatenate([combined_mesh, ground_geom.stl])
    else:
        final_mesh = combined_mesh
    
    # Create UDGeom
    geom = UDGeom(stl=final_mesh)
    
    return geom


if __name__ == "__main__":
    print("Geometry generation module loaded successfully")
    print("\nAvailable functions:")
    print("  - create_flat_surface(xsize, ysize, edgelength)")
    print("  - create_cubes(xsize, ysize, Hx, Hy, Hz, Cx, Cy, geom_option)")
    print("  - create_canyons(xsize, ysize, Hx, Hy, Hz, Cx, Cy, orientation)")
    print("  - create_realistic(xsize, ysize, building_configs)")
