"""
Geometry Tutorial - uDALES Python Tools

This example demonstrates how to work with geometry in uDALES using the
UDGeom class and its integration with UDBase.

Based on the MATLAB tutorial: udales-geometry-tutorial.md
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path to import udbase and udgeom
sys.path.insert(0, str(Path(__file__).parent.parent))

from udbase import UDBase
from udgeom import UDGeom

# ============================================================================
# Example 1: Load geometry standalone
# ============================================================================
print("=" * 70)
print("Example 1: Load and inspect geometry")
print("=" * 70)

# Load a geometry file
geom = UDGeom('../../../examples/101')
geom.load('geometry.001')

# Print geometry information
print(f"\nGeometry loaded:")
print(f"  Number of facets: {geom.n_faces}")
print(f"  Number of vertices: {geom.n_vertices}")
print(f"  Total surface area: {geom.total_area:.2f} m²")
print(f"  Watertight: {geom.is_watertight}")

# Get bounding box
bounds = geom.bounds
print(f"\nBounding box:")
print(f"  X: [{bounds[0, 0]:.2f}, {bounds[1, 0]:.2f}] m")
print(f"  Y: [{bounds[0, 1]:.2f}, {bounds[1, 1]:.2f}] m")
print(f"  Z: [{bounds[0, 2]:.2f}, {bounds[1, 2]:.2f}] m")

# ============================================================================
# Example 2: Access geometry properties
# ============================================================================
print("\n" + "=" * 70)
print("Example 2: Access geometry properties")
print("=" * 70)

# Get face centers
centers = geom.face_centers
print(f"\nFace centers shape: {centers.shape}")
print(f"First 3 face centers:")
print(centers[:3])

# Get face normals
normals = geom.face_normals
print(f"\nFace normals shape: {normals.shape}")
print(f"First 3 face normals:")
print(normals[:3])

# Get face areas
areas = geom.face_areas
print(f"\nFace areas:")
print(f"  Min area: {areas.min():.6f} m²")
print(f"  Max area: {areas.max():.6f} m²")
print(f"  Mean area: {areas.mean():.6f} m²")

# ============================================================================
# Example 3: Analyze surface orientation
# ============================================================================
print("\n" + "=" * 70)
print("Example 3: Analyze surface orientation")
print("=" * 70)

# Classify surfaces by orientation
# Horizontal: |nz| > 0.9
# Vertical: |nz| < 0.1
nz = normals[:, 2]

n_horizontal = np.sum(np.abs(nz) > 0.9)
n_vertical = np.sum(np.abs(nz) < 0.1)
n_sloped = geom.n_faces - n_horizontal - n_vertical

print(f"\nSurface classification:")
print(f"  Horizontal facets: {n_horizontal} ({100*n_horizontal/geom.n_faces:.1f}%)")
print(f"  Vertical facets: {n_vertical} ({100*n_vertical/geom.n_faces:.1f}%)")
print(f"  Sloped facets: {n_sloped} ({100*n_sloped/geom.n_faces:.1f}%)")

# Ground vs building facets
n_ground = np.sum(centers[:, 2] <= 0)
n_building = np.sum(centers[:, 2] > 0)

print(f"\nSurface location:")
print(f"  Ground facets (z ≤ 0): {n_ground} ({100*n_ground/geom.n_faces:.1f}%)")
print(f"  Building facets (z > 0): {n_building} ({100*n_building/geom.n_faces:.1f}%)")

# ============================================================================
# Example 4: Visualize geometry (basic)
# ============================================================================
print("\n" + "=" * 70)
print("Example 4: Visualize geometry")
print("=" * 70)

print("\nDisplaying geometry with colored buildings and normals...")
print("(Close the plot window to continue)")

try:
    # Show with default settings: colored buildings, normals displayed
    geom.show()
except Exception as e:
    print(f"Visualization failed: {e}")
    print("This may happen if running in a headless environment.")

# ============================================================================
# Example 5: Visualize geometry (custom options)
# ============================================================================
print("\n" + "=" * 70)
print("Example 5: Visualize geometry (custom options)")
print("=" * 70)

print("\nDisplaying geometry without building colors (faster for large meshes)...")
print("(Close the plot window to continue)")

try:
    # Show without building colors, no normals, larger figure
    geom.show(color_buildings=False, show_normals=False, figsize=(12, 10))
except Exception as e:
    print(f"Visualization failed: {e}")

# ============================================================================
# Example 6: Load geometry via UDBase
# ============================================================================
print("\n" + "=" * 70)
print("Example 6: Load geometry via UDBase")
print("=" * 70)

# Create UDBase object (automatically loads geometry if stl_file in namoptions)
sim = UDBase(101, '../../../examples/101')

if sim.geom is not None:
    print(f"\nGeometry loaded via UDBase:")
    print(f"  Facets: {sim.geom.n_faces}")
    print(f"  Vertices: {sim.geom.n_vertices}")
    print(f"  Area: {sim.geom.total_area:.2f} m²")
    
    # Access geometry through simulation object
    print(f"\nGeometry bounds:")
    print(f"  X: [{sim.geom.bounds[0, 0]:.2f}, {sim.geom.bounds[1, 0]:.2f}] m")
    print(f"  Y: [{sim.geom.bounds[0, 1]:.2f}, {sim.geom.bounds[1, 1]:.2f}] m")
    print(f"  Z: [{sim.geom.bounds[0, 2]:.2f}, {sim.geom.bounds[1, 2]:.2f}] m")
else:
    print("\nNo geometry loaded (stl_file may not be in namoptions)")

# ============================================================================
# Example 7: Save modified geometry
# ============================================================================
print("\n" + "=" * 70)
print("Example 7: Save geometry")
print("=" * 70)

# You can save geometry to a new file
# (This would save to the geometry path)
# geom.save('modified_geometry.stl')

print("\nTo save geometry:")
print("  geom.save('output_filename.stl')")
print("  Saves to the geometry path directory")

# ============================================================================
# Example 8: Create geometry from trimesh
# ============================================================================
print("\n" + "=" * 70)
print("Example 8: Create geometry from trimesh object")
print("=" * 70)

try:
    import trimesh
    
    # Create a simple box geometry
    box = trimesh.creation.box(extents=[10, 10, 5])
    
    # Create UDGeom from trimesh object
    box_geom = UDGeom(stl=box)
    
    print(f"\nCreated box geometry:")
    print(f"  Facets: {box_geom.n_faces}")
    print(f"  Vertices: {box_geom.n_vertices}")
    print(f"  Volume: {box_geom.volume:.2f} m³")
    print(f"  Watertight: {box_geom.is_watertight}")
    
    # Visualize the box
    print("\nDisplaying box geometry...")
    print("(Close the plot window to continue)")
    box_geom.show()
    
except ImportError:
    print("\ntrimesh not available - skipping box creation example")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY: Geometry Operations")
print("=" * 70)

print("""
Key functionality demonstrated:

1. Load STL files:
   geom = UDGeom('path/to/sim')
   geom.load('geometry.001')

2. Access properties:
   - geom.n_faces, geom.n_vertices
   - geom.bounds, geom.total_area, geom.volume
   - geom.face_centers, geom.face_normals, geom.face_areas
   - geom.is_watertight

3. Visualize geometry:
   geom.show()  # Default: colored buildings + normals
   geom.show(color_buildings=False)  # Faster for large meshes
   geom.show(show_normals=False, figsize=(12, 10))  # Custom

4. Save geometry:
   geom.save('output.stl')

5. Integration with UDBase:
   sim = UDBase(101, 'path')
   if sim.geom is not None:
       sim.geom.show()

6. Create from trimesh:
   box = trimesh.creation.box(...)
   geom = UDGeom(stl=box)

The UDGeom class provides a simple interface to trimesh functionality
while maintaining compatibility with uDALES workflows.
""")

print("\nGeometry tutorial complete!")
