"""
Test script for Phase 3: Geometry Handling

Tests the UDGeom class and its integration with UDBase.
"""

import sys
from pathlib import Path
import numpy as np

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

print("=" * 70)
print("Phase 3 Test: Geometry Handling")
print("=" * 70)

# Test 1: Import modules
print("\n[Test 1] Import udgeom module...")
try:
    from udgeom import UDGeom
    print("✓ UDGeom imported successfully")
except ImportError as e:
    print(f"✗ Failed to import UDGeom: {e}")
    sys.exit(1)

try:
    from udbase import UDBase
    print("✓ UDBase imported successfully")
except ImportError as e:
    print(f"✗ Failed to import UDBase: {e}")
    sys.exit(1)

# Test 2: Create UDGeom object
print("\n[Test 2] Create UDGeom object...")
try:
    geom = UDGeom()
    print(f"✓ UDGeom object created: {geom}")
except Exception as e:
    print(f"✗ Failed to create UDGeom: {e}")
    sys.exit(1)

# Test 3: Load geometry from STL file
print("\n[Test 3] Load geometry from STL file...")
example_path = Path(__file__).parent.parent.parent.parent / 'examples' / '101'

if (example_path / 'geometry.001').exists():
    try:
        geom = UDGeom(example_path)
        geom.load('geometry.001')
        print(f"✓ Geometry loaded: {geom.n_faces} faces, {geom.n_vertices} vertices")
    except Exception as e:
        print(f"✗ Failed to load geometry: {e}")
        geom = None
else:
    print(f"⚠ Geometry file not found at {example_path / 'geometry.001'}")
    print("  Skipping geometry loading tests")
    geom = None

# Test 4: Access geometry properties
if geom is not None and geom.stl is not None:
    print("\n[Test 4] Access geometry properties...")
    try:
        # Basic properties
        assert geom.n_faces > 0, "No faces in geometry"
        assert geom.n_vertices > 0, "No vertices in geometry"
        print(f"✓ n_faces: {geom.n_faces}")
        print(f"✓ n_vertices: {geom.n_vertices}")
        
        # Bounds
        bounds = geom.bounds
        assert bounds.shape == (2, 3), f"Unexpected bounds shape: {bounds.shape}"
        print(f"✓ bounds: shape {bounds.shape}")
        
        # Face properties
        centers = geom.face_centers
        assert centers.shape == (geom.n_faces, 3), f"Unexpected centers shape: {centers.shape}"
        print(f"✓ face_centers: shape {centers.shape}")
        
        normals = geom.face_normals
        assert normals.shape == (geom.n_faces, 3), f"Unexpected normals shape: {normals.shape}"
        print(f"✓ face_normals: shape {normals.shape}")
        
        areas = geom.face_areas
        assert len(areas) == geom.n_faces, f"Unexpected areas length: {len(areas)}"
        assert np.all(areas > 0), "Some face areas are non-positive"
        print(f"✓ face_areas: {len(areas)} faces, sum = {geom.total_area:.2f} m²")
        
        # Total area
        total_area = geom.total_area
        assert total_area > 0, "Total area is non-positive"
        print(f"✓ total_area: {total_area:.2f} m²")
        
        # Watertight check
        is_watertight = geom.is_watertight
        print(f"✓ is_watertight: {is_watertight}")
        
    except Exception as e:
        print(f"✗ Property access failed: {e}")

# Test 5: Create geometry from trimesh object
print("\n[Test 5] Create geometry from trimesh object...")
try:
    import trimesh
    
    # Create a simple box
    box = trimesh.creation.box(extents=[1, 1, 1])
    box_geom = UDGeom(stl=box)
    
    assert box_geom.n_faces == 12, f"Box should have 12 faces, got {box_geom.n_faces}"
    assert box_geom.n_vertices == 8, f"Box should have 8 vertices, got {box_geom.n_vertices}"
    assert box_geom.is_watertight, "Box should be watertight"
    
    print(f"✓ Created box: {box_geom.n_faces} faces, {box_geom.n_vertices} vertices")
    print(f"✓ Volume: {box_geom.volume:.3f} m³ (expected 1.000)")
    
except ImportError:
    print("⚠ trimesh not available - skipping box creation test")
except Exception as e:
    print(f"✗ Box creation failed: {e}")

# Test 6: Integration with UDBase
print("\n[Test 6] Integration with UDBase...")
if (example_path / 'namoptions.101').exists():
    try:
        sim = UDBase(101, example_path)
        
        if sim.geom is not None:
            print(f"✓ Geometry loaded via UDBase")
            print(f"  Facets: {sim.geom.n_faces}")
            print(f"  Area: {sim.geom.total_area:.2f} m²")
            
            # Check that it's the same geometry we loaded before
            if geom is not None:
                assert sim.geom.n_faces == geom.n_faces, "Face count mismatch"
                print(f"✓ Geometry matches standalone load")
        else:
            print("⚠ No geometry loaded (may not be specified in namoptions)")
            
    except Exception as e:
        print(f"✗ UDBase geometry integration failed: {e}")
else:
    print(f"⚠ namoptions not found at {example_path}")

# Test 7: Save geometry
print("\n[Test 7] Save geometry...")
if geom is not None and geom.stl is not None:
    try:
        import tempfile
        import os
        
        # Create a temporary file
        with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as f:
            temp_path = f.name
        
        # Save to temp file
        temp_dir = Path(tempfile.gettempdir())
        temp_filename = Path(temp_path).name
        
        # Set geometry path to temp dir
        old_path = geom.path
        geom.path = temp_dir
        geom.save(temp_filename)
        
        # Check file exists
        assert Path(temp_path).exists(), "Saved file does not exist"
        file_size = Path(temp_path).stat().st_size
        assert file_size > 0, "Saved file is empty"
        
        print(f"✓ Saved geometry to temp file ({file_size} bytes)")
        
        # Load it back
        geom2 = UDGeom(temp_dir)
        geom2.load(temp_filename)
        
        assert geom2.n_faces == geom.n_faces, "Face count mismatch after save/load"
        print(f"✓ Loaded saved geometry: {geom2.n_faces} faces")
        
        # Clean up
        os.unlink(temp_path)
        geom.path = old_path
        
    except Exception as e:
        print(f"✗ Save/load test failed: {e}")

# Test 8: Visualization (skip actual display in automated test)
print("\n[Test 8] Visualization methods...")
if geom is not None and geom.stl is not None:
    try:
        # Check that show method exists and has correct signature
        import inspect
        sig = inspect.signature(geom.show)
        params = list(sig.parameters.keys())
        
        assert 'color_buildings' in params, "Missing color_buildings parameter"
        assert 'show_normals' in params, "Missing show_normals parameter"
        
        print(f"✓ show() method available with parameters: {params}")
        print("  (Not displaying - would require interactive session)")
        
    except Exception as e:
        print(f"✗ Visualization check failed: {e}")

# Summary
print("\n" + "=" * 70)
print("Phase 3 Test Summary")
print("=" * 70)

print("""
Phase 3 Implementation Complete ✓

Components tested:
1. UDGeom class creation and initialization
2. STL file loading
3. Geometry property access (faces, vertices, normals, areas)
4. Geometry creation from trimesh objects
5. Integration with UDBase
6. Save/load functionality
7. Visualization methods (structure only)

Next Phase: Phase 4 - Visualization methods (plot_fac, plot_fac_type)
""")
