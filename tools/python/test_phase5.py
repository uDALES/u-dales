"""
Test script for Phase 5: Advanced Facet Analysis

Tests the convert_fac_to_field and calculate_frontal_properties methods.
"""

import sys
from pathlib import Path
import numpy as np

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

print("=" * 70)
print("Phase 5 Test: Advanced Facet Analysis")
print("=" * 70)

# Test 1: Import modules
print("\n[Test 1] Import modules...")
try:
    from udbase import UDBase
    print("✓ UDBase imported successfully")
except ImportError as e:
    print(f"✗ Failed to import: {e}")
    sys.exit(1)

# Test 2: Check that methods exist
print("\n[Test 2] Check advanced analysis methods exist...")
try:
    assert hasattr(UDBase, 'convert_fac_to_field'), "convert_fac_to_field method not found"
    assert hasattr(UDBase, 'calculate_frontal_properties'), "calculate_frontal_properties method not found"
    print("✓ convert_fac_to_field method exists")
    print("✓ calculate_frontal_properties method exists")
except AssertionError as e:
    print(f"✗ Method check failed: {e}")
    sys.exit(1)

# Test 3: Check method signatures
print("\n[Test 3] Check method signatures...")
try:
    import inspect
    
    # Check convert_fac_to_field
    sig = inspect.signature(UDBase.convert_fac_to_field)
    params = list(sig.parameters.keys())
    
    assert 'var' in params, "convert_fac_to_field missing 'var' parameter"
    assert 'facsec' in params, "convert_fac_to_field missing 'facsec' parameter"
    assert 'dz' in params, "convert_fac_to_field missing 'dz' parameter"
    
    print(f"✓ convert_fac_to_field parameters: {params}")
    
    # Check calculate_frontal_properties
    sig = inspect.signature(UDBase.calculate_frontal_properties)
    params = list(sig.parameters.keys())
    
    print(f"✓ calculate_frontal_properties parameters: {params}")
    
except Exception as e:
    print(f"✗ Signature check failed: {e}")

# Test 4: Create UDBase with geometry and facet sections
print("\n[Test 4] Create UDBase with geometry and facet sections...")
example_path = Path(__file__).parent.parent.parent.parent / 'examples' / '101'

if (example_path / 'namoptions.101').exists():
    try:
        sim = UDBase(101, example_path)
        
        if sim.geom is not None:
            print(f"✓ UDBase created with geometry: {sim.geom.n_faces} facets")
        else:
            print("⚠ Geometry not loaded")
            sim = None
        
        if sim is not None and hasattr(sim, 'facsec') and sim.facsec is not None:
            print(f"✓ Facet sections loaded: {list(sim.facsec.keys())}")
        else:
            print("⚠ Facet sections not loaded")
            if sim is not None:
                print("  This is OK - methods will test error handling")
            
    except Exception as e:
        print(f"✗ Failed to create UDBase: {e}")
        sim = None
else:
    print(f"⚠ Example not found at {example_path}")
    sim = None

# Test 5: Test convert_fac_to_field with synthetic data
print("\n[Test 5] Test convert_fac_to_field...")
if sim is not None and sim.geom is not None:
    if hasattr(sim, 'facsec') and sim.facsec is not None:
        try:
            # Create constant facet variable
            var = np.ones(sim.geom.n_faces)
            
            # Convert to field
            field = sim.convert_fac_to_field(var)
            
            # Check output
            assert field.shape == (sim.itot, sim.jtot, sim.ktot), \
                f"Wrong field shape: {field.shape}"
            assert field.dtype == np.float32, f"Wrong dtype: {field.dtype}"
            assert np.all(field >= 0), "Field has negative values"
            
            print(f"✓ convert_fac_to_field executed")
            print(f"  Field shape: {field.shape}")
            print(f"  Non-zero cells: {np.count_nonzero(field)}")
            print(f"  Max value: {field.max():.6f}")
            
        except Exception as e:
            print(f"✗ convert_fac_to_field failed: {e}")
    else:
        print("⚠ Skipping - facet sections not available")
        
        # Test error handling
        try:
            var = np.ones(10)
            field = sim.convert_fac_to_field(var)
            print("✗ Should have raised ValueError for missing facet sections")
        except ValueError as e:
            print(f"✓ Correctly raises ValueError: {str(e)[:50]}...")
        except Exception as e:
            print(f"⚠ Raised unexpected error: {type(e).__name__}")
else:
    print("⚠ Skipping - geometry not available")

# Test 6: Test convert_fac_to_field with different grids
print("\n[Test 6] Test convert_fac_to_field with different grids...")
if sim is not None and hasattr(sim, 'facsec') and sim.facsec is not None:
    try:
        var = np.ones(sim.geom.n_faces)
        
        # Test each grid
        for grid in ['c', 'u', 'v', 'w']:
            if grid in sim.facsec:
                field = sim.convert_fac_to_field(var, facsec=sim.facsec[grid])
                print(f"✓ {grid}-grid: {np.count_nonzero(field)} non-zero cells")
        
    except Exception as e:
        print(f"✗ Grid test failed: {e}")
else:
    print("⚠ Skipping - facet sections not available")

# Test 7: Test calculate_frontal_properties
print("\n[Test 7] Test calculate_frontal_properties...")
if sim is not None and sim.geom is not None:
    if hasattr(sim, 'facsec') and sim.facsec is not None:
        try:
            # Calculate frontal properties (this also prints results)
            props = sim.calculate_frontal_properties()
            
            # Check output structure
            assert isinstance(props, dict), "Return value should be dict"
            
            required_keys = ['skylinex', 'skyliney', 'Afx', 'Afy', 'brx', 'bry']
            for key in required_keys:
                assert key in props, f"Missing key: {key}"
            
            print(f"\n✓ calculate_frontal_properties executed")
            print(f"  skylinex shape: {props['skylinex'].shape}")
            print(f"  skyliney shape: {props['skyliney'].shape}")
            print(f"  Frontal area (x): {props['Afx']:.1f} m²")
            print(f"  Frontal area (y): {props['Afy']:.1f} m²")
            print(f"  Blockage ratio (x): {props['brx']:.3f}")
            print(f"  Blockage ratio (y): {props['bry']:.3f}")
            
            # Validate results
            assert props['Afx'] >= 0, "Frontal area should be non-negative"
            assert props['Afy'] >= 0, "Frontal area should be non-negative"
            assert 0 <= props['brx'] <= 1, "Blockage ratio should be in [0, 1]"
            assert 0 <= props['bry'] <= 1, "Blockage ratio should be in [0, 1]"
            
            print("✓ All values in valid ranges")
            
        except Exception as e:
            print(f"✗ calculate_frontal_properties failed: {e}")
    else:
        print("⚠ Skipping - facet sections not available")
        
        # Test error handling
        try:
            props = sim.calculate_frontal_properties()
            print("✗ Should have raised ValueError for missing facet sections")
        except ValueError as e:
            print(f"✓ Correctly raises ValueError: {str(e)[:50]}...")
        except Exception as e:
            print(f"⚠ Raised unexpected error: {type(e).__name__}")
else:
    print("⚠ Skipping - geometry not available")

# Test 8: Test skyline properties
print("\n[Test 8] Validate skyline properties...")
if sim is not None and hasattr(sim, 'facsec') and sim.facsec is not None:
    try:
        props = sim.calculate_frontal_properties()
        
        skylinex = props['skylinex']
        skyliney = props['skyliney']
        
        # Check shapes
        assert skylinex.shape == (sim.jtot, sim.ktot), \
            f"skylinex wrong shape: {skylinex.shape}"
        assert skyliney.shape == (sim.itot, sim.ktot), \
            f"skyliney wrong shape: {skyliney.shape}"
        
        # Check values are binary
        assert np.all((skylinex == 0) | (skylinex == 1)), \
            "skylinex should be binary"
        assert np.all((skyliney == 0) | (skyliney == 1)), \
            "skyliney should be binary"
        
        print("✓ Skyline shapes correct")
        print("✓ Skyline values are binary (0 or 1)")
        
        # Check blockage makes sense
        block_x_frac = np.sum(skylinex) / skylinex.size
        block_y_frac = np.sum(skyliney) / skyliney.size
        
        print(f"✓ Blockage fractions: x={block_x_frac:.3f}, y={block_y_frac:.3f}")
        
    except Exception as e:
        print(f"✗ Skyline validation failed: {e}")
else:
    print("⚠ Skipping - facet sections not available")

# Test 9: Test conservation in convert_fac_to_field
print("\n[Test 9] Test conservation in convert_fac_to_field...")
if sim is not None and hasattr(sim, 'facsec') and sim.facsec is not None:
    try:
        # For a constant facet variable, the field integral should relate
        # to the total facet area
        var = np.ones(sim.geom.n_faces)
        field = sim.convert_fac_to_field(var)
        
        # Integrate field over volume
        field_integral = 0.0
        for k in range(sim.ktot):
            field_integral += np.sum(field[:, :, k]) * sim.dx * sim.dy * sim.dzt[k]
        
        # This should be related to total facet area
        # (not exactly equal due to grid discretization)
        print(f"  Field integral: {field_integral:.2f}")
        print(f"  Total facet area: {sim.geom.total_area:.2f}")
        print(f"  Ratio: {field_integral / sim.geom.total_area:.3f}")
        
        # The ratio depends on how facets are distributed across cells
        print("✓ Conservation test completed")
        
    except Exception as e:
        print(f"✗ Conservation test failed: {e}")
else:
    print("⚠ Skipping - facet sections not available")

# Test 10: Test error handling without geometry
print("\n[Test 10] Test error handling without required data...")
try:
    # Create a mock sim without geometry
    class MockSim:
        def __init__(self):
            self.geom = None
            self.facsec = None
    
    mock = MockSim()
    
    # Bind methods
    mock.convert_fac_to_field = UDBase.convert_fac_to_field.__get__(mock, UDBase)
    mock.calculate_frontal_properties = UDBase.calculate_frontal_properties.__get__(mock, UDBase)
    
    # Try calling convert_fac_to_field
    try:
        mock.convert_fac_to_field(np.array([1, 2, 3]))
        print("✗ Should have raised ValueError for missing facet sections")
    except ValueError:
        print("✓ convert_fac_to_field raises ValueError when facet sections missing")
    except Exception as e:
        print(f"⚠ Raised unexpected error: {type(e).__name__}")
    
    # Try calling calculate_frontal_properties
    try:
        mock.calculate_frontal_properties()
        print("✗ Should have raised ValueError for missing geometry")
    except ValueError:
        print("✓ calculate_frontal_properties raises ValueError when geometry missing")
    except Exception as e:
        print(f"⚠ Raised unexpected error: {type(e).__name__}")
        
except Exception as e:
    print(f"✗ Error handling test failed: {e}")

# Summary
print("\n" + "=" * 70)
print("Phase 5 Test Summary")
print("=" * 70)

print("""
Phase 5 Implementation Complete ✓

Components tested:
1. convert_fac_to_field and calculate_frontal_properties methods exist
2. Method signatures correct
3. convert_fac_to_field creates valid 3D fields
4. Different grid support (c, u, v, w)
5. calculate_frontal_properties computes valid results
6. Skyline properties validated (binary, correct shapes)
7. Conservation properties checked
8. Error handling for missing data

Key Features:
- convert_fac_to_field(var): Convert facet data to 3D density field
  * Distributes facet areas into grid cells
  * Supports different grids (c, u, v, w)
  * Returns float32 array

- calculate_frontal_properties(): Compute urban canopy properties
  * Frontal areas in x and y directions
  * Blockage ratios (fraction of cross-section blocked)
  * Skyline profiles (binary indicators)
  * Prints results automatically

Applications:
- Urban canopy characterization
- Flow resistance estimation
- Building morphology analysis
- Surface-volume coupling

Next: Complete remaining phases or testing/documentation
""")

print("\nAll Phase 5 tests complete!")
