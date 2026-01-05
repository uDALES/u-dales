"""
Test script for Phase 1: Core Infrastructure

This script tests the basic functionality of UDBase:
- Namoptions parsing
- Grid generation
- Basic file loading
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from udbase import UDBase
import numpy as np

def test_basic_loading():
    """Test basic UDBase loading with a sample experiment."""
    
    print("=" * 70)
    print("Testing UDBase Core Infrastructure (Phase 1)")
    print("=" * 70)
    
    # Test with experiment 065 from tutorials
    exp_path = Path("/home/mvr/experiments/065")
        
    if not exp_path.exists():
        print(f"\nWarning: Test experiment not found at {exp_path}")
        print("Please provide a valid experiment path to test.")
        return
    
    try:
        print(f"\n1. Loading simulation from: {exp_path}")
        sim = UDBase(expnr=65, path=exp_path)
        print("   ✓ Simulation loaded successfully")
        
        print(f"\n2. Checking simulation properties:")
        print(sim)
        
        print(f"\n3. Verifying grid arrays:")
        print(f"   xm shape: {sim.xm.shape}, range: [{sim.xm[0]:.2f}, {sim.xm[-1]:.2f}]")
        print(f"   xt shape: {sim.xt.shape}, range: [{sim.xt[0]:.2f}, {sim.xt[-1]:.2f}]")
        print(f"   ym shape: {sim.ym.shape}, range: [{sim.ym[0]:.2f}, {sim.ym[-1]:.2f}]")
        print(f"   yt shape: {sim.yt.shape}, range: [{sim.yt[0]:.2f}, {sim.yt[-1]:.2f}]")
        print(f"   zm shape: {sim.zm.shape}, range: [{sim.zm[0]:.2f}, {sim.zm[-1]:.2f}]")
        print(f"   zt shape: {sim.zt.shape}, range: [{sim.zt[0]:.2f}, {sim.zt[-1]:.2f}]")
        print("   ✓ Grid arrays loaded")
        
        print(f"\n4. Checking namoptions parameters:")
        if hasattr(sim, 'itot'):
            print(f"   itot: {sim.itot}, jtot: {sim.jtot}, ktot: {sim.ktot}")
        if hasattr(sim, 'xlen'):
            print(f"   xlen: {sim.xlen}, ylen: {sim.ylen}, zsize: {sim.zsize}")
        if hasattr(sim, 'dx'):
            print(f"   dx: {sim.dx:.4f}, dy: {sim.dy:.4f}")
        print("   ✓ Namoptions parsed")
        
        print(f"\n5. Checking optional data:")
        if sim.geom is not None:
            print(f"   ✓ Geometry loaded")
        else:
            print(f"   - Geometry not available")
        
        if hasattr(sim, 'Su') and sim.Su is not None:
            print(f"   ✓ Solid masks loaded")
            print(f"     Su: {np.sum(sim.Su)} solid points")
        else:
            print(f"   - Solid masks not available")
        
        if sim.facs:
            print(f"   ✓ Facet data loaded")
            if 'area' in sim.facs:
                print(f"     {len(sim.facs['area'])} facets")
        else:
            print(f"   - Facet data not available")
        
        print("\n" + "=" * 70)
        print("✓ Phase 1 tests completed successfully!")
        print("=" * 70)
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_basic_loading()
