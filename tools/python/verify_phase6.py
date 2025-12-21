"""
Quick verification script for Phase 6 geometry generation module.
Tests basic imports and function signatures without requiring example data.
"""

import sys
from pathlib import Path

# Add tools/python to path
sys.path.insert(0, str(Path(__file__).parent))

print("=" * 60)
print("Phase 6 Quick Verification")
print("=" * 60 + "\n")

# Test 1: Import check
print("Test 1: Module import...")
try:
    from geometry_generation import (
        create_flat_surface,
        create_cubes,
        create_canyons,
        create_realistic
    )
    print("✓ Module imported successfully\n")
except ImportError as e:
    print(f"✗ Import failed: {e}\n")
    sys.exit(1)

# Test 2: Function signatures
print("Test 2: Function signatures...")
import inspect

functions = {
    'create_flat_surface': create_flat_surface,
    'create_cubes': create_cubes,
    'create_canyons': create_canyons,
    'create_realistic': create_realistic
}

for name, func in functions.items():
    sig = inspect.signature(func)
    params = list(sig.parameters.keys())
    print(f"✓ {name}({', '.join(params)})")

print("\n" + "=" * 60)
print("Phase 6 Module Structure Complete ✓")
print("=" * 60)
print("\nNote: Full tests require numpy, trimesh, and other dependencies.")
print("Install with: pip install -r requirements.txt")
