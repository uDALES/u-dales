#!/usr/bin/env python3
"""
Simple test script to verify cross-platform configuration parsing.
Tests udprep_in.py without requiring full uDALES dependencies.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add tools/python to path
script_dir = Path(__file__).resolve().parent
tools_python = script_dir.parent
if str(tools_python) not in sys.path:
    sys.path.insert(0, str(tools_python))

# Test 1: Import the configuration module
print("=" * 70)
print("Test 1: Importing udprep_in module...")
try:
    from udprep.udprep_init import (
        parse_shell_config,
        read_iexpnr_from_namoptions,
        setup_paths_from_config,
    )
    print("✓ Successfully imported udprep_in functions")
except ImportError as e:
    print(f"✗ Import failed: {e}")
    sys.exit(1)

# Test 2: Test default path handling
print("\n" + "=" * 70)
print("Test 2: Testing default path configuration...")
if len(sys.argv) > 1:
    print(f"Config directory provided: {sys.argv[1]}")
    try:
        config_dir = Path(sys.argv[1]).resolve()
        expdir, expnr = setup_paths_from_config(config_dir)
        print(f"✓ Successfully parsed configuration:")
        print(f"  - expdir: {expdir}")
        print(f"  - expnr: {expnr}")
        
        # Verify paths are Path objects
        assert isinstance(expdir, Path), "expdir should be a Path object"
        assert isinstance(expnr, str), "expnr should be a string"
        print("✓ All path types are correct")
        
    except Exception as e:
        print(f"✗ Configuration parsing failed: {e}")
        sys.exit(1)
else:
    print("No config directory provided - using hardcoded Windows path as default")
    expnr = "526"
    # Derive paths from script location
    script_dir = Path(__file__).resolve().parent
    tools_python = script_dir.parent
    udales_root = tools_python.parent
    expdir = (udales_root.parent / "experiments" / expnr).resolve()
    
    print(f"  - expdir: {expdir}")
    print(f"  - expnr: {expnr}")
    
    # On WSL, this path won't exist but pathlib should still handle it
    if not expdir.exists():
        print(f"  ⚠ Path doesn't exist (expected on WSL without config)")

# Test 3: Platform detection
print("\n" + "=" * 70)
print("Test 3: Platform information...")
import platform
print(f"  - Platform: {platform.system()}")
print(f"  - Platform release: {platform.release()}")
print(f"  - Python version: {platform.python_version()}")
print(f"  - Path separator: '{Path('/').as_posix()}' vs '{Path('\\\\').as_posix()}'")

# Test 4: Path handling
print("\n" + "=" * 70)
print("Test 4: Cross-platform path handling...")
test_paths = [
    "C:/Users/test/file.txt",
    "/mnt/c/Users/test/file.txt",
    "relative/path/file.txt",
]
for path_str in test_paths:
    p = Path(path_str)
    print(f"  - Input: '{path_str}'")
    print(f"    → Path object: {p}")
    print(f"    → POSIX: {p.as_posix()}")

print("\n" + "=" * 70)
print("✓ All basic tests passed!")
print("=" * 70)
