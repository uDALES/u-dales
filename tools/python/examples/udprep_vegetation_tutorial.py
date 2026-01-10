from __future__ import annotations

"""
User-facing walkthrough for vegetation preprocessing.

This tutorial keeps the short structure of the original script while adding
context for each step. It demonstrates:
- Loading the test vegetation case using UDPrep/UDBase.
- Loading a block file for vegetation placement.
- Loading an STL file for tree geometry.
"""

from pathlib import Path
import sys

udales_path = Path(__file__).resolve().parents[3]
PYTOOLS = udales_path / "tools" / "python"
if str(PYTOOLS) not in sys.path:
    sys.path.insert(0, str(PYTOOLS))

from udprep import UDPrep

# Test case included with uDALES
expnr = "525"
expdir = udales_path / "tests" / "tests_tree_input"

# Instantiate UDPrep (loads UDBase + geometry by default).
prep = UDPrep(expnr, expdir)

# Load vegetation block file (tree bounding boxes)
print("Load block file...")
fig = prep.vegetation.load_block('trees.inp.525')

prep.vegetation.save()  # Save veg.inp/veg_params.inp

fig.show(renderer="browser")

# Load STL tree geometry for visualization / preprocessing
print("Load STL file...")
stlfile = expdir / "veg_spheres_525.stl"
fig = prep.vegetation.load_stl(stlfile)
fig.show(renderer="browser")
