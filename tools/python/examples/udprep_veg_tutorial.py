from __future__ import annotations

from pathlib import Path
import sys

udales_path = Path(__file__).resolve().parents[3]
PYTOOLS = udales_path / "tools" / "python"
if str(PYTOOLS) not in sys.path:
    sys.path.insert(0, str(PYTOOLS))

from udprep import UDPrep

expnr = "525"
expdir = udales_path / "tests" / "tests_tree_input"
stlfile = expdir / "veg_spheres_525.stl"

prep = UDPrep(expnr, expdir)

print("Load block file...")
fig = prep.vegetation.load_block()
fig.show(renderer="browser")

print("Load STL file...")
fig = prep.vegetation.load_stl(stlfile)
fig.show(renderer="browser")
