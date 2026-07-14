from __future__ import annotations

"""
User-facing walkthrough for vegetation preprocessing.

This tutorial keeps the short structure of the original script while adding
context for each step. It demonstrates:
- Loading the test vegetation case using UDPrep/UDBase.
- Loading a block file (tree bounding boxes) and writing the sparse
  veg.inp/veg_params.inp files the solver consumes.
- Loading an STL tree volume and voxelising it onto the grid.
- Visualising the resulting vegetation points on top of the geometry.

Notes
-----
- ``load_block``/``load_stl`` populate ``prep.sim.veg`` (a dict of sparse
  points + per-point parameters) and update ``ntrees``; they do not return a
  figure. ``load_stl`` additionally returns a ready-to-show scene of the result.
- Visualisation goes through ``prep.sim.plot_veg`` and uses the package default
  rendering backend (PyVista). Set ``SHOW = False`` (or ``UDALES_SHOW=0``) to run
  the preprocessing without opening any interactive window.
"""

from pathlib import Path
import os
import sys

udales_path = Path(__file__).resolve().parents[3]
PYTOOLS = udales_path / "tools" / "python"
if str(PYTOOLS) not in sys.path:
    sys.path.insert(0, str(PYTOOLS))

from udprep import UDPrep  # noqa: E402

# Open interactive 3-D windows for the plots. Defaults to True for users; set the
# environment variable UDALES_SHOW=0 to run headless (e.g. in CI/verification).
SHOW = os.environ.get("UDALES_SHOW", "1") != "0"

# Test case included with uDALES
expnr = "525"
expdir = udales_path / "tests" / "cases" / expnr

# Instantiate UDPrep (loads UDBase + geometry by default).
prep = UDPrep(expnr, expdir)

# -----------------------------------------------------------------------------
# Block file: tree bounding boxes -> sparse vegetation points
# -----------------------------------------------------------------------------
print("Load block file...")
prep.vegetation.load_block("trees.inp.525")  # populates prep.sim.veg
prep.vegetation.save()  # writes veg.inp.<expnr> / veg_params.inp.<expnr>

# Visualise the vegetation points from the block file on top of the geometry.
if SHOW:
    prep.sim.plot_veg(show=True)

# -----------------------------------------------------------------------------
# STL file: voxelise a closed tree volume onto the grid
# -----------------------------------------------------------------------------
# This overwrites prep.sim.veg with the points found inside the STL volume and
# returns a scene (a PyVista Plotter or Plotly Figure, depending on the backend).
print("Load STL file...")
stlfile = expdir / "veg_spheres_525.stl"
fig = prep.vegetation.load_stl(stlfile)
if SHOW and fig is not None:
    fig.show()
