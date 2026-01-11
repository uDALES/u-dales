from __future__ import annotations

import time
from pathlib import Path

import numpy as np

start = time.perf_counter()

import sys

# Add the uDALES python path
udbase_path = Path("C:/Users/mvr/OneDrive - Imperial College London/codes/uDALES/u-dales").resolve()
tools_path = (udbase_path / "tools" / "python").resolve()
if tools_path not in sys.path:
    sys.path.insert(0, str(tools_path))

expnr = "065"
expdir = (udbase_path.parents[0] / "experiments" / expnr).resolve()

from udprep import UDPrep  # noqa: E402

elapsed = time.perf_counter() - start
print(f"loading libraries runtime: {elapsed:.3f} s")

start = time.perf_counter()
prep = UDPrep(expnr, expdir)
sim = prep.sim
elapsed = time.perf_counter() - start
print(f"UDbase startup runtime: {elapsed:.3f} s")

maxD = 250.0

vf, svf, paths = prep.radiation.calc_view_factors(maxD=maxD)

print(f"VS3 path: {paths['vs3']}")
print(f"VF path: {paths['vf']}")
print(f"SVF path: {paths['svf']}")
if paths["vfsparse"] is not None:
    print(f"VF sparse path: {paths['vfsparse']}")
print(f"Nonzero VF entries: {vf.nnz}")
