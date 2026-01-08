from __future__ import annotations

import sys
from pathlib import Path

# Add the uDALES python path
udbase_path = Path("C:/Users/mvr/OneDrive - Imperial College London/codes/uDALES/u-dales").resolve()
tools_path = (udbase_path / "tools" / "python").resolve()
if str(tools_path) not in sys.path:
    sys.path.insert(0, str(tools_path))

expnr = "526"
expdir = (udbase_path.parents[0] / "experiments" / expnr).resolve()

from udbase import UDBase  # noqa: E402
from udprep import UDPrep  # noqa: E402

sim = UDBase(expnr, expdir)
prep = UDPrep(sim)

print("UDPrep initialized")
print(f"expnr: {prep.expnr}")
print(f"path: {prep.path}")
print(prep)

