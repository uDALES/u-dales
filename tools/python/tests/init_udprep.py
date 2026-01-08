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

# Instantiate UDPrep either as ...
sim = UDBase(expnr, expdir)
prep = UDPrep(sim)

# ... or directly from experiment directory:
prep = UDPrep(expnr, expdir)

# Summary of the derived preprocessing configuration.
print(prep)

# Example access patterns.
print(f"ibm dx: {prep.ibm.dx}")
print(f"ic thl0: {prep.ic.thl0}")
print(f"vegetation ltrees: {prep.vegetation.ltrees}")
print(f"scalars nsv: {prep.scalars.nsv}")
print(f"seb facT: {prep.seb.facT}")

# Example 1: override a parameter, then run a single IC step.
# (Here we adjust the initial potential temperature before writing profiles.)
prep.ic.thl0 = 290.0
prep.ic.generate_prof()
prep.ic.write_prof()

# Example 2: adjust vegetation switch and run vegetation-only preprocessing.
prep.vegetation.ltrees = 1
prep.vegetation.run_all()

# Example 3: run a single section's run_all explicitly.
prep.ibm.run_all()

# Example 4: record parameter changes for a section.
# This will eventually diff against defaults and write only changed params.
prep.ic.write_changed_params()

# Example 3: run the full preprocessing chain.
prep.run_all()