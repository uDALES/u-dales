from __future__ import annotations

"""
User-facing walkthrough for UDPrep.

This tutorial demonstrates the current workflow:
- Load a case with UDPrep.
- Inspect section values.
- Override a few section parameters.
- Run section-level preprocessing.
- Show and write changed parameters.
- Run the full preprocessing chain.
"""

import sys
from pathlib import Path

# Resolve uDALES python tools from the script location.
udales_path = Path(__file__).resolve().parents[3]
pytools = udales_path / "tools" / "python"
if str(pytools) not in sys.path:
    sys.path.insert(0, str(pytools))

from udprep import UDPrep  # noqa: E402

# Tutorial case directory.
expnr = "526"
expdir = udales_path.parents[0] / "experiments" / expnr

print("Initializing UDPrep tutorial...")
prep = UDPrep(expnr, expdir)

print("-------------------------------------------------------------------")
print("Summary of the derived preprocessing configuration")
print(prep)

print("-------------------------------------------------------------------")
print("... or for an individual section")
print(prep.ibm)

print("-------------------------------------------------------------------")
print("Accessing section parameters...")
print(f"ibm dx: {prep.ibm.dx}")
print(f"ic thl0: {prep.ic.thl0}")
print(f"vegetation ltrees: {prep.vegetation.ltrees}")
print(f"scalars nsv: {prep.scalars.nsv}")
print(f"seb facT: {prep.seb.facT}")

print("-------------------------------------------------------------------")
print("Example 1: override IC parameter and write profiles...")
# Update an IC parameter and regenerate profile inputs.
prep.ic.thl0 = 290.0
prep.ic.generate_prof()
prep.ic.write_prof()

print("-------------------------------------------------------------------")
print("Example 2: enable vegetation and run vegetation preprocessing...")
prep.vegetation.ltrees = 1
prep.vegetation.run_all()
prep.vegetation.save()

print("-------------------------------------------------------------------")
print("Example 3: run IBM section preprocessing...")
prep.ibm.run_all()

print("-------------------------------------------------------------------")
print("Example 4: write changed params for IC section...")
# Write changed parameters from one section back to namelist files.
prep.ic.write_changed_params()

print("-------------------------------------------------------------------")
print("Example 4b: show changed params for all sections...")
prep.show_changed_params()

print("-------------------------------------------------------------------")
print("Example 4c: write changed params for all sections...")
prep.write_changed_params()

print("-------------------------------------------------------------------")
print("Example 5: run full preprocessing chain (recommended one-shot call)...")
# run_all applies current UDPrep gating, so optional sections execute only when enabled.
prep.run_all(force=True)
