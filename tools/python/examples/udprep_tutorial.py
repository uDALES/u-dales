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

from udprep import UDPrep  # noqa: E402

print("Initializing UDPrep demo...")
prep = UDPrep(expnr, expdir)

print("-------------------------------------------------------------------")
print("Summary of the derived preprocessing configuration")
print(prep)

print('-------------------------------------------------------------------')
print(' ... or for individual sections')
print(prep.ibm)

print('-------------------------------------------------------------------')
print("Accessing section parameters...")
# Example access patterns.
print(f"ibm dx: {prep.ibm.dx}")
print(f"ic thl0: {prep.ic.thl0}")
print(f"vegetation ltrees: {prep.vegetation.ltrees}")
print(f"scalars nsv: {prep.scalars.nsv}")
print(f"seb facT: {prep.seb.facT}")

print('-------------------------------------------------------------------')
print("Example 1: override IC parameter and write profiles...")
# (Here we adjust the initial potential temperature before writing profiles.)
prep.ic.thl0 = 290.0
prep.ic.generate_prof()
prep.ic.write_prof()

print('-------------------------------------------------------------------')
print("Example 2: enable vegetation and run vegetation preprocessing...")
prep.vegetation.ltrees = 1
prep.vegetation.run_all()

print('-------------------------------------------------------------------')
print("Example 3: run IBM section preprocessing...")
prep.ibm.run_all()

print('-------------------------------------------------------------------')
print("Example 4: write changed params for IC section...")
# This will eventually actually write the changed params.
prep.ic.write_changed_params()

print('-------------------------------------------------------------------')
print("Example 4b: show changed params for all sections...")
prep.show_changed_params()

print('-------------------------------------------------------------------')
print("Example 4c: write changed params for all sections...")
prep.write_changed_params()

print('-------------------------------------------------------------------')
print("Example 5: run full preprocessing chain...")
prep.run_all()
