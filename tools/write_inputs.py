#%% Import libraries
from __future__ import annotations

import sys
from pathlib import Path


#%% Add u-dales/tools/python to path for imports (this file is in u-dales/tools/)
script_dir = Path(__file__).resolve().parent
tools_python = script_dir / "python"
if str(tools_python) not in sys.path:
    sys.path.insert(0, str(tools_python))


#%% Key Inputs: Setup experiment directory and number
if len(sys.argv) > 1:
    # Set from command-line argument
    expdir = Path(sys.argv[1]).resolve()
else:
    # Default behavior: use hardcoded paths
    expnr = "999"
    # Derive udales root from script location (tools/ -> u-dales/)
    udales_root = script_dir.parent
    expdir = (udales_root.parent / "experiments" / expnr).resolve()


#%% Initialize UDPrep
from udprep import UDPrep  # noqa: E402

print("Initializing UDPrep...")
sys.stdout.flush()

# Lightweight approach (default) - only reads namoptions + STL
prep = UDPrep(expdir)

# Alternative: Use legacy UDBase mode
# prep = UDPrep(expdir, use_udbase=True)


#%% Display derived configuration
print("-------------------------------------------------------------------")
print("Summary of the derived preprocessing configuration")
print(prep)

print('-------------------------------------------------------------------')
print(' ... or for individual sections')
print(prep.grid)

print('-------------------------------------------------------------------')
print("Accessing section parameters...")
# Example access patterns.
print(f"grid dx: {prep.grid.dx}")
print(f"grid itot: {prep.grid.itot}")
print(f"ic thl0: {prep.ic.thl0}")
print(f"vegetation ltrees: {prep.vegetation.ltrees}")
print(f"scalars nsv: {prep.scalars.nsv}")
print(f"seb facT: {prep.seb.facT}")

print('-------------------------------------------------------------------')
print("Example 0: run grid generation...")
prep.grid.run_all()

print('-------------------------------------------------------------------')
print("Example 1: override IC parameter and write profiles...")
# (Here we adjust the initial potential temperature before writing profiles.)
prep.ic.thl0 = 290.0
prep.ic.generate_prof()
prep.ic.write_prof()

# print('-------------------------------------------------------------------')
# print("Example 2: enable vegetation and run vegetation preprocessing...")
# prep.vegetation.ltrees = 1
# prep.vegetation.run_all()

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
