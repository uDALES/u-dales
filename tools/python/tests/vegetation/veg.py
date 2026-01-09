from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/rds/general/user/mvr/home/udales/u-dales/tools/python/")

from udbase import UDBase  # noqa: E402
from udprep import UDPrep  # noqa: E402

sim_id = 525
base_path = Path("/rds/general/user/mvr/home/udales/u-dales/tests/tests_tree_input")

# Instantiate UDBase without loading STL geometry
prep = UDPrep(expnr=sim_id, path=base_path, load_geometry=False)
out_file = prep.vegetation.load_block()
print(f"Conversion complete: {out_file}")
#return out_file
