from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/rds/general/user/mvr/home/udales/u-dales/tools/python/")

from udbase import UDBase  # noqa: E402
from udprep import convert_block_to_sparse  # noqa: E402

sim_id = 525
base_path = Path("/rds/general/user/mvr/home/udales/u-dales/tests/tests_tree_input")

# Instantiate UDBase without loading STL geometry
sim = UDBase(expnr=sim_id, path=base_path, load_geometry=False)
if sim.trees is None:
    raise RuntimeError(f"trees.inp.{sim_id} not loaded from {base_path}")

out_file = convert_block_to_sparse(sim)
print(f"Conversion complete: {out_file}")
#return out_file
