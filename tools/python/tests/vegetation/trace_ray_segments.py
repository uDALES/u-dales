from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

udbase_path = Path(
    r"C:/Users/mvr/OneDrive - Imperial College London/codes/uDALES/u-dales"
).resolve()
tools_path = (udbase_path / "tools" / "python").resolve()
if str(tools_path) not in sys.path:
    sys.path.insert(0, str(tools_path))

from udbase import UDBase  # noqa: E402
from udprep.directshortwave import trace_ray_segments  # noqa: E402


def main() -> None:
    base = Path(__file__).parent
    detail_path = base / "expected_ray_detail.csv"
    if not detail_path.exists():
        raise SystemExit(f"Missing {detail_path}. Run expected_ray_analysis.py first.")
    rows = np.genfromtxt(detail_path, delimiter=",", names=True, dtype=None, encoding=None)
    if rows.size == 0:
        raise SystemExit("expected_ray_detail.csv is empty.")

    expnr = "064"
    expdir = (udbase_path.parents[0] / "experiments" / expnr).resolve
    sim = UDBase(expnr, expdir)

    for label in ("wall_low", "wall_high", "roof"):
        sel = rows["label"] == label
        if not np.any(sel):
            continue
        row = rows[sel][0]
        origin = np.array([row["start_x"], row["start_y"], row["start_z"]], dtype=float)
        direction = np.array([row["dir_x"], row["dir_y"], row["dir_z"]], dtype=float)
        segs = trace_ray_segments(sim, origin, direction, max_steps=4000, extend_bounds=True)
        out_path = base / f"ray_segments_{label}.csv"
        header = "i,j,k,t,t_next,t_max_x,t_max_y,t_max_z,ds,inside"
        np.savetxt(out_path, segs, delimiter=",", header=header, comments="")
        print(f"Wrote {out_path} ({segs.shape[0]} rows)")


if __name__ == "__main__":
    main()
