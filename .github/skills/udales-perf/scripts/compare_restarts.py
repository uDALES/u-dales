#!/usr/bin/env python3
"""Compare two uDALES restart-file sets (Class-B verification, issue #330).

Usage: compare_restarts.py <dirA> <dirB>

Reads every per-rank initd* file in dirA and its namesake in dirB
(Fortran unformatted sequential, double precision; record layout per
modsave.f90::writerestartfiles) and reports the relative L-inf difference
per field. Exit code 1 if any field exceeds TOL.
"""
import sys
import glob
import os
import numpy as np
from scipy.io import FortranFile

TOL = 1e-10  # relative L-inf tolerance for flag-only (Class B) changes over ~10 steps

# record order in writerestartfiles (after mindist and wall, which are static geometry)
FIELDS = ["mindist", "wall", "u0", "v0", "w0", "pres0", "thl0", "e120",
          "ekm", "qt0", "ql0", "ql0h"]


def read_records(path):
    recs = []
    with FortranFile(path, "r") as f:
        for _ in FIELDS:
            recs.append(f.read_reals(dtype=np.float64))
        recs.append(f.read_reals(dtype=np.float64))  # timee, dt
    return recs


def main(dira, dirb):
    files = sorted(os.path.basename(p) for p in glob.glob(os.path.join(dira, "initd*")))
    if not files:
        sys.exit(f"no initd* files in {dira}")
    partner_missing = [f for f in files if not os.path.exists(os.path.join(dirb, f))]
    if partner_missing:
        sys.exit(f"missing in {dirb}: {partner_missing[:3]} ...")

    worst = {name: 0.0 for name in FIELDS}
    for fname in files:
        ra = read_records(os.path.join(dira, fname))
        rb = read_records(os.path.join(dirb, fname))
        for name, a, b in zip(FIELDS, ra, rb):
            denom = max(np.abs(a).max(), np.abs(b).max(), 1e-300)
            diff = np.abs(a - b).max() / denom
            worst[name] = max(worst[name], diff)

    print(f"{'field':8s}  rel Linf (worst over {len(files)} ranks)")
    fail = False
    for name in FIELDS:
        status = ""
        if worst[name] > TOL:
            status = "  <-- EXCEEDS TOL" if name not in ("mindist", "wall") else "  <-- GEOMETRY DIFFERS?!"
            fail = True
        print(f"{name:8s}  {worst[name]:.3e}{status}")
    print(f"\nTOL = {TOL:g}: {'FAIL' if fail else 'PASS'}")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
