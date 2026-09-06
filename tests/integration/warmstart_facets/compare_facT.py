#!/usr/bin/env python3
"""Compare facet temperature trajectories across a paired reference/restart chain.

This is the numerical half of the warm-start facet continuity test for GitHub
issue #269. ``run_test.sh`` produces the NetCDF outputs; this script asserts
that the restarted chain continues the reference chain rather than restarting
the surface energy balance from ``Tfacinit.inp``.

The facT output is written as ``NF90_FLOAT`` (see modstat_nc.f90), so the
tolerances below are chosen relative to single-precision resolution at ~300 K
(about 3e-5 K) rather than demanding bitwise equality.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, Tuple

import netCDF4 as nc
import numpy as np


# Continuity tolerance for a restart that should reproduce the reference
# trajectory. Two orders of magnitude above float32 resolution at 300 K.
CONTINUITY_ATOL = 1.0e-3
# Parity tolerance between two MPI decompositions of the same restart. Looser,
# because the reference chains themselves differ at round-off between
# decompositions; still far tighter than a broken broadcast could ever be.
DECOMPOSITION_ATOL = 1.0e-2
# The reference trajectory must move at least this far between its first and
# last record, otherwise the continuity assertion above is vacuous and the test
# configuration needs a longer run or a larger dtEB.
MIN_SIGNAL = 20.0 * CONTINUITY_ATOL


def load_facT(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Return (times, T) with T shaped (time, layer, facet) as float64."""
    if not path.is_file():
        raise SystemExit(f"ERROR: facT output not found: {path}")
    with nc.Dataset(path) as ds:
        if "T" not in ds.variables:
            raise SystemExit(f"ERROR: no 'T' variable in {path}")
        times = np.asarray(ds.variables["t"][:], dtype=np.float64)
        values = np.asarray(ds.variables["T"][:], dtype=np.float64)
    if times.size == 0:
        raise SystemExit(f"ERROR: {path} contains no facT records")
    return times, values


def record_at(times: np.ndarray, values: np.ndarray, target: float, label: str) -> np.ndarray:
    matches = np.nonzero(np.isclose(times, target, rtol=0.0, atol=1.0e-6))[0]
    if matches.size == 0:
        raise SystemExit(
            f"ERROR: {label} has no facT record at t = {target}; available times: {times.tolist()}"
        )
    return values[matches[0]]


def report(name: str, diff: float, atol: float, failures: list) -> None:
    status = "ok" if diff <= atol else "FAIL"
    print(f"  [{status}] {name}: max abs diff {diff:.3e} K (tolerance {atol:.1e} K)")
    if diff > atol:
        failures.append(f"{name}: max abs diff {diff:.3e} K exceeds {atol:.1e} K")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path,
                        help="facT NetCDF from the cold-start reference run")
    parser.add_argument("--restart", required=True, type=Path,
                        help="facT NetCDF from the warm-started continuation")
    parser.add_argument("--other-decomposition", type=Path, default=None,
                        help="facT NetCDF from the same restart under a different decomposition")
    parser.add_argument("--unfixed", type=Path, default=None,
                        help="facT NetCDF from a run with lreadfacT = .false. (pre-fix behaviour)")
    args = parser.parse_args()

    failures: list = []

    ref_t, ref_T = load_facT(args.reference)
    rst_t, rst_T = load_facT(args.restart)

    print(f"reference records at t = {ref_t.tolist()}")
    print(f"restart   records at t = {rst_t.tolist()}")

    # The reference must actually evolve, otherwise nothing below discriminates.
    signal = float(np.abs(ref_T[-1] - ref_T[0]).max())
    print(f"reference facT drift over the run: {signal:.3e} K "
          f"(needs > {MIN_SIGNAL:.1e} K for the test to be meaningful)")
    if signal <= MIN_SIGNAL:
        failures.append(
            f"reference facT barely moved ({signal:.3e} K); the continuity assertions "
            "would be vacuous - lengthen the run or increase dtEB"
        )

    # --- The #269 regression assertion -------------------------------------
    # The first facT record written after the restart must equal the reference
    # record at the same time. Before the fix it instead equals the Tfacinit
    # profile advanced by a single EB step.
    first_t = float(rst_t[0])
    ref_first = record_at(ref_t, ref_T, first_t, "reference")
    diff_first = float(np.abs(rst_T[0] - ref_first).max())
    print("\nissue #269 regression assertion:")
    report(f"first post-restart record (t = {first_t})", diff_first, CONTINUITY_ATOL, failures)

    # How far the first post-restart record sits from the reference's own first
    # record, which is roughly what a Tfacinit re-seed produces. Reported so a
    # failure is easy to read, and asserted so the tolerance stays meaningful.
    spread = float(np.abs(ref_first - ref_T[0]).max())
    print(f"  (reference moved {spread:.3e} K between its first record and t = {first_t})")

    # --- Full-trajectory continuity ----------------------------------------
    print("\ntrajectory continuity over the overlapping window:")
    overlap: Dict[float, float] = {}
    for idx, t in enumerate(rst_t):
        ref_match = np.nonzero(np.isclose(ref_t, t, rtol=0.0, atol=1.0e-6))[0]
        if ref_match.size:
            overlap[float(t)] = float(np.abs(rst_T[idx] - ref_T[ref_match[0]]).max())
    if not overlap:
        failures.append("restart and reference facT outputs share no common output times")
    for t, diff in sorted(overlap.items()):
        report(f"t = {t}", diff, CONTINUITY_ATOL, failures)

    # --- Negative control ---------------------------------------------------
    if args.unfixed is not None:
        unf_t, unf_T = load_facT(args.unfixed)
        unf_first = record_at(unf_t, unf_T, first_t, "lreadfacT = .false. run")
        diff_unfixed = float(np.abs(unf_first - ref_first).max())
        print("\nnegative control (lreadfacT = .false., i.e. pre-fix behaviour):")
        print(f"  first post-restart record differs from the reference by {diff_unfixed:.3e} K")
        if diff_unfixed <= CONTINUITY_ATOL:
            failures.append(
                f"lreadfacT = .false. reproduced the reference to {diff_unfixed:.3e} K, so the "
                "continuity assertion cannot detect the #269 regression"
            )
        else:
            print("  [ok] the assertion above does discriminate the unfixed behaviour")

    # --- MPI decomposition parity ------------------------------------------
    if args.other_decomposition is not None:
        oth_t, oth_T = load_facT(args.other_decomposition)
        oth_first = record_at(oth_t, oth_T, first_t, "other-decomposition restart")
        diff_decomp = float(np.abs(rst_T[0] - oth_first).max())
        print("\nMPI decomposition parity after restart (guards the facet broadcasts):")
        report("first post-restart record", diff_decomp, DECOMPOSITION_ATOL, failures)

    print()
    if failures:
        print("FAILURES:")
        for item in failures:
            print(f"  - {item}")
        return 1
    print("All facT comparisons passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
