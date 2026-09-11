#!/usr/bin/env python3
"""Regression test for the driver-writer record-continuation fix on `nesting`
(src/moddriver.f90, src/modstartup.f90).

Before the fix, a `&DRIVER idriver = 1` writer that was warm-started always
began its inflow-plane record again at index 1 (overwriting the start of the
record it had already partly written), and `trestart` was unconditionally
replaced by `(tdriverstart + (driverstore-1)*dtdriver)` -- so a record longer
than one job could never be written across a restart. The fix: (a) a namelist
`trestart` shorter than the record length is honoured; (b) on a warm start at
or after the record start, if `tdriver_000.<expnr>` already exists, its record
count is resumed (`driver_records_written`) and the nominal dump clock keeps
counting from the cold `tdriverstart`, logging "Driver record continued: N
records already written".

This test builds a small, periodic (no IBM), fixed-timestep case -- based on
tests/integration/nesting/namoptions.1012's grid/BC and prof.inp.901, see
namoptions.903 for what was changed and why -- and runs the solver three
times:

  A (one job):  the whole driverstore record in a single run, trestart huge.
  B1 (job 1):   a cold start that stops partway through the record, with a
                namelist trestart short enough that a restart file is written
                mid-record.
  B2 (job 2):   a warm start from that restart file, running past the end of
                the record.

then asserts that B's tdriver_/udriver_/vdriver_/wdriver_ files are identical
to A's: with a fixed timestep, deterministic (seeded) dynamics and a bitwise
restart round-trip, the two must produce the same numbers whether or not the
run was interrupted and resumed mid-record. If that ever turns out not to be
bit-exact here, the comparison falls back to a numeric max-abs-diff against
TIGHT_TOLERANCE below -- chosen deliberately tight, so a real behavioural
regression still fails loudly while reporting the actual discrepancy instead
of just "not equal".

Environment: same as tests/integration/nesting/_launch.py (UDALES_BUILD,
UDALES_RUNTIME_MODULES, UDALES_MPIEXEC, UDALES_REQUIRE_LAUNCHER, TMPDIR).
"""

from __future__ import annotations

import array
import filecmp
import os
import re
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List

TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parents[2]

# Reuse the nesting suite's launch helper (module stack, mpiexec discovery,
# launcher skip/require logic) rather than forking a copy of it.
sys.path.insert(0, str(REPO_ROOT / "tests" / "integration" / "nesting"))
import _launch as launch  # noqa: E402

EXPNR = "903"
UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))

# Real is promoted to 8 bytes (-r8 / -fdefault-real-8, see CMakeLists.txt), so
# each direct-access record in tdriver_*.<expnr> (one real, no header) is 8
# bytes.
BYTES_PER_RECORD = 8

# Namelist timings (all in the &DRIVER/&RUN sense: seconds of simulated time).
DTMAX = 0.05          # fixed timestep (ladaptive = .false.)
TDRIVERSTART = 0.25   # a few steps in
DTDRIVER = 0.05       # == dtmax: one driver record dumped per step
DRIVERSTORE = 40
RECORD_END = TDRIVERSTART + (DRIVERSTORE - 1) * DTDRIVER  # 2.20: last record's time
TAIL_STEPS = 2        # steps run past the record end, both in A and in B2

REFERENCE_RUNTIME = RECORD_END + TAIL_STEPS * DTMAX  # 2.30

# B1 stops mid-record with a namelist trestart shorter than the record: the
# restart file is only written when timee >= tnextrestart at rk3step 3 (see
# writerestartfiles in src/modsave.f90), so B1's runtime must run a little
# past SEG1_TRESTART for the write to actually happen.
SEG1_TRESTART = 1.20
SEG1_RUNTIME = 1.25

# B2 warm-starts from that restart (btime == SEG1_TRESTART, since dt is fixed
# and lands exactly on it) and runs to the same absolute end time as A.
SEG2_RUNTIME = REFERENCE_RUNTIME - SEG1_TRESTART

# Fallback tolerance if the byte-identical comparison ever fails: see the
# module docstring. 0 is what this test actually expects.
TIGHT_TOLERANCE = 1.0e-9


def _read_doubles(path: Path) -> array.array:
    data = path.read_bytes()
    n = len(data) // 8
    a = array.array("d")
    a.frombytes(data[: n * 8])
    return a


def _max_abs_diff(path_a: Path, path_b: Path) -> float:
    a = _read_doubles(path_a)
    b = _read_doubles(path_b)
    n = min(len(a), len(b))
    if n == 0:
        return float("nan")
    return max(abs(x - y) for x, y in zip(a[:n], b[:n]))


def _tail(label: str, text: str) -> str:
    return launch.tail(label, text)


class DriverRestartTest(unittest.TestCase):
    """The driver writer must continue its record across a warm start."""

    @classmethod
    def setUpClass(cls) -> None:
        if not UDALES_BUILD.is_file():
            raise RuntimeError(f"u-dales executable not found at {UDALES_BUILD}")
        launch.require_launcher()
        cls.workdir = tempfile.TemporaryDirectory(prefix="udales-driverrestart-",
                                                    dir=launch.scratch_dir())
        cls.base_dir = Path(cls.workdir.name)

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "workdir"):
            cls.workdir.cleanup()

    def _make_run_dir(self, name: str) -> Path:
        d = self.base_dir / name
        d.mkdir(parents=True)
        shutil.copy2(TEST_DIR / f"prof.inp.{EXPNR}", d)
        shutil.copy2(TEST_DIR / f"lscale.inp.{EXPNR}", d)
        return d

    def _write_namelist(self, run_dir: Path, name: str, edits: Dict[str, str]) -> str:
        text = (TEST_DIR / f"namoptions.{EXPNR}").read_text(encoding="utf-8")
        for key, value in edits.items():
            text, n = re.subn(
                rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$", lambda m: m.group(1) + value, text
            )
            if n == 0:
                raise RuntimeError(f"setting '{key}' not found in namoptions.{EXPNR}")
        (run_dir / name).write_text(text, encoding="utf-8")
        return name

    def _run(self, run_dir: Path, name: str, nprocs: int):
        return launch.run(run_dir, name, nprocs, UDALES_BUILD)

    def _record_count(self, path: Path) -> int:
        return path.stat().st_size // BYTES_PER_RECORD

    def _run_protocol(self, label: str, nprocx: int, nprocy: int) -> None:
        nprocs = nprocx * nprocy
        dir_a = self._make_run_dir(f"{label}-A")
        dir_b = self._make_run_dir(f"{label}-B")
        decomp = {"nprocx": str(nprocx), "nprocy": str(nprocy)}

        # --- A: the whole record, one job ------------------------------- #
        name = self._write_namelist(dir_a, f"namoptions.{EXPNR}.A", {
            **decomp,
            "runtime": f"{REFERENCE_RUNTIME:.10g}",
            "trestart": "1.e9",
            "lwarmstart": ".false.",
        })
        done_a = self._run(dir_a, name, nprocs)
        if done_a.returncode != 0:
            self.fail(f"{label}: reference run A exited {done_a.returncode}\n"
                      + _tail("stdout", done_a.stdout) + "\n" + _tail("stderr", done_a.stderr))

        # --- B1: cold start, stopped mid-record with a short trestart --- #
        name = self._write_namelist(dir_b, f"namoptions.{EXPNR}.B1", {
            **decomp,
            "runtime": f"{SEG1_RUNTIME:.10g}",
            "trestart": f"{SEG1_TRESTART:.10g}",
            "lwarmstart": ".false.",
        })
        done_b1 = self._run(dir_b, name, nprocs)
        out_b1 = (done_b1.stdout or "") + (done_b1.stderr or "")
        if done_b1.returncode != 0:
            self.fail(f"{label}: segment 1 (B1) exited {done_b1.returncode}\n" + _tail("output", out_b1))

        self.assertNotIn(
            "trestart mentioned in namoptions", out_b1,
            f"{label}: B1's namelist trestart ({SEG1_TRESTART}) was overridden -- "
            f"it is shorter than the record and should have been honoured\n" + _tail("output", out_b1),
        )

        restarts = sorted(dir_b.glob(f"initd????????_???_???.{EXPNR}"))
        if len(restarts) != nprocs:
            self.fail(f"{label}: expected {nprocs} restart file(s) (one per rank) after B1, "
                      f"found {[p.name for p in restarts]}\n" + _tail("output", out_b1))
        startfile = restarts[0].name

        record_count_before = self._record_count(dir_b / f"tdriver_000.{EXPNR}")
        self.assertGreater(record_count_before, 0,
                            f"{label}: no driver records on disk after B1")
        self.assertLess(record_count_before, DRIVERSTORE,
                         f"{label}: B1 already wrote the whole record ({record_count_before}); "
                         "the restart would not be mid-record any more, tighten SEG1_TRESTART/SEG1_RUNTIME")

        # --- B2: warm start from that restart, past the record end ------ #
        name = self._write_namelist(dir_b, f"namoptions.{EXPNR}.B2", {
            **decomp,
            "runtime": f"{SEG2_RUNTIME:.10g}",
            "trestart": "1.e9",
            "lwarmstart": ".true.",
            "startfile": f"'{startfile}'",
        })
        done_b2 = self._run(dir_b, name, nprocs)
        out_b2 = (done_b2.stdout or "") + (done_b2.stderr or "")
        if done_b2.returncode != 0:
            self.fail(f"{label}: segment 2 (B2) exited {done_b2.returncode}\n" + _tail("output", out_b2))

        m = re.search(r"Driver record continued:\s*(\d+)\s*records already written", out_b2)
        if not m:
            self.fail(f"{label}: expected the 'Driver record continued:' log line in B2 -- "
                      + _tail("output", out_b2))
        self.assertEqual(
            int(m.group(1)), record_count_before,
            f"{label}: B2 reported continuing from {m.group(1)} records, "
            f"but {record_count_before} were on disk when it started",
        )

        # --- Byte parity between A and B for every driver file ---------- #
        # tdriver_*.<expnr> is written only by rank column 0 (writedriverfile
        # gates it on driverid==0, since the time stamps are the same for
        # every column); udriver_/vdriver_/wdriver_ are per rank column.
        mismatches: List[str] = []
        driver_files = [("t", 0)] + [(v, r) for v in "uvw" for r in range(nprocy)]
        for var, r in driver_files:
            fname = f"{var}driver_{r:03d}.{EXPNR}"
            fa, fb = dir_a / fname, dir_b / fname
            self.assertTrue(fa.is_file(), f"{label}: missing {fa}")
            self.assertTrue(fb.is_file(), f"{label}: missing {fb}")
            if not filecmp.cmp(fa, fb, shallow=False):
                diff = _max_abs_diff(fa, fb)
                if diff > TIGHT_TOLERANCE:
                    mismatches.append(f"{fname}: max abs diff {diff:.3e} > tolerance {TIGHT_TOLERANCE:g}")
        if mismatches:
            self.fail(f"{label}: A (one job) and B (restarted) driver files differ:\n  - "
                      + "\n  - ".join(mismatches))

        # --- Record count from file size equals driverstore in both ----- #
        for tag, d in (("A", dir_a), ("B", dir_b)):
            n = self._record_count(d / f"tdriver_000.{EXPNR}")
            self.assertEqual(n, DRIVERSTORE,
                              f"{label}: {tag} has {n} driver records on disk, "
                              f"expected driverstore = {DRIVERSTORE}")

        # --- B's time records strictly increase (the old bug reset/duped) #
        times = _read_doubles(dir_b / f"tdriver_000.{EXPNR}")
        self.assertEqual(len(times), DRIVERSTORE)
        for i in range(1, len(times)):
            self.assertGreater(
                times[i], times[i - 1],
                f"{label}: driver time record {i} ({times[i]}) does not strictly "
                f"increase over record {i - 1} ({times[i - 1]}) -- a duplicated or reset record "
                "is exactly the pre-fix bug",
            )

    def test_one_rank(self) -> None:
        self._run_protocol("1x1", nprocx=1, nprocy=1)

    def test_two_ranks_in_y(self) -> None:
        """Cheap to add: the plane files are written one per rank column in y."""
        self._run_protocol("1x2", nprocx=1, nprocy=2)


if __name__ == "__main__":
    unittest.main()
