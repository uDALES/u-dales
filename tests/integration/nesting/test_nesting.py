#!/usr/bin/env python3
"""Driver for the in-solver nesting unit tests (runmodes 1006-1011).

Builds the fixtures with the production Python writer
(``tools/python/udprep/nesting.py``), then runs each nesting runmode on
1x1, 2x1, 1x2 and 2x2 ranks and asserts the exit code. The abort cases
(U13, U22, U27, and the schema-2 cases) are separate invocations that must exit
non-zero with a specific message, because the routine under test calls
``nest_abort``.

See README.md in this directory for what each runmode covers.

Environment (see ``_launch.py`` for the full list and the defaults):
  UDALES_BUILD            path to the u-dales executable
                          (default build/release/u-dales)
  UDALES_RUNTIME_MODULES  module stack loaded before the run
  UDALES_MPIEXEC          MPI launcher (then MPIEXEC, then PATH)
  UDALES_REQUIRE_LAUNCHER =1: an unusable launcher fails instead of skipping
  UDALES_ABORT_EXIT_CODE  exit code the abort cases must return (default 1)
  TMPDIR                  where the run directories go
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Tuple

TEST_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TEST_DIR))

import _launch as launch  # noqa: E402

REPO_ROOT = launch.REPO_ROOT
EXPNR = "901"

UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))

#: What ``mpiexec`` returns when a rank calls ``nest_abort``.  Measured 1 for
#: both ``MPI_Abort(comm, 1)`` and ``stop 1`` under Intel MPI 2021.2 and
#: Open MPI 4.1.5, on 1 and 2 ranks.  The code an ``MPI_Abort`` propagates is
#: the launcher's business, so a launcher that maps it differently can say so
#: here; whatever the value, it must be non-zero and must not be a signal
#: death, and the message must be there.
ABORT_EXIT_CODE = int(os.environ.get("UDALES_ABORT_EXIT_CODE", "1"))

#: Output that means the process crashed rather than aborted on purpose.  An
#: abort test that accepted these would pass on a segfault in the code path it
#: is meant to be checking.
CRASH_SIGNATURES = (
    "forrtl: severe", "forrtl: error", "Segmentation fault", "SIGSEGV", "SIGFPE",
    "Floating point exception", "floating point exception", "Backtrace for this error",
    "Program received signal", "exited on signal", "BAD TERMINATION",
)

CONFIGS: Dict[str, Tuple[int, int]] = {
    "serial": (1, 1),
    "x_split": (2, 1),
    "y_split": (1, 2),
    "xy_split": (2, 2),
}

RUNMODES = {
    1006: "weights (U1-U7)",
    1007: "geometry (U8-U14)",
    1008: "io (U15-U19, U21-U22, U44, U46-U47)",
    1009: "flux (U23-U28, U35-U39, U48)",
    1010: "update (U29-U34, U45, U49)",
    1011: "cold-start init from the parent (U40-U43)",
}

#: Cases that must abort. Each entry is
#: (label, runmode, namelist edits, a string the output must contain).
ABORT_CASES = [
    (
        "U13 building-free enforcement",
        1007,
        {"nest_lparentgeom": ".false."},
        "solid points found inside the relaxation zone",
    ),
    (
        "U22 header validation: schema",
        1008,
        {"nestfile": f"'bad_schema.{EXPNR}.nc'"},
        "mismatch in udales_nesting_schema",
    ),
    (
        "U22 header validation: itot",
        1008,
        {"nestfile": f"'bad_itot.{EXPNR}.nc'"},
        "mismatch in itot",
    ),
    (
        "U22 header validation: xlen",
        1008,
        {"nestfile": f"'bad_xlen.{EXPNR}.nc'"},
        "mismatch in xlen",
    ),
    (
        "U22 header validation: zf",
        1008,
        {"nestfile": f"'bad_zf.{EXPNR}.nc'"},
        "mismatch in zf",
    ),
    (
        "U22 header validation: stagger tag",
        1008,
        {"nestfile": f"'bad_stagger.{EXPNR}.nc'"},
        "stagger: file =",
    ),
    (
        "U27 flux assertion fires",
        1009,
        {"nestfile": f"'assertfire.{EXPNR}.nc'"},
        "not flux balanced",
    ),
    # design 10.6 item 3: nest_lfluxcheckall is what catches a writer whose
    # stored residual does not describe the data it stored.
    (
        "flux: a lying stored residual is caught by nest_lfluxcheckall",
        1009,
        {"nestfile": f"'assertfire_lying.{EXPNR}.nc'",
         "nest_lfluxcheckall": ".true."},
        "not flux balanced",
    ),
    # ... and the reader falls back to the recompute on its own when the fluid
    # lateral area the file advertises is not this run's.
    (
        "flux: an area mismatch forces the full recompute",
        1009,
        {"nestfile": f"'assertfire_area.{EXPNR}.nc'"},
        "does not match this run",
    ),
    # design 10.6 item 4: the cold-start switch and its three failure modes.
    (
        "init: switch on but the file has no initial-condition block",
        1011,
        {"nestfile": f"'nesting_analytic.{EXPNR}.nc'"},
        "carries no initial-condition block",
    ),
    (
        "init: initial-condition block at the wrong shape",
        1011,
        {"nestfile": f"'bad_initdims.{EXPNR}.nc'"},
        "mismatch in u_init",
    ),
    (
        "init: initial-condition block at the wrong stagger",
        1011,
        {"nestfile": f"'bad_initstag.{EXPNR}.nc'"},
        "u_init stagger: file =",
    ),
    # review 2026-09-06 items F3-F5. F5: the stored residual is clean, so the
    # cheap init check passes and the per-level slab validator must catch the
    # value. F4: nest_lendabort (the default) refuses at init a run that would
    # outlast its parent record, rather than freezing the boundary hours in.
    # F3: a nested direction must impose both of its faces; checkinitvalues
    # stops with `stop 1`, whose message survives only on rank 0 under Intel
    # MPI -- the abort cases run on 1 rank for that reason among others.
    (
        "F5 poisoned slab",
        1009,
        {"nestfile": f"'assertfire_nan.{EXPNR}.nc'"},
        "non-finite or fill",
    ),
    (
        "F4 run outlasts the record",
        1008,
        {"runtime": "1000."},
        "extends past the end of the parent record",
    ),
    (
        "F3 partial nest_lateral",
        1007,
        {"nest_lateral": ".true.,.false.,.true.,.true."},
        "needs both x faces",
    ),
    # A legal partial face set: x nested, y periodic, so lface = (T,T,F,F).
    # The stored flux_residual and fluid_lateral_area are sums over all FOUR
    # faces (the spec defines them that way and the writer has no notion of a
    # subset), so they certify nothing about the imposed pair -- zero over four
    # faces does not imply zero over west/east.  The cheap check used to accept
    # that four-face zero, pass initialisation, and abort on the first timestep
    # with a large Phi and no explanation.  It must now be caught at init, over
    # the faces actually imposed, and say why rerunning the correction will not
    # help.
    (
        "partial face set is not certified by the four-face residual",
        1009,
        {"BCym": "1"},
        "imposes only some of the four lateral faces",
    ),
]

#: Cases that must succeed and whose output must contain a given string.
#: Each entry is (label, runmode, namelist edits, needle).
OUTPUT_CASES = [
    ("U14 narrow-zone warning", 1007, {}, "WARNING zone is only"),
    ("U14 wide-zone warning", 1007, {}, "WARNING zone occupies"),
    ("U13 warning half", 1007, {}, "solid points inside the relaxation zone (allowed"),
    # design 10.6 item 3: the cheap path must say so, and a schema 1 file must
    # say why it is falling back rather than doing it silently.
    ("flux: the cheap check reports itself", 1009, {},
     "flux balanced (from the stored residual, no slab read)"),
    ("flux: a schema 1 file warns and recomputes", 1009, {},
     "predates schema 2"),
    ("flux: the recompute reports itself", 1009, {},
     "flux balanced (recomputed from the boundary slabs)"),
    # design 10.6 item 4.
    ("init: the cold start says where it came from", 1011, {},
     "cold start initialised from the parent"),
    ("init: a warm start says it is ignoring the switch", 1011, {},
     "the restart file wins"),
    # review 2026-09-06 items F4 and F2: the frozen-boundary warning under
    # nest_lendabort = .false. (U44), and the solid-face mask count (U45).
    ("F4 freeze warning", 1008, {}, "the boundary now freezes on that level"),
    ("F2 masked faces", 1010, {}, "masked west faces = 16 (0 still carried a value)"),
]


# --------------------------------------------------------------------------- #
# Running the solver
# --------------------------------------------------------------------------- #


def _run(run_dir: Path, namelist: str, nprocs: int) -> subprocess.CompletedProcess:
    return launch.run(run_dir, namelist, nprocs, UDALES_BUILD)


def _write_namelist(run_dir: Path, runmode: int, name: str, edits: Dict[str, str]) -> str:
    text = (TEST_DIR / f"namoptions.{runmode}").read_text(encoding="utf-8")
    for key, value in edits.items():
        text, n = re.subn(
            rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$", lambda m: m.group(1) + value, text
        )
        if n == 0:
            raise RuntimeError(f"setting '{key}' not found in namoptions.{runmode}")
    (run_dir / name).write_text(text, encoding="utf-8")
    return name


def _tail(label: str, text: str, limit: int = 60) -> str:
    return launch.tail(label, text, limit)


def crash_signature(output: str) -> str:
    """The first crash marker found in ``output``, or the empty string."""
    for needle in CRASH_SIGNATURES:
        if needle in output:
            return needle
    return ""


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #


class NestingUnitRunmodes(unittest.TestCase):
    """Every nesting runmode, on every decomposition, must exit 0."""

    @classmethod
    def setUpClass(cls) -> None:
        if not UDALES_BUILD.is_file():
            raise RuntimeError(f"u-dales executable not found at {UDALES_BUILD}")
        launch.require_launcher()

        cls.workdir = tempfile.TemporaryDirectory(prefix="udales-nesting-",
                                                  dir=launch.scratch_dir())
        cls.run_dir = Path(cls.workdir.name)
        shutil.copy2(TEST_DIR / f"prof.inp.{EXPNR}", cls.run_dir)

        # Fixtures come from the production writer, so the writer and the
        # Fortran reader cannot drift apart (U15/U16).
        import make_fixtures  # noqa: E402

        make_fixtures.write_all(cls.run_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "workdir"):
            cls.workdir.cleanup()

    def _run_case(self, runmode: int, label: str, nprocx: int, nprocy: int) -> subprocess.CompletedProcess:
        name = _write_namelist(
            self.run_dir,
            runmode,
            f"namoptions.{runmode}.{label}",
            {"nprocx": str(nprocx), "nprocy": str(nprocy)},
        )
        return _run(self.run_dir, name, nprocx * nprocy)

    def _assert_runmode(self, runmode: int, label: str, nprocx: int, nprocy: int) -> None:
        done = self._run_case(runmode, label, nprocx, nprocy)
        if done.returncode != 0:
            self.fail(
                f"runmode {runmode} ({RUNMODES[runmode]}) on {label} "
                f"({nprocx}x{nprocy}) exited {done.returncode}\n"
                + _tail("stdout", done.stdout)
                + "\n"
                + _tail("stderr", done.stderr)
            )

    # -- one method per runmode, serial: the unit layer ------------------- #

    def test_runmode_1006_weights(self) -> None:
        self._assert_runmode(1006, "serial", 1, 1)

    def test_runmode_1007_geometry(self) -> None:
        self._assert_runmode(1007, "serial", 1, 1)

    def test_runmode_1008_io(self) -> None:
        self._assert_runmode(1008, "serial", 1, 1)

    def test_runmode_1009_flux(self) -> None:
        self._assert_runmode(1009, "serial", 1, 1)

    def test_runmode_1010_update(self) -> None:
        self._assert_runmode(1010, "serial", 1, 1)

    def test_runmode_1011_init(self) -> None:
        self._assert_runmode(1011, "serial", 1, 1)

    # -- every runmode on every decomposition: the integration layer ------ #

    def test_runmodes_on_all_decompositions(self) -> None:
        failures: List[str] = []
        for runmode, what in sorted(RUNMODES.items()):
            for label, (nprocx, nprocy) in CONFIGS.items():
                done = self._run_case(runmode, label, nprocx, nprocy)
                if done.returncode != 0:
                    failures.append(
                        f"runmode {runmode} ({what}) on {label} "
                        f"({nprocx}x{nprocy}) exited {done.returncode}\n"
                        + _tail("stdout", done.stdout)
                        + "\n"
                        + _tail("stderr", done.stderr)
                    )
        if failures:
            self.fail("\n\n".join(failures))

    def test_abort_cases(self) -> None:
        """Each abort case must stop *deliberately*: the documented exit code,
        the documented message, and no sign of a crash.

        "Any non-zero exit containing the needle" would also accept a run that
        printed the message and then died of a bounds error on the way out;
        the exit code and the crash signatures are what separate the two.
        """
        failures: List[str] = []
        codes: Dict[str, int] = {}
        for n, (label, runmode, edits, needle) in enumerate(ABORT_CASES):
            edits = dict(edits)
            edits.update({"nprocx": "1", "nprocy": "1"})
            name = _write_namelist(
                self.run_dir,
                runmode,
                f"namoptions.{runmode}.abort{n}",
                edits,
            )
            done = _run(self.run_dir, name, 1)
            output = (done.stdout or "") + (done.stderr or "")
            codes[label] = done.returncode
            if done.returncode == 0:
                failures.append(f"{label}: expected a non-zero exit, got 0")
                continue
            problems = []
            if needle not in output:
                problems.append(f"aborted but without the expected message {needle!r}")
            if done.returncode >= 128:
                problems.append(f"exit code {done.returncode} means a signal death, "
                                "not a deliberate abort")
            elif done.returncode != ABORT_EXIT_CODE:
                problems.append(f"exit code {done.returncode}, expected {ABORT_EXIT_CODE} "
                                "(UDALES_ABORT_EXIT_CODE) -- if this launcher maps "
                                "MPI_Abort's code differently, say so in the environment")
            marker = crash_signature(output)
            if marker:
                problems.append(f"the output carries a crash signature ({marker!r}), "
                                "so this is not the abort path being tested")
            if problems:
                failures.append(f"{label}:\n  - " + "\n  - ".join(problems)
                                + "\n" + _tail("output", output))
        print(f"\n[abort] exit codes: {sorted(set(codes.values()))} over "
              f"{len(codes)} cases (expected {ABORT_EXIT_CODE})", flush=True)
        if failures:
            self.fail("\n\n".join(failures))

    def test_reported_messages(self) -> None:
        """Each message case must both exit 0 and print its message."""
        failures: List[str] = []
        cache: Dict[str, Tuple[int, str]] = {}
        for n, (label, runmode, edits, needle) in enumerate(OUTPUT_CASES):
            key = f"{runmode}:{sorted(edits.items())}"
            if key not in cache:
                all_edits = dict(edits)
                all_edits.update({"nprocx": "1", "nprocy": "1"})
                name = _write_namelist(
                    self.run_dir, runmode, f"namoptions.{runmode}.msg{n}", all_edits
                )
                done = _run(self.run_dir, name, 1)
                cache[key] = (done.returncode, (done.stdout or "") + (done.stderr or ""))
            returncode, output = cache[key]
            if returncode != 0:
                failures.append(f"{label}: runmode {runmode} exited {returncode}\n"
                                + _tail("output", output))
            elif needle not in output:
                failures.append(f"{label}: runmode {runmode} did not report {needle!r}")
        if failures:
            self.fail("\n\n".join(failures))


if __name__ == "__main__":
    unittest.main()
