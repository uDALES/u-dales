#!/usr/bin/env python3
"""Driver for the in-solver nesting unit tests (runmodes 1006-1011).

Builds the fixtures with the production Python writer
(``tools/python/udprep/nesting.py``), then runs each nesting runmode on
1x1, 2x1, 1x2 and 2x2 ranks and asserts the exit code. The abort cases
(U13, U22, U27, and the schema-2 cases) are separate invocations that must exit
non-zero with a specific message, because the routine under test calls
``stop 1``.

See README.md in this directory for what each runmode covers.

Environment:
  UDALES_BUILD            path to the u-dales executable
                          (default build/release/u-dales)
  UDALES_RUNTIME_MODULES  module stack loaded before the run
  MPIEXEC                 MPI launcher
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
from typing import Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[3]
TEST_DIR = Path(__file__).resolve().parent
EXPNR = "901"

UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0",
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
    1008: "io (U15-U22)",
    1009: "flux (U23-U28, U35-U39)",
    1010: "update (U29-U34)",
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
    stripped = (text or "").strip()
    if not stripped:
        return f"{label}: <empty>"
    lines = stripped.splitlines()
    return f"{label} (last {min(len(lines), limit)} lines):\n" + "\n".join(lines[-limit:])


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
        failures: List[str] = []
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
            if done.returncode == 0:
                failures.append(f"{label}: expected a non-zero exit, got 0")
            elif needle not in output:
                failures.append(
                    f"{label}: aborted but without the expected message {needle!r}\n"
                    + _tail("output", output)
                )
        if failures:
            self.fail("\n\n".join(failures))

    def test_reported_messages(self) -> None:
        failures: List[str] = []
        cache: Dict[str, str] = {}
        for n, (label, runmode, edits, needle) in enumerate(OUTPUT_CASES):
            key = f"{runmode}:{sorted(edits.items())}"
            if key not in cache:
                all_edits = dict(edits)
                all_edits.update({"nprocx": "1", "nprocy": "1"})
                name = _write_namelist(
                    self.run_dir, runmode, f"namoptions.{runmode}.msg{n}", all_edits
                )
                done = _run(self.run_dir, name, 1)
                cache[key] = (done.stdout or "") + (done.stderr or "")
            if needle not in cache[key]:
                failures.append(f"{label}: runmode {runmode} did not report {needle!r}")
        if failures:
            self.fail("\n\n".join(failures))


if __name__ == "__main__":
    unittest.main()
