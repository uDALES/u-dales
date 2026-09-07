#!/usr/bin/env python3
"""Driver for the nesting integration matrix I1-I8.

Implements ``docs/udales-nesting-design.md`` section 10.3.  Where the unit
runmodes (``test_nesting.py``, U1-U34) check one mechanism at a time through
the public test hooks, these tests run the **whole solver** on a nested case
and check the properties the composed scheme is supposed to have.

Cases are built by ``make_case_fixtures.py`` from analytic parent fields whose
discrete divergence and net boundary flux are known exactly, written through
the production writer ``tools/python/udprep/nesting.py``.

| Test | What it establishes |
|---|---|
| I1 | `lnesting = .false.` leaves an existing case bit-identical to the pre-branch binary, built here from `origin/master` |
| I2 | the imposed face values survive `poisson` + `tstep_integrate` (design finding F2) |
| I3 | a uniform parent is preserved exactly, with `Phi`, `divmax`, `divtot` and `p` at round-off |
| I4 | a solenoidal parent needs no correction; `|grad p|` is second order in `h` |
| I5 | I3 and I4 give the same fields on 1x1, 2x1, 1x2 and 2x2; and a ramp + `tau > 0` + a cube on 2x2 (`TestI5CubeParity2x2`, in CI) |
| I6 | 100 steps == 50 + restart + 50, bitwise, on and off a parent interval boundary; and the same on 2x2 with a cube (`TestI6CubeRestartParity2x2`) |
| I7 | a difference confined to the interior stays out of the zone, and by how much |
| I8 | facet stresses (`tau_x`, `tau_y`, `tau_z`, `pres`) on a building at the zone edge, nested versus not |
| I9 | a cold start with `nest_linitfromparent` starts *at* the parent, divergence free |
| I10 | `BCtopm_pressure` (design case B) keeps the flux assertion on, and it is correct |

Environment (see ``_launch.py`` and ``_baseline.py`` for the full lists):
  UDALES_BUILD            path to the u-dales executable
                          (default build/release/u-dales)
  UDALES_BASELINE_REF     git ref I1 builds its pre-branch baseline from
                          (default origin/master, configured like UDALES_BUILD)
  UDALES_BASELINE         a ready-made pre-branch executable for I1, which
                          skips that build
  UDALES_RUNTIME_MODULES  module stack loaded before the run
  UDALES_MPIEXEC          MPI launcher (then MPIEXEC, then PATH)
  UDALES_REQUIRE_LAUNCHER =1: an unusable launcher fails instead of skipping
  TMPDIR                  where the run directories go
  UDALES_NESTING_KEEP     if set, run directories are kept, not deleted
"""

from __future__ import annotations

import math
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

TEST_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TEST_DIR))

import _launch as launch  # noqa: E402

REPO_ROOT = launch.REPO_ROOT
sys.path.insert(0, str(REPO_ROOT / "tools" / "python"))

import _baseline  # noqa: E402
import make_case_fixtures as mcf  # noqa: E402

EXPNR = mcf.EXPNR

UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))

#: Round-off yardsticks.  ``EPS`` is one double ulp; the field tolerances are
#: expressed as small multiples of it so that a failure means "not round-off",
#: not "slightly over an arbitrary number".
EPS = np.finfo(np.float64).eps
ROUNDOFF = 50.0 * EPS          # ~1.1e-14: a handful of ulps of an O(1) field
PARITY_TOL = 1.0e-9            # design 10.3 I5, matching processor_boundaries


# --------------------------------------------------------------------------- #
# Running the solver
# --------------------------------------------------------------------------- #


def run_solver(run_dir: Path, nprocs: int = 1, namelist: str = f"namoptions.{EXPNR}",
               executable: Optional[Path] = None) -> subprocess.CompletedProcess:
    return launch.run(run_dir, namelist, nprocs, executable or UDALES_BUILD,
                      quiet_intel_diagnostics=True)


def _tail(label: str, text: str, limit: int = 40) -> str:
    return launch.tail(label, text, limit)


def check_ok(case: unittest.TestCase, done: subprocess.CompletedProcess, what: str) -> str:
    if done.returncode != 0:
        case.fail(f"{what} exited {done.returncode}\n"
                  + _tail("stdout", done.stdout) + "\n" + _tail("stderr", done.stderr))
    return (done.stdout or "") + (done.stderr or "")


# --------------------------------------------------------------------------- #
# Reading what the solver reported
# --------------------------------------------------------------------------- #

_STATS_PATTERNS = {
    "t": re.compile(r"modnesting: t\s*=\s*(\S+)"),
    "phi": re.compile(r"modnesting: Phi \(norm\)\s*=\s*(\S+)"),
    "misfit": re.compile(r"modnesting: zone misfit rms \[m/s\]\s*=\s*(\S+)"),
}
_PHILID = re.compile(
    r"modnesting: Phi lid\s*=\s*(\S+)\s+closed faces\s*=\s*(\S+)"
)
_GRADP = re.compile(
    r"modnesting: \|grad p\| zone\s*=\s*(\S+)\s+interior\s*=\s*(\S+)\s+ratio\s*=\s*(\S+)"
)
_ENERGY = re.compile(
    r"modnesting: energy injected guard\s*=\s*(\S+)\s+relaxation\s*=\s*(\S+)"
)
_DIV = re.compile(r"divmax, divtot\s*=\s*(\S+)\s+(\S+)")


def parse_nesting_stats(output: str) -> List[Dict[str, float]]:
    """One record per ``nesting_stats`` block, in order of appearance.

    ``nesting_stats`` is called once per RK substep, so the last record of a run
    is the third substep of the last step.  Comparisons across runs must always
    use the same substep, which in practice means the last one.
    """
    records: List[Dict[str, float]] = []
    current: Dict[str, float] = {}
    for line in output.splitlines():
        for key, pattern in _STATS_PATTERNS.items():
            m = pattern.search(line)
            if m:
                if key == "t" and current:
                    records.append(current)
                    current = {}
                current[key] = float(m.group(1))
        m = _PHILID.search(line)
        if m:
            current["phi_lid"] = float(m.group(1))
            current["phi_closed"] = float(m.group(2))
        m = _GRADP.search(line)
        if m:
            current["gradp_zone"] = float(m.group(1))
            current["gradp_interior"] = float(m.group(2))
            current["gradp_ratio"] = float(m.group(3))
        m = _ENERGY.search(line)
        if m:
            current["e_guard"] = float(m.group(1))
            current["e_relax"] = float(m.group(2))
    if current:
        records.append(current)
    return records


def parse_divergence(output: str) -> List[Tuple[float, float]]:
    return [(float(m.group(1)), float(m.group(2))) for m in _DIV.finditer(output)]


# --------------------------------------------------------------------------- #
# Restart files: the only double-precision view of the solver's fields
# --------------------------------------------------------------------------- #
#
# All NetCDF output in uDALES is written as NF90_FLOAT (modstat_nc.f90), which
# cannot resolve a 1e-9 parity claim, let alone a round-off one.  The `initd`
# restart file is unformatted double precision, so every field assertion below
# reads that instead.  Record order is modsave.f90:88-100.

_REC_MINDIST, _REC_WALL = 0, 1
RESTART_FIELDS = {
    "u0": 2, "v0": 3, "w0": 4, "pres0": 5, "thl0": 6,
    "e120": 7, "ekm": 8, "qt0": 9, "ql0": 10, "ql0h": 11,
}
_REC_TIME = 12


def read_fortran_records(path: Path) -> List[bytes]:
    """Split an ifort sequential-unformatted file into its records."""
    data = path.read_bytes()
    records: List[bytes] = []
    off = 0
    while off < len(data):
        (n,) = struct.unpack("<i", data[off:off + 4])
        if n < 0:
            raise RuntimeError(f"{path}: subrecord markers are not supported")
        records.append(data[off + 4:off + 4 + n])
        off += 8 + n
    return records


def latest_restart(run_dir: Path, px: int, py: int, expnr: str = EXPNR) -> Path:
    files = sorted(run_dir.glob(f"initd*_{px:03d}_{py:03d}.{expnr}"))
    if not files:
        raise RuntimeError(f"no restart file for rank ({px},{py}) in {run_dir}")
    return files[-1]


def read_restart_fields(run_dir: Path, spec: mcf.CaseSpec, nprocx: int, nprocy: int,
                        names: Sequence[str] = ("u0", "v0", "w0", "pres0"),
                        expnr: str = EXPNR) -> Dict[str, np.ndarray]:
    """Stitch the per-rank restart files into global arrays.

    Returns, for each requested field, an array of shape
    ``(ktot + 1, jtot + 2, itot + 2)`` laid out like the solver's own
    ``(kb:ke+kh, jb-1:je+1, ib-1:ie+1)`` block: index ``[k, j+1, i+1]`` is the
    solver's ``(ib+i, jb+j, kb+k)``.  The two extra rows and columns therefore
    hold the *domain* ghost planes on the outer ranks -- which is where the
    east face ``u(ie+1)`` and the north face ``v(je+1)`` live, and I2 needs
    them.  Interior processor ghosts are simply overwritten by their owners.
    """
    imax, jmax = spec.itot // nprocx, spec.jtot // nprocy
    nx, ny, nz = imax + 2, jmax + 2, spec.ktot + 1
    out = {name: np.full((nz, spec.jtot + 2, spec.itot + 2), np.nan) for name in names}
    times: List[np.ndarray] = []
    for px in range(nprocx):
        for py in range(nprocy):
            records = read_fortran_records(latest_restart(run_dir, px, py, expnr))
            times.append(np.frombuffer(records[_REC_TIME], dtype="<f8"))
            for name in names:
                block = np.frombuffer(records[RESTART_FIELDS[name]],
                                      dtype="<f8").reshape((nz, ny, nx))
                i0, j0 = px * imax, py * jmax
                # interior first, then the outer ghost planes only where they
                # are domain boundaries rather than processor boundaries
                out[name][:, j0 + 1:j0 + 1 + jmax, i0 + 1:i0 + 1 + imax] = block[:, 1:-1, 1:-1]
                if px == 0:
                    out[name][:, j0 + 1:j0 + 1 + jmax, 0] = block[:, 1:-1, 0]
                if px == nprocx - 1:
                    out[name][:, j0 + 1:j0 + 1 + jmax, -1] = block[:, 1:-1, -1]
                if py == 0:
                    out[name][:, 0, i0 + 1:i0 + 1 + imax] = block[:, 0, 1:-1]
                if py == nprocy - 1:
                    out[name][:, -1, i0 + 1:i0 + 1 + imax] = block[:, -1, 1:-1]
    out["timee"] = np.array([t[0] for t in times])
    out["dt"] = np.array([t[1] for t in times])
    return out


def interior(field: np.ndarray) -> np.ndarray:
    """Drop the ghost ring, leaving ``(ktot+1, jtot, itot)``."""
    return field[:, 1:-1, 1:-1]


# --------------------------------------------------------------------------- #
# Base fixture
# --------------------------------------------------------------------------- #


class _NestingCase(unittest.TestCase):
    """Shared setup: a temporary root and the usual availability checks."""

    @classmethod
    def setUpClass(cls) -> None:
        if not UDALES_BUILD.is_file():
            raise RuntimeError(f"u-dales executable not found at {UDALES_BUILD}")
        launch.require_launcher()
        scratch = launch.scratch_dir()
        if os.environ.get("UDALES_NESTING_KEEP"):
            cls._tmp = None
            cls.root = Path(tempfile.mkdtemp(prefix="udales-nesting-int-", dir=scratch))
            print(f"\n[nesting] run directories kept in {cls.root}", flush=True)
        else:
            cls._tmp = tempfile.TemporaryDirectory(prefix="udales-nesting-int-", dir=scratch)
            cls.root = Path(cls._tmp.name)

    @classmethod
    def tearDownClass(cls) -> None:
        if getattr(cls, "_tmp", None) is not None:
            cls._tmp.cleanup()

    def make_run(self, name: str) -> Path:
        run_dir = self.root / name
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True)
        return run_dir


# --------------------------------------------------------------------------- #
# I1 -- no-op guarantee
# --------------------------------------------------------------------------- #
#
# Design section 10.3 asks for "bitwise identical to the pre-branch binary".
# That is achievable on a small case and NOT achievable on a large one, and the
# reason is not the nesting branch: uDALES plans all of its FFTs with
# FFTW_MEASURE (src/modpois.f90:110-191, 2decomp-fft/src/fft_fftw3.f90:26),
# which picks the transform algorithm by timing it at run time.  Two runs of
# *the same* binary on tests/cases/526 therefore differ at the 1e-16..1e-15
# relative level, with or without this branch.  So I1 is split:
#
#   TestI1NoOpSmallCase -- a case small enough that FFTW always picks the same
#       plan.  Self-reproducibility is verified first, then bitwise identity
#       against the baseline is REQUIRED.
#   TestI1NoOpExistingCase -- tests/cases/526 (IBM, trees, statistics, 4 ranks).
#       The baseline's own run-to-run spread is measured, and the branch is
#       required to be no further from the baseline than the baseline is from
#       itself.


class _I1Base(_NestingCase):
    """Shared plumbing: run the same directory twice, once per executable.

    The baseline is built here, by ``_baseline.ensure_baseline``, from
    ``UDALES_BASELINE_REF`` (default ``origin/master``) with the compiler,
    build type and library paths read back from the build under test -- so a
    Debug branch build is compared against a Debug baseline and a Release one
    against a Release one, and the test no longer depends on an uncommitted
    binary.  ``UDALES_BASELINE`` overrides it with a ready-made executable.
    """

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        try:
            cls.baseline = _baseline.ensure_baseline(UDALES_BUILD)
        except _baseline.BaselineError as exc:
            raise RuntimeError(f"I1 has no baseline to compare against: {exc}") from exc
        cls._runs: Dict[str, Tuple[Path, str]] = {}

    @staticmethod
    def _filtered_stdout(text: str) -> List[str]:
        """Drop everything that legitimately varies between two identical runs.

        gfortran's runtime diagnostics ("At line N of file /abs/path/x.f90")
        embed the source path, and the baseline is compiled from a different
        tree by construction, so the path is reduced to the file name.
        """
        drop = re.compile(r"Time of Day|CPU time|Elapsed|wall|WALL|unit =|^\s*$")
        path = re.compile(r"(At line \d+ of file )\S*/")
        return [path.sub(r"\1", ln.rstrip()) for ln in text.splitlines() if not drop.search(ln)]

    @staticmethod
    def _split_divergence(lines: List[str]) -> Tuple[List[str], List[Tuple[float, float]]]:
        kept, divs = [], []
        for ln in lines:
            m = _DIV.search(ln)
            if m:
                divs.append((float(m.group(1)), float(m.group(2))))
            else:
                kept.append(ln)
        return kept, divs

    @classmethod
    def _spread(cls, a: List[Tuple[float, float]], b: List[Tuple[float, float]]) -> float:
        worst = 0.0
        for (am, at), (bm, bt) in zip(a, b):
            for x, y in ((am, bm), (at, bt)):
                worst = max(worst, abs(x - y) / max(abs(x), abs(y), 1.0e-30))
        return worst

    @classmethod
    def _field_spread(cls, x: Path, y: Path) -> Dict[str, float]:
        """Largest relative difference per restart record."""
        ra, rb = read_fortran_records(x), read_fortran_records(y)
        out = {}
        for name, rec in RESTART_FIELDS.items():
            u = np.frombuffer(ra[rec], dtype="<f8")
            v = np.frombuffer(rb[rec], dtype="<f8")
            scale = max(float(np.max(np.abs(u))), 1.0e-30)
            out[name] = float(np.max(np.abs(u - v))) / scale
        return out


class TestI1NoOpSmallCase(_I1Base):
    """A cheap periodic case: the branch must be no further from the baseline
    than the baseline is from itself.

    Same grid, timestep and randomisation as the I6 control -- periodic
    laterals, `lnesting = .false.`.  This case was assumed small enough to be
    bitwise reproducible run to run and compared with a plain byte-equality
    assertion; in CI (gfortran, ubuntu-latest Release) two runs of the same
    baseline executable on the same input differed at the 1e-12 level in
    `pres0`, which rules out anything in this branch -- both those runs used
    the pre-branch `origin/master` binary.  The cause is the same one
    `TestI1NoOpExistingCase` below already documents: `ipoiss` here is the
    FFTW-based solver, and `FFTW_MEASURE` benchmarks planner variants against
    wall-clock time in `modpois.f90`, so even this small a transform is not
    guaranteed bit-identical on a noisy runner.  Follow the same "no worse
    than the baseline's own run-to-run spread" comparison used there instead
    of assuming exact reproducibility.
    """

    SPEC = mcf.ZONED
    DT = 0.125
    NSTEPS = 20

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.dirs: Dict[str, Path] = {}
        cls.out: Dict[str, str] = {}
        cls.error: Optional[str] = None
        for label, exe in (("base_a", cls.baseline), ("base_b", cls.baseline),
                           ("head", UDALES_BUILD)):
            run_dir = cls.root / f"i1_small_{label}"
            mcf.write_case(
                run_dir, cls.SPEC, "uniform", times=(0.0, 1.0e6), uprof=1.0, U=1.0,
                edits={"dtmax": f"{cls.DT:.10g}",
                       "runtime": f"{cls.DT * cls.NSTEPS:.10g}",
                       "trestart": f"{cls.DT * cls.NSTEPS:.10g}",
                       "randu": "0.05", "krand": str(cls.SPEC.ktot),
                       "lrandomize": ".true.",
                       "BCxm": "1", "BCym": "1", "lnesting": ".false."})
            done = run_solver(run_dir, executable=exe)
            if done.returncode != 0:
                cls.error = (f"I1 small case with {exe.name} exited {done.returncode}\n"
                             + _tail("stdout", done.stdout) + "\n"
                             + _tail("stderr", done.stderr))
                return
            cls.dirs[label] = run_dir
            cls.out[label] = (done.stdout or "") + (done.stderr or "")

    def setUp(self) -> None:
        if self.error:
            self.fail(self.error)

    def test_baseline_self_noise(self) -> None:
        a = latest_restart(self.dirs["base_a"], 0, 0)
        b = latest_restart(self.dirs["base_b"], 0, 0)
        bitwise = a.read_bytes() == b.read_bytes()
        spread = self._field_spread(a, b)
        worst = max(spread.values())
        self.__class__.self_noise = worst
        print(f"\n[I1] small case baseline vs baseline: bitwise {bitwise}, "
              f"worst relative field difference {worst:.3e}", flush=True)
        for name, value in sorted(spread.items(), key=lambda kv: -kv[1])[:4]:
            print(f"[I1]   {name}: {value:.3e}", flush=True)
        self.assertLess(worst, 1.0e-9,
                        "the baseline's own run-to-run spread is larger than round-off")

    def test_branch_matches_baseline(self) -> None:
        base = latest_restart(self.dirs["base_a"], 0, 0)
        head = latest_restart(self.dirs["head"], 0, 0)
        reference = getattr(self.__class__, "self_noise", None)
        if reference is None:
            self.test_baseline_self_noise()
            reference = self.self_noise
        spread = self._field_spread(base, head)
        for name, value in spread.items():
            print(f"[I1] small case {name}: relative diff {value:.3e}", flush=True)
        worst = max(spread.values())
        self.assertLessEqual(
            worst, max(10.0 * reference, 1.0e-13),
            "with lnesting = .false. the branch moved the small case further than "
            "the baseline moves against itself")

    def test_stdout_matches_baseline(self) -> None:
        base_lines, base_div = self._split_divergence(self._filtered_stdout(self.out["base_a"]))
        head_lines, head_div = self._split_divergence(self._filtered_stdout(self.out["head"]))
        self_lines, self_div = self._split_divergence(self._filtered_stdout(self.out["base_b"]))
        self.assertEqual(len(base_lines), len(head_lines),
                         "the branch prints a different number of lines")
        for n, (a, b) in enumerate(zip(base_lines, head_lines)):
            self.assertEqual(a, b, f"stdout line {n} differs:\n  base: {a}\n  head: {b}")
        self.assertEqual(len(base_div), len(head_div))
        reference = max(self._spread(base_div, self_div), 1.0e-3)
        against = self._spread(base_div, head_div)
        print(f"[I1] small case divmax/divtot: baseline-vs-branch spread {against:.3e}, "
              f"baseline-vs-baseline {self._spread(base_div, self_div):.3e}", flush=True)
        self.assertLessEqual(
            against, max(3.0 * reference, 2.0),
            "divmax/divtot moved further than the baseline's own run-to-run spread")


class TestI1NoOpExistingCase(_I1Base):
    """`tests/cases/526` unchanged: IBM, trees, statistics, on 4 ranks.

    This case is NOT bitwise reproducible run to run -- FFTW_MEASURE, see the
    note above -- so the assertion is the strongest one that is well posed:
    the branch must be no further from the baseline than the baseline is from
    itself.  `test_baseline_self_noise` measures that reference spread and
    prints it, and it also fails if the case turns out to be reproducible after
    all, at which point this test should demand bitwise identity instead.

    `ladaptive` is forced off so that a round-off difference cannot feed back
    into a different timestep and turn a 1e-16 difference into a 1e-3 one.
    """

    CASE_ID = 526
    NPROCS = 4

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.case_source = REPO_ROOT / "tests" / "cases" / str(cls.CASE_ID)
        if not cls.case_source.is_dir():
            raise unittest.SkipTest(f"case {cls.CASE_ID} not found")
        cls.dirs: Dict[str, Path] = {}
        cls.out: Dict[str, str] = {}
        cls.error: Optional[str] = None
        for label, exe in (("base_a", cls.baseline), ("base_b", cls.baseline),
                           ("head", UDALES_BUILD)):
            run_dir = cls.root / f"i1_case_{label}"
            shutil.copytree(cls.case_source, run_dir)
            nml = run_dir / f"namoptions.{cls.CASE_ID}"
            text = nml.read_text(encoding="utf-8")
            for key, value in (("runtime", "0.02"), ("trestart", "0.01"),
                               ("ladaptive", ".false."), ("dtmax", "0.01")):
                text = re.sub(rf"(?m)^(\s*{key}\s*=\s*).*$", rf"\g<1>{value}", text)
            nml.write_text(text, encoding="utf-8")
            done = run_solver(run_dir, nprocs=cls.NPROCS,
                              namelist=f"namoptions.{cls.CASE_ID}", executable=exe)
            if done.returncode != 0:
                cls.error = (f"case {cls.CASE_ID} with {exe.name} exited "
                             f"{done.returncode}\n" + _tail("stdout", done.stdout)
                             + "\n" + _tail("stderr", done.stderr))
                return
            cls.dirs[label] = run_dir
            cls.out[label] = (done.stdout or "") + (done.stderr or "")

    def setUp(self) -> None:
        if self.error:
            self.fail(self.error)

    def _restarts(self, label: str) -> List[Path]:
        return sorted(self.dirs[label].glob(f"initd*.{self.CASE_ID}"))

    def test_baseline_self_noise(self) -> None:
        a, b = self._restarts("base_a"), self._restarts("base_b")
        self.assertTrue(a, "the baseline run wrote no restart file")
        bitwise = all(x.read_bytes() == y.read_bytes() for x, y in zip(a, b))
        spread = self._field_spread(a[-1], b[-1])
        worst = max(spread.values())
        self.__class__.self_noise = worst
        print(f"\n[I1] case {self.CASE_ID} baseline vs baseline: "
              f"bitwise {bitwise}, worst relative field difference {worst:.3e}", flush=True)
        for name, value in sorted(spread.items(), key=lambda kv: -kv[1])[:4]:
            print(f"[I1]   {name}: {value:.3e}", flush=True)
        _, da = self._split_divergence(self._filtered_stdout(self.out["base_a"]))
        _, db = self._split_divergence(self._filtered_stdout(self.out["base_b"]))
        print(f"[I1] baseline vs baseline divmax/divtot relative spread: "
              f"{self._spread(da, db):.3e}", flush=True)
        self.assertFalse(
            bitwise,
            f"case {self.CASE_ID} turned out to be bitwise reproducible; the comparison "
            "against the branch should be tightened from 'no worse than self-noise' to "
            "'bitwise'")
        # Loose on purpose: `pres0` is an accumulated pressure whose values sit
        # near zero, so normalising by its own peak inflates the ratio.  The
        # point of the bound is only that the spread is round-off and not a real
        # difference; the branch is judged against the measured value, not this.
        self.assertLess(worst, 1.0e-9,
                        "the baseline's own run-to-run spread is larger than round-off")

    def test_branch_is_within_the_baseline_self_noise(self) -> None:
        base, head = self._restarts("base_a"), self._restarts("head")
        self.assertEqual([f.name for f in base], [f.name for f in head],
                         "the branch wrote a different set of restart files")
        reference = getattr(self.__class__, "self_noise", None)
        if reference is None:
            self.test_baseline_self_noise()
            reference = self.self_noise
        spread = self._field_spread(base[-1], head[-1])
        worst = max(spread.values())
        print(f"[I1] case {self.CASE_ID} baseline vs branch: worst relative field "
              f"difference {worst:.3e} (baseline self-noise {reference:.3e})", flush=True)
        for name, value in sorted(spread.items(), key=lambda kv: -kv[1])[:4]:
            print(f"[I1]   {name}: {value:.3e}", flush=True)
        self.assertLessEqual(
            worst, max(10.0 * reference, 1.0e-13),
            f"with lnesting = .false. the branch moved case {self.CASE_ID} further than "
            "the baseline moves against itself")

    def test_stdout_is_identical_apart_from_divergence(self) -> None:
        base_lines, base_div = self._split_divergence(self._filtered_stdout(self.out["base_a"]))
        head_lines, head_div = self._split_divergence(self._filtered_stdout(self.out["head"]))
        self_lines, self_div = self._split_divergence(self._filtered_stdout(self.out["base_b"]))
        self.assertEqual(len(base_lines), len(head_lines),
                         "the branch prints a different number of lines")
        for n, (a, b) in enumerate(zip(base_lines, head_lines)):
            self.assertEqual(a, b, f"stdout line {n} differs:\n  base: {a}\n  head: {b}")
        self.assertEqual(len(base_div), len(head_div))
        reference = max(self._spread(base_div, self_div), 1.0e-3)
        against = self._spread(base_div, head_div)
        print(f"[I1] divmax/divtot: baseline-vs-branch spread {against:.3e}, "
              f"baseline-vs-baseline {self._spread(base_div, self_div):.3e}", flush=True)
        self.assertLessEqual(
            against, max(3.0 * reference, 2.0),
            "divmax/divtot moved further than the baseline's own run-to-run spread")



# --------------------------------------------------------------------------- #
# I2 -- face imposition survives projection (design finding F2)
# --------------------------------------------------------------------------- #


class TestI2FaceImposition(_NestingCase):
    """The four lateral faces must still carry the parent value after `poisson`.

    The parent is `face_forced`: identical on the two `x` boundary faces (so
    `Phi = 0` exactly, and the flux assertion is not what is being tested) but
    strongly divergent inside, so `poisson` has real work to do.  If the
    pressure gradient leaked into the boundary-normal velocity the face values
    would move, which is exactly what design finding F2 says must not happen.
    """

    SPEC = mcf.BASE

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.run_dir = cls.root / "i2"
        cls.data = mcf.write_case(
            cls.run_dir, cls.SPEC, "face_forced", times=(0.0, 1.0e6), uprof=1.0,
            edits={"runtime": "0.5", "dtmax": "0.5", "trestart": "0.5"},
        )
        done = run_solver(cls.run_dir)
        cls.output = (done.stdout or "") + (done.stderr or "")
        cls.returncode = done.returncode

    def setUp(self) -> None:
        if self.returncode != 0:
            self.fail("I2 run failed\n" + _tail("output", self.output))

    def test_faces_match_the_parent_to_roundoff(self) -> None:
        spec = self.SPEC
        g = spec.grid()
        u_ref, v_ref, _ = mcf.FIELDS["face_forced"](g, 0.0)
        fields = read_restart_fields(self.run_dir, spec, 1, 1)
        # solver index [k, j+1, i+1] <-> (ib+i, jb+j, kb+k); u(ib) is i = 0 and
        # u(ie+1) is the outer ghost column, i.e. index -1.
        u = fields["u0"]
        v = fields["v0"]
        checks = {
            "u west (i = ib)":   (u[:spec.ktot, 1:-1, 1].T,  u_ref[0, :, :]),
            "u east (i = ie+1)": (u[:spec.ktot, 1:-1, -1].T, u_ref[spec.itot, :, :]),
            "v south (j = jb)":  (v[:spec.ktot, 1, 1:-1].T,  v_ref[:, 0, :]),
            "v north (j = je+1)": (v[:spec.ktot, -1, 1:-1].T, v_ref[:, spec.jtot, :]),
        }
        failures = []
        for label, (got, want) in checks.items():
            err = float(np.max(np.abs(got - want)))
            print(f"[I2] {label}: max|u - u_parent| = {err:.3e}", flush=True)
            if err > ROUNDOFF:
                failures.append(f"{label}: max abs error {err:.3e} > {ROUNDOFF:.1e}")
        if failures:
            self.fail("imposed face values did not survive the projection:\n- "
                      + "\n- ".join(failures))

    def test_the_projection_actually_did_something(self) -> None:
        """Guard against a vacuous pass: the field must be genuinely divergent."""
        g = self.SPEC.grid()
        div = mcf.discrete_divergence(g, mcf.FIELDS["face_forced"](g, 0.0))
        self.assertGreater(float(np.abs(div).max()), 1.0e-3,
                           "the I2 parent field is not divergent, so I2 proves nothing")
        stats = parse_nesting_stats(self.output)
        self.assertTrue(stats, "no nesting_stats output")
        self.assertGreater(stats[-1]["gradp_zone"], 1.0e-6,
                           "|grad p| is at round-off, so the projection did no work")


# --------------------------------------------------------------------------- #
# I3 -- uniform flow
# --------------------------------------------------------------------------- #


class TestI3UniformFlow(_NestingCase):
    """A constant `(U, 0, 0)` parent must be preserved exactly.

    The most informative single test in the matrix: a uniform field is a fixed
    point of advection, of the subgrid model on a free-slip box, and of the
    projection, and it carries `Phi = 0` identically.  Anything that moves is
    the nesting scheme moving it.
    """

    SPEC = mcf.BASE
    U = 1.0

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.run_dir = cls.root / "i3"
        mcf.write_case(cls.run_dir, cls.SPEC, "uniform", times=(0.0, 1.0e6),
                       uprof=cls.U, U=cls.U,
                       edits={"runtime": "5.", "dtmax": "0.5", "trestart": "5."})
        done = run_solver(cls.run_dir)
        cls.output = (done.stdout or "") + (done.stderr or "")
        cls.returncode = done.returncode

    def setUp(self) -> None:
        if self.returncode != 0:
            self.fail("I3 run failed\n" + _tail("output", self.output))

    def test_field_is_preserved(self) -> None:
        f = read_restart_fields(self.run_dir, self.SPEC, 1, 1)
        du = float(np.max(np.abs(interior(f["u0"]) - self.U)))
        dv = float(np.max(np.abs(interior(f["v0"]))))
        dw = float(np.max(np.abs(interior(f["w0"]))))
        dp = float(np.max(np.abs(interior(f["pres0"]))))
        print(f"\n[I3] max|u-U| = {du:.3e}  max|v| = {dv:.3e}  max|w| = {dw:.3e}  "
              f"max|pres0| = {dp:.3e}  ({du / EPS:.1f} ulp)", flush=True)
        self.assertLessEqual(du, ROUNDOFF, "u drifted from the imposed uniform value")
        self.assertLessEqual(dv, ROUNDOFF, "v was generated from nothing")
        self.assertLessEqual(dw, ROUNDOFF, "w was generated from nothing")
        self.assertLessEqual(dp, ROUNDOFF, "the accumulated pressure is not zero")

    def test_diagnostics_are_at_roundoff(self) -> None:
        stats = parse_nesting_stats(self.output)
        divs = parse_divergence(self.output)
        self.assertTrue(stats and divs, "no diagnostics in the output")
        worst_phi = max(abs(s["phi"]) for s in stats)
        worst_mis = max(s["misfit"] for s in stats)
        worst_gp = max(s["gradp_zone"] for s in stats)
        worst_divmax = max(abs(a) for a, _ in divs)
        worst_divtot = max(abs(b) for _, b in divs)
        print(f"[I3] max|Phi| = {worst_phi:.3e}  max misfit = {worst_mis:.3e}  "
              f"max |grad p|_zone = {worst_gp:.3e}  "
              f"max divmax = {worst_divmax:.3e}  max divtot = {worst_divtot:.3e}",
              flush=True)
        self.assertLessEqual(worst_phi, 1.0e-12, "the net boundary flux is not zero")
        self.assertLessEqual(worst_mis, ROUNDOFF, "the zone does not match the parent")
        self.assertLessEqual(worst_gp, ROUNDOFF, "the projection is correcting something")
        self.assertLessEqual(worst_divmax, ROUNDOFF, "divmax is not at round-off")
        self.assertLessEqual(worst_divtot, ROUNDOFF, "divtot is not at round-off")


# --------------------------------------------------------------------------- #
# I4 -- manufactured solenoidal field
# --------------------------------------------------------------------------- #


class TestI4ManufacturedSolenoidal(_NestingCase):
    """A solenoidal parent must leave the projection nothing to do.

    Two halves, because the field design section 10.3 names turns out to be a
    *stronger* case than the design expected:

    * `test_taylor_green_needs_no_correction` runs exactly that field.  On a
      staggered grid its discrete divergence cancels identically (see
      `test_taylor_green_divergence_is_exact`), so `|grad p|` sits at round-off
      on every grid and there is no convergence rate to measure -- the
      second-order claim is not observable for this field, and asserting it
      would be asserting a property of the round-off.
    * `test_mixed_mode_gradp_is_second_order` therefore runs a second
      manufactured field with the same construction (a stream function, so
      analytically solenoidal) but different wavenumbers in `x` and `y`, whose
      discrete divergence really is `O(h^2)`.  That is the field on which the
      second-order convergence of `||grad p||` can be, and is, measured.
    """

    SPECS = mcf.CONVERGENCE
    RUN_EDITS = {"runtime": "1.", "dtmax": "0.25", "trestart": "1."}

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.runs: Dict[Tuple[str, int], Tuple[Path, str, int]] = {}
        for field in ("taylor_green", "mixed_mode"):
            for spec in cls.SPECS:
                run_dir = cls.root / f"i4_{field}_{spec.itot}"
                mcf.write_case(run_dir, spec, field, times=(0.0, 1.0e6),
                               edits=dict(cls.RUN_EDITS))
                done = run_solver(run_dir)
                cls.runs[(field, spec.itot)] = (
                    run_dir, (done.stdout or "") + (done.stderr or ""), done.returncode
                )

    def _stats(self, field: str, spec: mcf.CaseSpec) -> Dict[str, float]:
        run_dir, output, rc = self.runs[(field, spec.itot)]
        if rc != 0:
            self.fail(f"I4 {field} on {spec.itot}^3 failed\n" + _tail("output", output))
        stats = parse_nesting_stats(output)
        self.assertTrue(stats, f"no nesting_stats output for {field} {spec.itot}")
        return stats[-1]

    def test_taylor_green_divergence_is_exact(self) -> None:
        """Why the second-order claim is not measurable for the design's field."""
        for spec in self.SPECS:
            g = spec.grid()
            div = mcf.discrete_divergence(g, mcf.FIELDS["taylor_green"](g, 0.0))
            peak = float(np.abs(div).max())
            print(f"[I4] taylor_green h = {spec.dx:g}: max|div| = {peak:.3e}", flush=True)
            self.assertLessEqual(peak, 1.0e-13,
                                 "the Taylor-Green field is no longer discretely solenoidal")

    def test_taylor_green_needs_no_correction(self) -> None:
        for spec in self.SPECS:
            s = self._stats("taylor_green", spec)
            divs = parse_divergence(self.runs[("taylor_green", spec.itot)][1])
            worst_div = max(abs(a) for a, _ in divs) if divs else 0.0
            print(f"[I4] taylor_green h = {spec.dx:g}: |grad p|_zone = "
                  f"{s['gradp_zone']:.3e}  misfit = {s['misfit']:.3e}  "
                  f"divmax = {worst_div:.3e}  Phi = {s['phi']:.3e}", flush=True)
            self.assertLessEqual(s["gradp_zone"], 1.0e-13,
                                 f"h = {spec.dx}: the projection corrected a solenoidal field")
            self.assertLessEqual(worst_div, 1.0e-13, f"h = {spec.dx}: divmax is not at round-off")
            self.assertLessEqual(abs(s["phi"]), 1.0e-12, f"h = {spec.dx}: Phi is not zero")

    def test_mixed_mode_divergence_is_second_order(self) -> None:
        peaks = []
        for spec in self.SPECS:
            g = spec.grid()
            div = mcf.discrete_divergence(g, mcf.FIELDS["mixed_mode"](g, 0.0))
            peaks.append(float(np.abs(div).max()))
        orders = [math.log2(peaks[n] / peaks[n + 1]) for n in range(len(peaks) - 1)]
        print(f"[I4] mixed_mode max|div| = {['%.3e' % p for p in peaks]}  "
              f"orders = {['%.2f' % o for o in orders]}", flush=True)
        for o in orders:
            self.assertGreater(o, 1.8, "the mixed-mode field is not second order in h")

    def test_mixed_mode_gradp_is_second_order(self) -> None:
        gradp = [self._stats("mixed_mode", spec)["gradp_zone"] for spec in self.SPECS]
        orders = [math.log2(gradp[n] / gradp[n + 1]) for n in range(len(gradp) - 1)]
        print(f"[I4] mixed_mode |grad p|_zone = {['%.3e' % g for g in gradp]}  "
              f"orders = {['%.2f' % o for o in orders]}", flush=True)
        for o in orders:
            self.assertGreater(o, 1.8,
                               "||grad p|| does not converge at second order")
            self.assertLess(o, 2.4, "||grad p|| converges faster than second order")


# --------------------------------------------------------------------------- #
# I5 -- decomposition parity
# --------------------------------------------------------------------------- #


CONFIGS: Dict[str, Tuple[int, int]] = {
    "serial": (1, 1), "x_split": (2, 1), "y_split": (1, 2), "xy_split": (2, 2),
}


class TestI5DecompositionParity(_NestingCase):
    """I3 and I4 must give the same fields on 1x1, 2x1, 1x2 and 2x2.

    Mirrors `tests/integration/processor_boundaries/`, but compares the
    double-precision restart state rather than the NF90_FLOAT statistics
    output, so the 1e-9 tolerance of design section 10.3 is meaningful.
    """

    SPEC = mcf.BASE
    FIELDS = ("uniform", "taylor_green")

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.fields: Dict[Tuple[str, str], Dict[str, np.ndarray]] = {}
        cls.failures: List[str] = []
        for field in cls.FIELDS:
            for label, (nprocx, nprocy) in CONFIGS.items():
                run_dir = cls.root / f"i5_{field}_{label}"
                mcf.write_case(
                    run_dir, cls.SPEC, field, times=(0.0, 1.0e6),
                    uprof=1.0 if field == "uniform" else 0.0,
                    edits={"runtime": "2.", "dtmax": "0.5", "trestart": "2.",
                           "nprocx": str(nprocx), "nprocy": str(nprocy)},
                )
                done = run_solver(run_dir, nprocs=nprocx * nprocy)
                if done.returncode != 0:
                    cls.failures.append(
                        f"{field} on {label} exited {done.returncode}\n"
                        + _tail("stdout", done.stdout) + "\n" + _tail("stderr", done.stderr))
                    continue
                cls.fields[(field, label)] = read_restart_fields(
                    run_dir, cls.SPEC, nprocx, nprocy)

    def setUp(self) -> None:
        if self.failures:
            self.fail("\n\n".join(self.failures))

    def _compare(self, field: str) -> None:
        reference = self.fields[(field, "serial")]
        problems = []
        for label in CONFIGS:
            if label == "serial":
                continue
            candidate = self.fields[(field, label)]
            for name in ("u0", "v0", "w0", "pres0"):
                diff = np.abs(interior(candidate[name]) - interior(reference[name]))
                worst = float(np.nanmax(diff))
                idx = np.unravel_index(int(np.nanargmax(diff)), diff.shape)
                print(f"[I5] {field} {label} {name}: max abs diff {worst:.3e} at {idx}",
                      flush=True)
                if worst > PARITY_TOL:
                    problems.append(f"{field} {label} {name}: {worst:.3e} at {idx}")
        if problems:
            self.fail("decomposition parity failed:\n- " + "\n- ".join(problems))

    def test_uniform_flow_parity(self) -> None:
        self._compare("uniform")

    def test_taylor_green_parity(self) -> None:
        self._compare("taylor_green")


# --------------------------------------------------------------------------- #
# I6 -- restart parity
# --------------------------------------------------------------------------- #


class TestI6RestartParity(_NestingCase):
    """`N` steps must equal `N/2` + restart + `N/2`, bitwise.

    Run with a genuinely unsteady parent (`unsteady_uniform`) and a real
    guard + ramp zone, so the time buffer, the ramp and the interior dynamics
    all matter.  Two restart points are tested: one in the middle of a parent
    interval and one exactly on an interval boundary, which is the case design
    section 9.5 singles out.

    Every time here is an exact binary fraction (`dt = 0.125 s`, parent levels
    every `1.25 s`), so `timeleft`, `timee` and `tnextrestart` accumulate
    without rounding and both legs take exactly the intended number of steps.

    `test_control_restart_without_nesting` runs the *same* case with periodic
    laterals and `lnesting = .false.`.  Without it a failure here could not be
    attributed: uDALES's restart machinery has to be bitwise on its own before
    the nesting state can be blamed for anything.

    HISTORY: these tests found a real defect in `src/program.f90` --
    `nesting_init` ran before `readinitfiles`, which is what assigns `timee`,
    so on a warm start the parent buffer and the target were positioned at
    `t = 0` instead of at the restart time, and every nested case trapped on a
    Debug build where the undefined `timee` is a signalling NaN.  The call
    order was fixed; these tests are now bitwise green.  See "One defect, two
    symptoms" in README.md.
    """

    SPEC = mcf.ZONED
    DT = 0.125
    NSTEPS = 100
    PARENT_DT = 1.25
    #: step 50 -> t = 6.25 s, which is parent level 5 exactly;
    #: step 37 -> t = 4.625 s, strictly inside the interval [3.75, 5.0).
    SPLIT_ON_BOUNDARY = 50
    SPLIT_MID_INTERVAL = 37
    U0, AMP, PERIOD = 1.0, 0.25, 8.0

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.times = np.arange(0.0, cls.DT * cls.NSTEPS + 2 * cls.PARENT_DT, cls.PARENT_DT)
        cls._cache: Dict[str, Path] = {}

    # -- case construction ------------------------------------------------ #

    def _write(self, run_dir: Path, nested: bool, extra: Dict[str, str]) -> None:
        edits = {
            "dtmax": f"{self.DT:.10g}",
            "randu": "0.05",
            "krand": str(self.SPEC.ktot),
            "lrandomize": ".true.",
        }
        if not nested:
            edits.update({"BCxm": "1", "BCym": "1", "lnesting": ".false."})
        edits.update(extra)
        mcf.write_case(run_dir, self.SPEC, "unsteady_uniform", times=self.times,
                       uprof=self.U0, U0=self.U0, A=self.AMP, T=self.PERIOD,
                       edits=edits)

    def _leg(self, name: str, nested: bool, nsteps: int,
             startfile: Optional[str] = None) -> Tuple[Path, int, float]:
        """Run one leg and return its directory, its final `ntrun` and `timee`."""
        run_dir = self.root / name
        if run_dir.exists():
            shutil.rmtree(run_dir)
        extra = {"runtime": f"{self.DT * nsteps:.10g}",
                 "trestart": f"{self.DT * nsteps:.10g}"}
        if startfile is not None:
            extra.update({"lwarmstart": ".true.", "startfile": f"'{startfile}'"})
        self._write(run_dir, nested, extra)
        if startfile is not None:
            source = self.root / name.replace("_second", "_first")
            for f in list(source.glob("initd*")) + list(source.glob("inits*")):
                shutil.copy2(f, run_dir / f.name)
        check_ok(self, run_solver(run_dir), f"I6 leg {name}")
        latest = latest_restart(run_dir, 0, 0)
        ntrun = int(latest.name[5:13])
        timee = float(np.frombuffer(read_fortran_records(latest)[_REC_TIME], dtype="<f8")[0])
        return run_dir, ntrun, timee

    def _continuous(self, nested: bool) -> Path:
        key = f"continuous_{'nested' if nested else 'plain'}"
        if key not in self._cache:
            run_dir, ntrun, timee = self._leg(f"i6_{key}", nested, self.NSTEPS)
            self.assertEqual(ntrun, self.NSTEPS,
                             "the reference run did not take the intended number of steps")
            self.assertAlmostEqual(timee, self.DT * self.NSTEPS, places=12)
            self._cache[key] = run_dir
        return self._cache[key]

    def _split(self, label: str, nested: bool, nsplit: int) -> Path:
        first, ntrun, timee = self._leg(f"i6_{label}_first", nested, nsplit)
        self.assertEqual(ntrun, nsplit, f"{label}: the first leg stopped at step {ntrun}")
        self.assertAlmostEqual(timee, self.DT * nsplit, places=12)
        second, ntrun2, timee2 = self._leg(
            f"i6_{label}_second", nested, self.NSTEPS - nsplit,
            startfile=f"initd{nsplit:08d}_000_000.{EXPNR}")
        self.assertEqual(ntrun2, self.NSTEPS,
                         f"{label}: the second leg stopped at step {ntrun2}")
        self.assertAlmostEqual(timee2, self.DT * self.NSTEPS, places=12)
        return second

    # -- comparison ------------------------------------------------------- #

    def _compare(self, label: str, cont: Path, split: Path) -> List[str]:
        a = read_fortran_records(latest_restart(cont, 0, 0))
        b = read_fortran_records(latest_restart(split, 0, 0))
        problems = []
        for name, rec in RESTART_FIELDS.items():
            x = np.frombuffer(a[rec], dtype="<f8")
            y = np.frombuffer(b[rec], dtype="<f8")
            worst = float(np.max(np.abs(x - y)))
            same = a[rec] == b[rec]
            print(f"[I6] {label} {name}: bitwise {'yes' if same else 'NO':3s}  "
                  f"max abs diff {worst:.3e}", flush=True)
            if not same:
                problems.append(f"{name}: max abs diff {worst:.3e}")
        return problems

    def _penetration(self, label: str, cont: Path, split: Path) -> None:
        """Where the difference sits, so a failure says *why* as well as *that*."""
        fa = read_restart_fields(cont, self.SPEC, 1, 1)
        fb = read_restart_fields(split, self.SPEC, 1, 1)
        du = np.abs(interior(fa["u0"]) - interior(fb["u0"]))
        n = self.SPEC.itot
        idx = np.arange(n)
        di = np.minimum(idx, n - 1 - idx).astype(float)
        dj = np.minimum(np.arange(self.SPEC.jtot),
                        self.SPEC.jtot - 1 - np.arange(self.SPEC.jtot)).astype(float)
        dist = np.minimum(np.broadcast_to(di[None, None, :], du.shape),
                          np.broadcast_to(dj[None, :, None], du.shape))
        n_zone = int(round((self.SPEC.guardwidth + self.SPEC.zonewidth) / self.SPEC.dx))
        by_d = [(d, float(du[dist == d].max())) for d in range(int(dist.max()) + 1)]
        print(f"[I6] {label} max|du| by inward distance: "
              + "  ".join(f"{d}:{v:.2e}" for d, v in by_d[:n_zone + 3]), flush=True)

    # -- the tests -------------------------------------------------------- #

    def test_control_restart_without_nesting(self) -> None:
        """The same case, periodic and unnested, must restart bitwise."""
        cont = self._continuous(nested=False)
        split = self._split("control", nested=False, nsplit=self.SPLIT_ON_BOUNDARY)
        problems = self._compare("control (no nesting)", cont, split)
        if problems:
            self.fail("uDALES does not restart bitwise even WITHOUT nesting, so I6 "
                      "cannot attribute anything to the nesting scheme:\n- "
                      + "\n- ".join(problems))

    def test_nesting_init_positions_the_buffer_at_the_restart_time(self) -> None:
        """Root-cause isolator for the two bitwise tests below.

        `nesting_init` must run *after* `readinitfiles`, which is what reads
        `timee` back from the restart file; otherwise `set_interval(timee)` /
        `eval_target(timee)` position the parent buffer and the interpolated
        target at `t = 0` rather than at the restart time.  The line
        `nesting_init` prints says which interval it chose, so a regression in
        the call order is visible directly rather than only through the two
        bitwise tests below.
        """
        nsplit = self.SPLIT_MID_INTERVAL
        self._continuous(nested=True)
        self._split("probe", nested=True, nsplit=nsplit)
        run_dir = self.root / f"i6_probe_second"
        done = run_solver(run_dir)
        output = check_ok(self, done, "I6 buffer-position probe")
        m = re.search(r"buffer at interval\s+(\d+)", output)
        self.assertIsNotNone(m, "nesting_init did not report its buffer interval")
        got = int(m.group(1))
        t_restart = self.DT * nsplit
        expected = int(np.searchsorted(self.times, t_restart, side="right"))
        print(f"\n[I6] warm start at t = {t_restart:g} s (parent levels "
              f"{list(self.times)}): nesting_init reports interval {got}, "
              f"expected {expected}", flush=True)
        self.assertEqual(
            got, expected,
            "nesting_init positioned the parent buffer at the wrong time on a warm "
            "start; see this test's docstring for the mechanism")

    def test_restart_with_a_time_constant_parent(self) -> None:
        """The same restart with a parent that does not vary in time.

        Isolates the mechanism: if the target is the same at `t = 0` as at the
        restart time, mispositioning the buffer cannot matter, and the restart
        should come back at round-off.  Whatever is left over is a *second*,
        much smaller effect and is reported rather than swept up.
        """
        nsplit = self.SPLIT_MID_INTERVAL
        runs = {}
        for name, nsteps, startfile in (("const_cont", self.NSTEPS, None),
                                        ("const_first", nsplit, None),
                                        ("const_second", self.NSTEPS - nsplit,
                                         f"initd{nsplit:08d}_000_000.{EXPNR}")):
            run_dir = self.root / f"i6_{name}"
            if run_dir.exists():
                shutil.rmtree(run_dir)
            edits = {"dtmax": f"{self.DT:.10g}", "randu": "0.05",
                     "krand": str(self.SPEC.ktot), "lrandomize": ".true.",
                     "runtime": f"{self.DT * nsteps:.10g}",
                     "trestart": f"{self.DT * nsteps:.10g}"}
            if startfile:
                edits.update({"lwarmstart": ".true.", "startfile": f"'{startfile}'"})
            mcf.write_case(run_dir, self.SPEC, "uniform", times=self.times,
                           uprof=self.U0, U=self.U0, edits=edits)
            if startfile:
                for f in (self.root / "i6_const_first").glob("initd*"):
                    shutil.copy2(f, run_dir / f.name)
            check_ok(self, run_solver(run_dir), f"I6 {name}")
            runs[name] = run_dir
        a = read_fortran_records(latest_restart(runs["const_cont"], 0, 0))
        b = read_fortran_records(latest_restart(runs["const_second"], 0, 0))
        worst = 0.0
        for name, rec in RESTART_FIELDS.items():
            x = np.frombuffer(a[rec], dtype="<f8")
            y = np.frombuffer(b[rec], dtype="<f8")
            worst = max(worst, float(np.max(np.abs(x - y))))
        print(f"[I6] time-constant parent: max abs restart difference {worst:.3e} "
              f"({worst / EPS:.1f} ulp)", flush=True)
        self.assertLessEqual(
            worst, ROUNDOFF,
            "even with a time-constant parent the nested restart differs by more than "
            "round-off, so there is a second effect beyond the buffer position")

    def test_restart_mid_parent_interval(self) -> None:
        cont = self._continuous(nested=True)
        split = self._split("mid", nested=True, nsplit=self.SPLIT_MID_INTERVAL)
        problems = self._compare("mid-interval", cont, split)
        self._penetration("mid-interval", cont, split)
        if problems:
            self.fail("I6 mid-interval: the nested restart is not bitwise:\n- "
                      + "\n- ".join(problems))

    def test_restart_on_parent_interval_boundary(self) -> None:
        cont = self._continuous(nested=True)
        split = self._split("boundary", nested=True, nsplit=self.SPLIT_ON_BOUNDARY)
        problems = self._compare("on-boundary", cont, split)
        self._penetration("on-boundary", cont, split)
        if problems:
            self.fail("I6 on a parent interval boundary: the nested restart is not "
                      "bitwise:\n- " + "\n- ".join(problems))


class TestI7ZoneIsolation(_NestingCase):
    """Two runs differing **only in the interior**: how far does that reach?

    Construction, in three stages:

    1. a short spin-up from a cold start, which writes an `initd` restart file;
    2. that restart file is copied and, in the copy, a *solenoidal* velocity
       perturbation is added whose support ends exactly at the inner edge of
       the relaxation zone.  It is built from a stream function
       `psi = A Wx(x) Wy(y)` with `u += dpsi/dy`, `v -= dpsi/dx`, where `Wx` and
       `Wy` are raised-cosine windows that are identically zero within
       `L_imp + L_rel` of a lateral face -- so `u` and `v` are untouched
       everywhere in the zone, bit for bit;
    3. both files are warm-started for the same short run, short enough that
       the mean flow carries the perturbation less than one cell.

    Whatever difference then appears inside the zone did not get there by
    advection.  It got there through the projection, which is elliptic and
    therefore instantaneous and global: that is concern C1 of design section 7,
    and this test measures it rather than assuming it away.

    The sharpest number is the difference **on the imposed boundary faces**
    themselves.  `nesting_bcpup` sets those to the parent value identically in
    both runs, so the only thing that can move them is the pressure increment.
    """

    SPEC = mcf.ZONED
    U = 1.0
    DT = 0.125
    SPINUP_STEPS = 16          # 2 s
    COMPARE_STEPS = 8          # 1 s: the mean flow advances 1 m = 1 cell
    AMPLITUDE = 0.05           # m/s, peak perturbation velocity

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.n_guard = int(round(cls.SPEC.guardwidth / cls.SPEC.dx))
        cls.n_zone = int(round((cls.SPEC.guardwidth + cls.SPEC.zonewidth) / cls.SPEC.dx))
        cls.outputs: Dict[str, str] = {}
        cls.fields: Dict[str, Dict[str, np.ndarray]] = {}
        cls.failures: List[str] = []
        try:
            cls._build()
        except Exception as exc:                       # noqa: BLE001
            cls.failures.append(f"I7 setup failed: {exc}")

    # -- the perturbation ------------------------------------------------- #

    @classmethod
    def _window(cls, s: np.ndarray) -> np.ndarray:
        """Raised cosine: exactly 0 for `s <= L_imp + L_rel`, 1 well inside."""
        d0 = cls.SPEC.guardwidth + cls.SPEC.zonewidth
        dw = 3.0 * cls.SPEC.dx
        xi = np.clip((s - d0) / dw, 0.0, 1.0)
        return 0.5 * (1.0 - np.cos(np.pi * xi))

    @classmethod
    def _dwindow(cls, s: np.ndarray, sign: np.ndarray) -> np.ndarray:
        d0 = cls.SPEC.guardwidth + cls.SPEC.zonewidth
        dw = 3.0 * cls.SPEC.dx
        xi = (s - d0) / dw
        inside = (xi > 0.0) & (xi < 1.0)
        return np.where(inside, 0.5 * np.pi / dw * np.sin(np.pi * np.clip(xi, 0.0, 1.0)), 0.0) * sign

    @classmethod
    def _perturb(cls, src: Path, dst: Path) -> Tuple[float, float]:
        """Copy `src` to `dst`, adding the interior-only solenoidal blob.

        Returns `(peak |du|, peak |du| anywhere in the zone)`; the second is a
        self-check and must be exactly zero.
        """
        spec = cls.SPEC
        records = read_fortran_records(src)
        nx, ny, nz = spec.itot + 2, spec.jtot + 2, spec.ktot + 1
        # Array column c is the solver's i = ib - 1 + c, row r is j = jb - 1 + r.
        ci = np.arange(nx)
        rj = np.arange(ny)
        xu, xv = (ci - 1.0) * spec.dx, (ci - 0.5) * spec.dx
        yu, yv = (rj - 0.5) * spec.dy, (rj - 1.0) * spec.dy

        def sd(c: np.ndarray, length: float) -> Tuple[np.ndarray, np.ndarray]:
            """Inward distance from the nearer face, and the inward direction."""
            return np.minimum(c, length - c), np.where(c < 0.5 * length, 1.0, -1.0)

        su_x, nu_x = sd(xu, spec.xlen)
        su_y, nu_y = sd(yu, spec.ylen)
        sv_x, nv_x = sd(xv, spec.xlen)
        sv_y, nv_y = sd(yv, spec.ylen)
        # u = A Wx(x) Wy'(y);  v = -A Wx'(x) Wy(y)
        du = cls.AMPLITUDE * np.outer(cls._dwindow(su_y, nu_y), cls._window(su_x))
        dv = -cls.AMPLITUDE * np.outer(cls._window(sv_y), cls._dwindow(sv_x, nv_x))
        scale = max(np.abs(du).max(), np.abs(dv).max())
        du, dv = du / scale * cls.AMPLITUDE, dv / scale * cls.AMPLITUDE

        out = list(records)
        for rec, delta in ((RESTART_FIELDS["u0"], du), (RESTART_FIELDS["v0"], dv)):
            arr = np.frombuffer(records[rec], dtype="<f8").reshape((nz, ny, nx)).copy()
            arr += delta[None, :, :]
            out[rec] = arr.tobytes()
        with dst.open("wb") as fh:
            for rec in out:
                head = struct.pack("<i", len(rec))
                fh.write(head + rec + head)

        # self-check: nothing may have moved inside the zone
        zone_x = np.minimum(ci, spec.itot + 1 - ci) <= cls.n_zone
        zone_y = np.minimum(rj, spec.jtot + 1 - rj) <= cls.n_zone
        in_zone = zone_x[None, :] | zone_y[:, None]
        leak = max(float(np.abs(du[in_zone]).max()), float(np.abs(dv[in_zone]).max()))
        return max(float(np.abs(du).max()), float(np.abs(dv).max())), leak

    # -- the runs --------------------------------------------------------- #

    def _edits(self, nsteps: int, extra: Optional[Dict[str, str]] = None) -> Dict[str, str]:
        edits = {"dtmax": f"{self.DT:.10g}",
                 "runtime": f"{self.DT * nsteps:.10g}",
                 "trestart": f"{self.DT * nsteps:.10g}",
                 "irandom": "43", "randu": "0.1",
                 "krand": str(self.SPEC.ktot), "lrandomize": ".true."}
        edits.update(extra or {})
        return edits

    @classmethod
    def _build(cls) -> None:
        self = cls  # the helpers below are classmethods in all but name
        spin = cls.root / "i7_spinup"
        mcf.write_case(spin, cls.SPEC, "uniform", times=(0.0, 1.0e6), uprof=cls.U, U=cls.U,
                       edits={"dtmax": f"{cls.DT:.10g}",
                              "runtime": f"{cls.DT * cls.SPINUP_STEPS:.10g}",
                              "trestart": f"{cls.DT * cls.SPINUP_STEPS:.10g}",
                              "irandom": "43", "randu": "0.1",
                              "krand": str(cls.SPEC.ktot), "lrandomize": ".true."})
        done = run_solver(spin)
        if done.returncode != 0:
            raise RuntimeError("spin-up failed\n" + _tail("stdout", done.stdout))
        base = latest_restart(spin, 0, 0)
        cls.start_name = base.name

        peak, leak = cls._perturb(base, cls.root / "i7_perturbed.bin")
        cls.perturbation_peak, cls.perturbation_leak = peak, leak

        for label in ("a", "b"):
            run_dir = cls.root / f"i7_{label}"
            mcf.write_case(run_dir, cls.SPEC, "uniform", times=(0.0, 1.0e6),
                           uprof=cls.U, U=cls.U,
                           edits={"dtmax": f"{cls.DT:.10g}",
                                  "runtime": f"{cls.DT * cls.COMPARE_STEPS:.10g}",
                                  "trestart": f"{cls.DT * cls.COMPARE_STEPS:.10g}",
                                  "lwarmstart": ".true.",
                                  "startfile": f"'{base.name}'"})
            if label == "a":
                shutil.copy2(base, run_dir / base.name)
            else:
                shutil.copy2(cls.root / "i7_perturbed.bin", run_dir / base.name)
            done = run_solver(run_dir)
            if done.returncode != 0:
                raise RuntimeError(f"I7 run {label} failed\n" + _tail("stdout", done.stdout))
            cls.outputs[label] = (done.stdout or "") + (done.stderr or "")
            cls.fields[label] = read_restart_fields(run_dir, cls.SPEC, 1, 1)

    def setUp(self) -> None:
        if self.failures:
            self.fail("\n\n".join(self.failures))

    # -- assertions ------------------------------------------------------- #

    def test_the_perturbation_really_is_interior_only(self) -> None:
        print(f"\n[I7] interior perturbation: peak {self.perturbation_peak:.3e} m/s, "
              f"peak inside the zone {self.perturbation_leak:.3e} m/s", flush=True)
        self.assertEqual(self.perturbation_leak, 0.0,
                         "the initial perturbation is not confined to the interior")

    def _difference(self, name: str) -> np.ndarray:
        return np.abs(self.fields["a"][name] - self.fields["b"][name])

    def test_penetration_profile(self) -> None:
        spec = self.SPEC
        du = interior(self._difference("u0"))
        idx_i = np.arange(spec.itot)
        idx_j = np.arange(spec.jtot)
        di = np.minimum(idx_i, spec.itot - 1 - idx_i).astype(float)
        dj = np.minimum(idx_j, spec.jtot - 1 - idx_j).astype(float)
        dist = np.minimum(np.broadcast_to(di[None, None, :], du.shape),
                          np.broadcast_to(dj[None, :, None], du.shape))
        print("[I7] max|u_a - u_b| by inward distance from the lateral boundary "
              f"(guard 0..{self.n_guard - 1}, ramp {self.n_guard}..{self.n_zone - 1}):",
              flush=True)
        profile = []
        for d in range(int(dist.max()) + 1):
            value = float(du[dist == d].max())
            profile.append((d, value))
            where = "guard" if d < self.n_guard else ("ramp" if d < self.n_zone else "interior")
            print(f"[I7]   d = {d:2d} ({where:8s}): {value:.3e}", flush=True)
        guard = max(v for d, v in profile if d < self.n_guard)
        ramp = max(v for d, v in profile if self.n_guard <= d < self.n_zone)
        inner = max(v for d, v in profile if d >= self.n_zone)
        print(f"[I7] peaks: guard {guard:.3e}  ramp {ramp:.3e}  interior {inner:.3e}  "
              f"guard/interior = {guard / inner:.3e}", flush=True)
        self.assertGreater(inner, 0.1 * self.AMPLITUDE,
                           "the two runs barely differ, so I7 proves nothing")
        self.assertLess(guard, 0.05 * inner,
                        "the interior difference has reached the guard strip")
        self.assertLess(ramp, inner, "the ramp is not attenuating the difference")

    def test_imposed_faces_move_only_through_the_pressure(self) -> None:
        """The C1 number: the imposed faces are identical in both runs by
        construction, so any difference there is the global pressure response."""
        spec = self.SPEC
        du = self._difference("u0")
        dv = self._difference("v0")
        faces = {
            "u west (i = ib)":    float(np.max(du[:spec.ktot, 1:-1, 1])),
            "u east (i = ie+1)":  float(np.max(du[:spec.ktot, 1:-1, -1])),
            "v south (j = jb)":   float(np.max(dv[:spec.ktot, 1, 1:-1])),
            "v north (j = je+1)": float(np.max(dv[:spec.ktot, -1, 1:-1])),
        }
        inner = float(np.max(interior(du)))
        for label, value in faces.items():
            print(f"[I7] {label}: max|u_a - u_b| = {value:.3e}  "
                  f"({value / inner:.2e} of the interior difference)", flush=True)
        worst = max(faces.values())
        self.assertLess(worst, 0.02 * inner,
                        "the imposed boundary faces move with the interior")

    def test_reported_zone_diagnostics(self) -> None:
        for label in ("a", "b"):
            stats = parse_nesting_stats(self.outputs[label])
            self.assertTrue(stats, f"no nesting_stats output for run {label}")
            s = stats[-1]
            print(f"[I7] run {label}: |grad p| zone = {s['gradp_zone']:.3e}  "
                  f"interior = {s['gradp_interior']:.3e}  ratio = {s['gradp_ratio']:.3f}  "
                  f"E_guard = {s['e_guard']:.3e}  E_relax = {s['e_relax']:.3e}", flush=True)
            self.assertGreater(s["gradp_interior"], 0.0,
                               "the interior |grad p| diagnostic is dead")
            self.assertNotEqual(s["e_relax"], 0.0,
                                "the relaxation ramp injected no energy, so it is inactive")
        a = parse_nesting_stats(self.outputs["a"])[-1]
        b = parse_nesting_stats(self.outputs["b"])[-1]
        print(f"[I7] |grad p| ratio (zone/interior): run a {a['gradp_ratio']:.3f}, "
              f"run b {b['gradp_ratio']:.3f}", flush=True)


# --------------------------------------------------------------------------- #
# The committed cube case, tests/cases/064, used by I8 and the 2x2 tests
# --------------------------------------------------------------------------- #

CASE_064 = "064"
CASE_064_DIR = REPO_ROOT / "tests" / "cases" / CASE_064

#: 3 m guard + 20 m ramp = 23 m, tau = 4 s; the cube's windward face is at
#: x = 24 m, one cell beyond the inner edge of the zone.
CUBE_SPEC = mcf.CaseSpec(itot=64, jtot=64, ktot=64, xlen=64.0, ylen=64.0, zsize=64.0,
                         nzone=24, guardwidth=3.0, zonewidth=20.0, tau=4.0)


def prepare_case_064(run_dir: Path, nested: bool, runtime: float, dt: float,
                     nprocx: int = 1, nprocy: int = 1, trestart: float = 1.0e9,
                     extra: Optional[Dict[str, str]] = None,
                     field: str = "uniform", times: Sequence[float] = (0.0, 1.0e6),
                     **field_kwargs) -> Path:
    """Copy ``tests/cases/064`` into ``run_dir`` as a nested (or periodic) run.

    Energy balance, temperature and moisture are switched off so the run is
    purely mechanical and takes seconds.  The nested run is driven by the
    parent; the periodic reference by the volume-flow-rate controller at the
    same 1 m/s.  Keys absent from ``namoptions.064`` (``lwarmstart``,
    ``startfile``) are inserted into ``&RUN``.
    """
    if not CASE_064_DIR.is_dir():
        raise unittest.SkipTest(f"case {CASE_064} not found")
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(CASE_064_DIR, run_dir)
    nml = run_dir / f"namoptions.{CASE_064}"
    text = nml.read_text(encoding="utf-8")
    settings = {
        "runtime": f"{runtime:.10g}",
        "dtmax": f"{dt:.10g}",
        "trestart": f"{trestart:.10g}",
        "ladaptive": ".false.",
        "nprocx": str(nprocx), "nprocy": str(nprocy),
        "lEB": ".false.", "ltempeq": ".false.", "lmoist": ".false.",
        "lbuoyancy": ".false.", "lxytdump": ".false.",
        "luvolflowr": ".false." if nested else ".true.",
    }
    settings.update(extra or {})
    for key, value in settings.items():
        text, n = re.subn(rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$",
                          lambda m: m.group(1) + value, text)
        if n == 0:
            text = text.replace("&RUN\n", f"&RUN\n{key} = {value}\n", 1)
    text = text.replace("&WALLS\n", "&WALLS\nlwritefac = .true.\n"
                        f"dtfac = {runtime:.10g}\n")
    text = text.replace("&BC\n", "&BC\n"
                        f"BCxm = {4 if nested else 1}\n"
                        f"BCym = {3 if nested else 1}\nBCtopm = 1\n")
    spec = CUBE_SPEC
    text += (f"\n&NESTING\nlnesting = {'.true.' if nested else '.false.'}\n"
             f"nest_guardwidth = {spec.guardwidth:.10g}\n"
             f"nest_zonewidth = {spec.zonewidth:.10g}\n"
             f"nest_tau = {spec.tau:.10g}\n"
             "nest_shape = 1\nnest_timeinterp = 2\nnest_nwall = 1\n"
             "nest_lparentgeom = .false.\n"
             "nest_fluxtol = 1.e-10\nnest_lfluxassert = .true.\n/\n")
    nml.write_text(text, encoding="utf-8")
    if nested:
        data = mcf.build_nesting_data(spec, field, list(times), **field_kwargs)
        from udprep.nesting import write_nesting_file
        write_nesting_file(run_dir / f"nesting.inp.{CASE_064}.nc", data, override=True)
    return run_dir


# --------------------------------------------------------------------------- #
# I8 -- IBM interaction
# --------------------------------------------------------------------------- #


class TestI8IbmInteraction(_NestingCase):
    """A building right at the inner edge of the zone, nested versus not.

    `tests/cases/064` is a single 6 m cube in a 64 x 64 x 64 m box, with its
    windward face at `x = 24 m`, so a `3 + 20 = 23 m` zone puts the inner edge
    exactly one cell upstream of the building -- design section 10.3's
    "buildings adjacent to the zone edge", and the tightest arrangement the
    building-free rule of design section 5 allows.  `nest_lparentgeom` is left
    at `.false.` on purpose: if any part of the geometry did fall inside the
    zone, `nesting_init` would abort and this test would say so.

    The reference is the same case, same initial condition, driven periodically
    at the same volume flow rate.  This is a **modelling** comparison, not an
    exactness one: the two runs have genuinely different boundary conditions,
    so the tolerance is a stated one and the measured value is printed.  What
    would be a real failure is a large difference concentrated on the first
    building row, which would mean the zone is contaminating the geometry it
    was supposed to stand off from.

    Energy balance, temperature and moisture are switched off so the comparison
    is purely mechanical and the case runs in a few seconds.
    """

    CASE_ID = 64
    CASE_NAME = CASE_064
    SPEC = CUBE_SPEC
    RUNTIME, DT = 4.0, 0.25
    #: fraction of the reference stress the two runs may differ by
    TOLERANCE = 0.10
    #: which stress components are non-trivial on which faces of the cube.
    #: `tau_<n>` is identically zero on a face whose normal is `n` (the normal
    #: load is `pres`), so each face is judged on the components it carries.
    FACES = {
        "windward (-x)": (lambda n: n[:, 0] < -0.5, ("tau_y", "tau_z", "pres")),
        "south side (-y)": (lambda n: n[:, 1] < -0.5, ("tau_x", "tau_z", "pres")),
        "north side (+y)": (lambda n: n[:, 1] > 0.5, ("tau_x", "tau_z", "pres")),
    }
    REPORT_ONLY = {"leeward (+x)": (lambda n: n[:, 0] > 0.5, ("tau_y", "tau_z", "pres"))}

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.case_source = CASE_064_DIR
        cls.error: Optional[str] = None
        cls.out: Dict[str, str] = {}
        cls.dirs: Dict[str, Path] = {}
        try:
            cls._build()
        except Exception as exc:                       # noqa: BLE001
            cls.error = f"I8 setup failed: {exc}"

    @classmethod
    def _prepare(cls, label: str, nested: bool) -> Path:
        return prepare_case_064(cls.root / f"i8_{label}", nested,
                                runtime=cls.RUNTIME, dt=cls.DT)

    @classmethod
    def _build(cls) -> None:
        for label, nested in (("nested", True), ("reference", False)):
            run_dir = cls._prepare(label, nested)
            done = run_solver(run_dir, nprocs=1, namelist=f"namoptions.{cls.CASE_NAME}")
            if done.returncode != 0:
                raise RuntimeError(f"the {label} run exited {done.returncode}\n"
                                   + _tail("stdout", done.stdout) + "\n"
                                   + _tail("stderr", done.stderr))
            cls.dirs[label] = run_dir
            cls.out[label] = (done.stdout or "") + (done.stderr or "")

    def setUp(self) -> None:
        if self.error:
            self.fail(self.error)

    # -- facets ----------------------------------------------------------- #

    def _normals(self) -> np.ndarray:
        rows = []
        with (self.dirs["nested"] / f"facets.inp.{self.CASE_NAME}").open() as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 4 or parts[0].startswith("#"):
                    continue
                rows.append([float(x) for x in parts[1:4]])
        return np.array(rows)

    def _facet_stress(self, label: str) -> Dict[str, np.ndarray]:
        import netCDF4 as nc
        with nc.Dataset(self.dirs[label] / f"fac.{self.CASE_NAME}.nc") as ds:
            return {v: np.asarray(ds.variables[v][:])[-1].astype(np.float64)
                    for v in ("tau_x", "tau_y", "tau_z", "pres")}

    def test_the_zone_is_building_free(self) -> None:
        """`nest_lparentgeom = .false.` means init aborts if it is not."""
        out = self.out["nested"]
        self.assertIn("modnesting: zone points", out, "nesting never initialised")
        self.assertNotIn("solid points found inside the relaxation zone", out)
        self.assertNotIn("solid points inside the relaxation zone (allowed", out)
        margin = 24.0 - (self.SPEC.guardwidth + self.SPEC.zonewidth)
        print(f"\n[I8] zone inner edge is {margin:g} m ({margin / self.SPEC.dx:g} cells) "
              "upstream of the building", flush=True)
        self.assertLessEqual(margin, 2.0 * self.SPEC.dx,
                             "the building is not adjacent to the zone edge, so I8 "
                             "is not testing what it claims to")

    def test_facet_stresses_match_the_reference(self) -> None:
        """Every non-trivial stress component on the cube's windward and side
        faces, nested versus periodic, within the stated tolerance.

        `tau_x` -- the streamwise shear, which is the stress the flow past a
        cube is mostly about -- lives on the side faces and the roof, not on
        the windward face where its normal load is `pres`; each face is
        therefore compared on the components it actually carries, and a
        component whose reference is identically zero is a failure of the test
        rather than a free pass.
        """
        normals = self._normals()
        nested = self._facet_stress("nested")
        reference = self._facet_stress("reference")
        failures = []
        groups = dict(self.FACES)
        groups.update(self.REPORT_ONLY)
        for face, (select, components) in groups.items():
            idx = np.where(select(normals))[0]
            self.assertGreater(idx.size, 0, f"no facets found for the {face} face")
            for name in components:
                a, b = nested[name][idx], reference[name][idx]
                scale = float(np.max(np.abs(b)))
                if scale <= 0.0:
                    failures.append(f"{face} {name}: the reference stress is identically "
                                    "zero, so this component cannot be compared here")
                    continue
                rel = float(np.max(np.abs(a - b))) / scale
                judged = face in self.FACES
                print(f"[I8] {face:16s} {name:5s}: nested {a.mean():+.4e}  "
                      f"reference {b.mean():+.4e}  relative diff {rel:.3e}"
                      f"{'' if judged else '   (reported, not judged)'}", flush=True)
                if judged and rel > self.TOLERANCE:
                    failures.append(f"{face} {name}: {rel:.3e} > {self.TOLERANCE:.2f}")
        if failures:
            self.fail("facet stresses on the building differ from the no-nesting "
                      "reference by more than the stated tolerance:\n- "
                      + "\n- ".join(failures))

    def test_reported_zone_diagnostics_with_buildings(self) -> None:
        stats = parse_nesting_stats(self.out["nested"])
        self.assertTrue(stats, "no nesting_stats output from the nested run")
        s = stats[-1]
        print(f"[I8] |grad p| zone = {s['gradp_zone']:.3e}  "
              f"interior = {s['gradp_interior']:.3e}  ratio = {s['gradp_ratio']:.3f}  "
              f"misfit = {s['misfit']:.3e}  "
              f"E_guard = {s['e_guard']:.3e}  E_relax = {s['e_relax']:.3e}", flush=True)
        divs = parse_divergence(self.out["nested"])
        worst = max(abs(a) for a, _ in divs) if divs else 0.0
        print(f"[I8] nested run max divmax = {worst:.3e}, |Phi| = {abs(s['phi']):.3e}",
              flush=True)
        self.assertLessEqual(worst, 1.0e-12,
                             "the nested IBM run does not solve to a divergence-free field")
        self.assertLessEqual(abs(s["phi"]), 1.0e-10, "the net boundary flux is not zero")
        self.assertLess(s["gradp_ratio"], 2.0,
                        "the projection works harder in the zone than in the interior")



# --------------------------------------------------------------------------- #
# I5 on 2x2 with a real zone and a cube -- the multi-rank case that runs in CI
# --------------------------------------------------------------------------- #


class TestI5CubeParity2x2(_NestingCase):
    """`tests/cases/064` on 1x1 and 2x2 with a ramp, `tau > 0` and a cube.

    `TestI5DecompositionParity` runs with `W == 1` everywhere, `tau = 0` and
    no buildings, so it never exercises the ramp weights, the relaxation
    factor or the IBM masking on a decomposed zone.  This does: a 3 m guard,
    a 20 m ramp, `tau = 4 s`, and a cube whose windward face sits one cell
    inside the zone's inner edge -- on one rank and on four, oversubscribed
    on a 2-4 core runner.  The fields must agree to `PARITY_TOL`.

    Sized for CI: `NSTEPS` steps of a 64^3 case take ~15 s (1x1) and ~11 s
    (2x2) on a gfortran Debug build, about half of it initialisation.
    """

    NSTEPS, DT = 4, 0.25
    LAYOUTS = {"serial": (1, 1), "xy_split": (2, 2)}

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.fields: Dict[str, Dict[str, np.ndarray]] = {}
        cls.out: Dict[str, str] = {}
        cls.dirs: Dict[str, Path] = {}
        cls.failures: List[str] = []
        runtime = cls.NSTEPS * cls.DT
        for label, (nprocx, nprocy) in cls.LAYOUTS.items():
            cls.dirs[label] = run_dir = prepare_case_064(
                cls.root / f"i5cube_{label}", nested=True, runtime=runtime, dt=cls.DT,
                nprocx=nprocx, nprocy=nprocy, trestart=runtime,
                extra={"lrandomize": ".false."})
            done = run_solver(run_dir, nprocs=nprocx * nprocy,
                              namelist=f"namoptions.{CASE_064}")
            cls.out[label] = (done.stdout or "") + (done.stderr or "")
            if done.returncode != 0:
                cls.failures.append(f"{label} ({nprocx}x{nprocy}) exited {done.returncode}\n"
                                    + _tail("output", cls.out[label]))
                continue
            cls.fields[label] = read_restart_fields(run_dir, CUBE_SPEC, nprocx, nprocy,
                                                    expnr=CASE_064)

    def setUp(self) -> None:
        if self.failures:
            self.fail("\n\n".join(self.failures))

    def test_fields_agree_on_2x2(self) -> None:
        reference, candidate = self.fields["serial"], self.fields["xy_split"]
        problems = []
        for name in ("u0", "v0", "w0", "pres0"):
            diff = np.abs(interior(candidate[name]) - interior(reference[name]))
            worst = float(np.nanmax(diff))
            idx = np.unravel_index(int(np.nanargmax(diff)), diff.shape)
            print(f"[I5 cube] 2x2 vs 1x1 {name}: max abs diff {worst:.3e} at {idx}",
                  flush=True)
            if worst > PARITY_TOL:
                problems.append(f"{name}: {worst:.3e} at {idx}")
        if problems:
            self.fail("decomposition parity with a cube and a ramp failed:\n- "
                      + "\n- ".join(problems))

    def test_the_case_has_the_ramp_and_the_cube(self) -> None:
        """Guard against a vacuous pass: the ramp relaxed and the cube is there."""
        for label in self.LAYOUTS:
            stats = parse_nesting_stats(self.out[label])
            self.assertTrue(stats, f"no nesting_stats output on {label}")
            s = stats[-1]
            print(f"[I5 cube] {label}: E_relax = {s['e_relax']:.3e}  "
                  f"E_guard = {s['e_guard']:.3e}  misfit = {s['misfit']:.3e}", flush=True)
            self.assertNotEqual(s["e_relax"], 0.0,
                                f"the relaxation ramp injected no energy on {label}")
            # The cube: its facets (fac.<expnr>.nc, lwritefac) carry a pressure
            # load only if the immersed boundary was read and the flow hit it.
            import netCDF4 as nc
            with nc.Dataset(self.dirs[label] / f"fac.{CASE_064}.nc") as ds:
                pres = np.asarray(ds.variables["pres"][:])[-1]
            print(f"[I5 cube] {label}: {pres.size} facets, max |pres| = "
                  f"{float(np.max(np.abs(pres))):.3e}", flush=True)
            self.assertGreater(float(np.max(np.abs(pres))), 1.0e-3,
                               f"no facet feels the flow on {label}, so there is no cube")


# --------------------------------------------------------------------------- #
# I6 on 2x2 with a cube
# --------------------------------------------------------------------------- #


class TestI6CubeRestartParity2x2(_NestingCase):
    """`N` steps == `N/2` + restart + `N/2`, bitwise, on four ranks with a cube.

    `TestI6RestartParity` establishes the restart on one rank with a uniform
    field and no buildings.  This repeats the two split points -- one strictly
    inside a parent interval, one exactly on a parent level -- on 2x2 ranks of
    `tests/cases/064`, so the per-rank restart files, the decomposed parent
    buffer and the IBM state all have to come back exactly.

    `dt = 0.25 s` and parent levels every `1.25 s` are exact binary
    fractions, so both legs take exactly the intended number of steps.
    """

    DT = 0.25
    NSTEPS = 8
    PARENT_DT = 1.25
    #: step 5 -> t = 1.25 s = parent level 1 exactly; step 3 -> t = 0.75 s.
    SPLIT_ON_BOUNDARY = 5
    SPLIT_MID_INTERVAL = 3
    NPROCX, NPROCY = 2, 2
    U0, AMP, PERIOD = 1.0, 0.25, 8.0

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.times = np.arange(0.0, cls.DT * cls.NSTEPS + 2 * cls.PARENT_DT, cls.PARENT_DT)
        cls._cache: Dict[str, Path] = {}

    def _leg(self, name: str, nsteps: int, startfile: Optional[str] = None) -> Path:
        extra = {"lrandomize": ".false."}
        if startfile is not None:
            extra.update({"lwarmstart": ".true.", "startfile": f"'{startfile}'"})
        run_dir = prepare_case_064(
            self.root / name, nested=True, runtime=self.DT * nsteps, dt=self.DT,
            nprocx=self.NPROCX, nprocy=self.NPROCY, trestart=self.DT * nsteps,
            extra=extra, field="unsteady_uniform", times=self.times,
            U0=self.U0, A=self.AMP, T=self.PERIOD)
        if startfile is not None:
            source = self.root / name.replace("_second", "_first")
            for f in list(source.glob("initd*")) + list(source.glob("inits*")):
                shutil.copy2(f, run_dir / f.name)
        check_ok(self, run_solver(run_dir, nprocs=self.NPROCX * self.NPROCY,
                                  namelist=f"namoptions.{CASE_064}"), f"I6 cube leg {name}")
        latest = latest_restart(run_dir, 0, 0, expnr=CASE_064)
        ntrun = int(latest.name[5:13])
        self.assertEqual(ntrun, nsteps if startfile is None else self.NSTEPS,
                         f"{name}: stopped at step {ntrun}")
        return run_dir

    def _continuous(self) -> Path:
        if "continuous" not in self._cache:
            self._cache["continuous"] = self._leg("i6cube_continuous", self.NSTEPS)
        return self._cache["continuous"]

    def _split(self, label: str, nsplit: int) -> Path:
        self._leg(f"i6cube_{label}_first", nsplit)
        return self._leg(f"i6cube_{label}_second", self.NSTEPS - nsplit,
                         startfile=f"initd{nsplit:08d}_000_000.{CASE_064}")

    def _compare(self, label: str, cont: Path, split: Path) -> List[str]:
        problems = []
        for px in range(self.NPROCX):
            for py in range(self.NPROCY):
                a = read_fortran_records(latest_restart(cont, px, py, expnr=CASE_064))
                b = read_fortran_records(latest_restart(split, px, py, expnr=CASE_064))
                for name, rec in RESTART_FIELDS.items():
                    x = np.frombuffer(a[rec], dtype="<f8")
                    y = np.frombuffer(b[rec], dtype="<f8")
                    worst = float(np.max(np.abs(x - y)))
                    same = a[rec] == b[rec]
                    if not same:
                        problems.append(f"rank ({px},{py}) {name}: max abs diff {worst:.3e}")
        print(f"[I6 cube] {label}: {len(RESTART_FIELDS) * self.NPROCX * self.NPROCY} "
              f"restart records over 4 ranks, {len(problems)} differ", flush=True)
        return problems

    def test_restart_mid_parent_interval(self) -> None:
        cont = self._continuous()
        split = self._split("mid", self.SPLIT_MID_INTERVAL)
        problems = self._compare("mid-interval", cont, split)
        if problems:
            self.fail("I6 on 2x2 with a cube, mid-interval: the restart is not bitwise:\n- "
                      + "\n- ".join(problems))

    def test_restart_on_parent_interval_boundary(self) -> None:
        cont = self._continuous()
        split = self._split("boundary", self.SPLIT_ON_BOUNDARY)
        problems = self._compare("on-boundary", cont, split)
        if problems:
            self.fail("I6 on 2x2 with a cube, on a parent level: the restart is not "
                      "bitwise:\n- " + "\n- ".join(problems))


# --------------------------------------------------------------------------- #
# I9 -- cold start from the parent (design section 10.6 item 4)
# --------------------------------------------------------------------------- #


class TestI9ColdStartFromParent(_NestingCase):
    """A cold start must be able to begin *at* the parent, not at `prof.inp`.

    The case is the ZONED one, with a real guard strip and ramp, so most of the
    domain is **not** imposed: what the interior holds is what the initial
    condition put there and nothing else.

    How "reproduces the parent to round-off at `t = 0`" is established.  A
    restart file can only be written *after* a step, so the field at `t = 0` is
    not directly readable from a run.  It is pinned instead in two independent
    ways:

    * runmode 1011 (U40) reads `u0`/`um`, `v0`/`vm`, `w0`/`wm` straight after
      `nesting_init` and compares every point against the stored block -- a
      non-separable analytic 3-D field, on 1x1, 2x1, 1x2 and 2x2.  Measured
      max error 1.1e-16.
    * here, end to end: a run whose interior comes from the parent block is
      required to be **bitwise identical** to a run whose interior comes from a
      `prof.inp` carrying the same profile.  Two runs of this solver cannot
      agree bit for bit unless they started from the same bits.

    The parent is `uniform`, which the second construction needs (`prof.inp` is
    a profile, so only a horizontally uniform field can be expressed both ways).
    Note that it is *not* a fixed point of the run: `closurebc` gives the ground
    a no-slip molecular viscosity regardless of `BCbotm`
    (`src/modboundary.f90:466`), which decelerates the lowest layer by
    ~3e-5 m/s over the 2 s of the run.  That is physics, it is identical in both
    runs, and it is why this test compares two runs rather than comparing one
    run against `U`.
    """

    SPEC = mcf.ZONED
    U = 1.0
    #: (label, uprof, nest_linitfromparent)
    RUNS = (("from_parent", 0.0, ".true."),
            ("from_prof_equal", 1.0, ".false."),
            ("from_prof_zero", 0.0, ".false."))

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.out: Dict[str, str] = {}
        cls.dirs: Dict[str, Path] = {}
        cls.data = None
        for label, uprof, switch in cls.RUNS:
            run_dir = cls.root / f"i9_{label}"
            data = mcf.write_case(
                run_dir, cls.SPEC, "uniform", times=(0.0, 1.0e6),
                uprof=uprof, U=cls.U, initial=True,
                edits={"runtime": "2.", "dtmax": "0.5", "trestart": "2.",
                       "nest_linitfromparent": switch},
            )
            cls.data = cls.data or data
            done = run_solver(run_dir)
            cls.dirs[label] = run_dir
            cls.out[label] = (done.stdout or "") + (done.stderr or "")
            if done.returncode != 0:
                raise RuntimeError(f"I9 run '{label}' failed\n"
                                   + _tail("output", cls.out[label]))

    def test_the_stored_block_is_solenoidal(self) -> None:
        """What the solver reads must already be divergence free on the child grid.

        Checked on the fixture itself, before any solver is involved, so a
        failure here is the writer's projection and not the reader's.
        """
        block = self.data.initial_fields
        div = mcf.discrete_divergence(
            self.SPEC.grid(), (block["u"], block["v"], block["w"])
        )
        peak = float(np.max(np.abs(div)))
        print(f"\n[I9] max |div| of the stored initial condition = {peak:.3e}", flush=True)
        self.assertLessEqual(peak, ROUNDOFF, "the writer stored a divergent initial condition")

    def test_it_is_bitwise_an_equivalent_prof_inp_start(self) -> None:
        a = read_fortran_records(latest_restart(self.dirs["from_parent"], 0, 0))
        b = read_fortran_records(latest_restart(self.dirs["from_prof_equal"], 0, 0))
        self.assertEqual(len(a), len(b), "the two runs wrote different restart records")
        differing = [name for name, rec in RESTART_FIELDS.items() if a[rec] != b[rec]]
        worst = 0.0
        for rec in RESTART_FIELDS.values():
            x = np.frombuffer(a[rec], dtype="<f8")
            y = np.frombuffer(b[rec], dtype="<f8")
            worst = max(worst, float(np.max(np.abs(x - y))))
        print(f"[I9] parent-initialised vs equivalent prof.inp start: "
              f"max abs restart difference {worst:.3e}, differing records {differing}",
              flush=True)
        self.assertEqual(differing, [],
                         "the cold start from the parent block is not bit-for-bit the "
                         "same state as a prof.inp start carrying the same field")

    def test_the_first_projection_is_already_clean(self) -> None:
        divs = parse_divergence(self.out["from_parent"])
        self.assertTrue(divs, "the run reported no divergence diagnostics")
        first, worst = divs[0], max(abs(a) for a, _ in divs)
        stats = parse_nesting_stats(self.out["from_parent"])
        control = parse_nesting_stats(self.out["from_prof_zero"])
        self.assertTrue(stats and control, "no nesting_stats output")
        gp, gp0 = stats[0]["gradp_interior"], control[0]["gradp_interior"]
        mis, mis0 = stats[0]["misfit"], control[0]["misfit"]
        print(f"[I9] first divmax = {first[0]:.3e}, divtot = {first[1]:.3e}; "
              f"worst divmax over the run = {worst:.3e}", flush=True)
        print(f"[I9] first substep |grad p|_interior = {gp:.3e} (from prof.inp: "
              f"{gp0:.3e}); zone misfit {mis:.3e} (from prof.inp: {mis0:.3e})", flush=True)
        self.assertLessEqual(abs(first[0]), ROUNDOFF,
                             "the field is not divergence free after the first projection")
        self.assertLessEqual(worst, ROUNDOFF)
        # The initial condition is solenoidal AND consistent with the imposed
        # boundary, so the first projection has nothing of its own to remove.
        self.assertLess(gp, 1.0e-3 * max(gp0, 1.0e-30),
                        "the first projection worked as hard as it does from prof.inp, "
                        "so the initial condition was not consistent with the boundary")
        # The misfit is measured after a substep has been taken, so it cannot be
        # at round-off -- the ground's molecular no-slip has already moved the
        # lowest layer by ~1e-6 m/s (see this class's docstring). It is reported
        # rather than asserted; what is asserted is the ratio to the control.
        self.assertLess(mis, 1.0e-3 * max(mis0, 1.0e-30),
                        "the zone starts as far from the parent as a prof.inp start does")

    def test_the_switch_is_what_did_it(self) -> None:
        """The control: without the switch the interior starts from `prof.inp`.

        Without this the tests above could pass on a case where `prof.inp`
        happened to agree with the parent, which would make them vacuous.
        """
        f = read_restart_fields(self.dirs["from_prof_zero"], self.SPEC, 1, 1)
        du = float(np.max(np.abs(interior(f["u0"]) - self.U)))
        dp = float(np.max(np.abs(interior(f["pres0"]))))
        g = read_restart_fields(self.dirs["from_parent"], self.SPEC, 1, 1)
        du_on = float(np.max(np.abs(interior(g["u0"]) - self.U)))
        print(f"[I9] with the switch off: max|u-U| = {du:.3e}  max|pres0| = {dp:.3e};  "
              f"with it on: max|u-U| = {du_on:.3e}", flush=True)
        self.assertGreater(du, 1.0e-3,
                           "prof.inp and the parent agree, so I9 proves nothing")
        self.assertLess(du_on, 1.0e-3 * du,
                        "the switch made no appreciable difference")
        self.assertIn("cold start initialised from the parent", self.out["from_parent"])
        self.assertNotIn("cold start initialised from the parent",
                         self.out["from_prof_zero"])


# --------------------------------------------------------------------------- #
# I10 -- the leaky lid (design case B, section 10.6 item 5)
# --------------------------------------------------------------------------- #


class TestI10LeakyLid(_NestingCase):
    """`BCtopm_pressure` must run with the flux assertion **on**, and correctly.

    Design section 3.2 case B: `bcpup` sets the predicted lid velocity from the
    accumulated pressure and `tderive` adds the matching increment, which is
    exactly the Dirichlet-in-the-mean-mode row the solver pins (F3).  The
    projection there is complete for any net flux, so `Phi_total = 0` is not a
    requirement -- the lid flux is the child breathing against its reservoir.
    What must still vanish is the flux through the faces the scheme controls.

    Making the lid actually breathe takes some care, because the feedback is
    autonomous and starts from rest: the lid velocity is driven by the
    horizontal mean of the accumulated pressure at `k = ke`, that mean is
    proportional to `-Phi_total` through the pin, and with a flux-balanced
    parent and a cold start both are identically zero for ever.  So the case is
    warm-started from a restart file whose `pres0` carries a uniform offset --
    a child arriving with a column-pressure excess, which is precisely the mass
    excess case B exists to let out.

    The test then requires: the run completes with the assertion on; the lid
    flux is large enough that the old six-face assertion **would** have fired
    (otherwise nothing was fixed); the closed faces stay balanced to round-off;
    and the post-projection divergence stays at round-off, which is case B's
    substantive claim.
    """

    SPEC = mcf.BASE
    DT = 0.125
    SPINUP_STEPS = 8
    COMPARE_STEPS = 8
    #: uniform offset added to `pres0` in the restart file [m2 s-2]
    POFFSET = 0.5

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        spin = cls.root / "i10_spinup"
        mcf.write_case(spin, cls.SPEC, "face_forced", times=(0.0, 1.0e6), uprof=1.0,
                       edits={"dtmax": f"{cls.DT:.10g}",
                              "runtime": f"{cls.DT * cls.SPINUP_STEPS:.10g}",
                              "trestart": f"{cls.DT * cls.SPINUP_STEPS:.10g}",
                              "BCtopm": "3"})
        done = run_solver(spin)
        if done.returncode != 0:
            raise RuntimeError("I10 spin-up failed\n"
                               + _tail("output", (done.stdout or "") + (done.stderr or "")))
        source = latest_restart(spin, 0, 0)

        cls.out: Dict[str, str] = {}
        cls.rc: Dict[str, int] = {}
        for label, offset in (("breathing", cls.POFFSET), ("control", 0.0)):
            run_dir = cls.root / f"i10_{label}"
            mcf.write_case(run_dir, cls.SPEC, "face_forced", times=(0.0, 1.0e6), uprof=1.0,
                           edits={"dtmax": f"{cls.DT:.10g}",
                                  "runtime": f"{cls.DT * (cls.SPINUP_STEPS + cls.COMPARE_STEPS):.10g}",
                                  "trestart": "1.e9",
                                  "BCtopm": "3",
                                  "lwarmstart": ".true.",
                                  "startfile": f"'{source.name}'"})
            cls._offset_pres0(source, run_dir / source.name, offset)
            done = run_solver(run_dir)
            cls.out[label] = (done.stdout or "") + (done.stderr or "")
            cls.rc[label] = done.returncode

    @classmethod
    def _offset_pres0(cls, src: Path, dst: Path, offset: float) -> None:
        """Copy a restart file, adding a uniform offset to the `pres0` record."""
        records = list(read_fortran_records(src))
        rec = RESTART_FIELDS["pres0"]
        arr = np.frombuffer(records[rec], dtype="<f8").copy() + offset
        records[rec] = arr.tobytes()
        with dst.open("wb") as fh:
            for r in records:
                head = struct.pack("<i", len(r))
                fh.write(head + r + head)

    def test_the_run_completes_with_the_assertion_on(self) -> None:
        if self.rc["breathing"] != 0:
            self.fail("a BCtopm_pressure run still trips the flux assertion\n"
                      + _tail("output", self.out["breathing"]))
        self.assertNotIn("boundary flux residual out of tolerance", self.out["breathing"])

    def test_the_lid_breathes_and_the_old_assertion_would_have_fired(self) -> None:
        stats = parse_nesting_stats(self.out["breathing"])
        self.assertTrue(stats, "no nesting_stats output")
        self.assertTrue(all("phi_lid" in s for s in stats),
                        "nesting_stats does not report the lid flux")
        worst_lid = max(abs(s["phi_lid"]) for s in stats)
        worst_all = max(abs(s["phi"]) for s in stats)
        worst_closed = max(abs(s["phi_closed"]) for s in stats)
        control = parse_nesting_stats(self.out["control"])
        control_all = max(abs(s["phi"]) for s in control) if control else 0.0
        print(f"\n[I10] with a column-pressure offset: max |Phi| (six faces) = "
              f"{worst_all:.3e}   max |Phi_lid| = {worst_lid:.3e}   "
              f"max |Phi_closed| = {worst_closed:.3e}", flush=True)
        print(f"[I10] without the offset (control): max |Phi| = {control_all:.3e}",
              flush=True)
        self.assertGreater(worst_lid, 1.0e-6,
                           "the lid carries no flux, so this case does not exercise case B")
        # This is the number the pre-fix assertion tested, and it is what used to
        # force nest_lfluxassert = .false. on every case B run.
        self.assertGreater(worst_all, 1.0e-10,
                           "the six-face residual is within tolerance, so the old "
                           "assertion would not have fired and nothing was fixed")
        self.assertLessEqual(worst_closed, 1.0e-10,
                             "the imposed faces are not flux balanced")

    def test_the_projection_is_still_complete(self) -> None:
        """Case B's substantive claim: the lid realises the flux the pressure
        implies, so the post-projection divergence stays at round-off -- unlike
        case A, where an unbalanced flux leaves a source in the top cell layer."""
        divs = parse_divergence(self.out["breathing"])
        self.assertTrue(divs, "no divergence diagnostics")
        worst = max(abs(a) for a, _ in divs)
        worst_tot = max(abs(b) for _, b in divs)
        print(f"[I10] max divmax = {worst:.3e}, max divtot = {worst_tot:.3e}", flush=True)
        self.assertLessEqual(worst, 1.0e-12,
                             "case B does not leave a divergence-free field")


if __name__ == "__main__":
    unittest.main(verbosity=2)
