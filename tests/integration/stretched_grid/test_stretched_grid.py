"""Stretched vertical grid: preprocessing and solver integration test.

The committed energy-balance case 064 is preprocessed again with an
exponentially stretched vertical grid and then run for a few time steps. No
other committed case uses ``lzstretch``, so this is the test that exercises the
stretched-grid code paths end to end:

* the grid generator and the profile written to ``prof.inp``;
* the IBM preprocessing, whose facet sections must cover every facet exactly
  once whatever the vertical spacing of the cells they are cut by;
* the solver start-up (wall-function reconstruction points, Poisson guard) and
  a short run with the IBM wall functions, the energy balance and the periodic
  energy-balance correction active.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Optional

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
CASE_ID = "064"
CASE_SOURCE = REPO_ROOT / "tests" / "cases" / CASE_ID
CLEANSIM = REPO_ROOT / "tools" / "cleansim.sh"
WRITE_INPUTS = REPO_ROOT / "tools" / "write_inputs.py"
UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0",
)
# decomposition of the solver run; the macOS CI runner has three CPUs and uses 1x2
NPROCX = int(os.environ.get("UDALES_TEST_NPROCX", "2"))
NPROCY = int(os.environ.get("UDALES_TEST_NPROCY", "2"))
DZLIN = 0.5
HLIN = 8.0
# (block, key, value): linear spacing DZLIN up to HLIN, exponential stretch above
NAMELIST_OVERRIDES = (
    ("RUN", "runtime", "20."),
    ("RUN", "nprocx", str(NPROCX)),
    ("RUN", "nprocy", str(NPROCY)),
    ("INP", "lzstretch", ".true."),
    ("INP", "lstretchexp", ".true."),
    ("INP", "stretchconst", "0.01"),
    ("INP", "dzlin", str(DZLIN)),
    ("INP", "hlin", str(HLIN)),
    ("INP", "ishortwave", "1"),
)
# facet_sections_*.txt carry four decimals per section; a facet is cut into up
# to a few hundred sections, so sums agree to a few 1e-3 m2 at best.
AREA_ATOL = 5.0e-3


def _set_namoption(path: Path, block: str, key: str, value: str) -> None:
    """Set ``key = value`` in namelist ``&block``, replacing an existing entry."""
    lines = path.read_text(encoding="ascii").splitlines(keepends=True)
    start = next(
        (i for i, line in enumerate(lines)
         if line.strip().lower().startswith("&") and line.strip()[1:].strip().lower() == block.lower()),
        None,
    )
    if start is None:
        raise ValueError(f"Namelist block &{block} not found in {path}")
    end = next(i for i in range(start + 1, len(lines)) if lines[i].strip().startswith("/"))
    for i in range(start + 1, end):
        if "=" in lines[i] and lines[i].split("=", 1)[0].strip().lower() == key.lower():
            lines[i] = f"{key} = {value}\n"
            break
    else:
        lines.insert(end, f"{key} = {value}\n")
    path.write_text("".join(lines), encoding="ascii")


def _facet_section_sums(case_dir: Path, grid: str, nfcts: int) -> np.ndarray:
    """Total section area per facet from ``facet_sections_<grid>.txt``."""
    data = np.loadtxt(case_dir / f"facet_sections_{grid}.txt", comments="#", ndmin=2)
    sums = np.zeros(nfcts)
    np.add.at(sums, data[:, 0].astype(int) - 1, data[:, 1])
    return sums


# launcher resolution and the one-rank-per-CPU pre-flight check shared by the runners
sys.path.insert(0, str(REPO_ROOT / "tests" / "integration"))
import mpi_env  # noqa: E402

_mpi_exec_and_args = mpi_env.mpi_exec_and_args


def _shell_prefix() -> str:
    module_setup = ""
    if Path("/etc/profile.d/modules.sh").is_file():
        module_setup = "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    return f"{module_setup}if command -v module >/dev/null 2>&1; then module load {RUNTIME_MODULES}; fi && "


def _launcher_unavailable_reason() -> Optional[str]:
    mpiexec, extra_args = _mpi_exec_and_args()
    probe = subprocess.run(
        ["bash", "-lc", f"{_shell_prefix()}'{mpiexec}' {extra_args} -n 1 /bin/true"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if probe.returncode == 0:
        return None
    return (probe.stderr or "").strip() or (probe.stdout or "").strip() or (
        f"mpiexec probe failed with exit code {probe.returncode}"
    )


class TestStretchedGrid(unittest.TestCase):
    """Preprocess case 064 on a stretched grid and run the solver on it."""

    @classmethod
    def setUpClass(cls) -> None:
        if shutil.which("bash") is None:
            raise unittest.SkipTest("bash is required to run tools/cleansim.sh")
        cls.workdir = tempfile.TemporaryDirectory(prefix="udales-stretched-grid-")
        cls.case_dir = Path(cls.workdir.name) / CASE_ID
        shutil.copytree(CASE_SOURCE, cls.case_dir)
        namelist = cls.case_dir / f"namoptions.{CASE_ID}"
        for block, key, value in NAMELIST_OVERRIDES:
            _set_namoption(namelist, block, key, value)

        clean = subprocess.run(
            ["bash", str(CLEANSIM), str(cls.case_dir), "--preprocess-only", "--delete"],
            cwd=REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if clean.returncode != 0:
            raise RuntimeError(f"cleansim failed for the stretched case:\n{clean.stdout}")

        env = dict(os.environ)
        env.setdefault("PREPROC_NCPU", "8")
        result = subprocess.run(
            [sys.executable, str(WRITE_INPUTS), str(cls.case_dir)],
            cwd=REPO_ROOT,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Python preprocessing failed on the stretched grid:\n{result.stdout}")

        cls.facet_area = np.loadtxt(cls.case_dir / f"facetarea.inp.{CASE_ID}", skiprows=1).reshape(-1)
        cls.nfcts = cls.facet_area.size

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "workdir"):
            cls.workdir.cleanup()

    def test_profile_describes_the_requested_stretched_grid(self) -> None:
        zf = np.loadtxt(self.case_dir / f"prof.inp.{CASE_ID}", skiprows=2)[:, 0]
        zh = np.zeros(zf.size + 1)
        for k, centre in enumerate(zf):  # the solver's face recursion
            zh[k + 1] = 2.0 * centre - zh[k]
        dz = np.diff(zh)

        self.assertTrue(np.all(dz > 0.0))
        np.testing.assert_allclose(dz[: int(round(HLIN / DZLIN))], DZLIN, rtol=0, atol=1e-8)
        self.assertGreater(dz.max() / dz.min(), 1.5)
        self.assertTrue(np.all(np.diff(dz) >= -1e-8), "spacing must not decrease with height")
        self.assertAlmostEqual(zh[-1], 64.0, places=6)

    def test_scalar_grid_facet_sections_cover_every_facet_exactly(self) -> None:
        sums = _facet_section_sums(self.case_dir, "c", self.nfcts)
        np.testing.assert_allclose(sums, self.facet_area, rtol=0, atol=AREA_ATOL)

    def test_u_and_v_grid_facet_sections_match_the_uniform_grid(self) -> None:
        # The area of a facet seen by the u or v cells does not depend on how
        # the cells are spaced vertically, so it must equal the committed
        # uniform-grid preprocessing of the same geometry.
        for grid in ("u", "v"):
            with self.subTest(grid=grid):
                stretched = _facet_section_sums(self.case_dir, grid, self.nfcts)
                uniform = _facet_section_sums(CASE_SOURCE, grid, self.nfcts)
                np.testing.assert_allclose(stretched, uniform, rtol=0, atol=AREA_ATOL)

    def test_solver_runs_on_the_stretched_grid(self) -> None:
        if not UDALES_BUILD.is_file():
            self.skipTest(f"u-dales executable not found at {UDALES_BUILD}")
        launcher_issue = _launcher_unavailable_reason()
        if launcher_issue is not None:
            self.skipTest(f"MPI launcher is not usable in this environment: {launcher_issue}")

        try:
            mpi_env.preflight(NPROCX * NPROCY)
        except RuntimeError as exc:
            self.fail(str(exc))
        mpiexec, extra_args = _mpi_exec_and_args()
        command = (
            f"{_shell_prefix()}export HDF5_USE_FILE_LOCKING=FALSE && cd '{self.case_dir}' && "
            f"'{mpiexec}' {extra_args} -n {NPROCX * NPROCY} '{UDALES_BUILD}' namoptions.{CASE_ID}"
        )
        completed = subprocess.run(
            ["bash", "-lc", command], cwd=REPO_ROOT, check=False, capture_output=True, text=True
        )
        log = completed.stdout + completed.stderr
        self.assertEqual(completed.returncode, 0, msg="\n".join(log.splitlines()[-40:]))
        self.assertNotRegex(log, r"(?i)\bnan\b")

        # The pressure solve treats the vertical direction by Gaussian
        # elimination with the true spacings, so the projected velocity must be
        # divergence free to round-off.
        divergences = re.findall(r"divmax, divtot =\s*(\S+)\s+(\S+)", log)
        self.assertTrue(divergences, "no divergence diagnostics found in the solver log")
        self.assertLess(max(abs(float(divmax)) for divmax, _ in divergences), 1.0e-10)

        import netCDF4 as nc

        with nc.Dataset(self.case_dir / f"facEB.{CASE_ID}.nc") as ds:
            heat_flux = np.ma.filled(ds["hf"][:], np.nan)
        self.assertGreater(heat_flux.size, 0)
        self.assertTrue(np.all(np.isfinite(heat_flux)))


if __name__ == "__main__":
    unittest.main()
