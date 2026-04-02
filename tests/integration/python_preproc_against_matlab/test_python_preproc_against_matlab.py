import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
CASE = "100"
CASE_SOURCE = REPO_ROOT / "tests" / "cases" / CASE
EXPECTED_OUTPUTS = [
    f"lscale.inp.{CASE}",
    f"prof.inp.{CASE}",
    f"factypes.inp.{CASE}",
    f"facets.inp.{CASE}",
    f"facetarea.inp.{CASE}",
    "solid_c.txt",
    "solid_u.txt",
    "solid_v.txt",
    "solid_w.txt",
    "fluid_boundary_c.txt",
    "fluid_boundary_u.txt",
    "fluid_boundary_v.txt",
    "fluid_boundary_w.txt",
    "facet_sections_c.txt",
    "facet_sections_u.txt",
    "facet_sections_v.txt",
    "facet_sections_w.txt",
]


def _load_numeric_table(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    return np.atleast_2d(data) if np.ndim(data) else np.asarray([[data]])


class TestPythonPreprocAgainstMatlab(unittest.TestCase):
    """Compare the Python preprocessing chain against MATLAB on no-tree case 100."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory(prefix="udales-preproc-parity-")
        self.addCleanup(self.temp_dir.cleanup)
        temp_root = Path(self.temp_dir.name)

        self.matlab_root = temp_root / "matlab"
        self.python_root = temp_root / "python"
        self.matlab_case = self.matlab_root / CASE
        self.python_case = self.python_root / CASE

        shutil.copytree(CASE_SOURCE, self.matlab_case)
        shutil.copytree(CASE_SOURCE, self.python_case)

        for case_dir in (self.matlab_case, self.python_case):
            for relpath in EXPECTED_OUTPUTS:
                target = case_dir / relpath
                if target.exists():
                    target.unlink()

    def _run_matlab_preprocessing(self) -> None:
        matlab = shutil.which("matlab")
        if matlab is None:
            raise unittest.SkipTest("matlab not found")

        env = os.environ.copy()
        env["DA_EXPDIR"] = str(self.matlab_root)
        env["DA_TOOLSDIR"] = str(REPO_ROOT / "tools")
        env["MATLAB_USE_USERWORK"] = "0"

        cmd = [
            matlab,
            "-nodesktop",
            "-noFigureWindows",
            "-nosplash",
            "-nodisplay",
            "-r",
            f"expnr={CASE}; write_inputs; quit",
        ]
        result = subprocess.run(
            cmd,
            cwd=REPO_ROOT / "tools",
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode != 0:
            self.fail(f"MATLAB preprocessing failed:\n{result.stdout}")

    def _run_python_preprocessing(self) -> None:
        cmd = [sys.executable, str(REPO_ROOT / "tools" / "write_inputs.py"), str(self.python_case)]
        result = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode != 0:
            self.fail(f"Python preprocessing failed:\n{result.stdout}")

    def _assert_same_numeric_file(self, relpath: str) -> None:
        matlab_path = self.matlab_case / relpath
        python_path = self.python_case / relpath

        self.assertTrue(matlab_path.exists(), f"MATLAB output missing: {relpath}")
        self.assertTrue(python_path.exists(), f"Python output missing: {relpath}")

        matlab_data = _load_numeric_table(matlab_path)
        python_data = _load_numeric_table(python_path)

        self.assertEqual(
            matlab_data.shape,
            python_data.shape,
            f"shape mismatch for {relpath}: MATLAB {matlab_data.shape}, Python {python_data.shape}",
        )
        np.testing.assert_allclose(
            python_data,
            matlab_data,
            rtol=0.0,
            atol=1.0e-10,
            err_msg=f"content mismatch for {relpath}",
        )

    def test_case_100_matches_matlab_outputs(self) -> None:
        self._run_matlab_preprocessing()
        self._run_python_preprocessing()

        for relpath in EXPECTED_OUTPUTS:
            with self.subTest(output=relpath):
                self._assert_same_numeric_file(relpath)


if __name__ == "__main__":
    unittest.main()
