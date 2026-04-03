import shutil
import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
CASE_SOURCE = REPO_ROOT / "tests" / "cases" / "100"
COMPARE_TOOL = REPO_ROOT / "tools" / "compare_preprocessing.py"


class TestPythonPreprocAgainstMatlab(unittest.TestCase):
    """Compare the Python preprocessing chain against MATLAB on no-tree case 100."""

    def test_case_100_matches_matlab_outputs(self) -> None:
        matlab = shutil.which("matlab")
        if matlab is None:
            raise unittest.SkipTest("matlab not found")

        cmd = [sys.executable, str(COMPARE_TOOL), "--case-dir", str(CASE_SOURCE)]
        result = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode != 0:
            self.fail(f"Preprocessing comparison failed:\n{result.stdout}")


if __name__ == "__main__":
    unittest.main()
