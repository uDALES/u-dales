import json
import sys
import unittest
from pathlib import Path
from typing import Dict

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
PYTHON_DIR = REPO_ROOT / "tools" / "python"

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udbase import UDBase


DATA_DIR = REPO_ROOT / "tests" / "integration" / "udbase_against_matlab" / "data"


def _load_reference(case: str) -> Dict:
    return json.loads((DATA_DIR / f"{case}.json").read_text())


class TestUDBaseAgainstMatlab(unittest.TestCase):
    cases = ("064", "101")

    def test_reference_files_exist(self):
        self.assertTrue((DATA_DIR / "manifest.json").exists())
        for case in self.cases:
            self.assertTrue((DATA_DIR / f"{case}.json").exists(), case)

    def test_load_facsec_c_matches_matlab_after_index_normalization(self):
        for case in self.cases:
            with self.subTest(case=case):
                sim = UDBase(case, REPO_ROOT / "tests" / "cases" / case, load_geometry=True)
                ref = _load_reference(case)["facsec_c"]
                facsec = sim.load_facsec("c")

                np.testing.assert_array_equal(facsec["facid"], np.asarray(ref["facid"], dtype=int) - 1)
                np.testing.assert_allclose(facsec["area"], np.asarray(ref["area"], dtype=float), atol=1e-12)
                np.testing.assert_array_equal(facsec["locs"], np.asarray(ref["locs"], dtype=int) - 1)
                np.testing.assert_allclose(
                    facsec["distance"],
                    np.asarray(ref["distance"], dtype=float),
                    atol=1e-12,
                )

    def test_calculate_frontal_properties_matches_matlab(self):
        for case in self.cases:
            with self.subTest(case=case):
                sim = UDBase(case, REPO_ROOT / "tests" / "cases" / case, load_geometry=True)
                ref = _load_reference(case)["frontal"]
                frontal = sim.calculate_frontal_properties()

                np.testing.assert_allclose(
                    frontal["skylinex"],
                    np.asarray(ref["skylinex"], dtype=float),
                    atol=1e-12,
                )
                np.testing.assert_allclose(
                    frontal["skyliney"],
                    np.asarray(ref["skyliney"], dtype=float),
                    atol=1e-12,
                )
                self.assertAlmostEqual(frontal["Afx"], float(ref["Afx"]), places=12)
                self.assertAlmostEqual(frontal["Afy"], float(ref["Afy"]), places=12)
                self.assertAlmostEqual(frontal["brx"], float(ref["brx"]), places=12)
                self.assertAlmostEqual(frontal["bry"], float(ref["bry"]), places=12)


if __name__ == "__main__":
    unittest.main()
