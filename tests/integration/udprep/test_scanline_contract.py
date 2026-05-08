from __future__ import annotations

import importlib
import tempfile
import unittest
from pathlib import Path

import numpy as np

from tools.python.tests._common import REPO_ROOT
from udbase import UDBase

_IMPORT_ERROR = None
try:
    from udprep.directshortwave import DirectShortwaveSolver
except ImportError as exc:
    DirectShortwaveSolver = None
    _IMPORT_ERROR = exc


def _missing_prereqs():
    missing = []
    if _IMPORT_ERROR is not None:
        missing.append(f"python directshortwave import failed: {_IMPORT_ERROR}")
    elif importlib.util.find_spec("udprep.directshortwave_f2py") is None:
        missing.append("compiled directshortwave_f2py wrapper not found")
    return missing


class TestScanlineContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        missing = _missing_prereqs()
        if missing:
            raise unittest.SkipTest("; ".join(missing))
        cls.case_dir = REPO_ROOT / "tests" / "cases" / "064"
        cls.sim = UDBase("064", cls.case_dir, load_geometry=True)
        cls.solver = DirectShortwaveSolver(cls.sim, method="scanline", ray_density=12.0, ray_jitter=0.0)
        cls.nsun = np.array([0.0, 0.0, 1.0], dtype=float)

    def test_scanline_inputs_follow_matlab_contract(self) -> None:
        inputs = self.solver._build_scanline_inputs(self.nsun, 800.0, resolution=0.1)
        mesh = self.sim.geom.stl

        np.testing.assert_array_equal(inputs.faces_1based, np.asarray(mesh.faces, dtype=np.int32) + 1)
        np.testing.assert_allclose(inputs.facet_points, np.asarray(self.sim.geom.face_incenters, dtype=float))
        np.testing.assert_allclose(inputs.face_normals, np.asarray(mesh.face_normals, dtype=float))
        np.testing.assert_allclose(inputs.vertices, np.asarray(mesh.vertices, dtype=float))
        np.testing.assert_allclose(inputs.nsun, self.nsun)
        self.assertEqual(inputs.irradiance, 800.0)
        self.assertEqual(inputs.resolution, 0.1)

    def test_scanline_inputs_are_fortran_ready(self) -> None:
        inputs = self.solver._build_scanline_inputs(self.nsun, 800.0, resolution=0.1)
        self.assertTrue(inputs.faces_1based.flags["F_CONTIGUOUS"])
        self.assertTrue(inputs.facet_points.flags["F_CONTIGUOUS"])
        self.assertTrue(inputs.face_normals.flags["F_CONTIGUOUS"])
        self.assertTrue(inputs.vertices.flags["F_CONTIGUOUS"])
        self.assertTrue(inputs.nsun.flags["F_CONTIGUOUS"])
        self.assertEqual(inputs.faces_1based.dtype, np.int32)

    def test_legacy_serialization_matches_canonical_inputs(self) -> None:
        inputs = self.solver._build_scanline_inputs(self.nsun, 800.0, resolution=0.1)
        with tempfile.TemporaryDirectory(prefix="udales-scanline-contract-") as td:
            td_path = Path(td)
            self.solver._write_scanline_legacy_inputs(td_path, inputs)

            vertices = np.loadtxt(td_path / "vertices.txt")
            faces = np.loadtxt(td_path / "faces.txt")
            with (td_path / "info_directShortwave.txt").open("r", encoding="ascii") as f:
                info_lines = tuple(f.readlines())

        np.testing.assert_allclose(vertices, inputs.vertices)
        np.testing.assert_allclose(faces, self.solver._scanline_face_rows(inputs))
        self.assertEqual(info_lines, self.solver._scanline_info_lines(inputs))


if __name__ == "__main__":
    unittest.main()
