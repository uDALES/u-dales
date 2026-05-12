from __future__ import annotations

import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from scipy import sparse
import trimesh

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR  # noqa: E402

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udgeom.view3d import compute_svf, read_view3d_output, stl_to_view3d, write_vfsparse  # noqa: E402


class TestView3DUtils(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)

    def test_compute_svf_clamps_negative_values(self) -> None:
        vf = sparse.csr_matrix(
            np.array(
                [
                    [0.1, 0.2, 0.3],
                    [0.6, 0.5, 0.0],
                    [0.0, 0.0, 0.0],
                ],
                dtype=float,
            )
        )

        svf = compute_svf(vf)

        np.testing.assert_allclose(svf, np.array([0.4, 0.0, 1.0], dtype=float))

    def test_write_vfsparse_filters_and_sorts_entries(self) -> None:
        vf = sparse.coo_matrix(
            (
                np.array([0.4, 1.0e-8, 0.2, 0.5], dtype=float),
                (
                    np.array([2, 0, 1, 0], dtype=int),
                    np.array([1, 2, 0, 1], dtype=int),
                ),
            ),
            shape=(3, 3),
        )
        out_path = self.workdir / "vfsparse.inp.101"

        write_vfsparse(out_path, vf, threshold=1.0e-6)

        lines = out_path.read_text(encoding="ascii").splitlines()
        self.assertEqual(lines, ["1 2 0.500000", "2 1 0.200000", "3 2 0.400000"])

    def test_read_view3d_output_sparse_text_converts_one_based_indices(self) -> None:
        out_path = self.workdir / "vf_sparse.txt"
        out_path.write_text("1 2 0.5\n3 1 0.25\n", encoding="ascii")

        vf = read_view3d_output(out_path, nfacets=3, outformat=2, one_based=True)

        expected = np.array(
            [
                [0.0, 0.5, 0.0],
                [0.0, 0.0, 0.0],
                [0.25, 0.0, 0.0],
            ],
            dtype=float,
        )
        np.testing.assert_allclose(vf.toarray(), expected)

    def test_stl_to_view3d_writes_expected_header_and_surface_records(self) -> None:
        mesh = trimesh.Trimesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ],
                dtype=float,
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
            process=False,
        )
        stl_path = self.workdir / "triangle.stl"
        vs3_path = self.workdir / "triangle.vs3"
        mesh.export(stl_path)

        result = stl_to_view3d(stl_path, vs3_path, outformat=2, maxD=50.0, row=3, col=4)

        self.assertEqual(result, vs3_path)
        content = vs3_path.read_bytes().decode("ascii")
        self.assertIn("T\r\n", content)
        self.assertIn("C out=2 maxD=50.0 row=3 col=4\r\n", content)
        self.assertIn("V    1 0.000000 0.000000 0.000000\r\n", content)
        self.assertIn("S    1      1      2      3      0      0      0      0      1f\r\n", content)
        self.assertTrue(content.endswith("End of Data\r\n"))


if __name__ == "__main__":
    unittest.main()
