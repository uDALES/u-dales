from __future__ import annotations

import json
import os
import sys
import unittest
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from unittest import mock

import numpy as np
from scipy import sparse
import trimesh

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR  # noqa: E402

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udgeom.view3d import (  # noqa: E402
    ViewFactorRepairError,
    ViewFactorRepairLimits,
    ViewFactorValidationError,
    compute_svf,
    count_sparse_entries,
    default_view3d_config_path,
    inspect_view_factors,
    load_view3d_runtime_env,
    read_view3d_output,
    repair_view_factors,
    resolve_view3d_exe,
    stl_to_view3d,
    validate_view_factors,
    write_vf,
    write_view_factor_repair_report,
    write_vfsparse,
)

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

    def test_validate_view_factors_accepts_closed_sparse_matrix(self) -> None:
        vf = sparse.csr_matrix(
            np.array([[0.0, 0.25], [0.4, 0.0]], dtype=float)
        )

        validate_view_factors(vf, compute_svf(vf))

    def test_validate_view_factors_rejects_factor_above_one(self) -> None:
        vf = sparse.csr_matrix(np.array([[0.0, 1.2], [0.0, 0.0]], dtype=float))

        with self.assertRaisesRegex(ValueError, "matrix entries above"):
            validate_view_factors(vf, compute_svf(vf))

    def test_validate_view_factors_rejects_row_closure_error(self) -> None:
        vf = sparse.csr_matrix(np.array([[0.0, 0.8], [0.2, 0.0]], dtype=float))
        svf = np.array([0.3, 0.8], dtype=float)

        with self.assertRaisesRegex(ValueError, "rows fail"):
            validate_view_factors(vf, svf)

    def test_inspect_view_factors_marks_closure_only_failure(self) -> None:
        vf = sparse.csr_matrix(np.array([[0.0, 0.8], [0.2, 0.0]]))
        report = inspect_view_factors(vf, np.array([0.3, 0.8]))

        self.assertFalse(report.is_valid)
        self.assertTrue(report.closure_only)
        self.assertEqual(report.bad_closure_rows, 1)

    def test_repair_view_factors_preserves_reciprocity_and_closure(self) -> None:
        areas = np.array([1.0, 2.0, 4.0])
        exchange = np.array(
            [
                [0.0, 0.8, 0.5],
                [0.8, 0.0, 0.4],
                [0.5, 0.4, 0.0],
            ]
        )
        vf = sparse.csr_matrix(exchange / areas[:, None])
        limits = ViewFactorRepairLimits(
            max_overfull_area_fraction=1.0,
            max_exchange_area_reduction_fraction=1.0,
        )

        repaired, sky, report = repair_view_factors(vf, areas, limits=limits)

        repaired_exchange = sparse.diags(areas) @ repaired
        np.testing.assert_allclose(
            repaired_exchange.toarray(), repaired_exchange.toarray().T, atol=1.0e-12
        )
        np.testing.assert_allclose(
            np.asarray(repaired.sum(axis=1)).reshape(-1) + sky,
            np.ones(3),
            atol=1.0e-12,
        )
        self.assertTrue(report.repaired)
        self.assertEqual(report.overfull_rows, 1)
        self.assertGreater(report.scaled_entries, 0)
        self.assertLessEqual(report.max_row_sum_after, 1.0)
        self.assertAlmostEqual(report.reciprocity_l1_relative_after, 0.0)

    def test_repair_view_factors_rejects_excessive_area_correction(self) -> None:
        vf = sparse.csr_matrix(
            np.array(
                [
                    [0.0, 0.7, 0.6],
                    [0.7, 0.0, 0.6],
                    [0.6, 0.6, 0.0],
                ]
            )
        )

        with self.assertRaisesRegex(ViewFactorRepairError, "facet area"):
            repair_view_factors(vf, np.ones(3))

    def test_repair_view_factors_rejects_individual_factor_above_one(self) -> None:
        vf = sparse.csr_matrix(np.array([[0.0, 1.2], [0.2, 0.0]]))

        with self.assertRaisesRegex(ViewFactorRepairError, "individual view factors"):
            repair_view_factors(vf, np.ones(2))

    def test_view_factor_repair_limits_reject_invalid_values(self) -> None:
        invalid = (
            {"max_row_sum": 1.0},
            {"max_overfull_area_fraction": -0.1},
            {"max_exchange_area_reduction_fraction": 1.1},
            {"max_reciprocity_l1_relative": 0.0},
            {"max_row_sum": np.inf},
        )

        for kwargs in invalid:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    ViewFactorRepairLimits(**kwargs)

    def test_repair_view_factors_does_not_treat_roundoff_as_material(self) -> None:
        vf = sparse.csr_matrix(np.array([[0.0, 1.0005], [1.0005, 0.0]]))

        unchanged, sky, report = repair_view_factors(vf, np.ones(2))

        self.assertFalse(report.repaired)
        self.assertEqual(report.overfull_rows, 2)
        self.assertEqual(report.materially_overfull_rows, 0)
        self.assertEqual(report.materially_overfull_area_fraction, 0.0)
        np.testing.assert_allclose(unchanged.toarray(), vf.toarray())
        validate_view_factors(unchanged, sky, areas=np.ones(2))

    def test_repair_view_factors_honours_custom_reciprocity_limit(self) -> None:
        vf = sparse.csr_matrix(np.array([[0.0, 0.5], [0.4, 0.0]]))
        limits = ViewFactorRepairLimits(max_reciprocity_l1_relative=0.5)

        unchanged, sky, report = repair_view_factors(
            vf, np.ones(2), limits=limits
        )

        self.assertFalse(report.repaired)
        np.testing.assert_allclose(unchanged.toarray(), vf.toarray())
        validate_view_factors(
            unchanged, sky, areas=np.ones(2), reciprocity_tolerance=0.5
        )

    def test_write_view_factor_repair_report_is_machine_readable(self) -> None:
        areas = np.ones(3)
        vf = sparse.csr_matrix(
            np.array(
                ([0.0, 0.7, 0.6], [0.7, 0.0, 0.2], [0.6, 0.2, 0.0])
            )
        )
        limits = ViewFactorRepairLimits(
            max_overfull_area_fraction=1.0,
            max_exchange_area_reduction_fraction=1.0,
        )
        _, _, report = repair_view_factors(vf, areas, limits=limits)
        out_path = self.workdir / "view3d_repair.101.json"

        write_view_factor_repair_report(out_path, report)
        payload = json.loads(out_path.read_text(encoding="ascii"))

        self.assertEqual(payload["algorithm"], "reciprocal-open-domain-v1")
        self.assertEqual(payload["materially_overfull_rows"], 1)
        self.assertGreater(payload["scaled_entries"], 0)

    def test_strict_validation_rejects_material_reciprocity_error(self) -> None:
        vf = sparse.csr_matrix(np.array([[0.0, 0.8], [0.2, 0.0]]))
        sky = compute_svf(vf)

        with self.assertRaisesRegex(ViewFactorValidationError, "reciprocity"):
            validate_view_factors(vf, sky, areas=np.ones(2))

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
        self.assertEqual(
            lines,
            ["1 2 0.50000000", "2 1 0.20000000", "3 2 0.40000000"],
        )

    def test_repaired_sparse_round_trip_remains_physical(self) -> None:
        areas = np.array([1.0, 2.0, 4.0])
        exchange = np.array(
            [
                [0.0, 0.8, 0.5],
                [0.8, 0.0, 0.4],
                [0.5, 0.4, 0.0],
            ]
        )
        vf = sparse.csr_matrix(exchange / areas[:, None])
        limits = ViewFactorRepairLimits(
            max_overfull_area_fraction=1.0,
            max_exchange_area_reduction_fraction=1.0,
        )
        repaired, _, _ = repair_view_factors(vf, areas, limits=limits)
        out_path = self.workdir / "vfsparse.inp.101"

        write_vfsparse(out_path, repaired, threshold=0.0)
        reread = read_view3d_output(out_path, nfacets=3, outformat=2)
        sky = compute_svf(reread)

        validate_view_factors(reread, sky, areas=areas)
        np.testing.assert_allclose(
            np.asarray(reread.sum(axis=1)).reshape(-1) + sky,
            np.ones(3),
            atol=1.0e-7,
        )

    def test_count_sparse_entries_ignores_blank_lines(self) -> None:
        out_path = self.workdir / "vfsparse.inp.101"
        out_path.write_text(
            "1 2 0.500000\n\n  \n2 1 0.200000\n",
            encoding="ascii",
        )

        self.assertEqual(count_sparse_entries(out_path), 2)

    def test_write_vf_uses_matlab_compatible_orientation(self) -> None:
        try:
            import netCDF4
        except ImportError:
            self.skipTest("netCDF4 not available")

        vf = np.array([[0.0, 0.25, 0.5], [0.75, 0.0, 1.0], [0.125, 0.375, 0.0]], dtype=float)
        out_path = self.workdir / "vf.nc.inp.101"

        write_vf(out_path, vf)

        with netCDF4.Dataset(out_path) as ds:
            stored = np.asarray(ds.variables["view factor"][:])
        np.testing.assert_allclose(stored, vf.T)

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

    def test_read_view3d_output_sparse_text_accepts_single_zero_based_entry(self) -> None:
        out_path = self.workdir / "vf_sparse_zero_based.txt"
        out_path.write_text("0 2 0.75\n", encoding="ascii")

        vf = read_view3d_output(out_path, nfacets=3, outformat=2, one_based=False)

        expected = np.array(
            [
                [0.0, 0.0, 0.75],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ],
            dtype=float,
        )
        np.testing.assert_allclose(vf.toarray(), expected)

    def test_read_view3d_binary_matches_matlab_reshape_transpose(self) -> None:
        out_path = self.workdir / "vf.bin"
        expected = np.array(
            [
                [0.0, 0.12, 0.34],
                [0.56, 0.0, 0.78],
                [0.91, 0.23, 0.0],
            ],
            dtype=np.float32,
        )
        header_and_area = np.zeros(8 + expected.shape[0], dtype=np.float32)
        with out_path.open("wb") as f:
            header_and_area.tofile(f)
            expected.ravel(order="C").tofile(f)

        vf = read_view3d_output(out_path, nfacets=expected.shape[0], outformat=1)

        np.testing.assert_allclose(vf.toarray(), expected)
        self.assertEqual(vf.dtype, np.dtype(float))

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

    @unittest.skipIf(shutil.which("bash") is None, "bash is required for View3D shell config sourcing")
    def test_load_view3d_runtime_env_sources_shell_config(self) -> None:
        config = self.workdir / "view3d_config.sh"
        config.write_text(
            "\n".join(
                [
                    "export VIEW3D_NUM_THREADS=7",
                    "OMP_NUM_THREADS=7",
                    'export VIEW3D_MAX_DENSE_MATRIX_GIB="${VIEW3D_MAX_DENSE_MATRIX_GIB:-112}"',
                    "echo ignored stdout",
                ]
            )
            + "\n",
            encoding="ascii",
        )

        env, path = load_view3d_runtime_env(
            base_env={"PATH": os.environ.get("PATH", ""), "VIEW3D_MAX_DENSE_MATRIX_GIB": "128"},
            config_path=config,
        )

        self.assertEqual(path, config.resolve())
        self.assertEqual(env["VIEW3D_NUM_THREADS"], "7")
        self.assertEqual(env["OMP_NUM_THREADS"], "7")
        self.assertEqual(env["VIEW3D_MAX_DENSE_MATRIX_GIB"], "128")

    @unittest.skipIf(shutil.which("bash") is None, "bash is required for View3D shell config sourcing")
    def test_default_view3d_config_derives_dense_limit_from_preproc_mem(self) -> None:
        env, path = load_view3d_runtime_env(
            base_env={"PATH": os.environ.get("PATH", ""), "PREPROC_NCPU": "8", "PREPROC_MEM": "128gb"},
            config_path=default_view3d_config_path(),
        )

        self.assertEqual(path, default_view3d_config_path().resolve())
        self.assertEqual(env["VIEW3D_MAX_DENSE_MATRIX_GIB"], "112")

    @unittest.skipIf(shutil.which("bash") is None, "bash is required for View3D shell config sourcing")
    def test_default_view3d_config_uses_low_preproc_mem_as_dense_limit(self) -> None:
        env, _ = load_view3d_runtime_env(
            base_env={"PATH": os.environ.get("PATH", ""), "PREPROC_NCPU": "8", "PREPROC_MEM": "8gb"},
            config_path=default_view3d_config_path(),
        )

        self.assertEqual(env["VIEW3D_MAX_DENSE_MATRIX_GIB"], "8")

    @unittest.skipIf(shutil.which("bash") is None, "bash is required for View3D shell config sourcing")
    def test_default_view3d_config_rejects_invalid_preproc_mem(self) -> None:
        for value in ("128", "abcgb"):
            with self.subTest(PREPROC_MEM=value):
                with self.assertRaisesRegex(RuntimeError, "PREPROC_MEM must be set like 128gb"):
                    load_view3d_runtime_env(
                        base_env={"PATH": os.environ.get("PATH", ""), "PREPROC_NCPU": "8", "PREPROC_MEM": value},
                        config_path=default_view3d_config_path(),
                    )

    @unittest.skipIf(shutil.which("bash") is None, "bash is required for View3D shell config sourcing")
    def test_default_view3d_config_uses_preproc_ncpu_when_available(self) -> None:
        env, _ = load_view3d_runtime_env(
            base_env={"PATH": os.environ.get("PATH", ""), "PREPROC_NCPU": "6"},
            config_path=default_view3d_config_path(),
        )

        self.assertEqual(env["VIEW3D_NUM_THREADS"], "6")
        self.assertEqual(env["OMP_NUM_THREADS"], "6")

    @unittest.skipIf(shutil.which("bash") is None, "bash is required for View3D shell config sourcing")
    def test_default_view3d_config_requires_preproc_ncpu(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "PREPROC_NCPU must be set"):
            load_view3d_runtime_env(
                base_env={"PATH": os.environ.get("PATH", "")},
                config_path=default_view3d_config_path(),
            )

    @unittest.skipIf(shutil.which("bash") is None, "bash is required for View3D shell config sourcing")
    def test_default_view3d_config_preserves_explicit_dense_limit(self) -> None:
        env, _ = load_view3d_runtime_env(
            base_env={
                "PATH": os.environ.get("PATH", ""),
                "PREPROC_NCPU": "8",
                "PREPROC_MEM": "128gb",
                "VIEW3D_MAX_DENSE_MATRIX_GIB": "96",
            },
            config_path=default_view3d_config_path(),
        )

        self.assertEqual(env["VIEW3D_MAX_DENSE_MATRIX_GIB"], "96")

    @unittest.skipIf(shutil.which("bash") is None, "bash is required for View3D shell config sourcing")
    def test_resolve_view3d_exe_uses_configured_executable(self) -> None:
        config = self.workdir / "view3d_config.sh"
        custom_exe = self.workdir / "custom-view3d"
        config.write_text(f'export VIEW3D_EXE="{custom_exe}"\n', encoding="ascii")

        with mock.patch.dict(os.environ, {"VIEW3D_CONFIG": str(config)}):
            self.assertEqual(resolve_view3d_exe(), custom_exe.resolve())


if __name__ == "__main__":
    unittest.main()
