from __future__ import annotations

import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
PYTHON_DIR = REPO_ROOT / "tools" / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udbase import UDBase
from udprep.directshortwave import DirectShortwaveSolver


class _DirectShortwaveCaseMixin:
    RAY_DENSITY = 12.0
    AZIMUTH_DEG = 20.0
    ELEVATION_DEG = 15.0
    IRRADIANCE = 800.0

    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = REPO_ROOT
        cls.case_dir = cls.repo_root / "tests" / "cases" / cls.CASE
        cls.sim = UDBase(cls.CASE, cls.case_dir, load_geometry=True)
        cls.areas = np.asarray(cls.sim.geom.stl.area_faces, dtype=float)

        az = np.deg2rad(cls.AZIMUTH_DEG)
        el = np.deg2rad(cls.ELEVATION_DEG)
        cls.nsun = np.array(
            [np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)],
            dtype=float,
        )

    def _cell_min(self) -> float:
        return min(self.sim.dx, self.sim.dy, float(np.min(self.sim.dzt)))

    def _scanline_spacing(self, ray_density: float) -> float:
        return 0.25 * self._cell_min() / ray_density

    def _ray_spacing(self, ray_density: float) -> float:
        return self._cell_min() / ray_density

    def _run_python_method(
        self,
        method: str,
        ray_density: float,
        *,
        veg_data=None,
    ) -> tuple[np.ndarray, np.ndarray, dict[str, float], float]:
        solver = DirectShortwaveSolver(
            self.sim,
            method=method,
            ray_density=ray_density,
            ray_jitter=0.0,
            veg_data=veg_data,
        )
        t0 = time.perf_counter()
        sdir, veg, budget = solver.compute(
            nsun=self.nsun,
            irradiance=self.IRRADIANCE,
            periodic_xy=False,
        )
        elapsed = time.perf_counter() - t0
        return sdir, veg, budget, elapsed

    def _compare_fields(self, candidate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
        diff = candidate - reference
        common = (reference > 1.0e-6) & (candidate > 1.0e-6)
        rel = np.abs(diff[common]) / np.maximum(reference[common], 1.0)
        return {
            "energy_ratio": float(np.sum(candidate * self.areas) / np.sum(reference * self.areas)),
            "mean_abs": float(np.mean(np.abs(diff))),
            "p95_rel": float(np.quantile(rel, 0.95)),
            "corr": float(np.corrcoef(reference[common], candidate[common])[0, 1]),
        }

    def _assert_budget_conserved(self, budget: dict[str, float], *, rel_tol: float) -> None:
        self.assertIn("in", budget)
        self.assertIn("fac", budget)
        self.assertIn("veg", budget)
        self.assertIn("out", budget)
        closure = budget["fac"] + budget["veg"] + budget["out"]
        delta = rel_tol * max(abs(budget["in"]), 1.0)
        self.assertAlmostEqual(budget["in"], closure, delta=delta)


class TestDirectShortwaveReferenceIntegration(_DirectShortwaveCaseMixin, unittest.TestCase):
    """Reference integration test on the no-tree case `100`."""

    CASE = "100"

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.build_dir = tempfile.TemporaryDirectory(prefix="udales-directshortwave-")
        cls.fortran_exe = Path(cls.build_dir.name) / "DS.exe"
        subprocess.run(
            ["gfortran", "-O3", str(cls.repo_root / "tools" / "SEB" / "directShortwave.f90"), "-o", str(cls.fortran_exe)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "fortran_exe") and cls.fortran_exe.exists():
            cls.fortran_exe.unlink()
        if hasattr(cls, "build_dir"):
            cls.build_dir.cleanup()

    def _run_fortran_scanline(self, ray_density: float) -> tuple[np.ndarray, float]:
        mesh = self.sim.geom.stl
        resolution = self._scanline_spacing(ray_density)
        faces = np.asarray(mesh.faces, dtype=np.int32) + 1
        incenter = np.asarray(mesh.triangles_center, dtype=float)
        face_normal = np.asarray(mesh.face_normals, dtype=float)
        vertices = np.asarray(mesh.vertices, dtype=float)

        with tempfile.TemporaryDirectory(prefix="udales-directshortwave-case-") as td:
            td_path = Path(td)
            np.savetxt(td_path / "vertices.txt", vertices, fmt="%15.10f %15.10f %15.10f")
            face_rows = np.hstack([faces, incenter, face_normal])
            np.savetxt(
                td_path / "faces.txt",
                face_rows,
                fmt="%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f",
            )
            with (td_path / "info_directShortwave.txt").open("w", encoding="ascii") as f:
                f.write(f"{len(faces):8d} {len(vertices):8d}\n")
                f.write(f"{self.nsun[0]:15.10f} {self.nsun[1]:15.10f} {self.nsun[2]:15.10f}\n")
                f.write(f"{self.IRRADIANCE:15.10f}\n")
                f.write(f"{resolution:15.10f}\n")

            t0 = time.perf_counter()
            subprocess.run(
                [str(self.fortran_exe)],
                cwd=td_path,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            elapsed = time.perf_counter() - t0
            return np.loadtxt(td_path / "Sdir.txt"), elapsed

    def test_scanline_matches_legacy_fortran_reference(self) -> None:
        reference, fortran_seconds = self._run_fortran_scanline(self.RAY_DENSITY)
        scanline, _, budget, f2py_seconds = self._run_python_method("scanline", self.RAY_DENSITY)
        metrics = self._compare_fields(scanline, reference)
        budget_energy_ratio = float(np.sum(reference * self.areas) / budget["fac"])

        print(
            "\nscanline benchmark:",
            {
                "case": self.CASE,
                "ray_density": self.RAY_DENSITY,
                "scanline_spacing": round(self._scanline_spacing(self.RAY_DENSITY), 6),
                "fortran_seconds": round(fortran_seconds, 3),
                "f2py_seconds": round(f2py_seconds, 3),
                "speedup_f2py_vs_fortran": round(fortran_seconds / f2py_seconds, 3),
            },
        )

        self.assertLess(fortran_seconds, 10.0)
        self.assertLess(f2py_seconds, 10.0)
        self.assertGreater(metrics["corr"], 0.999999)
        self.assertLess(metrics["mean_abs"], 0.01)
        self.assertLess(metrics["p95_rel"], 1.0e-3)
        self.assertAlmostEqual(metrics["energy_ratio"], 1.0, delta=5.0e-4)
        self.assertAlmostEqual(budget_energy_ratio, 1.0, delta=5.0e-4)

    def test_python_methods_track_legacy_reference(self) -> None:
        reference, _ = self._run_fortran_scanline(self.RAY_DENSITY)

        method_metrics: dict[str, dict[str, float]] = {}
        for method in ("facsec", "moller"):
            candidate, _, budget, elapsed = self._run_python_method(method, self.RAY_DENSITY)
            metrics = self._compare_fields(candidate, reference)
            metrics["seconds"] = elapsed
            metrics["budget_energy_ratio"] = float(np.sum(reference * self.areas) / budget["fac"])
            method_metrics[method] = metrics
            self._assert_budget_conserved(budget, rel_tol=5.0e-3)

        print(
            "\npython method benchmark:",
            {
                method: {
                    "ray_spacing": round(self._ray_spacing(self.RAY_DENSITY), 6),
                    "seconds": round(metrics["seconds"], 3),
                    "energy_ratio": round(metrics["energy_ratio"], 6),
                    "corr": round(metrics["corr"], 6),
                    "p95_rel": round(metrics["p95_rel"], 6),
                }
                for method, metrics in method_metrics.items()
            },
        )

        facsec = method_metrics["facsec"]
        moller = method_metrics["moller"]

        self.assertLess(facsec["seconds"], 120.0)
        self.assertLess(moller["seconds"], 30.0)

        self.assertGreater(facsec["corr"], 0.99)
        self.assertLess(facsec["mean_abs"], 10.0)
        self.assertLess(facsec["p95_rel"], 0.6)
        self.assertAlmostEqual(facsec["energy_ratio"], 1.0, delta=5.0e-3)
        self.assertAlmostEqual(facsec["budget_energy_ratio"], 1.0, delta=5.0e-3)

        self.assertGreater(moller["corr"], 0.999)
        self.assertLess(moller["mean_abs"], 1.0)
        self.assertLess(moller["p95_rel"], 0.1)
        self.assertAlmostEqual(moller["energy_ratio"], 1.0, delta=2.0e-3)
        self.assertAlmostEqual(moller["budget_energy_ratio"], 1.0, delta=2.0e-3)

        self.assertLess(moller["mean_abs"], facsec["mean_abs"])
        self.assertLess(moller["p95_rel"], facsec["p95_rel"])
        self.assertGreater(moller["corr"], facsec["corr"])


class TestDirectShortwaveTreesIntegration(_DirectShortwaveCaseMixin, unittest.TestCase):
    """Tree-enabled integration test on flat terrain case `525`."""

    CASE = "525"

    def test_trees_reduce_direct_flux_and_absorb_energy(self) -> None:
        for method in ("facsec", "moller"):
            sdir_clear, _, bud_clear, _ = self._run_python_method(method, self.RAY_DENSITY, veg_data=None)
            sdir_trees, veg_trees, bud_trees, _ = self._run_python_method(method, self.RAY_DENSITY, veg_data=self.sim.veg)

            self.assertTrue(np.all(np.isfinite(sdir_clear)))
            self.assertTrue(np.all(np.isfinite(sdir_trees)))
            self.assertGreater(np.count_nonzero(sdir_clear > 0.0), 0)
            self.assertGreater(np.count_nonzero(sdir_trees > 0.0), 0)
            self.assertGreater(len(veg_trees), 0)
            self.assertGreater(bud_trees["veg"], 0.0)
            self.assertLess(bud_trees["fac"], bud_clear["fac"])
            self._assert_budget_conserved(bud_clear, rel_tol=5.0e-3)
            self._assert_budget_conserved(bud_trees, rel_tol=5.0e-3)

    def test_tree_methods_agree_on_flat_terrain_case(self) -> None:
        facsec, _, bud_facsec, facsec_seconds = self._run_python_method(
            "facsec",
            self.RAY_DENSITY,
            veg_data=self.sim.veg,
        )
        moller, _, bud_moller, moller_seconds = self._run_python_method(
            "moller",
            self.RAY_DENSITY,
            veg_data=self.sim.veg,
        )
        metrics = self._compare_fields(facsec, moller)

        print(
            "\ntree benchmark:",
            {
                "case": self.CASE,
                "ray_density": self.RAY_DENSITY,
                "ray_spacing": round(self._ray_spacing(self.RAY_DENSITY), 6),
                "facsec_seconds": round(facsec_seconds, 3),
                "moller_seconds": round(moller_seconds, 3),
                "facsec_fac": round(bud_facsec["fac"], 3),
                "moller_fac": round(bud_moller["fac"], 3),
                "veg_absorption": round(bud_facsec["veg"], 3),
                "corr": round(metrics["corr"], 6),
                "p95_rel": round(metrics["p95_rel"], 6),
            },
        )

        self.assertAlmostEqual(bud_facsec["veg"], bud_moller["veg"], delta=1.0e-9)
        self.assertGreater(metrics["corr"], 0.999)
        self.assertLess(metrics["mean_abs"], 1.0)
        self.assertLess(metrics["p95_rel"], 0.01)
        self.assertAlmostEqual(metrics["energy_ratio"], 1.0, delta=5.0e-3)


if __name__ == "__main__":
    unittest.main()
