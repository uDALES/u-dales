from __future__ import annotations

import subprocess
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np

try:
    from _common import REPO_ROOT
except ModuleNotFoundError:
    from tools.python.tests._common import REPO_ROOT
from udbase import UDBase
from udprep.directshortwave import DirectShortwaveSolver


class TestDirectShortwaveIntegration(unittest.TestCase):
    CASE = "100"
    RAY_DENSITY = 12.0

    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = REPO_ROOT
        cls.case_dir = cls.repo_root.parent / "experiments" / cls.CASE
        cls.sim = UDBase(cls.CASE, cls.case_dir, load_geometry=True)
        cls.areas = np.asarray(cls.sim.geom.stl.area_faces, dtype=float)

        azimuth_deg = 20.0
        elevation_deg = 15.0
        az = np.deg2rad(azimuth_deg)
        el = np.deg2rad(elevation_deg)
        cls.nsun = np.array(
            [np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)],
            dtype=float,
        )
        cls.irradiance = 800.0

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
        resolution = 0.25 * min(self.sim.dx, self.sim.dy, float(np.min(self.sim.dzt))) / ray_density
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
                f.write(f"{self.irradiance:15.10f}\n")
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

    def _run_python_method(self, method: str, ray_density: float) -> tuple[np.ndarray, dict[str, float], float]:
        solver = DirectShortwaveSolver(
            self.sim,
            method=method,
            ray_density=ray_density,
            ray_jitter=0.0,
        )
        t0 = time.perf_counter()
        sdir, _, budget = solver.compute(
            nsun=self.nsun,
            irradiance=self.irradiance,
            periodic_xy=False,
        )
        elapsed = time.perf_counter() - t0
        return sdir, budget, elapsed

    def _compare_to_reference(
        self,
        candidate: np.ndarray,
        *,
        reference: np.ndarray,
        budget: dict[str, float],
        elapsed_seconds: float,
    ) -> dict[str, float]:
        diff = candidate - reference
        common = (reference > 1.0e-6) & (candidate > 1.0e-6)
        rel = np.abs(diff[common]) / np.maximum(reference[common], 1.0)
        return {
            "seconds": elapsed_seconds,
            "energy_ratio": float(np.sum(candidate * self.areas) / np.sum(reference * self.areas)),
            "budget_energy_ratio": float(np.sum(reference * self.areas) / budget["fac"]),
            "max_abs": float(np.max(np.abs(diff))),
            "mean_abs": float(np.mean(np.abs(diff))),
            "median_rel": float(np.median(rel)),
            "p95_rel": float(np.quantile(rel, 0.95)),
            "corr": float(np.corrcoef(reference[common], candidate[common])[0, 1]),
            "common_lit": int(np.count_nonzero(common)),
            "lit": int(np.count_nonzero(candidate > 1.0e-6)),
        }

    def test_scanline_matches_legacy_fortran_reference(self) -> None:
        reference, fortran_seconds = self._run_fortran_scanline(self.RAY_DENSITY)
        scanline, budget, f2py_seconds = self._run_python_method("scanline", self.RAY_DENSITY)
        metrics = self._compare_to_reference(
            scanline,
            reference=reference,
            budget=budget,
            elapsed_seconds=f2py_seconds,
        )

        print(
            "\nscanline benchmark:",
            {
                "case": self.CASE,
                "ray_density": self.RAY_DENSITY,
                "fortran_seconds": round(fortran_seconds, 3),
                "f2py_seconds": round(f2py_seconds, 3),
                "speedup_f2py_vs_fortran": round(fortran_seconds / f2py_seconds, 3),
            },
        )

        self.assertLess(fortran_seconds, 10.0)
        self.assertLess(f2py_seconds, 10.0)
        self.assertGreater(metrics["corr"], 0.999999)
        self.assertLess(metrics["max_abs"], 0.01)
        self.assertLess(metrics["mean_abs"], 0.01)
        self.assertLess(metrics["p95_rel"], 1.0e-3)
        self.assertAlmostEqual(metrics["energy_ratio"], 1.0, delta=5.0e-4)
        self.assertAlmostEqual(metrics["budget_energy_ratio"], 1.0, delta=5.0e-4)

    def test_python_methods_track_legacy_reference(self) -> None:
        reference, _ = self._run_fortran_scanline(self.RAY_DENSITY)

        method_metrics: dict[str, dict[str, float]] = {}
        for method in ("facsec", "moller"):
            candidate, budget, elapsed = self._run_python_method(method, self.RAY_DENSITY)
            method_metrics[method] = self._compare_to_reference(
                candidate,
                reference=reference,
                budget=budget,
                elapsed_seconds=elapsed,
            )

        print(
            "\npython method benchmark:",
            {
                method: {
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


if __name__ == "__main__":
    unittest.main()
