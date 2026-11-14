"""Exercise the water-facet (walltype -21) energy-balance behaviour.

Self-contained on examples/201: ground facets west of x = 64 m are retyped to
the water band, a 27-column (nfaclyrs = 5) factypes row with an effective
conductivity of 1000 W/m/K is appended, and short MPI runs verify that

 * the water column stays isothermal at the anchor temperature (waterT, or
   flrT when waterT is not set), while ordinary facets keep the bldT anchor;
 * water facets produce a negative (evaporative) latent heat flux ef through
   the saturated, resistance-free open-water branch, while non-vegetated,
   non-water facets produce none.

Skips when the executable, an MPI launcher, or trimesh are unavailable.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import f90nml
import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / "examples" / "201"
EXECUTABLE = Path(os.environ.get("UDALES_EXE", ROOT / "bin" / "u-dales")).resolve()

CASE = 201
FLRT = 295.0
BLDT = 301.0
WATER_LAMBDA = 1000.0
WATER_C = 4.18e6
RUNTIME = 30.0


def _water_facet_ids() -> np.ndarray:
    """Ground-level upward facets west of x=64 m, in STL (= facets.inp) order."""
    import trimesh

    mesh = trimesh.load(FIXTURE / "geom.201.STL", process=False)
    normals = np.asarray(mesh.face_normals)
    centres = np.asarray(mesh.triangles_center)
    ground = (normals[:, 2] > 0.99) & (np.abs(centres[:, 2]) < 1e-3)
    return np.flatnonzero(ground & (centres[:, 0] < 64.0))


class WaterFacetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not EXECUTABLE.is_file():
            raise unittest.SkipTest(f"uDALES executable not found: {EXECUTABLE}")
        cls.mpiexec = shutil.which(os.environ.get("MPIEXEC", "mpiexec"))
        if cls.mpiexec is None:
            raise unittest.SkipTest("MPI launcher not found")
        try:
            cls.water = _water_facet_ids()
        except ImportError:
            raise unittest.SkipTest("trimesh not available")
        cls.types = np.loadtxt(FIXTURE / f"facets.inp.{CASE}", skiprows=1)[:, 0].astype(int)

    def _build_case(self, case_dir: Path, *, watert: float | None) -> None:
        case_dir.mkdir(parents=True)
        for source in FIXTURE.iterdir():
            if source.name in (f"namoptions.{CASE}", f"facets.inp.{CASE}",
                               f"factypes.inp.{CASE}", "config.sh", "info.txt"):
                continue
            shutil.copy(source, case_dir / source.name)

        lines = (FIXTURE / f"facets.inp.{CASE}").read_text().splitlines()
        for i in self.water:
            parts = lines[i + 1].split()
            parts[0] = "-21"
            lines[i + 1] = "  " + "  ".join(parts)
        (case_dir / f"facets.inp.{CASE}").write_text("\n".join(lines) + "\n")

        # nfaclyrs = 5 here, so the row must be in the 27-column layout:
        # id lGR z0 z0h al em d1..d5 C1..C5 l1..l5 k1..k6
        k = WATER_LAMBDA / WATER_C
        row = ("     -21    0   0.003  0.00003    0.06    0.95" + "    0.30" * 5
               + f"  {WATER_C:14.0f}" * 5 + f"  {WATER_LAMBDA:12.4f}" * 5
               + f"  {k:14.8e}" * 6)
        text = (FIXTURE / f"factypes.inp.{CASE}").read_text().rstrip("\n")
        (case_dir / f"factypes.inp.{CASE}").write_text(text + "\n" + row + "\n")

        namelist = f90nml.read(FIXTURE / f"namoptions.{CASE}")
        namelist["run"].update(runtime=RUNTIME, trestart=1.0e6)
        output = namelist["output"]
        output.pop("tfielddump", None)
        output.update(lfielddump=False, ltdump=False, lxytdump=False,
                      tinstantdump=1.0e6)
        eb = namelist["energybalance"]
        eb["flrt"] = FLRT
        if watert is not None:
            eb["watert"] = watert
        namelist.write(case_dir / f"namoptions.{CASE}", force=True)

    def _run(self, case_dir: Path) -> None:
        result = subprocess.run(
            [self.mpiexec, "-n", "4", str(EXECUTABLE), f"namoptions.{CASE}"],
            cwd=case_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=1800,
        )
        tail = "\n".join(result.stdout.splitlines()[-30:])
        self.assertEqual(result.returncode, 0, f"solver failed:\n{tail}")

    def _facet_temperature(self, case_dir: Path) -> np.ndarray:
        with Dataset(case_dir / f"facT.{CASE}.nc") as ds:
            ds.set_auto_mask(False)
            return np.asarray(ds.variables["T"][-1])  # (layers, facets) at t_end

    def test_water_column_anchored_at_watert_and_evaporating(self):
        watert = 296.5
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = Path(tmp) / "case"
            self._build_case(case_dir, watert=watert)
            self._run(case_dir)

            T = self._facet_temperature(case_dir)
            water = np.zeros(len(self.types), dtype=bool)
            water[self.water] = True
            # column isothermal at the anchor; deep node exactly the anchor
            self.assertLess(np.abs(T[-1, :][water] - watert).max(), 1e-6)
            self.assertLess(np.abs(T[:, water] - watert).max(), 0.5)
            # ordinary facets keep the building anchor
            self.assertLess(np.abs(T[-1, :][~water] - BLDT).max(), 1e-6)

            with Dataset(case_dir / f"facEB.{CASE}.nc") as ds:
                ds.set_auto_mask(False)
                ef = np.asarray(ds.variables["ef"][-1])
            ef_water = ef[self.water]
            # evaporation: energy leaves the facet (uDALES sign: into facet > 0)
            self.assertLess(ef_water.mean(), -1.0)
            self.assertGreater(ef_water.min(), -2000.0)
            # non-water, non-vegetated facets have no latent flux
            self.assertEqual(float(np.abs(ef[~water]).max()), 0.0)

    def test_default_anchor_falls_back_to_flrt(self):
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = Path(tmp) / "case"
            self._build_case(case_dir, watert=None)
            self._run(case_dir)
            T = self._facet_temperature(case_dir)
            water = np.zeros(len(self.types), dtype=bool)
            water[self.water] = True
            self.assertLess(np.abs(T[-1, :][water] - FLRT).max(), 1e-6)
            self.assertLess(np.abs(T[:, water] - FLRT).max(), 0.5)


if __name__ == "__main__":
    unittest.main()
