import sys
import unittest
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from _common import PYTHON_DIR

from udprep.udprep_vegetation import VegetationSection, vegetation_block_to_veg as module_veg_fn  # noqa: E402
from udbase import UDBase  # noqa: E402

class DummySim:
    # Reuse the real vegetation loader/cache so the provenance-flag regression
    # test exercises UDBase.load_veg's actual cache condition.
    load_veg = UDBase.load_veg
    _load_sparse_file = UDBase._load_sparse_file

    def __init__(self, path: Path, expnr="777"):
        self.path = path
        self.expnr = expnr
        self.itot = 8
        self.jtot = 8
        self.ktot = 6
        self.ntrees = 0

class TestVegetationSection(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        self.sim = DummySim(self.workdir)
        (self.workdir / f"namoptions.{self.sim.expnr}").write_text(
            "&RUN\niexpnr = 777\n/\n", encoding="ascii"
        )

    def _make_section(self, **overrides):
        values = {
            "treesfile": f"trees.inp.{self.sim.expnr}",
            "lad": 1.0,
            "cd": 0.3,
            "ud": 0.0002,
            "dec": 0.0,
            "lsize": 0.15,
            "r_s": 0.0,
            "ltrees": True,
            "ltreesfile": True,
        }
        values.update(overrides)
        return VegetationSection("vegetation", values, sim=self.sim, defaults={})

    def _make_veg(self, points_0based, **param_overrides):
        """Build a minimal veg dict with uniform params (0-based points)."""
        n = len(points_0based)
        params = {
            "id":    np.arange(1, n + 1, dtype=int),
            "lad":   np.ones(n, dtype=float),
            "cd":    np.full(n, 0.3),
            "ud":    np.full(n, 0.0002),
            "dec":   np.zeros(n),
            "lsize": np.full(n, 0.15),
            "r_s":   np.zeros(n),
        }
        params.update(param_overrides)
        return {"points": np.array(points_0based, dtype=int), "params": params}

    def test_load_block_populates_sparse_points_without_geometry(self):
        (self.workdir / "trees.inp.777").write_text("1 2 3 4 1 2\n", encoding="ascii")
        section = self._make_section()

        section.load_block()

        self.assertEqual(section.ntrees, 8)
        np.testing.assert_array_equal(section.veg["points"][0], [0, 2, 0])
        self.assertEqual(section.veg["params"]["id"][0], 1)
        self.assertTrue(np.all(section.veg["params"]["lad"] == 1.0))

    def test_section_built_veg_survives_cached_load_veg(self):
        """load_block must set the veg provenance flag so a cached load_veg
        before save() returns the section-built veg (cache hit) rather than
        re-reading a nonexistent veg.inp and replacing it with an empty dict."""
        (self.workdir / "trees.inp.777").write_text("1 2 3 4 1 2\n", encoding="ascii")
        section = self._make_section()

        section.load_block()  # no save() -> veg.inp.777 is not on disk
        self.assertGreater(len(self.sim.veg["points"]), 0)

        veg = self.sim.load_veg(zero_based=True, cache=True)

        # Cache hit: the populated section-built veg must be returned intact.
        self.assertGreater(len(veg["points"]), 0)
        self.assertGreater(len(self.sim.veg["points"]), 0)

    def test_vegetation_block_to_veg_writes_sparse_files(self):
        (self.workdir / "trees.inp.777").write_text("1 2 1 1 1 1\n", encoding="ascii")
        section = self._make_section()

        out = section.vegetation_block_to_veg()

        veg_path = self.workdir / "veg.inp.777"
        params_path = self.workdir / "veg_params.inp.777"
        self.assertEqual(out["veg"], veg_path)
        self.assertEqual(out["params"], params_path)
        self.assertTrue(veg_path.exists())
        self.assertTrue(params_path.exists())
        veg_lines = veg_path.read_text(encoding="ascii").splitlines()
        params_lines = params_path.read_text(encoding="ascii").splitlines()
        self.assertEqual(veg_lines[0], "# position (i,j,k)")
        self.assertEqual(params_lines[0], "# id lad cd ud dec lsize r_s")
        self.assertEqual(len(veg_lines), 3)
        self.assertEqual(len(params_lines), 3)

    def test_load_block_rejects_out_of_bounds_block(self):
        (self.workdir / "trees.inp.777").write_text("1 9 1 1 1 1\n", encoding="ascii")
        section = self._make_section()
        with self.assertRaises(ValueError):
            section.load_block()

    def test_load_block_rejects_missing_tree_file(self):
        section = self._make_section()
        with self.assertRaises(FileNotFoundError):
            section.load_block()

    def test_load_block_accepts_legacy_seven_column_format(self):
        (self.workdir / "trees.inp.777").write_text("5 2 1 3 3 2 1\n", encoding="ascii")
        section = self._make_section()
        section.load_block()
        self.assertEqual(section.ntrees, 4)
        self.assertTrue(np.all(section.veg["params"]["id"] == 1))

    # ------------------------------------------------------------------
    # save() — dedicated tests
    # ------------------------------------------------------------------

    def test_save_writes_1based_indices_to_veg_file(self):
        """0-based point (4, 6, 16) must be written as (5, 7, 17) in veg.inp."""
        section = self._make_section()
        self.sim.veg = self._make_veg([(4, 6, 16)])

        section.save()

        lines = (self.workdir / "veg.inp.777").read_text(encoding="ascii").splitlines()
        self.assertEqual(lines[0], "# position (i,j,k)")
        parts = lines[1].split()
        self.assertEqual(parts, ["5", "7", "17"])

    def test_save_writes_correct_params_to_params_file(self):
        """veg_params.inp must contain the correct id, lad, and uniform scalar params."""
        section = self._make_section()
        self.sim.veg = self._make_veg([(0, 0, 0)])

        section.save()

        lines = (self.workdir / "veg_params.inp.777").read_text(encoding="ascii").splitlines()
        self.assertEqual(lines[0], "# id lad cd ud dec lsize r_s")
        parts = lines[1].split()
        self.assertEqual(int(parts[0]), 1)
        self.assertAlmostEqual(float(parts[1]), 1.0)   # lad
        self.assertAlmostEqual(float(parts[2]), 0.3)   # cd

    def test_save_warns_when_param_varies_per_point(self):
        """_coerce_scalar must emit RuntimeWarning when a param differs across points."""
        section = self._make_section()
        self.sim.veg = self._make_veg([(0, 0, 0), (1, 1, 1)], cd=np.array([0.3, 0.6]))

        with self.assertWarns(RuntimeWarning) as ctx:
            section.save()

        self.assertIn("cd", str(ctx.warning))

    def test_save_uses_first_value_when_param_varies(self):
        """When a param is not uniform, the first value must be written for all rows."""
        section = self._make_section()
        self.sim.veg = self._make_veg([(0, 0, 0), (1, 1, 1)], cd=np.array([0.3, 0.9]))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            section.save()

        lines = (self.workdir / "veg_params.inp.777").read_text(encoding="ascii").splitlines()
        cd_row1 = float(lines[1].split()[2])
        cd_row2 = float(lines[2].split()[2])
        self.assertAlmostEqual(cd_row1, 0.3)   # first value used
        self.assertAlmostEqual(cd_row2, 0.3)   # same scalar written for both rows

    def test_load_stl_plot_false_is_side_effect_free(self):
        # P28: load_stl(plot=False) must load vegetation without touching a
        # rendering backend, returning the veg dict; plot=True keeps the
        # opt-in visualisation.
        try:
            import trimesh
        except ImportError:
            self.skipTest("trimesh not available")
        from unittest.mock import MagicMock

        mesh = trimesh.creation.box(extents=(4.0, 4.0, 4.0))  # watertight, [-2,2]^3
        mesh.export(str(self.workdir / "cube.stl"))

        self.sim.dx = 1.0
        self.sim.dy = 1.0
        self.sim.dzt = np.array([1.0, 1.0, 1.0])
        self.sim.xt = np.array([-1.0, 0.0, 1.0])
        self.sim.yt = np.array([-1.0, 0.0, 1.0])
        self.sim.zt = np.array([-1.0, 0.0, 1.0])
        self.sim.vis = MagicMock()

        section = self._make_section(treesfile="cube.stl")

        veg = section.load_stl(plot=False)
        self.assertIsInstance(veg, dict)
        self.assertIn("points", veg)
        self.assertGreater(len(veg["points"]), 0)
        self.sim.vis.plot_veg.assert_not_called()

        fig = section.load_stl(plot=True)
        self.sim.vis.plot_veg.assert_called_once()
        self.assertIs(fig, self.sim.vis.plot_veg.return_value)

    def test_save_raises_on_missing_param(self):
        """Omitting a required param key must raise ValueError."""
        section = self._make_section()
        veg = self._make_veg([(0, 0, 0)])
        del veg["params"]["cd"]
        self.sim.veg = veg

        with self.assertRaises(ValueError):
            section.save()

class TestModuleLevelVegetationBlockToVeg(unittest.TestCase):
    """Tests for the module-level vegetation_block_to_veg shim function."""

    def _fake_udprep_cls(self, return_value):
        """Return a minimal UDPrep stand-in whose vegetation method returns *return_value*."""
        rv = return_value

        class FakeUDPrep:
            def __init__(self, sim=None):
                from unittest.mock import MagicMock
                self.vegetation = MagicMock()
                self.vegetation.vegetation_block_to_veg.return_value = rv

        return FakeUDPrep

    def test_wraps_non_udprep_in_udprep_and_delegates(self):
        """A non-UDPrep argument must be wrapped in UDPrep before delegation."""
        expected = {"veg": Path("veg.inp.777"), "params": Path("veg_params.inp.777")}
        FakeUDPrep = self._fake_udprep_cls(expected)

        with patch("udprep.udprep.UDPrep", FakeUDPrep):
            result = module_veg_fn(object(), filename="custom.inp")

        self.assertEqual(result, expected)

    def test_uses_udprep_instance_directly_without_rewrapping(self):
        """An existing UDPrep instance must be used as-is; no new UDPrep is constructed."""
        expected = {"veg": Path("veg.inp.777"), "params": Path("veg_params.inp.777")}
        FakeUDPrep = self._fake_udprep_cls(expected)
        fake_prep = FakeUDPrep()   # already a FakeUDPrep instance

        # Replace __init__ AFTER fake_prep is created so we can detect any new construction.
        new_construction_calls = []
        original_init = FakeUDPrep.__init__

        def tracking_init(self_inner, sim=None):
            if sim is not None:
                new_construction_calls.append(sim)
            original_init(self_inner, sim)

        FakeUDPrep.__init__ = tracking_init

        with patch("udprep.udprep.UDPrep", FakeUDPrep):
            result = module_veg_fn(fake_prep)

        self.assertEqual(result, expected)
        self.assertEqual(
            new_construction_calls, [],
            "UDPrep must not be re-constructed when an existing UDPrep instance is passed",
        )

if __name__ == "__main__":
    unittest.main()
