import io
import os
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udbase import UDBase  # noqa: E402
from udvis import UDVis  # noqa: E402


class TestUDBaseCore(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        (self.workdir / "namoptions.001").write_text(
            "\n".join(
                [
                    "&DOMAIN",
                    " itot = 4",
                    " jtot = 3",
                    " ktot = 2",
                    " xlen = 40.0",
                    " ylen = 30.0",
                    " zsize = 20.0",
                    "/",
                ]
            )
            + "\n",
            encoding="ascii",
        )

    def test_constructor_initializes_visualization_facade(self):
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        self.assertIsInstance(sim.vis, UDVis)
        self.assertIs(sim.vis.sim, sim)

    def test_constructor_expands_user_home_in_path(self):
        fake_home = self.workdir / "home"
        expdir = fake_home / "simulation" / "udtest" / "experiments" / "001"
        expdir.mkdir(parents=True)
        (expdir / "namoptions.001").write_text(
            "\n".join(
                [
                    "&DOMAIN",
                    " itot = 4",
                    " jtot = 3",
                    " ktot = 2",
                    " xlen = 40.0",
                    " ylen = 30.0",
                    " zsize = 20.0",
                    "/",
                ]
            )
            + "\n",
            encoding="ascii",
        )

        with patch.dict(os.environ, {"HOME": str(fake_home)}):
            sim = UDBase(
                "1",
                "~/simulation/udtest/experiments/001",
                load_geometry=False,
                suppress_load_warnings=True,
            )

        self.assertEqual(sim.path, expdir)

    def test_constructor_sets_scalar_defaults_when_namoptions_omits_scalars(self):
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        self.assertEqual(sim.nsv, 0)
        self.assertFalse(sim.lscasrc)
        self.assertFalse(sim.lscasrcl)
        self.assertEqual(sim.nscasrc, 0)
        self.assertEqual(sim.nscasrcl, 0)

    def test_load_scalar_sources_reads_point_and_line_sources_by_scalar_index(self):
        (self.workdir / "namoptions.001").write_text(
            "\n".join(
                [
                    "&DOMAIN",
                    " itot = 4",
                    " jtot = 3",
                    " ktot = 2",
                    " xlen = 40.0",
                    " ylen = 30.0",
                    " zsize = 20.0",
                    "/",
                    "&SCALARS",
                    " nsv = 2",
                    " lscasrc = .true.",
                    " lscasrcl = .true.",
                    " nscasrc = 2",
                    " nscasrcl = 1",
                    "/",
                ]
            )
            + "\n",
            encoding="ascii",
        )
        (self.workdir / "scalarsourcep.inp.1.001").write_text(
            "# Scalar point source data\n"
            "#xS yS zS SS sigS\n"
            "1 2 3 4 5\n"
            "6 7 8 9 10\n",
            encoding="ascii",
        )
        (self.workdir / "scalarsourcep.inp.2.001").write_text(
            "# Scalar point source data\n"
            "#xS yS zS SS sigS\n"
            "11 12 13 14 15\n"
            "16 17 18 19 20\n",
            encoding="ascii",
        )
        (self.workdir / "scalarsourcel.inp.1.001").write_text(
            "# Scalar line source data\n"
            "#xSb ySb zSb xSe ySe zSe SS sigS\n"
            "1 2 3 4 5 6 7 8\n",
            encoding="ascii",
        )
        (self.workdir / "scalarsourcel.inp.2.001").write_text(
            "# Scalar line source data\n"
            "#xSb ySb zSb xSe ySe zSe SS sigS\n"
            "9 10 11 12 13 14 15 16\n",
            encoding="ascii",
        )
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        sources = sim.load_scalar_sources()

        self.assertEqual(sorted(sources["point"]), [1, 2])
        self.assertEqual(sorted(sources["line"]), [1, 2])
        np.testing.assert_allclose(
            sources["point"][1],
            np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]], dtype=float),
        )
        np.testing.assert_allclose(
            sources["line"][2],
            np.array([[9, 10, 11, 12, 13, 14, 15, 16]], dtype=float),
        )

    def test_load_scalar_sources_warns_on_missing_file_name(self):
        (self.workdir / "namoptions.001").write_text(
            "\n".join(
                [
                    "&DOMAIN",
                    " itot = 4",
                    " jtot = 3",
                    " ktot = 2",
                    " xlen = 40.0",
                    " ylen = 30.0",
                    " zsize = 20.0",
                    "/",
                    "&SCALARS",
                    " nsv = 1",
                    " lscasrc = .true.",
                    " nscasrc = 1",
                    "/",
                ]
            )
            + "\n",
            encoding="ascii",
        )
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=False)

        stderr = io.StringIO()
        with redirect_stderr(stderr):
            sources = sim.load_scalar_sources()

        self.assertEqual(sources, {"point": {}, "line": {}})
        self.assertIn(
            "scalarsourcep.inp.1.001 not found.",
            stderr.getvalue(),
        )

    def test_load_scalar_sources_warns_on_load_error(self):
        (self.workdir / "namoptions.001").write_text(
            "\n".join(
                [
                    "&DOMAIN",
                    " itot = 4",
                    " jtot = 3",
                    " ktot = 2",
                    " xlen = 40.0",
                    " ylen = 30.0",
                    " zsize = 20.0",
                    "/",
                    "&SCALARS",
                    " nsv = 1",
                    " lscasrc = .true.",
                    " nscasrc = 1",
                    "/",
                ]
            )
            + "\n",
            encoding="ascii",
        )
        (self.workdir / "scalarsourcep.inp.1.001").write_text(
            "# Scalar point source data\n"
            "#xS yS zS SS sigS\n"
            "not-a-number\n",
            encoding="ascii",
        )
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        with self.assertWarnsRegex(UserWarning, "Error loading scalarsourcep.inp.1.001"):
            sources = sim.load_scalar_sources()

        self.assertEqual(sources, {"point": {}, "line": {}})

    def test_load_facsec_applies_one_based_to_zero_based_mapping(self):
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)
        sim.ffacet_sections = "facet_sections"
        sim.ffluid_boundary = "fluid_boundary"
        (self.workdir / "facet_sections_c.txt").write_text(
            "\n".join(
                [
                    "# facid area fluid_idx distance",
                    "1 2.5 2 0.10",
                    "3 1.0 1 0.20",
                ]
            )
            + "\n",
            encoding="ascii",
        )
        (self.workdir / "fluid_boundary_c.txt").write_text(
            "\n".join(
                [
                    "# i j k",
                    "1 2 1",
                    "4 3 2",
                ]
            )
            + "\n",
            encoding="ascii",
        )

        facsec = sim.load_facsec("c")

        np.testing.assert_array_equal(facsec["facid"], np.array([0, 2]))
        np.testing.assert_allclose(facsec["area"], np.array([2.5, 1.0]))
        np.testing.assert_array_equal(facsec["locs"], np.array([[3, 2, 1], [0, 1, 0]]))
        np.testing.assert_allclose(facsec["distance"], np.array([0.10, 0.20]))

    def test_convert_fac_to_field_accumulates_cell_density(self):
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)
        facsec = {
            "facid": np.array([0, 1, 0]),
            "area": np.array([2.0, 1.0, 3.0]),
            "locs": np.array([[0, 0, 0], [0, 0, 1], [0, 0, 0]]),
        }
        var = np.array([10.0, 20.0])

        field = sim.convert_fac_to_field(var, facsec=facsec)

        self.assertEqual(field.shape, (4, 3, 2))
        self.assertAlmostEqual(field[0, 0, 0], 0.05)
        self.assertAlmostEqual(field[0, 0, 1], 0.02)
        self.assertTrue(np.all(field[1:, :, :] == 0.0))


class TestTreeLoading(unittest.TestCase):
    """Golden-file coverage for ``UDBase._load_tree_data`` / ``trees.inp.<expnr>``.

    Pins the current format so future format drift is caught here rather than in
    the tutorials: ``#`` header line(s), then rows of 1-based ``il iu jl ju kl ku``
    bounding-box indices, stored 0-based as ``(n_trees, 6)``.
    """

    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        (self.workdir / "namoptions.001").write_text(
            "&DOMAIN\n itot = 4\n jtot = 3\n ktot = 2\n"
            " xlen = 40.0\n ylen = 30.0\n zsize = 20.0\n/\n",
            encoding="ascii",
        )

    def _sim_with_trees(self, text):
        (self.workdir / "trees.inp.001").write_text(text, encoding="ascii")
        return UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

    def test_current_format_two_headers_six_columns_stored_zero_based(self):
        sim = self._sim_with_trees(
            "# Trees data\n"
            "# tree_n\t il\t iu\t jl\t ju\t kl\t ku\n"
            " 129\t 173\t 106\t 150\t   5\t  17\n"
        )
        self.assertTrue(sim._lftrees)
        self.assertEqual(sim.trees.shape, (1, 6))
        np.testing.assert_array_equal(sim.trees, [[128, 172, 105, 149, 4, 16]])

    def test_multiple_trees(self):
        sim = self._sim_with_trees("#h\n#h\n1 2 3 4 5 6\n7 8 9 10 11 12\n")
        np.testing.assert_array_equal(
            sim.trees, [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11]]
        )

    def test_optional_leading_tree_n_column_is_dropped(self):
        sim = self._sim_with_trees("#h\n#h\n1 129 173 106 150 5 17\n")
        np.testing.assert_array_equal(sim.trees, [[128, 172, 105, 149, 4, 16]])

    def test_robust_to_header_row_count(self):
        # comment-based skipping, not a fixed skiprows=2.
        sim = self._sim_with_trees("# one header line only\n 129 173 106 150 5 17\n")
        np.testing.assert_array_equal(sim.trees, [[128, 172, 105, 149, 4, 16]])

    def test_empty_or_header_only_file_gives_empty_array(self):
        sim = self._sim_with_trees("# Trees data\n# header only, no trees\n")
        self.assertTrue(sim._lftrees)
        self.assertEqual(sim.trees.shape, (0, 6))

    def test_missing_file_sets_flag_false_and_trees_none(self):
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)
        self.assertFalse(sim._lftrees)
        self.assertIsNone(sim.trees)


class TestReadMatrix(unittest.TestCase):
    """`read_matrix` is a @staticmethod, so `sim.read_matrix(path, n)` must not
    pass the instance as the filename (regression: it lacked @staticmethod)."""

    def test_read_matrix_called_on_instance_skips_header(self):
        with TemporaryDirectory() as d:
            p = Path(d) / "data.txt"
            p.write_text("header row\n1 2 3\n4 5 6\n", encoding="ascii")
            # Access through an instance, exactly as udprep_forcing does.
            sim = UDBase.__new__(UDBase)
            arr = sim.read_matrix(p, 1)
        np.testing.assert_array_equal(arr, [[1, 2, 3], [4, 5, 6]])
