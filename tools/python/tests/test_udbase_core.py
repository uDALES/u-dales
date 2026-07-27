import io
import os
import sys
import unittest
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from _common import PYTHON_DIR

from udbase import UDBase  # noqa: E402
from udvis import UDVis  # noqa: E402

from exceptions import DataFormatError  # noqa: E402

_DOMAIN_NAMOPTIONS = "\n".join(
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

def _factypes_row(nfaclyrs, wallid=1):
    """Build one valid factypes data row: 6 header + d(k) + C(k) + l(k) + k(k+1)."""
    k = nfaclyrs
    fields = [str(wallid), "0", "0.1", "0.01", "0.2", "0.9"]
    fields += ["0.1"] * k          # d: layer thicknesses
    fields += ["1000"] * k         # C: heat capacity
    fields += ["1.0"] * k          # l: conductivity
    fields += ["0.4"] * (k + 1)    # k: conductivity (k+1 columns)
    return " ".join(fields)

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

    def test_module_docstring_is_real_docstring(self):
        # C12: the module docstring must precede `from __future__ import ...`
        # so it becomes the actual __doc__ instead of a discarded expression.
        import udbase

        self.assertIsNotNone(udbase.__doc__)
        self.assertIn("Post-Processing", udbase.__doc__)

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

        # Path.expanduser() reads HOME on POSIX but USERPROFILE on Windows, so
        # patch both to make the test platform-neutral.
        with patch.dict(os.environ, {"HOME": str(fake_home), "USERPROFILE": str(fake_home)}):
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

        with self.assertWarnsRegex(UserWarning, "scalarsourcep.inp.1.001 not found."):
            sources = sim.load_scalar_sources()

        self.assertEqual(sources, {"point": {}, "line": {}})

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

    def test_repr_does_not_crash_when_geometry_not_loaded(self):
        # C1: with load_geometry=False, self.geom must still be initialized so
        # repr(sim) (which reads self.geom) does not raise AttributeError.
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        self.assertIsNone(sim.geom)
        text = repr(sim)  # must not raise
        self.assertIn("UDBase(expnr='001')", text)

    def test_repr_is_concise_and_describe_holds_full_dump(self):
        # C24: __repr__ is a short summary that points to describe(); the full
        # per-attribute dump (private flags, every scalar/array) lives in
        # describe().
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        r = repr(sim)
        self.assertIn("UDBase(expnr='001')", r)
        self.assertIn("describe()", r)
        self.assertNotIn("_lfprof", r)  # private flags not dumped in repr

        d = sim.describe()
        self.assertIn("UDBase(expnr='001')", d)
        self.assertIn("_lfprof", d)  # full dump includes private flags
        self.assertIn("itot", d)
        # describe() is strictly the longer of the two
        self.assertGreater(len(d.splitlines()), len(r.splitlines()))

    def test_load_facsec_rejects_invalid_var(self):
        # C29: load_facsec validates its grid designator like load_slice does,
        # instead of raising an obscure FileNotFoundError for an impossible path.
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)
        with self.assertRaises(ValueError):
            sim.load_facsec("x")

    def test_facsec_missing_error_message_names_real_files(self):
        # C22: the "requires facet section data" message must reference the real
        # filenames (facet_sections_(u,v,w,c).txt), not the nonexistent
        # facet_sections_(u,v,w,c).001.
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)
        del sim.facsec  # force the missing-facsec branch (facsec is normally {})

        with self.assertRaises(ValueError) as ctx:
            sim.convert_fac_to_field(np.array([1.0]))

        msg = str(ctx.exception)
        self.assertIn("facet_sections_(u,v,w,c).txt", msg)
        self.assertIn("fluid_boundary_(u,v,w,c).txt", msg)
        self.assertNotIn(".001", msg)

    def test_load_veg_cache_respects_zero_based_flag(self):
        # C4: the cache must not return 0-based data for a zero_based=False call.
        (self.workdir / "veg.inp.001").write_text(
            "# i j k\n2 3 4\n", encoding="ascii"
        )
        (self.workdir / "veg_params.inp.001").write_text(
            "# id lad cd ud dec lsize r_s\n1 0.1 0.2 0.3 0.4 0.5 0.6\n",
            encoding="ascii",
        )
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        veg0 = sim.load_veg(zero_based=True)
        veg1 = sim.load_veg(zero_based=False)

        np.testing.assert_array_equal(veg1["points"] - veg0["points"], 1)

    def test_single_row_facet_files_load_with_correct_shapes(self):
        # C7: one-facet files must not collapse to 1-D/0-D and raise IndexError
        # (which was swallowed into a spurious "missing facet data" warning).
        (self.workdir / "facetarea.inp.001").write_text(
            "# area\n2.5\n", encoding="ascii"
        )
        (self.workdir / "facets.inp.001").write_text(
            "# typeid nx ny nz\n1 0.0 0.0 1.0\n", encoding="ascii"
        )
        # A valid 3-layer factypes row: 6 + 4*3 + 1 = 19 columns
        # (6 header + d(3) + C(3) + l(3) + k(4)).
        (self.workdir / "factypes.inp.001").write_text(
            "# header line 1\n# header line 2\n# header line 3\n"
            + _factypes_row(3)
            + "\n",
            encoding="ascii",
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        self.assertEqual([str(w.message) for w in caught], [])
        self.assertEqual(sim.facs["area"].shape, (1,))
        self.assertEqual(sim.facs["typeid"].shape, (1,))
        self.assertEqual(sim.facs["normals"].shape, (1, 3))
        self.assertEqual(sim.factypes["id"].shape, (1,))

    # ----- C2: factypes column-width validation -----

    def _write_namoptions(self, extra_lines=()):
        text = _DOMAIN_NAMOPTIONS
        if extra_lines:
            text += "\n" + "\n".join(extra_lines)
        (self.workdir / "namoptions.001").write_text(text + "\n", encoding="ascii")

    def _write_factypes(self, rows):
        (self.workdir / "factypes.inp.001").write_text(
            "# header line 1\n# header line 2\n# header line 3\n"
            + "\n".join(rows)
            + "\n",
            encoding="ascii",
        )

    def test_factypes_non_default_nfaclyrs_parses_end_to_end(self):
        # (a) A consistent file with a non-default layer count (5) parses, and the
        # property blocks are sliced to the right widths.
        self._write_namoptions(["&WALLS", " nfaclyrs = 5", "/"])
        self._write_factypes([_factypes_row(5, wallid=1), _factypes_row(5, wallid=2)])

        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        self.assertEqual(sim.nfaclyrs, 5)
        self.assertEqual(sim.factypes["d"].shape, (2, 5))
        self.assertEqual(sim.factypes["C"].shape, (2, 5))
        self.assertEqual(sim.factypes["lam"].shape, (2, 6))  # k + 1 columns

    def test_factypes_width_mismatch_raises_dataformaterror(self):
        # (b) nfaclyrs says 5 but the file is a 3-layer (19-col) file: the layer
        # count disagrees, so material properties would be mis-sliced -> raise.
        self._write_namoptions(["&WALLS", " nfaclyrs = 5", "/"])
        self._write_factypes([_factypes_row(3)])  # 19 columns, not 27

        with self.assertRaises(DataFormatError) as ctx:
            UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        msg = str(ctx.exception)
        self.assertIn("factypes.inp.001", msg)
        self.assertIn("nfaclyrs=5", msg)
        self.assertIn("27", msg)  # expected column count
        self.assertIn("19", msg)  # actual column count

    def test_factypes_empty_file_raises_dataformaterror(self):
        # A header-only factypes file has no data rows; np.loadtxt returns an
        # empty array whose "width" is meaningless. Raise a clear "no data rows"
        # error instead of a misleading width-mismatch message.
        self._write_namoptions(["&WALLS", " nfaclyrs = 3", "/"])
        (self.workdir / "factypes.inp.001").write_text(
            "# header line 1\n# header line 2\n# header line 3\n",
            encoding="ascii",
        )

        with self.assertRaises(DataFormatError) as ctx:
            UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        msg = str(ctx.exception)
        self.assertIn("factypes.inp.001", msg)
        self.assertIn("no data rows", msg)

    def test_factypes_default_guess_with_matching_three_layer_file(self):
        # (c) nfaclyrs absent from namoptions -> guessed 3; a matching 19-col file
        # still loads.
        self._write_namoptions()  # no nfaclyrs
        self._write_factypes([_factypes_row(3)])  # 19 columns

        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        self.assertEqual(sim.nfaclyrs, 3)
        self.assertEqual(sim.factypes["d"].shape, (1, 3))
        self.assertEqual(sim.factypes["lam"].shape, (1, 4))

    # ----- C3: prof.inp grid parse errors must not silently degrade -----

    def test_missing_prof_inp_warns_and_uses_uniform_grid(self):
        # (a) No prof.inp present: warn, fall back to a uniform grid, stay usable.
        with self.assertWarnsRegex(UserWarning, r"prof\.inp\.001 not found"):
            sim = UDBase("1", self.workdir, load_geometry=False)

        self.assertFalse(sim._lfprof)
        # uniform z-grid for ktot=2, zsize=20 -> centres at 5, 15
        np.testing.assert_allclose(sim.zt, [5.0, 15.0])

    def test_corrupt_prof_inp_raises_dataformaterror(self):
        # (b) prof.inp present but unparseable: raise, do not silently substitute
        # a uniform grid.
        (self.workdir / "prof.inp.001").write_text(
            "# z thl qt u v tke\nnot a number here\n", encoding="ascii"
        )

        with self.assertRaises(DataFormatError) as ctx:
            UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        self.assertIn("prof.inp.001", str(ctx.exception))

    def test_valid_stretched_prof_inp_loads(self):
        # (c) a valid stretched prof.inp is read and used verbatim.
        (self.workdir / "prof.inp.001").write_text(
            "# stretched profile\n# z thl qt u v tke\n"
            "2.0 288 0 1 0 0\n14.0 288 0 1 0 0\n",
            encoding="ascii",
        )

        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        self.assertTrue(sim._lfprof)
        np.testing.assert_allclose(sim.zt, [2.0, 14.0])

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

class TestGridHelpers(unittest.TestCase):
    """Pure grid-coordinate maths (udgrid)."""

    def test_horizontal_grid_edges_and_centres(self):
        import udgrid

        xm, xt, ym, yt = udgrid.horizontal_grid(3, 2.0, 2, 5.0)
        np.testing.assert_allclose(xm, [0, 2, 4])
        np.testing.assert_allclose(xt, [1, 3, 5])
        np.testing.assert_allclose(ym, [0, 5])
        np.testing.assert_allclose(yt, [2.5, 7.5])

    def test_uniform_z_grid(self):
        import udgrid

        zm, zt, dzt = udgrid.uniform_z_grid(20.0, 4)
        np.testing.assert_allclose(zm, [0, 5, 10, 15])
        np.testing.assert_allclose(zt, [2.5, 7.5, 12.5, 17.5])
        np.testing.assert_allclose(dzt, [5, 5, 5, 5])

    def test_z_grid_from_profile_matches_centres(self):
        import udgrid

        zt = np.array([2.5, 7.5, 12.5, 17.5])
        zm, dzt = udgrid.z_grid_from_profile(zt, 20.0)
        np.testing.assert_allclose(zm, [0, 5, 10, 15])
        np.testing.assert_allclose(dzt, [5, 5, 5, 5])

class TestNamoptionsParsing(unittest.TestCase):
    """Pure namelist parsing (udconfig)."""

    def test_parse_value_types(self):
        import udconfig

        self.assertIs(udconfig.parse_value(".true."), True)
        self.assertIs(udconfig.parse_value(".FALSE."), False)
        self.assertEqual(udconfig.parse_value("42"), 42)
        self.assertEqual(udconfig.parse_value("1.5"), 1.5)
        self.assertEqual(udconfig.parse_value("1e3"), 1000.0)
        self.assertEqual(udconfig.parse_value("'geom.stl'"), "geom.stl")

    def test_parse_namoptions_skips_headers_and_comments_and_applies_defaults(self):
        import udconfig

        with TemporaryDirectory() as d:
            p = Path(d) / "namoptions.001"
            p.write_text(
                "&RUN\n! a comment\n itot = 4\n xlen = 40.0  ! inline\n"
                " ltest = .true.\n stl_file = 'geom.stl'\n/\n",
                encoding="ascii",
            )
            vals = udconfig.parse_namoptions(p)
        self.assertEqual(vals["itot"], 4)
        self.assertEqual(vals["xlen"], 40.0)
        self.assertIs(vals["ltest"], True)
        self.assertEqual(vals["stl_file"], "geom.stl")
        self.assertEqual(vals["nsv"], 0)  # scalar default present

class TestNcDataHandle(unittest.TestCase):
    """_load_ncdata must release the NetCDF file handle: load-and-close for a
    requested variable, load-into-memory for browse. On Windows an open handle
    blocks os.remove, so deleting the file afterwards catches a leak."""

    @staticmethod
    def _make_nc(path):
        import xarray as xr

        ds = xr.Dataset(
            {"a": (("time", "x"), np.arange(6.0).reshape(3, 2))},
            coords={"time": [0, 1, 2], "x": [0, 1]},
        )
        ds.to_netcdf(path)
        ds.close()

    def test_requested_variable_releases_handle(self):
        with TemporaryDirectory() as d:
            p = Path(d) / "f.nc"
            self._make_nc(p)
            sim = UDBase.__new__(UDBase)
            arr = sim._load_ncdata(p, "a")
            self.assertEqual(arr.shape, (2, 3))  # transposed (time, x) -> (x, time)
            os.remove(p)  # PermissionError on Windows if the handle leaked

    def test_browse_releases_handle(self):
        import io
        from contextlib import redirect_stdout

        with TemporaryDirectory() as d:
            p = Path(d) / "f.nc"
            self._make_nc(p)
            sim = UDBase.__new__(UDBase)
            with redirect_stdout(io.StringIO()):
                ds = sim._load_ncdata(p, None)
            self.assertIn("a", ds)
            os.remove(p)  # deletable => no lingering handle
