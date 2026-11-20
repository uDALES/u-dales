"""Warm-start facet-output merging without loading whole Paris arrays."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import types
import unittest

from netCDF4 import Dataset
import numpy as np

from udcomf import UDComf
from udcomf.sources import load_facets, load_statistics
from udcomf.udcomf_io import merge_facet_outputs, merge_statistics_outputs


def _write_segment(directory, case, times, *, facets=2, with_temperature=False):
    directory.mkdir(exist_ok=True)
    with Dataset(directory / f"facEB.{case}.nc", "w") as ds:
        ds.createDimension("time", None)
        ds.createDimension("fct", facets)
        ds.createVariable("fct", "i4", ("fct",))
        ds.createVariable("t", "f4", ("time",))[:] = times
        for name in ("LWout", "netsw", "LWin", "hf"):
            var = ds.createVariable(name, "f4", ("time", "fct"), fill_value=-999.)
            var[:] = np.asarray(times)[:, None] + np.arange(facets)[None, :]
            var.units = "W m-2"
    if with_temperature:
        with Dataset(directory / f"facT.{case}.nc", "w") as ds:
            ds.createDimension("time", None)
            ds.createDimension("fct", facets)
            ds.createVariable("t", "f4", ("time",))[:] = times


class TestFacetOutputMerge(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.destination = self.root / "305"
        self.destination.mkdir()
        (self.destination / "namoptions.305").write_text(
            "&WALLS\n nfcts = 2\n nfaclyrs = 1\n/\n", encoding="ascii"
        )

    def test_merges_only_t_and_lwout_in_time_order(self):
        first = self.root / "312"
        _write_segment(first, 312, [0, 300])
        _write_segment(self.destination, 305, [650, 950])   # the case's own records
        sim = types.SimpleNamespace(path=self.destination, expnr="305")

        result = UDComf(sim).prepare_facet_outputs([self.destination, first])

        self.assertEqual(result.path.name, "merged_facEB.305.nc")
        self.assertEqual(result.source_cases, ("312", "305"))
        self.assertEqual(result.record_count, 4)
        self.assertEqual(result.boundary_gaps, (350.,))
        with Dataset(result.path) as ds:
            self.assertEqual(set(ds.variables), {"t", "LWout"})
            self.assertEqual(set(ds.dimensions), {"time", "fct"})
            np.testing.assert_array_equal(ds["t"][:], [0, 300, 650, 950])
            np.testing.assert_array_equal(ds["LWout"][:, 1], [1, 301, 651, 951])
            self.assertEqual(ds["LWout"].units, "W m-2")
            self.assertEqual(json.loads(ds.getncattr("udcomf_merged_cases")), ["312", "305"])
            self.assertEqual(json.loads(ds.getncattr("udcomf_boundary_gaps_seconds")), [350.])
        # the sources are untouched, and the merged file is readable as a facet source
        with Dataset(self.destination / "facEB.305.nc") as ds:
            self.assertIn("netsw", ds.variables)
        np.testing.assert_array_equal(load_facets(sim, "merged_facEB", "t"), [0, 300, 650, 950])
        self.assertEqual(load_facets(sim, "merged_facEB", "LWout", time_index=2).shape, (2,))
        # merged files are rebuilt on every call
        again = merge_facet_outputs(self.destination, [first])
        self.assertEqual(again.record_count, 2)

    def test_does_not_need_fact(self):
        first = self.root / "312"
        _write_segment(first, 312, [0, 300], with_temperature=False)
        result = merge_facet_outputs(self.destination, [first])
        self.assertFalse((self.destination / "facT.305.nc").exists())
        self.assertEqual(result.record_count, 2)

    def test_rejects_overlapping_times_before_writing(self):
        first = self.root / "312"
        second = self.root / "313"
        _write_segment(first, 312, [0, 300])
        _write_segment(second, 313, [300, 600])
        with self.assertRaisesRegex(ValueError, "overlap"):
            merge_facet_outputs(self.destination, [first, second])
        self.assertFalse((self.destination / "merged_facEB.305.nc").exists())

    def test_rejects_incompatible_facet_count(self):
        first = self.root / "312"
        _write_segment(first, 312, [0, 300], facets=3)
        with self.assertRaisesRegex(ValueError, "Facet dimension"):
            merge_facet_outputs(self.destination, [first])

    def test_rejects_incompatible_variable_schema(self):
        first = self.root / "312"
        second = self.root / "313"
        _write_segment(first, 312, [0, 300])
        _write_segment(second, 313, [600, 900])
        with Dataset(second / "facEB.313.nc", "a") as ds:
            ds.renameVariable("LWout", "LWout_old")
            ds.createVariable("LWout", "f8", ("time", "fct"))[:] = 1.0
        with self.assertRaisesRegex(ValueError, "schema"):
            merge_facet_outputs(self.destination, [first, second])
        self.assertFalse((self.destination / "merged_facEB.305.nc").exists())

    def test_merged_file_cannot_be_a_source(self):
        first = self.root / "312"
        _write_segment(first, 312, [0, 300])
        result = merge_facet_outputs(self.destination, [first])
        with self.assertRaises(ValueError):
            merge_facet_outputs(self.destination, [first, result.path])


def _write_stats_segment(directory, case, times, *, levels=(0.75, 2.25, 5.25, 9.75),
                         heights=(1.1, 1.5), name="stats_tkslice", tstatsdump=900.):
    directory.mkdir(exist_ok=True)
    (directory / f"namoptions.{case}").write_text(
        f"&OUTPUT\n tstatsdump={tstatsdump}\n/\n", encoding="ascii"
    )
    times = np.asarray(times, dtype=float)
    with Dataset(directory / f"{name}.{case}.nc", "w") as ds:
        ds.createDimension("time", None)
        for dim, length in (("zt", len(levels)), ("yt", 2), ("xt", 3), ("receptor_height", len(heights))):
            ds.createDimension(dim, length)
        ds.createVariable("time", "f8", ("time",))[:] = times
        ds.createVariable("zt", "f4", ("zt",))[:] = levels
        ds.createVariable("yt", "f4", ("yt",))[:] = [5., 15.]
        ds.createVariable("xt", "f4", ("xt",))[:] = [5., 15., 25.]
        ds.createVariable("receptor_height", "f4", ("receptor_height",))[:] = heights
        for field in ("tha", "pabs", "qt", "ql"):
            var = ds.createVariable(field, "f4", ("time", "zt", "yt", "xt"), fill_value=-999.)
            var.units = "K"
            var[:] = (times[:, None, None, None] + 10 * np.arange(len(levels))[None, :, None, None]
                      + np.arange(2)[None, None, :, None] * 0.1 + np.arange(3)[None, None, None, :] * 0.01)
        ds.createVariable("ws_local", "f4", ("time", "receptor_height", "yt", "xt"))[:] = (
            times[:, None, None, None] + np.arange(len(heights))[None, :, None, None]
        )
        ds.createVariable("ws_10", "f4", ("time", "yt", "xt"))[:] = times[:, None, None] * 2
        ds.createVariable("u", "f4", ("time", "zt", "yt", "xt"))[:] = 7.  # not needed; must be dropped
        ds.title = "fixture"


class TestStatisticsOutputMerge(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.destination = self.root / "305"
        self.destination.mkdir()
        (self.destination / "namoptions.305").write_text(
            "&OUTPUT\n receptor_height=1.1\n nreceptor_heights=2\n receptor_heights=1.1,1.5\n"
            " tstatsdump=900.\n tstatsgap=2700.\n/\n", encoding="ascii"
        )
        (self.destination / "prof.inp.305").write_text(
            "# z thl\n0.75 300\n2.25 300\n5.25 300\n9.75 300\n", encoding="ascii"
        )

    def test_merges_subset_levels_in_time_order_and_records_segments(self):
        _write_stats_segment(self.root / "312", 312, [3600., 7200.])
        _write_stats_segment(self.root / "313", 313, [14400., 18000.])
        result = merge_statistics_outputs(
            self.destination, [self.root / "313", self.root / "312"], source="stats_tkslice")
        self.assertEqual(result.source_cases, ("312", "313"))
        self.assertEqual(result.record_count, 4)
        self.assertEqual(result.selected_levels, (0, 1))
        self.assertEqual(result.boundary_gaps, (14400. - 900. - 7200.,))
        with Dataset(result.path) as ds:
            ds.set_auto_mask(False)
            self.assertEqual(result.path.name, "merged_stats_tkslice.305.nc")
            np.testing.assert_array_equal(ds["time"][:], [3600., 7200., 14400., 18000.])
            np.testing.assert_allclose(ds["zt"][:], [0.75, 2.25])
            self.assertNotIn("u", ds.variables)
            self.assertEqual(ds["tha"].dimensions, ("time", "zt", "yt", "xt"))
            self.assertEqual(ds["tha"].units, "K")
            np.testing.assert_allclose(ds["tha"][2, 1, 1, 2], 14400. + 10 + 0.1 + 0.02, rtol=1e-6)
            np.testing.assert_allclose(ds["ws_local"][3, 1, 0, 0], 18001.)
            np.testing.assert_allclose(ds["ws_10"][0, 0, 0], 7200.)
            self.assertEqual(json.loads(ds.getncattr("udcomf_merged_cases")), ["312", "313"])
            self.assertEqual(json.loads(ds.getncattr("udcomf_segment_first_times")), [3600., 14400.])
            self.assertEqual(json.loads(ds.getncattr("udcomf_selected_scalar_levels")), [1, 2])
            self.assertEqual(ds.getncattr("udcomf_averaging_window_s"), 900.)
            self.assertEqual(ds.title, "fixture")
        # merged files are rebuilt on every call
        again = merge_statistics_outputs(self.destination, [self.root / "312"], source="stats_tkslice")
        self.assertEqual(again.record_count, 2)

    def test_reads_only_the_requested_file_type_without_fallback(self):
        _write_stats_segment(self.root / "312", 312, [3600.])                   # stats_tkslice only
        _write_stats_segment(self.root / "313", 313, [7200.], name="stats_t")  # stats_t only
        with self.assertRaisesRegex(FileNotFoundError, "stats_t file is missing for case 312"):
            merge_statistics_outputs(self.destination, [self.root / "312", self.root / "313"])
        with self.assertRaisesRegex(FileNotFoundError, "stats_tkslice file is missing for case 313"):
            merge_statistics_outputs(self.destination, [self.root / "312", self.root / "313"],
                                     source="stats_tkslice")
        result = merge_statistics_outputs(self.destination, [self.root / "313"])
        self.assertEqual(result.source_files[0].name, "stats_t.313.nc")
        self.assertEqual(result.path.name, "merged_stats_tkslice.305.nc")
        with self.assertRaises(ValueError):
            merge_statistics_outputs(self.destination, [self.root / "313"], source="stats_x")

    def test_rejects_overlapping_windows_and_mismatched_tstatsdump(self):
        kw = dict(source="stats_tkslice")
        _write_stats_segment(self.root / "312", 312, [3600., 7200.])
        _write_stats_segment(self.root / "313", 313, [7900., 8800.])
        with self.assertRaisesRegex(ValueError, "windows overlap"):
            merge_statistics_outputs(self.destination, [self.root / "312", self.root / "313"], **kw)
        _write_stats_segment(self.root / "314", 314, [14400.], tstatsdump=600.)
        with self.assertRaisesRegex(ValueError, "tstatsdump differs"):
            merge_statistics_outputs(self.destination, [self.root / "312", self.root / "314"], **kw)
        self.assertFalse((self.destination / "merged_stats_tkslice.305.nc").exists())

    def test_rejects_missing_bracketing_level_and_merged_file_as_source(self):
        kw = dict(source="stats_tkslice")
        _write_stats_segment(self.root / "312", 312, [3600.], levels=(2.25, 5.25))
        with self.assertRaisesRegex(ValueError, "lacks the scalar level"):
            merge_statistics_outputs(self.destination, [self.root / "312"], **kw)
        _write_stats_segment(self.destination, 305, [3600.])
        # the case's own statistics may be a segment; the merged output may not
        result = merge_statistics_outputs(self.destination, [self.destination], **kw)
        self.assertEqual(result.source_cases, ("305",))
        with self.assertRaises(ValueError):   # merged_* is not a stats_tkslice source name
            merge_statistics_outputs(self.destination, [result.path], **kw)

    def test_facade_entry_point_and_merged_source_loader(self):
        _write_stats_segment(self.root / "312", 312, [3600.])
        sim = types.SimpleNamespace(path=self.destination, expnr="305")
        result = UDComf(sim).prepare_statistics_outputs([self.root / "312"], source="stats_tkslice")
        self.assertEqual(result.record_count, 1)
        self.assertEqual(result.path.name, "merged_stats_tkslice.305.nc")
        np.testing.assert_array_equal(load_statistics(sim, "merged_stats_tkslice", "time"), [3600.])
        self.assertEqual(load_statistics(sim, "merged_stats_tkslice", "ws_10", time_index=0).shape, (3, 2))
        with self.assertRaises(ValueError):
            load_statistics(sim, "stats_x", "time")
