"""Warm-start facet-output merging without loading whole Paris arrays."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import types
import unittest

from netCDF4 import Dataset
import numpy as np

from udcomf import UDComf
from udcomf.udcomf_io import merge_facet_outputs


def _write_segment(directory, case, times, *, facets=2, layers=2):
    directory.mkdir(exist_ok=True)
    for kind in ("facEB", "facT"):
        with Dataset(directory / f"{kind}.{case}.nc", "w") as ds:
            ds.createDimension("time", None)
            ds.createDimension("fct", facets)
            if kind == "facT":
                ds.createDimension("lyr", layers)
            ds.createVariable("fct", "i4", ("fct",))
            ds.createVariable("t", "f4", ("time",))[:] = times
            if kind == "facEB":
                for name in ("LWout", "netsw"):
                    ds.createVariable(name, "f4", ("time", "fct"), fill_value=-999.)[:] = (
                        np.asarray(times)[:, None] + np.arange(facets)[None, :]
                    )
            else:
                ds.createVariable("lyr", "i4", ("lyr",))
                for name in ("T", "dTdz"):
                    ds.createVariable(name, "f4", ("time", "lyr", "fct"))[:] = (
                        np.asarray(times)[:, None, None]
                        + np.arange(layers)[None, :, None]
                        + np.arange(facets)[None, None, :]
                    )


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

    def test_merges_directories_and_copied_files_in_time_order(self):
        first = self.root / "312"
        _write_segment(first, 312, [0, 300])
        _write_segment(self.destination, 313, [650, 950])
        copied_eb = self.destination / "facEB.313.nc"
        sim = types.SimpleNamespace(path=self.destination)

        result = UDComf(sim).prepare_facet_outputs([copied_eb, first])

        self.assertEqual(result.source_cases, ("312", "313"))
        self.assertEqual(result.record_count, 4)
        self.assertEqual(result.boundary_gaps, (350.,))
        for kind in ("facEB", "facT"):
            with Dataset(self.destination / f"{kind}.305.nc") as ds:
                np.testing.assert_array_equal(ds.variables["t"][:], [0, 300, 650, 950])
                self.assertEqual(json.loads(ds.getncattr("udcomf_merged_cases")), ["312", "313"])
                self.assertEqual(len(ds.dimensions["fct"]), 2)
                np.testing.assert_array_equal(ds.variables["fct"][:].data,
                                              [-2147483647, -2147483647])
                if kind == "facEB":
                    np.testing.assert_array_equal(ds.variables["LWout"][2, :], [650, 651])
                else:
                    np.testing.assert_array_equal(ds.variables["T"][3, :, :],
                                                  [[950, 951], [951, 952]])
        with Dataset(copied_eb) as source:
            self.assertEqual(len(source.dimensions["time"]), 2)
        with self.assertRaises(FileExistsError):
            merge_facet_outputs(self.destination, [first, copied_eb])

    def test_rejects_overlapping_times_before_writing(self):
        first = self.root / "312"
        second = self.root / "313"
        _write_segment(first, 312, [0, 300])
        _write_segment(second, 313, [300, 600])
        with self.assertRaisesRegex(ValueError, "overlap"):
            merge_facet_outputs(self.destination, [first, second])
        self.assertFalse((self.destination / "facEB.305.nc").exists())

    def test_rejects_incompatible_facet_count(self):
        first = self.root / "312"
        _write_segment(first, 312, [0], facets=3)
        with self.assertRaisesRegex(ValueError, "Facet or layer dimensions"):
            merge_facet_outputs(self.destination, [first])

    def test_rejects_incompatible_variable_schema(self):
        first = self.root / "312"
        second = self.root / "313"
        _write_segment(first, 312, [0])
        _write_segment(second, 313, [300])
        with Dataset(second / "facEB.313.nc", "a") as ds:
            ds.createVariable("extra", "f4", ("time", "fct"))[:] = [[1, 2]]
        with self.assertRaisesRegex(ValueError, "Incompatible eb_path NetCDF schema"):
            merge_facet_outputs(self.destination, [first, second])
        self.assertFalse((self.destination / "facEB.305.nc").exists())

    def test_rejects_unpaired_temperature_times(self):
        first = self.root / "312"
        _write_segment(first, 312, [0, 300])
        with Dataset(first / "facT.312.nc", "a") as ds:
            ds.variables["t"][1] = 301
        with self.assertRaisesRegex(ValueError, "timestamps must match"):
            merge_facet_outputs(self.destination, [first])


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
        from udcomf.udcomf_io import merge_statistics_outputs
        _write_stats_segment(self.root / "312", 312, [3600., 7200.])
        _write_stats_segment(self.root / "313", 313, [14400., 18000.], name="stats_t")
        result = merge_statistics_outputs(self.destination, [self.root / "313", self.root / "312"])
        self.assertEqual(result.source_cases, ("312", "313"))
        self.assertEqual(result.record_count, 4)
        self.assertEqual(result.selected_levels, (0, 1))
        self.assertEqual(result.boundary_gaps, (14400. - 900. - 7200.,))
        with Dataset(result.path) as ds:
            ds.set_auto_mask(False)
            self.assertEqual(result.path.name, "stats_tkslice.305.nc")
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
        with self.assertRaises(FileExistsError):
            merge_statistics_outputs(self.destination, [self.root / "312", self.root / "313"])

    def test_rejects_overlapping_windows_and_mismatched_tstatsdump(self):
        from udcomf.udcomf_io import merge_statistics_outputs
        _write_stats_segment(self.root / "312", 312, [3600., 7200.])
        _write_stats_segment(self.root / "313", 313, [7900., 8800.])
        with self.assertRaisesRegex(ValueError, "windows overlap"):
            merge_statistics_outputs(self.destination, [self.root / "312", self.root / "313"])
        _write_stats_segment(self.root / "314", 314, [14400.], tstatsdump=600.)
        with self.assertRaisesRegex(ValueError, "tstatsdump differs"):
            merge_statistics_outputs(self.destination, [self.root / "312", self.root / "314"])
        self.assertFalse((self.destination / "stats_tkslice.305.nc").exists())

    def test_rejects_missing_bracketing_level_and_destination_as_source(self):
        from udcomf.udcomf_io import merge_statistics_outputs
        _write_stats_segment(self.root / "312", 312, [3600.], levels=(2.25, 5.25))
        with self.assertRaisesRegex(ValueError, "lacks the scalar level"):
            merge_statistics_outputs(self.destination, [self.root / "312"])
        _write_stats_segment(self.destination, 305, [3600.])
        with self.assertRaisesRegex(ValueError, "cannot also be a warm-start source"):
            merge_statistics_outputs(self.destination, [self.destination / "stats_tkslice.305.nc"])

    def test_facade_entry_point(self):
        _write_stats_segment(self.root / "312", 312, [3600.])
        sim = types.SimpleNamespace(path=self.destination)
        result = UDComf(sim).prepare_statistics_outputs([self.root / "312"])
        self.assertEqual(result.record_count, 1)
