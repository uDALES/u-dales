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
