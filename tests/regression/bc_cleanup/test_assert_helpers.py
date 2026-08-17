#!/usr/bin/env python3

"""Standalone unit tests for the fill-value-agnostic sentinel/NaN assertion
helpers in run_test.py (`_var_fill_value`, `_marker_mask`,
`_assert_sentinels_only_in_solid_slabs`, `_assert_no_sentinel_or_nan`).

These exercise the helpers directly against tiny, hand-built NetCDF fixtures
(no solver build, no case fixtures under tests/cases) so the fill-value
plumbing can be verified in isolation from the full regression harness that
run_test.py otherwise requires (MPI, a built u-dales executable, git
worktrees).

Written against stdlib unittest (the CI python environment has no pytest);
sub-cases over the fill values use subTest. Run with either:

    python -m unittest discover -s tests/regression/bc_cleanup -p "test_*.py"
    python -m pytest tests/regression/bc_cleanup/test_assert_helpers.py -q

run_test.py is a script, not a package (no __init__.py in this directory),
so it is loaded via importlib.util.spec_from_file_location rather than a
normal import.
"""

import importlib.util
import math
import tempfile
import unittest
from pathlib import Path

import netCDF4 as nc
import numpy as np

RUN_TEST_PATH = Path(__file__).resolve().parent / "run_test.py"


def _load_run_test():
    spec = importlib.util.spec_from_file_location("bc_cleanup_run_test", RUN_TEST_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


rt = _load_run_test()

# Fill values under test: the current Fortran marker (-999.0), a plausible
# alternative finite marker, and a signalling-NaN-style marker -- the whole
# point of keying on `_FillValue` is that the harness must not care which of
# these the file declares.
FILL_VALUES = (-999.0, -1.0e30, float("nan"))


def _fill_id(value: float) -> str:
    return "nan" if math.isnan(value) else str(value)


def _write_namoptions(run_dir: Path, case: str, itot: int, jtot: int) -> None:
    # _fully_solid_levels only regex-matches "key = digits" lines; it does not
    # parse namelist sections, so this minimal text is sufficient.
    (run_dir / f"namoptions.{case}").write_text(
        f"&DOMAIN\nitot = {itot}\njtot = {jtot}\n/\n", encoding="utf-8"
    )


def _write_solid_c(run_dir: Path, itot: int, jtot: int, solid_k_levels) -> None:
    """Fabricate solid_c.txt marking every (i,j) cell solid at each k in
    solid_k_levels -- i.e. those slabs are fully solid, matching the format
    _fully_solid_levels expects: '# position (i,j,k)' header, then 1-based
    i j k rows, one per solid cell centre.
    """
    lines = ["# position (i,j,k)"]
    for k in solid_k_levels:
        for i in range(1, itot + 1):
            for j in range(1, jtot + 1):
                lines.append(f"{i:5d}{j:5d}{k:5d}")
    (run_dir / "solid_c.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_nc_var(
    path: Path,
    *,
    n_zt: int,
    values: np.ndarray,
    fill,
    var_dims=("time", "zt"),
    extra_dims=None,
) -> None:
    """Create a minimal NetCDF file with coordinate vars (time, zt) and a
    single float32 data variable 'thl'.

    `values` must already include the leading time=1 axis, e.g. shape
    (1, n_zt) for var_dims=("time", "zt").

    `fill` is the value to declare as the variable's `_FillValue` (may be
    float('nan')); pass fill=None to create the variable WITHOUT a
    `_FillValue` attribute at all (exercises the _var_fill_value fallback).
    """
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("time", None)
        ds.createDimension("zt", n_zt)
        for dname, size in (extra_dims or {}).items():
            ds.createDimension(dname, size)
        time_var = ds.createVariable("time", "f8", ("time",))
        time_var[:] = [0.0]
        zt_var = ds.createVariable("zt", "f4", ("zt",))
        zt_var[:] = np.arange(1, n_zt + 1, dtype="f4")
        kwargs = {"fill_value": fill} if fill is not None else {"fill_value": False}
        var = ds.createVariable("thl", "f4", var_dims, **kwargs)
        var[:] = values.astype("f4")


def _spec_for_dumps() -> "rt.CaseSpec":
    return rt.CaseSpec(case="unit", nc_pattern="tdump.*.unit.nc", fields=None, abs_tol=None, default_atol=0.0)


class TmpDirTestCase(unittest.TestCase):
    """Base: a fresh temporary directory per test, mirroring pytest's tmp_path."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp_path = Path(self._tmp.name)


class TestVarFillValue(TmpDirTestCase):
    def test_reads_declared_attribute(self):
        path = self.tmp_path / "one.nc"
        _write_nc_var(path, n_zt=3, values=np.zeros((1, 3)), fill=-1.0e30)
        with nc.Dataset(path) as ds:
            got = rt._var_fill_value(ds.variables["thl"])
        # the attribute is stored as float32, so compare with a relative tolerance
        self.assertTrue(math.isclose(got, -1.0e30, rel_tol=1e-6))

    def test_reads_nan_attribute(self):
        path = self.tmp_path / "one.nc"
        _write_nc_var(path, n_zt=3, values=np.zeros((1, 3)), fill=float("nan"))
        with nc.Dataset(path) as ds:
            self.assertTrue(math.isnan(rt._var_fill_value(ds.variables["thl"])))

    def test_falls_back_to_sentinel_when_attribute_absent(self):
        path = self.tmp_path / "one.nc"
        _write_nc_var(path, n_zt=3, values=np.zeros((1, 3)), fill=None)
        with nc.Dataset(path) as ds:
            var = ds.variables["thl"]
            self.assertFalse(hasattr(var, "_FillValue"))
            self.assertEqual(rt._var_fill_value(var), rt.SENTINEL_VALUE)


class TestMarkerMask(unittest.TestCase):
    def test_matches_only_the_fill_value(self):
        for fill in FILL_VALUES:
            with self.subTest(fill=_fill_id(fill)):
                arr = np.array([1.0, fill, 3.0, fill])
                mask = rt._marker_mask(arr, fill)
                np.testing.assert_array_equal(mask, [False, True, False, True])


class TestSentinelsOnlyInSolidSlabs(TmpDirTestCase):
    def test_passes_when_marker_at_solid_levels_only(self):
        for fill in FILL_VALUES:
            with self.subTest(fill=_fill_id(fill)):
                itot = jtot = 2
                n_zt = 4
                _write_namoptions(self.tmp_path, "unit", itot, jtot)
                _write_solid_c(self.tmp_path, itot, jtot, solid_k_levels=(1, 2))

                # zt indices are 0-based; k=1,2 (solid) -> zt index 0,1.
                # k=3,4 (fluid) -> zt index 2,3, which must carry ordinary
                # finite data, not the marker.
                values = np.array([[fill, fill, 5.0, 6.0]])
                _write_nc_var(self.tmp_path / "xytdump.000.000.unit.nc", n_zt=n_zt, values=values, fill=fill)

                rt._assert_sentinels_only_in_solid_slabs(self.tmp_path, _spec_for_dumps())  # must not raise
                (self.tmp_path / "xytdump.000.000.unit.nc").unlink()

    def test_raises_on_fluid_level_marker(self):
        for fill in FILL_VALUES:
            with self.subTest(fill=_fill_id(fill)):
                itot = jtot = 2
                n_zt = 4
                _write_namoptions(self.tmp_path, "unit", itot, jtot)
                _write_solid_c(self.tmp_path, itot, jtot, solid_k_levels=(1, 2))

                # Marker at zt index 2 (k=3) is a fluid level -- must be flagged.
                values = np.array([[fill, fill, fill, 6.0]])
                _write_nc_var(self.tmp_path / "xytdump.000.000.unit.nc", n_zt=n_zt, values=values, fill=fill)

                with self.assertRaisesRegex(RuntimeError, "fluid"):
                    rt._assert_sentinels_only_in_solid_slabs(self.tmp_path, _spec_for_dumps())
                (self.tmp_path / "xytdump.000.000.unit.nc").unlink()

    def test_raises_for_tdump_with_any_marker(self):
        for fill in FILL_VALUES:
            with self.subTest(fill=_fill_id(fill)):
                itot = jtot = 2
                n_zt = 4
                _write_namoptions(self.tmp_path, "unit", itot, jtot)
                _write_solid_c(self.tmp_path, itot, jtot, solid_k_levels=(1, 2))

                # tdump (not xytdump) must carry no marker at all, even at a
                # genuinely fully-solid level.
                values = np.array([[fill, 2.0, 3.0, 4.0]])
                _write_nc_var(self.tmp_path / "tdump.000.000.unit.nc", n_zt=n_zt, values=values, fill=fill)

                with self.assertRaisesRegex(RuntimeError, "nodata marker"):
                    rt._assert_sentinels_only_in_solid_slabs(self.tmp_path, _spec_for_dumps())
                (self.tmp_path / "tdump.000.000.unit.nc").unlink()


class TestNoSentinelOrNan(TmpDirTestCase):
    def test_passes_on_clean_file(self):
        for fill in FILL_VALUES:
            with self.subTest(fill=_fill_id(fill)):
                n_zt = 3
                values = np.array([[[1.0], [2.0], [3.0]]])  # (time, zt, y), no marker present
                path = self.tmp_path / "tdump.000.000.unit.nc"
                _write_nc_var(
                    path,
                    n_zt=n_zt,
                    values=values,
                    fill=fill,
                    var_dims=("time", "zt", "y"),
                    extra_dims={"y": 1},
                )
                rt._assert_no_sentinel_or_nan(self.tmp_path, _spec_for_dumps())  # must not raise
                path.unlink()

    def test_raises_on_fill_marker(self):
        for fill in FILL_VALUES:
            with self.subTest(fill=_fill_id(fill)):
                n_zt = 3
                values = np.array([[[1.0], [fill], [3.0]]])
                path = self.tmp_path / "tdump.000.000.unit.nc"
                _write_nc_var(
                    path,
                    n_zt=n_zt,
                    values=values,
                    fill=fill,
                    var_dims=("time", "zt", "y"),
                    extra_dims={"y": 1},
                )
                with self.assertRaisesRegex(RuntimeError, "sentinel|NaN|nodata"):
                    rt._assert_no_sentinel_or_nan(self.tmp_path, _spec_for_dumps())
                path.unlink()

    def test_raises_on_stray_nan_when_fill_is_finite(self):
        # Fill value is finite (-999.0): a NaN elsewhere is never the declared
        # marker and must always be reported as an error in its own right.
        n_zt = 3
        fill = -999.0
        values = np.array([[[1.0], [float("nan")], [3.0]]])
        _write_nc_var(
            self.tmp_path / "tdump.000.000.unit.nc",
            n_zt=n_zt,
            values=values,
            fill=fill,
            var_dims=("time", "zt", "y"),
            extra_dims={"y": 1},
        )
        with self.assertRaisesRegex(RuntimeError, "NaN"):
            rt._assert_no_sentinel_or_nan(self.tmp_path, _spec_for_dumps())


if __name__ == "__main__":
    unittest.main()
