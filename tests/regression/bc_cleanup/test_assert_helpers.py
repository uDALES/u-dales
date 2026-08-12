#!/usr/bin/env python3

"""Standalone unit tests for the fill-value-agnostic sentinel/NaN assertion
helpers in run_test.py (`_var_fill_value`, `_marker_mask`,
`_assert_sentinels_only_in_solid_slabs`, `_assert_no_sentinel_or_nan`).

These exercise the helpers directly against tiny, hand-built NetCDF fixtures
(no solver build, no case fixtures under tests/cases) so the fill-value
plumbing can be verified in isolation from the full regression harness that
run_test.py otherwise requires (MPI, a built u-dales executable, git
worktrees). Run with:

    python -m pytest tests/regression/bc_cleanup/test_assert_helpers.py -q

run_test.py is a script, not a package (no __init__.py in this directory),
so it is loaded via importlib.util.spec_from_file_location rather than a
normal import.
"""

import importlib.util
import math
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pytest

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


# ---------------------------------------------------------------------------
# _var_fill_value
# ---------------------------------------------------------------------------


def test_var_fill_value_reads_declared_attribute(tmp_path):
    path = tmp_path / "one.nc"
    _write_nc_var(path, n_zt=3, values=np.zeros((1, 3)), fill=-1.0e30)
    with nc.Dataset(path) as ds:
        assert rt._var_fill_value(ds.variables["thl"]) == pytest.approx(-1.0e30)


def test_var_fill_value_reads_nan_attribute(tmp_path):
    path = tmp_path / "one.nc"
    _write_nc_var(path, n_zt=3, values=np.zeros((1, 3)), fill=float("nan"))
    with nc.Dataset(path) as ds:
        assert math.isnan(rt._var_fill_value(ds.variables["thl"]))


def test_var_fill_value_falls_back_to_sentinel_when_attribute_absent(tmp_path):
    path = tmp_path / "one.nc"
    _write_nc_var(path, n_zt=3, values=np.zeros((1, 3)), fill=None)
    with nc.Dataset(path) as ds:
        var = ds.variables["thl"]
        assert not hasattr(var, "_FillValue")
        assert rt._var_fill_value(var) == rt.SENTINEL_VALUE


# ---------------------------------------------------------------------------
# _marker_mask
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fill", FILL_VALUES, ids=_fill_id)
def test_marker_mask_matches_only_the_fill_value(fill):
    arr = np.array([1.0, fill, 3.0, fill])
    mask = rt._marker_mask(arr, fill)
    np.testing.assert_array_equal(mask, [False, True, False, True])


# ---------------------------------------------------------------------------
# _assert_sentinels_only_in_solid_slabs
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fill", FILL_VALUES, ids=_fill_id)
def test_sentinels_only_in_solid_slabs_passes_when_marker_at_solid_levels_only(tmp_path, fill):
    itot = jtot = 2
    n_zt = 4
    _write_namoptions(tmp_path, "unit", itot, jtot)
    _write_solid_c(tmp_path, itot, jtot, solid_k_levels=(1, 2))

    # zt indices are 0-based; k=1,2 (solid) -> zt index 0,1. k=3,4 (fluid) ->
    # zt index 2,3, which must carry ordinary finite data, not the marker.
    values = np.array([[fill, fill, 5.0, 6.0]])
    _write_nc_var(tmp_path / "xytdump.000.000.unit.nc", n_zt=n_zt, values=values, fill=fill)

    spec = rt.CaseSpec(case="unit", nc_pattern="tdump.*.unit.nc", fields=None, abs_tol=None, default_atol=0.0)
    rt._assert_sentinels_only_in_solid_slabs(tmp_path, spec)  # must not raise


@pytest.mark.parametrize("fill", FILL_VALUES, ids=_fill_id)
def test_sentinels_only_in_solid_slabs_raises_on_fluid_level_marker(tmp_path, fill):
    itot = jtot = 2
    n_zt = 4
    _write_namoptions(tmp_path, "unit", itot, jtot)
    _write_solid_c(tmp_path, itot, jtot, solid_k_levels=(1, 2))

    # Marker at zt index 2 (k=3) is a fluid level -- must be flagged.
    values = np.array([[fill, fill, fill, 6.0]])
    _write_nc_var(tmp_path / "xytdump.000.000.unit.nc", n_zt=n_zt, values=values, fill=fill)

    spec = rt.CaseSpec(case="unit", nc_pattern="tdump.*.unit.nc", fields=None, abs_tol=None, default_atol=0.0)
    with pytest.raises(RuntimeError, match="fluid"):
        rt._assert_sentinels_only_in_solid_slabs(tmp_path, spec)


@pytest.mark.parametrize("fill", FILL_VALUES, ids=_fill_id)
def test_sentinels_only_in_solid_slabs_raises_for_tdump_with_any_marker(tmp_path, fill):
    itot = jtot = 2
    n_zt = 4
    _write_namoptions(tmp_path, "unit", itot, jtot)
    _write_solid_c(tmp_path, itot, jtot, solid_k_levels=(1, 2))

    # tdump (not xytdump) must carry no marker at all, even at a genuinely
    # fully-solid level.
    values = np.array([[fill, 2.0, 3.0, 4.0]])
    _write_nc_var(tmp_path / "tdump.000.000.unit.nc", n_zt=n_zt, values=values, fill=fill)

    spec = rt.CaseSpec(case="unit", nc_pattern="tdump.*.unit.nc", fields=None, abs_tol=None, default_atol=0.0)
    with pytest.raises(RuntimeError, match="nodata marker"):
        rt._assert_sentinels_only_in_solid_slabs(tmp_path, spec)


# ---------------------------------------------------------------------------
# _assert_no_sentinel_or_nan
# ---------------------------------------------------------------------------


def _spec_for_dumps() -> "rt.CaseSpec":
    return rt.CaseSpec(case="unit", nc_pattern="tdump.*.unit.nc", fields=None, abs_tol=None, default_atol=0.0)


@pytest.mark.parametrize("fill", FILL_VALUES, ids=_fill_id)
def test_no_sentinel_or_nan_passes_on_clean_file(tmp_path, fill):
    n_zt = 3
    values = np.array([[[1.0], [2.0], [3.0]]])  # (time, zt, y), no marker present
    _write_nc_var(
        tmp_path / "tdump.000.000.unit.nc",
        n_zt=n_zt,
        values=values,
        fill=fill,
        var_dims=("time", "zt", "y"),
        extra_dims={"y": 1},
    )
    rt._assert_no_sentinel_or_nan(tmp_path, _spec_for_dumps())  # must not raise


@pytest.mark.parametrize("fill", FILL_VALUES, ids=_fill_id)
def test_no_sentinel_or_nan_raises_on_fill_marker(tmp_path, fill):
    n_zt = 3
    values = np.array([[[1.0], [fill], [3.0]]])
    _write_nc_var(
        tmp_path / "tdump.000.000.unit.nc",
        n_zt=n_zt,
        values=values,
        fill=fill,
        var_dims=("time", "zt", "y"),
        extra_dims={"y": 1},
    )
    with pytest.raises(RuntimeError, match="sentinel|NaN|nodata"):
        rt._assert_no_sentinel_or_nan(tmp_path, _spec_for_dumps())


def test_no_sentinel_or_nan_raises_on_stray_nan_when_fill_is_finite(tmp_path):
    # Fill value is finite (-999.0): a NaN elsewhere is never the declared
    # marker and must always be reported as an error in its own right.
    n_zt = 3
    fill = -999.0
    values = np.array([[[1.0], [float("nan")], [3.0]]])
    _write_nc_var(
        tmp_path / "tdump.000.000.unit.nc",
        n_zt=n_zt,
        values=values,
        fill=fill,
        var_dims=("time", "zt", "y"),
        extra_dims={"y": 1},
    )
    with pytest.raises(RuntimeError, match="NaN"):
        rt._assert_no_sentinel_or_nan(tmp_path, _spec_for_dumps())
