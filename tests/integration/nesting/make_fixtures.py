#!/usr/bin/env python3
"""Generate the nesting unit-test fixtures (docs/udales-nesting-design.md section 10.1).

Everything the in-solver tests read is produced here, through the *production*
writer in tools/python/udprep/nesting.py, so the Python writer and the Fortran
reader cannot drift apart (tests U15/U16 close that seam).

Files written into the target directory:

  nesting_analytic.<expnr>.nc    exactly analytic_field(), NOT divergence
                                 corrected, so every stored value equals the
                                 analytic function to the last bit. Used by
                                 U15-U21, U23-U25, U28-U34 (with nest_fluxtol
                                 raised so the init-time flux check passes).
  nesting_corrected.<expnr>.nc   the same field with the offline divergence
                                 correction applied (U26).
  assertfire.<expnr>.nc          a copy of the uncorrected file, used with the
                                 production tolerance so the flux assertion
                                 must fire (U27).
  nesting_nonlinear.<expnr>.nc   analytic in space, NON-linear in time (the
                                 stored level n holds the analytic field at
                                 pseudo-time s_n = t_n^2 / t_max). Needed by
                                 U18/U19: a field linear in t is reproduced
                                 exactly by every interpolant, so it cannot
                                 distinguish linear from Hermite.
  nesting_initial.<expnr>.nc     the analytic field with the schema 2 FULL-DOMAIN
                                 initial-condition block, uncorrected so the block
                                 is exactly analytic_field() (runmode 1011,
                                 U35-U37 and the cold-start abort cases).
  nesting_v1.<expnr>.nc          the corrected field written as SCHEMA 1: no
                                 flux_residual, no initial condition. The
                                 backwards-compatibility fixture (U38).
  assertfire_lying.<expnr>.nc    schema 2, data NOT flux balanced, but with
                                 flux_residual overwritten with zeros: the cheap
                                 init check believes it, nest_lfluxcheckall does
                                 not (U39 and an abort case).
  assertfire_area.<expnr>.nc     the same, plus a corrupted fluid_lateral_area, so
                                 the reader must distrust the stored residual and
                                 fall back to the full recompute -- which aborts.
  bad_schema.<expnr>.nc          header-validation fixtures for U22, each with
  bad_itot.<expnr>.nc            exactly one corrupted header field.
  bad_xlen.<expnr>.nc
  bad_zf.<expnr>.nc
  bad_stagger.<expnr>.nc
  bad_initdims.<expnr>.nc        an initial-condition block at the wrong shape and
  bad_initstag.<expnr>.nc        at the wrong stagger (runmode 1011 abort cases).
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np

from udprep.nesting import (
    ANALYTIC_COEFFS,
    NestGrid,
    NestingData,
    analytic_initial_fields,
    analytic_slabs,
    apply_divergence_correction,
    net_volume_flux,
    write_nesting_file,
)

# Fixture geometry. Kept small so the runmodes finish in a second, but large
# enough that a 2x2 decomposition of a 12-cell zone is non-trivial.
EXPNR = "901"
ITOT = 32
JTOT = 32
KTOT = 16
XLEN = 32.0
YLEN = 32.0
ZSIZE = 16.0
NZONE = 12
TIMES = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0])


def grid() -> NestGrid:
    return NestGrid.uniform(ITOT, JTOT, KTOT, XLEN, YLEN, ZSIZE)


def _data(times, slab_times=None, corrected=False, initial=False) -> NestingData:
    g = grid()
    slab_times = times if slab_times is None else slab_times
    data = NestingData(
        grid=g,
        nzone=NZONE,
        times=np.asarray(times, dtype=np.float64),
        slabs=analytic_slabs(g, NZONE, slab_times),
        initial_fields=(analytic_initial_fields(g, float(np.asarray(slab_times)[0]))
                        if initial else None),
        parent_model="analytic",
        parent_dx=float(XLEN / ITOT),
        parent_dt=float(np.min(np.diff(times))) if len(times) > 1 else 0.0,
        child_origin_x=0.0,
        child_origin_y=0.0,
    )
    if corrected:
        apply_divergence_correction(data)
    else:
        data.net_volume_flux = net_volume_flux(data)
        data.flux_residual = data.net_volume_flux.copy()
    return data


def _corrupt(src: Path, dst: Path, **edits) -> None:
    """Copy ``src`` to ``dst`` and rewrite the named header items."""
    import netCDF4 as nc

    shutil.copy2(src, dst)
    with nc.Dataset(dst, "a") as ds:
        for key, value in edits.items():
            if key in ("zf", "flux_residual"):
                ds.variables[key][:] = value
            elif key == "stagger":
                var, tag = value
                ds.variables[var].stagger = tag
            else:
                ds.setncattr(key, value)


def _copy_with_dims(src: Path, dst: Path, varname: str, dims) -> None:
    """Copy a NetCDF file verbatim, redefining one variable's dimensions.

    netCDF4 cannot change a variable's shape in place, so the wrong-shape
    fixture has to be rebuilt.  Everything except ``varname`` is copied as it
    stands, so the file differs from a good one in exactly one respect.
    """
    import netCDF4 as nc

    with nc.Dataset(src, "r") as a, nc.Dataset(dst, "w", format="NETCDF4") as b:
        for name, dim in a.dimensions.items():
            b.createDimension(name, None if dim.isunlimited() else len(dim))
        for name, var in a.variables.items():
            newdims = tuple(dims) if name == varname else var.dimensions
            out = b.createVariable(name, var.datatype, newdims)
            out.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            if name == varname:
                shape = tuple(len(b.dimensions[d]) for d in newdims)
                values = np.asarray(var[:], dtype=np.float64)
                out[:] = np.resize(values, shape)
            else:
                out[:] = var[:]
        b.setncatts({k: a.getncattr(k) for k in a.ncattrs()})


def write_all(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    analytic = outdir / f"nesting_analytic.{EXPNR}.nc"
    write_nesting_file(analytic, _data(TIMES), override=True)

    write_nesting_file(
        outdir / f"nesting_corrected.{EXPNR}.nc", _data(TIMES, corrected=True), override=True
    )

    shutil.copy2(analytic, outdir / f"assertfire.{EXPNR}.nc")

    # Non-linear in time: stored level n is the analytic field at s_n, but the
    # file says the level is at t_n. The Fortran side knows both, so it can
    # predict every node value exactly while d(target)/dt is genuinely
    # discontinuous under linear interpolation.
    pseudo = TIMES**2 / TIMES[-1]
    write_nesting_file(
        outdir / f"nesting_nonlinear.{EXPNR}.nc",
        _data(TIMES, slab_times=pseudo),
        override=True,
    )

    # --- schema 2: the full-domain initial condition (runmode 1011) ---
    initial = outdir / f"nesting_initial.{EXPNR}.nc"
    write_nesting_file(initial, _data(TIMES, initial=True), override=True)

    # --- schema 1: the backwards-compatibility fixture (U38) ---
    write_nesting_file(
        outdir / f"nesting_v1.{EXPNR}.nc", _data(TIMES, corrected=True),
        override=True, schema=1,
    )

    # --- a writer that lies about its own residual (U39 and two abort cases) ---
    for name, area in ((f"assertfire_lying.{EXPNR}.nc", None),
                       (f"assertfire_area.{EXPNR}.nc", 1.0)):
        edits = {"flux_residual": np.zeros(len(TIMES))}
        if area is not None:
            edits["fluid_lateral_area"] = np.float64(area)
        _corrupt(analytic, outdir / name, **edits)

    zf_bad = grid().zf.copy()
    zf_bad[KTOT // 2] += 0.25
    _corrupt(analytic, outdir / f"bad_schema.{EXPNR}.nc", udales_nesting_schema=np.int32(7))
    _corrupt(analytic, outdir / f"bad_itot.{EXPNR}.nc", itot=np.int32(ITOT + 1))
    _corrupt(analytic, outdir / f"bad_xlen.{EXPNR}.nc", xlen=np.float64(XLEN * 1.5))
    _corrupt(analytic, outdir / f"bad_zf.{EXPNR}.nc", zf=zf_bad)
    _corrupt(analytic, outdir / f"bad_stagger.{EXPNR}.nc", stagger=("u_west", "xf yf zf"))
    _corrupt(initial, outdir / f"bad_initstag.{EXPNR}.nc", stagger=("u_init", "xf yf zf"))
    _copy_with_dims(initial, outdir / f"bad_initdims.{EXPNR}.nc",
                    "u_init", ("xf", "yf", "zf"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="directory to write the fixtures into")
    args = parser.parse_args()
    write_all(args.outdir)
    print(f"nesting fixtures written to {args.outdir}")
    print(f"  coefficients: {ANALYTIC_COEFFS}")


if __name__ == "__main__":
    main()
