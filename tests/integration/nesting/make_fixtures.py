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
  nesting_maskwest.<expnr>.nc    the analytic field, divergence corrected over the
                                 FLUID lateral faces only, with the west faces
                                 j = 5..8, k = 3..6 (1-based) declared solid: the
                                 footprint of the 2-cell-deep box U45 stands on
                                 the west boundary. The solid faces keep their
                                 (non-zero) analytic values, i.e. the parent does
                                 NOT resolve the child's building -- the case in
                                 which the solver must mask them (review F2).
  nesting_solenoidal.<expnr>.nc  the analytic field projected onto the discretely
                                 solenoidal subspace of the child grid (w = 0 on
                                 floor and lid), cut into slabs the way
                                 make_child_case cuts a parent dump, constant in
                                 time (U48, review F9).
  assertfire_nan.<expnr>.nc      the corrected file with one NaN in u_west at time
                                 level 2. The cheap init-time flux check passes
                                 (the stored residual is clean), so the abort has
                                 to come from the slab validator when the level is
                                 read (F5 abort case, runmode 1009).
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np

from udprep.nesting import (
    ANALYTIC_COEFFS,
    FaceMasks,
    NestGrid,
    NestingData,
    analytic_initial_fields,
    analytic_slabs,
    apply_divergence_correction,
    net_volume_flux,
    slabs_from_fields,
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
# West-face footprint of the U45 test box (global, 1-based, inclusive): the
# j and k extent of nest_set_solid_box(1, 2, 5, 8, 3, 6) in src/tests.f90.
MASK_J = (5, 8)
MASK_K = (3, 6)


def grid() -> NestGrid:
    return NestGrid.uniform(ITOT, JTOT, KTOT, XLEN, YLEN, ZSIZE)


def _data(times, slab_times=None, corrected=False, initial=False, slabs=None,
          masks=None, parent_model="analytic") -> NestingData:
    """The analytic field as a NestingData, optionally corrected.

    ``slabs`` replaces the analytic slabs (the solenoidal fixture cuts its own);
    ``masks`` are the lateral fluid masks the correction is restricted to
    (recorded in the data, so the writer stores the matching residual and
    fluid area).
    """
    g = grid()
    slab_times = times if slab_times is None else slab_times
    data = NestingData(
        grid=g,
        nzone=NZONE,
        times=np.asarray(times, dtype=np.float64),
        slabs=analytic_slabs(g, NZONE, slab_times) if slabs is None else slabs,
        initial_fields=(analytic_initial_fields(g, float(np.asarray(slab_times)[0]))
                        if initial else None),
        parent_model=parent_model,
        parent_dx=float(XLEN / ITOT),
        parent_dt=float(np.min(np.diff(times))) if len(times) > 1 else 0.0,
        child_origin_x=0.0,
        child_origin_y=0.0,
    )
    if corrected:
        apply_divergence_correction(data, masks=masks)
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


def _poison(src: Path, dst: Path) -> None:
    """Copy ``src`` to ``dst`` and put one NaN into ``u_west`` at time level 2.

    The stored ``flux_residual`` stays clean, so the init-time cheap check
    passes and only the per-level slab validator can catch the value.
    """
    import netCDF4 as nc

    shutil.copy2(src, dst)
    with nc.Dataset(dst, "a") as ds:
        ds.variables["u_west"][1, 0, 0, 0] = np.nan


def write_all(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    analytic = outdir / f"nesting_analytic.{EXPNR}.nc"
    write_nesting_file(analytic, _data(TIMES), override=True)

    corrected = outdir / f"nesting_corrected.{EXPNR}.nc"
    write_nesting_file(corrected, _data(TIMES, corrected=True), override=True)

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

    # --- U45 (review F2): the west faces inside a box the parent does not resolve ---
    # Corrected over the fluid faces only; the solid faces keep their analytic
    # values, so a solver that failed to mask them would inject their flux.
    west = np.ones((JTOT, KTOT), dtype=bool)
    west[MASK_J[0] - 1:MASK_J[1], MASK_K[0] - 1:MASK_K[1]] = False
    write_nesting_file(
        outdir / f"nesting_maskwest.{EXPNR}.nc",
        _data(TIMES, corrected=True, masks=FaceMasks(west=west)),
        override=True,
    )

    # --- U48 (review F9): slabs cut from a discretely solenoidal full field ---
    # apply_divergence_correction on a data set WITH an initial-condition block
    # syncs the block to the corrected faces and projects it (w = 0 at floor
    # and lid); the slabs are then cut from that field the way make_child_case
    # cuts a parent dump, and held constant in time.
    base = _data(TIMES, corrected=True, initial=True)
    f = base.initial_fields
    level = slabs_from_fields(grid(), NZONE, f["u"], f["v"], f["w"])
    slabs = {k: np.repeat(v[None], len(TIMES), axis=0) for k, v in level.items()}
    write_nesting_file(
        outdir / f"nesting_solenoidal.{EXPNR}.nc",
        _data(TIMES, corrected=True, slabs=slabs, parent_model="analytic-projected"),
        override=True,
    )

    # --- F5: a NaN the cheap init check cannot see ---
    _poison(corrected, outdir / f"assertfire_nan.{EXPNR}.nc")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="directory to write the fixtures into")
    args = parser.parse_args()
    write_all(args.outdir)
    print(f"nesting fixtures written to {args.outdir}")
    print(f"  coefficients: {ANALYTIC_COEFFS}")


if __name__ == "__main__":
    main()
