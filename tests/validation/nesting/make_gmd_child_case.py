#!/usr/bin/env python3
"""Build the GMD 2024 indoor-outdoor nested child case (plan:
``~/udales/gmd-indoor-plan-2026-09-10.md``, "Child box" and "Decisions").

Parent: experiment 567, the paper-resolution run (1024 x 512 x 128 cells on
3.42 x 1.80 x 0.96 m, stretched vertical), driven by the synthetic-inflow
precursor and dumping its ``&NESTPARENT`` band around the child box
i = 108..323, j = 171..341 (216 x 171 parent cells, full depth).  This module
builds the child case that box drives: the enclosure of
``indoor_geometry.py`` (GMD 17, 6277, 2024, section 4.2) at the paper's wall
(20 mm) or thinner, at the parent's own horizontal resolution (``r = 1``,
slabs cut -- the V1 no-op standard) or twice as fine (``r = 2``, slabs
interpolated).  The vertical is never refined: the child keeps the parent's
own stretched ``zh`` exactly, which is also why ``&INPS``'s stretching
parameters are copied verbatim rather than re-derived (see
:func:`child_sections`) -- ``UDPrep`` reproduces the same vertical grid from
them for the child's own preprocessing.

Two-stage build, deliberately split so the first stage is testable without
the parent's simulation having produced anything:

* :func:`build_case` (steps 1-3: grid, geometry, namelist) reads only the
  parent's ALREADY-WRITTEN ``namoptions.567``/``prof.inp.567`` -- static
  inputs that exist as soon as the parent case is built, not run output.
* :func:`build_from_parent_run` (steps 4-5: preprocessing, the nesting file)
  additionally reads the parent's ``nesting.out.*.567.nc`` band files (and,
  once, ``nesting.out.init.*.567.nc``), which only exist after the parent has
  actually run.

Usage
-----
    python make_gmd_child_case.py <parent_dir> <outdir> --expnr 601 \\
        [--refine 1|2] [--wall 0.020] [--case-only]
"""

from __future__ import annotations

import argparse
import json
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np

import indoor_geometry as ig
from caselib import (
    NestParent,
    check_finite_slabs,
    load_solid_mask,
    run_preprocessing,
    write_lscale,
    write_namoptions,
    write_prof,
)
from udprep.nesting import (
    DEFAULT_PROLONGATION,
    PROLONGATIONS,
    NestGrid,
    NestingWriter,
    check_alignment,
    face_masks_from_ibm,
    initial_fields_from_parent,
    slabs_from_fields,
    slabs_from_parent,
    stagger_masks_from_ibm,
)

# --------------------------------------------------------------------------- #
# The parent's own geometry (paper grid, experiment 567) -- fixed facts of the
# published case, not read off any file, since the horizontal grid is uniform
# and its size is part of the case's identity.  Only the STRETCHED vertical is
# read back (from prof.inp.567, see parent_grid()), because it cannot be
# reproduced analytically here without re-running UDPrep's own stretching.
# --------------------------------------------------------------------------- #

PARENT_ITOT = 1024
PARENT_JTOT = 512
PARENT_KTOT = 128
PARENT_XLEN = 3.42
PARENT_YLEN = 1.80
PARENT_EXPNR = "567"

#: Child box, in parent cell indices -- plan "Child box": i = 108..323 (216
#: cells), j = 171..341 (171 cells), full depth.  Parametrised with these as
#: defaults so a different box can be tried without editing this module.
CHILD_I0 = 108
CHILD_J0 = 171
CHILD_ITOT_R1 = 216
CHILD_JTOT_R1 = 171

#: Guard + ramp, in CHILD cells, independent of the refinement ratio (plan:
#: "Zone: guard 4 + ramp 8 child cells").
GUARD_CELLS = 4
ZONE_CELLS = 8
NZONE_CHILD = GUARD_CELLS + ZONE_CELLS

#: Relaxation time in TIMESTEPS.  The design (docs/udales-nesting-design.md,
#: "tau = n_tau dt") makes tau a multiple of the timestep, and the campaign's
#: PRODUCTION preset ran tau = 1.0 s at dtmax = 0.5 s, i.e. n_tau = 2, for
#: V1-V6.  The number that transfers between cases is n_tau, not the second:
#: this case's timestep is 0.4-0.8 ms and its zone is 40 mm wide, crossed by
#: the flow in ~8 ms, so a literal 1.0 s would leave the zone essentially
#: unforced.
N_TAU_STEPS = 2


def nest_tau(dtmax_child: float) -> float:
    """``nest_tau`` [s] for a child with timestep ``dtmax_child``: n_tau = 2."""
    return N_TAU_STEPS * float(dtmax_child)


def parent_dx(itot: int = PARENT_ITOT, xlen: float = PARENT_XLEN) -> float:
    """The parent's own (uniform) x spacing [m]."""
    return xlen / itot


def parent_dy(jtot: int = PARENT_JTOT, ylen: float = PARENT_YLEN) -> float:
    """The parent's own (uniform) y spacing [m]."""
    return ylen / jtot


# --------------------------------------------------------------------------- #
# Namelist parsing -- just enough of &SECTION ... key = value ... / to clone
# and patch namoptions.567.  Not a general Fortran namelist parser: scalar
# assignments only, which is everything namoptions.567 uses.
# --------------------------------------------------------------------------- #

_BOOL_TRUE = {".true.", ".t.", "t", "true"}
_BOOL_FALSE = {".false.", ".f.", "false"}


def _parse_nml_scalar(raw: str) -> Any:
    s = raw.strip()
    low = s.lower()
    if low in _BOOL_TRUE:
        return True
    if low in _BOOL_FALSE:
        return False
    if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
        return s[1:-1]
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        # An unquoted string (namoptions.567 writes 'stl_file = indoor_object_final.stl'
        # with no quotes) -- caselib.write_namoptions's own _fmt quotes any str it is
        # given, so round-tripping through here still produces a valid namelist.
        return s


def read_namoptions_sections(path: Path) -> "OrderedDict[str, OrderedDict[str, Any]]":
    """Parse ``namoptions.<expnr>`` into ``{SECTION: {key: typed value}}``, in file order.

    Good enough for a namelist that :func:`child_sections` is about to rewrite
    with :func:`caselib.write_namoptions` (whose ``_fmt`` needs typed Python
    values, not raw text) -- not a general namelist reader.
    """
    sections: "OrderedDict[str, OrderedDict[str, Any]]" = OrderedDict()
    current: Optional[str] = None
    for line in Path(path).read_text(encoding="ascii", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("!"):
            continue
        if stripped.startswith("&"):
            current = stripped[1:].strip().upper()
            sections[current] = OrderedDict()
            continue
        if stripped == "/":
            current = None
            continue
        if current is None or "=" not in stripped:
            continue
        key, _, value = stripped.partition("=")
        sections[current][key.strip()] = _parse_nml_scalar(value.split("!", 1)[0].strip())
    return sections


def read_parent_zf(parent_dir: Path, expnr: str = PARENT_EXPNR) -> np.ndarray:
    """The parent's own vertical cell centres -- column 1 of ``prof.inp.<expnr>``.

    ``prof.inp`` *defines* the vertical grid (``caselib.write_prof``); reading
    it back is how this module learns the parent's stretched ``zf`` without
    re-running UDPrep's stretching construction.
    """
    path = Path(parent_dir) / f"prof.inp.{expnr}"
    data = np.loadtxt(path, comments="#", ndmin=2)
    return np.asarray(data[:, 0], dtype=np.float64)


def zh_from_zf(zf: np.ndarray) -> np.ndarray:
    """Reconstruct cell faces from centres: ``zh[0] = 0``, ``zh[k+1] = 2 zf[k] - zh[k]``.

    The exact inverse of how a cell centre sits at the midpoint of its two
    faces; applied here because ``prof.inp`` stores only ``zf``.
    """
    zf = np.asarray(zf, dtype=np.float64)
    zh = np.empty(zf.size + 1, dtype=np.float64)
    zh[0] = 0.0
    for k in range(zf.size):
        zh[k + 1] = 2.0 * zf[k] - zh[k]
    return zh


def parent_grid(parent_dir: Path, expnr: str = PARENT_EXPNR) -> NestGrid:
    """The paper-resolution parent's own grid, full domain, absolute coordinates."""
    dx, dy = parent_dx(), parent_dy()
    xh = np.arange(PARENT_ITOT + 1, dtype=np.float64) * dx
    yh = np.arange(PARENT_JTOT + 1, dtype=np.float64) * dy
    zf = read_parent_zf(parent_dir, expnr)
    if zf.size != PARENT_KTOT:
        raise ValueError(
            f"prof.inp.{expnr} has {zf.size} vertical levels, expected {PARENT_KTOT}"
        )
    zh = zh_from_zf(zf)
    return NestGrid.from_faces(xh, yh, zh)


# --------------------------------------------------------------------------- #
# The child box and grid
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ChildBox:
    """The child's footprint, in the parent's own (absolute) coordinates."""

    refine: int
    i0: int = CHILD_I0
    j0: int = CHILD_J0
    itot_r1: int = CHILD_ITOT_R1
    jtot_r1: int = CHILD_JTOT_R1

    def __post_init__(self) -> None:
        if self.refine not in (1, 2):
            raise ValueError(f"refine must be 1 or 2 (horizontal only), got {self.refine}")

    @property
    def itot(self) -> int:
        return self.itot_r1 * self.refine

    @property
    def jtot(self) -> int:
        return self.jtot_r1 * self.refine

    @property
    def x0(self) -> float:
        return self.i0 * parent_dx()

    @property
    def y0(self) -> float:
        return self.j0 * parent_dy()

    @property
    def xlen(self) -> float:
        return self.itot_r1 * parent_dx()

    @property
    def ylen(self) -> float:
        return self.jtot_r1 * parent_dy()


def child_grid(box: ChildBox, parent: NestGrid) -> NestGrid:
    """The child's own grid: refined horizontally, the parent's own ``zh`` vertically.

    Raises :class:`udprep.nesting.NestingAlignmentError` if the box does not
    nest in the parent's cells (:func:`udprep.nesting.check_alignment`).
    """
    r = box.refine
    xh = box.x0 + np.arange(box.itot + 1, dtype=np.float64) * (parent_dx() / r)
    yh = box.y0 + np.arange(box.jtot + 1, dtype=np.float64) * (parent_dy() / r)
    grid = NestGrid.from_faces(xh, yh, parent.zh.copy())
    check_alignment(parent, grid)
    return grid


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #


def build_child_stl(casedir: Path, nr: str, box: ChildBox, t: float) -> "tuple[Path, float]":
    """Write the enclosure STL in the CHILD's own coordinates.

    Refuses (``indoor_geometry.WallTooThin``) a wall the child's horizontal
    grid cannot carry (``check_wall(t, dx/r, dy/r)``).  The vertical rule is
    the parent's own dz (not refined, and not enforced here as a hard limit):
    returns the roof's own cell count for the caller to record and warn about
    if it falls below the paper's 3-cell rule, rather than raising.
    """
    dx_child = parent_dx() / box.refine
    dy_child = parent_dy() / box.refine
    t = ig.check_wall(t, dx_child, dy_child)
    stl_path = casedir / f"geom.{nr}.stl"
    ig.write_stl(stl_path, t, ground=True, domain=(box.xlen, box.ylen),
                 origin=(box.x0, box.y0))
    return stl_path, t


# --------------------------------------------------------------------------- #
# Namelist
# --------------------------------------------------------------------------- #


def child_sections(
    parent_sections: "OrderedDict[str, OrderedDict[str, Any]]",
    nr: str, box: ChildBox, t: float, dtmax_child: float,
    runtime: float, trestart: float, nprocx: int, nprocy: int, stl_name: str,
) -> "OrderedDict[str, OrderedDict[str, Any]]":
    """The child's namoptions: the parent's own sections, changed only where a child needs.

    Everything not named below -- ``&PHYSICS`` (``ps``, ``igrw_damp``),
    ``&DYNAMICS`` (``ipoiss = 0``, already what ``lnesting`` requires),
    ``&NAMSUBGRID``, ``&OUTPUT``, and ``&BC``'s ``BCtopm``/roughness/surface
    keys (``z0``, ``z0h``, ``wtsurf``, ``wqsurf``, ``thls``) -- is the
    parent's own value, untouched.
    """
    sec: "OrderedDict[str, OrderedDict[str, Any]]" = OrderedDict(
        (name, OrderedDict(fields)) for name, fields in parent_sections.items()
    )

    run = sec["RUN"]
    run["iexpnr"] = int(nr) if str(nr).isdigit() else nr
    run["dtmax"] = float(dtmax_child)
    run["runtime"] = float(runtime)
    run["trestart"] = float(trestart)
    run["nprocx"] = int(nprocx)
    run["nprocy"] = int(nprocy)

    dom = sec["DOMAIN"]
    dom["itot"] = int(box.itot)
    dom["jtot"] = int(box.jtot)
    dom["ktot"] = PARENT_KTOT
    dom["xlen"] = float(box.xlen)
    dom["ylen"] = float(box.ylen)

    bc = sec["BC"]
    bc["BCxm"] = 4   # BCxm_nesting
    bc["BCym"] = 3   # BCym_nesting
    # BCtopm and the roughness/surface keys are whatever the parent had --
    # left alone by the clone above.

    # &DRIVER: this child is driven by nesting, not by inflow planes.
    sec["DRIVER"] = OrderedDict([("idriver", 0)])

    # &WALLS: keep iwallmom (and any of its siblings the parent set); drop the
    # solid/boundary/facet-section counts, which are the PARENT's -- wrong for
    # this case's own subdomain and STL. Preprocessing (run_preprocessing)
    # recomputes and writes the child's own counts, as caselib.run_preprocessing
    # / make_child_case.build's WALLS template already assumes.
    walls_src = sec.get("WALLS", OrderedDict())
    sec["WALLS"] = OrderedDict(
        (k, v) for k, v in walls_src.items()
        if not k.startswith(("nfcts", "nsolpts", "nbndpts", "nfctsecs"))
    )

    # &INPS: only the geometry changes; zsize/lzstretch/hlin/dzlin/lstretchexp/
    # stretchconst/u0/diag_neighbs stay the parent's own, so UDPrep reproduces
    # the same vertical grid for the child.
    inps = sec["INPS"]
    inps["stl_file"] = stl_name
    inps["stl_ground"] = True

    # &NESTPARENT: this is a child, not a parent.
    sec.pop("NESTPARENT", None)

    dx_child = parent_dx() / box.refine
    sec["NESTING"] = OrderedDict([
        ("lnesting", True),
        ("nest_guardwidth", GUARD_CELLS * dx_child),
        ("nest_zonewidth", ZONE_CELLS * dx_child),
        ("nest_tau", nest_tau(dtmax_child)),
        ("nest_shape", 1),
        ("nest_timeinterp", 2),
        ("nest_nwall", 1),
        ("nest_lparentgeom", bool(abs(float(t) - ig.PAPER_WALL) < 1.0e-9)),
        ("nest_linitfromparent", True),
        ("nest_lendabort", True),
    ])
    return sec


# --------------------------------------------------------------------------- #
# Steps 1-3: grid, geometry, namelist -- testable without the parent's run
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ChildCase:
    """Everything :func:`build_from_parent_run` needs, produced by :func:`build_case`."""

    casedir: Path
    box: ChildBox
    t: float
    nr: str
    parent_grid: NestGrid
    child_grid: NestGrid
    dtmax: float
    nprocx: int
    nprocy: int
    prolongation: str
    roof_cells: float


def build_case(
    parent_dir: Path, outdir: Path, expnr: Any, refine: int, wall: float = ig.PAPER_WALL,
    *, runtime: float = 150.0, trestart: float = 12.5,
    nprocx: int = 8, nprocy: int = 9,
    i0: int = CHILD_I0, j0: int = CHILD_J0,
    itot_r1: int = CHILD_ITOT_R1, jtot_r1: int = CHILD_JTOT_R1,
    prolongation: str = DEFAULT_PROLONGATION,
) -> ChildCase:
    """Steps 1-3: grid, geometry and namelist.

    Reads only the parent's already-written ``namoptions.567``/``prof.inp.567``
    (static inputs, not run output), so this can be exercised before the
    parent has actually run.
    """
    if prolongation not in PROLONGATIONS:
        raise ValueError(f"unknown prolongation {prolongation!r}; expected one of {PROLONGATIONS}")
    nr = str(expnr)
    parent_dir = Path(parent_dir)
    parent = parent_grid(parent_dir)
    box = ChildBox(refine=refine, i0=i0, j0=j0, itot_r1=itot_r1, jtot_r1=jtot_r1)
    child = child_grid(box, parent)  # raises NestingAlignmentError if not nested

    casedir = Path(outdir) / nr
    casedir.mkdir(parents=True, exist_ok=True)

    parent_sections = read_namoptions_sections(parent_dir / f"namoptions.{PARENT_EXPNR}")

    _, t = build_child_stl(casedir, nr, box, wall)
    dzlin = float(parent_sections["INPS"]["dzlin"])
    roof_cells = t / dzlin
    if roof_cells < 3.0:
        print(f"[make_gmd_child_case] WARNING: the roof spans {roof_cells:.1f} cells "
              f"(t = {t * 1.0e3:g} mm over the parent's dz = {dzlin * 1.0e3:g} mm); "
              "the paper's own 3-cell rule is not met vertically -- accepted, not refused, "
              "since the vertical is not refined")

    dtmax_child = float(parent_sections["RUN"]["dtmax"]) / refine
    sections = child_sections(parent_sections, nr, box, t, dtmax_child,
                              runtime, trestart, nprocx, nprocy, f"geom.{nr}.stl")
    write_namoptions(
        casedir / f"namoptions.{nr}", sections,
        header=[
            "GMD 2024 indoor-outdoor nested child, experiment 567 parent",
            "generated by tests/validation/nesting/make_gmd_child_case.py; do not hand-edit",
            f"parent {parent_dir}, refine r = {refine}, wall {t * 1.0e3:g} mm",
        ],
    )

    zf = read_parent_zf(parent_dir)
    prof_data = np.loadtxt(parent_dir / f"prof.inp.{PARENT_EXPNR}", comments="#", ndmin=2)
    write_prof(casedir / f"prof.inp.{nr}", zf,
               u=prof_data[:, 3], v=prof_data[:, 4], e12=prof_data[:, 5],
               thl=float(prof_data[0, 1]), qt=float(prof_data[0, 2]),
               comment="cut from the parent's prof.inp.567 (same z levels; vertical ratio 1)")
    write_lscale(casedir / f"lscale.inp.{nr}", zf,
                 comment="cut from the parent's lscale.inp.567 (all zero; forcing is the "
                         "inflow/nesting, not dpdx)")

    return ChildCase(casedir=casedir, box=box, t=t, nr=nr, parent_grid=parent,
                     child_grid=child, dtmax=dtmax_child, nprocx=nprocx, nprocy=nprocy,
                     prolongation=prolongation, roof_cells=roof_cells)


# --------------------------------------------------------------------------- #
# Steps 4-5: preprocessing and the nesting file -- need the parent's own run
# --------------------------------------------------------------------------- #


def build_from_parent_run(parent_dir: Path, case: ChildCase, *,
                          ibm_backend: str = "auto") -> Path:
    """Steps 4-5: preprocessing and the nesting file.

    Needs the parent's ``nesting.out.*.567.nc`` band files (and, once,
    ``nesting.out.init.*.567.nc``) to already exist -- i.e. the parent must
    have actually run. Call after :func:`build_case`.
    """
    parent_dir = Path(parent_dir)
    casedir = case.casedir
    box = case.box
    r = box.refine
    nr = case.nr
    parent = case.parent_grid
    child = case.child_grid
    prolongation = case.prolongation
    dx_child = parent_dx() / r

    run_preprocessing(casedir, ibm_backend=ibm_backend)
    child_fluid = load_solid_mask(
        casedir, (child.itot, child.jtot, child.ktot))
    child_masks = face_masks_from_ibm(child_fluid)

    pgrid = NestGrid.from_faces(
        box.x0 + np.arange(box.itot_r1 + 1, dtype=np.float64) * parent_dx(),
        box.y0 + np.arange(box.jtot_r1 + 1, dtype=np.float64) * parent_dy(),
        parent.zh,
    )
    fluid_window = load_solid_mask(
        parent_dir, (PARENT_ITOT, PARENT_JTOT, PARENT_KTOT)
    )[box.i0:box.i0 + box.itot_r1, box.j0:box.j0 + box.jtot_r1, :PARENT_KTOT]
    parent_masks = stagger_masks_from_ibm(fluid_window)

    dump = NestParent(parent_dir, PARENT_EXPNR, parent_dx())
    try:
        times = np.asarray(dump.times, dtype=np.float64)
        if times.size == 0:
            raise RuntimeError(f"{parent_dir}: the parent wrote no nestparent band levels")
        t_offset = float(times[0])
        parent_dt = (float(np.median(np.diff(times))) if times.size > 1
                    else float(dump.tnestparent))

        span_chunk = min(box.itot_r1 * r // case.nprocx, box.jtot_r1 * r // case.nprocy)
        nestfile = casedir / f"nesting.inp.{nr}.nc"
        writer = NestingWriter(
            nestfile, child, NZONE_CHILD,
            masks=child_masks,
            parent_model=f"udales:{PARENT_EXPNR}:gmd-indoor-outdoor r={r}",
            parent_dx=parent_dx(),
            parent_dy=parent_dy(),
            parent_dz=float(np.min(np.diff(parent.zh))),
            parent_dt=parent_dt,
            child_origin_x=box.x0,
            child_origin_y=box.y0,
            child_dt=case.dtmax,
            span_chunk=span_chunk,
        )
        phi_after = np.empty(times.size, dtype=np.float64)
        has_initial = False
        with writer:
            for n in range(times.size):
                pu, pv, pw = dump.read_level(n)
                cu, cv, cw = dump.child_block(pu, pv, pw, box.i0, box.j0,
                                              box.itot_r1, box.jtot_r1)
                initial_fields = None
                if n == 0:
                    if not dump.has_init:
                        raise RuntimeError(
                            f"{parent_dir}: no nesting.out.init.*.{PARENT_EXPNR}.nc; the "
                            "parent must run with nestparent_linit = .true. (already set "
                            "in namoptions.567)"
                        )
                    iu, iv, iw = dump.read_init()
                    if abs(float(dump.init_time) - t_offset) > 1.0e-6:
                        raise RuntimeError(
                            f"nestparent_init is stamped t = {dump.init_time} but the "
                            f"first band level is t = {times[0]}"
                        )
                    if r == 1:
                        initial_fields = {"u": iu.copy(), "v": iv.copy(), "w": iw.copy()}
                    else:
                        initial_fields = initial_fields_from_parent(
                            pgrid, iu, iv, iw, child, parent_masks=parent_masks,
                            prolongation=prolongation,
                        )
                    has_initial = True
                if r == 1:
                    level = slabs_from_fields(child, NZONE_CHILD, cu, cv, cw)
                else:
                    level = slabs_from_parent(pgrid, cu, cv, cw, child=child,
                                              nzone=NZONE_CHILD, parent_masks=parent_masks,
                                              prolongation=prolongation)
                check_finite_slabs(level, dump.nzone, parent_dx())
                result = writer.append_level(times[n] - t_offset, level,
                                             initial_fields=initial_fields)
                phi_after[n] = result["flux_residual"]
        diagnostics = writer.diagnostics
    finally:
        dump.close()

    manifest = {
        "case": "gmd2024-indoor-outdoor",
        "parent_dir": str(parent_dir.resolve()),
        "parent_expnr": PARENT_EXPNR,
        "child_expnr": nr,
        "refine": r,
        "wall_thickness_m": case.t,
        "nest_lparentgeom": bool(abs(case.t - ig.PAPER_WALL) < 1.0e-9),
        "roof_cells": case.roof_cells,
        "box": {"i0": box.i0, "j0": box.j0, "x0": box.x0, "y0": box.y0,
                "xlen": box.xlen, "ylen": box.ylen},
        "nzone_child_cells": NZONE_CHILD,
        "nest_guardwidth_m": GUARD_CELLS * dx_child,
        "nest_zonewidth_m": ZONE_CELLS * dx_child,
        "nest_tau_s": nest_tau(case.dtmax),
        "n_tau_steps": N_TAU_STEPS,
        "prolongation": prolongation if r > 1 else None,
        "n_parent_levels": int(times.size),
        "t_offset": t_offset,
        "runtime": float(times[-1] - t_offset),
        "parent_dt_median": parent_dt,
        "child_dtmax": case.dtmax,
        "nprocx": case.nprocx,
        "nprocy": case.nprocy,
        "writer_diagnostics": diagnostics,
        "has_initial_condition": has_initial,
        "flux_residual_after_correction_max_abs": float(np.max(np.abs(phi_after))),
    }
    (casedir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n",
                                           encoding="ascii")
    return casedir


def build(parent_dir: Path, outdir: Path, expnr: Any, refine: int, wall: float = ig.PAPER_WALL,
          *, ibm_backend: str = "auto", **kwargs: Any) -> Path:
    """Convenience wrapper: :func:`build_case` then :func:`build_from_parent_run`."""
    case = build_case(parent_dir, outdir, expnr, refine, wall, **kwargs)
    return build_from_parent_run(parent_dir, case, ibm_backend=ibm_backend)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("parent_dir", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--expnr", required=True)
    parser.add_argument("--refine", type=int, choices=(1, 2), default=1)
    parser.add_argument("--wall", type=float, default=ig.PAPER_WALL,
                        help="wall thickness [m] (default: the paper's 0.020 m)")
    parser.add_argument("--runtime", type=float, default=150.0)
    parser.add_argument("--trestart", type=float, default=12.5)
    parser.add_argument("--nprocx", type=int, default=8)
    parser.add_argument("--nprocy", type=int, default=9)
    parser.add_argument("--ibm-backend", default="auto")
    parser.add_argument("--prolongation", default=DEFAULT_PROLONGATION, choices=PROLONGATIONS)
    parser.add_argument("--case-only", action="store_true",
                        help="stop after the grid/geometry/namelist (steps 1-3); skip "
                             "preprocessing and the nesting file, which need the parent's "
                             "own run to have produced its nesting.out.* band files")
    args = parser.parse_args()

    case = build_case(args.parent_dir, args.outdir, args.expnr, args.refine, args.wall,
                      runtime=args.runtime, trestart=args.trestart,
                      nprocx=args.nprocx, nprocy=args.nprocy, prolongation=args.prolongation)
    print(f"[make_gmd_child_case] child case files written to {case.casedir} "
          f"(r = {args.refine}, wall = {case.t * 1.0e3:g} mm, roof = {case.roof_cells:.1f} cells)")
    if args.case_only:
        print("[make_gmd_child_case] --case-only: stopping before preprocessing/nesting "
              "file (needs the parent's own run)")
        return

    build_from_parent_run(args.parent_dir, case, ibm_backend=args.ibm_backend)
    manifest = json.loads((case.casedir / "manifest.json").read_text())
    print(f"[make_gmd_child_case] nesting file written; {manifest['n_parent_levels']} parent "
          f"levels, child runtime {manifest['runtime']:.3f} s")
    print(f"  Phi after correction  max |Phi| = "
          f"{manifest['flux_residual_after_correction_max_abs']:.3e} m3/s")


if __name__ == "__main__":
    main()
