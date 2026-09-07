#!/usr/bin/env python3
"""Shared machinery for the V1 Big Brother nesting validation.

Nothing here is specific to one preset: everything takes a
:class:`config.Preset` and does the same thing at either size, which is what
makes the tiny smoke test an actual test of the production pipeline.

Contents
--------
``write_namoptions``      render a namelist file from an ordered dict of dicts
``cube_geometry``         the aligned cube array, via ``udgeom.create_cubes``
``write_prof``/``_lscale``  the 1-D input profiles (prof.inp defines the z grid)
``run_preprocessing``     the repo's standard ``UDPrep`` path for the IBM inputs
``FieldDump``             reader for the per-rank ``fielddump.XXX.YYY.<nr>.nc``
``NestDump``              reader for the per-rank ``nestdump[_init].XXX.YYY.<nr>.nc``
``load_solid_mask``       the IBM solid cell-centre mask, from ``solid_c.txt``
``cell_centred``          face-staggered (u, v, w) -> co-located cell centres
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
PYTOOLS = REPO_ROOT / "tools" / "python"
if str(PYTOOLS) not in sys.path:
    sys.path.insert(0, str(PYTOOLS))


# --------------------------------------------------------------------------- #
# Namelists
# --------------------------------------------------------------------------- #


def _fmt(value) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.10g}"
    if isinstance(value, str):
        return value if value.startswith("'") else f"'{value}'"
    raise TypeError(f"cannot render {value!r} into a namelist")


def write_namoptions(path: Path, sections: "OrderedDict[str, Dict[str, object]]",
                     header: Sequence[str] = ()) -> Path:
    """Write ``namoptions.<expnr>`` from ``{section: {key: value}}``."""
    lines: List[str] = [f"! {line}" for line in header]
    for name, entries in sections.items():
        lines.append("")
        lines.append(f"&{name}")
        width = max((len(k) for k in entries), default=1)
        for key, value in entries.items():
            lines.append(f"{key:<{width}} = {_fmt(value)}")
        lines.append("/")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def read_namoption(path: Path, key: str) -> Optional[str]:
    """Return the raw text of ``key`` in a namelist file, or None."""
    pattern = re.compile(rf"(?im)^\s*{re.escape(key)}\s*=\s*(.*?)\s*$")
    match = pattern.search(path.read_text(encoding="ascii", errors="ignore"))
    return match.group(1) if match else None


# --------------------------------------------------------------------------- #
# Geometry and 1-D inputs
# --------------------------------------------------------------------------- #


def cube_mesh(preset, centres: np.ndarray, xsize: float, ysize: float):
    """Build a cube array with an **arbitrary** layout, plus the ground.

    ``udgeom.create_cubes`` only generates the unbroken regular array, and V1's
    default geometry is that array with the cubes over the child's relaxation
    zone removed (a plaza), so the two private helpers ``create_cubes`` itself
    uses are called directly with an explicit list of centres.  This is
    deliberately the *same* pair of calls in the *same* order that
    ``create_cubes(..., 'AC')`` makes, so the tessellation is identical; the
    coupling is turned into a checked invariant by
    ``test_v1_tiny.test_uniform_layout_reproduces_create_cubes``, which requires
    the mesh built here from the full layout to match ``create_cubes`` vertex
    for vertex.
    """
    from udgeom.geometry_generation import (
        _generate_ground_matlab_style,
        _operate_unit_cube,
    )
    from udgeom.udgeom import UDGeom

    centres = np.asarray(centres, dtype=float).reshape(-1, 2)
    if centres.size == 0:
        raise ValueError("a cube layout with no cubes at all is not a useful case")
    divisions = int(round(preset.building_width / preset.edgelength))
    shifts = np.column_stack(
        [centres[:, 0], centres[:, 1],
         np.full(len(centres), 0.5 * preset.building_height)]
    )
    cubes = _operate_unit_cube(
        [preset.building_width, preset.building_width, preset.building_height],
        shifts, divisions,
    )
    ground, _ = _generate_ground_matlab_style(cubes, xsize, ysize, preset.edgelength)
    return UDGeom(stl=ground)


def cube_geometry(preset, x0: float, y0: float, xsize: float, ysize: float,
                  out_stl: Path, centres: Optional[np.ndarray] = None) -> Path:
    """Write the STL for the window ``[x0, x0+xsize] x [y0, y0+ysize]``.

    Parent and child go through this same call with different windows.  The
    default layout is :meth:`config.Preset.cube_centres_in`, the parent's array
    restricted to the window, so there is no second independent description of
    the geometry that could drift.

    ``centres`` overrides it, in **window** coordinates.  The child uses that,
    passing :meth:`config.Preset.child_cube_centres`, because the child's layout
    is allowed to differ from the parent's restriction: with
    ``clear_child_zone`` it drops the cubes that would fall in its guard + ramp
    band, so the zone sits over open ground while the parent keeps them and goes
    on imprinting them through the imposed velocity field (design section 9.4).

    Nothing is committed: the STL is regenerated on every run.
    """
    if centres is None:
        centres = preset.cube_centres_in(x0, y0, xsize, ysize)
    geom = cube_mesh(preset, centres, xsize, ysize)
    out_stl.parent.mkdir(parents=True, exist_ok=True)
    geom.save(str(out_stl))
    return out_stl


def write_prof(path: Path, zf: np.ndarray, u: np.ndarray, v: np.ndarray,
               e12: np.ndarray, thl: float = 288.0, qt: float = 0.0,
               comment: str = "") -> Path:
    """``prof.inp.<expnr>``.  This file also *defines* the vertical grid.

    Column order is the one ``modstartup::readinitfiles`` reads:
    ``z thl qt u v e12``.
    """
    zf = np.asarray(zf, dtype=float)
    u = np.broadcast_to(np.asarray(u, dtype=float), zf.shape)
    v = np.broadcast_to(np.asarray(v, dtype=float), zf.shape)
    e12 = np.broadcast_to(np.asarray(e12, dtype=float), zf.shape)
    with path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write(f"# {comment}\n")
        fh.write("# z thl qt u v e12\n")
        for k, z in enumerate(zf):
            fh.write(f"{z:20.12f} {thl:16.9f} {qt:16.9f} "
                     f"{u[k]:16.9f} {v[k]:16.9f} {e12[k]:16.9f}\n")
    return path


def write_lscale(path: Path, zf: np.ndarray, comment: str = "") -> Path:
    """``lscale.inp.<expnr>``, all zeros.

    **DO NOT "simplify" this to ``UDPrep.forcing.run_all()``.**  The mean
    pressure gradient is supplied through ``dpdx`` in ``&PHYSICS`` and the
    ``pgx`` column here is deliberately zero, because ``modstartup`` forms

        dpdxl(k) = om23_gs*vg(k) - pgx(k) - dpdx        (src/modstartup.f90:2231)

    from the ``lscale.inp`` column **and** the namelist value.
    ``UDPrep.forcing.generate_lscale`` writes ``ls[:, 3] = dpdx`` into that
    column when no forcing switch is set
    (``tools/python/udprep/udprep_forcing.py``), so a case that also sets
    ``dpdx`` in ``namoptions`` and runs the standard forcing preprocessing is
    forced **twice**.  That is a real bug in the preprocessing, tracked
    separately from the nesting work (see README.md, "Finding N3"); this harness
    writes its own ``prof.inp``/``lscale.inp`` to stay out of its way, and uses
    ``UDPrep`` only for the IBM/geometry section.  Whoever fixes the doubling
    should delete this note, not this function.
    """
    zf = np.asarray(zf, dtype=float)
    with path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write(f"# {comment}\n")
        fh.write("# z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls dthlrad\n")
        for z in zf:
            fh.write(f"{z:20.12f}" + "".join(f" {0.0:16.9f}" for _ in range(9)) + "\n")
    return path


# --------------------------------------------------------------------------- #
# Preprocessing
# --------------------------------------------------------------------------- #


def run_preprocessing(casedir: Path, ibm_backend: str = "auto") -> None:
    """Run the repo's standard preprocessing chain on ``casedir``.

    Only the IBM section is exercised: ``prof.inp``/``lscale.inp`` are written
    by :func:`write_prof`/:func:`write_lscale` above rather than by
    ``UDPrep.forcing``, for the ``pgx``/``dpdx`` reason documented there.

    ``ibm_backend='auto'`` prefers the f2py extension and falls back to the
    prebuilt ``IBM_preproc`` executable, so the harness runs on a checkout where
    ``tools/build_preprocessing.sh`` has produced only one of the two.
    """
    from udprep import UDPrep

    prep = UDPrep(str(casedir))
    if ibm_backend != "auto":
        prep.ibm.run_all(backend=ibm_backend)
        return
    try:
        prep.ibm.run_all(backend="f2py")
    except RuntimeError as exc:
        if "ibm_preproc_f2py" not in str(exc):
            raise
        print(f"[caselib] {exc}\n[caselib] falling back to the IBM_preproc executable")
        prep.ibm.run_all(backend="legacy")


# --------------------------------------------------------------------------- #
# Reading uDALES field dumps
# --------------------------------------------------------------------------- #


class FieldDump:
    """Assemble the per-rank ``fielddump.XXX.YYY.<expnr>.nc`` into global fields.

    ``modfielddump`` writes ``u0(ib:ie, jb:je, kb:ke)``, so the returned arrays
    carry ``u`` at ``xh(1..itot)``, ``v`` at ``yh(1..jtot)`` and ``w`` at
    ``zh(1..ktot)``: the upper face of each direction (``u(itot+1)``,
    ``v(jtot+1)``, ``w(ktot+1)``) is **not** in the file.  For a periodic parent
    the first two are the periodic images of index 1 and the third is zero at
    the rigid lid; :meth:`child_block` uses that, and :func:`cell_centred`
    simply drops the last cell in each direction so that parent and child are
    reduced identically.
    """

    def __init__(self, rundir: Path, expnr: str, dx: float):
        self.rundir = Path(rundir)
        self.expnr = str(expnr)
        self.dx = float(dx)
        self.files = sorted(self.rundir.glob(f"fielddump.???.???.{self.expnr}.nc"))
        if not self.files:
            raise FileNotFoundError(
                f"no fielddump.???.???.{self.expnr}.nc under {self.rundir}"
            )
        self._offsets: Dict[Path, Tuple[int, int]] = {}
        self._shapes: Dict[Path, Tuple[int, int, int]] = {}
        self._open: Dict[Path, Any] = {}
        self.times: np.ndarray = np.zeros(0)
        self._probe()

    # The production parent writes ~600 levels across 64 rank files and the
    # pipeline walks the record three times (slab cut, parent statistics, child
    # statistics).  Re-opening every file for every level would be ~10^5 opens on
    # a shared Lustre filesystem, which AGENTS.md warns about specifically, so
    # the handles are kept.
    def _dataset(self, path: Path):
        from netCDF4 import Dataset

        ds = self._open.get(path)
        if ds is None:
            ds = Dataset(path, "r")
            self._open[path] = ds
        return ds

    def close(self) -> None:
        for ds in self._open.values():
            try:
                ds.close()
            except Exception:
                pass
        self._open.clear()

    def __enter__(self) -> "FieldDump":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self) -> None:  # pragma: no cover - best effort
        self.close()

    def _probe(self) -> None:
        from netCDF4 import Dataset

        itot = jtot = ktot = 0
        times = None
        for path in self.files:
            with Dataset(path, "r") as ds:
                xm = np.asarray(ds.variables["xm"][:], dtype=float)
                ym = np.asarray(ds.variables["ym"][:], dtype=float)
                zm = np.asarray(ds.variables["zm"][:], dtype=float)
                t = np.asarray(ds.variables["time"][:], dtype=float)
            i0 = int(round(xm[0] / self.dx))
            j0 = int(round(ym[0] / self.dx))
            if abs(zm[0]) > 1.0e-9:
                raise ValueError(f"{path}: zm does not start at 0")
            self._offsets[path] = (i0, j0)
            self._shapes[path] = (xm.size, ym.size, zm.size)
            itot = max(itot, i0 + xm.size)
            jtot = max(jtot, j0 + ym.size)
            ktot = max(ktot, zm.size)
            if times is None or t.size < times.size:
                times = t
        self.itot, self.jtot, self.ktot = itot, jtot, ktot
        self.times = times if times is not None else np.zeros(0)

    @property
    def ntime(self) -> int:
        return int(self.times.size)

    def read_level(self, n: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Global ``(u, v, w)`` at time level ``n``, shape ``(itot, jtot, ktot)``."""
        u = np.empty((self.itot, self.jtot, self.ktot), dtype=np.float64)
        v = np.empty_like(u)
        w = np.empty_like(u)
        for path in self.files:
            i0, j0 = self._offsets[path]
            ni, nj, nk = self._shapes[path]
            ds = self._dataset(path)
            for name, out in (("u", u), ("v", v), ("w", w)):
                # stored (time, z, y, x); transpose to (x, y, z)
                block = np.asarray(ds.variables[name][n, :, :, :], dtype=np.float64)
                out[i0:i0 + ni, j0:j0 + nj, 0:nk] = block.transpose(2, 1, 0)
        return u, v, w

    def child_block(self, u: np.ndarray, v: np.ndarray, w: np.ndarray,
                    i0: int, j0: int, ni: int, nj: int
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Cut a child-sized staggered block out of global parent fields.

        Returns ``u`` of shape ``(ni + 1, nj, ktot)``, ``v`` of ``(ni, nj + 1,
        ktot)`` and ``w`` of ``(ni, nj, ktot + 1)`` -- the full staggered set the
        nesting writer wants.  ``w`` at the lid is zero (``BCtopm_freeslip``);
        the extra ``u``/``v`` faces come from inside the parent, so no periodic
        wrap is involved as long as the window is not flush with the parent's
        east/north edge.
        """
        if i0 + ni + 1 > self.itot or j0 + nj + 1 > self.jtot:
            raise ValueError(
                "the child window touches the parent's east/north face; the extra "
                "staggered face would need the periodic image, which this cut does "
                "not do"
            )
        uc = u[i0:i0 + ni + 1, j0:j0 + nj, :].copy()
        vc = v[i0:i0 + ni, j0:j0 + nj + 1, :].copy()
        wc = np.zeros((ni, nj, self.ktot + 1), dtype=np.float64)
        wc[:, :, :self.ktot] = w[i0:i0 + ni, j0:j0 + nj, :]
        return uc, vc, wc


def cell_centred(u: np.ndarray, v: np.ndarray, w: np.ndarray
                 ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Co-locate face-staggered velocities at cell centres.

    ``u``, ``v``, ``w`` come in as the raw :meth:`FieldDump.read_level` blocks,
    i.e. shape ``(N, M, K)`` with the *upper* face of each direction missing.
    The result therefore has shape ``(N-1, M-1, K-1)``: the last cell in each
    direction is dropped rather than closed with a boundary value, so that the
    parent sub-region and the child are reduced by exactly the same rule.  The
    cells lost sit on the domain faces, deep inside the guard strip, and are
    excluded from every comparison anyway.
    """
    uc = 0.5 * (u[:-1, :-1, :-1] + u[1:, :-1, :-1])
    vc = 0.5 * (v[:-1, :-1, :-1] + v[:-1, 1:, :-1])
    wc = 0.5 * (w[:-1, :-1, :-1] + w[:-1, :-1, 1:])
    return uc, vc, wc


def load_solid_mask(casedir: Path, shape: Tuple[int, int, int]) -> np.ndarray:
    """Boolean *fluid* mask at cell centres, from ``solid_c.txt``.

    ``True`` where the cell is fluid.  ``solid_c.txt`` holds 1-based ``i j k``
    triples of solid cell centres.
    """
    path = Path(casedir) / "solid_c.txt"
    fluid = np.ones(shape, dtype=bool)
    if not path.exists():
        return fluid
    idx = np.loadtxt(path, comments="#", dtype=int, ndmin=2)
    if idx.size == 0:
        return fluid
    i, j, k = idx[:, 0] - 1, idx[:, 1] - 1, idx[:, 2] - 1
    keep = (i < shape[0]) & (j < shape[1]) & (k < shape[2])
    fluid[i[keep], j[keep], k[keep]] = False
    return fluid


# --------------------------------------------------------------------------- #
# Coarsening a field dump (V0, refinement)
# --------------------------------------------------------------------------- #


def coarsen_staggered(u: np.ndarray, v: np.ndarray, w: np.ndarray, factor: int
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flux-conservative coarsening of one raw :meth:`FieldDump.read_level` set.

    Each coarse face is *co-planar with* a fine face -- every ``factor``-th one
    -- and takes the area-weighted mean of the ``factor x factor`` fine faces it
    contains.  With uniform spacing that is the plain mean, so the coarse face
    flux is exactly the sum of the fine face fluxes it replaces.  Two things
    follow, and both are what make this the right filter for V0 rather than a
    convenience:

    * the coarse field is discretely solenoidal wherever the fine one is (a
      coarse cell's net flux is the sum of the net fluxes of the ``factor**3``
      fine cells inside it), so the coarse parent hands the prolongation of
      design section 1.3 exactly the input it assumes;
    * it is the *exact left inverse* of that prolongation's tangential half --
      piecewise-constant distribution followed by block averaging is the
      identity -- so a round trip through coarsening and interpolation changes
      only what the normal-direction linear interpolation changes.  Any
      difference V0 measures between a filtered-parent child and V1 is therefore
      attributable to the parent's filter scale, not to a mismatch between the
      two operators.

    Inputs are in the raw dump convention (the upper face of each direction is
    missing), and so is the output.  ``factor`` must divide every dimension.
    """
    factor = int(factor)
    if factor < 1:
        raise ValueError(f"coarsening factor must be >= 1, got {factor}")
    if factor == 1:
        return u, v, w
    n = np.asarray(u).shape
    for size in n:
        if size % factor:
            raise ValueError(
                f"field dump shape {n} is not divisible by the coarsening factor {factor}"
            )
    ni, nj, nk = (s // factor for s in n)
    f = factor
    # u sits on xh: keep every f-th x-face, average over the y and z blocks.
    uc = np.asarray(u)[::f].reshape(ni, nj, f, nk, f).mean(axis=(2, 4))
    vc = np.asarray(v)[:, ::f].reshape(ni, f, nj, nk, f).mean(axis=(1, 4))
    wc = np.asarray(w)[:, :, ::f].reshape(ni, f, nj, f, nk).mean(axis=(1, 3))
    return uc, vc, wc


def coarsen_fluid_mask(fluid: np.ndarray, factor: int) -> np.ndarray:
    """Block-AND a cell-centred fluid mask: a coarse cell is fluid iff all of its is.

    Exact for the cube arrays used here, whose faces are aligned to every grid
    in the suite; conservative (it never calls a partly solid coarse cell fluid)
    for anything else.
    """
    factor = int(factor)
    if factor == 1:
        return fluid
    ni, nj, nk = (s // factor for s in fluid.shape)
    f = factor
    return fluid[:ni * f, :nj * f, :nk * f].reshape(ni, f, nj, f, nk, f).all(axis=(1, 3, 5))


class CoarsenedFieldDump:
    """A :class:`FieldDump` presented on a grid ``factor`` times coarser.

    Same interface as ``FieldDump`` for everything the child builder uses --
    ``times``, ``itot``/``jtot``/``ktot``, ``read_level``, ``child_block``,
    ``close`` -- so a coarse view and a genuinely coarse run are interchangeable
    at the call site.  That is what lets V0's two arms differ only in which
    object is constructed.
    """

    def __init__(self, dump: FieldDump, factor: int):
        self.dump = dump
        self.factor = int(factor)
        if self.factor < 1:
            raise ValueError(f"coarsening factor must be >= 1, got {factor}")
        for label, size in (("itot", dump.itot), ("jtot", dump.jtot),
                            ("ktot", dump.ktot)):
            if size % self.factor:
                raise ValueError(
                    f"{label} = {size} is not divisible by the coarsening factor "
                    f"{self.factor}"
                )
        self.itot = dump.itot // self.factor
        self.jtot = dump.jtot // self.factor
        self.ktot = dump.ktot // self.factor
        self.dx = dump.dx * self.factor
        self.times = dump.times
        self.rundir = dump.rundir
        self.expnr = dump.expnr

    @property
    def ntime(self) -> int:
        return int(self.times.size)

    def read_level(self, n: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return coarsen_staggered(*self.dump.read_level(n), factor=self.factor)

    def child_block(self, u, v, w, i0: int, j0: int, ni: int, nj: int):
        return FieldDump.child_block(self, u, v, w, i0, j0, ni, nj)

    def close(self) -> None:
        self.dump.close()

    def __enter__(self) -> "CoarsenedFieldDump":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


# --------------------------------------------------------------------------- #
# Reading the parent-side zone dump (&NESTDUMP, src/modnestdump.f90)
# --------------------------------------------------------------------------- #


class NestDump:
    """Assemble the per-rank ``nestdump.XXX.YYY.<expnr>.nc`` into the child box.

    The parent wrote, per rank, the intersection of its subdomain with each of
    the four strips of the band (``nestdump_nzone`` parent cells inside each
    lateral face of the child box) and, once, its part of the whole box
    (``nestdump_init.XXX.YYY.<expnr>.nc``).  Every block carries its global
    index range as attributes, so the box is assembled here without knowing the
    decomposition; overlapping blocks (the strips meet at the corners, ranks
    share their upper face through the halo) carry identical values.

    Same interface as :class:`FieldDump` for what the child builder uses --
    ``times``, ``itot``/``jtot``/``ktot`` (of the **parent**), ``dx``,
    ``read_level``, ``child_block``, ``close`` -- with one difference:
    :meth:`read_level` returns the **box**, already staggered the way
    ``child_block`` would cut it, and ``NaN`` wherever the parent did not write
    (the interior).  A slab cut or prolongation that reaches outside the band
    therefore comes out non-finite, which ``make_child_case`` checks for, rather
    than silently taking zeros.  :meth:`read_init` returns the whole box.

    ``u`` is at ``xh(i0..i0+ni)`` (``ni + 1`` faces), ``v`` at ``yh(j0..j0+nj)``,
    ``w`` at ``zh(1..ktot+1)`` with the lid value included -- the complete
    staggered set the nesting writer wants.
    """

    def __init__(self, rundir: Path, expnr: str, dx: float):
        self.rundir = Path(rundir)
        self.expnr = str(expnr)
        self.dx = float(dx)
        self.files = sorted(self.rundir.glob(f"nestdump.???.???.{self.expnr}.nc"))
        self.init_files = sorted(self.rundir.glob(f"nestdump_init.???.???.{self.expnr}.nc"))
        if not self.files:
            raise FileNotFoundError(
                f"no nestdump.???.???.{self.expnr}.nc under {self.rundir}"
            )
        self._open: Dict[Path, Any] = {}
        #: per band file: list of (suffix, i1, i2, j1, j2), 0-based inclusive cells
        self._pieces: Dict[Path, List[Tuple[str, int, int, int, int]]] = {}
        self.times: np.ndarray = np.zeros(0)
        self._probe()

    def _dataset(self, path: Path):
        from netCDF4 import Dataset

        ds = self._open.get(path)
        if ds is None:
            ds = Dataset(path, "r")
            self._open[path] = ds
        return ds

    def close(self) -> None:
        for ds in self._open.values():
            try:
                ds.close()
            except Exception:
                pass
        self._open.clear()

    def __enter__(self) -> "NestDump":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self) -> None:  # pragma: no cover - best effort
        self.close()

    @staticmethod
    def _block_range(var) -> Tuple[int, int, int, int]:
        """0-based inclusive cell range of one block, from its attributes."""
        return (int(var.i_start) - 1, int(var.i_end) - 1,
                int(var.j_start) - 1, int(var.j_end) - 1)

    def _probe(self) -> None:
        from netCDF4 import Dataset

        times = None
        header = None
        for path in self.files:
            with Dataset(path, "r") as ds:
                attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
                if header is None:
                    header = attrs
                    if abs(float(attrs["dx"]) - self.dx) > 1.0e-9 * self.dx:
                        raise ValueError(
                            f"{path}: parent dx = {attrs['dx']} but {self.dx} was expected")
                else:
                    for key in ("itot", "jtot", "ktot", "box_i_start", "box_i_end",
                                "box_j_start", "box_j_end", "nzone"):
                        if attrs[key] != header[key]:
                            raise ValueError(f"{path}: {key} differs between rank files")
                pieces = []
                for face in ("west", "east", "south", "north"):
                    name = f"u_{face}"
                    if name in ds.variables:
                        pieces.append((f"_{face}",) + self._block_range(ds.variables[name]))
                self._pieces[path] = pieces
                t = np.asarray(ds.variables["time"][:], dtype=float)
            if times is None or t.size < times.size:
                times = t
        assert header is not None
        self.itot = int(header["itot"])
        self.jtot = int(header["jtot"])
        self.ktot = int(header["ktot"])
        self.i0 = int(header["box_i_start"]) - 1
        self.j0 = int(header["box_j_start"]) - 1
        self.ni = int(header["box_i_end"]) - self.i0
        self.nj = int(header["box_j_end"]) - self.j0
        self.nzone = int(header["nzone"])
        self.tnestdump = float(header["tnestdump"])
        self.times = times if times is not None else np.zeros(0)

    @property
    def ntime(self) -> int:
        return int(self.times.size)

    @property
    def has_init(self) -> bool:
        return bool(self.init_files)

    def _empty_box(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        ni, nj, nk = self.ni, self.nj, self.ktot
        u = np.full((ni + 1, nj, nk), np.nan, dtype=np.float64)
        v = np.full((ni, nj + 1, nk), np.nan, dtype=np.float64)
        w = np.full((ni, nj, nk + 1), np.nan, dtype=np.float64)
        return u, v, w

    def _place(self, out, ds, name: str, i1: int, i2: int, j1: int, j2: int,
               n: Optional[int]) -> None:
        """Copy block ``name`` of ``ds`` into the box arrays ``out = (u, v, w)``."""
        u, v, w = out
        ia, ja = i1 - self.i0, j1 - self.j0
        ni, nj = i2 - i1 + 1, j2 - j1 + 1
        for comp, arr, di, dj in (("u", u, 1, 0), ("v", v, 0, 1), ("w", w, 0, 0)):
            var = ds.variables[f"{comp}{name}"]
            # stored (time, z, y, x) or (z, y, x); to (x, y, z)
            block = np.asarray(var[n, :, :, :] if n is not None else var[:, :, :],
                               dtype=np.float64).transpose(2, 1, 0)
            if ia < 0 or ja < 0 or ia + ni + di > arr.shape[0] or ja + nj + dj > arr.shape[1]:
                raise ValueError(f"block {comp}{name} at cells {i1}..{i2} x {j1}..{j2} "
                                 f"falls outside the box")
            arr[ia:ia + ni + di, ja:ja + nj + dj, :] = block

    def read_level(self, n: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Box ``(u, v, w)`` at time level ``n``; ``NaN`` outside the band."""
        out = self._empty_box()
        for path in self.files:
            ds = self._dataset(path)
            for suffix, i1, i2, j1, j2 in self._pieces[path]:
                self._place(out, ds, suffix, i1, i2, j1, j2, n)
        return out

    def read_init(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The whole box from the ``nestdump_init`` files, and its time stamp."""
        from netCDF4 import Dataset

        if not self.init_files:
            raise FileNotFoundError(
                f"no nestdump_init.???.???.{self.expnr}.nc under {self.rundir}; "
                "the parent ran with nestdump_linit = .false.")
        out = self._empty_box()
        stamp = None
        for path in self.init_files:
            with Dataset(path, "r") as ds:
                i1, i2, j1, j2 = self._block_range(ds.variables["u"])
                self._place(out, ds, "", i1, i2, j1, j2, None)
                t = float(ds.variables["time"][...])
                if stamp is None:
                    stamp = t
                elif t != stamp:
                    raise ValueError(f"{path}: init block time {t} != {stamp}")
        if not all(np.isfinite(a).all() for a in out):
            raise ValueError("the nestdump_init files do not cover the whole box")
        self.init_time = stamp
        return out

    def child_block(self, u: np.ndarray, v: np.ndarray, w: np.ndarray,
                    i0: int, j0: int, ni: int, nj: int
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The box *is* the child block; only check that it is the one asked for."""
        if (i0, j0, ni, nj) != (self.i0, self.j0, self.ni, self.nj):
            raise ValueError(
                f"the child window (cells {i0}..{i0 + ni - 1} x {j0}..{j0 + nj - 1}) is "
                f"not the box the parent dumped ({self.i0}..{self.i0 + self.ni - 1} x "
                f"{self.j0}..{self.j0 + self.nj - 1}); rebuild the parent with the "
                "matching &NESTDUMP box"
            )
        return u, v, w


def check_finite_slabs(slabs: Dict[str, np.ndarray], band_cells: int, dx: float) -> None:
    """Refuse a slab set with non-finite entries -- the band was too thin.

    With :class:`NestDump` the parent's interior is ``NaN``, so a zone slab (or
    a prolongation stencil) that reached past the band shows up here, naming
    the slab, rather than as zeros in the nesting file.
    """
    bad = {name: int((~np.isfinite(arr)).sum()) for name, arr in slabs.items()
           if not np.isfinite(arr).all()}
    if bad:
        raise ValueError(
            f"non-finite values in slabs {sorted(bad)} ({bad}): the parent's nestdump "
            f"band ({band_cells} cells of {dx:g} m) is too thin for this child's zone "
            "and prolongation stencil; increase nestdump_nzone in the parent"
        )


# --------------------------------------------------------------------------- #
# Running the solver
# --------------------------------------------------------------------------- #

DEFAULT_RUNTIME_MODULES = (
    "tools/prod intel/2021a netCDF/4.8.0-iimpi-2021a "
    "netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a"
)


def solver_binary() -> Path:
    build = os.environ.get("UDALES_BUILD")
    if build:
        return Path(build)
    return REPO_ROOT / "build" / "release" / "u-dales"


def run_solver(casedir: Path, namelist: str, nprocs: int, logfile: Path,
               timeout: Optional[float] = None) -> None:
    """Run uDALES in ``casedir`` under ``mpiexec``, streaming stdout to ``logfile``.

    ``namelist`` is the file name passed as ``argv[1]`` -- the solver takes the
    namelist path from the command line, so a directory can hold several
    (``namoptions.<nr>`` and ``namoptions_spinup.<nr>``) and the experiment
    number comes from ``iexpnr`` inside the file.

    The module stack has to be loaded in the *same* shell as the launch (module
    state does not survive between commands on CX3), so the run goes through
    ``bash -lc`` with an explicit ``module load``.  Override the stack with
    ``UDALES_RUNTIME_MODULES`` and the launcher with ``MPIEXEC``.
    """
    binary = solver_binary()
    if not binary.exists():
        raise FileNotFoundError(f"solver binary not found at {binary}")
    modules = os.environ.get("UDALES_RUNTIME_MODULES", DEFAULT_RUNTIME_MODULES)
    mpiexec = os.environ.get("MPIEXEC", "mpiexec")
    script = (
        "set -e\n"
        + (f"module purge && module load {modules}\n" if modules else "")
        + "export HDF5_USE_FILE_LOCKING=FALSE\n"
        + f"cd {casedir}\n"
        + f"{mpiexec} -n {nprocs} {binary} {namelist}\n"
    )
    with logfile.open("w", encoding="utf-8", errors="replace") as fh:
        completed = subprocess.run(
            ["bash", "-lc", script], stdout=fh, stderr=subprocess.STDOUT,
            timeout=timeout,
        )
    if completed.returncode != 0:
        tail = "\n".join(logfile.read_text(errors="replace").splitlines()[-40:])
        raise RuntimeError(
            f"uDALES failed in {casedir} (exit {completed.returncode}); last lines:\n{tail}"
        )
