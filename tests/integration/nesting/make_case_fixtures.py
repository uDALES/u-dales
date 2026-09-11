#!/usr/bin/env python3
"""Build the run directories for the nesting integration matrix I1-I8.

Companion to :mod:`make_fixtures`, which builds the fixtures for the in-solver
unit runmodes (U1-U34).  This module builds *whole cases* -- a nesting input
file, ``prof.inp``, ``lscale.inp`` and a ``namoptions`` -- so the solver can be
run end to end against a parent field whose exact properties are known
analytically.  See ``docs/udales-nesting-design.md`` section 10.3 for the
matrix and ``README.md`` in this directory for what each case is for.

Every nesting file is written through the *production* writer
``tools/python/udprep/nesting.py`` (``slabs_from_fields`` +
``apply_divergence_correction`` + ``write_nesting_file``), so the writer and
the Fortran reader cannot drift apart, exactly as for the unit fixtures.

The parent fields are defined on the **child** grid (refinement ratio 1), which
is what isolates the solver: any error seen in a run is the solver's, not the
interpolation's.  ``tools/python/tests/test_nesting.py`` (P1-P11) covers the
interpolation separately.

Analytic fields
---------------

``uniform``
    ``u = U``, ``v = w = 0``.  Discretely solenoidal, zero net boundary flux,
    and a fixed point of every term in the momentum equation on a free-slip
    box.  Test I3.

``face_forced``
    ``u = U0 + A cos(2 pi x/Lx) cos(2 pi y/Ly) sin(2 pi z/Lz)``, ``v = w = 0``.
    Chosen so the two ``x`` boundary faces carry *identical* (and non-uniform)
    data -- so ``Phi = 0`` exactly -- while ``du/dx != 0`` in the interior, so
    the projection has real work to do.  Test I2: the face values must survive
    that projection.

``taylor_green``
    ``u = A sin(2 pi x/L) cos(2 pi y/L)``, ``v = -A cos(2 pi x/L) sin(2 pi y/L)``,
    ``w = 0`` -- the field design section 10.3 names for I4.  On a staggered
    (MAC) grid its *discrete* divergence vanishes identically, not merely to
    ``O(h^2)``; see the note in ``README.md``.

``mixed_mode``
    From the stream function ``psi = A sin(2 pi x/L) sin(4 pi y/L)``:
    ``u = dpsi/dy``, ``v = -dpsi/dx``, ``w = 0``.  Analytically solenoidal, but
    because the two directions carry *different* wavenumbers the discrete
    divergence is genuinely ``O(h^2)``.  This is the field the second-order
    convergence half of I4 needs.

``unsteady_uniform``
    ``u = U0 + A sin(2 pi t/T)``, ``v = w = 0``: uniform in space, non-trivial
    in time, so the time buffer and its restart repositioning matter.  Test I6.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from udprep.nesting import (
    NestGrid,
    NestingData,
    apply_divergence_correction,
    discrete_divergence as _nesting_divergence,
    initial_fields_from_fields,
    net_volume_flux,
    slabs_from_fields,
    write_nesting_file,
)

TEST_DIR = Path(__file__).resolve().parent
EXPNR = "902"

#: (u, v, w) on the child grid, at each component's own stagger.
Fields = Tuple[np.ndarray, np.ndarray, np.ndarray]


# --------------------------------------------------------------------------- #
# Case geometry
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CaseSpec:
    """Grid and zone geometry of one integration case."""

    itot: int = 32
    jtot: int = 32
    ktot: int = 16
    xlen: float = 32.0
    ylen: float = 32.0
    zsize: float = 16.0
    #: zone thickness stored in the file, in child cells
    nzone: int = 16
    #: L_imp [m]
    guardwidth: float = 16.0
    #: L_rel [m]
    zonewidth: float = 0.0
    #: relaxation time [s]; <= 0 means Dirichlet wherever W > 0
    tau: float = 0.0

    @property
    def dx(self) -> float:
        return self.xlen / self.itot

    @property
    def dy(self) -> float:
        return self.ylen / self.jtot

    @property
    def dz(self) -> float:
        return self.zsize / self.ktot

    def grid(self) -> NestGrid:
        return NestGrid.uniform(
            self.itot, self.jtot, self.ktot, self.xlen, self.ylen, self.zsize
        )


def full_zone_spec(itot: int, jtot: int, ktot: int,
                   xlen: float = 32.0, ylen: float = 32.0, zsize: float = 16.0) -> CaseSpec:
    """A spec with ``W == 1`` everywhere: a pure guard strip covering the box.

    ``L_imp`` is half the domain, so every point is within ``L_imp`` of a
    lateral face and the bounded union of design section 1.1 gives ``W = 1``.
    ``nzone`` is half the horizontal size, which is the smallest value for
    which the west/east (south/north) slabs together cover every child index.
    """
    if itot % 2 or jtot % 2:
        raise ValueError("full_zone_spec needs an even itot and jtot")
    if abs(xlen / itot - ylen / jtot) > 1.0e-12:
        raise ValueError("full_zone_spec assumes dx == dy")
    return CaseSpec(
        itot=itot, jtot=jtot, ktot=ktot,
        xlen=xlen, ylen=ylen, zsize=zsize,
        nzone=min(itot, jtot) // 2,
        guardwidth=0.5 * min(xlen, ylen),
        zonewidth=0.0,
        tau=0.0,
    )


#: I2/I3/I5 baseline: 32 x 32 x 16 cells over 32 x 32 x 16 m, W == 1 everywhere.
BASE = full_zone_spec(32, 32, 16)

#: I4 convergence ladder: the same box at h = 2, 1 and 0.5 m.
CONVERGENCE = [
    full_zone_spec(16, 16, 8),
    full_zone_spec(32, 32, 16),
    full_zone_spec(64, 64, 32),
]

#: I6/I7 realistic zone: 3 m guard + 8 m ramp, tau = 4 s (design section 1.4).
ZONED = CaseSpec(
    itot=32, jtot=32, ktot=16,
    xlen=32.0, ylen=32.0, zsize=16.0,
    nzone=12, guardwidth=3.0, zonewidth=8.0, tau=4.0,
)


# --------------------------------------------------------------------------- #
# Analytic parent fields, sampled at each component's own stagger
# --------------------------------------------------------------------------- #


def _mesh(g: NestGrid, component: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = g.component_coords(component, 0)
    y = g.component_coords(component, 1)
    z = g.component_coords(component, 2)
    return np.meshgrid(x, y, z, indexing="ij")


def uniform(g: NestGrid, t: float = 0.0, U: float = 1.0) -> Fields:
    u = np.full(g.component_shape("u"), float(U))
    v = np.zeros(g.component_shape("v"))
    w = np.zeros(g.component_shape("w"))
    return u, v, w


def unsteady_uniform(g: NestGrid, t: float, U0: float = 1.0,
                     A: float = 0.25, T: float = 20.0) -> Fields:
    return uniform(g, U=U0 + A * math.sin(2.0 * math.pi * t / T))


def face_forced(g: NestGrid, t: float = 0.0, U0: float = 1.0, A: float = 0.2) -> Fields:
    """Non-uniform on the x faces, identical on both of them, divergent inside."""
    x, y, z = _mesh(g, "u")
    u = U0 + A * (np.cos(2.0 * np.pi * x / g.xlen)
                  * np.cos(2.0 * np.pi * y / g.ylen)
                  * np.sin(2.0 * np.pi * z / g.zsize))
    v = np.zeros(g.component_shape("v"))
    w = np.zeros(g.component_shape("w"))
    return u, v, w


def taylor_green(g: NestGrid, t: float = 0.0, A: float = 1.0) -> Fields:
    """The field of design section 10.3, I4.  ``L`` is the (square) box side."""
    L = g.xlen
    xu, yu, _ = _mesh(g, "u")
    xv, yv, _ = _mesh(g, "v")
    u = A * np.sin(2.0 * np.pi * xu / L) * np.cos(2.0 * np.pi * yu / L)
    v = -A * np.cos(2.0 * np.pi * xv / L) * np.sin(2.0 * np.pi * yv / L)
    w = np.zeros(g.component_shape("w"))
    return u, v, w


def mixed_mode(g: NestGrid, t: float = 0.0, A: float = 1.0) -> Fields:
    """``u = dpsi/dy``, ``v = -dpsi/dx`` for ``psi = A sin(2pi x/L) sin(4pi y/L)``.

    Solenoidal analytically; its *discrete* divergence is ``O(h^2)`` because the
    two directions carry different wavenumbers, so it does not benefit from the
    exact cancellation :func:`taylor_green` enjoys on a MAC grid.
    """
    L = g.xlen
    kx = 2.0 * np.pi / L
    ky = 4.0 * np.pi / L
    xu, yu, _ = _mesh(g, "u")
    xv, yv, _ = _mesh(g, "v")
    u = A * ky * np.sin(kx * xu) * np.cos(ky * yu)
    v = -A * kx * np.cos(kx * xv) * np.sin(ky * yv)
    w = np.zeros(g.component_shape("w"))
    return u, v, w


FIELDS: Dict[str, Callable[..., Fields]] = {
    "uniform": uniform,
    "unsteady_uniform": unsteady_uniform,
    "face_forced": face_forced,
    "taylor_green": taylor_green,
    "mixed_mode": mixed_mode,
}


def discrete_divergence(g: NestGrid, fields: Fields) -> np.ndarray:
    """Cell-centred discrete divergence, the same operator ``fillps`` applies.

    Used by the fixture self-check below and by the tests, so that a claim
    about a manufactured field is verified rather than asserted.  It is the
    production operator of ``udprep.nesting``, so the fixture check and the
    writer's own projection cannot disagree about what "divergence" means.
    """
    return _nesting_divergence(g, *fields)


# --------------------------------------------------------------------------- #
# Writing a case
# --------------------------------------------------------------------------- #


def build_nesting_data(spec: CaseSpec, fieldname: str, times: Sequence[float],
                       correct: bool = True, initial: bool = False,
                       **kwargs) -> NestingData:
    """Cut the twelve zone slabs out of an analytic child field, per time level.

    With ``initial=True`` the field at the first time is also carried whole, as
    the schema 2 initial-condition block, so a cold start can be run with
    ``nest_linitfromparent``.  The correction then projects it (design section
    10.6 item 4); with ``correct=False`` it is stored exactly as built.
    """
    g = spec.grid()
    fn = FIELDS[fieldname]
    times = np.asarray(times, dtype=np.float64).reshape(-1)
    per_time = [slabs_from_fields(g, spec.nzone, *fn(g, float(t), **kwargs)) for t in times]
    slabs = {name: np.stack([s[name] for s in per_time], axis=0) for name in per_time[0]}
    data = NestingData(
        grid=g,
        nzone=spec.nzone,
        times=times,
        slabs=slabs,
        initial_fields=(initial_fields_from_fields(g, *fn(g, float(times[0]), **kwargs))
                        if initial else None),
        parent_model=f"analytic:{fieldname}",
        parent_dx=spec.dx,
        parent_dt=float(np.min(np.diff(times))) if times.size > 1 else 0.0,
        child_origin_x=0.0,
        child_origin_y=0.0,
    )
    if correct:
        apply_divergence_correction(data)
    else:
        data.net_volume_flux = net_volume_flux(data)
    return data


def write_profiles(rundir: Path, spec: CaseSpec, uprof: float = 0.0,
                   vprof: float = 0.0, thl: float = 288.0,
                   e12: float = 0.0, expnr: str = EXPNR) -> None:
    """``prof.inp`` (which also *defines* the vertical grid) and ``lscale.inp``.

    Column order is the one ``modstartup::readinitfiles`` actually reads:
    ``z thl qt u v e12`` -- not the order the header comment of some existing
    cases suggests.
    """
    zf = spec.grid().zf
    with (rundir / f"prof.inp.{expnr}").open("w", encoding="utf-8") as fh:
        fh.write("# nesting integration case: uniform vertical grid\n")
        fh.write("# z thl qt u v e12\n")
        for z in zf:
            fh.write(f"{z:18.11f} {thl:16.9f} {0.0:16.9f} "
                     f"{uprof:16.9f} {vprof:16.9f} {e12:16.9f}\n")
    with (rundir / f"lscale.inp.{expnr}").open("w", encoding="utf-8") as fh:
        fh.write("# nesting integration case: no large-scale forcing\n")
        fh.write("# z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls dthlrad\n")
        for z in zf:
            fh.write(f"{z:18.11f}" + "".join(f" {0.0:16.9f}" for _ in range(9)) + "\n")


def render_namoptions(spec: CaseSpec, edits: Optional[Dict[str, str]] = None,
                      template: Optional[Path] = None) -> str:
    """Fill the shared ``namoptions.902`` template for one case."""
    template = template or (TEST_DIR / f"namoptions.{EXPNR}")
    text = template.read_text(encoding="utf-8")
    settings: Dict[str, str] = {
        "itot": str(spec.itot),
        "jtot": str(spec.jtot),
        "ktot": str(spec.ktot),
        "xlen": f"{spec.xlen:.10g}",
        "ylen": f"{spec.ylen:.10g}",
        "zsize": f"{spec.zsize:.10g}",
        "nest_guardwidth": f"{spec.guardwidth:.10g}",
        "nest_zonewidth": f"{spec.zonewidth:.10g}",
        "nest_tau": f"{spec.tau:.10g}",
    }
    settings.update(edits or {})
    for key, value in settings.items():
        text, n = re.subn(
            rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$", lambda m: m.group(1) + value, text
        )
        if n == 0:
            raise RuntimeError(f"setting '{key}' not found in {template}")
    return text


def write_case(rundir: Path, spec: CaseSpec, fieldname: str,
               times: Sequence[float] = (0.0, 1.0e6),
               edits: Optional[Dict[str, str]] = None,
               uprof: float = 0.0, correct: bool = True, initial: bool = False,
               expnr: str = EXPNR, **kwargs) -> NestingData:
    """Write a complete, runnable case directory and return the parent data."""
    rundir.mkdir(parents=True, exist_ok=True)
    data = build_nesting_data(spec, fieldname, times, correct=correct,
                              initial=initial, **kwargs)
    write_nesting_file(rundir / f"nesting.inp.{expnr}.nc", data, override=True)
    write_profiles(rundir, spec, uprof=uprof, expnr=expnr)
    (rundir / f"namoptions.{expnr}").write_text(
        render_namoptions(spec, edits), encoding="utf-8"
    )
    return data


# --------------------------------------------------------------------------- #
# Self-check: what the manufactured fields actually do, discretely
# --------------------------------------------------------------------------- #


def report(spec: CaseSpec) -> List[str]:
    """One line per field: its peak discrete divergence and net boundary flux."""
    g = spec.grid()
    lines = [f"grid {g.itot}x{g.jtot}x{g.ktot}  h = {spec.dx:g} m"]
    for name in ("uniform", "face_forced", "taylor_green", "mixed_mode"):
        f = FIELDS[name](g, 0.0)
        div = discrete_divergence(g, f)
        data = build_nesting_data(spec, name, [0.0], correct=False)
        phi = float(net_volume_flux(data)[0])
        lines.append(f"  {name:16s} max|div| = {np.abs(div).max():10.3e}   Phi = {phi:10.3e}")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, nargs="?",
                        help="write the I3 baseline case here (a smoke test)")
    args = parser.parse_args()
    for spec in CONVERGENCE:
        for line in report(spec):
            print(line)
    if args.outdir is not None:
        write_case(args.outdir, BASE, "uniform", uprof=1.0)
        print(f"I3 baseline case written to {args.outdir}")


if __name__ == "__main__":
    main()
