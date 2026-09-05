#!/usr/bin/env python3
"""Presets for the V1 "Big Brother" nesting validation experiment.

One :class:`Preset` fixes *every* number the experiment depends on, so the
tiny smoke-test version and the production version differ only in the values
here -- never in the code path.  ``docs/udales-nesting-design.md`` section 10.4
row V1 is what this implements; ``README.md`` in this directory explains the
protocol and the acceptance criterion.

The parent is a periodic neutral urban LES on an aligned cube array, driven by
a **fixed** mean pressure gradient ``dpdx = ustar**2 / zsize``.  The child is
the central sub-region of that parent, at the same resolution and carrying the
same buildings, driven at its lateral boundaries by slabs cut from the parent's
own field dumps and by the same ``dpdx``.

Why a fixed ``dpdx`` and not ``luvolflowr``: the two runs must share the same
momentum source *exactly*.  A volume-flow-rate controller adjusts a body force
in time to hold a target flux, and the two runs would then be driven
differently -- the child's controller would mask precisely the boundary error
the experiment is trying to measure.  With a fixed pressure gradient the source
is identical by construction and the friction velocity ``ustar`` is known a
priori, which also gives the error metrics a natural, run-independent scale.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np


@dataclass(frozen=True)
class Preset:
    """Every parameter of one V1 experiment.

    Lengths are metres, times seconds.  ``dx = dy = dz`` throughout: the point
    of V1 is a refinement ratio of exactly 1, so parent and child share the grid
    and the interpolation contributes nothing (P1-P11 cover it separately).
    """

    name: str

    # -- parent grid ------------------------------------------------------- #
    itot: int
    jtot: int
    ktot: int
    dx: float

    # -- geometry: an aligned array of cubes ------------------------------- #
    building_height: float
    building_width: float
    street_width: float
    #: target facet edge length handed to the cube-array generator
    edgelength: float
    #: ``"plaza"``   -- the cubes whose footprint would fall in (or within
    #:                 ``nwall`` cells of) the child's guard + ramp band are
    #:                 removed from the **parent**, so the child's zone sits
    #:                 over open ground.  This is the configuration the scheme
    #:                 was designed for (design section 5, "no buildings
    #:                 anywhere W > 0"), and it lets ``nest_lparentgeom`` be
    #:                 ``.false.`` so that ``nesting_init`` *asserts* the rule
    #:                 rather than warning about it.
    #: ``"uniform"`` -- the unbroken array, which puts buildings inside the
    #:                 zone.  Legal only for self-nesting (the parent resolves
    #:                 the same geometry) and needs ``nest_lparentgeom = .true.``
    #:                 It is the harder case and is kept available, but it is
    #:                 not the V1 default.
    geometry: str

    # -- child window, in parent columns ----------------------------------- #
    child_itot: int
    child_jtot: int

    # -- relaxation zone (design section 1.1, 1.4) ------------------------- #
    guardwidth: float   #: L_imp [m]
    zonewidth: float    #: L_rel [m]
    tau: float          #: relaxation time [s]
    nzone: int          #: zone thickness stored in the nesting file, in cells
    nwall: int          #: nest_nwall, wall erosion in cells
    #: nest_timeinterp.  1 = linear, 2 = monotone cubic Hermite.  **Use 1.**
    #: Mode 2's Fritsch-Carlson slope limiter is a nonlinear function of the
    #: four buffered levels, so the interpolated boundary field is not a fixed
    #: linear combination of levels whose net flux is individually zero, and
    #: the divergence compatibility of design section 3.1 is lost between
    #: parent levels.  Measured on the tiny preset: mode 2 gives a normalised
    #: flux residual of 5.2e-5 (nest_fluxtol is 1e-10, so nest_lfluxassert
    #: aborts the run) and, with the assertion off, divmax = 2.1e-4 and
    #: divtot = 1.7; mode 1 gives Phi = 8.9e-16 and divmax = 4.0e-16.
    #: See README.md, "Finding N1".
    timeinterp: int
    #: cold-start the child from the parent's own field at times[0]
    #: (nest_linitfromparent + the schema-2 u_init/v_init/w_init block)

    # -- forcing ----------------------------------------------------------- #
    ustar: float        #: target friction velocity; sets dpdx = ustar^2/zsize
    u0: float           #: initial uniform streamwise velocity [m/s]
    tke0: float         #: initial e12 in prof.inp

    # -- schedule ---------------------------------------------------------- #
    spinup: float       #: parent spin-up before the production window [s]
    production: float   #: parent production window, dumping fields [s]
    dtdump: float       #: field dump interval during the production window [s]
    child_spinup: float #: child time discarded before statistics start [s]

    # -- runtime ----------------------------------------------------------- #
    nprocx: int
    nprocy: int
    child_nprocx: int
    child_nprocy: int
    dtmax: float

    # -- analysis ---------------------------------------------------------- #
    #: heights [m] at which streamwise spectra are taken (below / at / above h)
    spectra_heights: Tuple[float, ...] = (8.0, 16.0, 32.0)
    #: keep every n-th dumped level when accumulating statistics
    stride: int = 1
    init_from_parent: bool = True
    #: experiment numbers
    parent_expnr: str = "903"
    child_expnr: str = "904"

    # -- derived ----------------------------------------------------------- #

    @property
    def dy(self) -> float:
        return self.dx

    @property
    def dz(self) -> float:
        return self.dx

    @property
    def xlen(self) -> float:
        return self.itot * self.dx

    @property
    def ylen(self) -> float:
        return self.jtot * self.dx

    @property
    def zsize(self) -> float:
        return self.ktot * self.dx

    @property
    def child_ktot(self) -> int:
        return self.ktot

    @property
    def child_xlen(self) -> float:
        return self.child_itot * self.dx

    @property
    def child_ylen(self) -> float:
        return self.child_jtot * self.dx

    @property
    def period(self) -> float:
        """Streamwise/spanwise period of the cube array [m]."""
        return self.building_width + self.street_width

    @property
    def child_i0(self) -> int:
        """0-based parent column index of the child's west face."""
        return (self.itot - self.child_itot) // 2

    @property
    def child_j0(self) -> int:
        """0-based parent column index of the child's south face."""
        return (self.jtot - self.child_jtot) // 2

    @property
    def child_origin(self) -> Tuple[float, float]:
        return (self.child_i0 * self.dx, self.child_j0 * self.dx)

    @property
    def dpdx(self) -> float:
        """Body force [m/s^2] that gives the requested friction velocity."""
        return self.ustar ** 2 / self.zsize

    @property
    def zone_cells(self) -> int:
        """Total zone thickness in child cells, L_imp + L_rel."""
        return int(np.ceil((self.guardwidth + self.zonewidth) / self.dx - 1.0e-9))

    @property
    def t_start(self) -> float:
        """Parent time at which the production window (and the child) begins."""
        return self.spinup

    @property
    def t_end(self) -> float:
        return self.spinup + self.production

    @property
    def optical_depth(self) -> float:
        """Design section 1.4(b) absorption depth D at the bulk velocity ``u0``.

        A disturbance leaving the domain is attenuated by exp(-D); a reflection
        makes the round trip, so its residual amplitude is ~exp(-2 D).  D >~ 2.3
        gives 99 % suppression.
        """
        return (self.guardwidth + 0.5 * self.zonewidth) / (self.u0 * self.tau)

    # -- the cube array ------------------------------------------------------ #

    @property
    def building_free_zone(self) -> bool:
        """True when the parent geometry keeps the child's zone clear."""
        return self.geometry == "plaza"

    @property
    def zone_clearance(self) -> float:
        """How far from a lateral face a building must stay, in metres.

        ``L_imp + L_rel`` is where ``W`` becomes zero; ``nest_nwall`` cells more
        because the weights are additionally eroded away from any solid point
        (design section 5), so a building just outside the ramp would still eat
        zone points.
        """
        return self.guardwidth + self.zonewidth + self.nwall * self.dx

    def _full_cube_centres(self) -> np.ndarray:
        """Centres of the unbroken aligned array, in parent metres.

        Same formula and same ordering as ``udgeom.create_cubes(..., 'AC')``:
        ``c = i (C + H) - H/2 - C/2`` for ``i = 1..N``, ``i`` outer, ``j`` inner.
        ``test_v1_tiny`` asserts that the mesh built from this reproduces
        ``create_cubes`` exactly, so the reimplementation cannot drift.
        """
        h, c = self.building_width, self.street_width
        nx = int(round(self.xlen / (h + c)))
        ny = int(round(self.ylen / (h + c)))
        out = [(i * (c + h) - h / 2.0 - c / 2.0, j * (c + h) - h / 2.0 - c / 2.0)
               for i in range(1, nx + 1) for j in range(1, ny + 1)]
        return np.asarray(out, dtype=float)

    def _child_box(self) -> Tuple[float, float, float, float]:
        x0, y0 = self.child_origin
        return x0, y0, x0 + self.child_xlen, y0 + self.child_ylen

    def _interior_box(self) -> Tuple[float, float, float, float]:
        x0, y0, x1, y1 = self._child_box()
        d = self.zone_clearance
        return x0 + d, y0 + d, x1 - d, y1 - d

    def cube_centres(self) -> np.ndarray:
        """Centres of the cubes the **parent** actually carries, in parent metres.

        For ``geometry = "uniform"`` that is the whole array.  For ``"plaza"``
        every cube is required to be either entirely inside the child's interior
        or entirely outside the child, so the guard + ramp band -- and a
        ``nest_nwall`` margin around it -- is open ground.  Anything else is
        removed.
        """
        centres = self._full_cube_centres()
        if self.geometry == "uniform":
            return centres
        if self.geometry != "plaza":
            raise ValueError(f"unknown geometry {self.geometry!r}")
        half = 0.5 * self.building_width
        cx0, cy0, cx1, cy1 = self._child_box()
        ix0, iy0, ix1, iy1 = self._interior_box()
        keep = []
        for cx, cy in centres:
            x0, x1, y0, y1 = cx - half, cx + half, cy - half, cy + half
            overlaps_child = (x1 > cx0 and x0 < cx1 and y1 > cy0 and y0 < cy1)
            inside_interior = (x0 >= ix0 and x1 <= ix1 and y0 >= iy0 and y1 <= iy1)
            if (not overlaps_child) or inside_interior:
                keep.append((cx, cy))
        return np.asarray(keep, dtype=float).reshape(-1, 2)

    def cube_centres_in(self, x0: float, y0: float, xsize: float, ysize: float
                        ) -> np.ndarray:
        """The same cubes restricted to a window, in that window's own metres.

        A cube is kept when its whole footprint lies in the window, so cutting
        the child out of the parent reproduces the parent's array exactly -- the
        smoke test checks that cell by cell against the IBM solid masks.
        """
        half = 0.5 * self.building_width
        out = []
        for cx, cy in self.cube_centres():
            if (cx - half >= x0 - 1.0e-9 and cx + half <= x0 + xsize + 1.0e-9
                    and cy - half >= y0 - 1.0e-9 and cy + half <= y0 + ysize + 1.0e-9):
                out.append((cx - x0, cy - y0))
        return np.asarray(out, dtype=float).reshape(-1, 2)

    @property
    def n_cubes_removed(self) -> int:
        return len(self._full_cube_centres()) - len(self.cube_centres())

    # -- self-consistency --------------------------------------------------- #

    def validate(self) -> None:
        """Fail loudly on a preset that cannot mean what it says."""
        errors = []
        p = self.period
        for label, size in (("xlen", self.xlen), ("ylen", self.ylen),
                            ("child_xlen", self.child_xlen),
                            ("child_ylen", self.child_ylen)):
            if abs(size / p - round(size / p)) > 1.0e-9:
                errors.append(f"{label} = {size} is not a whole number of {p} m cube periods")
        # The child must start on a period boundary, or its cube array would be
        # out of phase with the parent's and the two geometries would differ.
        for label, origin in zip(("child_origin_x", "child_origin_y"), self.child_origin):
            if abs(origin / p - round(origin / p)) > 1.0e-9:
                errors.append(
                    f"{label} = {origin} m is not a multiple of the {p} m cube period; "
                    "the child's buildings would not line up with the parent's"
                )
        if self.nzone < self.zone_cells:
            errors.append(
                f"nzone = {self.nzone} is thinner than the {self.zone_cells}-cell zone "
                f"(L_imp + L_rel = {self.guardwidth + self.zonewidth} m)"
            )
        if 2 * self.zone_cells >= min(self.child_itot, self.child_jtot):
            errors.append("the zone covers the whole child; there is no interior to compare")
        for n, tot, label in ((self.nprocx, self.itot, "nprocx/itot"),
                              (self.nprocy, self.jtot, "nprocy/jtot"),
                              (self.child_nprocx, self.child_itot, "child nprocx/itot"),
                              (self.child_nprocy, self.child_jtot, "child nprocy/jtot")):
            if tot % n:
                errors.append(f"{label}: {tot} is not divisible by {n}")
        if self.child_spinup >= self.production:
            errors.append("child_spinup leaves no statistics window")
        if self.building_height % self.dx:
            errors.append("building_height is not a whole number of cells")
        if self.timeinterp not in (1, 2):
            errors.append("timeinterp must be 1 (linear) or 2 (Hermite)")
        if self.geometry not in ("plaza", "uniform"):
            errors.append(f"geometry must be 'plaza' or 'uniform', got {self.geometry!r}")
        elif self.geometry == "plaza":
            # The point of the plaza is that the assertion in nesting_init
            # (nest_lparentgeom = .false.) can be switched on, so check here that
            # it will pass rather than discovering it three stages later.
            half = 0.5 * self.building_width
            ix0, iy0, ix1, iy1 = self._interior_box()
            inside = 0
            for cx, cy in self.cube_centres():
                if (cx - half >= ix0 - 1e-9 and cx + half <= ix1 + 1e-9
                        and cy - half >= iy0 - 1e-9 and cy + half <= iy1 + 1e-9):
                    inside += 1
                elif (cx + half > self.child_origin[0] and cx - half < self.child_origin[0] + self.child_xlen
                      and cy + half > self.child_origin[1] and cy - half < self.child_origin[1] + self.child_ylen):
                    errors.append(
                        f"a cube at ({cx:g}, {cy:g}) m survives inside the child but not "
                        "inside its interior; the zone would not be building-free"
                    )
            if inside == 0:
                errors.append(
                    "the plaza layout leaves no buildings at all in the child interior; "
                    "the child would be an empty box and V1 would test nothing"
                )
        if errors:
            raise ValueError(
                f"preset '{self.name}' is inconsistent:\n  " + "\n  ".join(errors)
            )

    def summary(self) -> str:
        u = self.ustar
        lines = [
            f"preset               {self.name}",
            f"parent               {self.itot} x {self.jtot} x {self.ktot} cells, "
            f"{self.xlen:g} x {self.ylen:g} x {self.zsize:g} m, dx = {self.dx:g} m",
            f"child                {self.child_itot} x {self.child_jtot} x {self.child_ktot} cells, "
            f"origin ({self.child_origin[0]:g}, {self.child_origin[1]:g}) m",
            f"buildings            {self.building_width:g} m cubes, h = {self.building_height:g} m "
            f"({self.building_height / self.dx:g} cells), {self.street_width:g} m streets, "
            f"period {self.period:g} m",
            f"geometry             '{self.geometry}': {len(self.cube_centres())} cubes "
            f"({self.n_cubes_removed} removed of {len(self._full_cube_centres())}), "
            f"zone {'building-free' if self.building_free_zone else 'CONTAINS buildings'} "
            f"(clearance {self.zone_clearance:g} m), "
            f"nest_lparentgeom = {'.false.' if self.building_free_zone else '.true.'}",
            f"zone                 L_imp = {self.guardwidth:g} m "
            f"({self.guardwidth / self.dx:g} cells), L_rel = {self.zonewidth:g} m "
            f"({self.zonewidth / self.dx:g} cells), tau = {self.tau:g} s, nzone = {self.nzone}",
            f"                     total {self.zone_cells} cells per side, "
            f"optical depth D = {self.optical_depth:.2f} (want >= 2.3)",
            f"interior             {self.child_itot - 2 * self.zone_cells} x "
            f"{self.child_jtot - 2 * self.zone_cells} cells = "
            f"{(self.child_itot - 2 * self.zone_cells) * self.dx:g} x "
            f"{(self.child_jtot - 2 * self.zone_cells) * self.dy:g} m = "
            f"{(self.child_itot - 2 * self.zone_cells) * self.dx / self.building_height:.0f}h x "
            f"{(self.child_jtot - 2 * self.zone_cells) * self.dy / self.building_height:.0f}h",
            f"forcing              dpdx = {self.dpdx:.4e} m/s^2 -> ustar = {u:g} m/s",
            f"time interpolation   nest_timeinterp = {self.timeinterp} "
            f"({'linear' if self.timeinterp == 1 else 'monotone Hermite'})",
            f"child init           {'from the parent block' if self.init_from_parent else 'from prof.inp'}",
            f"schedule             spin-up {self.spinup:g} s, production "
            f"[{self.t_start:g}, {self.t_end:g}] s, dump every {self.dtdump:g} s",
            f"                     child spin-up {self.child_spinup:g} s, statistics over "
            f"{self.production - self.child_spinup:g} s",
            f"ranks                parent {self.nprocx} x {self.nprocy}, "
            f"child {self.child_nprocx} x {self.child_nprocy}",
        ]
        return "\n".join(lines)


#: Production experiment: design section 10.4 V1 as costed in the brief.
PRODUCTION = Preset(
    name="production",
    itot=256, jtot=256, ktot=64, dx=2.0,
    building_height=16.0, building_width=16.0, street_width=16.0, edgelength=16.0,
    geometry="plaza",
    child_itot=128, child_jtot=128,
    # design section 1.4: N_imp = 3 cells, N_rel = 9 cells -> 12 cells total,
    # the recommended default.  L_rel = 18 m >= max(8 dx, h) = 16 m.
    guardwidth=6.0, zonewidth=18.0, tau=1.0, nzone=12, nwall=1, timeinterp=1,
    ustar=0.4, u0=3.0, tke0=0.1,
    spinup=3600.0, production=1800.0, dtdump=3.0, child_spinup=300.0,
    nprocx=8, nprocy=8, child_nprocx=8, child_nprocy=8,
    dtmax=0.5,
    spectra_heights=(8.0, 16.0, 32.0),
    stride=1,
)

#: V1 rerun sized by the statistics rather than by the compute.
#:
#: The first production run (3600 s spin-up, 1800 s window) left the parent's
#: own half-window spread at 13.6 % above z/h = 2, which is the same size as
#: the child's TKE deficit there -- so it could not say whether the deficit was
#: real.  It also finished with the interior TKE still climbing, i.e. the
#: spin-up had not reached equilibrium, and with the child still decorrelating
#: from the parent through the first ~700 s of the averaging window.
#:
#: This preset triples the spin-up, takes the window out to 10800 s (a 6.8x
#: longer statistics window once the longer child spin-up is discarded, so the
#: sampling spread should fall by ~2.6x to about 5 %), and discards the
#: decorrelation transient.  Everything else is identical to PRODUCTION so the
#: two runs stay comparable -- in particular dtdump stays at 3 s, since the
#: parent dump cadence is part of what is under test.
#:
#: Costed from the measured rates of run 3990472 on 64 ranks: parent 2.36 s
#: simulated per s wall while spinning up and 1.96 while dumping, child 6.5x
#: real time.  Parent 76 + 92 min, child 28 min, builds and analysis ~20 min
#: -> about 3.8 h against an 8 h walltime.  Field dumps are ~180 GB.
CONVERGED = Preset(
    name="converged",
    itot=256, jtot=256, ktot=64, dx=2.0,
    building_height=16.0, building_width=16.0, street_width=16.0, edgelength=16.0,
    geometry="plaza",
    child_itot=128, child_jtot=128,
    guardwidth=6.0, zonewidth=18.0, tau=1.0, nzone=12, nwall=1, timeinterp=1,
    ustar=0.4, u0=3.0, tke0=0.1,
    spinup=10800.0, production=10800.0, dtdump=3.0, child_spinup=600.0,
    nprocx=8, nprocy=8, child_nprocx=8, child_nprocy=8,
    dtmax=0.5,
    spectra_heights=(8.0, 16.0, 32.0),
    stride=1,
)

#: Smoke test: the identical pipeline at a size that runs on a login node in
#: minutes.  Everything that differs is a number in this object.
TINY = Preset(
    name="tiny",
    itot=96, jtot=96, ktot=32, dx=2.0,
    building_height=16.0, building_width=16.0, street_width=16.0, edgelength=16.0,
    geometry="plaza",
    child_itot=64, child_jtot=64,
    guardwidth=6.0, zonewidth=8.0, tau=1.0, nzone=7, nwall=1, timeinterp=1,
    ustar=0.4, u0=3.0, tke0=0.1,
    spinup=120.0, production=120.0, dtdump=3.0, child_spinup=40.0,
    nprocx=2, nprocy=2, child_nprocx=2, child_nprocy=2,
    dtmax=0.5,
    spectra_heights=(8.0, 16.0, 32.0),
    stride=1,
)

PRESETS: Dict[str, Preset] = {p.name: p for p in (TINY, PRODUCTION, CONVERGED)}


def get_preset(name: str) -> Preset:
    try:
        preset = PRESETS[name]
    except KeyError:
        raise SystemExit(
            f"unknown preset {name!r}; choose one of {', '.join(sorted(PRESETS))}"
        ) from None
    preset.validate()
    return preset


if __name__ == "__main__":
    for name in sorted(PRESETS):
        print(get_preset(name).summary())
        print()
