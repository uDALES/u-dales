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

from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class PlazaWindow:
    """The child window whose relaxation zone carved the parent's plaza.

    The parent geometry is generated once, **for one child**, by removing every
    cube that would fall in that child's guard + ramp band (plus a
    ``nest_nwall`` margin).  A later experiment that reuses the same parent
    field dumps -- which is the whole point of V2 -- must therefore describe the
    parent with the window that carved it, not with its own child window, or
    the cube layout it regenerates is not the one the parent actually ran and
    the child's buildings would not be the parent's.

    Lengths are parent metres.
    """

    x0: float
    y0: float
    xsize: float
    ysize: float
    clearance: float

    @property
    def box(self) -> Tuple[float, float, float, float]:
        return (self.x0, self.y0, self.x0 + self.xsize, self.y0 + self.ysize)

    @property
    def interior_box(self) -> Tuple[float, float, float, float]:
        d = self.clearance
        x0, y0, x1, y1 = self.box
        return (x0 + d, y0 + d, x1 - d, y1 - d)

    def describe(self) -> str:
        return (f"window ({self.x0:g}, {self.y0:g}) + {self.xsize:g} x "
                f"{self.ysize:g} m, clearance {self.clearance:g} m")


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
    #: The child window that carved the parent's plaza, when that window is not
    #: this preset's own child.  ``None`` -- the V1 case -- means "this preset
    #: defines the parent geometry itself".  V2 sets it to the V1 child window,
    #: so every sweep point regenerates the *same* 228-cube parent layout and
    #: the V1 parent's field dumps can be reused unchanged.  ``Sweep.validate``
    #: turns that into a checked invariant.
    plaza: Optional[PlazaWindow] = None
    #: Remove from the **child** the cubes that would fall in its guard + ramp
    #: band, leaving them in the parent.
    #:
    #: Parent and child geometry are not required to match (design section 9.4;
    #: V3 and V4 exist precisely to vary them), and a building-free relaxation
    #: zone is the configuration the scheme was designed for.  So rather than
    #: letting a child inherit whatever cubes the parent happens to have near
    #: its boundary, the child simply clears its own zone.  The parent's
    #: buildings still reach the child: their wakes are in the velocity field
    #: that is imposed on the boundary.  What the child does not carry is solid
    #: cells where ``W > 0``, so ``nest_lparentgeom`` can stay ``.false.`` and
    #: ``nesting_init`` asserts the rule instead of warning about it.
    #:
    #: ``False`` -- the V1 default -- leaves the child's layout as the parent's
    #: restriction, which is right when the plaza was carved for this child and
    #: is what makes the child an exact sub-model of the parent.
    clear_child_zone: bool = False
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
        """True when no cube of **the child's own** layout intrudes into its zone.

        Computed, not declared, and computed from the geometry the child
        actually carries.  With ``clear_child_zone`` the answer is always true
        by construction: whatever the parent has near the child's boundary, the
        child does not carry it.  Without it, the answer depends on the parent
        -- true for V1, where the plaza was carved for this child, and false for
        a smaller or wider-zoned child that inherits the parent's cubes.
        """
        return self.clear_child_zone or len(self.cubes_in_zone()) == 0

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

    @property
    def plaza_window(self) -> PlazaWindow:
        """The window the **parent's** plaza was carved for.

        Defaults to this preset's own child, which is what V1 does.  V2 points
        override it with the V1 child window so that they describe the parent
        that is already on disk.
        """
        if self.plaza is not None:
            return self.plaza
        x0, y0 = self.child_origin
        return PlazaWindow(x0, y0, self.child_xlen, self.child_ylen,
                           self.zone_clearance)

    @property
    def owns_parent_geometry(self) -> bool:
        """True when the parent's plaza was carved for *this* child."""
        return self.plaza is None

    def _interior_box(self) -> Tuple[float, float, float, float]:
        """The building-free interior the **parent's** plaza was cut to leave."""
        return self.plaza_window.interior_box

    def _zone_box(self) -> Tuple[float, float, float, float]:
        """Inside this child's own guard + ramp + ``nwall`` band.

        Identical to :meth:`_interior_box` when this preset owns the parent
        geometry; different -- and the whole point -- when it does not.
        """
        x0, y0, x1, y1 = self._child_box()
        d = self.zone_clearance
        return x0 + d, y0 + d, x1 - d, y1 - d

    def _analysis_interior_box(self) -> Tuple[float, float, float, float]:
        """The box ``analyse.interior_indices`` compares over (no ``nwall``)."""
        x0, y0, x1, y1 = self._child_box()
        d = self.guardwidth + self.zonewidth
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
        cx0, cy0, cx1, cy1 = self.plaza_window.box
        ix0, iy0, ix1, iy1 = self.plaza_window.interior_box
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

    def cubes_in_zone(self) -> np.ndarray:
        """Cubes of the parent's layout that intrude into *this* child's zone.

        A cube counts when it overlaps the child at all but is not wholly
        inside the child's guard + ramp + ``nest_nwall`` band -- exactly the
        condition ``nesting_init`` turns into an abort when
        ``nest_lparentgeom = .false.``.
        """
        half = 0.5 * self.building_width
        cx0, cy0, cx1, cy1 = self._child_box()
        ix0, iy0, ix1, iy1 = self._zone_box()
        out = []
        for cx, cy in self.cube_centres():
            x0, x1, y0, y1 = cx - half, cx + half, cy - half, cy + half
            overlaps = (x1 > cx0 + 1.0e-9 and x0 < cx1 - 1.0e-9
                        and y1 > cy0 + 1.0e-9 and y0 < cy1 - 1.0e-9)
            inside = (x0 >= ix0 - 1.0e-9 and x1 <= ix1 + 1.0e-9
                      and y0 >= iy0 - 1.0e-9 and y1 <= iy1 + 1.0e-9)
            if overlaps and not inside:
                out.append((cx, cy))
        return np.asarray(out, dtype=float).reshape(-1, 2)

    def child_cube_centres(self) -> np.ndarray:
        """Centres of the cubes the **child** carries, in child metres.

        The parent's layout restricted to the child window, minus -- when
        ``clear_child_zone`` is set -- the cubes that would fall in the child's
        guard + ramp band.  This is the one description of the child's geometry;
        ``make_child_case`` builds the STL from it and nothing else.
        """
        kept = self.cube_centres_in(self.child_origin[0], self.child_origin[1],
                                    self.child_xlen, self.child_ylen)
        if not self.clear_child_zone or kept.size == 0:
            return kept
        drop = self.cubes_in_zone()
        if drop.size == 0:
            return kept
        x0, y0 = self.child_origin
        dropped = {(round(cx - x0, 9), round(cy - y0, 9)) for cx, cy in drop}
        out = [(cx, cy) for cx, cy in kept
               if (round(cx, 9), round(cy, 9)) not in dropped]
        return np.asarray(out, dtype=float).reshape(-1, 2)

    def child_cubes_removed(self) -> np.ndarray:
        """Cubes the parent has inside the child window that the child drops."""
        return (self.cubes_in_zone() if self.clear_child_zone
                else np.zeros((0, 2), dtype=float))

    @property
    def n_child_cubes_removed(self) -> int:
        return len(self.child_cubes_removed())

    def removed_cubes_reaching_the_interior(self) -> np.ndarray:
        """Removed cubes whose footprint reaches the region compared.

        Clearing the child's zone makes the child stop being an exact sub-model
        of the parent -- but only where it is allowed to: inside the band, where
        the solution is imposed and no criterion is applied.  A removed cube
        that also reached the *analysis interior* would break that, because the
        statistics there would then be taken over two different geometries.  So
        it is checked rather than hoped for; ``validate`` refuses a preset where
        this is non-empty.
        """
        half = 0.5 * self.building_width
        ix0, iy0, ix1, iy1 = self._analysis_interior_box()
        out = [(cx, cy) for cx, cy in self.child_cubes_removed()
               if (cx + half > ix0 + 1.0e-9 and cx - half < ix1 - 1.0e-9
                   and cy + half > iy0 + 1.0e-9 and cy - half < iy1 - 1.0e-9)]
        return np.asarray(out, dtype=float).reshape(-1, 2)

    def cubes_in_analysis_interior(self) -> np.ndarray:
        """Cubes **the child carries** overlapping the region compared."""
        half = 0.5 * self.building_width
        ix0, iy0, ix1, iy1 = self._analysis_interior_box()
        x0, y0 = self.child_origin
        out = [(cx + x0, cy + y0) for cx, cy in self.child_cube_centres()
               if (cx + x0 + half > ix0 + 1.0e-9 and cx + x0 - half < ix1 - 1.0e-9
                   and cy + y0 + half > iy0 + 1.0e-9 and cy + y0 - half < iy1 - 1.0e-9)]
        return np.asarray(out, dtype=float).reshape(-1, 2)

    @property
    def building_clearance_available(self) -> float:
        """Widest zone this child could have and still sit over open ground [m].

        The smallest distance from any of the four lateral faces to the nearest
        cube face, over the cubes the **parent** actually carries.  A zone whose
        ``zone_clearance`` exceeds it necessarily contains buildings.  ``inf``
        when the child holds no cubes at all.
        """
        half = 0.5 * self.building_width
        cx0, cy0, cx1, cy1 = self._child_box()
        best = float("inf")
        for cx, cy in self.cube_centres():
            x0, x1, y0, y1 = cx - half, cx + half, cy - half, cy + half
            if not (x1 > cx0 + 1.0e-9 and x0 < cx1 - 1.0e-9
                    and y1 > cy0 + 1.0e-9 and y0 < cy1 - 1.0e-9):
                continue
            best = min(best, x0 - cx0, cx1 - x1, y0 - cy0, cy1 - y1)
        return best

    @property
    def interior_cells(self) -> int:
        """Cells per side outside the guard + ramp; the free fetch, in cells."""
        return min(self.child_itot, self.child_jtot) - 2 * self.zone_cells

    @property
    def interior_extent_m(self) -> float:
        return self.interior_cells * self.dx

    @property
    def interior_extent_h(self) -> float:
        return self.interior_extent_m / self.building_height

    @property
    def zone_fraction(self) -> float:
        """``(L_imp + L_rel)`` as a fraction of the shorter child side.

        ``nesting_init`` prints a warning above 0.15 (``modnesting.f90:198``).
        At the small end of the V2b sweep it fires, and that is expected: the
        warning is about how much of the domain the zone eats, not about the
        zone being wrong.
        """
        return (self.guardwidth + self.zonewidth) / min(self.child_xlen,
                                                        self.child_ylen)

    @property
    def zone_fraction_warns(self) -> bool:
        return self.zone_fraction > 0.15

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
        if self.nzone > min(self.child_itot, self.child_jtot):
            errors.append(
                f"nzone = {self.nzone} exceeds the child domain "
                f"({self.child_itot} x {self.child_jtot} cells); the writer refuses it"
            )
        if self.geometry not in ("plaza", "uniform"):
            errors.append(f"geometry must be 'plaza' or 'uniform', got {self.geometry!r}")
        for cx, cy in self.removed_cubes_reaching_the_interior():
            errors.append(
                f"clearing the child's zone would remove a cube at ({cx:g}, {cy:g}) m "
                "whose footprint also reaches the analysis interior; the statistics "
                "would be taken over two different geometries"
            )
        if self.geometry == "plaza" and not self.owns_parent_geometry \
                and not self.clear_child_zone and self.cubes_in_zone().size:
            # Not an error -- 'uniform' means to do exactly this -- but on an
            # inherited plaza it is almost always an oversight, so say so.
            print(f"[config] note: preset '{self.name}' inherits {len(self.cubes_in_zone())} "
                  "of the parent's cubes into its relaxation zone and does not clear "
                  "them; it will need nest_lparentgeom = .true.")
        if self.geometry == "plaza" and self.owns_parent_geometry:
            # When the plaza was carved for *this* child, a cube in the zone is
            # a bug in the preset: the whole point of the plaza is that
            # nest_lparentgeom = .false. can be switched on, so check here that
            # the assertion will pass rather than discovering it three stages
            # later.  When the plaza was carved for a *different* child -- a V2
            # point reusing the V1 parent -- a cube in the zone is not a bug but
            # a measured property of the configuration, reported by
            # `building_free_zone` and honoured by nest_lparentgeom.
            for cx, cy in self.cubes_in_zone():
                errors.append(
                    f"a cube at ({cx:g}, {cy:g}) m survives inside the child but not "
                    "inside its interior; the zone would not be building-free"
                )
        if len(self.cubes_in_analysis_interior()) == 0:
            errors.append(
                "no buildings anywhere in the region the statistics are taken over; "
                "the child would be an empty box and the experiment would test nothing"
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
            f"plaza {'own child' if self.owns_parent_geometry else self.plaza_window.describe()}",
            f"child geometry       {len(self.child_cube_centres())} cubes"
            + (f", {self.n_child_cubes_removed} of the parent's cleared from the "
               f"child's zone (child != parent sub-model inside the band)"
               if self.n_child_cubes_removed else
               " = the parent's restriction, cube for cube"),
            f"zone geometry        {'building-free' if self.building_free_zone else 'CONTAINS buildings'}, "
            f"clearance needed {self.zone_clearance:g} m of "
            f"{self.building_clearance_available:g} m the parent leaves, "
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
            f"{(self.child_itot - 2 * self.zone_cells) * self.dx / self.building_height:.2f}h x "
            f"{(self.child_jtot - 2 * self.zone_cells) * self.dy / self.building_height:.2f}h, "
            f"{len(self.cubes_in_analysis_interior())} cubes",
            f"zone fraction        {100 * self.zone_fraction:.1f} % of the shorter side "
            f"({'nesting_init WARNS above 15 %' if self.zone_fraction_warns else 'no warning'})",
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

# --------------------------------------------------------------------------- #
# V2 -- the falsification sweep (design section 10.4 row V2)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SweepPoint:
    """One child of a sweep, and which arm(s) of the table it belongs to."""

    key: str
    #: ``"zone"`` (V2a, zone width at fixed child size), ``"size"`` (V2b, child
    #: size at fixed zone), or both -- the reference point is shared.
    arms: Tuple[str, ...]
    preset: Preset
    #: ``True`` when this child has already been run and is to be reused rather
    #: than repeated.  The driver takes its case directory from ``--reuse-dir``.
    reuse: bool = False
    note: str = ""

    @property
    def expnr(self) -> str:
        return self.preset.child_expnr


@dataclass(frozen=True)
class Sweep:
    """A set of children driven by **one** parent run, and how to read them.

    Every point shares the parent grid, geometry, forcing and schedule of
    ``parent``; :meth:`validate` checks that rather than trusting it, because
    the whole economy of V2 rests on the parent's field dumps being reusable.
    """

    name: str
    parent: Preset
    points: Tuple[SweepPoint, ...]
    #: Side, in child cells, of the central block over which **every** point is
    #: additionally compared.  Interiors shrink faster than domains, so the
    #: per-point interiors are not the same region; this one is, which
    #: separates "smaller measurement window" from "shorter fetch".
    common_block_cells: int = 0

    def point(self, key: str) -> SweepPoint:
        for p in self.points:
            if p.key == key:
                return p
        raise KeyError(f"no sweep point {key!r} in '{self.name}'; "
                       f"have {', '.join(p.key for p in self.points)}")

    def arm(self, arm: str) -> List[SweepPoint]:
        return [p for p in self.points if arm in p.arms]

    @property
    def to_run(self) -> List[SweepPoint]:
        return [p for p in self.points if not p.reuse]

    def validate(self) -> None:
        """Fail loudly on a sweep whose points cannot share one parent run."""
        errors: List[str] = []
        base = self.parent
        base.validate()
        seen_keys, seen_expnr = set(), {}
        for pt in self.points:
            q = pt.preset
            try:
                q.validate()
            except ValueError as exc:
                errors.append(str(exc))
                continue
            if pt.key in seen_keys:
                errors.append(f"duplicate sweep key {pt.key!r}")
            seen_keys.add(pt.key)
            if q.child_expnr in seen_expnr and seen_expnr[q.child_expnr] != pt.key:
                errors.append(
                    f"points {seen_expnr[q.child_expnr]!r} and {pt.key!r} share "
                    f"child_expnr {q.child_expnr}; their case directories would collide"
                )
            seen_expnr[q.child_expnr] = pt.key
            # -- what makes the parent reusable ---------------------------- #
            for field in ("itot", "jtot", "ktot", "dx", "building_height",
                          "building_width", "street_width", "edgelength",
                          "geometry", "ustar", "u0", "spinup", "production",
                          "dtdump", "child_spinup", "nprocx", "nprocy", "dtmax",
                          "parent_expnr", "stride", "tau", "guardwidth",
                          "nwall", "timeinterp", "init_from_parent"):
                if getattr(q, field) != getattr(base, field):
                    errors.append(
                        f"point {pt.key!r}: {field} = {getattr(q, field)!r} differs from "
                        f"the parent preset's {getattr(base, field)!r}; the parent's "
                        "dumps would not describe this child"
                    )
            if q.plaza_window != base.plaza_window:
                errors.append(
                    f"point {pt.key!r}: plaza {q.plaza_window} does not match the "
                    f"parent's {base.plaza_window}; it would regenerate a different "
                    "cube layout from the one the parent ran"
                )
            if len(q.cube_centres()) != len(base.cube_centres()):
                errors.append(
                    f"point {pt.key!r}: {len(q.cube_centres())} cubes against the "
                    f"parent's {len(base.cube_centres())}"
                )
        if self.common_block_cells:
            for pt in self.points:
                if self.common_block_cells > pt.preset.interior_cells:
                    errors.append(
                        f"common_block_cells = {self.common_block_cells} does not fit "
                        f"inside point {pt.key!r}'s {pt.preset.interior_cells}-cell interior"
                    )
        if not any(p.arms.count("zone") for p in self.points):
            errors.append("the sweep has no 'zone' arm")
        if not any(p.arms.count("size") for p in self.points):
            errors.append("the sweep has no 'size' arm")
        if errors:
            raise ValueError(
                f"sweep '{self.name}' is inconsistent:\n  " + "\n  ".join(errors)
            )

    def summary(self) -> str:
        h = self.parent.building_height
        lines = [
            f"sweep '{self.name}': {len(self.points)} points "
            f"({len(self.to_run)} to run, {len(self.points) - len(self.to_run)} reused), "
            f"one parent '{self.parent.name}' ({self.parent.parent_expnr})",
            f"common comparison block {self.common_block_cells} cells = "
            f"{self.common_block_cells * self.parent.dx / h:.2f}h",
            "",
            f"{'key':10s} {'arms':11s} {'nr':4s} {'child':9s} {'N_imp+N_rel':12s} "
            f"{'nzone':6s} {'interior':16s} {'zone':22s} {'run':6s}",
        ]
        for pt in self.points:
            q = pt.preset
            lines.append(
                f"{pt.key:10s} {'+'.join(pt.arms):11s} {q.child_expnr:4s} "
                f"{q.child_itot:3d}x{q.child_jtot:<5d} "
                f"{int(q.guardwidth / q.dx):3d}+{int(q.zonewidth / q.dx):<8d} "
                f"{q.nzone:<6d} "
                f"{q.interior_cells:3d} cells {q.interior_extent_h:5.2f}h  "
                f"{'clear' if q.building_free_zone else 'BUILDINGS':9s} "
                f"{100 * q.zone_fraction:4.1f}%  "
                f"cut {q.n_child_cubes_removed:<3d}  "
                f"{'reuse' if pt.reuse else 'run':6s}"
            )
        return "\n".join(lines)


def _sweep_child(base: Preset, *, name: str, child_expnr: str,
                 zonewidth: Optional[float] = None,
                 child_cells: Optional[int] = None) -> Preset:
    """One sweep point, expressed as a delta on ``base``.

    ``base`` also fixes the parent, so ``plaza`` is pinned to *its* window: the
    generated cube layout is the one the parent on disk actually ran, whatever
    this child's own zone would have asked for.  ``nzone`` follows the zone
    width, because the nesting file has to store at least as many cells as the
    weights are nonzero over.
    """
    zonewidth = base.zonewidth if zonewidth is None else float(zonewidth)
    n = base.child_itot if child_cells is None else int(child_cells)
    trial = replace(base, name=name, child_expnr=child_expnr,
                    zonewidth=zonewidth, child_itot=n, child_jtot=n,
                    plaza=base.plaza_window, clear_child_zone=True)
    return replace(trial, nzone=trial.zone_cells)


def _v2_sweep(base: Preset, name: str,
              zone_cells_ramp: Sequence[int],
              child_sizes: Sequence[int]) -> Sweep:
    """Build the V2 sweep as deltas on ``base``, which is also its parent.

    ``base`` is the reference point: it sits in **both** arms, and it is the
    child that has already been run, so it is marked ``reuse``.
    """
    dx = base.dx
    ref_rel = int(round(base.zonewidth / dx))
    ref_size = base.child_itot
    points = [SweepPoint(
        key="ref", arms=("zone", "size"), preset=base, reuse=True,
        note=f"the V1 child: N_rel = {ref_rel} cells at {ref_size} x {ref_size}")]
    expnr = int(base.child_expnr)
    for nrel in zone_cells_ramp:
        if nrel == ref_rel:
            continue
        expnr += 1
        points.append(SweepPoint(
            key=f"nrel{nrel}", arms=("zone",),
            preset=_sweep_child(base, name=f"{name}-nrel{nrel}",
                                child_expnr=str(expnr), zonewidth=nrel * dx),
            note=f"N_rel = {nrel} cells at the reference child size"))
    for size in child_sizes:
        if size == ref_size:
            continue
        expnr += 1
        points.append(SweepPoint(
            key=f"size{size}", arms=("size",),
            preset=_sweep_child(base, name=f"{name}-size{size}",
                                child_expnr=str(expnr), child_cells=size),
            note=f"{size} x {size} child at the reference zone"))
    common = min(p.preset.interior_cells for p in points)
    return Sweep(name=name, parent=base, points=tuple(points),
                 common_block_cells=common)


#: **V2 -- the falsification sweep.**  Design section 10.4 row V2, reframed by
#: section 10.5 as a test of the *fetch* interpretation of the V1 TKE deficit:
#:
#:   P1  if the deficit is fetch-limited, the zone width should barely move it;
#:   P2  the child domain size should move it a lot.
#:
#: Both arms hang off the V1 ``converged`` child, which is therefore the shared
#: reference point and is **reused**, not repeated -- so every point is driven
#: by numerically identical parent forcing.
#:
#: Why the zone ramp stops at 16 cells and not the 20 of the design table.  The
#: parent on disk carries the V1 plaza, which leaves 40 m of open ground inside
#: each lateral face of the 128-cell child.  A zone needs
#: ``L_imp + L_rel + nest_nwall*dx`` of it, so ``N_rel = 16`` (38 m + one cell
#: of margin) is the widest ramp that still sits over open ground.  ``N_rel =
#: 20`` would need 48 m: it would have to either re-run the parent -- which
#: breaks the identical-forcing property that makes this comparison clean -- or
#: put buildings in the relaxation zone, which changes the boundary treatment
#: and so confounds exactly the variable under test.  4 -> 16 is still a factor
#: of four in ``L_rel`` and 7 -> 19 cells in total zone thickness, which is
#: ample lever for P1.
#:
#: The size arm keeps a building-free zone the same way every other point does:
#: by clearing it in the **child**.  Parent and child geometry are not required
#: to match (design section 9.4), so the 64- and 96-cell children simply do not
#: carry the cubes that would fall in their guard + ramp band, while the parent
#: keeps them and goes on imprinting them on the child through the imposed
#: velocity field -- their wakes are in the flow that arrives at the boundary.
#: Every point therefore runs ``nest_lparentgeom = .false.`` and child size is
#: the only variable moving along the arm, which is what makes P2 a sharp test
#: rather than corroboration.  ``Preset.removed_cubes_reaching_the_interior``
#: turns "the two geometries differ only inside the band" into a checked
#: invariant; ``validate`` refuses a preset where a cleared cube would also
#: reach the region the statistics are taken over.
V2 = _v2_sweep(CONVERGED, "v2", zone_cells_ramp=(4, 9, 12, 16),
               child_sizes=(64, 96, 128))

#: Parent and reference child for the tiny sweep.  Deliberately **not**
#: ``TINY``.
#:
#: The size arm has to clear cubes out of the smaller children's zones, and for
#: a cleared cube to stay out of the analysis interior -- the invariant
#: ``removed_cubes_reaching_the_interior`` enforces -- the zone has to be at
#: least as deep as a cube is far from the child's face.  With a 32 m array and
#: child origins on the period, that distance is always 24 m: cube centres sit
#: at 16 mod 32 and faces at 0 mod 32, so every cube occupies 8-24 m in from
#: each face.  ``CONVERGED``'s 24 m zone clears it exactly; ``TINY``'s 14 m one
#: does not, and a 32-cell tiny child would have had its interior geometry
#: changed by the clearing.
#:
#: So the tiny sweep carries the **production zone on a tiny domain** rather
#: than a tiny zone, which is the right thing for a smoke test anyway: it is the
#: production code path that wants exercising.  ``TINY`` itself is untouched, so
#: ``test_v1_tiny`` is unaffected.
TINY_SWEEP = Preset(
    name="tiny-sweep",
    itot=128, jtot=128, ktot=32, dx=2.0,
    building_height=16.0, building_width=16.0, street_width=16.0, edgelength=16.0,
    geometry="plaza",
    child_itot=96, child_jtot=96,
    guardwidth=6.0, zonewidth=18.0, tau=1.0, nzone=12, nwall=1, timeinterp=1,
    ustar=0.4, u0=3.0, tke0=0.1,
    spinup=120.0, production=120.0, dtdump=3.0, child_spinup=40.0,
    nprocx=2, nprocy=2, child_nprocx=2, child_nprocy=2,
    dtmax=0.5,
    spectra_heights=(8.0, 16.0, 32.0),
    stride=1,
)

#: The same sweep a few minutes instead of a few hours, exercising the identical
#: code path.  Its size arm really does clear cubes out of the child (the
#: production one clears 12 and 20; this one clears some too), so the "child is
#: not the parent's restriction" path is covered, and both its narrow-domain and
#: wide-zone points make the 15 % zone-fraction warning fire -- all before the
#: production job is submitted.
V2_TINY = _v2_sweep(TINY_SWEEP, "v2-tiny", zone_cells_ramp=(4, 9, 12),
                    child_sizes=(64, 96))

SWEEPS: Dict[str, Sweep] = {s.name: s for s in (V2_TINY, V2)}


def get_sweep(name: str) -> Sweep:
    try:
        sweep = SWEEPS[name]
    except KeyError:
        raise SystemExit(
            f"unknown sweep {name!r}; choose one of {', '.join(sorted(SWEEPS))}"
        ) from None
    sweep.validate()
    return sweep


PRESETS: Dict[str, Preset] = {p.name: p
                              for p in (TINY, TINY_SWEEP, PRODUCTION, CONVERGED)}
PRESETS.update({pt.preset.name: pt.preset
                for sweep in SWEEPS.values() for pt in sweep.points})


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
    import argparse

    ap = argparse.ArgumentParser(description="print the presets and the sweeps")
    ap.add_argument("--sweep", default=None,
                    help="print one sweep's table instead of every preset")
    ns = ap.parse_args()
    if ns.sweep:
        sweep = get_sweep(ns.sweep)
        print(sweep.summary())
        print()
        for pt in sweep.points:
            print(f"--- {pt.key} ({'reused' if pt.reuse else 'to run'}): {pt.note}")
            print(pt.preset.summary())
            print()
    else:
        for name in sorted(PRESETS):
            print(get_preset(name).summary())
            print()
        for name in sorted(SWEEPS):
            print(get_sweep(name).summary())
            print()
