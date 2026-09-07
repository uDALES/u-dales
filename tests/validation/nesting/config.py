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

from collections import OrderedDict
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
    #: nest_timeinterp.  1 = linear, 2 = cubic Hermite with Catmull-Rom slopes
    #: (``modnesting.f90``, ``nest_interp_time``).  Mode 2 is **unlimited**: the
    #: interpolant is a fixed linear combination of the four buffered levels,
    #: so a set of levels whose net flux is individually zero interpolates to
    #: zero net flux and the divergence compatibility of design section 3.1
    #: survives between parent levels.  That was not true of the earlier
    #: Fritsch-Carlson *monotone* Hermite (README.md, "Finding N1": its
    #: nonlinear slope limiter gave Phi = 5.2e-5 and tripped nest_lfluxassert);
    #: the limiter is gone, and 2 is now a legitimate choice.  V1/V2 ran 1; C0
    #: compares the two at the same cadence.
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
    #: Cadence of the boundary data handed to the child [s].  ``make_child_case``
    #: reads every ``cadence / dtdump``-th parent dump level and never touches
    #: the rest, so the nesting file -- and the I/O of building it -- shrink
    #: with it.  Must be a positive whole multiple of ``dtdump``.  ``None``
    #: means ``dtdump``, i.e. every level, which is what everything before C0
    #: ran at.  C0 sweeps it: by Taylor's hypothesis a boundary sampled every
    #: ``cadence`` seconds and interpolated linearly carries nothing below the
    #: wavelength ``2 U cadence``, and the dump Courant number
    #: ``C_dump = U cadence / dx`` (:meth:`dump_courant`) says whether that is
    #: below the parent's own filter scale (``C_dump <= 2``) or not.
    cadence: Optional[float] = None
    #: Field-dump interval of the **child** [s]; ``None`` means ``dtdump``.  A
    #: parent dumped at 0.5 s (C0b) does not need children dumping at 0.5 s for
    #: statistics that were always taken at 3 s; :attr:`analysis_parent_stride`
    #: keeps parent and child analysed at the same sampling whatever the two
    #: dump intervals are.  Must be a positive whole multiple of ``dtdump``.
    child_dtdump: Optional[float] = None
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
    #: What the parent writes during its production window, and therefore
    #: what ``make_child_case`` can build the child from:
    #:
    #: ``"fielddump"`` -- full-domain 3-D dumps of ``u, v, w`` every ``dtdump``
    #:                    (``&OUTPUT lfielddump``), the V1/V2/C0 path.  Needed
    #:                    whenever the parent is also the *reference* the child
    #:                    is compared against, or a V0 driver is box-filtered.
    #: ``"nestdump"``  -- only the child's band and one initial block
    #:                    (``&NESTDUMP``, ``src/modnestdump.f90``, plan item D1):
    #:                    what a driving parent needs to write at a fine cadence.
    #: ``"both"``      -- the two side by side, at the same ``dtdump``; the
    #:                    bit-identity check of ``test_nestdump_tiny`` uses it.
    parent_output: str = "fielddump"
    #: experiment numbers
    parent_expnr: str = "903"
    child_expnr: str = "904"

    def __post_init__(self) -> None:
        # Frozen, so the two defaults that depend on another field are filled
        # in here.  ``dataclasses.replace`` passes the *resolved* values on, so
        # a preset derived from another with a new ``dtdump`` has to set these
        # explicitly (C0_FINE does).
        if self.cadence is None:
            object.__setattr__(self, "cadence", float(self.dtdump))
        if self.child_dtdump is None:
            object.__setattr__(self, "child_dtdump", float(self.dtdump))

    # -- derived ----------------------------------------------------------- #

    @property
    def cadence_stride(self) -> int:
        """Parent dump levels per boundary level: ``round(cadence / dtdump)``."""
        return max(1, int(round(self.cadence / self.dtdump)))

    @property
    def child_dump_stride(self) -> int:
        """Parent dump levels per child dump level."""
        return max(1, int(round(self.child_dtdump / self.dtdump)))

    @property
    def analysis_parent_stride(self) -> int:
        """Stride ``analyse`` applies to the **parent** dumps.

        ``stride`` applies to the child's; the parent's is scaled by the ratio
        of the two dump intervals so that both runs are sampled at
        ``stride * child_dtdump`` seconds.  Equal to ``stride`` whenever the
        two dump intervals agree, which is every preset before C0b.
        """
        return self.stride * self.child_dump_stride

    def dump_courant(self, u: float) -> float:
        """Dump Courant number ``C_dump = u * cadence / dx`` at wind speed ``u``.

        The boundary is sampled every ``cadence`` seconds; with Taylor's
        hypothesis that removes every wavelength below ``2 u cadence`` from the
        imposed field.  Nothing the parent resolved is lost when that is at or
        below its own filter scale ``4 dx``, i.e. ``C_dump <= 2``.  Evaluate it
        at the largest wind in the zone, not the bulk value.
        """
        return u * self.cadence / self.dx

    @property
    def c_dump_u0(self) -> float:
        """:meth:`dump_courant` at the preset's initial bulk wind ``u0``.

        A single-number label for a preset, not the criterion: ``u0`` is the
        initial uniform velocity, and the mean wind aloft in the converged run
        is above it (3.6 m/s at z/h = 2 against ``u0 = 3``).
        """
        return self.dump_courant(self.u0)

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


    # -- the parent-side zone dump (&NESTDUMP, src/modnestdump.f90) --------- #

    @property
    def writes_fielddump(self) -> bool:
        return self.parent_output in ("fielddump", "both")

    @property
    def writes_nestdump(self) -> bool:
        return self.parent_output in ("nestdump", "both")

    def nestdump_nzone(self, child: Optional["Preset"] = None) -> int:
        """Band thickness the parent has to dump, in **parent** cells.

        The child's stored zone (``child.nzone`` cells of ``child.dx``) or its
        guard + ramp, whichever is wider, expressed on this preset's grid and
        rounded up, **plus one cell**: the linear tangential reconstruction of
        ``udprep.nesting.conservative_interpolate`` takes the slope of the
        outermost zone cell from its inner neighbour, so a refined child needs
        one parent cell more than the zone itself.  ``child`` defaults to this
        preset's own child (ratio 1), for which the margin is simply spare.
        """
        c = self if child is None else child
        width = max(c.nzone * c.dx, c.guardwidth + c.zonewidth)
        return int(np.ceil(width / self.dx - 1.0e-9)) + 1

    def nestdump_sections(self, child: Optional["Preset"] = None
                          ) -> "OrderedDict[str, object]":
        """The ``&NESTDUMP`` block for this parent, boxed on ``child``'s window.

        ``child`` (default: this preset's own child) must sit on the same
        physical window -- a V0 driver is a ``dataclasses.replace`` of the fine
        preset with a coarser mesh, so it does -- and the box must land on this
        grid's faces, which the solver checks again and refuses otherwise.
        """
        c = self if child is None else child
        x0, y0 = c.child_origin
        for label, value in (("x0", x0), ("y0", y0),
                             ("xsize", c.child_xlen), ("ysize", c.child_ylen)):
            if abs(value / self.dx - round(value / self.dx)) > 1.0e-9:
                raise ValueError(
                    f"nestdump {label} = {value} m is not on the {self.dx} m parent grid")
        return OrderedDict([
            ("lnestdump", True),
            ("tnestdump", float(self.dtdump)),
            ("nestdump_x0", float(x0)),
            ("nestdump_y0", float(y0)),
            ("nestdump_xsize", float(c.child_xlen)),
            ("nestdump_ysize", float(c.child_ylen)),
            ("nestdump_nzone", int(self.nestdump_nzone(c))),
            ("nestdump_linit", True),
        ])

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
            errors.append("timeinterp must be 1 (linear) or 2 (Catmull-Rom Hermite)")
        for label, value in (("cadence", self.cadence), ("child_dtdump", self.child_dtdump)):
            if not (value > 0):
                errors.append(f"{label} = {value} s is not positive")
            elif abs(value / self.dtdump - round(value / self.dtdump)) > 1.0e-9 \
                    or round(value / self.dtdump) < 1:
                errors.append(
                    f"{label} = {value} s is not a whole multiple of dtdump = "
                    f"{self.dtdump} s; the parent dumps cannot be subsampled to it"
                )
        if self.stride < 1:
            errors.append(f"stride = {self.stride} must be at least 1")
        if self.nzone > min(self.child_itot, self.child_jtot):
            errors.append(
                f"nzone = {self.nzone} exceeds the child domain "
                f"({self.child_itot} x {self.child_jtot} cells); the writer refuses it"
            )
        if self.geometry not in ("plaza", "uniform"):
            errors.append(f"geometry must be 'plaza' or 'uniform', got {self.geometry!r}")
        if self.parent_output not in ("fielddump", "nestdump", "both"):
            errors.append(
                f"parent_output must be 'fielddump', 'nestdump' or 'both', got "
                f"{self.parent_output!r}")
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
            f"({'linear' if self.timeinterp == 1 else 'Catmull-Rom cubic Hermite, unlimited'})",
            f"child init           {'from the parent block' if self.init_from_parent else 'from prof.inp'}",
            f"schedule             spin-up {self.spinup:g} s, production "
            f"[{self.t_start:g}, {self.t_end:g}] s, dump every {self.dtdump:g} s",
            f"                     child spin-up {self.child_spinup:g} s, statistics over "
            f"{self.production - self.child_spinup:g} s",
            f"boundary cadence     {self.cadence:g} s = every "
            f"{self.cadence_stride}{'st' if self.cadence_stride == 1 else 'nd' if self.cadence_stride == 2 else 'rd' if self.cadence_stride == 3 else 'th'} "
            f"parent dump level ({self.production / self.cadence:.0f} levels); "
            f"C_dump = u0 * cadence / dx = {self.c_dump_u0:.2f} at u0 = {self.u0:g} m/s "
            f"({'<= 2, nothing resolved is lost' if self.c_dump_u0 <= 2.0 else '> 2, the parent-resolved band below 2 U cadence is lost at the boundary'})",
            f"child dumps          every {self.child_dtdump:g} s; analysis samples parent "
            f"every {self.analysis_parent_stride} level(s), child every {self.stride}",
            f"ranks                parent {self.nprocx} x {self.nprocy}, "
            f"child {self.child_nprocx} x {self.child_nprocy}",
        ]
        return "\n".join(lines)

    def describe(self) -> str:
        """Alias of :meth:`summary`: the human-readable block, one line per fact."""
        return self.summary()


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
    #: size at fixed zone), ``"cadence"`` (C0, boundary-data cadence at the
    #: reference child), or several -- the reference point is shared.
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
    #: The arms this sweep is read along; each must have at least one point.
    #: ``sweep_summary.ARMS`` says how each is tabulated and plotted.
    arms: Tuple[str, ...] = ("zone", "size")

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
            # Not in the list, because they are the child's business and the
            # parent's dumps describe the child whatever they are: timeinterp
            # (how the child interpolates between the levels it is given),
            # cadence (which of the parent's levels it is given) and
            # child_dtdump.  C0 sweeps the first two.
            for field in ("itot", "jtot", "ktot", "dx", "building_height",
                          "building_width", "street_width", "edgelength",
                          "geometry", "ustar", "u0", "spinup", "production",
                          "dtdump", "child_spinup", "nprocx", "nprocy", "dtmax",
                          "parent_expnr", "stride", "tau", "guardwidth",
                          "nwall", "init_from_parent"):
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
        if not self.arms:
            errors.append("the sweep declares no arms")
        for arm in self.arms:
            if not any(arm in p.arms for p in self.points):
                errors.append(f"the sweep has no {arm!r} arm")
        for pt in self.points:
            for arm in pt.arms:
                if arm not in self.arms:
                    errors.append(f"point {pt.key!r} is in arm {arm!r}, which the "
                                  f"sweep does not declare ({', '.join(self.arms)})")
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
            f"{'nzone':6s} {'interior':16s} {'zone':22s} "
            f"{'cadence':9s} {'C_dump':7s} {'interp':7s} {'run':6s}",
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
                f"{q.cadence:<5g} s  {q.c_dump_u0:<7.2f} "
                f"{'linear' if q.timeinterp == 1 else 'CR':7s} "
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

# --------------------------------------------------------------------------- #
# C0 -- the cadence discriminator (nesting-plan-2026-09-06.md section 1)
# --------------------------------------------------------------------------- #
#
# Experiment numbers in use across this directory, so a new sweep does not
# collide with a run already on disk (case directories are named by them):
#
#   903 / 904        V1 parent / child (config.PRODUCTION, CONVERGED, TINY, ...)
#   905-909          V2 children (config.V2, V2_TINY)
#   911, 912         V0 coarse driving parents; 921-924 V0 children (config.V0)
#   920-922, 930-931, 940-954   V3 / V4 geometry cases (presets_geometry.py)
#   960-969          C0: 960 the fine-cadence parent (C0_FINE); 961-963 the C0a
#                    children (cad6, cad9, cr3); 964-969 the C0b children
#                    (0.5, 1, 1.5, 3, 6, 9 s)
#   970-975          C0c: the Catmull-Rom ladder off the same fine parent
#                    (0.5, 1, 1.5, 3, 6, 9 s, nest_timeinterp = 2)
#
# 910 and 913-919 are free but sit between V0's two blocks; C0 takes the next
# clear decade instead.


def _c0_sweep(base: Preset, name: str, *, cadences: Sequence[float],
              cr_cadences: Sequence[float], expnrs: Sequence[str],
              reuse_reference: bool) -> Sweep:
    """A cadence ladder of children off **one** parent, each child otherwise
    ``base``'s child exactly -- same window, zone, geometry, forcing, init.

    ``reuse_reference``: ``base`` is a child that has already been run (V1's)
    and is the sweep's ``ref`` point; a point at ``base.cadence`` is then not
    generated again.  Otherwise every point runs, and the point at
    ``base.cadence`` *is* ``base`` (so it carries ``base.child_expnr``).
    ``cr_cadences`` adds Catmull-Rom (``timeinterp = 2``) points; everything
    else is linear.  ``expnrs`` are consumed in order by the generated points.
    """
    nrs = iter(expnrs)
    points: List[SweepPoint] = []
    if reuse_reference:
        points.append(SweepPoint(
            key="ref", arms=("cadence",), preset=base, reuse=True,
            note=f"the V1 child: {base.cadence:g} s cadence, linear "
                 f"(C_dump = {base.c_dump_u0:.2f} at u0)"))
    for c in cadences:
        key = f"cad{c:g}"
        if abs(c - base.cadence) < 1.0e-9:
            if reuse_reference:
                continue
            preset = base
        else:
            preset = replace(base, name=f"{name}-{key}", child_expnr=next(nrs),
                             cadence=float(c), timeinterp=1)
        points.append(SweepPoint(
            key=key, arms=("cadence",), preset=preset,
            note=f"{c:g} s cadence, linear (C_dump = {preset.c_dump_u0:.2f} at u0), "
                 f"every {preset.cadence_stride} parent level(s)"))
    for c in cr_cadences:
        key = f"cr{c:g}"
        preset = replace(base, name=f"{name}-{key}", child_expnr=next(nrs),
                         cadence=float(c), timeinterp=2)
        points.append(SweepPoint(
            key=key, arms=("cadence",), preset=preset,
            note=f"{c:g} s cadence, Catmull-Rom cubic in time "
                 f"(C_dump = {preset.c_dump_u0:.2f} at u0)"))
    common = min(p.preset.interior_cells for p in points)
    return Sweep(name=name, parent=base, points=tuple(points),
                 common_block_cells=common, arms=("cadence",))


#: **C0a -- coarser cadences from the existing 3 s dumps.**  The review of
#: 2026-09-06 attributes the V1 TKE deficit above the canopy not to fetch but
#: to the 3 s boundary cadence: sampled every ``cadence`` seconds and
#: interpolated linearly, the boundary carries nothing below ``2 U cadence``
#: (21.5 m at z/h = 2), which is exactly the band the child was short of.  The
#: cheap half of the test subsamples the converged parent's dumps to 6 s and
#: 9 s -- if the cadence causes the deficit these must be *worse* than V1 --
#: and adds the unlimited Catmull-Rom interpolant at 3 s, which cannot restore
#: a band the samples do not contain and so should leave 8-16 m unchanged.
#: Three children through the V2 machinery, off the same 903 bundle, against
#: the same reused 904 reference.  Pre-registered predictions in README.md.
C0 = _c0_sweep(CONVERGED, "c0", cadences=(6.0, 9.0), cr_cadences=(3.0,),
               expnrs=("961", "962", "963"), reuse_reference=True)

#: **C0b -- the fine-cadence ladder.**  The converged parent left its
#: end-of-spin-up restart (``initd00031204_*.903``, t = 10800 s).  This preset
#: warm-starts it for 2400 s dumping every 0.5 s (dt is about 0.38 s, so every
#: 1-2 steps; 4800 levels), and the ladder below slices those dumps to 0.5, 1,
#: 1.5, 3, 6 and 9 s.  One parent realisation drives all six, so the
#: comparison is paired, and the 3 s point cross-checks V1 and C0a.  The
#: children keep the 600 s discard and get an 1800 s window, enough for the
#: band ratios (reproduced to 0.004 between the 1491 s and 10191 s windows)
#: though not for the profile deficit -- which is why the band ratios are the
#: primary metric.  The child dumps every 3 s as V1's did; the analysis samples
#: the parent every 6th level to match (``analysis_parent_stride``).
#:
#: Memory: the slab cut holds the whole nesting file, 9.8 MB per level for the
#: 128^2 child at nzone = 12, so the 0.5 s point is 47 GB of slabs.  The
#: writer stores each slab variable separately and the flux correction works
#: in place, so the peak is that plus one slab (V1 measured 39 GB for a 35 GB
#: file); mem=128gb holds it.
C0_FINE = replace(
    CONVERGED, name="c0-fine", parent_expnr="960", child_expnr="964",
    production=2400.0, dtdump=0.5, child_spinup=600.0,
    cadence=0.5, child_dtdump=3.0, stride=1,
)

C0B = _c0_sweep(C0_FINE, "c0b", cadences=(0.5, 1.0, 1.5, 3.0, 6.0, 9.0),
                cr_cadences=(), expnrs=("965", "966", "967", "968", "969"),
                reuse_reference=False)

#: **C0c -- the same ladder with the Catmull-Rom interpolant.**  C0a found
#: that ``nest_timeinterp = 2`` at 3 s halves the deficit V1 measured (-9.9 %
#: to -5.7 %; 8-16 m ratio 0.83 to 0.89), which the pre-registered prediction
#: had not allowed for: the better-supplied 16-64 m band feeds the cascade
#: that rebuilds 8-16 m.  This sweep runs the six C0b cadences again with the
#: cubic interpolant off the same 960 dumps (symlink the 960 case directory
#: into the run directory and the parent is not re-run), so that the design
#: can state the operating curve for the interpolant it recommends.
C0C = _c0_sweep(C0_FINE, "c0c", cadences=(),
                cr_cadences=(0.5, 1.0, 1.5, 3.0, 6.0, 9.0),
                expnrs=("970", "971", "972", "973", "974", "975"),
                reuse_reference=False)

#: The same two sweeps in minutes, on the ``tiny`` parent (dtdump 3 s, 120 s
#: window).  ``c0-tiny`` reuses the tiny V1 child as ``ref`` exactly as the
#: production sweep reuses 904; ``c0b-tiny`` warm-starts a 0.5 s parent from
#: the tiny spin-up's restart and runs three children off it.
C0_TINY = _c0_sweep(TINY, "c0-tiny", cadences=(6.0, 9.0), cr_cadences=(3.0,),
                    expnrs=("961", "962", "963"), reuse_reference=True)

C0_FINE_TINY = replace(
    TINY, name="c0-fine-tiny", parent_expnr="960", child_expnr="964",
    production=120.0, dtdump=0.5, child_spinup=40.0,
    cadence=0.5, child_dtdump=3.0, stride=1,
)

C0B_TINY = _c0_sweep(C0_FINE_TINY, "c0b-tiny", cadences=(0.5, 1.5, 3.0),
                     cr_cadences=(), expnrs=("965", "966"),
                     reuse_reference=False)

SWEEPS: Dict[str, Sweep] = {s.name: s for s in (V2_TINY, V2, C0_TINY, C0,
                                                 C0B_TINY, C0B, C0C)}


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
                              for p in (TINY, TINY_SWEEP, PRODUCTION, CONVERGED,
                                        C0_FINE, C0_FINE_TINY)}
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


# --------------------------------------------------------------------------- #
# V0 -- validation at refinement ratio > 1 (design section 10.4 row V0)
# --------------------------------------------------------------------------- #
#
# V1 and V2 both run parent and child on the SAME grid, so the conservative
# interpolation in tools/python/udprep/nesting.py has never carried a running
# simulation -- only its own unit tests (P1-P17).  V0 is the first end-to-end
# test of refinement AND the first end-to-end test of that interpolation, and
# the configuration below is arranged so that when something moves it is
# possible to say which of the two moved it.


@dataclass(frozen=True)
class RefinedPoint:
    """One refined child, and where the data that drives it comes from.

    A refined experiment lives on **two** grids, so it needs two
    :class:`Preset` objects and both of them are real:

    ``child``
        the **fine** side: the child's own grid, geometry, zone, forcing and
        schedule -- and, through ``parent_expnr``, the fine reference run the
        child is *measured against*.  It is an ordinary ratio-1 preset, so
        ``analyse.run`` compares child and reference on identical grids by
        exactly the code V1 used, and the V1 result is literally the ``r = 1``
        row of the same table.
    ``driver``
        the **coarse** side: the grid the boundary data lives on.  Same
        physical domain, same cubes, same forcing, spacing ``refine`` times
        coarser.  For the ``coarse`` arm it is a case that is actually built
        and run; for the ``filtered`` arm nothing is run and it serves only to
        describe the grid the fine reference's own dumps are filtered onto.

    Separating "what drives the child" from "what the child is compared to" is
    the whole design of V0.  The reference is the **fine** run in both arms --
    comparing a 2 m child against a 4 m or 8 m parent would score the child
    down for resolving turbulence its parent cannot represent, which is not an
    error but the point of nesting.

    The two arms measure different things and must be reported separately:

    ``filtered``
        the driving data is the fine reference's own field, box-filtered onto
        the coarse grid.  The filter is flux-conservative, so a discretely
        solenoidal fine field gives a discretely solenoidal coarse one, and the
        coarse field is a *perfect* coarse parent -- it knows exactly what the
        fine run was doing at the scales it can represent.  This is an
        idealisation: it hands the child filtered fine-scale information a
        genuinely coarse LES would never have had.  What it isolates is the
        prolongation and the parent's filter scale, and nothing else -- it is
        V1 with one variable changed, driven by the same realisation of the
        same turbulence, so it is paired with the reference exactly as V1 was.
    ``coarse``
        the driving data comes from a genuinely coarse LES of the same domain.
        This is the real use case.  It is *not* paired: the coarse run is an
        independent realisation, so its eddies are not the reference's eddies
        and only statistics can be compared.  It also carries the coarse run's
        own biases -- a 16 m cube is 4 cells wide at ``r = 2`` and 2 at
        ``r = 4``, and its drag will not be the 8-cell version's.  A mean-flow
        error measured here is therefore the parent's error plus the nesting's,
        and the two are separated by the driving parent's own profile, which
        ``make_child_case`` records while it cuts the slabs.

    The difference between the arms at the same ``refine`` is the quantity of
    interest: how much the child suffers from its parent genuinely not knowing
    the small scales, as opposed to from the interpolation.
    """

    key: str
    #: ``"filtered"`` or ``"coarse"``; see the class docstring.
    arm: str
    #: Spatial refinement ratio ``dx_parent / dx_child``.  Integer, and at most
    #: ``udprep.nesting.MAX_SPATIAL_REFINEMENT`` (4), which the writer enforces.
    refine: int
    driver: Preset
    child: Preset
    note: str = ""

    # -- derived ------------------------------------------------------------ #

    @property
    def expnr(self) -> str:
        """Experiment number of the child."""
        return self.child.child_expnr

    @property
    def reference_expnr(self) -> str:
        """Experiment number of the fine run the child is compared against."""
        return self.child.parent_expnr

    @property
    def driver_expnr(self) -> str:
        return self.driver.parent_expnr

    @property
    def runs_driver(self) -> bool:
        """True when the driving parent is a case of its own that must be run."""
        return self.arm == "coarse"

    @property
    def coarsen(self) -> int:
        """Factor the *driving directory's* dumps are filtered by before use.

        ``refine`` for the ``filtered`` arm, where the driving directory holds
        the fine reference's dumps; 1 for the ``coarse`` arm, where they are
        already on the coarse grid.
        """
        return self.refine if self.arm == "filtered" else 1

    @property
    def parent_nyquist_wavelength(self) -> float:
        """The parent's own filter scale, ``2 dx_parent`` [m].

        The wavelength either side of which the child is doing two different
        jobs: below it the child must *generate* structure the parent never
        resolved, above it the child is reproducing structure the parent had.
        Section 10.5's deficit sits in the 8-64 m band, which straddles this at
        both ratios (16 m at ``r = 2``, 32 m at ``r = 4``), which is why the
        band ratios have to be split here and not merely quoted.
        """
        return 2.0 * self.driver.dx

    @property
    def resolves_the_ramp(self) -> bool:
        """Design section 1.4(c): ``L_rel >= 2 dx_parent``.

        The transition has to be resolved in the *parent's* terms, not only in
        the child's.  Reported rather than enforced -- a preset that violates it
        is a legitimate (if unfavourable) configuration, and the tiny smoke test
        does violate it at ``r = 4``.
        """
        return self.child.zonewidth >= 2.0 * self.driver.dx - 1.0e-9

    def describe(self) -> str:
        c, d = self.child, self.driver
        return (f"{self.key:14s} r = {self.refine} ({self.arm:8s}) "
                f"parent {d.itot:3d}x{d.jtot:3d}x{d.ktot:3d} @ {d.dx:g} m"
                f"{' (filtered from ' + self.reference_expnr + ')' if not self.runs_driver else ' (run as ' + d.parent_expnr + ')'}"
                f"  ->  child {c.child_itot}x{c.child_jtot}x{c.child_ktot} @ "
                f"{c.dx:g} m ({c.child_expnr})")

    # -- self-consistency --------------------------------------------------- #

    def validate(self) -> None:
        errors: List[str] = []
        c, d = self.child, self.driver
        c.validate()
        d.validate()
        if self.arm not in ("filtered", "coarse"):
            errors.append(f"arm must be 'filtered' or 'coarse', got {self.arm!r}")
        if self.refine < 1 or self.refine != int(self.refine):
            errors.append(f"refine must be a positive integer, got {self.refine!r}")
        # udprep.nesting.MAX_SPATIAL_REFINEMENT; not imported, to keep config.py
        # free of solver-tooling imports, but it is the same number and
        # write_nesting_file enforces it.
        if self.refine > 4:
            errors.append(f"refine = {self.refine} exceeds the writer's validated "
                          "maximum spatial refinement of 4")
        if abs(d.dx - self.refine * c.dx) > 1.0e-9:
            errors.append(f"driver dx = {d.dx} is not {self.refine} x the child's {c.dx}")
        # Same physical box, or the driving field does not cover the child.
        for label, a, b in (("xlen", d.xlen, c.xlen), ("ylen", d.ylen, c.ylen),
                            ("zsize", d.zsize, c.zsize)):
            if abs(a - b) > 1.0e-9:
                errors.append(f"driver {label} = {a} does not match the reference's {b}")
        # Same child window, in metres.
        for label, a, b in (("origin x", d.child_origin[0], c.child_origin[0]),
                            ("origin y", d.child_origin[1], c.child_origin[1]),
                            ("xlen", d.child_xlen, c.child_xlen),
                            ("ylen", d.child_ylen, c.child_ylen)):
            if abs(a - b) > 1.0e-9:
                errors.append(f"the driver's child window {label} = {a} m does not "
                              f"match the child's {b} m")
        # Same geometry.  The plaza is carved in metres, so the two grids must
        # agree cube for cube or the child's buildings are not the parent's.
        if d.plaza_window != c.plaza_window:
            errors.append(f"driver plaza {d.plaza_window.describe()} does not match "
                          f"the reference's {c.plaza_window.describe()}")
        if len(d.cube_centres()) != len(c.cube_centres()):
            errors.append(f"driver carries {len(d.cube_centres())} cubes against the "
                          f"reference's {len(c.cube_centres())}")
        elif d.cube_centres().size and not np.allclose(
                np.sort(d.cube_centres(), axis=0), np.sort(c.cube_centres(), axis=0)):
            errors.append("driver and reference cube layouts differ")
        # Same forcing and the same schedule: the child has to see one momentum
        # source and one time axis, whichever arm drives it.
        for field in ("building_height", "building_width", "street_width",
                      "geometry", "ustar", "u0", "tke0", "spinup", "production",
                      "dtdump", "child_spinup", "dtmax"):
            if getattr(d, field) != getattr(c, field):
                errors.append(f"driver {field} = {getattr(d, field)!r} differs from the "
                              f"reference's {getattr(c, field)!r}")
        if abs(d.dpdx - c.dpdx) > 1.0e-15:
            errors.append(f"driver dpdx = {d.dpdx} differs from the reference's {c.dpdx}")
        if self.arm == "coarse" and d.parent_expnr == c.parent_expnr:
            errors.append(f"the coarse driving parent and the fine reference share "
                          f"expnr {d.parent_expnr}; their case directories would collide")
        if self.refine == 1:
            errors.append("a V0 point at refine = 1 is V1; use the V1 preset instead")
        if not self.resolves_the_ramp:
            # Not an error: it is a real, reportable property of the point.
            print(f"[config] note: point '{self.key}' has L_rel = {c.zonewidth:g} m "
                  f"< 2 dx_parent = {2 * d.dx:g} m, so the relaxation ramp is not "
                  "resolved in the parent's own terms (design section 1.4c)")
        if errors:
            raise ValueError(
                f"refinement point '{self.key}' is inconsistent:\n  " + "\n  ".join(errors)
            )


@dataclass(frozen=True)
class RefinementSuite:
    """The V0 points, and the one fine reference run they are all measured against."""

    name: str
    #: The fine, unnested run that supplies the truth -- and, for the
    #: ``filtered`` arm, the field that is filtered to drive the child.  For the
    #: production suite this is the V1 ``converged`` parent already on disk.
    reference: Preset
    points: Tuple[RefinedPoint, ...]

    def point(self, key: str) -> RefinedPoint:
        for p in self.points:
            if p.key == key:
                return p
        raise KeyError(f"no V0 point {key!r} in '{self.name}'; "
                       f"have {', '.join(p.key for p in self.points)}")

    def arm(self, arm: str) -> List[RefinedPoint]:
        return [p for p in self.points if p.arm == arm]

    def ratio(self, refine: int) -> List[RefinedPoint]:
        return [p for p in self.points if p.refine == refine]

    @property
    def refinements(self) -> Tuple[int, ...]:
        return tuple(sorted({p.refine for p in self.points}))

    @property
    def drivers_to_run(self) -> List[Preset]:
        """The coarse parent cases that have to be built and run, one per ratio."""
        out: List[Preset] = []
        seen = set()
        for p in self.points:
            if p.runs_driver and p.driver.parent_expnr not in seen:
                seen.add(p.driver.parent_expnr)
                out.append(p.driver)
        return out

    def validate(self) -> None:
        errors: List[str] = []
        self.reference.validate()
        seen_keys, seen_expnr, drivers = set(), {}, {}
        for pt in self.points:
            try:
                pt.validate()
            except ValueError as exc:
                errors.append(str(exc))
                continue
            if pt.key in seen_keys:
                errors.append(f"duplicate point key {pt.key!r}")
            seen_keys.add(pt.key)
            if pt.expnr in seen_expnr:
                errors.append(f"points {seen_expnr[pt.expnr]!r} and {pt.key!r} share "
                              f"child expnr {pt.expnr}; their directories would collide")
            seen_expnr[pt.expnr] = pt.key
            # Every child must be the reference's child, or `analyse.run` cannot
            # compare it against the reference's own dumps.
            for field in ("itot", "jtot", "ktot", "dx", "building_height",
                          "building_width", "street_width", "edgelength",
                          "geometry", "child_itot", "child_jtot", "ustar", "u0",
                          "spinup", "production", "dtdump", "child_spinup",
                          "guardwidth", "zonewidth", "tau", "nzone", "nwall",
                          "timeinterp", "init_from_parent", "clear_child_zone",
                          "child_nprocx", "child_nprocy", "dtmax", "stride",
                          "parent_expnr"):
                if getattr(pt.child, field) != getattr(self.reference, field):
                    errors.append(
                        f"point {pt.key!r}: child {field} = "
                        f"{getattr(pt.child, field)!r} differs from the reference's "
                        f"{getattr(self.reference, field)!r}; the reference run's dumps "
                        "would not describe this child"
                    )
            if pt.child.plaza_window != self.reference.plaza_window:
                errors.append(f"point {pt.key!r}: child plaza does not match the "
                              "reference's; it would carry different buildings")
            # One coarse grid per ratio, shared by both arms, so the filtered and
            # coarse arms of a ratio differ ONLY in where the data came from.
            prev = drivers.get(pt.refine)
            if prev is not None:
                for field in ("itot", "jtot", "ktot", "dx", "nprocx", "nprocy"):
                    if getattr(prev, field) != getattr(pt.driver, field):
                        errors.append(
                            f"point {pt.key!r}: driver {field} differs from the other "
                            f"r = {pt.refine} point's; the two arms would not be "
                            "comparing like with like"
                        )
            drivers[pt.refine] = pt.driver
        for refine in sorted(drivers):
            if not any(p.refine == refine and p.arm == "filtered" for p in self.points):
                errors.append(f"r = {refine} has no 'filtered' arm")
            if not any(p.refine == refine and p.arm == "coarse" for p in self.points):
                errors.append(f"r = {refine} has no 'coarse' arm")
        if errors:
            raise ValueError(
                f"refinement suite '{self.name}' is inconsistent:\n  "
                + "\n  ".join(errors)
            )

    def summary(self) -> str:
        r = self.reference
        lines = [
            f"suite '{self.name}': {len(self.points)} points at r = "
            f"{', '.join(str(x) for x in self.refinements)}, arms "
            f"{', '.join(sorted({p.arm for p in self.points}))}",
            f"reference (truth)    '{r.name}' ({r.parent_expnr}): "
            f"{r.itot} x {r.jtot} x {r.ktot} @ {r.dx:g} m, "
            f"{r.xlen:g} x {r.ylen:g} x {r.zsize:g} m",
            f"child (every point)  {r.child_itot} x {r.child_jtot} x {r.child_ktot} @ "
            f"{r.dx:g} m, origin ({r.child_origin[0]:g}, {r.child_origin[1]:g}) m, "
            f"zone {r.zone_cells} cells = {r.guardwidth + r.zonewidth:g} m",
            "coarse parents to run "
            + (", ".join("%s @ %g m" % (p.parent_expnr, p.dx)
                         for p in self.drivers_to_run) or "none"),
            "",
        ]
        for pt in self.points:
            lines.append("  " + pt.describe())
            lines.append(f"                 parent Nyquist wavelength "
                         f"{pt.parent_nyquist_wavelength:g} m; L_rel = "
                         f"{pt.child.zonewidth:g} m "
                         f"({'>=' if pt.resolves_the_ramp else '<'} 2 dx_parent = "
                         f"{2 * pt.driver.dx:g} m); {pt.note}")
        return "\n".join(lines)


def _v0_driver(base: Preset, refine: int, *, expnr: str,
               nprocx: int, nprocy: int) -> Preset:
    """The coarse parent grid for one refinement ratio, as a delta on ``base``.

    Same physical domain, same cubes, same forcing and the same schedule; only
    the mesh is coarser.  ``plaza`` is pinned to ``base``'s window so that the
    coarse run regenerates the layout the fine reference actually ran -- the
    default would compute the plaza from the coarse grid's own ``nest_nwall``
    margin and carve a slightly larger one, and the child's buildings would
    then not be the parent's.
    """
    return replace(
        base,
        name=f"{base.name}-parent-r{refine}",
        itot=base.itot // refine, jtot=base.jtot // refine, ktot=base.ktot // refine,
        dx=base.dx * refine,
        child_itot=base.child_itot // refine, child_jtot=base.child_jtot // refine,
        plaza=base.plaza_window,
        parent_expnr=expnr, nprocx=nprocx, nprocy=nprocy,
    )


def _v0_suite(base: Preset, name: str, *,
              ranks: Dict[int, Tuple[int, int]],
              driver_expnr: Dict[int, str],
              child_expnr: Dict[Tuple[int, str], str]) -> RefinementSuite:
    """Build a V0 suite: both arms at every ratio in ``ranks``, one grid per ratio."""
    points: List[RefinedPoint] = []
    for refine in sorted(ranks):
        nprocx, nprocy = ranks[refine]
        driver = _v0_driver(base, refine, expnr=driver_expnr[refine],
                            nprocx=nprocx, nprocy=nprocy)
        for arm in ("filtered", "coarse"):
            nr = child_expnr[(refine, arm)]
            child = replace(base, name=f"{name}-r{refine}-{arm}", child_expnr=nr)
            points.append(RefinedPoint(
                key=f"r{refine}-{arm}", arm=arm, refine=refine,
                driver=driver, child=child,
                note=("boundary data box-filtered from the fine reference's own dumps: "
                      "a perfect coarse parent, paired with the truth"
                      if arm == "filtered" else
                      f"boundary data from a genuine {driver.dx:g} m LES "
                      f"({driver.parent_expnr}): the real use case, unpaired")))
    return RefinementSuite(name=name, reference=base, points=tuple(points))


#: **V0 -- the refinement validation.**
#:
#: The child is *exactly* the V1 converged child -- 128 x 128 x 64 cells at 2 m
#: over the same 256 x 256 x 128 m box, the same 36 cubes, the same 3 + 9 cell
#: zone, the same forcing and the same 10 800 s window -- and the reference it
#: is compared against is the V1 converged parent, which is already on disk.
#: Only the parent's mesh changes.  So the V1 result is the ``r = 1`` row of
#: this table, measured by the same code over the same window, and every number
#: here is directly readable against section 10.5.
#:
#: The two coarse grids are 128^2 x 32 at 4 m and 64^2 x 16 at 8 m: one eighth
#: and one sixty-fourth of the fine parent's cells, so the coarse arm's two
#: parent runs are cheap next to the children they drive.  At 4 m a 16 m cube is
#: 4 cells wide and at 8 m it is 2, which is a real limitation of the coarse
#: parent and is the point -- V0 asks what a child can recover from a parent
#: like that, not from a good one.
#:
#: L_rel = 18 m clears 2 dx_parent at both ratios (8 m and 16 m), so the
#: relaxation ramp is resolved in the parent's own terms as design section
#: 1.4(c) requires; nothing about the zone changes between V1 and V0.
V0 = _v0_suite(
    CONVERGED, "v0",
    ranks={2: (8, 8), 4: (4, 4)},
    driver_expnr={2: "911", 4: "912"},
    child_expnr={(2, "filtered"): "921", (4, "filtered"): "922",
                 (2, "coarse"): "923", (4, "coarse"): "924"},
)

#: The same suite in minutes rather than hours, on the ``tiny`` fine reference,
#: which this one has to *run* as well (the production suite reuses V1's).
#:
#: Two things are deliberately unlike production and are reported rather than
#: hidden: ``tiny``'s L_rel = 8 m is below ``2 dx_parent`` at both ratios, so
#: the ramp is not resolved in the parent's terms, and at r = 4 the coarse
#: parent is 24 x 24 x 8 cells with 2-cell cubes.  Neither matters for a smoke
#: test -- what is being exercised is the code path -- and both would matter a
#: great deal for a physical claim, which this preset does not make.
V0_TINY = _v0_suite(
    TINY, "v0-tiny",
    ranks={2: (2, 2), 4: (2, 2)},
    driver_expnr={2: "911", 4: "912"},
    child_expnr={(2, "filtered"): "921", (4, "filtered"): "922",
                 (2, "coarse"): "923", (4, "coarse"): "924"},
)

SUITES: Dict[str, RefinementSuite] = {s.name: s for s in (V0_TINY, V0)}


def get_suite(name: str) -> RefinementSuite:
    try:
        suite = SUITES[name]
    except KeyError:
        raise SystemExit(
            f"unknown refinement suite {name!r}; choose one of {', '.join(sorted(SUITES))}"
        ) from None
    suite.validate()
    return suite


PRESETS.update({p.name: p for s in SUITES.values()
                for pt in s.points for p in (pt.driver, pt.child)})


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="print the presets, sweeps and suites")
    ap.add_argument("--sweep", default=None,
                    help="print one sweep's table instead of every preset")
    ap.add_argument("--suite", default=None,
                    help="print one refinement suite's table instead of every preset")
    ns = ap.parse_args()
    if ns.suite:
        suite = get_suite(ns.suite)
        print(suite.summary())
        print()
        for pt in suite.points:
            print(f"--- {pt.key} ({pt.arm}, r = {pt.refine}): {pt.note}")
            print(pt.child.summary())
            print()
        for driver in suite.drivers_to_run:
            print(f"--- coarse driving parent {driver.parent_expnr}")
            print(driver.summary())
            print()
    elif ns.sweep:
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
        for name in sorted(SUITES):
            print(get_suite(name).summary())
            print()
