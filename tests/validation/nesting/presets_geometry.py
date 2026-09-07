#!/usr/bin/env python3
"""Presets for **V3** and **V4** -- the experiments where parent and child
geometry differ (``docs/udales-nesting-design.md`` section 9.4, section 10.4
rows V3 and V4).

This module imports the shared machinery from :mod:`config` and registers its
own presets there; it deliberately adds nothing to ``config.py`` itself, which
V0 is editing in parallel.

What the two experiments are
----------------------------

**V3 -- parent without buildings.**  The parent resolves no geometry at all: a
flat, periodic, neutral channel.  The child carries a cube canopy that starts
immediately at the inner edge of its relaxation zone.  Section 9.4 predicts
that an internal boundary layer must develop from that point, and that a
building-free *standoff* between the zone and the first building row is
"actively counterproductive" because the flow would then adjust twice.  V3
measures the **adjustment length** and compares standoffs of 0, 5, 15 and 40
cells so that the claim can fail.

**V4 -- different parent geometry.**  The parent carries a *staggered* cube
array; the child carries the *aligned* array of V1, at V1's size, zone,
forcing and schedule.  Everything is V1's except the parent's layout, so the
V1 ``converged`` child is the reference and the answer is a difference from it.

Three things this module does **not** do, on purpose
----------------------------------------------------

1. **It does not carve a plaza.**  A building-free zone is a constraint on the
   *child* alone (section 9.4, and the ``clear_child_zone`` work in
   ``config.py``): only the child's solid mask has to be clear where ``W > 0``,
   which is what ``nest_lparentgeom = .false.`` asserts.  Every child here
   generates its own layout inside its own clear box, so the parent is never
   modified for the child's benefit.  For V3 the parent has no buildings at
   all; for V4 the parent keeps its staggered array intact, zone included.
2. **It does not require the child's lattice to be in phase with the
   parent's.**  For V3 there is no parent lattice to be in phase with, which is
   exactly the freedom that lets the standoff be varied in cells rather than in
   whole cube periods.
3. **It does not treat the parent sub-region as truth.**  In both experiments
   the parent has different buildings from the child inside the compared
   region, so ``analyse.py``'s child-versus-parent numbers -- criterion A
   included -- measure the *geometry difference*, not a nesting error.  The
   references are elsewhere: a periodic equilibrium run for V3, the V1
   ``converged`` child for V4.  See ``analyse_geometry.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Tuple

import numpy as np

import config
from config import Preset

# --------------------------------------------------------------------------- #
# Geometry primitives
# --------------------------------------------------------------------------- #


def aligned_centres(xlen: float, ylen: float, width: float, street: float
                    ) -> np.ndarray:
    """Centres of the unbroken aligned array on ``[0, xlen] x [0, ylen]``.

    Same formula and ordering as ``udgeom.create_cubes(..., 'AC')`` and as
    ``config.Preset._full_cube_centres``: ``c = i (C + H) - H/2 - C/2``.
    """
    p = width + street
    nx, ny = int(round(xlen / p)), int(round(ylen / p))
    out = [(i * p - 0.5 * width - 0.5 * street, j * p - 0.5 * width - 0.5 * street)
           for i in range(1, nx + 1) for j in range(1, ny + 1)]
    return np.asarray(out, dtype=float).reshape(-1, 2)


def staggered_centres(xlen: float, ylen: float, width: float, street: float
                      ) -> np.ndarray:
    """Centres of the staggered array, exactly as ``create_cubes(..., 'SC')``.

    Every second column (constant ``x``) is displaced half a period in ``y``,
    so a cube sits behind the *gap* of the column upstream of it -- the
    standard staggered cube array, at the same plan area density as the aligned
    one.  ``create_cubes`` writes the displaced columns as ``Ny + 1`` cubes at
    ``y = 0, p, 2p, ..., ylen`` and clips the two at the periodic faces into
    halves that join across the boundary; both are returned here as nominal
    whole-cube centres, because the only thing this list is used for is
    restricting the array to a *child* window (``cube_centres_in`` keeps a cube
    only when its whole footprint is inside), where the clipped pair never
    reaches.  The parent's STL comes from ``create_cubes`` itself, not from
    this list, so the clipping is done by the generator that owns it.
    """
    p = width + street
    nx, ny = int(round(xlen / p)), int(round(ylen / p))
    out: List[Tuple[float, float]] = []
    for i in range(1, nx + 1):
        x = i * p - 0.5 * width - 0.5 * street
        if i % 2 == 0:
            out.append((x, 0.0))
            out.extend((x, j * p) for j in range(1, ny + 1))
        else:
            out.extend((x, j * p - 0.5 * width - 0.5 * street)
                       for j in range(1, ny + 1))
    return np.asarray(out, dtype=float).reshape(-1, 2)


LAYOUTS = {"none": None, "aligned": aligned_centres, "staggered": staggered_centres}


# --------------------------------------------------------------------------- #
# The preset
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class GeoPreset(Preset):
    """A :class:`config.Preset` whose parent and child layouts are independent.

    ``config.Preset`` describes one array, carried by the parent, of which the
    child's is a restriction (optionally with its zone cleared).  V3 and V4 need
    the two to be genuinely different, so the layout methods are overridden.

    **Coordinate convention.**  ``cube_centres`` stays in *parent* metres, as in
    the base class.  Everything about the child's own array --
    ``child_cube_centres``, ``cubes_in_zone``, ``child_cubes_removed``,
    ``cubes_in_analysis_interior``, ``removed_cubes_reaching_the_interior`` --
    is in *child* metres.  The base class mixes the two (its child array is a
    restriction of the parent's, so both are meaningful); here the child's array
    has no parent counterpart, and child metres are the only frame in which it
    can be described.  Every consumer of those methods either takes a length
    (``len``) or is overridden below.
    """

    #: ``"reference"``  a periodic equilibrium run, with no child at all;
    #: ``"parent"``     the run whose dumps drive a child;
    #: ``"child"``      a nested child.
    role: str = "child"
    #: The layout the **parent** resolves: ``"none"`` (V3 -- a flat parent),
    #: ``"aligned"`` or ``"staggered"`` (V4).
    parent_layout: str = "aligned"
    #: The layout the **child** resolves: ``"aligned"`` or ``"staggered"``.
    child_layout: str = "aligned"
    #: ``"parent"``   the child's lattice is in phase with the parent's global
    #:                lattice, so where the two layouts agree they line up
    #:                cube for cube (V4, and V1's arrangement);
    #: ``"standoff"`` the lattice is placed by :attr:`standoff_cells` instead --
    #:                the first building face sits ``standoff_cells`` cells
    #:                inside the child's clear box, and the spanwise rows are
    #:                centred.  V3 uses this; it is only meaningful when the
    #:                parent has no lattice to be in phase with.
    child_phase: str = "parent"
    #: Building-free standoff between the inner edge of the zone (strictly, the
    #: ``nest_nwall``-eroded clear box) and the first building face, in child
    #: cells.  ``0`` is what section 9.4 recommends.
    standoff_cells: int = 0
    #: When set, the periodic run is driven by a volume-flow-rate controller at
    #: this bulk velocity instead of by ``dpdx``.  Only meaningful for
    #: ``role != "child"``.  ``None`` means fixed ``dpdx``, as V1 uses.
    uflowrate: Optional[float] = None
    #: The child's geometry is never the parent's restriction here -- that is
    #: the whole experiment -- so this stays ``True``.  ``make_child_case`` uses
    #: it to decide whether the seed ``prof.inp`` profile, which it accumulates
    #: from the *parent's* field over the child window, has to be masked with the
    #: parent's solid cells as well as the child's.  Without it a cell that is
    #: fluid in the child and solid in the parent would fold the parent's
    #: near-zero in-building velocity into the seed.  It matters little (the
    #: child cold-starts from the parent's 3-D block anyway) but it is free and
    #: it is right.
    clear_child_zone: bool = True
    #: Fraction of the canopy's spanwise extent, centred, over which the
    #: headline streamwise statistics are taken.  The lateral internal boundary
    #: layers spreading in from the two spanwise zone edges contaminate the
    #: outer part of the patch, increasingly so with fetch; 0.5 keeps the
    #: central half.  ``analyse_geometry`` reports the full width as well, so
    #: the choice is checkable rather than assumed.
    y_core_fraction: float = 0.5

    # -- the parent's array -------------------------------------------------- #

    def cube_centres(self) -> np.ndarray:
        """Centres of the cubes the **parent** carries, in parent metres."""
        fn = LAYOUTS.get(self.parent_layout)
        if self.parent_layout not in LAYOUTS:
            raise ValueError(f"unknown parent_layout {self.parent_layout!r}")
        if fn is None:
            return np.zeros((0, 2), dtype=float)
        return fn(self.xlen, self.ylen, self.building_width, self.street_width)

    @property
    def parent_is_flat(self) -> bool:
        return self.parent_layout == "none"

    @property
    def n_cubes_removed(self) -> int:
        """Nothing is ever removed from the parent here -- that is the point."""
        return 0

    # -- the child's array --------------------------------------------------- #

    def _child_clear_box(self) -> Tuple[float, float, float, float]:
        """The child box inset by ``zone_clearance``, in **child** metres.

        ``L_imp + L_rel`` is where ``W`` reaches zero and ``nest_nwall`` cells
        more because the weights are additionally eroded away from any solid
        point, so this is the region in which a building may legally stand.
        """
        c = self.zone_clearance
        return c, c, self.child_xlen - c, self.child_ylen - c

    def _child_lattice(self) -> np.ndarray:
        """The child's nominal lattice before the clear box is applied.

        ``child_phase = "parent"`` regenerates the parent's global lattice and
        shifts it into child coordinates, so a child whose layout matches its
        parent's lines up with it cube for cube.  ``"standoff"`` places the
        lattice by hand: the first building face ``standoff_cells`` cells inside
        the clear box, the spanwise rows centred in what is left.
        """
        fn = LAYOUTS.get(self.child_layout)
        if fn is None:
            raise ValueError(
                f"child_layout {self.child_layout!r}: a child with no buildings "
                "is not one of these experiments"
            )
        if self.child_phase == "parent":
            centres = fn(self.xlen, self.ylen, self.building_width,
                         self.street_width)
            x0, y0 = self.child_origin
            local = np.column_stack([centres[:, 0] - x0, centres[:, 1] - y0])
            # Only cubes whose whole footprint is in the child window exist for
            # the child at all; the rest are the parent's business.  Restricting
            # here rather than later keeps ``child_cubes_removed`` meaning "in
            # the window but not in the clear box" -- the cubes whose absence
            # could confound the comparison -- rather than "everywhere else in
            # the parent".
            half = 0.5 * self.building_width
            keep = ((local[:, 0] - half >= -1.0e-9)
                    & (local[:, 0] + half <= self.child_xlen + 1.0e-9)
                    & (local[:, 1] - half >= -1.0e-9)
                    & (local[:, 1] + half <= self.child_ylen + 1.0e-9))
            return local[keep].reshape(-1, 2)
        if self.child_phase != "standoff":
            raise ValueError(f"unknown child_phase {self.child_phase!r}")
        if self.child_layout != "aligned":
            raise ValueError("the standoff phase is only defined for an aligned array")
        p, half = self.period, 0.5 * self.building_width
        cx0, cy0, cx1, cy1 = self._child_clear_box()
        x_first = cx0 + self.standoff_cells * self.dx
        nx = int(np.floor((cx1 - x_first) / p + 1.0e-9))
        ny = int(np.floor((cy1 - cy0) / p + 1.0e-9))
        # Centre the spanwise rows in the span the clear box leaves: the two
        # lateral zones are symmetric, so the canopy should be too.
        y_first = cy0 + 0.5 * ((cy1 - cy0) - ny * p)
        out = [(x_first + k * p + half, y_first + m * p + half)
               for k in range(max(nx, 0)) for m in range(max(ny, 0))]
        return np.asarray(out, dtype=float).reshape(-1, 2)

    def _inside_clear_box(self, centres: np.ndarray) -> np.ndarray:
        half = 0.5 * self.building_width
        cx0, cy0, cx1, cy1 = self._child_clear_box()
        keep = []
        for cx, cy in centres:
            if (cx - half >= cx0 - 1.0e-9 and cx + half <= cx1 + 1.0e-9
                    and cy - half >= cy0 - 1.0e-9 and cy + half <= cy1 + 1.0e-9):
                keep.append((cx, cy))
        return np.asarray(keep, dtype=float).reshape(-1, 2)

    def child_cube_centres(self) -> np.ndarray:
        """Centres of the cubes the **child** carries, in child metres.

        The nominal lattice minus every cube that would not fit entirely inside
        the clear box.  That single rule gives a building-free zone by
        construction at every point of both experiments, so
        ``nest_lparentgeom = .false.`` and ``nesting_init`` *asserts* the design
        section 5 rule rather than warning about it.
        """
        return self._inside_clear_box(self._child_lattice())

    def child_cubes_removed(self) -> np.ndarray:
        """Lattice cubes the child does not carry, in child metres."""
        kept = {(round(x, 9), round(y, 9)) for x, y in self.child_cube_centres()}
        out = [(x, y) for x, y in self._child_lattice()
               if (round(x, 9), round(y, 9)) not in kept]
        return np.asarray(out, dtype=float).reshape(-1, 2)

    def cubes_in_zone(self) -> np.ndarray:
        """Cubes of the **child's own** layout that intrude into its zone.

        Empty by construction -- :meth:`child_cube_centres` keeps only cubes
        inside the clear box -- and computed rather than asserted so that
        :attr:`config.Preset.building_free_zone`, which is what decides
        ``nest_lparentgeom``, stays a measurement.
        """
        half = 0.5 * self.building_width
        cx0, cy0, cx1, cy1 = self._child_clear_box()
        out = [(cx, cy) for cx, cy in self.child_cube_centres()
               if not (cx - half >= cx0 - 1.0e-9 and cx + half <= cx1 + 1.0e-9
                       and cy - half >= cy0 - 1.0e-9 and cy + half <= cy1 + 1.0e-9)]
        return np.asarray(out, dtype=float).reshape(-1, 2)

    def _analysis_interior_box_child(self) -> Tuple[float, float, float, float]:
        """What ``analyse.interior_indices`` compares over, in child metres."""
        d = self.guardwidth + self.zonewidth
        return d, d, self.child_xlen - d, self.child_ylen - d

    def cubes_in_analysis_interior(self) -> np.ndarray:
        """Cubes **the child carries** overlapping the region compared."""
        half = 0.5 * self.building_width
        ix0, iy0, ix1, iy1 = self._analysis_interior_box_child()
        out = [(cx, cy) for cx, cy in self.child_cube_centres()
               if (cx + half > ix0 + 1.0e-9 and cx - half < ix1 - 1.0e-9
                   and cy + half > iy0 + 1.0e-9 and cy - half < iy1 - 1.0e-9)]
        return np.asarray(out, dtype=float).reshape(-1, 2)

    def removed_cubes_reaching_the_interior(self) -> np.ndarray:
        """Dropped cubes whose footprint would have reached the compared region.

        Always empty here: a cube is dropped only when it does not fit inside
        the clear box, which is ``nest_nwall`` cells *deeper* than the analysis
        interior, so a dropped cube can never poke into the interior.  Checked
        rather than argued, exactly as ``config.Preset`` checks it.
        """
        half = 0.5 * self.building_width
        ix0, iy0, ix1, iy1 = self._analysis_interior_box_child()
        out = [(cx, cy) for cx, cy in self.child_cubes_removed()
               if (cx + half > ix0 + 1.0e-9 and cx - half < ix1 - 1.0e-9
                   and cy + half > iy0 + 1.0e-9 and cy - half < iy1 - 1.0e-9)]
        return np.asarray(out, dtype=float).reshape(-1, 2)

    # -- the canopy patch, as the streamwise analysis needs it ---------------- #

    @property
    def canopy_x_range(self) -> Optional[Tuple[float, float]]:
        """``(first face, last face)`` of the canopy in child metres."""
        c = self.child_cube_centres()
        if c.size == 0:
            return None
        half = 0.5 * self.building_width
        return float(c[:, 0].min() - half), float(c[:, 0].max() + half)

    @property
    def canopy_y_range(self) -> Optional[Tuple[float, float]]:
        c = self.child_cube_centres()
        if c.size == 0:
            return None
        half = 0.5 * self.building_width
        return float(c[:, 1].min() - half), float(c[:, 1].max() + half)

    @property
    def y_core_range(self) -> Optional[Tuple[float, float]]:
        """The central part of the canopy the headline statistics are taken over."""
        span = self.canopy_y_range
        if span is None:
            return None
        y0, y1 = span
        mid, halfwidth = 0.5 * (y0 + y1), 0.5 * self.y_core_fraction * (y1 - y0)
        return mid - halfwidth, mid + halfwidth

    @property
    def n_rows(self) -> int:
        """Number of streamwise building rows in the canopy.

        One row is one spanwise line of cubes at a common ``x``; the streamwise
        statistics are reduced onto the period-wide block that contains it.
        """
        c = self.child_cube_centres()
        if c.size == 0:
            return 0
        return int(np.unique(np.round(c[:, 0], 6)).size)

    def row_edges(self) -> np.ndarray:
        """Leading face of each streamwise row block, in child metres.

        A row block is one cube period: the cube and the street behind it.  The
        streamwise statistics are reduced onto these blocks rather than onto
        bare ``x`` slabs, because ``<u>(x)`` in a cube array oscillates strongly
        within a period (in front of / on top of / behind a cube) and the
        adjustment is the trend *through* that oscillation.
        """
        span = self.canopy_x_range
        if span is None:
            return np.zeros(0)
        return span[0] + self.period * np.arange(self.n_rows + 1)

    @property
    def standoff_m(self) -> float:
        """Building-free standoff beyond the clear box, in metres."""
        span = self.canopy_x_range
        if span is None:
            return float("nan")
        return span[0] - self.zone_clearance

    @property
    def first_row_fetch_m(self) -> float:
        """Distance from the **inner edge of the zone** to the first building face.

        This is the quantity section 9.4's standoff argument is about: it is the
        domain the layout spends before the canopy starts.  It is
        ``nest_nwall * dx`` larger than :attr:`standoff_m`, because even a
        "zero standoff" layout has to keep the wall-erosion margin clear.
        """
        span = self.canopy_x_range
        if span is None:
            return float("nan")
        return span[0] - (self.guardwidth + self.zonewidth)

    # -- consistency --------------------------------------------------------- #

    def validate(self) -> None:
        errors: List[str] = []
        p = self.period
        if self.role not in ("reference", "parent", "child"):
            errors.append(f"unknown role {self.role!r}")
        if self.parent_layout not in LAYOUTS:
            errors.append(f"unknown parent_layout {self.parent_layout!r}")
        if self.timeinterp not in (1, 2):
            errors.append("timeinterp must be 1 (linear) or 2 (Hermite)")
        if self.building_height % self.dx:
            errors.append("building_height is not a whole number of cells")
        if self.child_spinup >= self.production:
            errors.append("child_spinup leaves no statistics window")
        if self.parent_output not in ("fielddump", "nestdump", "both"):
            errors.append(f"parent_output must be 'fielddump', 'nestdump' or 'both', "
                          f"got {self.parent_output!r}")
        if self.fielddump_dtdump is not None and not (self.fielddump_dtdump > 0):
            errors.append(f"fielddump_dtdump = {self.fielddump_dtdump} s is not positive")
        for label, value in (("cadence", self.cadence), ("child_dtdump", self.child_dtdump)):
            if not (value > 0):
                errors.append(f"{label} = {value} s is not positive")
            elif abs(value / self.dtdump - round(value / self.dtdump)) > 1.0e-9 \
                    or round(value / self.dtdump) < 1:
                errors.append(
                    f"{label} = {value} s is not a whole multiple of dtdump = "
                    f"{self.dtdump} s; the parent dumps cannot be subsampled to it")
        for n, tot, label in ((self.nprocx, self.itot, "nprocx/itot"),
                              (self.nprocy, self.jtot, "nprocy/jtot")):
            if tot % n:
                errors.append(f"{label}: {tot} is not divisible by {n}")
        if not self.parent_is_flat:
            for label, size in (("xlen", self.xlen), ("ylen", self.ylen)):
                if abs(size / p - round(size / p)) > 1.0e-9:
                    errors.append(f"{label} = {size} is not a whole number of "
                                  f"{p} m cube periods")
        if self.role != "reference":
            for n, tot, label in ((self.child_nprocx, self.child_itot,
                                   "child nprocx/itot"),
                                  (self.child_nprocy, self.child_jtot,
                                   "child nprocy/jtot")):
                if tot % n:
                    errors.append(f"{label}: {tot} is not divisible by {n}")
            # FieldDump.child_block needs one cell of parent beyond the child's
            # east/north face for the extra staggered velocity plane.
            if self.child_i0 + self.child_itot + 1 > self.itot:
                errors.append(
                    f"the child window ({self.child_i0} + {self.child_itot}) reaches "
                    f"the parent's east face ({self.itot}); the extra staggered face "
                    "would need the periodic image, which the slab cut does not do")
            if self.child_j0 + self.child_jtot + 1 > self.jtot:
                errors.append(
                    f"the child window ({self.child_j0} + {self.child_jtot}) reaches "
                    f"the parent's north face ({self.jtot})")
            if self.nzone < self.zone_cells:
                errors.append(
                    f"nzone = {self.nzone} is thinner than the {self.zone_cells}-cell "
                    f"zone (L_imp + L_rel = {self.guardwidth + self.zonewidth} m)")
            if self.nzone > min(self.child_itot, self.child_jtot):
                errors.append(f"nzone = {self.nzone} exceeds the child domain")
            if 2 * self.zone_cells >= min(self.child_itot, self.child_jtot):
                errors.append("the zone covers the whole child; no interior is left")
            if self.standoff_cells < 0:
                errors.append("standoff_cells must not be negative")
            if len(self.child_cube_centres()) == 0:
                errors.append(
                    "the child's clear box holds no cubes at all; there is no canopy "
                    "for the flow to adjust to")
            elif self.n_rows < 2:
                errors.append(
                    f"the child's canopy is {self.n_rows} row(s) deep; the streamwise "
                    "adjustment cannot be resolved with fewer than two")
            if len(self.cubes_in_zone()):
                errors.append("a cube of the child's own layout intrudes into its "
                              "zone; nest_lparentgeom = .false. would abort")
            for cx, cy in self.removed_cubes_reaching_the_interior():
                errors.append(
                    f"a cube dropped at ({cx:g}, {cy:g}) m would also have reached the "
                    "analysis interior; the statistics would span two geometries")
            if len(self.cubes_in_analysis_interior()) == 0:
                errors.append("no buildings in the region the statistics are taken "
                              "over; the experiment would test nothing")
        else:
            for label, size in (("itot*dx", self.xlen), ("jtot*dx", self.ylen)):
                if abs(size / p - round(size / p)) > 1.0e-9:
                    errors.append(f"reference {label} = {size} m is not a whole "
                                  f"number of {p} m periods, so it is not periodic")
        if errors:
            raise ValueError(
                f"preset '{self.name}' is inconsistent:\n  " + "\n  ".join(errors))

    # -- reporting ----------------------------------------------------------- #

    @property
    def forcing_description(self) -> str:
        """How the run this preset describes is driven.

        A *child* is always driven by the fixed ``dpdx``
        (``make_child_case.child_sections`` writes ``luvolflowr = .false.``),
        whatever ``uflowrate`` says: a child preset inherits the field from the
        parent it was derived from and never uses it.
        """
        if self.uflowrate is not None and self.role != "child":
            return (f"volume flow rate, uflowrate = {self.uflowrate:g} m/s "
                    f"(dpdx = 0)")
        return f"dpdx = {self.dpdx:.4e} m/s^2 -> ustar = {self.ustar:g} m/s"

    def summary(self) -> str:
        lines = [
            f"preset               {self.name}  [{self.role}]",
            f"parent grid          {self.itot} x {self.jtot} x {self.ktot} cells, "
            f"{self.xlen:g} x {self.ylen:g} x {self.zsize:g} m, dx = {self.dx:g} m",
            f"parent geometry      '{self.parent_layout}': "
            f"{len(self.cube_centres())} cubes",
            f"forcing              {self.forcing_description}",
            f"schedule             spin-up {self.spinup:g} s, production "
            f"[{self.t_start:g}, {self.t_end:g}] s, dtdump {self.dtdump:g} s",
            f"parent output        {self.parent_output}"
            + (f" (&OUTPUT every {self.fielddump_interval:g} s, &NESTDUMP every "
               f"{self.dtdump:g} s)" if self.parent_output == "both"
               and abs(self.fielddump_interval - self.dtdump) > 1.0e-9 else ""),
            f"ranks                {self.nprocx} x {self.nprocy}",
        ]
        if self.role != "reference":
            span_x = self.canopy_x_range
            span_y = self.canopy_y_range
            lines += [
                f"child grid           {self.child_itot} x {self.child_jtot} x "
                f"{self.child_ktot} cells, {self.child_xlen:g} x "
                f"{self.child_ylen:g} m, origin "
                f"({self.child_origin[0]:g}, {self.child_origin[1]:g}) m",
                f"child geometry       '{self.child_layout}' phase "
                f"'{self.child_phase}': {len(self.child_cube_centres())} cubes, "
                f"{self.n_rows} streamwise rows, "
                f"{len(self.child_cubes_removed())} dropped for not fitting the "
                f"clear box",
                f"canopy               x [{span_x[0]:g}, {span_x[1]:g}] m, "
                f"y [{span_y[0]:g}, {span_y[1]:g}] m; standoff "
                f"{self.standoff_cells} cells = {self.standoff_m:g} m beyond the "
                f"clear box, {self.first_row_fetch_m:g} m beyond the ramp",
                f"zone                 L_imp = {self.guardwidth:g} m, L_rel = "
                f"{self.zonewidth:g} m, tau = {self.tau:g} s, nzone = {self.nzone}, "
                f"total {self.zone_cells} cells, D = {self.optical_depth:.2f}",
                f"zone geometry        "
                f"{'building-free' if self.building_free_zone else 'CONTAINS BUILDINGS'}"
                f", nest_lparentgeom = "
                f"{'.false.' if self.building_free_zone else '.true.'}",
                f"interior             streamwise "
                f"{self.child_itot - 2 * self.zone_cells} cells = "
                f"{(self.child_xlen - 2 * (self.guardwidth + self.zonewidth)) / self.building_height:.2f} h, "
                f"spanwise {self.child_jtot - 2 * self.zone_cells} cells = "
                f"{(self.child_ylen - 2 * (self.guardwidth + self.zonewidth)) / self.building_height:.2f} h, "
                f"{len(self.cubes_in_analysis_interior())} cubes; zone fraction "
                f"{100 * self.zone_fraction:.1f} %"
                f"{' (nesting_init WARNS)' if self.zone_fraction_warns else ''}",
                f"child spin-up        {self.child_spinup:g} s discarded, statistics "
                f"over {self.production - self.child_spinup:g} s",
                f"child ranks          {self.child_nprocx} x {self.child_nprocy}",
                f"boundary cadence     {self.cadence:g} s = every "
                f"{self.cadence_stride} parent dump level; "
                f"C_dump = u0 * cadence / dx = {self.c_dump_u0:.2f} at u0 = "
                f"{self.u0:g} m/s "
                f"({'<= 2, nothing resolved is lost' if self.c_dump_u0 <= 2.0 else '> 2, the parent-resolved band below 2 U cadence is lost at the boundary'})",
                f"time interpolation   nest_timeinterp = {self.timeinterp} "
                f"({'linear' if self.timeinterp == 1 else 'Catmull-Rom cubic Hermite, unlimited'})",
                f"child dumps          every {self.child_dtdump:g} s",
            ]
        return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Experiment containers
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PeriodicRun:
    """A periodic run: the equilibrium reference, or a child's parent."""

    key: str
    preset: GeoPreset
    note: str = ""
    #: key of the periodic run whose measured bulk velocity sets this one's
    #: ``uflowrate``.  ``None`` -- fixed ``dpdx``.
    bulk_from: Optional[str] = None


@dataclass(frozen=True)
class ChildRun:
    key: str
    preset: GeoPreset
    parent_key: str
    note: str = ""
    #: run this point unless ``--only`` says otherwise
    default: bool = True


@dataclass(frozen=True)
class Experiment:
    """One of V3 / V4, at production or tiny size."""

    name: str
    kind: str                      # "v3" or "v4"
    periodic: Tuple[PeriodicRun, ...]
    children: Tuple[ChildRun, ...]
    #: V3: which periodic run is the canopy's equilibrium reference.
    equilibrium_key: Optional[str] = None
    #: V4: which child is the matched-geometry baseline.  ``None`` means it is
    #: supplied from outside (the V1 ``converged`` child on disk).
    baseline_child_key: Optional[str] = None
    #: V4: the preset describing that external baseline child, when there is one.
    external_baseline: Optional[Preset] = None
    headline: str = ""

    def periodic_run(self, key: str) -> PeriodicRun:
        for r in self.periodic:
            if r.key == key:
                return r
        raise KeyError(f"no periodic run {key!r} in '{self.name}'")

    def child(self, key: str) -> ChildRun:
        for c in self.children:
            if c.key == key:
                return c
        raise KeyError(f"no child {key!r} in '{self.name}'; have "
                       f"{', '.join(c.key for c in self.children)}")

    @property
    def default_children(self) -> List[ChildRun]:
        return [c for c in self.children if c.default]

    def validate(self) -> None:
        errors: List[str] = []
        keys = [r.key for r in self.periodic]
        if len(set(keys)) != len(keys):
            errors.append(f"duplicate periodic keys {keys}")
        expnrs: Dict[str, str] = {}
        for r in self.periodic:
            try:
                r.preset.validate()
            except ValueError as exc:
                errors.append(str(exc))
            if r.bulk_from is not None and r.bulk_from not in keys:
                errors.append(f"{r.key!r}: bulk_from {r.bulk_from!r} is not a run here")
            if r.bulk_from is not None and r.preset.uflowrate is None:
                errors.append(f"{r.key!r}: bulk_from is set but uflowrate is None, "
                              "so the measured bulk would be ignored")
            if r.preset.parent_expnr in expnrs:
                errors.append(f"expnr {r.preset.parent_expnr} used twice")
            expnrs[r.preset.parent_expnr] = r.key
        ckeys = [c.key for c in self.children]
        if len(set(ckeys)) != len(ckeys):
            errors.append(f"duplicate child keys {ckeys}")
        for c in self.children:
            try:
                c.preset.validate()
            except ValueError as exc:
                errors.append(str(exc))
                continue
            if c.parent_key not in keys:
                errors.append(f"child {c.key!r}: no periodic run {c.parent_key!r}")
                continue
            parent = self.periodic_run(c.parent_key).preset
            # What makes this child's slab cut describe that parent's dumps.
            for f in ("itot", "jtot", "ktot", "dx", "spinup", "production",
                      "dtdump", "parent_expnr", "parent_layout"):
                if getattr(c.preset, f) != getattr(parent, f):
                    errors.append(
                        f"child {c.key!r}: {f} = {getattr(c.preset, f)!r} differs from "
                        f"parent {c.parent_key!r}'s {getattr(parent, f)!r}")
            if c.preset.child_expnr in expnrs:
                errors.append(f"expnr {c.preset.child_expnr} used twice "
                              f"({expnrs[c.preset.child_expnr]} and {c.key})")
            expnrs[c.preset.child_expnr] = c.key
        if self.kind == "v3" and self.equilibrium_key not in keys:
            errors.append("a V3 experiment needs an equilibrium reference run")
        if self.kind == "v4":
            if self.baseline_child_key is None and self.external_baseline is None:
                errors.append("a V4 experiment needs a matched-geometry baseline, "
                              "internal or external")
            if self.baseline_child_key is not None \
                    and self.baseline_child_key not in ckeys:
                errors.append(f"baseline child {self.baseline_child_key!r} is not here")
        if errors:
            raise ValueError(
                f"experiment '{self.name}' is inconsistent:\n  " + "\n  ".join(errors))

    def summary(self) -> str:
        lines = [f"experiment '{self.name}' ({self.kind}): {self.headline}", ""]
        for r in self.periodic:
            f = (f"uflowrate from '{r.bulk_from}'" if r.bulk_from
                 else ("uflowrate fixed" if r.preset.uflowrate is not None
                       else "dpdx"))
            lines.append(f"  periodic {r.key:12s} nr {r.preset.parent_expnr}  "
                         f"{r.preset.itot:4d}x{r.preset.jtot:<4d}x{r.preset.ktot:<3d} "
                         f"layout {r.preset.parent_layout:9s} {f:22s} {r.note}")
        for c in self.children:
            p = c.preset
            lines.append(
                f"  child    {c.key:12s} nr {p.child_expnr}  "
                f"{p.child_itot:4d}x{p.child_jtot:<4d}x{p.child_ktot:<3d} "
                f"layout {p.child_layout + '/' + p.child_phase:18s} "
                f"standoff {p.standoff_cells:2d} cells = {p.standoff_m:5.1f} m  "
                f"{p.n_rows:2d} rows  on '{c.parent_key}'"
                f"{'' if c.default else '  [not run by default]'}")
        return "\n".join(lines)


# --------------------------------------------------------------------------- #
# V3 -- parent without buildings
# --------------------------------------------------------------------------- #

#: Everything V3 shares with V1: 2 m cells, 16 m cubes on a 32 m period, a
#: neutral rigid-lid channel 128 m deep, u* = 0.4 m/s.  Holding these fixed is
#: what makes the V3 canopy the *same* canopy V1 and V2 measured.
#: ``dtdump`` is set per preset below, not here: a "reference" run needs
#: nothing finer than the historical 3 s (nobody drives a child from it), while
#: a "parent" run's ``dtdump`` is the ``&NESTDUMP`` cadence a child is driven
#: from and has to satisfy ``C_dump <= 2`` (nesting-plan-2026-09-06.md section
#: 1, "V3, V4"; docs/udales-nesting-design.md section 10.5's operating rule).
#: ``timeinterp = 2`` (Catmull-Rom, unlimited) is C0c's recommended default.
_V3_COMMON = dict(
    ktot=64, dx=2.0,
    building_height=16.0, building_width=16.0, street_width=16.0, edgelength=16.0,
    geometry="uniform",
    guardwidth=6.0, zonewidth=18.0, tau=1.0, nzone=12, nwall=1, timeinterp=2,
    ustar=0.4, tke0=0.1, dtmax=0.5,
    spectra_heights=(8.0, 16.0, 32.0), stride=1,
)

#: The ``&NESTDUMP`` cadence every V3/V4 parent writes its boundary band at,
#: and the child's ``cadence``: with the mean wind at the domain top ~5.5 m/s
#: (V1's converged run, z/h = 2) and ``dx = 2`` m, ``C_dump <= 2`` needs
#: ``dtdump <= 2 * 2 / 5.5 = 0.73`` s; 0.5 s matches C0's own fine parent
#: (``config.C0_FINE``) and gives ``C_dump = 1.4`` at 5.5 m/s (0.82 at the
#: preset's own ``u0 = 3`` m/s, the number ``summary()`` prints).  The parent's
#: own ``&OUTPUT`` field dump -- read back by ``periodic-stats`` for the bulk
#: velocity and canopy statistics ``run_geometry.py`` needs, never by the
#: child -- stays at the old, statistically-sufficient 3 s
#: (``fielddump_dtdump``), so ``parent_output = "both"`` does not multiply the
#: full-domain dump volume by 6x for no benefit.
_NESTDUMP_CADENCE = dict(dtdump=0.5, cadence=0.5, child_dtdump=3.0,
                         fielddump_dtdump=3.0, parent_output="both")

#: **The equilibrium reference.**  A periodic 6 x 6 array of the child's cubes,
#: driven by the same fixed ``dpdx`` the child is driven by, so its statistics
#: are what "in equilibrium with its own canopy" *means* for the V3 child.
#: Without it the adjustment length could only ever be self-referential -- "the
#: profile stops changing" -- which cannot distinguish an adjusted canopy from a
#: canopy that has stopped adjusting short of equilibrium.
V3_REFERENCE = GeoPreset(
    name="v3-reference", role="reference",
    parent_layout="aligned", child_layout="aligned",
    itot=96, jtot=96,
    child_itot=64, child_jtot=64,          # nominal; a reference has no child
    u0=3.0,
    spinup=3600.0, production=3600.0, child_spinup=600.0,
    nprocx=8, nprocy=8, child_nprocx=8, child_nprocy=8,
    parent_expnr="920", child_expnr="920",
    dtdump=3.0,  # no child is ever driven from this run; the old cadence suffices
    **_V3_COMMON,
)

#: **The parent.**  Flat: no buildings anywhere, the mesoscale-parent case.
#:
#: Forcing.  A flat surface at the ground roughness of the cases
#: (``factypes`` type 1, z0 = 0.05 m) driven by the canopy's own ``dpdx`` would
#: run at roughly twice the canopy's equilibrium bulk velocity, because the
#: canopy drag it is missing is most of the drag.  That is a *bulk momentum*
#: mismatch, and it is not the one section 9.4 is about: the child's interior
#: would then be decelerating everywhere and no near-surface equilibrium would
#: exist to adjust to.  So the parent is driven by a volume-flow-rate
#: controller at the reference run's measured bulk velocity instead
#: (``bulk_from="reference"``, which fills ``uflowrate`` in), leaving the child
#: in global momentum balance and the mismatch confined to the *shape* of the
#: near-surface profile -- a log layer down to the ground where the child has a
#: canopy.  That is exactly the imposed profile section 9.4 says will be out of
#: equilibrium, and the adjustment length is then a property of the internal
#: boundary layer rather than of a bulk imbalance.
#:
#: The consequence to keep in mind when reading the result: the flat parent
#: carries a *weaker* surface stress than the canopy does (u* about 0.18 m/s
#: against 0.4), so it also delivers less turbulence to the boundary than a
#: matched parent would.  Both are reported.
#:
#: The domain is 640 x 320 m = 5 x 2.5 lid depths, which the child's
#: 512 x 256 m window sits inside with one cell of margin for the staggered
#: face.  ``uflowrate`` is a placeholder here; ``run_geometry`` overwrites it
#: with the reference run's measurement and records which value was used.
V3_PARENT = GeoPreset(
    name="v3-parent", role="parent",
    parent_layout="none", child_layout="aligned", child_phase="standoff",
    itot=320, jtot=160,
    child_itot=256, child_jtot=128,
    # ``uflowrate`` is the parent's forcing; the children inherit the field by
    # ``replace`` but never use it -- ``make_child_case.child_sections`` always
    # writes ``luvolflowr = .false.`` and the fixed ``dpdx``, which is what puts
    # the child's interior in balance with its own canopy.
    u0=3.0, uflowrate=3.0,
    # 900 s discarded of a 3600 s window leaves 900 samples.  The child cold
    # starts from the parent's field, which already carries the right bulk and
    # the right large scales, so what has to develop is the canopy and its
    # internal boundary layer; 900 s is about 4.4 flow-throughs of the 512 m
    # domain at the equilibrium bulk velocity.  The child's own bulk trace is
    # emitted (``time_series`` in ``v3_metrics.json``) so the choice can be
    # checked against the data rather than trusted.
    spinup=9000.0, production=3600.0, child_spinup=900.0,
    nprocx=8, nprocy=8, child_nprocx=8, child_nprocy=8,
    parent_expnr="921", child_expnr="922",
    **_NESTDUMP_CADENCE, **_V3_COMMON,
)


def _v3_child(standoff: int, expnr: str, base: GeoPreset = V3_PARENT) -> GeoPreset:
    return replace(base, name=f"{base.name.rsplit('-', 1)[0]}-standoff{standoff}",
                   role="child", standoff_cells=standoff, child_expnr=expnr,
                   uflowrate=base.uflowrate)


#: Standoffs 0, 5 and 15 cells are the design table's own comparison
#: (section 10.4, row V3).  40 cells = 80 m = 5 h is added because the first
#: three are all short compared with any plausible adjustment scale -- 0, 0.6
#: and 1.9 building heights -- so a sweep of just those three could fail to
#: separate the hypotheses for want of lever rather than because the claim is
#: right.  A standoff that section 9.4's mechanism has room to act over is what
#: makes the claim falsifiable.
V3_STANDOFFS: Tuple[int, ...] = (0, 5, 15, 40)

#: **The cleared-parent-cubes arm** (nesting-plan-2026-09-06.md section 0, "New
#: from V2").  V2 found that a child which clears cubes out of its *own* zone
#: fails criterion A by 0.05-0.07 u* in the canopy layer, because the wakes
#: those removed cubes would have shed inside the zone are missing from the
#: inflow -- the mirror image of V3's own question, seen from the other side.
#: This arm isolates it: **the same child construction as** ``standoff0`` --
#: ``child_phase = "standoff"``, ``standoff_cells = 0`` -- so the child's own
#: lattice is placed from the clear box exactly as before and never depends on
#: the parent at all; only the periodic run it is driven from changes, from
#: flat (``V3_PARENT``) to V1's own aligned canopy filling the same domain
#: (``V3_PARENT_CUBES``).  ``_v3_child`` builds it unchanged -- no new
#: mechanism, just a different ``base``.  (``child_phase = "parent"`` was
#: tried first and rejected: it regenerates the *global* lattice and clips it
#: to the clear box, which is anchored to the parent's coordinate origin, not
#: to the clear box the way ``"standoff"`` is, so it places an entirely
#: different, disjoint set of cubes -- confirmed by comparing the two child
#: layouts directly, not assumed.  That would have compared two different
#: canopies as well as two different parents.)  Driven by the same fixed
#: ``dpdx``: a real canopy carries its own drag, unlike the flat parent, so no
#: volume-flow calibration against the reference is needed.
V3_PARENT_CUBES = GeoPreset(
    name="v3-parent-cubes", role="parent",
    parent_layout="aligned", child_layout="aligned", child_phase="standoff",
    itot=320, jtot=160,
    child_itot=256, child_jtot=128,
    u0=3.0,
    spinup=3600.0, production=3600.0, child_spinup=900.0,
    nprocx=8, nprocy=8, child_nprocx=8, child_nprocy=8,
    parent_expnr="926", child_expnr="927",
    **_NESTDUMP_CADENCE, **_V3_COMMON,
)

V3_CLEARED = replace(_v3_child(0, "927", base=V3_PARENT_CUBES),
                    name="v3-cleared-parent-cubes")

V3 = Experiment(
    name="v3", kind="v3",
    headline="parent without buildings: the adjustment length, and whether a "
             "standoff shortens or lengthens it",
    periodic=(
        PeriodicRun("reference", V3_REFERENCE,
                    "periodic 6x6 cube array at the child's dpdx -- the "
                    "equilibrium the child's canopy is measured against"),
        PeriodicRun("parent", V3_PARENT,
                    "flat, volume-flow forced at the reference's bulk velocity",
                    bulk_from="reference"),
        PeriodicRun("parent-cubes", V3_PARENT_CUBES,
                    "V1's own aligned canopy where the child's zone sits, "
                    "fixed dpdx -- the 'parent had cubes there' arm"),
    ),
    children=tuple(
        ChildRun(f"standoff{s}", _v3_child(s, str(922 + n)), "parent",
                 f"buildings start {s} cells past the clear box")
        for n, s in enumerate(V3_STANDOFFS)) + (
        ChildRun("cleared-parent-cubes", V3_CLEARED, "parent-cubes",
                 "standoff 0, but the parent had cubes in the child's zone and "
                 "the child clears them (nest_lparentgeom = .false.); compare "
                 "against 'standoff0', where the parent never had cubes there"),
    ),
    equilibrium_key="reference",
)


# --------------------------------------------------------------------------- #
# V4 -- different parent geometry
# --------------------------------------------------------------------------- #

#: V1's ``converged`` configuration, to the digit, except that the parent's
#: cube array is **staggered** rather than aligned.  Same grid, same cubes, same
#: plan area density, same forcing, same schedule, same zone, same child window,
#: same ranks.  That is the whole design: the only thing that moves is the
#: layout the parent resolves, so the V1 result is the reference and the answer
#: is a difference from it.
#: ``dtdump = 0.5`` s (``_NESTDUMP_CADENCE``) is the ``&NESTDUMP`` cadence the
#: child's boundary needs (``C_dump <= 2``, see ``_V3_COMMON``'s note); the
#: parent's own ``&OUTPUT`` field dump, read back by ``periodic-stats`` for its
#: canopy statistics, stays at the historical 3 s (``fielddump_dtdump``).
_V4_COMMON = dict(
    itot=256, jtot=256, ktot=64, dx=2.0,
    building_height=16.0, building_width=16.0, street_width=16.0, edgelength=16.0,
    geometry="uniform",
    child_itot=128, child_jtot=128,
    guardwidth=6.0, zonewidth=18.0, tau=1.0, nzone=12, nwall=1, timeinterp=2,
    ustar=0.4, u0=3.0, tke0=0.1,
    spinup=10800.0, production=10800.0, child_spinup=600.0,
    nprocx=8, nprocy=8, child_nprocx=8, child_nprocy=8, dtmax=0.5,
    spectra_heights=(8.0, 16.0, 32.0), stride=1,
    parent_layout="staggered", parent_expnr="930",
    **_NESTDUMP_CADENCE,
)

V4_PARENT = GeoPreset(name="v4-parent", role="parent",
                      child_layout="aligned", child_expnr="931", **_V4_COMMON)

#: The experiment: the child carries V1's **aligned** array while its parent
#: carries the staggered one.  ``child_phase = "parent"`` regenerates the global
#: aligned lattice, so this child's cubes are V1's child's cubes -- the tiny
#: test asserts that cube for cube against ``config.CONVERGED``.
V4_MISMATCH = GeoPreset(name="v4-mismatch", role="child",
                        child_layout="aligned", child_expnr="931", **_V4_COMMON)

#: **There is deliberately no matched-geometry control at production size, and
#: the reason is a geometric impossibility rather than a budget.**
#:
#: The obvious control is the same child carrying its parent's *staggered*
#: array -- V1 repeated on this realisation.  It cannot have a building-free
#: zone.  Write ``c`` for a cube centre in child metres and ``L`` for the child
#: extent: a cube is dropped from the clear box when ``c < 34`` or ``c > L-34``
#: (26 m of clearance plus the 8 m half-width), and it also reaches the analysis
#: interior when ``16 < c < L-16``.  So the residues mod the 32 m period that a
#: dropped cube may *not* occupy span two 18 m windows -- 36 m of a 32 m period
#: -- and the two column families of a staggered array sit exactly half a period
#: apart, so one of them always lands in a blocked window whatever the child's
#: origin or size.  Clearing the zone would therefore remove buildings from the
#: region the statistics are taken over, which
#: ``removed_cubes_reaching_the_interior`` refuses and should refuse.
#:
#: The alternatives are to keep the parent's cubes inside the zone and run
#: ``nest_lparentgeom = .true.`` -- legal for self-nesting, but then the control
#: no longer shares V1's boundary treatment and stops being a control -- or to
#: widen the zone, which changes the variable under test.  The tiny experiment
#: does carry a matched control (a 64-cell child in a 96-cell parent clears
#: cleanly), so the code path is exercised; at production size the comparison
#: that answers V4 is the mismatched child against V1's, and that needs no
#: control on this parent.
V4 = Experiment(
    name="v4", kind="v4",
    headline="different parent geometry: what a mismatched parent layout costs "
             "in the interior, as a difference from V1",
    periodic=(PeriodicRun("parent", V4_PARENT,
                          "V1's parent with a staggered array instead of an "
                          "aligned one; everything else identical"),),
    children=(
        ChildRun("mismatch", V4_MISMATCH, "parent",
                 "V1's child (aligned) driven by a staggered parent"),
    ),
    external_baseline=config.CONVERGED,
)


# --------------------------------------------------------------------------- #
# Tiny versions -- the identical code path in minutes
# --------------------------------------------------------------------------- #

#: ``dtdump`` is set per preset below (a "reference" run keeps 3 s; a "parent"
#: run gets the ``_NESTDUMP_CADENCE`` fine cadence), exactly as in the
#: production dicts above.
_TINY_BASE = dict(
    ktot=32, dx=2.0,
    building_height=16.0, building_width=16.0, street_width=16.0, edgelength=16.0,
    geometry="uniform",
    tau=1.0, nwall=1, timeinterp=2,
    ustar=0.4, u0=3.0, tke0=0.1, dtmax=0.5,
    spinup=40.0, production=90.0, child_spinup=30.0,
    nprocx=2, nprocy=2, child_nprocx=2, child_nprocy=2,
    spectra_heights=(8.0, 16.0), stride=1,
)

#: V3's tiny children build their canopy inside the clear box, so nothing is
#: ever dropped and a narrow zone costs nothing.
_TINY_V3 = dict(guardwidth=6.0, zonewidth=8.0, nzone=7, **_TINY_BASE)

#: V4's tiny children restrict the *parent's* lattice and then drop whatever
#: does not fit the clear box, so the zone has to be at least as deep as a cube
#: sits in from a child face -- 24 m in this 32 m array, since cube centres are
#: at 16 mod 32 and child faces at 0 mod 32.  The production zone (3 + 9 cells,
#: 26 m of clearance) clears it exactly and a tiny zone does not, which is the
#: same reason ``config.TINY_SWEEP`` carries the production zone.  Get this
#: wrong and ``removed_cubes_reaching_the_interior`` fires -- as it should.  All
#: five of V4-tiny's presets (both parents, both children, and the
#: never-run ``V4_TINY_MATCHED_IMPOSSIBLE``) are "parent"/"child" roles that
#: drive or are driven, so the fine ``&NESTDUMP`` cadence goes in here directly.
_TINY_V4 = dict(guardwidth=6.0, zonewidth=18.0, nzone=12,
               **_NESTDUMP_CADENCE, **_TINY_BASE)

V3_TINY_REFERENCE = GeoPreset(
    name="v3-tiny-reference", role="reference",
    parent_layout="aligned", child_layout="aligned",
    itot=48, jtot=48, child_itot=32, child_jtot=32,
    parent_expnr="940", child_expnr="940", dtdump=3.0, **_TINY_V3)

V3_TINY_PARENT = GeoPreset(
    name="v3-tiny-parent", role="parent",
    parent_layout="none", child_layout="aligned", child_phase="standoff",
    itot=96, jtot=64, child_itot=64, child_jtot=48,
    uflowrate=3.0,
    parent_expnr="941", child_expnr="942",
    **_NESTDUMP_CADENCE, **_TINY_V3)

#: The tiny cleared-parent-cubes arm, exactly as the production
#: ``V3_PARENT_CUBES``: same domain and zone as ``V3_TINY_PARENT`` (the
#: child's own ``"standoff"`` lattice depends on the clear box alone, not on
#: the parent), only ``parent_layout`` and the forcing differ.
V3_TINY_PARENT_CUBES = GeoPreset(
    name="v3-tiny-parent-cubes", role="parent",
    parent_layout="aligned", child_layout="aligned", child_phase="standoff",
    itot=96, jtot=64, child_itot=64, child_jtot=48,
    parent_expnr="945", child_expnr="946",
    **_NESTDUMP_CADENCE, **_TINY_V3)

V3_TINY_CLEARED = replace(_v3_child(0, "946", base=V3_TINY_PARENT_CUBES),
                         name="v3-tiny-cleared-parent-cubes")

V3_TINY = Experiment(
    name="v3-tiny", kind="v3",
    headline="V3 at a size that runs on a login node",
    periodic=(
        PeriodicRun("reference", V3_TINY_REFERENCE, "tiny periodic 3x3 array"),
        PeriodicRun("parent", V3_TINY_PARENT, "tiny flat parent",
                    bulk_from="reference"),
        PeriodicRun("parent-cubes", V3_TINY_PARENT_CUBES,
                    "tiny aligned canopy parent -- 'parent had cubes there'"),
    ),
    children=tuple(
        ChildRun(f"standoff{s}",
                 replace(V3_TINY_PARENT, name=f"v3-tiny-standoff{s}", role="child",
                         standoff_cells=s, child_expnr=str(942 + n)),
                 "parent", f"standoff {s} cells")
        for n, s in enumerate((0, 5, 15))) + (
        ChildRun("cleared-parent-cubes", V3_TINY_CLEARED, "parent-cubes",
                 "standoff 0, parent had cubes in the zone, child clears them"),
    ),
    equilibrium_key="reference",
)

_V4_TINY_COMMON = dict(
    itot=96, jtot=96, child_itot=64, child_jtot=64, **_TINY_V4)

#: The tiny V4 runs its own matched-geometry baseline -- an aligned parent
#: driving the same aligned child -- because the production baseline is the V1
#: ``converged`` child on disk and the smoke test must exercise the comparison
#: without it.  Like the production experiment it carries **no** staggered
#: control child: a staggered array admits no clean clear-box cut at any size
#: (see ``V4`` above), and ``V4_TINY_MATCHED_IMPOSSIBLE`` below turns that
#: argument into something the tiny test can fire.
V4_TINY = Experiment(
    name="v4-tiny", kind="v4",
    headline="V4 at a size that runs on a login node",
    periodic=(
        PeriodicRun("aligned", GeoPreset(
            name="v4-tiny-parent-aligned", role="parent", parent_layout="aligned",
            child_layout="aligned", parent_expnr="950", child_expnr="952",
            **_V4_TINY_COMMON), "the matched parent"),
        PeriodicRun("staggered", GeoPreset(
            name="v4-tiny-parent-staggered", role="parent",
            parent_layout="staggered", child_layout="aligned",
            parent_expnr="951", child_expnr="953",
            **_V4_TINY_COMMON), "the mismatched parent"),
    ),
    children=(
        ChildRun("baseline", GeoPreset(
            name="v4-tiny-baseline", role="child", parent_layout="aligned",
            child_layout="aligned", parent_expnr="950", child_expnr="952",
            **_V4_TINY_COMMON), "aligned",
            "aligned child on an aligned parent -- the V1 arrangement"),
        ChildRun("mismatch", GeoPreset(
            name="v4-tiny-mismatch", role="child", parent_layout="staggered",
            child_layout="aligned", parent_expnr="951", child_expnr="953",
            **_V4_TINY_COMMON), "staggered",
            "the same child on a staggered parent"),
    ),
    baseline_child_key="baseline",
)

#: A staggered control child, defined **only** so that the tiny test can assert
#: that ``validate`` refuses it.  See the note on ``V4`` for why no window and
#: no size admits one: the residues a dropped cube may not occupy cover more
#: than a full period, and a staggered array's two column families are exactly
#: half a period apart.
V4_TINY_MATCHED_IMPOSSIBLE = GeoPreset(
    name="v4-tiny-matched", role="child", parent_layout="staggered",
    child_layout="staggered", parent_expnr="951", child_expnr="954",
    **_V4_TINY_COMMON)


# --------------------------------------------------------------------------- #
# Registration
# --------------------------------------------------------------------------- #

EXPERIMENTS: Dict[str, Experiment] = {e.name: e for e in (V3, V4, V3_TINY, V4_TINY)}

#: Every preset this module defines, by name.
GEO_PRESETS: Dict[str, GeoPreset] = {}
for _exp in EXPERIMENTS.values():
    for _r in _exp.periodic:
        GEO_PRESETS[_r.preset.name] = _r.preset
    for _c in _exp.children:
        GEO_PRESETS[_c.preset.name] = _c.preset

# Register with config so that `config.get_preset(name)` -- which analyse.py's
# and make_child_case.py's command-line entry points use -- can find them.  This
# is an import-time side effect on purpose: it is the "registers itself" half of
# keeping V3/V4 out of config.py while still sharing its machinery.
config.PRESETS.update(GEO_PRESETS)


def get_experiment(name: str) -> Experiment:
    try:
        exp = EXPERIMENTS[name]
    except KeyError:
        raise SystemExit(
            f"unknown experiment {name!r}; choose one of "
            f"{', '.join(sorted(EXPERIMENTS))}") from None
    exp.validate()
    return exp


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="print the V3 and V4 presets")
    ap.add_argument("--experiment", default=None)
    ap.add_argument("--full", action="store_true",
                    help="print every preset in full, not just the table")
    ns = ap.parse_args()
    names = [ns.experiment] if ns.experiment else sorted(EXPERIMENTS)
    for name in names:
        exp = get_experiment(name)
        print(exp.summary())
        print()
        if ns.full:
            for r in exp.periodic:
                print(f"--- periodic {r.key}: {r.note}")
                print(r.preset.summary())
                print()
            for c in exp.children:
                print(f"--- child {c.key}: {c.note}")
                print(c.preset.summary())
                print()
