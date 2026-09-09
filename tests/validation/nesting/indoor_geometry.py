#!/usr/bin/env python3
"""The GMD 2024 indoor-outdoor enclosure, parametrised on wall thickness.

The case is section 4.2 of *A conservative immersed boundary method for the
multi-physics urban large-eddy simulation model uDALES v2.0* (GMD 17, 6277,
2024), experiment 567, inputs at doi:10.5281/zenodo.12510825.

**The invariant is the interior.**  The cavity is the physical quantity -- it
sets the ventilation dynamics -- and the walls were thickened *outward* to
whatever the grid could represent.  So:

    outer footprint = interior + 2 t
    outer height    = interior + t          (the floor is the ground)

and the window is a *duct through the wall*, so its depth is ``t`` while its
opening and sill are fixed relative to the interior.  A finer child can
therefore use a thinner, more realistic wall and keep exactly the same cavity,
which is what makes the coarse-parent / fine-child comparison fair: same
volume, same openings, and one difference -- the duct the jet passes through.

At ``t = 0.02`` this reproduces every structural plane of the published STL
(``indoor_object_final.stl``); see ``test_indoor_geometry.py``.  The published
mesh additionally splits coplanar faces at their midpoints (0.49, 0.60, 0.71 in
x; 0.90 in y; 0.08, 0.09 in z), which is tessellation and carries no geometry.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Sequence, Tuple

#: Interior cavity [m].  Invariant: the paper's "0.2 x 0.2 x 0.16 m enclosure".
INTERIOR_X = 0.200
INTERIOR_Y = 0.200
INTERIOR_Z = 0.160

#: Window opening [m], and its sill above the floor.  Positioned relative to
#: the INTERIOR: 0.054 m in from each interior side wall, 0.062 m up.
WINDOW_WIDTH = 0.092
WINDOW_HEIGHT = 0.036
WINDOW_SILL = 0.062

#: Interior centre in the published 3.42 x 1.80 m domain: 3 H_out from the
#: inlet with H_out = 0.18 m.
CENTRE_X = 0.600
CENTRE_Y = 0.900

#: The published wall thickness, forced by the 3.34 mm grid.
PAPER_WALL = 0.020

Box = Tuple[float, float, float, float, float, float]


def interior_bounds() -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
    """``((x0, x1), (y0, y1), (z0, z1))`` of the cavity -- independent of ``t``."""
    return ((CENTRE_X - INTERIOR_X / 2, CENTRE_X + INTERIOR_X / 2),
            (CENTRE_Y - INTERIOR_Y / 2, CENTRE_Y + INTERIOR_Y / 2),
            (0.0, INTERIOR_Z))


def interior_volume() -> float:
    """Cavity volume [m3].  Constant by construction, whatever ``t`` is."""
    return INTERIOR_X * INTERIOR_Y * INTERIOR_Z


def outer_bounds(t: float) -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
    (xi0, xi1), (yi0, yi1), (_, zi1) = interior_bounds()
    return ((xi0 - t, xi1 + t), (yi0 - t, yi1 + t), (0.0, zi1 + t))


def window_bounds(t: float, face: str) -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
    """The opening as a box.  ``face`` is 'windward' (-x) or 'leeward' (+x)."""
    (xi0, xi1), _, _ = interior_bounds()
    if face == "windward":
        xr = (xi0 - t, xi0)
    elif face == "leeward":
        xr = (xi1, xi1 + t)
    else:
        raise ValueError(f"face must be 'windward' or 'leeward', not {face!r}")
    return (xr,
            (CENTRE_Y - WINDOW_WIDTH / 2, CENTRE_Y + WINDOW_WIDTH / 2),
            (WINDOW_SILL, WINDOW_SILL + WINDOW_HEIGHT))


def boxes(t: float) -> List[Box]:
    """The shell as axis-aligned boxes: 4 walls (2 pierced), plus the roof.

    A pierced wall is decomposed into the four rectangles around its opening,
    so every piece stays an axis-aligned box and the result is watertight
    without triangulating a face with a hole in it.
    """
    if t <= 0.0:
        raise ValueError(f"wall thickness must be positive, got {t}")
    (xi0, xi1), (yi0, yi1), (_, zi1) = interior_bounds()
    (xo0, xo1), (yo0, yo1), (_, zo1) = outer_bounds(t)
    wy0, wy1 = CENTRE_Y - WINDOW_WIDTH / 2, CENTRE_Y + WINDOW_WIDTH / 2
    wz0, wz1 = WINDOW_SILL, WINDOW_SILL + WINDOW_HEIGHT

    out: List[Box] = []
    # The pierced walls span only the INTERIOR width; the solid side walls take
    # the full outer x range and so own the four corner columns.  Splitting it
    # this way keeps the boxes disjoint, which the shape does not care about but
    # any volume or mass accounting does -- overlapping them double-counts
    # 4 t^2 h of corner.
    for x0, x1 in ((xo0, xi0), (xi1, xo1)):
        out += [(x0, x1, yi0, wy0, 0.0, zi1),
                (x0, x1, wy1, yi1, 0.0, zi1),
                (x0, x1, wy0, wy1, 0.0, wz0),
                (x0, x1, wy0, wy1, wz1, zi1)]
    out += [(xo0, xo1, yo0, yi0, 0.0, zi1),            # solid side walls
            (xo0, xo1, yi1, yo1, 0.0, zi1),            # (including the corners)
            (xo0, xo1, yo0, yo1, zi1, zo1)]            # roof
    return out


def structural_planes(t: float) -> Tuple[List[float], List[float], List[float]]:
    """Sorted distinct box-face coordinates per axis."""
    b = boxes(t)
    return (sorted({round(v, 9) for x in b for v in x[0:2]}),
            sorted({round(v, 9) for x in b for v in x[2:4]}),
            sorted({round(v, 9) for x in b for v in x[4:6]}))


def solid_volume(t: float) -> float:
    """Total solid volume [m3] -- the boxes do not overlap, so this is a sum."""
    return sum((x1 - x0) * (y1 - y0) * (z1 - z0) for x0, x1, y0, y1, z0, z1 in boxes(t))


_QUADS = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
          (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))


def write_stl(path: Path, t: float, ground: bool = False,
              domain: Sequence[float] = (3.42, 1.80)) -> Path:
    """Write the shell (optionally with a ground plane) as an ascii STL."""
    path = Path(path)
    tris: List[Tuple[Tuple[float, float, float], ...]] = []
    for (x0, x1, y0, y1, z0, z1) in boxes(t):
        c = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
             (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
        for a, b, d, e in _QUADS:
            tris += [(c[a], c[b], c[d]), (c[a], c[d], c[e])]
    if ground:
        gx, gy = float(domain[0]), float(domain[1])
        g = [(0.0, 0.0, 0.0), (gx, 0.0, 0.0), (gx, gy, 0.0), (0.0, gy, 0.0)]
        tris += [(g[0], g[1], g[2]), (g[0], g[2], g[3])]
    with path.open("w", encoding="ascii") as f:
        f.write("solid enclosure\n")
        for tri in tris:
            f.write("facet normal 0 0 0\n  outer loop\n")
            for v in tri:
                f.write("    vertex %.6f %.6f %.6f\n" % v)
            f.write("  endloop\nendfacet\n")
        f.write("endsolid enclosure\n")
    return path
