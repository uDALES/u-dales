#!/usr/bin/env python3
"""``indoor_geometry`` must reproduce the published GMD 2024 enclosure.

The reference is the case's own ``indoor_object_final.stl`` (experiment 567,
doi:10.5281/zenodo.12510825).  If that file is reachable -- via
``UDALES_GMD2024_INPUTS`` pointing at the extracted ``indoor-outdoor``
directory -- the structural planes are compared against it directly.  The
plane values are also hard-coded, so the check still means something without
the download.

No solver, no MPI.
"""

from __future__ import annotations

import os
import struct
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import numpy as np

import indoor_geometry as ig

#: Structural planes of the published STL, read off it once (see
#: test_matches_the_published_stl for the live comparison).
PAPER_X = [0.480, 0.500, 0.700, 0.720]
PAPER_Y = [0.780, 0.800, 0.854, 0.946, 1.000, 1.020]
PAPER_Z = [0.000, 0.062, 0.098, 0.160, 0.180]

#: Coplanar-face midpoints the published mesh adds. Tessellation, not geometry.
PAPER_MIDPOINTS = {"x": [0.490, 0.600, 0.710], "y": [0.900], "z": [0.080, 0.090]}


def _published_stl() -> Path | None:
    root = os.environ.get("UDALES_GMD2024_INPUTS")
    if not root:
        return None
    p = Path(root) / "567" / "indoor_object_final.stl"
    return p if p.is_file() else None


def _stl_vertices(path: Path) -> np.ndarray:
    with path.open("rb") as f:
        f.read(80)
        (n,) = struct.unpack("<I", f.read(4))
        buf = bytearray(f.read(n * 50))
    tri = np.frombuffer(buf, count=n, dtype=np.dtype(
        [("normal", "<3f4"), ("v", "<3, 3f4"), ("attr", "<u2")]))
    return tri["v"]


def _overlaps_positively(a, b) -> bool:
    """Whether 1-D intervals ``a`` and ``b`` overlap by more than a point."""
    return min(a[1], b[1]) - max(a[0], b[0]) > 1e-12


def _box_intersects_region(box, region) -> bool:
    """Whether an axis-aligned box overlaps a region on all three axes at once
    -- i.e. with positive volume, not just a touching face."""
    x0, x1, y0, y1, z0, z1 = box
    xr, yr, zr = region
    return (_overlaps_positively((x0, x1), xr)
            and _overlaps_positively((y0, y1), yr)
            and _overlaps_positively((z0, z1), zr))


def _box_within(box, region) -> bool:
    """Whether an axis-aligned box lies inside a region (touching allowed)."""
    x0, x1, y0, y1, z0, z1 = box
    (rx0, rx1), (ry0, ry1), (rz0, rz1) = region
    eps = 1e-9
    return (rx0 - eps <= x0 and x1 <= rx1 + eps
            and ry0 - eps <= y0 and y1 <= ry1 + eps
            and rz0 - eps <= z0 and z1 <= rz1 + eps)


class TestTheInteriorIsInvariant(unittest.TestCase):
    """The cavity is the physical quantity and must not move with ``t``."""

    #: The paper's "0.2 x 0.2 x 0.16 m enclosure", centred as in the domain.
    EXPECTED_INTERIOR = ((0.5, 0.7), (0.8, 1.0), (0.0, 0.16))

    def _assert_bounds_almost_equal(self, got, want):
        for (g0, g1), (w0, w1) in zip(got, want):
            self.assertAlmostEqual(g0, w0, places=12)
            self.assertAlmostEqual(g1, w1, places=12)

    def test_volume_is_independent_of_wall_thickness(self):
        self.assertAlmostEqual(ig.interior_volume(), 0.2 * 0.2 * 0.16, places=12)
        interior = ig.interior_bounds()
        self._assert_bounds_almost_equal(interior, self.EXPECTED_INTERIOR)
        for t in (0.002, 0.007, 0.020, 0.040):
            with self.subTest(t=t):
                self.assertEqual(ig.interior_bounds(), interior,
                                  "the cavity must not move with the wall")
                self.assertAlmostEqual(ig.interior_volume(), 0.0064, places=12)

                # No wall/roof box may intrude into the open cavity, at any t.
                for box in ig.boxes(t):
                    self.assertFalse(
                        _box_intersects_region(box, interior),
                        f"a shell box intrudes into the cavity at t={t}")

                # Every box stays inside the outer envelope, which is inset
                # from the interior by exactly t on each side (and on top;
                # the floor is the ground, shared by interior and outer).
                outer = ig.outer_bounds(t)
                for box in ig.boxes(t):
                    self.assertTrue(
                        _box_within(box, outer),
                        f"a box escapes the outer envelope at t={t}")
                (xi0, xi1), (yi0, yi1), (_, zi1) = interior
                (xo0, xo1), (yo0, yo1), (_, zo1) = outer
                self.assertAlmostEqual(xi0 - xo0, t, places=12)
                self.assertAlmostEqual(xo1 - xi1, t, places=12)
                self.assertAlmostEqual(yi0 - yo0, t, places=12)
                self.assertAlmostEqual(yo1 - yi1, t, places=12)
                self.assertAlmostEqual(zo1 - zi1, t, places=12,
                                       msg="the floor is the ground: one wall of height")

    def test_the_outer_box_grows_with_the_wall(self):
        for t in (0.007, 0.020):
            (xo0, xo1), (yo0, yo1), (_, zo1) = ig.outer_bounds(t)
            (xi0, xi1), (yi0, yi1), (_, zi1) = ig.interior_bounds()
            self.assertAlmostEqual(xi0 - xo0, t, places=12)
            self.assertAlmostEqual(xo1 - xi1, t, places=12)
            self.assertAlmostEqual(yi0 - yo0, t, places=12)
            self.assertAlmostEqual(zo1 - zi1, t, places=12,
                                   msg="the floor is the ground: one wall of height")

    def test_the_window_is_a_duct_whose_depth_is_the_wall(self):
        for t in (0.007, 0.020):
            for face in ("windward", "leeward"):
                (x0, x1), (y0, y1), (z0, z1) = ig.window_bounds(t, face)
                self.assertAlmostEqual(x1 - x0, t, places=12, msg=f"{face} duct depth")
                self.assertAlmostEqual(y1 - y0, ig.WINDOW_WIDTH, places=12)
                self.assertAlmostEqual(z1 - z0, ig.WINDOW_HEIGHT, places=12)
                self.assertAlmostEqual(z0, ig.WINDOW_SILL, places=12)
                # 0.054 m in from each interior wall, per the paper
                (_, _), (yi0, yi1), _ = ig.interior_bounds()
                self.assertAlmostEqual(y0 - yi0, 0.054, places=12)
                self.assertAlmostEqual(yi1 - y1, 0.054, places=12)

    def test_an_unusable_wall_is_refused(self):
        for bad in (0.0, -0.01):
            with self.assertRaises(ValueError):
                ig.boxes(bad)


class TestReproducesThePaper(unittest.TestCase):

    def test_structural_planes_at_the_published_thickness(self):
        px, py, pz = ig.structural_planes(ig.PAPER_WALL)
        for got, want, ax in ((px, PAPER_X, "x"), (py, PAPER_Y, "y"), (pz, PAPER_Z, "z")):
            self.assertEqual([round(v, 6) for v in got], want, f"{ax} planes")

    def test_the_shell_does_not_overlap_itself(self):
        """solid_volume sums the boxes, so it is only right if they are disjoint."""
        t = ig.PAPER_WALL
        (xo0, xo1), (yo0, yo1), (_, zo1) = ig.outer_bounds(t)
        outer = (xo1 - xo0) * (yo1 - yo0) * zo1
        cavity = ig.interior_volume()
        self.assertAlmostEqual(ig.solid_volume(t) + cavity + self._window_volume(t),
                               outer, places=9)

    @staticmethod
    def _window_volume(t: float) -> float:
        return 2 * t * ig.WINDOW_WIDTH * ig.WINDOW_HEIGHT

    @unittest.skipIf(_published_stl() is None,
                     "set UDALES_GMD2024_INPUTS to the extracted indoor-outdoor directory")
    def test_matches_the_published_stl(self):
        v = _stl_vertices(_published_stl())
        house = v[v[:, :, 2].max(axis=1) > 1e-4].reshape(-1, 3)   # drop the ground
        px, py, pz = ig.structural_planes(ig.PAPER_WALL)
        for i, (ax, mine) in enumerate(zip("xyz", (px, py, pz))):
            theirs = sorted({round(float(c), 4) for c in house[:, i]})
            extra = [p for p in theirs
                     if not any(abs(p - m) < 1e-4 for m in mine)]
            self.assertEqual([round(p, 3) for p in extra], PAPER_MIDPOINTS[ax],
                             f"{ax}: unexpected planes in the published mesh")
            missing = [m for m in mine
                       if not any(abs(p - m) < 1e-4 for p in theirs)]
            self.assertEqual(missing, [], f"{ax}: planes the published mesh has that we lack")

    @unittest.skipIf(_published_stl() is None, "needs the published inputs")
    def test_every_extra_plane_is_a_midpoint(self):
        """Justifies ignoring them: each bisects a pair we do reproduce."""
        px, py, pz = ig.structural_planes(ig.PAPER_WALL)
        for ax, mine in zip("xyz", (px, py, pz)):
            for p in PAPER_MIDPOINTS[ax]:
                ok = any(abs((a + b) / 2 - p) < 1e-9
                         for a in mine for b in mine if b > a)
                self.assertTrue(ok, f"{ax}={p} is not the midpoint of any pair we have")


class TestTheThinWalledVariant(unittest.TestCase):
    """What the refined child changes, and what it must not."""

    THIN = 0.007

    def test_the_cavity_and_openings_are_unchanged(self):
        interior = ig.interior_bounds()
        for t in (0.002, 0.007, 0.020, 0.040):
            with self.subTest(t=t):
                self.assertEqual(ig.interior_bounds(), interior,
                                  "the cavity must not move with the wall")
                for face in ("windward", "leeward"):
                    _, y, z = ig.window_bounds(t, face)
                    _, paper_y, paper_z = ig.window_bounds(ig.PAPER_WALL, face)
                    self.assertEqual(y, paper_y, "opening moved in y")
                    self.assertEqual(z, paper_z, "opening moved in z")

    def test_only_the_duct_depth_and_the_outer_box_change(self):
        self.assertAlmostEqual(
            ig.window_bounds(ig.PAPER_WALL, "windward")[0][1]
            - ig.window_bounds(ig.PAPER_WALL, "windward")[0][0], 0.020, places=12)
        (x0, x1), _, _ = ig.window_bounds(self.THIN, "windward")
        self.assertAlmostEqual(x1 - x0, self.THIN, places=12)
        self.assertLess(ig.solid_volume(self.THIN), ig.solid_volume(ig.PAPER_WALL))


class TestWallAdmissibility(unittest.TestCase):
    """min_wall / check_wall / WallTooThin: the paper's own admissibility rule."""

    def test_min_wall_uses_the_coarser_spacing(self):
        self.assertAlmostEqual(ig.min_wall(0.01, 0.02), 3 * 0.02, places=12)
        self.assertAlmostEqual(ig.min_wall(0.02, 0.01), 3 * 0.02, places=12)

    def test_min_wall_honours_a_different_points(self):
        self.assertAlmostEqual(ig.min_wall(0.01, 0.02, points=5), 5 * 0.02, places=12)

    def test_check_wall_returns_the_thickness_when_admissible(self):
        dx = dy = 0.01
        got = ig.check_wall(0.05, dx, dy)
        self.assertEqual(got, 0.05)
        self.assertIsInstance(got, float)

    def test_check_wall_passes_at_the_exact_boundary(self):
        dx = dy = 0.01
        t = 3 * dx
        self.assertAlmostEqual(ig.check_wall(t, dx, dy), t, places=12)

    def test_check_wall_raises_when_too_thin(self):
        self.assertTrue(issubclass(ig.WallTooThin, ValueError))
        dx = dy = 0.004
        with self.assertRaises(ig.WallTooThin) as ctx:
            ig.check_wall(0.005, dx, dy)
        msg = str(ctx.exception)
        self.assertIn("5.00 mm", msg, "should report the offending wall thickness in mm")
        self.assertIn("12.00 mm", msg, "should report the required minimum in mm")

    def test_the_paper_wall_passes_but_the_thin_walled_variant_is_refused(self):
        """The 3-cell rule stands, and it is what picks the thin-walled variant's wall.

        On the paper's own grid (dx = dy = 0.01 / 3 m, i.e. 3.333 mm) the
        published 0.020 m wall spans 6 cells and passes. The refined child can
        carry a thinner wall than the parent -- that is the point of running
        it -- but the thin-walled variant must still be at least 3 child
        cells: 3 * 3.333 mm = 10 mm. So 10 mm is the floor on this grid, and
        the 7 mm figure quoted in commit 511ae566's message (2.1 cells) is
        refused; a 7 mm wall would need a finer child still.
        """
        dx = dy = 0.01 / 3
        self.assertEqual(ig.check_wall(ig.PAPER_WALL, dx, dy), ig.PAPER_WALL)
        self.assertAlmostEqual(ig.check_wall(0.010, dx, dy), 0.010, places=12)
        with self.assertRaises(ig.WallTooThin):
            ig.check_wall(0.007, dx, dy)


class TestStlOutput(unittest.TestCase):

    def test_it_writes_a_readable_ascii_stl(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            p = ig.write_stl(Path(tmp) / "shell.stl", ig.PAPER_WALL)
            text = p.read_text(encoding="ascii")
            self.assertTrue(text.startswith("solid "))
            self.assertTrue(text.rstrip().endswith("endsolid enclosure"))
            nfacet = text.count("facet normal")
            self.assertEqual(nfacet, 12 * len(ig.boxes(ig.PAPER_WALL)))

    def test_the_ground_plane_is_optional(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            bare = ig.write_stl(Path(tmp) / "a.stl", ig.PAPER_WALL, ground=False)
            withg = ig.write_stl(Path(tmp) / "b.stl", ig.PAPER_WALL, ground=True)
            self.assertEqual(withg.read_text().count("facet normal"),
                             bare.read_text().count("facet normal") + 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
