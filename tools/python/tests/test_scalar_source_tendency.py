"""Tests for the scalar-source-tendency kernel formulas.

Replicates the point-source and line-source Gaussian kernels from
initialize_scalar_source_tendency (src/modscalsource.f90) in Python and
verifies mathematical properties: peak location, Gaussian decay, 3-sigma
cutoff, symmetry, and linear scaling with source strength.

These tests run without a compiled solver and are auto-discovered by the
python-library CI suite.
"""

from __future__ import annotations

import math
import unittest

import numpy as np


def _point_source_tendency(
    nx: int,
    ny: int,
    nz: int,
    dx: float,
    dy: float,
    dz: float,
    x_src: float,
    y_src: float,
    z_src: float,
    strength: float,
    sigma: float,
) -> np.ndarray:
    """Gaussian point-source kernel matching initialize_scalar_source_tendency.

    Cell centres follow the single-rank Fortran convention:
      px = (i_0based + 0.5) * dx  (equivalent to (ig - 0.5)*dx with ig=i, ib=1, zstart=1)
    """
    tendency = np.zeros((nx, ny, nz), dtype=np.float64)
    sigma2 = sigma ** 2
    cutoff2 = 9.0 * sigma2
    dxi = 1.0 / dx
    dyi = 1.0 / dy

    for iz in range(nz):
        pz = (iz + 0.5) * dz
        dzfi = 1.0 / dz
        for iy in range(ny):
            py = (iy + 0.5) * dy
            for ix in range(nx):
                px = (ix + 0.5) * dx
                ra2 = (px - x_src) ** 2 + (py - y_src) ** 2 + (pz - z_src) ** 2
                if ra2 <= cutoff2:
                    tendency[ix, iy, iz] += dxi * dyi * dzfi * strength * math.exp(-ra2 / (2.0 * sigma2))
    return tendency


def _line_source_tendency(
    nx: int,
    ny: int,
    nz: int,
    dx: float,
    dy: float,
    dz: float,
    x0: float,
    y0: float,
    z0: float,
    x1: float,
    y1: float,
    z1: float,
    strength: float,
    sigma: float,
) -> np.ndarray:
    """Gaussian line-source kernel matching initialize_scalar_source_tendency."""
    tendency = np.zeros((nx, ny, nz), dtype=np.float64)
    sigma2 = sigma ** 2
    cutoff2 = 9.0 * sigma2
    sqrt_2pi = math.sqrt(2.0 * math.pi)
    lsx, lsy, lsz = x1 - x0, y1 - y0, z1 - z0
    line_length2 = lsx ** 2 + lsy ** 2 + lsz ** 2
    dxi = 1.0 / dx
    dyi = 1.0 / dy

    for iz in range(nz):
        pz = (iz + 0.5) * dz
        dzfi = 1.0 / dz
        for iy in range(ny):
            py = (iy + 0.5) * dy
            for ix in range(nx):
                px = (ix + 0.5) * dx
                vx, vy, vz = px - x0, py - y0, pz - z0
                t = (vx * lsx + vy * lsy + vz * lsz) / line_length2
                if t < 0.0:
                    ra2 = vx ** 2 + vy ** 2 + vz ** 2
                elif t > 1.0:
                    ra2 = (px - x1) ** 2 + (py - y1) ** 2 + (pz - z1) ** 2
                else:
                    ra2 = (px - (x0 + t * lsx)) ** 2 + \
                          (py - (y0 + t * lsy)) ** 2 + \
                          (pz - (z0 + t * lsz)) ** 2
                if ra2 <= cutoff2:
                    tendency[ix, iy, iz] += (
                        dxi * dyi * dzfi * sqrt_2pi * strength * sigma *
                        math.exp(-ra2 / (2.0 * sigma2)) *
                        math.erf(math.sqrt((cutoff2 - ra2) / (2.0 * sigma2)))
                    )
    return tendency


class TestPointSourceTendency(unittest.TestCase):
    """Mathematical-property tests for the Gaussian point-source kernel."""

    def setUp(self):
        # 8×8×8 uniform grid matching the GPU smoke-test case (thermo-scalar-wall)
        self.nx = self.ny = self.nz = 8
        self.dx = self.dy = self.dz = 1.0
        self.x_src, self.y_src, self.z_src = 4.0, 4.0, 1.5
        self.strength = 1.0
        self.sigma = 1.0
        self.t = _point_source_tendency(
            self.nx, self.ny, self.nz,
            self.dx, self.dy, self.dz,
            self.x_src, self.y_src, self.z_src,
            self.strength, self.sigma,
        )

    def test_peak_z_at_source(self):
        # z_src=1.5 is the centre of cell iz=1 (dz=1, centre at 1.5)
        peak = np.unravel_index(np.argmax(self.t), self.t.shape)
        self.assertEqual(peak[2], 1)

    def test_peak_xy_adjacent_to_source(self):
        # x_src=4.0 is equidistant from cell centres at 3.5 (ix=3) and 4.5 (ix=4)
        peak = np.unravel_index(np.argmax(self.t), self.t.shape)
        self.assertIn(peak[0], (3, 4))
        self.assertIn(peak[1], (3, 4))

    def test_peak_value_matches_gaussian(self):
        # Cell (3,3,1): centre (3.5, 3.5, 1.5), r²=(0.5²+0.5²+0²)=0.5
        expected = math.exp(-0.5 / (2.0 * self.sigma ** 2)) / (self.dx * self.dy * self.dz)
        self.assertAlmostEqual(float(self.t[3, 3, 1]), expected, places=12)
        self.assertGreater(expected, 0.0)

    def test_cutoff_applied_beyond_three_sigma(self):
        for iz in range(self.nz):
            pz = (iz + 0.5) * self.dz
            for iy in range(self.ny):
                py = (iy + 0.5) * self.dy
                for ix in range(self.nx):
                    px = (ix + 0.5) * self.dx
                    r = math.sqrt((px - self.x_src) ** 2 + (py - self.y_src) ** 2 + (pz - self.z_src) ** 2)
                    if r > 3.0 * self.sigma + 1e-10:
                        self.assertEqual(
                            self.t[ix, iy, iz], 0.0,
                            msg=f"Expected zero at ({ix},{iy},{iz}): r={r:.3f} > 3σ",
                        )

    def test_decay_with_distance(self):
        # Along the z-axis from the source: values must decrease as r grows
        ix, iy = 3, 3  # nearest xy cell
        t_z1 = float(self.t[ix, iy, 1])   # iz=1, pz=1.5 = z_src
        t_z2 = float(self.t[ix, iy, 2])   # iz=2, pz=2.5, r_z=1.0
        t_z3 = float(self.t[ix, iy, 3])   # iz=3, pz=3.5, r_z=2.0
        self.assertGreater(t_z1, t_z2)
        self.assertGreater(t_z2, t_z3)

    def test_xy_equidistant_cells_equal(self):
        # (3,3) and (4,4) both at r_xy=sqrt(0.5²+0.5²) from (4.0,4.0); must be equal
        iz = 1
        self.assertAlmostEqual(float(self.t[3, 3, iz]), float(self.t[4, 4, iz]), places=12)
        self.assertAlmostEqual(float(self.t[3, 4, iz]), float(self.t[4, 3, iz]), places=12)

    def test_all_non_negative(self):
        self.assertTrue(np.all(self.t >= 0.0))

    def test_strength_scales_linearly(self):
        t2 = _point_source_tendency(
            self.nx, self.ny, self.nz,
            self.dx, self.dy, self.dz,
            self.x_src, self.y_src, self.z_src,
            strength=2.0, sigma=self.sigma,
        )
        np.testing.assert_allclose(t2, 2.0 * self.t)

    def test_minimum_sigma_floor(self):
        # sigma below dx should still produce a non-zero tendency (no divide-by-zero)
        t_small = _point_source_tendency(
            self.nx, self.ny, self.nz,
            self.dx, self.dy, self.dz,
            self.x_src, self.y_src, self.z_src,
            strength=1.0, sigma=1.0,  # sigma == dx: the Fortran-validated minimum
        )
        self.assertGreater(float(np.max(t_small)), 0.0)

    def test_source_outside_domain_still_affects_border(self):
        # Source just outside the domain edge should contribute to edge cells via
        # the Gaussian tail if within 3σ
        t_edge = _point_source_tendency(
            4, 4, 4, 1.0, 1.0, 1.0,
            x_src=-0.1, y_src=2.0, z_src=2.0,
            strength=1.0, sigma=1.0,
        )
        # Cell (0,1,1): centre (0.5, 1.5, 1.5), r≈0.6 < 3σ → nonzero
        self.assertGreater(float(t_edge[0, 1, 1]), 0.0)


class TestLineSourceTendency(unittest.TestCase):
    """Mathematical-property tests for the Gaussian line-source kernel."""

    def setUp(self):
        self.nx = self.ny = self.nz = 8
        self.dx = self.dy = self.dz = 1.0
        # Vertical line at x=4, y=4 from z=2 to z=6
        self.t = _line_source_tendency(
            self.nx, self.ny, self.nz,
            self.dx, self.dy, self.dz,
            4.0, 4.0, 2.0, 4.0, 4.0, 6.0,
            strength=1.0, sigma=1.0,
        )

    def test_peak_xy_adjacent_to_line(self):
        peak = np.unravel_index(np.argmax(self.t), self.t.shape)
        self.assertIn(peak[0], (3, 4))
        self.assertIn(peak[1], (3, 4))

    def test_peak_z_within_line_span(self):
        # Line spans iz=2 (pz=2.5) to iz=5 (pz=5.5); peak must lie within
        peak = np.unravel_index(np.argmax(self.t), self.t.shape)
        self.assertIn(peak[2], (1, 2, 3, 4, 5))

    def test_all_non_negative(self):
        self.assertTrue(np.all(self.t >= 0.0))

    def test_values_symmetric_around_line(self):
        # Cells equidistant from the line axis must have equal tendencies
        # For the vertical line at x=4,y=4: cells (3,3,iz) and (4,4,iz) are equidistant
        for iz in range(self.nz):
            self.assertAlmostEqual(
                float(self.t[3, 3, iz]),
                float(self.t[4, 4, iz]),
                places=12,
            )

    def test_cutoff_applied_far_from_line(self):
        # Cells with perpendicular distance > 3σ from the line must be zero
        for iz in range(self.nz):
            pz = (iz + 0.5) * self.dz
            # Clamp projection to line segment [2.0, 6.0]
            t_proj = max(2.0, min(6.0, pz))
            for iy in range(self.ny):
                py = (iy + 0.5) * self.dy
                for ix in range(self.nx):
                    px = (ix + 0.5) * self.dx
                    r_perp = math.sqrt((px - 4.0) ** 2 + (py - 4.0) ** 2 + (pz - t_proj) ** 2)
                    if r_perp > 3.0 + 1e-10:
                        self.assertEqual(
                            self.t[ix, iy, iz], 0.0,
                            msg=f"Expected zero at ({ix},{iy},{iz}): r_perp={r_perp:.3f}",
                        )

    def test_strength_scales_linearly(self):
        t2 = _line_source_tendency(
            self.nx, self.ny, self.nz,
            self.dx, self.dy, self.dz,
            4.0, 4.0, 2.0, 4.0, 4.0, 6.0,
            strength=2.0, sigma=1.0,
        )
        np.testing.assert_allclose(t2, 2.0 * self.t)


if __name__ == "__main__":
    unittest.main()
