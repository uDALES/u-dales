"""Characterization tests for pure numeric helpers in the radiation section.

These pin the current behaviour of the IO-free compute helpers so they can later
be extracted from the IO/subprocess-heavy RadiationSection without changing
results (task 7 prerequisite).
"""
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from _common import PYTHON_DIR  # noqa: E402

from udprep.udprep_radiation import RadiationSection  # noqa: E402

class TestInterpMakima(unittest.TestCase):
    """Modified-Akima interpolation (RadiationSection._interp_makima)."""

    interp = staticmethod(RadiationSection._interp_makima)

    def test_passes_through_data_points(self):
        x = np.array([0.0, 1, 2, 3, 4])
        y = np.array([0.0, 1, 4, 9, 16])
        np.testing.assert_allclose(self.interp(x, y, x), y)

    def test_reproduces_linear_data_exactly(self):
        x = np.array([0.0, 1, 2])
        y = np.array([0.0, 2, 4])
        np.testing.assert_allclose(self.interp(x, y, np.array([0.5, 1.5])), [1.0, 3.0])

    def test_known_values_on_quadratic_samples(self):
        x = np.array([0.0, 1, 2, 3, 4])
        y = np.array([0.0, 1, 4, 9, 16])
        np.testing.assert_allclose(
            self.interp(x, y, np.array([0.5, 2.5])), [0.3125, 6.239583], rtol=1e-6
        )

    def test_two_points_falls_back_to_linear(self):
        # Regression: two points used to IndexError on m[1]. Now linear.
        x = np.array([0.0, 1.0])
        y = np.array([0.0, 2.0])
        np.testing.assert_allclose(
            self.interp(x, y, np.array([0.0, 0.5, 1.0, 2.0])),  # incl. extrapolation
            [0.0, 1.0, 2.0, 4.0],
        )

    def test_rejects_degenerate_input(self):
        with self.assertRaises(ValueError):
            self.interp(np.array([0.0]), np.array([0.0]), np.array([0.0]))
        with self.assertRaises(ValueError):
            self.interp(np.array([0.0, 0.0]), np.array([1.0, 2.0]), np.array([0.0]))

class TestNetShortwaveReflections(unittest.TestCase):
    """Iterative multi-bounce net-shortwave physics (calc_reflections_sw).

    Ground-truth values captured from the current implementation; also checks the
    physical invariants (no-reflection limit, energy split) so the numerics are
    pinned before the compute is separated from the section's IO (task 7).
    """

    @staticmethod
    def _reflect(*args, **kwargs):
        # self is unused by calc_reflections_sw; a bare instance is enough.
        section = object.__new__(RadiationSection)
        return section.calc_reflections_sw(*args, **kwargs)

    def test_no_view_factors_reduces_to_absorbed_direct_plus_diffuse(self):
        # With vf=0 there are no reflections: knet = (1-albedo)*(sdir + dsky*svf).
        knet = self._reflect(
            np.array([100.0, 200.0]), 0.0, np.zeros((2, 2)),
            np.array([1.0, 1.0]), np.array([0.5, 0.5]),
        )
        np.testing.assert_allclose(knet, [50.0, 100.0])

    def test_mutual_reflection_between_two_facets(self):
        vf = np.array([[0.0, 0.5], [0.5, 0.0]])
        knet = self._reflect(
            np.array([100.0, 0.0]), 0.0, vf,
            np.array([0.5, 0.5]), np.array([0.5, 0.5]),
        )
        np.testing.assert_allclose(knet, [53.320312, 13.28125], rtol=1e-6)

    def test_diffuse_sky_contribution(self):
        vf = np.array([[0.0, 0.5], [0.5, 0.0]])
        knet = self._reflect(
            np.array([50.0, 50.0]), 100.0, vf,
            np.array([0.3, 0.7]), np.array([0.2, 0.2]),
        )
        np.testing.assert_allclose(knet, [74.24, 103.36], rtol=1e-6)

    def test_zero_albedo_absorbs_all_incoming(self):
        # albedo=0 -> nothing reflected -> knet == sdir + dsky*svf, one iteration.
        knet = self._reflect(
            np.array([80.0, 120.0]), 50.0, np.array([[0.0, 1.0], [1.0, 0.0]]),
            np.array([0.5, 0.5]), np.array([0.0, 0.0]),
        )
        np.testing.assert_allclose(knet, [80.0 + 50.0 * 0.5, 120.0 + 50.0 * 0.5])

    def test_validation(self):
        with self.assertRaises(ValueError):
            self._reflect(np.array([1.0, 2.0]), 0.0, np.zeros((2, 2)),
                          np.array([1.0]), np.array([0.5, 0.5]))  # shape mismatch
        with self.assertRaises(ValueError):
            self._reflect(np.array([1.0]), 0.0, np.zeros((1, 1)),
                          np.array([1.0]), np.array([0.5]), tol=0.0)

    def test_section_wrapper_delegates_to_pure_function(self):
        from udprep import _radiation_compute as rc

        args = (np.array([100.0, 0.0]), 0.0, np.array([[0.0, 0.5], [0.5, 0.0]]),
                np.array([0.5, 0.5]), np.array([0.5, 0.5]))
        np.testing.assert_allclose(
            self._reflect(*args), rc.net_shortwave_reflections(*args)
        )

class TestNonScatteringShortwave(unittest.TestCase):
    def test_absorbed_direct_plus_diffuse(self):
        from udprep import _radiation_compute as rc

        knet = rc.net_shortwave_nonscattering(
            np.array([100.0, 200.0]), 50.0, np.array([0.5, 1.0]), np.array([0.2, 0.4])
        )
        # (1-albedo) * (sdir + dsky*fss)
        np.testing.assert_allclose(knet, [0.8 * (100 + 25), 0.6 * (200 + 50)])


class TestSolverSolidMaskHandling(unittest.TestCase):
    """P4: the solver must not silently mis-compute when the IBM solid mask is
    absent (Moller truncated the traced volume; facsec produced sdir~0)."""

    @staticmethod
    def _grid(ktot=10):
        from types import SimpleNamespace

        return SimpleNamespace(
            Sc=None, ktot=ktot, zm=np.arange(float(ktot)),
            zsize=float(ktot), dzt=np.ones(ktot),
        )

    def test_ktot_covers_mesh_extent_when_solid_absent(self):
        from types import SimpleNamespace

        from udprep.directshortwave import _compute_ktot_and_z_edges

        grid = self._grid()
        # Without any hint, no-solid + no-veg collapses to ktot=2 (the old bug).
        ktot_bare, *_ = _compute_ktot_and_z_edges(grid, None)
        self.assertEqual(ktot_bare, 2)
        # With a mesh whose top is at z=5.5, the traced volume now reaches it.
        mesh = SimpleNamespace(vertices=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 5.5]]))
        ktot, z_edges, z_max, _ = _compute_ktot_and_z_edges(grid, None, mesh=mesh)
        self.assertGreater(ktot, ktot_bare)
        self.assertGreaterEqual(float(z_edges[-1]), 5.5)

    def test_facsec_raises_without_solid_mask(self):
        import trimesh
        from types import SimpleNamespace

        from exceptions import RadiationError
        from udprep.directshortwave import DirectShortwaveSolver

        box = trimesh.creation.box(extents=(2.0, 2.0, 2.0))
        sim = SimpleNamespace(
            geom=SimpleNamespace(stl=box),
            Sc=None, dx=1.0, dy=1.0, itot=4, jtot=4, ktot=6,
            zm=np.arange(6.0), zsize=6.0, dzt=np.ones(6),
        )
        with self.assertRaises(RadiationError):
            DirectShortwaveSolver(sim, method="facsec")


class TestShortwaveStaleness(unittest.TestCase):
    """P5: cached Sdir/netsw must be re-derived when the inputs change, so a
    sun-position sweep does not keep returning the first run's result."""

    def _section(self, tmp, nfcts=3):
        from types import SimpleNamespace

        from udprep.udprep_radiation import RadiationSection

        sim = SimpleNamespace(
            path=Path(tmp), expnr="001", stl_file="geom.stl",
            geom=SimpleNamespace(stl=SimpleNamespace(
                faces=np.zeros((nfcts, 3), dtype=int),
                face_normals=np.zeros((nfcts, 3)),
            )),
        )
        sim.assign_prop_to_fac = lambda _p: np.zeros(nfcts)
        section = RadiationSection("radiation", {}, sim=sim, defaults={})
        section.view3d_out = 0
        section.lvfsparse = False
        section.maxD = 100.0
        section.Dsky = 50.0
        # Isolate the file-cache logic from the numba solver / View3D.
        section.calc_view_factors = lambda **kw: (object(), np.zeros(nfcts), {})
        section.calc_reflections_sw = lambda sdir, *a, **k: sdir * 2.0
        return section

    def test_recomputes_when_sun_changes(self):
        with TemporaryDirectory() as tmp:
            section = self._section(tmp)
            calls = []

            def fake_direct(nsun, irr, method=None, **kw):
                calls.append(np.asarray(nsun, float).copy())
                return np.full(3, float(nsun[2]) * irr), np.zeros(0), {}

            section.calc_direct_sw = fake_direct
            sdir1, k1, _ = section.calc_short_wave(np.array([0.0, 0.0, 1.0]), 100.0)
            sdir2, k2, _ = section.calc_short_wave(np.array([0.0, 1.0, 0.5]), 100.0)
            self.assertEqual(len(calls), 2)               # recomputed, not cached-stale
            self.assertFalse(np.allclose(sdir1, sdir2))   # reflects the new sun
            np.testing.assert_allclose(k2, sdir2 * 2.0)   # netsw re-derived too

    def test_reuses_when_inputs_identical(self):
        with TemporaryDirectory() as tmp:
            section = self._section(tmp)
            calls = []

            def fake_direct(nsun, irr, method=None, **kw):
                calls.append(1)
                return np.full(3, 42.0), np.zeros(0), {}

            section.calc_direct_sw = fake_direct
            section.calc_short_wave(np.array([0.0, 0.0, 1.0]), 100.0)
            section.calc_short_wave(np.array([0.0, 0.0, 1.0]), 100.0)
            self.assertEqual(len(calls), 1)               # identical inputs -> reuse Sdir.txt


if __name__ == "__main__":
    unittest.main()
