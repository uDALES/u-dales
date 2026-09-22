"""Small-scene tests for pedestrian receptor geometry and visibility."""

import types
import unittest

import numpy as np
import trimesh

from udcomf.udcomf_radiation import UDComfRadiation


def _case(mesh=None):
    solid = np.zeros((3, 2, 3), dtype=bool)
    solid[1, 1, 1] = True
    return types.SimpleNamespace(
        receptor_height=1.1,
        xt=np.array([0.0, 1.0, 2.0]),
        yt=np.array([0.0, 1.0]),
        zm=np.array([0.0, 1.0, 2.0]),
        zsize=3.0,
        Sc=solid,
        xazimuth=90.0,
        geom=types.SimpleNamespace(stl=mesh) if mesh is not None else None,
        facs={"typeid": np.ones(len(mesh.faces), dtype=int)} if mesh is not None else {},
    )


class TestReceptorGrid(unittest.TestCase):
    def test_uses_native_cell_centres_and_scalar_solid_mask(self):
        sim = _case()
        grid = UDComfRadiation(sim).receptor_grid()
        self.assertEqual(grid.valid.shape, (3, 2))
        self.assertFalse(grid.valid[1, 1])
        np.testing.assert_array_equal(grid.valid_indices(), [0, 1, 2, 4, 5])
        np.testing.assert_allclose(
            grid.points(np.array([0, 4, 5])),
            [[0.0, 0.0, 1.1], [2.0, 0.0, 1.1], [2.0, 1.0, 1.1]],
        )
        with self.assertRaisesRegex(ValueError, "solid or invalid"):
            grid.points(np.array([3]))

        sim.receptor_height = 2.1
        self.assertTrue(UDComfRadiation(sim).receptor_grid().valid.all())
        sim.Sc = None
        with self.assertRaisesRegex(ValueError, "solid_c.txt"):
            UDComfRadiation(sim).receptor_grid()

    def test_rejects_height_outside_grid(self):
        sim = _case()
        for height in (0.0, -1.0, np.nan, 3.0):
            with self.subTest(height=height), self.assertRaises(ValueError):
                sim.receptor_height = height
                UDComfRadiation(sim).receptor_grid()

    def test_true_north_uses_existing_solar_rotation_convention(self):
        sim = _case()
        normals = UDComfRadiation(sim).plane_normals()
        np.testing.assert_allclose(normals["upface"], [0, 0, 1], atol=1e-14)
        np.testing.assert_allclose(normals["downface"], [0, 0, -1], atol=1e-14)
        np.testing.assert_allclose(normals["northface"], [0, 1, 0], atol=1e-14)
        np.testing.assert_allclose(normals["eastface"], [1, 0, 0], atol=1e-14)
        np.testing.assert_allclose(normals["southface"], [0, -1, 0], atol=1e-14)
        np.testing.assert_allclose(normals["westface"], [-1, 0, 0], atol=1e-14)
        sim.xazimuth = 0.0
        normals = UDComfRadiation(sim).plane_normals()
        np.testing.assert_allclose(normals["northface"], [1, 0, 0], atol=1e-14)
        np.testing.assert_allclose(normals["eastface"], [0, -1, 0], atol=1e-14)


class TestVisibility(unittest.TestCase):
    @staticmethod
    def _two_walls():
        near = trimesh.creation.box(extents=(0.2, 2.0, 2.0))
        near.apply_translation((1.0, 0.0, 1.0))
        far = trimesh.creation.box(extents=(0.2, 2.0, 2.0))
        far.apply_translation((2.0, 0.0, 1.0))
        return trimesh.util.concatenate((near, far))

    def test_open_sky_and_ground_horizon(self):
        radiation = UDComfRadiation(_case())
        point = np.array([0.0, 0.0, 1.1])
        np.testing.assert_array_equal(
            radiation.sky_visibility(point, np.array([[0, 0, 1], [1, 0, 1], [0, 0, -1]])),
            [True, True, False],
        )
        self.assertTrue(radiation.direct_solar_visibility(point, 60.0, 0.0))
        self.assertFalse(radiation.direct_solar_visibility(point, 90.0, 0.0))
        self.assertFalse(radiation.direct_solar_visibility(point, 95.0, 0.0))

    def test_wall_blocks_sun_and_far_facets_but_not_clear_sky(self):
        mesh = self._two_walls()
        radiation = UDComfRadiation(_case(mesh))
        point = np.array([0.0, 0.0, 1.0])
        np.testing.assert_array_equal(
            radiation.sky_visibility(point, np.array([[0, 0, 1], [1, 0, 0.5], [-1, 0, 1]])),
            [True, False, True],
        )
        self.assertFalse(radiation.direct_solar_visibility(point, 60.0, 0.0))
        self.assertTrue(radiation.direct_solar_visibility(point, 60.0, 180.0))

        centers = mesh.triangles_center
        normals = mesh.face_normals
        near = np.flatnonzero((np.abs(centers[:, 0] - 0.9) < 1e-12) & (normals[:, 0] < -0.9))
        far = np.flatnonzero((np.abs(centers[:, 0] - 1.9) < 1e-12) & (normals[:, 0] < -0.9))
        backside = np.flatnonzero((np.abs(centers[:, 0] - 1.1) < 1e-12) & (normals[:, 0] > 0.9))
        self.assertTrue(len(near) and len(far) and len(backside))
        np.testing.assert_array_equal(
            radiation.facet_visibility(point, np.array([near[0], far[0], backside[0]])),
            [True, False, False],
        )
        sample_visibility = radiation.facet_sample_visibility(
            point, np.array([near[0], far[0]]),
            np.array([[0.6, 0.2, 0.2], [0.2, 0.6, 0.2], [0.2, 0.2, 0.6]]),
        )
        np.testing.assert_array_equal(sample_visibility, [[True] * 3, [False] * 3])
        with self.assertRaisesRegex(ValueError, "barycentric_samples"):
            radiation.facet_sample_visibility(point, np.array([near[0]]), np.array([[1, 1, 0]]))

    def test_requires_loaded_mesh_when_stl_is_configured(self):
        sim = _case()
        sim.stl_file = "city.stl"
        with self.assertRaisesRegex(ValueError, "Load the case geometry"):
            UDComfRadiation(sim).sky_visibility(np.array([0, 0, 1]), np.array([[0, 0, 1]]))

    def test_facet_hit_id_is_stable_at_large_map_coordinates(self):
        mesh = self._two_walls()
        mesh.apply_translation((1000.0, 2000.0, 0.0))
        radiation = UDComfRadiation(_case(mesh))
        centers = mesh.triangles_center
        normals = mesh.face_normals
        near = np.flatnonzero((np.abs(centers[:, 0] - 1000.9) < 1e-9) & (normals[:, 0] < -0.9))
        far = np.flatnonzero((np.abs(centers[:, 0] - 1001.9) < 1e-9) & (normals[:, 0] < -0.9))
        np.testing.assert_array_equal(
            radiation.facet_visibility(
                np.array([1000.0, 2000.0, 1.0]), np.array([near[0], far[0]])
            ),
            [True, False],
        )

    def test_rejects_facet_normals_with_inconsistent_order(self):
        mesh = self._two_walls()
        sim = _case(mesh)
        sim.facs["normals"] = mesh.face_normals.copy()
        sim.facs["normals"][0] *= -1
        with self.assertRaisesRegex(ValueError, "normals do not match"):
            UDComfRadiation(sim).facet_visibility(np.array([0, 0, 1]), np.array([0]))

    def test_one_facet_can_be_partly_hidden(self):
        target = trimesh.Trimesh(
            vertices=[[2, -1, 0], [2, 0, 2], [2, 1, 0]],
            faces=[[0, 1, 2]], process=False,
        )
        blocker = trimesh.Trimesh(
            vertices=[[1, -1, 0], [1, 0, 0], [1, 0, 2], [1, -1, 2]],
            faces=[[0, 1, 2], [0, 2, 3]], process=False,
        )
        radiation = UDComfRadiation(_case(trimesh.util.concatenate((target, blocker))))
        visible = radiation.facet_sample_visibility(
            np.array([0.0, 0.0, 1.0]), np.array([0]),
            np.array([[0.7, 0.2, 0.1], [0.1, 0.2, 0.7]]),
        )
        np.testing.assert_array_equal(visible, [[False, True]])
