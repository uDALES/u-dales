import sys
import unittest
import warnings
from pathlib import Path

import numpy as np
import trimesh

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR, REPO_ROOT, copy_case

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udgeom import UDGeom


DATA_DIR = REPO_ROOT / "tools" / "python" / "tests" / "data" / "udgeom_matlab"


class TestUDGeomApi(unittest.TestCase):
    @staticmethod
    def _wedge_mesh():
        vertices = np.array(
            [
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [0.0, 0.0, 2.0],
            ],
            dtype=float,
        )
        faces = np.array([[0, 1, 2], [0, 3, 1]], dtype=int)
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    def test_constructor_accepts_trimesh(self):
        mesh = trimesh.creation.box(extents=(1.0, 2.0, 3.0))
        geom = UDGeom(stl=mesh)

        self.assertIs(geom.stl, mesh)
        self.assertEqual(geom.n_faces, len(mesh.faces))
        self.assertEqual(geom.n_vertices, len(mesh.vertices))

    def test_constructor_rejects_invalid_stl(self):
        with self.assertRaises(TypeError):
            UDGeom(stl="not-a-mesh")

    def test_load_concatenates_trimesh_scene(self):
        temp_dir, case_dir = copy_case(DATA_DIR)
        self.addCleanup(temp_dir.cleanup)

        box_a = trimesh.creation.box(extents=(1.0, 1.0, 1.0))
        box_b = trimesh.creation.box(extents=(1.0, 1.0, 1.0))
        box_b.apply_translation((3.0, 0.0, 0.0))

        scene = trimesh.Scene()
        scene.add_geometry(box_a)
        scene.add_geometry(box_b)
        scene_path = case_dir / "scene.glb"
        scene.export(scene_path)

        geom = UDGeom(case_dir)
        geom.load("scene.glb")

        self.assertEqual(geom.n_faces, len(box_a.faces) + len(box_b.faces))
        self.assertEqual(geom.n_vertices, len(box_a.vertices) + len(box_b.vertices))

    def test_load_missing_file_raises(self):
        geom = UDGeom(DATA_DIR)
        with self.assertRaises(FileNotFoundError):
            geom.load("missing.stl")

    def test_save_without_geometry_raises(self):
        geom = UDGeom(DATA_DIR)
        with self.assertRaises(ValueError):
            geom.save("nope.stl")

    def test_save_round_trip_preserves_geometry(self):
        temp_dir, case_dir = copy_case(DATA_DIR)
        self.addCleanup(temp_dir.cleanup)

        geom = UDGeom(DATA_DIR)
        geom.load("single_box.stl")
        geom.path = case_dir
        geom.save("saved.stl")

        reloaded = UDGeom(case_dir)
        reloaded.load("saved.stl")

        np.testing.assert_allclose(reloaded.bounds, geom.bounds, atol=1e-12)
        np.testing.assert_allclose(reloaded.face_centers, geom.face_centers, atol=1e-12)
        self.assertEqual(reloaded.n_faces, geom.n_faces)
        self.assertEqual(reloaded.n_vertices, geom.n_vertices)

    def test_methods_without_geometry_follow_contract(self):
        geom = UDGeom()

        self.assertEqual(geom.n_faces, 0)
        self.assertEqual(geom.n_vertices, 0)
        np.testing.assert_array_equal(geom.bounds, np.array([[0, 0, 0], [0, 0, 0]]))
        np.testing.assert_array_equal(geom.face_centers, np.array([]))
        np.testing.assert_array_equal(geom.face_incenters, np.array([]))
        np.testing.assert_array_equal(geom.face_normals, np.array([]))
        np.testing.assert_array_equal(geom.face_areas, np.array([]))
        self.assertEqual(geom.total_area, 0.0)
        self.assertEqual(geom.volume, 0.0)
        self.assertFalse(geom.is_watertight)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self.assertEqual(geom.get_buildings(), [])
            self.assertEqual(len(caught), 1)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            np.testing.assert_array_equal(geom.get_face_to_building_map(), np.array([]))
            self.assertEqual(len(caught), 1)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self.assertEqual(geom.get_building_outlines(), [])
            self.assertEqual(len(caught), 1)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            np.testing.assert_array_equal(geom.get_outline(), np.empty((0, 2), dtype=int))
            self.assertEqual(len(caught), 1)

    def test_volume_warns_for_non_watertight_mesh(self):
        geom = UDGeom(DATA_DIR)
        geom.load("flat_ground.stl")

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            volume = geom.volume

        self.assertGreaterEqual(len(caught), 1)
        self.assertTrue(any("not watertight" in str(w.message).lower() for w in caught))
        self.assertIsInstance(volume, float)

    def test_building_related_results_are_cached(self):
        geom = UDGeom(DATA_DIR)
        geom.load("xie_castro_2008.stl")

        buildings_a = geom.get_buildings()
        buildings_b = geom.get_buildings()
        self.assertIs(buildings_a, buildings_b)

        map_a = geom.get_face_to_building_map()
        map_b = geom.get_face_to_building_map()
        self.assertIs(map_a, map_b)

        outlines_a = geom.calculate_outline2d()
        outlines_b = geom.calculate_outline2d()
        self.assertIs(outlines_a, outlines_b)

    def test_load_invalidates_caches(self):
        geom = UDGeom(DATA_DIR)
        geom.load("xie_castro_2008.stl")

        old_buildings = geom.get_buildings()
        old_map = geom.get_face_to_building_map()
        old_outline2d = geom.calculate_outline2d()

        geom.load("single_box.stl")

        self.assertIsNone(geom._buildings)
        self.assertIsNone(geom._face_to_building_map)
        self.assertIsNone(geom._outline2d)
        self.assertIsNot(old_buildings, geom.get_buildings())
        self.assertIsNot(old_map, geom.get_face_to_building_map())
        self.assertIsNot(old_outline2d, geom.calculate_outline2d())

    def test_outline_threshold_argument_is_supported(self):
        geom = UDGeom(stl=self._wedge_mesh())

        default_edges = np.asarray(sorted(map(tuple, geom.get_outline())), dtype=int)
        strict_edges = np.asarray(sorted(map(tuple, geom.get_outline(120.0))), dtype=int)
        self.assertFalse(np.array_equal(default_edges, strict_edges))

        default_building = geom.get_building_outlines()
        strict_building = geom.get_building_outlines(120.0)
        self.assertEqual(len(default_building), len(strict_building))
        for edges_a, edges_b in zip(default_building, strict_building):
            self.assertFalse(np.array_equal(np.asarray(edges_a), np.asarray(edges_b)))

    def test_repr_contains_geometry_summary(self):
        geom = UDGeom(DATA_DIR)
        self.assertIn("stl=None", repr(geom))

        geom.load("single_box.stl")
        text = repr(geom)
        self.assertIn("n_faces=12", text)
        self.assertIn("n_vertices=8", text)


if __name__ == "__main__":
    unittest.main()
