import json
import sys
import unittest
from pathlib import Path
from typing import Dict, List

import numpy as np

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR, REPO_ROOT

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udgeom import UDGeom


DATA_DIR = REPO_ROOT / "tools" / "python" / "tests" / "data" / "udgeom_matlab"


def _load_reference(name: str) -> Dict:
    return json.loads((DATA_DIR / f"{name}.json").read_text())


def _normalize_polygon(polygon: np.ndarray) -> np.ndarray:
    polygon = np.asarray(polygon, dtype=float)
    if polygon.size == 0:
        return polygon.reshape(0, 2)
    polygon = polygon[:, :2]
    if len(polygon) > 1 and np.allclose(polygon[0], polygon[-1]):
        polygon = polygon[:-1]
    order = np.lexsort((polygon[:, 1], polygon[:, 0]))
    return polygon[order]


def _building_summary(geom: UDGeom, face_map: np.ndarray) -> List[Dict]:
    buildings = geom.get_buildings()
    summaries = []
    for building_id, building in enumerate(buildings, start=1):
        face_indices = np.where(face_map == building_id)[0] + 1
        summaries.append(
            {
                "n_faces": int(len(building.faces)),
                "n_vertices": int(len(building.vertices)),
                "bounds": np.asarray(building.bounds, dtype=float),
                "original_face_indices": face_indices,
            }
        )
    return summaries


class TestUDGeomAgainstMatlab(unittest.TestCase):
    fixtures = ("flat_ground", "single_box", "two_boxes_with_ground", "xie_castro_2008")

    def _load_geom(self, name: str) -> UDGeom:
        geom = UDGeom(DATA_DIR)
        geom.load(f"{name}.stl")
        return geom

    def test_reference_files_exist(self):
        self.assertTrue((DATA_DIR / "manifest.json").exists())
        for name in self.fixtures:
            self.assertTrue((DATA_DIR / f"{name}.stl").exists(), name)
            self.assertTrue((DATA_DIR / f"{name}.json").exists(), name)

    def test_basic_mesh_properties_match_matlab(self):
        for name in self.fixtures:
            with self.subTest(name=name):
                geom = self._load_geom(name)
                ref = _load_reference(name)

                self.assertEqual(geom.n_faces, ref["n_faces"])
                self.assertEqual(geom.n_vertices, ref["n_vertices"])
                np.testing.assert_allclose(geom.bounds, ref["bounds"], atol=1e-12)
                np.testing.assert_allclose(geom.face_centers, ref["face_centers"], atol=1e-12)
                np.testing.assert_allclose(geom.face_incenters, ref["face_incenters"], atol=1e-12)
                np.testing.assert_allclose(geom.face_normals, ref["face_normals"], atol=1e-12)
                np.testing.assert_allclose(geom.face_areas, ref["face_areas"], atol=1e-12)
                self.assertAlmostEqual(geom.total_area, ref["total_area"], places=12)
                self.assertEqual(bool(geom.is_watertight), bool(ref["is_watertight"]))
                if ref["is_watertight"]:
                    self.assertAlmostEqual(geom.volume, ref["volume"], places=12)

    def test_building_partition_matches_matlab(self):
        for name in self.fixtures:
            with self.subTest(name=name):
                geom = self._load_geom(name)
                ref = _load_reference(name)
                face_map = geom.get_face_to_building_map()
                np.testing.assert_array_equal(face_map, ref["face_to_building_map"])

                summaries = _building_summary(geom, face_map)
                self.assertEqual(len(summaries), len(ref["buildings"]))
                for summary, ref_summary in zip(summaries, ref["buildings"]):
                    self.assertEqual(summary["n_faces"], ref_summary["n_faces"])
                    self.assertEqual(summary["n_vertices"], ref_summary["n_vertices"])
                    np.testing.assert_allclose(summary["bounds"], ref_summary["bounds"], atol=1e-12)
                    np.testing.assert_array_equal(
                        summary["original_face_indices"],
                        np.asarray(ref_summary["original_face_indices"], dtype=int),
                    )

    def test_outline2d_matches_matlab_after_normalization(self):
        for name in self.fixtures:
            with self.subTest(name=name):
                geom = self._load_geom(name)
                outlines = geom.calculate_outline2d()
                ref = _load_reference(name)

                self.assertEqual(len(outlines), len(ref["outline2d"]))
                for outline, ref_outline in zip(outlines, ref["outline2d"]):
                    np.testing.assert_allclose(
                        np.asarray(outline["centroid"], dtype=float)[:2],
                        np.asarray(ref_outline["centroid"], dtype=float)[:2],
                        atol=1e-12,
                    )
                    np.testing.assert_allclose(
                        _normalize_polygon(np.asarray(outline["polygon"], dtype=float)),
                        _normalize_polygon(np.asarray(ref_outline["polygon"], dtype=float)),
                        atol=1e-12,
                    )

    def test_outline_edges_match_matlab(self):
        for name in self.fixtures:
            with self.subTest(name=name):
                geom = self._load_geom(name)
                edges = np.asarray(sorted(map(tuple, geom.get_outline())), dtype=int)
                ref = _load_reference(name)
                ref_edges = np.asarray(sorted(map(tuple, ref["outline_edges"])), dtype=int) - 1
                np.testing.assert_array_equal(edges, ref_edges)

    def test_building_outline_edges_match_matlab(self):
        for name in self.fixtures:
            with self.subTest(name=name):
                geom = self._load_geom(name)
                outlines = geom.get_building_outlines()
                ref = _load_reference(name)
                self.assertEqual(len(outlines), len(ref["building_outlines"]))
                for outline, ref_outline in zip(outlines, ref["building_outlines"]):
                    outline_edges = np.asarray(sorted(map(tuple, outline)), dtype=int)
                    ref_edges = np.asarray(sorted(map(tuple, ref_outline)), dtype=int) - 1
                    np.testing.assert_array_equal(outline_edges, ref_edges)


if __name__ == "__main__":
    unittest.main()
