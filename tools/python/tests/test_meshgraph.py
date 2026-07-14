"""Unit tests for the shared face-adjacency graph primitive."""
import sys
import unittest
from pathlib import Path

import numpy as np

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR  # noqa: E402

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udgeom._meshgraph import (  # noqa: E402
    build_edge_face_adjacency,
    connected_components,
)

try:
    import trimesh  # noqa: E402

    HAVE_TRIMESH = True
except ImportError:
    HAVE_TRIMESH = False


class TestConnectedComponents(unittest.TestCase):
    @staticmethod
    def _adj(pairs, n):
        adj = {i: set() for i in range(n)}
        for a, b in pairs:
            adj[a].add(b)
            adj[b].add(a)
        return adj

    def test_two_components_in_seed_order(self):
        # 0-1-2 and 3-4 ; default seeds = ascending node order.
        adj = self._adj([(0, 1), (1, 2), (3, 4)], 5)
        comps = connected_components(adj)
        self.assertEqual([sorted(c.tolist()) for c in comps], [[0, 1, 2], [3, 4]])

    def test_isolated_node_is_its_own_component(self):
        adj = self._adj([(0, 1)], 3)  # node 2 isolated
        comps = connected_components(adj)
        self.assertEqual(len(comps), 2)
        self.assertIn([2], [sorted(c.tolist()) for c in comps])

    def test_seeds_limit_and_order_components(self):
        adj = self._adj([(0, 1), (3, 4)], 5)
        comps = connected_components(adj, seeds=[3, 0])  # start at the second group
        self.assertEqual(sorted(comps[0].tolist()), [3, 4])

    def test_keep_restricts_traversal(self):
        # chain 0-1-2-3, but keep excludes node 2 -> traversal stops at the gap.
        adj = self._adj([(0, 1), (1, 2), (2, 3)], 4)
        keep = np.array([True, True, False, True])
        comps = connected_components(adj, seeds=[0], keep=keep)
        self.assertEqual(sorted(comps[0].tolist()), [0, 1])  # 2 (and thus 3) not reached


class TestBuildEdgeFaceAdjacency(unittest.TestCase):
    def test_two_triangles_sharing_an_edge_are_adjacent(self):
        # Faces 0 and 1 share edge (1,2); faces 1 and 2 share edge (2,3).
        faces = np.array([[0, 1, 2], [1, 2, 3], [2, 3, 4]])
        adj = build_edge_face_adjacency(faces)
        self.assertEqual(adj[0], {1})
        self.assertEqual(adj[1], {0, 2})
        self.assertEqual(adj[2], {1})

    def test_non_manifold_edge_connects_all_sharing_faces(self):
        # Edge (0,1) is shared by three faces -> all mutually adjacent.
        faces = np.array([[0, 1, 2], [0, 1, 3], [0, 1, 4]])
        adj = build_edge_face_adjacency(faces)
        self.assertEqual(adj[0], {1, 2})
        self.assertEqual(adj[1], {0, 2})
        self.assertEqual(adj[2], {0, 1})

    def test_disjoint_triangles_have_no_edges(self):
        faces = np.array([[0, 1, 2], [3, 4, 5]])
        adj = build_edge_face_adjacency(faces)
        self.assertEqual(adj, {0: set(), 1: set()})

    def test_components_over_edge_adjacency(self):
        # Two edge-connected pairs -> two components.
        faces = np.array([[0, 1, 2], [1, 2, 3], [4, 5, 6], [5, 6, 7]])
        comps = connected_components(build_edge_face_adjacency(faces))
        self.assertEqual(
            sorted(sorted(c.tolist()) for c in comps), [[0, 1], [2, 3]]
        )


@unittest.skipUnless(HAVE_TRIMESH, "trimesh not installed")
class TestSplitBuildingsConsolidation(unittest.TestCase):
    """Both building splitters group two separated boxes into two buildings."""

    @staticmethod
    def _two_boxes():
        a = trimesh.creation.box(extents=(1, 1, 1))
        a.apply_translation((0, 0, 0.5))
        b = trimesh.creation.box(extents=(1, 1, 1))
        b.apply_translation((5, 5, 0.5))  # well separated in x+y
        return trimesh.util.concatenate([a, b])

    def test_split_buildings_function_finds_two(self):
        from udgeom.split_buildings import split_buildings

        buildings, face_map = split_buildings(self._two_boxes(), remove_ground=False)
        self.assertEqual(len(buildings), 2)
        self.assertEqual(set(np.unique(face_map).tolist()), {1, 2})

    def test_udgeom_split_finds_two_sorted_southwest_first(self):
        from udgeom.udgeom import UDGeom

        geom = UDGeom.__new__(UDGeom)
        geom.stl = self._two_boxes()
        geom._buildings = None
        geom._face_to_building_map = None
        geom._outline2d = None
        geom._outline3d = None
        geom._split_buildings()
        self.assertEqual(len(geom._buildings), 2)
        # Building ids 1 and 2 assigned; id 0 marks the z=0 box-bottom faces.
        self.assertTrue(
            {1, 2}.issubset(np.unique(geom._face_to_building_map).tolist())
        )
        # Building 1 is the south-west box (smaller x+y centroid).
        c1 = geom._buildings[0].vertices[:, :2].mean(axis=0).sum()
        c2 = geom._buildings[1].vertices[:, :2].mean(axis=0).sum()
        self.assertLess(c1, c2)


if __name__ == "__main__":
    unittest.main()
