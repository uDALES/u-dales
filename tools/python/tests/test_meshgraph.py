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

from udgeom._meshgraph import connected_components  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
