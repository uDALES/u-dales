import sys
import unittest
from pathlib import Path

import numpy as np
import trimesh

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udgeom import UDGeom


def _ring_prism_mesh():
    outer_bottom = np.array(
        [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 4.0, 0.0], [0.0, 4.0, 0.0]]
    )
    inner_bottom = np.array(
        [[1.0, 1.0, 0.0], [3.0, 1.0, 0.0], [3.0, 3.0, 0.0], [1.0, 3.0, 0.0]]
    )
    outer_top = outer_bottom + np.array([0.0, 0.0, 2.0])
    inner_top = inner_bottom + np.array([0.0, 0.0, 2.0])
    vertices = np.vstack([outer_bottom, inner_bottom, outer_top, inner_top])

    def quad(a, b, c, d):
        return [[a, b, c], [a, c, d]]

    faces = []

    # Outer walls.
    faces += quad(0, 1, 9, 8)
    faces += quad(1, 2, 10, 9)
    faces += quad(2, 3, 11, 10)
    faces += quad(3, 0, 8, 11)

    # Inner courtyard walls.
    faces += quad(4, 12, 13, 5)
    faces += quad(5, 13, 14, 6)
    faces += quad(6, 14, 15, 7)
    faces += quad(7, 15, 12, 4)

    # Top ring.
    faces += quad(8, 9, 13, 12)
    faces += quad(9, 10, 14, 13)
    faces += quad(10, 11, 15, 14)
    faces += quad(11, 8, 12, 15)

    return trimesh.Trimesh(vertices=vertices, faces=np.asarray(faces, dtype=int), process=False)


class TestUDGeomOutlineLoops(unittest.TestCase):
    def test_calculate_outline2d_keeps_multiple_boundary_loops(self):
        geom = UDGeom(stl=_ring_prism_mesh())

        outlines = geom.calculate_outline2d()

        self.assertEqual(len(outlines), 1)
        outline = outlines[0]
        polygon_xy = np.asarray(outline["polygon"])[:, :2]

        expected_vertices = {
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0),
            (1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0),
        }
        observed_vertices = {tuple(row) for row in np.unique(polygon_xy, axis=0)}

        self.assertEqual(observed_vertices, expected_vertices)
        np.testing.assert_allclose(outline["centroid"][:2], np.array([2.0, 2.0]), atol=1e-12)


if __name__ == "__main__":
    unittest.main()
