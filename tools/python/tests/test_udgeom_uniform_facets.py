"""Tests for udgeom.uniform_facets (facet-size uniformisation).

The fixtures are small synthetic meshes built by hand so that every expected
outcome (pass-through, split, merge, type preservation, conformity) can be
checked exactly.
"""
import unittest

import numpy as np
import trimesh

from _common import PYTHON_DIR  # noqa: F401  (adds tools/python to sys.path)

import udgeom
from udgeom import (
    FacetCountError,
    UniformFacetsResult,
    check,
    create_cubes,
    target_area_from_edge,
    target_edge_from_area,
    uniformize_facets,
)


def _quad(a, b, c, d):
    return [[a, b, c], [a, c, d]]


def _box(x0, y0, dx, dy, h):
    """Open-bottom box: four walls and a roof, two triangles each."""
    v = np.array(
        [
            [x0, y0, 0.0], [x0 + dx, y0, 0.0], [x0 + dx, y0 + dy, 0.0], [x0, y0 + dy, 0.0],
            [x0, y0, h], [x0 + dx, y0, h], [x0 + dx, y0 + dy, h], [x0, y0 + dy, h],
        ],
        dtype=float,
    )
    f = _quad(0, 1, 5, 4) + _quad(1, 2, 6, 5) + _quad(2, 3, 7, 6) + _quad(3, 0, 4, 7) + _quad(4, 5, 6, 7)
    return trimesh.Trimesh(vertices=v, faces=np.asarray(f), process=False)


def _wall_with_slivers(n_cells=8, roof_depth=100.0):
    """A 30 m x 30 m wall in the y = 0 plane made of ``2 * n_cells**2`` small
    triangles (7 m^2 for 8 cells), plus a roof strip of ``2 * n_cells`` large
    triangles (type 2) that shares the wall's top edge vertex by vertex."""
    xs = np.linspace(0.0, 30.0, n_cells + 1)
    zs = np.linspace(0.0, 30.0, n_cells + 1)
    gx, gz = np.meshgrid(xs, zs, indexing="ij")
    vertices = np.column_stack([gx.ravel(), np.zeros(gx.size), gz.ravel()])

    def vid(i, k):
        return i * (n_cells + 1) + k

    faces = []
    for i in range(n_cells):
        for k in range(n_cells):
            # Wall normal towards -y (outside).
            faces.append([vid(i, k), vid(i + 1, k + 1), vid(i, k + 1)])
            faces.append([vid(i, k), vid(i + 1, k), vid(i + 1, k + 1)])
    far0 = len(vertices)
    vertices = np.vstack([vertices, np.column_stack([xs, np.full(n_cells + 1, -roof_depth), np.full(n_cells + 1, 30.0)])])
    for i in range(n_cells):
        # Roof normal towards +z.
        faces.append([vid(i, n_cells), far0 + i, far0 + i + 1])
        faces.append([vid(i, n_cells), far0 + i + 1, vid(i + 1, n_cells)])
    n_wall = 2 * n_cells ** 2
    types = np.concatenate([np.ones(n_wall, dtype=int), np.full(2 * n_cells, 2, dtype=int)])
    mesh = trimesh.Trimesh(vertices=vertices, faces=np.asarray(faces), process=False)
    assert np.all(mesh.face_normals[:n_wall, 1] < 0)
    assert np.all(mesh.face_normals[n_wall:, 2] > 0)
    return mesh, types


def _edge_face_counts(mesh):
    faces = np.asarray(mesh.faces)
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    return counts


def _has_t_junction(mesh, tol=1e-9):
    """True when some vertex lies strictly inside an edge of a face it does not belong to."""
    vertices = np.asarray(mesh.vertices)
    faces = np.asarray(mesh.faces)
    edges = np.unique(np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1), axis=0)
    a = vertices[edges[:, 0]]
    b = vertices[edges[:, 1]]
    ab = b - a
    ab_len2 = np.einsum("ij,ij->i", ab, ab)
    for vid, p in enumerate(vertices):
        mask = (edges[:, 0] != vid) & (edges[:, 1] != vid)
        t = np.einsum("ij,ij->i", p - a[mask], ab[mask]) / ab_len2[mask]
        inside = (t > 1e-6) & (t < 1 - 1e-6)
        if not inside.any():
            continue
        foot = a[mask][inside] + t[inside, None] * ab[mask][inside]
        if (np.linalg.norm(foot - p, axis=1) <= tol).any():
            return True
    return False


class TestExports(unittest.TestCase):
    def test_package_exports(self):
        for name in ("uniformize_facets", "UniformFacetsResult", "FacetCountError",
                     "target_edge_from_area", "target_area_from_edge"):
            self.assertIn(name, udgeom.__all__)
            self.assertTrue(hasattr(udgeom, name))

    def test_target_conversions_round_trip(self):
        edge = 7.5
        self.assertAlmostEqual(target_edge_from_area(target_area_from_edge(edge)), edge)
        self.assertAlmostEqual(target_area_from_edge(1.0), np.sqrt(3.0) / 4.0)


class TestPassThrough(unittest.TestCase):
    def test_uniform_mesh_is_returned_unchanged(self):
        mesh = _box(0.0, 0.0, 10.0, 10.0, 10.0)  # ten 50 m^2 triangles
        result = uniformize_facets(mesh, target_edge=target_edge_from_area(50.0))
        self.assertIsInstance(result, UniformFacetsResult)
        self.assertEqual(result.report["n_patches_remeshed"], 0)
        self.assertEqual(len(result.mesh.faces), len(mesh.faces))
        np.testing.assert_array_equal(np.sort(result.face_map), np.arange(len(mesh.faces)))
        self.assertIsNone(result.face_types)
        self.assertAlmostEqual(result.report["area_out"], result.report["area_in"])

    def test_invalid_settings_raise(self):
        mesh = _box(0.0, 0.0, 10.0, 10.0, 10.0)
        with self.assertRaises(ValueError):
            uniformize_facets(mesh, size_band=0.5)
        with self.assertRaises(ValueError):
            uniformize_facets(mesh, allow_split=False, allow_merge=False)
        with self.assertRaises(ValueError):
            uniformize_facets(mesh, face_types=np.ones(3, dtype=int))


class TestSplit(unittest.TestCase):
    def test_large_box_faces_are_split_within_band_and_conforming(self):
        mesh = _box(0.0, 0.0, 30.0, 30.0, 30.0)  # 450 m^2 triangles
        result = uniformize_facets(mesh, target_edge=10.0, max_faces=1000)
        report = result.report
        self.assertEqual(report["n_faces_split"], 10)
        self.assertEqual(report["n_out_above_band"], 0)
        self.assertEqual(report["n_out_below_band"], 0)
        self.assertGreater(len(result.mesh.faces), len(mesh.faces))
        self.assertAlmostEqual(report["area_out"], report["area_in"], places=6)
        # No unused vertices, every edge used by one or two faces, no T-junctions.
        self.assertEqual(len(np.unique(result.mesh.faces)), len(result.mesh.vertices))
        self.assertTrue(np.all(_edge_face_counts(result.mesh) <= 2))
        self.assertFalse(_has_t_junction(result.mesh))
        # Normals keep pointing outward (roof up, walls outward).
        roof = np.isclose(result.mesh.triangles_center[:, 2], 30.0)
        self.assertTrue(np.all(result.mesh.face_normals[roof][:, 2] > 0.99))
        south = np.isclose(result.mesh.triangles_center[:, 1], 0.0)
        self.assertTrue(np.all(result.mesh.face_normals[south][:, 1] < -0.99))
        # Face map points into the input faces.
        self.assertEqual(result.face_map.min(), 0)
        self.assertEqual(result.face_map.max(), len(mesh.faces) - 1)

    def test_split_is_idempotent(self):
        mesh = _box(0.0, 0.0, 30.0, 30.0, 30.0)
        first = uniformize_facets(mesh, target_edge=10.0, max_faces=1000)
        second = uniformize_facets(first.mesh, target_edge=10.0)
        self.assertEqual(second.report["n_patches_remeshed"], 0)
        self.assertEqual(len(second.mesh.faces), len(first.mesh.faces))

    def test_allow_split_false_keeps_large_faces(self):
        mesh = _box(0.0, 0.0, 30.0, 30.0, 30.0)
        result = uniformize_facets(mesh, target_edge=10.0, allow_split=False)
        self.assertEqual(result.report["n_patches_remeshed"], 0)
        self.assertEqual(len(result.mesh.faces), len(mesh.faces))

    def test_count_guard(self):
        mesh = _box(0.0, 0.0, 30.0, 30.0, 30.0)
        with self.assertRaises(FacetCountError):
            uniformize_facets(mesh, target_edge=1.0, max_faces=100)
        with self.assertRaises(FacetCountError):
            # Prediction passes, the exact count must still be checked.
            uniformize_facets(mesh, target_edge=10.0, max_faces=11)


class TestMerge(unittest.TestCase):
    def test_slivers_merge_into_their_wall_and_types_are_preserved(self):
        mesh, types = _wall_with_slivers()
        result = uniformize_facets(mesh, face_types=types, target_edge=10.0)
        report = result.report
        self.assertEqual(report["n_patches"], 2)
        self.assertEqual(report["n_faces_merged"], 128)
        self.assertEqual(report["n_faces_split"], 16)
        self.assertLess(len(result.mesh.faces), 120)
        self.assertEqual(report["n_out_below_band"], 0)
        self.assertEqual(report["n_out_above_band"], 0)
        self.assertGreater(report["n_boundary_vertices_removed"], 0)
        self.assertAlmostEqual(report["area_out"], report["area_in"], places=6)
        # Types: wall faces are type 1, roof faces type 2, and the roof normal is unchanged.
        wall = np.isclose(result.mesh.triangles_center[:, 1], 0.0)
        self.assertTrue(np.all(result.face_types[wall] == 1))
        self.assertTrue(np.all(result.face_types[~wall] == 2))
        self.assertTrue(np.all(result.mesh.face_normals[~wall][:, 2] > 0.99))
        self.assertTrue(np.all(result.mesh.face_normals[wall][:, 1] < -0.99))
        self.assertEqual(len(np.unique(result.mesh.faces)), len(result.mesh.vertices))
        self.assertFalse(_has_t_junction(result.mesh))
        # Shared top edge: still exactly two faces per edge there.
        self.assertTrue(np.all(_edge_face_counts(result.mesh) <= 2))

    def test_faces_of_different_types_are_not_merged(self):
        mesh, types = _wall_with_slivers()
        # Give the lower half of the wall another type: two patches on one plane.
        lower = mesh.triangles_center[:, 2] < 15.0 - 1e-9
        types = types.copy()
        types[lower] = 3
        result = uniformize_facets(mesh, face_types=types, target_edge=10.0)
        self.assertEqual(result.report["n_patches"], 3)
        wall = np.isclose(result.mesh.triangles_center[:, 1], 0.0)
        z = result.mesh.triangles_center[:, 2]
        self.assertTrue(np.all(result.face_types[wall & (z < 15.0)] == 3))
        self.assertTrue(np.all(result.face_types[wall & (z > 15.0)] == 1))
        # The type boundary (z = 15) survives as a straight edge: no face straddles it.
        tri_z = result.mesh.vertices[result.mesh.faces][:, :, 2]
        straddle = wall & (tri_z.min(axis=1) < 15.0 - 1e-9) & (tri_z.max(axis=1) > 15.0 + 1e-9)
        self.assertFalse(straddle.any())

    def test_allow_merge_false_leaves_slivers(self):
        mesh, types = _wall_with_slivers()
        result = uniformize_facets(mesh, face_types=types, target_edge=10.0, allow_merge=False)
        self.assertEqual(result.report["n_faces_merged"], 0)
        self.assertEqual(result.report["n_small_isolated"], 128)
        self.assertGreaterEqual(len(result.mesh.faces), 128)

    def test_isolated_small_face_is_kept(self):
        # A single tiny triangle floating on its own plane cannot be merged.
        mesh = trimesh.util.concatenate([
            _box(0.0, 0.0, 10.0, 10.0, 10.0),
            trimesh.Trimesh(vertices=[[50, 50, 0], [50.1, 50, 0], [50, 50.1, 0]], faces=[[0, 1, 2]], process=False),
        ])
        result = uniformize_facets(mesh, target_edge=target_edge_from_area(50.0))
        self.assertEqual(result.report["n_small_isolated"], 1)
        self.assertEqual(len(result.mesh.faces), 11)

    def test_degenerate_faces_are_removed(self):
        box = _box(0.0, 0.0, 10.0, 10.0, 10.0)
        faces = np.vstack([box.faces, [[0, 1, 1]]])
        mesh = trimesh.Trimesh(vertices=box.vertices, faces=faces, process=False)
        result = uniformize_facets(mesh, target_edge=target_edge_from_area(50.0))
        self.assertEqual(result.report["n_degenerate_removed"], 1)
        self.assertEqual(len(result.mesh.faces), 10)


class TestWorkers(unittest.TestCase):
    """Worker processes must give exactly the serial result."""

    @staticmethod
    def _assert_same(a, b):
        np.testing.assert_array_equal(a.mesh.vertices, b.mesh.vertices)
        np.testing.assert_array_equal(a.mesh.faces, b.mesh.faces)
        np.testing.assert_array_equal(a.face_map, b.face_map)
        if a.face_types is not None:
            np.testing.assert_array_equal(a.face_types, b.face_types)
        for key in ("n_faces_out", "n_patches_remeshed", "n_patches_kept", "n_boundary_points_remeshed"):
            assert a.report[key] == b.report[key], key

    def test_workers_give_the_serial_result_on_a_small_mesh(self):
        import udgeom.uniform_facets as impl

        mesh, types = _wall_with_slivers()
        serial = uniformize_facets(mesh, face_types=types, target_edge=10.0)
        saved = impl._MIN_PATCHES_PER_WORKER
        impl._MIN_PATCHES_PER_WORKER = 1  # two patches: make the pool run anyway
        try:
            parallel = uniformize_facets(mesh, face_types=types, target_edge=10.0, workers=2)
        finally:
            impl._MIN_PATCHES_PER_WORKER = saved
        self._assert_same(parallel, serial)

    def test_workers_give_the_serial_result_on_many_patches(self):
        import udgeom.uniform_facets as impl

        mesh = create_cubes(240.0, 180.0, 30.0, 30.0, 30.0, 30.0, 30.0, "AC", 10.0).stl
        serial = uniformize_facets(mesh, target_edge=10.0)
        self.assertGreater(serial.report["n_patches_remeshed"], 16)
        saved = impl._MIN_PATCHES_PER_WORKER
        impl._MIN_PATCHES_PER_WORKER = 8  # enough patches for two workers
        try:
            parallel = uniformize_facets(mesh, target_edge=10.0, workers=2)
        finally:
            impl._MIN_PATCHES_PER_WORKER = saved
        self._assert_same(parallel, serial)

    def test_invalid_workers_raise(self):
        with self.assertRaises(ValueError):
            uniformize_facets(_box(0.0, 0.0, 10.0, 10.0, 10.0), workers=0)


class TestRealisticGeometry(unittest.TestCase):
    def test_cube_array_with_ground_stays_valid(self):
        geom = create_cubes(240.0, 180.0, 30.0, 30.0, 30.0, 30.0, 30.0, "AC", 10.0)
        mesh = geom.stl
        result = uniformize_facets(mesh)  # fixed-count default
        self.assertTrue(result.report["target_auto"])
        self.assertEqual(result.report["n_out_below_band"], 0)
        self.assertEqual(result.report["n_out_above_band"], 0)
        self.assertLess(len(result.mesh.faces), 1.5 * len(mesh.faces))
        report = check(result.mesh, require_single_component=False)
        self.assertTrue(report["valid"], report["issues"])
        self.assertFalse(_has_t_junction(result.mesh))
        again = uniformize_facets(result.mesh, target_edge=result.report["target_edge"])
        self.assertEqual(again.report["n_patches_remeshed"], 0)

    def test_udgeom_method_wrapper(self):
        geom = create_cubes(240.0, 180.0, 30.0, 30.0, 30.0, 30.0, 30.0, "AC", 10.0)
        n_before = geom.n_faces
        result = geom.uniformize_facets(target_edge=10.0, inplace=True)
        self.assertEqual(geom.n_faces, len(result.mesh.faces))
        self.assertNotEqual(geom.n_faces, n_before)


if __name__ == "__main__":
    unittest.main()
