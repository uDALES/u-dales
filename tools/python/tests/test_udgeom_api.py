import sys
import struct
import unittest
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import trimesh
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR, REPO_ROOT, copy_case

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udgeom import (
    UDGeom,
    add_ground,
    calculate_independent_surfaces,
    check,
    create_canyons,
    delete_ground,
    extrude_to_ground,
    find_internal_touching_wall_regions,
    find_nonmanifold_regions,
    find_unstitched_touching_regions,
    fix,
    identify_ground_faces,
    repair_adjacent_buildings,
    resolve_vertical_coplanar_overlaps,
    truncate_below_ground,
)
from udgeom.geometry_generation import _ground_footprint_union
from udgeom.split_buildings import split_buildings as geom_split_buildings


DATA_DIR = REPO_ROOT / "tools" / "python" / "tests" / "data" / "udgeom_matlab"
EXAMPLES_DIR = REPO_ROOT / "tools" / "python" / "examples"


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

    @staticmethod
    def _open_bottom_box_mesh():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 3.0, 0.0],
                [0.0, 3.0, 0.0],
                [0.0, 0.0, 4.0],
                [2.0, 0.0, 4.0],
                [2.0, 3.0, 4.0],
                [0.0, 3.0, 4.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 5], [0, 5, 4],
                [1, 2, 6], [1, 6, 5],
                [2, 3, 7], [2, 7, 6],
                [3, 0, 4], [3, 4, 7],
                [4, 5, 6], [4, 6, 7],
            ],
            dtype=int,
        )
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _ring_prism_mesh():
        outer_bottom = np.array(
            [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 4.0, 0.0], [0.0, 4.0, 0.0]],
            dtype=float,
        )
        inner_bottom = np.array(
            [[1.0, 1.0, 0.0], [3.0, 1.0, 0.0], [3.0, 3.0, 0.0], [1.0, 3.0, 0.0]],
            dtype=float,
        )
        outer_top = outer_bottom + np.array([0.0, 0.0, 2.0])
        inner_top = inner_bottom + np.array([0.0, 0.0, 2.0])
        vertices = np.vstack([outer_bottom, inner_bottom, outer_top, inner_top])

        def quad(a, b, c, d):
            return [[a, b, c], [a, c, d]]

        faces = []
        faces += quad(0, 1, 9, 8)
        faces += quad(1, 2, 10, 9)
        faces += quad(2, 3, 11, 10)
        faces += quad(3, 0, 8, 11)
        faces += quad(4, 12, 13, 5)
        faces += quad(5, 13, 14, 6)
        faces += quad(6, 14, 15, 7)
        faces += quad(7, 15, 12, 4)
        faces += quad(8, 9, 13, 12)
        faces += quad(9, 10, 14, 13)
        faces += quad(10, 11, 15, 14)
        faces += quad(11, 8, 12, 15)
        return trimesh.Trimesh(vertices=vertices, faces=np.asarray(faces, dtype=int), process=False)

    @staticmethod
    def _mesh_with_unused_vertex():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [5.0, 5.0, 5.0],
            ],
            dtype=float,
        )
        faces = np.array([[0, 1, 2]], dtype=int)
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _disconnected_two_triangle_mesh():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [5.0, 0.0, 0.0],
                [6.0, 0.0, 0.0],
                [5.0, 1.0, 0.0],
            ],
            dtype=float,
        )
        faces = np.array([[0, 1, 2], [3, 4, 5]], dtype=int)
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _mesh_with_duplicate_face():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            dtype=float,
        )
        faces = np.array([[0, 1, 2], [2, 1, 0]], dtype=int)
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _mesh_with_nonfinite_vertex():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [np.nan, 1.0, 0.0],
            ],
            dtype=float,
        )
        faces = np.array([[0, 1, 2]], dtype=int)
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _mesh_with_zero_area_face():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ],
            dtype=float,
        )
        faces = np.array([[0, 1, 2]], dtype=int)
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _mesh_with_downward_ground_face():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
            ],
            dtype=float,
        )
        faces = np.array([[0, 1, 2]], dtype=int)
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _mesh_with_topography_and_bridge():
        vertices = np.array(
            [
                [0.0, 0.0, 5.0],
                [2.0, 0.0, 5.0],
                [2.0, 2.0, 5.0],
                [0.0, 2.0, 5.0],
                [0.5, 0.5, 8.0],
                [1.5, 0.5, 8.0],
                [1.5, 1.5, 8.0],
                [0.5, 1.5, 8.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],
                [4, 5, 6],
                [4, 6, 7],
            ],
            dtype=int,
        )
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _mesh_with_nonmanifold_edge():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, -1.0, 0.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 2],
                [0, 3, 1],
                [0, 1, 4],
            ],
            dtype=int,
        )
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _mesh_with_unstitched_touching_edge():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, -1.0, 0.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 2],
                [3, 5, 4],
            ],
            dtype=int,
        )
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _mesh_with_unstitched_touching_edge_at_large_coordinates():
        mesh = TestUDGeomApi._mesh_with_unstitched_touching_edge().copy()
        mesh.vertices += np.array([1.0e5, 2.0e5, 0.0], dtype=float)
        return mesh

    @staticmethod
    def _mesh_with_unstitched_ground_building_t_junction():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],  # 0
                [1.0, 0.0, 0.0],  # 1 midpoint only on ground
                [2.0, 0.0, 0.0],  # 2
                [0.0, 1.0, 0.0],  # 3
                [1.0, 1.0, 0.0],  # 4
                [2.0, 1.0, 0.0],  # 5
                [0.0, 0.0, 1.0],  # 6
                [2.0, 0.0, 1.0],  # 7
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 4],
                [0, 4, 3],
                [1, 2, 5],
                [1, 5, 4],
                [0, 2, 7],
                [0, 7, 6],
            ],
            dtype=int,
        )
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _adjacent_boxes_different_heights_mesh():
        tall = trimesh.creation.box(extents=(1.0, 1.0, 10.0))
        tall.apply_translation((0.0, 0.0, 5.0))
        short = trimesh.creation.box(extents=(1.0, 1.0, 5.0))
        short.apply_translation((1.0, 0.0, 2.5))
        return trimesh.util.concatenate([tall, short])

    @staticmethod
    def _adjacent_buildings_problem_shell():
        def make_block(x0, x1, y0, y1, z1, nx, ny, nz, include_bottom=False):
            xs = np.linspace(x0, x1, nx + 1)
            ys = np.linspace(y0, y1, ny + 1)
            zs = np.linspace(0.0, z1, nz + 1)
            vertices = []
            index = {}

            def add_vertex(x, y, z):
                key = (float(x), float(y), float(z))
                if key not in index:
                    index[key] = len(vertices)
                    vertices.append([float(x), float(y), float(z)])
                return index[key]

            faces = []

            def add_rect(p00, p10, p11, p01):
                i00 = add_vertex(*p00)
                i10 = add_vertex(*p10)
                i11 = add_vertex(*p11)
                i01 = add_vertex(*p01)
                faces.append([i00, i10, i11])
                faces.append([i00, i11, i01])

            for i in range(nx):
                for j in range(ny):
                    if include_bottom:
                        add_rect((xs[i], ys[j], 0.0), (xs[i + 1], ys[j], 0.0), (xs[i + 1], ys[j + 1], 0.0), (xs[i], ys[j + 1], 0.0))
                    add_rect((xs[i], ys[j], z1), (xs[i], ys[j + 1], z1), (xs[i + 1], ys[j + 1], z1), (xs[i + 1], ys[j], z1))

            for j in range(ny):
                for k in range(nz):
                    add_rect((x0, ys[j], zs[k]), (x0, ys[j + 1], zs[k]), (x0, ys[j + 1], zs[k + 1]), (x0, ys[j], zs[k + 1]))
                    add_rect((x1, ys[j], zs[k]), (x1, ys[j], zs[k + 1]), (x1, ys[j + 1], zs[k + 1]), (x1, ys[j + 1], zs[k]))

            for i in range(nx):
                for k in range(nz):
                    add_rect((xs[i], y0, zs[k]), (xs[i + 1], y0, zs[k]), (xs[i + 1], y0, zs[k + 1]), (xs[i], y0, zs[k + 1]))
                    add_rect((xs[i], y1, zs[k]), (xs[i], y1, zs[k + 1]), (xs[i + 1], y1, zs[k + 1]), (xs[i + 1], y1, zs[k]))

            return trimesh.Trimesh(vertices=np.asarray(vertices, dtype=float), faces=np.asarray(faces, dtype=int), process=False)

        tall = make_block(60.0, 120.0, 30.0, 150.0, 60.0, nx=1, ny=4, nz=4, include_bottom=False)
        short = make_block(120.0, 180.0, 60.0, 120.0, 30.0, nx=1, ny=2, nz=2, include_bottom=False)
        return trimesh.util.concatenate([tall, short])

    @staticmethod
    def _adjacent_boxes_different_sizes_side_overlap_mesh():
        tall = trimesh.creation.box(extents=(1.0, 2.0, 10.0))
        tall.apply_translation((0.0, 0.0, 5.0))
        short = trimesh.creation.box(extents=(1.0, 1.0, 5.0))
        short.apply_translation((1.0, 0.5, 2.5))
        return trimesh.util.concatenate([tall, short])

    @staticmethod
    def _triangulated_adjacent_wall_mesh():
        def make_block(x0, x1, y0, y1, z1, nx, ny, nz):
            xs = np.linspace(x0, x1, nx + 1)
            ys = np.linspace(y0, y1, ny + 1)
            zs = np.linspace(0.0, z1, nz + 1)
            vertices = []
            index = {}

            def add_vertex(x, y, z):
                key = (float(x), float(y), float(z))
                if key not in index:
                    index[key] = len(vertices)
                    vertices.append([float(x), float(y), float(z)])
                return index[key]

            faces = []

            def add_rect(p00, p10, p11, p01):
                i00 = add_vertex(*p00)
                i10 = add_vertex(*p10)
                i11 = add_vertex(*p11)
                i01 = add_vertex(*p01)
                faces.append([i00, i10, i11])
                faces.append([i00, i11, i01])

            # bottom/top
            for i in range(nx):
                for j in range(ny):
                    add_rect((xs[i], ys[j], 0.0), (xs[i + 1], ys[j], 0.0), (xs[i + 1], ys[j + 1], 0.0), (xs[i], ys[j + 1], 0.0))
                    add_rect((xs[i], ys[j], z1), (xs[i], ys[j + 1], z1), (xs[i + 1], ys[j + 1], z1), (xs[i + 1], ys[j], z1))

            # x-normal walls
            for j in range(ny):
                for k in range(nz):
                    add_rect((x0, ys[j], zs[k]), (x0, ys[j + 1], zs[k]), (x0, ys[j + 1], zs[k + 1]), (x0, ys[j], zs[k + 1]))
                    add_rect((x1, ys[j], zs[k]), (x1, ys[j], zs[k + 1]), (x1, ys[j + 1], zs[k + 1]), (x1, ys[j + 1], zs[k]))

            # y-normal walls
            for i in range(nx):
                for k in range(nz):
                    add_rect((xs[i], y0, zs[k]), (xs[i + 1], y0, zs[k]), (xs[i + 1], y0, zs[k + 1]), (xs[i], y0, zs[k + 1]))
                    add_rect((xs[i], y1, zs[k]), (xs[i], y1, zs[k + 1]), (xs[i + 1], y1, zs[k + 1]), (xs[i + 1], y1, zs[k]))

            return trimesh.Trimesh(vertices=np.asarray(vertices, dtype=float), faces=np.asarray(faces, dtype=int), process=False)

        tall = make_block(-0.5, 0.5, -0.5, 0.5, 10.0, nx=1, ny=3, nz=5)
        short = make_block(0.5, 1.5, -0.5, 0.5, 5.0, nx=1, ny=3, nz=4)
        return trimesh.util.concatenate([tall, short])

    @staticmethod
    def _mesh_with_flat_ground_and_subground_building():
        xs = np.linspace(0.0, 4.0, 5)
        ys = np.linspace(0.0, 4.0, 5)
        ground_vertices = []
        ground_faces = []
        for y in ys:
            for x in xs:
                ground_vertices.append([float(x), float(y), 0.0])
        ground_vertices = np.asarray(ground_vertices, dtype=float)

        def vid(i, j):
            return j * len(xs) + i

        for j in range(len(ys) - 1):
            for i in range(len(xs) - 1):
                v00 = vid(i, j)
                v10 = vid(i + 1, j)
                v11 = vid(i + 1, j + 1)
                v01 = vid(i, j + 1)
                ground_faces.append([v00, v10, v11])
                ground_faces.append([v00, v11, v01])

        ground = trimesh.Trimesh(vertices=ground_vertices, faces=np.asarray(ground_faces, dtype=int), process=False)
        box = trimesh.creation.box(extents=(1.0, 1.0, 3.0))
        box.apply_translation((2.0, 2.0, 0.5))
        return trimesh.util.concatenate([ground, box])

    @staticmethod
    def _mesh_with_flat_ground_and_above_ground_building():
        xs = np.linspace(0.0, 4.0, 5)
        ys = np.linspace(0.0, 4.0, 5)
        ground_vertices = []
        ground_faces = []
        for y in ys:
            for x in xs:
                ground_vertices.append([float(x), float(y), 0.0])
        ground_vertices = np.asarray(ground_vertices, dtype=float)

        def vid(i, j):
            return j * len(xs) + i

        for j in range(len(ys) - 1):
            for i in range(len(xs) - 1):
                v00 = vid(i, j)
                v10 = vid(i + 1, j)
                v11 = vid(i + 1, j + 1)
                v01 = vid(i, j + 1)
                ground_faces.append([v00, v10, v11])
                ground_faces.append([v00, v11, v01])

        ground = trimesh.Trimesh(vertices=ground_vertices, faces=np.asarray(ground_faces, dtype=int), process=False)
        box = trimesh.creation.box(extents=(1.0, 1.0, 3.0))
        box.apply_translation((2.0, 2.0, 6.5))
        return trimesh.util.concatenate([ground, box])

    @staticmethod
    def _mesh_with_sloped_ground_and_building():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.2],
                [0.0, 2.0, 0.2],
                [0.5, 0.5, -0.5],
                [1.0, 0.5, -0.5],
                [1.0, 1.0, 1.5],
                [0.5, 1.0, 1.5],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],
                [4, 5, 6],
                [4, 6, 7],
                [4, 5, 1],
                [4, 1, 0],
            ],
            dtype=int,
        )
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    @staticmethod
    def _count_plane_face_signs(mesh, axis: int, plane_value: float, tol: float = 1.0e-8):
        vertices = np.asarray(mesh.vertices, dtype=float)
        faces = np.asarray(mesh.faces, dtype=int)
        normals = np.asarray(mesh.face_normals, dtype=float)
        plane_mask = np.all(np.isclose(vertices[faces][:, :, axis], plane_value, atol=tol), axis=1)
        pos = int(np.count_nonzero(plane_mask & (normals[:, axis] > 0.5)))
        neg = int(np.count_nonzero(plane_mask & (normals[:, axis] < -0.5)))
        return pos, neg

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

    def test_load_stl_honors_stored_normal_when_winding_disagrees(self):
        with TemporaryDirectory() as tmpdir:
            case_dir = Path(tmpdir)
            stl_path = case_dir / "stored_normal.stl"
            vertices = [
                (0.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (1.0, 0.0, 0.0),
            ]

            with stl_path.open("wb") as f:
                f.write(b"test stored normals".ljust(80, b" "))
                f.write(struct.pack("<I", 1))
                f.write(struct.pack("<3f", 0.0, 0.0, 1.0))
                for vertex in vertices:
                    f.write(struct.pack("<3f", *vertex))
                f.write(struct.pack("<H", 0))

            geom = UDGeom(case_dir)
            geom.load("stored_normal.stl")

            np.testing.assert_allclose(
                geom.face_normals,
                np.array([[0.0, 0.0, 1.0]]),
                atol=1e-12,
            )

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

    def test_face_incenters_fall_back_to_centroids_for_collapsed_faces(self):
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [2.0, 2.0, 3.0],
                [2.0, 2.0, 3.0],
                [2.0, 2.0, 3.0],
            ],
            dtype=float,
        )
        faces = np.array([[0, 1, 2], [3, 4, 5]], dtype=int)
        geom = UDGeom(stl=trimesh.Trimesh(vertices=vertices, faces=faces, process=False))

        with self.assertWarnsRegex(UserWarning, "fully collapsed"):
            incenters = geom.face_incenters

        expected_valid = np.array([1.0, 1.0, 0.0]) / (np.sqrt(2.0) + 2.0)
        np.testing.assert_allclose(incenters[0], expected_valid)
        np.testing.assert_allclose(incenters[1], [2.0, 2.0, 3.0])
        self.assertTrue(np.isfinite(incenters).all())

    def test_check_accepts_valid_box_mesh(self):
        geom = UDGeom(stl=trimesh.creation.box(extents=(1.0, 1.0, 1.0)))
        report = check(geom)

        self.assertTrue(report["valid"])
        self.assertEqual(report["n_unused_vertices"], 0)
        self.assertEqual(report["n_connected_components"], 1)
        self.assertEqual(report["n_independent_surfaces"], 1)

    def test_calculate_independent_surfaces_reports_face_partition(self):
        mesh = self._disconnected_two_triangle_mesh()
        result = calculate_independent_surfaces(mesh)

        self.assertEqual(result["n_surfaces"], 2)
        self.assertEqual(len(result["face_surface_ids"]), len(mesh.faces))
        self.assertEqual(sorted(np.unique(result["face_surface_ids"]).tolist()), [1, 2])
        self.assertEqual(sorted(item["n_faces"] for item in result["surfaces"]), [1, 1])

    def test_udgeom_calculate_independent_surfaces_delegates(self):
        geom = UDGeom(stl=self._disconnected_two_triangle_mesh())
        result = geom.calculate_independent_surfaces()

        self.assertEqual(result["n_surfaces"], 2)
        self.assertEqual(len(result["surfaces"]), 2)

    def test_check_flags_unused_vertices(self):
        mesh = self._mesh_with_unused_vertex()
        report = check(mesh)

        self.assertFalse(report["valid"])
        self.assertEqual(report["n_unused_vertices"], 1)
        self.assertTrue(any("unused vertices" in issue for issue in report["issues"]))

    def test_check_flags_disconnected_components(self):
        mesh = self._disconnected_two_triangle_mesh()
        report = check(mesh)

        self.assertFalse(report["valid"])
        self.assertEqual(report["n_connected_components"], 2)
        self.assertEqual(report["n_independent_surfaces"], 2)
        self.assertTrue(any("disconnected face components" in issue for issue in report["issues"]))

    def test_check_flags_duplicate_faces_as_internal_wall_proxy(self):
        mesh = self._mesh_with_duplicate_face()
        report = check(mesh)

        self.assertFalse(report["valid"])
        self.assertEqual(report["n_duplicate_faces"], 1)
        self.assertTrue(any("duplicate faces" in issue for issue in report["issues"]))

    def test_check_flags_nonfinite_vertices(self):
        mesh = self._mesh_with_nonfinite_vertex()
        report = check(mesh)

        self.assertFalse(report["valid"])
        self.assertEqual(report["n_nonfinite_vertices"], 1)
        self.assertTrue(any("non-finite vertices" in issue for issue in report["issues"]))

    def test_check_flags_zero_area_faces(self):
        mesh = self._mesh_with_zero_area_face()
        report = check(mesh)

        self.assertFalse(report["valid"])
        self.assertEqual(report["n_zero_area_faces"], 1)
        self.assertGreaterEqual(report["n_zero_normals"], 1)
        self.assertTrue(any("zero-area faces" in issue for issue in report["issues"]))

    def test_check_flags_downward_ground_faces(self):
        mesh = self._mesh_with_downward_ground_face()
        report = check(mesh)

        self.assertFalse(report["valid"])
        self.assertEqual(report["n_downward_ground_faces"], 1)
        self.assertTrue(any("downward-facing ground faces" in issue for issue in report["issues"]))

    def test_identify_ground_faces_prefers_lowest_large_horizontal_component(self):
        mesh = self._mesh_with_topography_and_bridge()
        mask = identify_ground_faces(mesh)

        np.testing.assert_array_equal(mask, np.array([True, True, False, False]))

    def test_identify_ground_faces_matches_z0_rule_for_flat_ground_case(self):
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 2.0],
                [1.0, 0.0, 2.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],
                [0, 1, 5],
                [0, 5, 4],
            ],
            dtype=int,
        )
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

        expected = np.all(np.isclose(vertices[faces][:, :, 2], 0.0), axis=1)
        actual = identify_ground_faces(mesh)

        np.testing.assert_array_equal(actual, expected)

    def test_udgeom_identify_ground_faces_delegates(self):
        geom = UDGeom(stl=self._mesh_with_topography_and_bridge())
        mask = geom.identify_ground_faces()

        np.testing.assert_array_equal(mask, np.array([True, True, False, False]))

    def test_find_nonmanifold_regions_reports_cluster(self):
        mesh = self._mesh_with_nonmanifold_edge()
        regions = find_nonmanifold_regions(mesh)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["n_edges"], 1)
        self.assertEqual(regions[0]["n_faces"], 3)
        self.assertEqual(regions[0]["edge_vertex_ids"], [(0, 1)])
        np.testing.assert_allclose(regions[0]["bbox"], np.array([[0.0, -1.0, 0.0], [1.0, 1.0, 1.0]]))

    def test_udgeom_find_nonmanifold_regions_delegates(self):
        geom = UDGeom(stl=self._mesh_with_nonmanifold_edge())
        regions = geom.find_nonmanifold_regions()

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["edge_vertex_ids"], [(0, 1)])

    def test_find_unstitched_touching_regions_reports_shared_geometric_edge(self):
        mesh = self._mesh_with_unstitched_touching_edge()
        regions = find_unstitched_touching_regions(mesh)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["n_edges"], 2)
        self.assertEqual(sorted(regions[0]["face_ids"]), [0, 1])

    def test_find_unstitched_touching_regions_handles_large_coordinates(self):
        mesh = self._mesh_with_unstitched_touching_edge_at_large_coordinates()
        regions = find_unstitched_touching_regions(mesh)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["n_edges"], 2)
        self.assertEqual(len(regions[0]["edge_points"]), 2)
        self.assertEqual(sorted(regions[0]["face_ids"]), [0, 1])

    def test_udgeom_find_unstitched_touching_regions_delegates(self):
        geom = UDGeom(stl=self._mesh_with_unstitched_touching_edge())
        regions = geom.find_unstitched_touching_regions()

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["n_edges"], 2)

    def test_find_internal_touching_wall_regions_reports_adjacent_building_overlap(self):
        mesh = self._adjacent_boxes_different_heights_mesh()
        regions = find_internal_touching_wall_regions(mesh)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["axis"], 0)
        self.assertAlmostEqual(regions[0]["plane_value"], 0.5, places=6)
        self.assertGreater(regions[0]["overlap_area"], 0.0)

    def test_udgeom_find_internal_touching_wall_regions_delegates(self):
        geom = UDGeom(stl=self._adjacent_boxes_different_heights_mesh())
        regions = geom.find_internal_touching_wall_regions()

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["axis"], 0)

    def test_check_includes_diagnostic_details(self):
        mesh = self._mesh_with_nonmanifold_edge()
        report = check(mesh, require_single_component=False)

        self.assertIn("details", report)
        self.assertEqual(len(report["details"]["nonmanifold_regions"]), 1)
        self.assertGreaterEqual(len(report["details"]["connected_components"]), 1)
        self.assertEqual(report["details"]["nonmanifold_regions"][0]["edge_vertex_ids"], [(0, 1)])

    def test_check_flags_internal_touching_walls(self):
        mesh = self._adjacent_boxes_different_heights_mesh()
        report = check(mesh, require_single_component=False)

        self.assertFalse(report["valid"])
        self.assertEqual(report["n_internal_touching_wall_regions"], 1)
        self.assertTrue(any("internal touching wall regions" in issue for issue in report["issues"]))
        self.assertEqual(len(report["details"]["internal_touching_wall_regions"]), 1)

    def test_check_flags_unstitched_touching_regions(self):
        mesh = self._mesh_with_unstitched_touching_edge()
        report = check(mesh, require_single_component=False)

        self.assertFalse(report["valid"])
        self.assertEqual(report["n_unstitched_touching_regions"], 1)
        self.assertTrue(any("unstitched touching regions" in issue for issue in report["issues"]))
        self.assertEqual(len(report["details"]["unstitched_touching_regions"]), 1)

    def test_check_flags_unstitched_ground_building_t_junction(self):
        mesh = self._mesh_with_unstitched_ground_building_t_junction()
        report = check(mesh, require_single_component=False)

        self.assertFalse(report["valid"])
        self.assertEqual(report["n_unstitched_touching_regions"], 1)
        self.assertTrue(any("unstitched touching regions" in issue for issue in report["issues"]))

    def test_udgeom_check_delegates_to_package_validator(self):
        geom = UDGeom(stl=trimesh.creation.box(extents=(1.0, 1.0, 1.0)))
        report = geom.check()

        self.assertTrue(report["valid"])

    def test_fix_removes_duplicate_and_zero_area_faces(self):
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 2],
                [2, 1, 0],
                [0, 1, 3],
            ],
            dtype=int,
        )
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

        fixed, report = fix(mesh)

        self.assertEqual(len(fixed.faces), 1)
        self.assertEqual(report["removed_duplicate_faces"], 1)
        self.assertEqual(report["removed_zero_area_faces"], 1)
        self.assertTrue(report["after"]["valid"])

    def test_fix_can_remove_tiny_disconnected_components(self):
        large = trimesh.creation.box(extents=(1.0, 1.0, 1.0))
        small_vertices = np.array(
            [
                [5.0, 0.0, 0.0],
                [6.0, 0.0, 0.0],
                [5.0, 1.0, 0.0],
            ],
            dtype=float,
        )
        small_faces = np.array([[0, 1, 2]], dtype=int)
        small = trimesh.Trimesh(vertices=small_vertices, faces=small_faces, process=False)
        mesh = trimesh.util.concatenate([large, small])

        fixed, report = fix(mesh, remove_small_components=True, min_component_faces=2, min_component_area_fraction=0.05)

        self.assertEqual(len(fixed.faces), len(large.faces))
        self.assertEqual(report["removed_small_component_faces"], 1)
        self.assertEqual(report["after"]["n_connected_components"], 1)

    def test_udgeom_fix_inplace_updates_mesh(self):
        mesh = self._mesh_with_duplicate_face()
        geom = UDGeom(stl=mesh)

        _, report = geom.fix(inplace=True)

        self.assertEqual(len(geom.stl.faces), 1)
        self.assertEqual(report["removed_duplicate_faces"], 1)

    def test_fix_can_resolve_adjacent_buildings_different_heights_overlap(self):
        mesh = self._adjacent_boxes_different_heights_mesh()
        pos_before, neg_before = self._count_plane_face_signs(mesh, axis=0, plane_value=0.5)
        self.assertGreater(pos_before, 0)
        self.assertGreater(neg_before, 0)

        fixed, report = fix(mesh, resolve_vertical_coplanar_overlaps=True)

        pos_after, neg_after = self._count_plane_face_signs(fixed, axis=0, plane_value=0.5)
        self.assertGreaterEqual(report["resolved_vertical_overlap_regions"], 1)
        self.assertGreater(pos_after, 0)
        self.assertEqual(neg_after, 0)

    def test_fix_can_resolve_partial_side_overlap_for_different_sizes(self):
        mesh = self._adjacent_boxes_different_sizes_side_overlap_mesh()
        pos_before, neg_before = self._count_plane_face_signs(mesh, axis=0, plane_value=0.5)
        self.assertGreater(pos_before, 0)
        self.assertGreater(neg_before, 0)

        fixed, report = fix(mesh, resolve_vertical_coplanar_overlaps=True)
        fixed_check = check(fixed)

        pos_after, neg_after = self._count_plane_face_signs(fixed, axis=0, plane_value=0.5)
        self.assertGreaterEqual(report["resolved_vertical_overlap_regions"], 1)
        self.assertGreater(pos_after, 0)
        self.assertEqual(neg_after, 0)
        self.assertTrue(fixed_check["valid"])
        self.assertEqual(fixed_check["n_nonmanifold_edges"], 0)
        self.assertEqual(fixed_check["n_connected_components"], 1)

    def test_fix_can_resolve_triangulated_internal_wall_overlap(self):
        mesh = self._triangulated_adjacent_wall_mesh()
        pos_before, neg_before = self._count_plane_face_signs(mesh, axis=0, plane_value=0.5)
        self.assertGreater(pos_before, 1)
        self.assertGreater(neg_before, 1)

        fixed, report = fix(mesh, resolve_vertical_coplanar_overlaps=True)
        fixed_check = check(fixed)

        pos_after, neg_after = self._count_plane_face_signs(fixed, axis=0, plane_value=0.5)
        self.assertGreaterEqual(report["resolved_vertical_overlap_regions"], 1)
        self.assertTrue((pos_after == 0 and neg_after > 0) or (neg_after == 0 and pos_after > 0))
        self.assertTrue(fixed_check["valid"])
        self.assertEqual(fixed_check["n_nonmanifold_edges"], 0)
        self.assertEqual(fixed_check["n_connected_components"], 1)

    def test_fix_can_weld_unstitched_ground_building_t_junction(self):
        mesh = self._mesh_with_unstitched_ground_building_t_junction()
        before = check(mesh, require_single_component=False)
        self.assertEqual(before["n_unstitched_touching_regions"], 1)

        fixed, report = fix(mesh, weld_touching_boundaries=True)
        after = check(fixed, require_single_component=False)

        self.assertGreater(report["welded_touching_boundary_faces"], 0)
        self.assertTrue(after["valid"])
        self.assertEqual(after["n_unstitched_touching_regions"], 0)
        self.assertEqual(after["n_nonmanifold_edges"], 0)

    def test_repair_adjacent_buildings_preserves_ground_and_returns_clean_geometry(self):
        shell = self._adjacent_buildings_problem_shell()
        before_geom = add_ground(shell, 240.0, 180.0, edgelength=30.0, preserve_existing_edges=True)
        before_report = check(before_geom, require_single_component=False)

        self.assertEqual(before_report["n_internal_touching_wall_regions"], 1)
        self.assertGreater(before_report["n_nonmanifold_edges"], 0)

        fixed_geom, fix_report = repair_adjacent_buildings(before_geom)
        self.assertIsInstance(fixed_geom, UDGeom)
        fixed_report = check(fixed_geom, require_single_component=False)

        self.assertEqual(fix_report["ground_face_count"], int(np.count_nonzero(identify_ground_faces(before_geom))))
        self.assertEqual(fix_report["resolved_vertical_overlap_regions"], 1)
        self.assertGreaterEqual(fix_report["welded_touching_boundary_faces"], 0)
        self.assertTrue(fixed_report["valid"])
        self.assertEqual(fixed_report["issues"], [])
        self.assertEqual(fixed_report["n_internal_touching_wall_regions"], 0)
        self.assertEqual(fixed_report["n_nonmanifold_edges"], 0)
        self.assertEqual(fixed_report["n_unstitched_touching_regions"], 0)

    def test_udgeom_repair_adjacent_buildings_delegates(self):
        shell = self._adjacent_buildings_problem_shell()
        geom = add_ground(shell, 240.0, 180.0, edgelength=30.0, preserve_existing_edges=True)

        cleaned, report = geom.repair_adjacent_buildings()
        self.assertIsInstance(cleaned, UDGeom)
        cleaned_report = check(cleaned, require_single_component=False)

        self.assertEqual(report["resolved_vertical_overlap_regions"], 1)
        self.assertTrue(cleaned_report["valid"])

    def test_truncate_below_ground_removes_subground_geometry_and_restitches(self):
        mesh = self._mesh_with_flat_ground_and_subground_building()

        cleaned, report = truncate_below_ground(mesh, edgelength=1.0, return_trimesh=True)
        cleaned_check = check(cleaned)
        ground_mask = identify_ground_faces(cleaned)
        centers = np.asarray(cleaned.triangles_center, dtype=float)
        footprint = Polygon([(1.5, 1.5), (2.5, 1.5), (2.5, 2.5), (1.5, 2.5)])
        ground_centers = centers[ground_mask]
        under_building = [
            (float(x), float(y))
            for x, y in ground_centers[:, :2]
            if footprint.contains(Point(float(x), float(y)))
        ]

        self.assertGreater(report["removed_below_ground_faces"] + report["clipped_faces"], 0)
        self.assertGreater(report["preserved_ground_faces"], 0)
        self.assertGreater(report["retriangulated_ground_faces"], 0)
        self.assertTrue(np.all(np.asarray(cleaned.vertices, dtype=float)[:, 2] >= -1.0e-9))
        self.assertTrue(cleaned_check["valid"])
        self.assertEqual(cleaned_check["n_internal_touching_wall_regions"], 0)
        self.assertEqual(cleaned_check["n_unstitched_touching_regions"], 0)
        self.assertEqual(under_building, [])

    def test_check_flags_non_ground_vertices_below_planar_ground(self):
        mesh = self._mesh_with_flat_ground_and_subground_building()

        report = check(mesh, require_single_component=False)

        self.assertFalse(report["valid"])
        self.assertGreater(report["n_below_ground_vertices"], 0)
        self.assertTrue(any("below planar ground" in issue for issue in report["issues"]))

    def test_check_does_not_flag_above_ground_case_as_below_ground(self):
        mesh = self._mesh_with_flat_ground_and_above_ground_building()

        report = check(mesh, require_single_component=False)

        self.assertEqual(report["n_below_ground_vertices"], 0)

    def test_check_summary_reports_multiple_independent_surfaces(self):
        mesh = self._mesh_with_flat_ground_and_above_ground_building()

        report = check(mesh, require_single_component=False)

        self.assertIn("independent surfaces: 2", report["summary"])
        self.assertIn("mesh has 2 independent surfaces", report["summary"])

    def test_truncate_below_ground_rejects_nonplanar_ground(self):
        mesh = self._mesh_with_sloped_ground_and_building()

        with self.assertRaises(NotImplementedError):
            truncate_below_ground(mesh, return_trimesh=True)

    def test_udgeom_truncate_below_ground_delegates(self):
        geom = UDGeom(stl=self._mesh_with_flat_ground_and_subground_building())

        cleaned, report = geom.truncate_below_ground(edgelength=1.0)
        self.assertIsInstance(cleaned, UDGeom)

        self.assertGreater(len(cleaned.stl.faces), 0)
        self.assertGreater(report["ground_face_count"], 0)
        self.assertGreater(report["preserved_ground_faces"], 0)

    def test_truncate_below_ground_is_noop_for_above_ground_case(self):
        mesh = self._mesh_with_flat_ground_and_above_ground_building()

        cleaned, report = truncate_below_ground(mesh)
        self.assertIsInstance(cleaned, UDGeom)
        cleaned_check = check(cleaned, require_single_component=False)

        self.assertEqual(report["removed_below_ground_faces"], 0)
        self.assertEqual(report["clipped_faces"], 0)
        self.assertTrue(cleaned_check["valid"])
        self.assertEqual(cleaned_check["n_unstitched_touching_regions"], 0)

    def test_extrude_to_ground_connects_above_ground_surface(self):
        mesh = self._mesh_with_flat_ground_and_above_ground_building()
        surface_result = calculate_independent_surfaces(mesh)
        centers = np.asarray(mesh.triangles_center, dtype=float)
        building_surface_id = None
        for surface in surface_result["surfaces"]:
            face_ids = np.asarray(surface["face_ids"], dtype=int)
            if float(np.min(centers[face_ids, 2])) > 0.0:
                building_surface_id = int(surface["surface_id"])
                break

        self.assertIsNotNone(building_surface_id)
        cleaned, report = extrude_to_ground(mesh, building_surface_id, return_trimesh=True)
        cleaned_check = check(cleaned, require_single_component=False)

        self.assertGreater(report["moved_vertices"], 0)
        self.assertTrue(cleaned_check["valid"])
        self.assertEqual(cleaned_check["n_unstitched_touching_regions"], 0)
        self.assertEqual(cleaned_check["n_nonmanifold_edges"], 0)
        self.assertEqual(cleaned_check["n_independent_surfaces"], 1)

    def test_udgeom_extrude_to_ground_delegates(self):
        mesh = self._mesh_with_flat_ground_and_above_ground_building()
        geom = UDGeom(stl=mesh)
        surface_result = geom.calculate_independent_surfaces()
        centers = np.asarray(geom.stl.triangles_center, dtype=float)
        building_surface_id = None
        for surface in surface_result["surfaces"]:
            face_ids = np.asarray(surface["face_ids"], dtype=int)
            if float(np.min(centers[face_ids, 2])) > 0.0:
                building_surface_id = int(surface["surface_id"])
                break

        cleaned, report = geom.extrude_to_ground(building_surface_id)
        self.assertIsInstance(cleaned, UDGeom)

        self.assertGreater(len(cleaned.stl.faces), 0)
        self.assertEqual(report["surface_id"], building_surface_id)

    def test_extrude_to_ground_rejects_unknown_surface_id(self):
        mesh = self._mesh_with_flat_ground_and_above_ground_building()

        with self.assertRaises(ValueError):
            extrude_to_ground(mesh, surface_id=999, return_trimesh=True)

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

    def test_show_forwards_to_vis_facade(self):
        class RecordingVis:
            def __init__(self):
                self.calls = []

            def show_geometry(self, **kwargs):
                self.calls.append(kwargs)
                return "show_geometry_result"

        geom = UDGeom.__new__(UDGeom)
        geom.vis = RecordingVis()

        result = UDGeom.show(
            geom,
            color_buildings=False,
            plot_quiver=False,
            normal_scale=0.5,
            show_edges=False,
            show_ground=True,
        )

        self.assertEqual(result, "show_geometry_result")
        self.assertEqual(
            geom.vis.calls,
            [
                {
                    "color_buildings": False,
                    "plot_quiver": False,
                    "normal_scale": 0.5,
                    "show_edges": False,
                    "show_ground": True,
                    "show": True,
                }
            ],
        )

    def test_show_outline_forwards_to_vis_facade(self):
        class RecordingVis:
            def __init__(self):
                self.calls = []

            def show_geometry_outline(self, **kwargs):
                self.calls.append(kwargs)
                return "show_geometry_outline_result"

        geom = UDGeom.__new__(UDGeom)
        geom.vis = RecordingVis()

        result = UDGeom.show_outline(geom, angle_threshold=30.0, show_ground=False)

        self.assertEqual(result, "show_geometry_outline_result")
        self.assertEqual(
            geom.vis.calls,
            [{"angle_threshold": 30.0, "show_ground": False, "show": True}],
        )

    def test_show_outline_forwards_show_flag_to_vis_facade(self):
        class RecordingVis:
            def __init__(self):
                self.calls = []

            def show_geometry_outline(self, **kwargs):
                self.calls.append(kwargs)
                return "show_geometry_outline_result"

        geom = UDGeom.__new__(UDGeom)
        geom.vis = RecordingVis()

        result = UDGeom.show_outline(geom, show=False)

        self.assertEqual(result, "show_geometry_outline_result")
        self.assertEqual(
            geom.vis.calls,
            [{"angle_threshold": 45.0, "show_ground": True, "show": False}],
        )

    def test_show_forwards_show_flag_to_vis_facade(self):
        class RecordingVis:
            def __init__(self):
                self.calls = []

            def show_geometry(self, **kwargs):
                self.calls.append(kwargs)
                return "show_geometry_result"

        geom = UDGeom.__new__(UDGeom)
        geom.vis = RecordingVis()

        # plot_quiver deliberately NOT passed: this pins the default (False).
        result = UDGeom.show(geom, show=False, show_edges=False)

        self.assertEqual(result, "show_geometry_result")
        self.assertEqual(
            geom.vis.calls,
            [{
                "color_buildings": True,
                "plot_quiver": False,
                "normal_scale": 0.2,
                "show_edges": False,
                "show_ground": True,
                "show": False,
            }],
        )

    def test_show_plot_quiver_defaults_false_in_all_signatures(self):
        import inspect

        from udvis.udbase_vis import UDVis

        self.assertIs(
            inspect.signature(UDGeom.show).parameters["plot_quiver"].default, False
        )
        self.assertIs(
            inspect.signature(UDVis.show_geometry).parameters["plot_quiver"].default,
            False,
        )
        self.assertIs(
            inspect.signature(UDVis.show_geometry_pyvista).parameters["plot_quiver"].default,
            False,
        )

    def test_show_routes_to_pyvista_backend(self):
        class RecordingVis:
            def __init__(self):
                self.calls = []

            def show_geometry_pyvista(self, **kwargs):
                self.calls.append(kwargs)
                return "pyvista_result"

        geom = UDGeom.__new__(UDGeom)
        geom.vis = RecordingVis()

        result = UDGeom.show(geom, pyvista=True, show=False)

        self.assertEqual(result, "pyvista_result")
        self.assertEqual(len(geom.vis.calls), 1)
        self.assertEqual(geom.vis.calls[0]["show"], False)

    def test_plot_independent_surfaces_forwards_to_plot_fac(self):
        class RecordingVis:
            def __init__(self):
                self.calls = []

            def plot_independent_surfaces(self, **kwargs):
                self.calls.append(kwargs)
                if kwargs.get("return_result", False):
                    return "surface-figure", {"n_surfaces": 2}
                return "surface-figure"

        geom = UDGeom.__new__(UDGeom)
        geom.vis = RecordingVis()

        fig = UDGeom.plot_independent_surfaces(geom, show=False)

        self.assertEqual(fig, "surface-figure")
        self.assertEqual(geom.vis.calls, [{"show": False, "return_result": False}])

    def test_plot_independent_surfaces_can_return_partition_result(self):
        class RecordingVis:
            def __init__(self):
                self.calls = []

            def plot_independent_surfaces(self, **kwargs):
                self.calls.append(kwargs)
                return "surface-figure", {"n_surfaces": 2}

        geom = UDGeom.__new__(UDGeom)
        geom.vis = RecordingVis()

        fig, result = UDGeom.plot_independent_surfaces(geom, show=False, return_result=True)

        self.assertEqual(fig, "surface-figure")
        self.assertEqual(result["n_surfaces"], 2)
        self.assertEqual(geom.vis.calls, [{"show": False, "return_result": True}])

    def test_add_ground_respects_building_footprint_constraints(self):
        building_mesh = self._open_bottom_box_mesh()
        building_mesh.apply_translation((1.0, 2.0, 0.0))
        mesh = add_ground(building_mesh, 10.0, 10.0, edgelength=2.0, preserve_existing_edges=True).stl

        footprint = Polygon([(1.0, 2.0), (3.0, 2.0), (3.0, 5.0), (1.0, 5.0)])
        crossing = 0
        inside = 0

        vertices = np.asarray(mesh.vertices, dtype=float)
        faces = np.asarray(mesh.faces, dtype=int)
        for face in faces:
            tri = vertices[face]
            if not np.allclose(tri[:, 2], 0.0):
                continue
            poly = Polygon(tri[:, :2])
            overlap = poly.intersection(footprint).area
            if overlap <= 1e-9:
                continue
            if footprint.covers(poly):
                inside += 1
            else:
                crossing += 1

        self.assertEqual(crossing, 0)
        self.assertEqual(inside, 0)

    def test_add_ground_contains_hole_for_open_bottom_box(self):
        building_mesh = self._open_bottom_box_mesh()
        building_mesh.apply_translation((1.0, 2.0, 0.0))
        mesh = add_ground(building_mesh, 10.0, 10.0, edgelength=2.0, preserve_existing_edges=True).stl

        vertices = np.asarray(mesh.vertices, dtype=float)
        faces = np.asarray(mesh.faces, dtype=int)
        ground_mask = np.all(np.isclose(vertices[faces][:, :, 2], 0.0), axis=1)

        polys = []
        for tri in vertices[faces[ground_mask]]:
            poly = Polygon(tri[:, :2])
            if poly.area > 1e-9:
                polys.append(poly)

        union = unary_union(polys)
        self.assertEqual(union.geom_type, "Polygon")
        self.assertEqual(len(union.interiors), 1)

    def test_add_ground_can_return_step1_debug(self):
        geom, debug = add_ground(
            DATA_DIR / "single_box.stl",
            10.0,
            10.0,
            shift=(1.0, 2.0, 0.0),
            edgelength=2.0,
            return_debug=True,
        )

        self.assertIsInstance(geom, UDGeom)
        self.assertIn("outside_region", debug)
        self.assertIn("grid_lines", debug)
        self.assertGreaterEqual(len(debug["grid_lines"]), 1)

    def test_add_ground_with_examples_udales_stl(self):
        geom = add_ground(
            EXAMPLES_DIR / "uDALES.stl",
            256.0,
            128.0,
            shift=(20.0, 0.0, 0.0),
            edgelength=16.0,
        )

        self.assertIsInstance(geom, UDGeom)
        self.assertGreater(geom.n_faces, 0)
        self.assertGreater(geom.n_vertices, 0)

    def test_udales_realistic_ground_identification_survives_sine_cosine_topography(self):
        geom = add_ground(
            EXAMPLES_DIR / "uDALES.stl",
            256.0,
            128.0,
            shift=(20.0, 0.0, 0.0),
            edgelength=16.0,
        )
        base_mesh = geom.stl.copy()
        base_ground_mask = identify_ground_faces(base_mesh)
        self.assertGreater(int(np.count_nonzero(base_ground_mask)), 0)

        base_buildings, _ = geom_split_buildings(base_mesh)
        self.assertGreater(len(base_buildings), 0)

        topo_mesh = base_mesh.copy()
        topo_vertices = np.asarray(topo_mesh.vertices, dtype=float).copy()
        topo_faces = np.asarray(topo_mesh.faces, dtype=int)
        ground_vertex_ids = np.unique(topo_faces[base_ground_mask].reshape(-1))

        x = topo_vertices[ground_vertex_ids, 0]
        y = topo_vertices[ground_vertex_ids, 1]
        topo_vertices[ground_vertex_ids, 2] += (
            1.5 * np.sin(2.0 * np.pi * x / 256.0)
            + 0.75 * np.cos(2.0 * np.pi * y / 128.0)
        )
        topo_mesh.vertices = topo_vertices

        topo_ground_mask = identify_ground_faces(topo_mesh)
        self.assertEqual(int(np.count_nonzero(topo_ground_mask)), int(np.count_nonzero(base_ground_mask)))

        filtered, mapping = delete_ground(topo_mesh)
        self.assertEqual(len(mapping), len(base_mesh.faces) - int(np.count_nonzero(base_ground_mask)))
        self.assertEqual(len(filtered.faces), len(mapping))

        topo_buildings, _ = geom_split_buildings(topo_mesh)
        self.assertEqual(len(topo_buildings), len(base_buildings))

    def test_create_canyons_ground_boundary_edges_are_connected_to_walls(self):
        xsize = 96.0
        ysize = 96.0
        B = 12.0
        W = 12.0
        H = 16.0
        shift = 20.0
        edgelength = 6.0

        geom = create_canyons(xsize, ysize, B, W, H, shift, edgelength, rotate90=False)
        mesh = geom.stl
        vertices = np.asarray(mesh.vertices, dtype=float)
        faces = np.asarray(mesh.faces, dtype=int)

        L = B + W
        Nx = int(round(xsize / (B + W)))
        Ny = int(round(ysize / L))
        x_max = xsize

        footprint_polys = []
        for i in range(Nx):
            for j in range(Ny):
                x0 = i * (B + W) + W / 2.0
                x1 = x0 + B
                y0 = j * L
                y1 = y0 + L
                rect = np.array(
                    [
                        [x0, y0],
                        [x1, y0],
                        [x1, y1],
                        [x0, y1],
                    ],
                    dtype=float,
                )
                mask = (rect[:, 0] > 0.0) & (rect[:, 0] < x_max)
                rect[mask, 0] += shift
                footprint_polys.append(Polygon(rect))
        footprint_boundary = unary_union(footprint_polys).boundary

        edge_to_faces = {}
        for face_id, face in enumerate(faces):
            for a, b in ((0, 1), (1, 2), (2, 0)):
                edge = tuple(sorted((int(face[a]), int(face[b]))))
                edge_to_faces.setdefault(edge, []).append(face_id)

        failures = []
        for edge, face_ids in edge_to_faces.items():
            p0 = vertices[edge[0]]
            p1 = vertices[edge[1]]
            if not (np.isclose(p0[2], 0.0) and np.isclose(p1[2], 0.0)):
                continue

            midpoint = Point(float((p0[0] + p1[0]) / 2.0), float((p0[1] + p1[1]) / 2.0))
            if footprint_boundary.distance(midpoint) > 1e-9:
                continue

            adjacent = vertices[faces[np.asarray(face_ids, dtype=int)]]
            has_ground_face = bool(np.any(np.all(np.isclose(adjacent[:, :, 2], 0.0), axis=1)))
            has_wall_face = bool(np.any(~np.all(np.isclose(adjacent[:, :, 2], 0.0), axis=1)))
            if not (has_ground_face and has_wall_face):
                failures.append((edge, face_ids))

        self.assertEqual(failures, [])

    def test_ground_footprint_union_preserves_internal_holes(self):
        mesh = self._ring_prism_mesh()
        footprint = _ground_footprint_union(mesh)

        self.assertEqual(footprint.geom_type, "Polygon")
        self.assertEqual(len(footprint.interiors), 1)
        self.assertAlmostEqual(footprint.area, 12.0, places=9)

    def test_delete_ground_matches_matlab_contract(self):
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 2.0],
                [1.0, 0.0, 2.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],
                [0, 1, 5],
                [0, 5, 4],
            ],
            dtype=int,
        )
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

        filtered, mapping = delete_ground(mesh)

        self.assertEqual(len(filtered.faces), 2)
        np.testing.assert_array_equal(mapping, np.array([2, 3], dtype=int))
        self.assertEqual(len(filtered.vertices), 4)
        self.assertTrue(np.all(filtered.vertices[:, 2] >= 0.0))

    def test_delete_ground_uses_identified_ground_not_absolute_z_zero(self):
        mesh = self._mesh_with_topography_and_bridge()

        filtered, mapping = delete_ground(mesh)

        np.testing.assert_array_equal(mapping, np.array([2, 3], dtype=int))
        self.assertEqual(len(filtered.faces), 2)
        self.assertTrue(np.allclose(filtered.vertices[:, 2], 8.0))

    def test_delete_ground_matches_old_z0_contract_for_flat_ground_case(self):
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 2.0],
                [1.0, 0.0, 2.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],
                [0, 1, 5],
                [0, 5, 4],
            ],
            dtype=int,
        )
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

        old_ground_mask = np.all(np.isclose(vertices[faces][:, :, 2], 0.0), axis=1)
        expected_mapping = np.where(~old_ground_mask)[0]

        filtered, mapping = delete_ground(mesh)

        np.testing.assert_array_equal(mapping, expected_mapping)
        self.assertEqual(len(filtered.faces), 2)

    def test_split_buildings_uses_delete_ground_mapping(self):
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [3.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [4.0, 1.0, 0.0],
                [3.0, 1.0, 0.0],
                [3.0, 0.0, 1.0],
                [4.0, 0.0, 1.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],
                [0, 1, 5],
                [0, 5, 4],
                [6, 7, 8],
                [6, 8, 9],
                [6, 7, 11],
                [6, 11, 10],
            ],
            dtype=int,
        )
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

        buildings, face_map = geom_split_buildings(mesh)

        self.assertEqual(len(buildings), 2)
        np.testing.assert_array_equal(face_map[:2], np.array([0, 0], dtype=int))
        self.assertTrue(np.all(face_map[2:4] == 1))
        self.assertTrue(np.all(face_map[6:8] == 2))


if __name__ == "__main__":
    unittest.main()
