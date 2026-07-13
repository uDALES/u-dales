import importlib.util
import io
import sys
import types
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udbase import UDBase  # noqa: E402
from udvis import UDVis  # noqa: E402


class RecordingVis:
    def __init__(self):
        self.calls = []

    def plot_fac(self, **kwargs):
        self.calls.append(("plot_fac", kwargs))
        return "plot_fac_result"

    def plot_2dmap(self, **kwargs):
        self.calls.append(("plot_2dmap", kwargs))
        return "plot_2dmap_result"

    def plot_scalar_source(self, **kwargs):
        self.calls.append(("plot_scalar_source", kwargs))
        return "plot_scalar_source_result"


class TestUDBaseVisualizationCompatibility(unittest.TestCase):
    def _assert_clean_plot_error(self, call, message):
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            result = call()

        text = stderr.getvalue().strip()
        self.assertIsNone(result)
        self.assertEqual(text, f"ERROR: {message}")
        self.assertEqual(len(text.splitlines()), 1)

    def test_udbase_plot_fac_forwards_to_vis_facade(self):
        sim = UDBase.__new__(UDBase)
        sim.vis = RecordingVis()

        result = sim.plot_fac(var="dummy_var", building_ids=[1, 2], show=False)

        self.assertEqual(result, "plot_fac_result")
        self.assertEqual(
            sim.vis.calls,
            [("plot_fac", {"var": "dummy_var", "building_ids": [1, 2], "show": False})],
        )

    def test_udbase_plot_2dmap_forwards_to_vis_facade(self):
        sim = UDBase.__new__(UDBase)
        sim.vis = RecordingVis()

        result = sim.plot_2dmap(val=1.0, labels=["A"], show=False)

        self.assertEqual(result, "plot_2dmap_result")
        self.assertEqual(
            sim.vis.calls,
            [("plot_2dmap", {"val": 1.0, "labels": ["A"], "show": False})],
        )

    def test_udbase_plot_scalar_source_forwards_to_vis_facade(self):
        sim = UDBase.__new__(UDBase)
        sim.vis = RecordingVis()

        result = sim.plot_scalar_source(
            scalar_sources={"point": {}, "line": {}},
            scalar_index=2,
            show=False,
        )

        self.assertEqual(result, "plot_scalar_source_result")
        self.assertEqual(
            sim.vis.calls,
            [
                (
                    "plot_scalar_source",
                    {
                        "scalar_sources": {"point": {}, "line": {}},
                        "scalar_index": 2,
                        "show": False,
                    },
                )
            ],
        )

    def test_udvis_plot_scalar_source_errors_when_source_files_missing(self):
        class DummySim:
            _lfgeom = True
            geom = object()

            def load_scalar_sources(self):
                return {"point": {}, "line": {}}

        vis = UDVis(DummySim())

        self._assert_clean_plot_error(
            lambda: vis.plot_scalar_source(show=False),
            "Scalar source data not found",
        )

    def test_udvis_plot_profiles_errors_when_prof_data_missing(self):
        class DummySim:
            expnr = "001"

            def load_prof(self):
                raise FileNotFoundError("missing")

        vis = UDVis(DummySim())

        self._assert_clean_plot_error(
            lambda: vis.plot_profiles(save=False, show=False),
            "prof.inp.001 data not found",
        )

    def test_udvis_plot_dz_variation_errors_when_prof_data_missing(self):
        class DummySim:
            expnr = "001"

            def load_prof(self):
                raise FileNotFoundError("missing")

        vis = UDVis(DummySim())

        self._assert_clean_plot_error(
            lambda: vis.plot_dz_variation(save=False, show=False),
            "prof.inp.001 data not found",
        )

    def test_udvis_plot_lscale_errors_when_lscale_data_missing(self):
        class DummySim:
            expnr = "001"

            def load_lscale(self):
                raise FileNotFoundError("missing")

        vis = UDVis(DummySim())

        self._assert_clean_plot_error(
            lambda: vis.plot_lscale(save=False, show=False),
            "lscale.inp.001 data not found",
        )

    def test_udvis_plot_solid_errors_when_solid_data_missing(self):
        class DummySim:
            _lfgeom = True
            geom = object()

        vis = UDVis(DummySim())

        self._assert_clean_plot_error(
            lambda: vis.plot_solid("c", show=False),
            "solid_c.txt data not found",
        )

    def test_udvis_plot_fluid_boundary_errors_when_data_missing(self):
        class DummySim:
            _lfgeom = True
            geom = object()

        with TemporaryDirectory() as tmpdir:
            sim = DummySim()
            sim.path = Path(tmpdir)
            vis = UDVis(sim)

            self._assert_clean_plot_error(
                lambda: vis.plot_fluid_boundary("c", show=False),
                "fluid_boundary_c.txt data not found",
            )

    def test_udvis_plot_veg_errors_when_vegetation_data_missing(self):
        class DummySim:
            _lfgeom = True
            geom = object()
            veg = None
            expnr = "001"

            def load_veg(self, cache=True):
                return None

        vis = UDVis(DummySim())

        self._assert_clean_plot_error(
            lambda: vis.plot_veg(show=False),
            "veg.inp.001 data not found",
        )


class TestUDVisRenderingHelpers(unittest.TestCase):
    @staticmethod
    def _make_mesh_data():
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ]
        )
        faces = np.array(
            [
                [0, 1, 2],
                [1, 3, 2],
            ],
            dtype=int,
        )
        return types.SimpleNamespace(vertices=vertices, faces=faces)

    @staticmethod
    def _fake_trimesh_module():
        fake_trimesh = types.ModuleType("trimesh")

        class FakeScene:
            def __init__(self):
                self.geometry = []

            def add_geometry(self, geometry):
                self.geometry.append(geometry)

        class FakeLine:
            def __init__(self, points):
                self.points = points

        class FakePath3D:
            def __init__(self, entities, vertices, colors):
                self.entities = entities
                self.vertices = vertices
                self.colors = colors

        fake_trimesh.Scene = FakeScene
        fake_trimesh.path = types.SimpleNamespace(
            entities=types.SimpleNamespace(Line=FakeLine),
            Path3D=FakePath3D,
        )
        return fake_trimesh

    @staticmethod
    def _fake_ipython_modules():
        fake_ipython = types.ModuleType("IPython")
        fake_display = types.ModuleType("IPython.display")
        fake_display.display = lambda *args, **kwargs: None
        return fake_ipython, fake_display

    def _make_vis(self, face_to_building=None):
        class DummyGeom:
            def __init__(self, mesh, building_map):
                self.stl = mesh
                self._face_to_building = np.asarray(building_map, dtype=int)

            def get_face_to_building_map(self):
                return self._face_to_building

            def _calculate_outline_edges(self, angle_threshold=45.0):
                return [(0, 1), (1, 3)]

        mesh = self._make_mesh_data()
        building_map = [1, 2] if face_to_building is None else face_to_building
        geom = DummyGeom(mesh, building_map)
        return UDVis(geom), mesh

    def test_collect_mesh_edges_returns_unique_sorted_edges(self):
        faces = np.array(
            [
                [0, 1, 2],
                [2, 1, 3],
            ],
            dtype=int,
        )

        edges = UDVis._collect_mesh_edges(faces)

        self.assertEqual(edges, [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)])

    def test_dz_from_zt_reconstructs_grid_spacing(self):
        zt = np.array([1.0, 3.0, 6.0, 10.0])

        dz = UDVis._dz_from_zt(zt)

        np.testing.assert_allclose(dz, [2.0, 2.0, 4.0, 4.0])

    def test_render_scene_routes_to_plotly_with_custom_edges(self):
        vis, mesh = self._make_vis()
        calls = {}

        def render_plotly(meshes, outline_edges, show=True):
            calls["meshes"] = meshes
            calls["outline_edges"] = outline_edges
            calls["show"] = show
            return "plotly-figure"

        vis._render_plotly = render_plotly
        fake_ipython, fake_display = self._fake_ipython_modules()

        with patch.dict(
            sys.modules,
            {
                "trimesh": self._fake_trimesh_module(),
                "IPython": fake_ipython,
                "IPython.display": fake_display,
            },
        ):
            result = vis._render_scene(mesh, custom_edges=[(0, 1), (1, 2)], show=False)

        self.assertEqual(result, "plotly-figure")
        self.assertEqual(calls["meshes"], [mesh])
        self.assertEqual(calls["outline_edges"], [(0, 1), (1, 2)])
        self.assertFalse(calls["show"])

    def test_render_scene_filters_outline_edges_by_building_id(self):
        vis, mesh = self._make_vis(face_to_building=[1, 2])
        calls = {}

        def render_plotly(meshes, outline_edges, show=True):
            calls["outline_edges"] = outline_edges
            return "plotly-figure"

        vis._render_plotly = render_plotly
        fake_ipython, fake_display = self._fake_ipython_modules()

        with patch.dict(
            sys.modules,
            {
                "trimesh": self._fake_trimesh_module(),
                "IPython": fake_ipython,
                "IPython.display": fake_display,
            },
        ):
            result = vis._render_scene(
                mesh,
                custom_edges=[(0, 1), (1, 3), (0, 3)],
                building_ids=np.array([2]),
                show=False,
            )

        self.assertEqual(result, "plotly-figure")
        self.assertEqual(calls["outline_edges"], [(1, 3)])


@unittest.skipIf(
    importlib.util.find_spec("matplotlib") is None,
    "matplotlib is required for dz-variation visualization tests",
)
class TestUDVisDzVariation(unittest.TestCase):
    def test_plot_dz_variation_uses_prof_zt_column(self):
        import matplotlib.pyplot as plt

        class DummySim:
            expnr = "001"
            path = Path(".")

            def load_prof(self):
                return np.array(
                    [
                        [1.0, 300.0],
                        [3.0, 301.0],
                        [6.0, 302.0],
                        [10.0, 303.0],
                    ]
                )

        vis = UDVis(DummySim())

        fig = vis.plot_dz_variation(save=False, show=False)

        ax = fig.axes[0]
        line = ax.lines[0]
        np.testing.assert_allclose(line.get_xdata(), [1, 2, 3, 4])
        np.testing.assert_allclose(line.get_ydata(), [2.0, 2.0, 4.0, 4.0])
        self.assertEqual(ax.get_title(), "dz variation")
        self.assertEqual(ax.get_xlabel(), "$k$")
        self.assertEqual(ax.get_ylabel(), "$dz$ [m]")
        plt.close(fig)


@unittest.skipIf(
    importlib.util.find_spec("matplotlib") is None or importlib.util.find_spec("trimesh") is None,
    "matplotlib and trimesh are required for mesh rendering helper tests",
)
class TestUDVisMeshRendering(unittest.TestCase):
    @staticmethod
    def _make_mesh():
        import trimesh

        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ]
        )
        faces = np.array(
            [
                [0, 1, 2],
                [1, 3, 2],
            ],
            dtype=int,
        )
        return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    def _make_vis(self, face_to_building=None):
        class DummyGeom:
            def __init__(self, mesh, building_map):
                self.stl = mesh
                self._face_to_building = np.asarray(building_map, dtype=int)

            def get_face_to_building_map(self):
                return self._face_to_building

            def _calculate_outline_edges(self, angle_threshold=45.0):
                return [(0, 1), (1, 3)]

        mesh = self._make_mesh()
        building_map = [1, 2] if face_to_building is None else face_to_building
        geom = DummyGeom(mesh, building_map)
        return UDVis(geom), mesh

    def test_create_colored_mesh_filters_faces_by_building_id(self):
        vis, _ = self._make_vis()

        mesh = vis._create_colored_mesh(np.array([10.0, 20.0]), building_ids=np.array([2]))

        np.testing.assert_array_equal(mesh.faces, np.array([[1, 3, 2]]))
        self.assertEqual(mesh.visual.face_colors.shape, (1, 4))


@unittest.skipIf(
    importlib.util.find_spec("plotly") is None or importlib.util.find_spec("trimesh") is None,
    "plotly and trimesh are required for scalar-source visualization tests",
)
class TestUDVisScalarSources(unittest.TestCase):
    def test_plot_scalar_source_uses_one_color_per_scalar_index(self):
        import plotly.graph_objects as go
        import trimesh

        class DummyGeom:
            pass

        class DummySim:
            _lfgeom = True
            lscasrc = True
            lscasrcl = True

            def load_scalar_sources(self):
                return {
                    "point": {
                        1: np.array([[1, 2, 3, 1, 0.5]], dtype=float),
                        2: np.array([[4, 5, 6, 1, 0.5]], dtype=float),
                    },
                    "line": {
                        1: np.array([[0, 0, 0, 1, 1, 1, 1, 0.5]], dtype=float),
                        2: np.array([[2, 2, 2, 3, 3, 3, 1, 0.5]], dtype=float),
                    },
                }

        geom = DummyGeom()
        geom.stl = trimesh.Trimesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        sim = DummySim()
        sim.geom = geom
        vis = UDVis(sim)
        vis._render_scene = lambda *args, **kwargs: go.Figure()

        fig = vis.plot_scalar_source(show=False)

        traces = {trace.name: trace for trace in fig.data}
        self.assertEqual(
            traces["scalar 1 point source"].marker.color,
            traces["scalar 1 line source"].line.color,
        )
        self.assertEqual(
            traces["scalar 2 point source"].marker.color,
            traces["scalar 2 line source"].line.color,
        )
        self.assertNotEqual(
            traces["scalar 1 point source"].marker.color,
            traces["scalar 2 point source"].marker.color,
        )


@unittest.skipUnless(
    importlib.util.find_spec("IPython")
    and importlib.util.find_spec("plotly")
    and importlib.util.find_spec("trimesh"),
    "notebook backend deps required",
)
class TestDisplayOnceSemantics(unittest.TestCase):
    """The plotting entry points must build the COMPLETE figure, display it at
    most once, and return None when they displayed it (so a bare notebook call
    cannot auto-display a duplicate)."""

    @staticmethod
    def _vis_with_box():
        import trimesh

        from udgeom import UDGeom

        geom = UDGeom(stl=trimesh.creation.box(extents=(1.0, 1.0, 1.0)))
        return UDVis(sim=geom)

    def test_show_geometry_show_true_displays_once_and_returns_none(self):
        import plotly.graph_objects as go

        vis = self._vis_with_box()
        with patch.object(go.Figure, "show", autospec=True) as mock_show:
            result = vis.show_geometry(show=True)
        self.assertIsNone(result)
        self.assertEqual(mock_show.call_count, 1)

    def test_show_geometry_show_false_returns_complete_figure_without_display(self):
        import plotly.graph_objects as go

        vis = self._vis_with_box()
        with patch.object(go.Figure, "show", autospec=True) as mock_show:
            fig = vis.show_geometry(show=False, plot_quiver=True)
        self.assertEqual(mock_show.call_count, 0)
        self.assertIsNotNone(fig)
        # Composition must be finished before return: title + quiver present.
        self.assertIn("Geometry:", fig.layout.title.text)
        self.assertIn("normals", [getattr(t, "name", None) for t in fig.data])

    def test_plot_fac_show_true_displays_once_and_returns_none(self):
        import plotly.graph_objects as go

        vis = self._vis_with_box()
        var = np.arange(vis.geom.n_faces, dtype=float)
        with patch.object(go.Figure, "show", autospec=True) as mock_show:
            result = vis.plot_fac(var, show=True)
        self.assertIsNone(result)
        self.assertEqual(mock_show.call_count, 1)

        with patch.object(go.Figure, "show", autospec=True) as mock_show:
            fig = vis.plot_fac(var, show=False)
        self.assertEqual(mock_show.call_count, 0)
        self.assertIsNotNone(fig)

    def test_plot_independent_surfaces_show_true_displays_decorated_figure(self):
        import plotly.graph_objects as go

        vis = self._vis_with_box()
        shown = []
        with patch.object(
            go.Figure, "show", autospec=True,
            side_effect=lambda self_fig, *a, **k: shown.append(self_fig),
        ):
            result = vis.plot_independent_surfaces(show=True)
        self.assertIsNone(result)
        self.assertEqual(len(shown), 1)
        # The DISPLAYED figure must carry the decoration (regression guard:
        # plot_fac(show=True) returning None used to skip it entirely).
        self.assertIn("Independent Surfaces", shown[0].layout.title.text)
        self.assertTrue(
            any(getattr(t, "mode", "") == "markers+text" for t in shown[0].data)
        )

    def test_plot_independent_surfaces_show_false_returns_decorated_figure(self):
        import plotly.graph_objects as go

        vis = self._vis_with_box()
        with patch.object(go.Figure, "show", autospec=True) as mock_show:
            fig, result = vis.plot_independent_surfaces(show=False, return_result=True)
        self.assertEqual(mock_show.call_count, 0)
        self.assertIsNotNone(fig)
        self.assertIn("Independent Surfaces", fig.layout.title.text)
        self.assertIn("n_surfaces", result)


if __name__ == "__main__":
    unittest.main()
