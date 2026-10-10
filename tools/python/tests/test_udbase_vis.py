import importlib.util
import os
import sys
import unittest
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

# PyVista *rendering* needs a GL context. It is exercised only when off-screen
# rendering is explicitly enabled (PYVISTA_OFF_SCREEN=1), i.e. in the dedicated
# GL-equipped CI job, so the GL-free matrix gate can install PyVista (the default
# backend) without needing xvfb/mesa. find_spec keeps it skipping when the
# backend isn't installed at all.
_PYVISTA_RENDER_OK = bool(
    importlib.util.find_spec("pyvista")
    and importlib.util.find_spec("trimesh")
    and os.environ.get("PYVISTA_OFF_SCREEN")
)

from _common import PYTHON_DIR

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
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = call()

        self.assertIsNone(result)
        messages = [str(w.message) for w in caught]
        self.assertEqual(messages, [message])  # exactly one, matching diagnostic
        self.assertEqual(len(messages[0].splitlines()), 1)
        formatted = warnings.formatwarning(
            caught[0].message,
            caught[0].category,
            caught[0].filename,
            caught[0].lineno,
            caught[0].line,
        )
        self.assertEqual(len(formatted.splitlines()), 1)

    def test_udbase_plot_fac_forwards_to_vis_facade(self):
        sim = UDBase.__new__(UDBase)
        sim.vis = RecordingVis()

        result = sim.plot_fac(var="dummy_var", building_ids=[1, 2], show=False)

        self.assertEqual(result, "plot_fac_result")
        self.assertEqual(
            sim.vis.calls,
            [("plot_fac", {"var": "dummy_var", "building_ids": [1, 2], "show": False, "backend": None})],
        )

    def test_udbase_plot_trees_forwards_to_plot_veg(self):
        # V1: plot_trees is a backward-compatible alias for plot_veg; it must
        # forward to the (existing) vis.plot_veg, not a nonexistent plot_trees.
        class RecordingVegVis:
            def __init__(self):
                self.calls = []

            def plot_veg(self, **kwargs):
                self.calls.append(("plot_veg", kwargs))
                return "plot_veg_result"

        sim = UDBase.__new__(UDBase)
        sim.vis = RecordingVegVis()

        result = sim.plot_trees(show=False)

        self.assertEqual(result, "plot_veg_result")
        self.assertEqual(sim.vis.calls, [("plot_veg", {"show": False})])

        # regression: the aliases must forward backend/color/opacity kwargs
        # instead of raising TypeError on cases exercising them (issue #337)
        sim.vis.calls.clear()
        sim.plot_trees(show=False, backend="plotly", opacity=0.5)
        self.assertEqual(sim.vis.calls, [
            ("plot_veg", {"show": False, "backend": "plotly", "opacity": 0.5}),
        ])

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
                        "backend": None,
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
        vis = UDVis(sim, backend="plotly")  # asserts on Plotly traces

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

        # These assert Plotly display-once semantics, so pin the Plotly backend
        # (the module default is now pyvista).
        geom = UDGeom(stl=trimesh.creation.box(extents=(1.0, 1.0, 1.0)))
        return UDVis(sim=geom, backend="plotly")

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

class TestBackendSelection(unittest.TestCase):
    """The geometry views expose a single method per view with a ``backend``
    selector, rather than parallel ``*_pyvista`` methods."""

    @staticmethod
    def _vis_with_box():
        import trimesh

        from udgeom import UDGeom

        return UDVis(sim=UDGeom(stl=trimesh.creation.box(extents=(2.0, 2.0, 1.0))))

    def test_render_scene_rejects_unknown_backend(self):
        from udvis.scene import Scene, render_scene

        with self.assertRaises(ValueError):
            render_scene(Scene(), backend="opengl")

    def test_no_parallel_pyvista_methods_remain(self):
        # The refactor collapses the parallel *_pyvista methods into a backend arg.
        for name in ("show_geometry_pyvista", "show_geometry_outline_pyvista", "plot_fac_pyvista"):
            self.assertFalse(hasattr(UDVis, name), f"{name} should have been removed")

    def test_constructor_backend_is_the_default_and_overridable(self):
        import trimesh

        from udgeom import UDGeom

        geom = UDGeom(stl=trimesh.creation.box(extents=(1.0, 1.0, 1.0)), backend="pyvista")
        self.assertEqual(geom.vis.backend, "pyvista")
        # backend=None on a call resolves to the instance default...
        self.assertEqual(geom.vis._resolve_backend(None), "pyvista")
        # ...but an explicit per-call backend overrides it.
        self.assertEqual(geom.vis._resolve_backend("plotly"), "plotly")
        # default when unspecified is the module default (pyvista)
        self.assertEqual(UDVis(sim=geom).backend, "pyvista")

    def test_pyvista_default_backend_routes_without_per_call_arg(self):
        import trimesh

        from udgeom import UDGeom

        geom = UDGeom(stl=trimesh.creation.box(extents=(1.0, 1.0, 1.0)), backend="pyvista")
        calls = {}
        geom.vis.__dict__  # ensure real instance
        import udvis.scene as scene_mod

        orig = scene_mod.render_scene

        def spy(scene, backend="plotly", show=True):
            calls["backend"] = backend
            return "sentinel"

        scene_mod.render_scene = spy
        # udbase_vis imported render_scene by name, patch there too
        import udvis.udbase_vis as vis_mod
        orig_vis = vis_mod.render_scene
        vis_mod.render_scene = spy
        try:
            result = geom.vis.show_geometry(show=False)  # no backend= passed
        finally:
            scene_mod.render_scene = orig
            vis_mod.render_scene = orig_vis
        self.assertEqual(result, "sentinel")
        self.assertEqual(calls["backend"], "pyvista")

    @unittest.skipUnless(
        _PYVISTA_RENDER_OK,
        "pyvista backend + PYVISTA_OFF_SCREEN required for GL rendering tests",
    )
    def test_pyvista_backend_builds_plotters_for_all_views(self):
        import numpy as _np
        import pyvista as pv

        pv.OFF_SCREEN = True
        vis = self._vis_with_box()
        var = _np.arange(vis.geom.n_faces, dtype=float)
        # show=False returns the backend-native object; building it exercises the
        # whole PyVista render path (meshes, glyphs, outline lines, axes, colorbar)
        # short of the GL rasterization that a headless CI cannot do.
        for label, plotter in [
            ("show_geometry", vis.show_geometry(show=False, plot_quiver=True, backend="pyvista")),
            ("show_geometry_outline", vis.show_geometry_outline(show=False, backend="pyvista")),
            ("plot_fac", vis.plot_fac(var, show=False, backend="pyvista")),
            ("plot_independent_surfaces", vis.plot_independent_surfaces(show=False, backend="pyvista")),
        ]:
            self.assertIsInstance(plotter, pv.Plotter, label)
            self.assertGreater(len(plotter.renderer.actors), 0, label)
            plotter.close()

    @unittest.skipUnless(
        importlib.util.find_spec("plotly") and importlib.util.find_spec("trimesh"),
        "plotly backend is optional and not installed",
    )
    def test_independent_surfaces_plotly_has_labels_and_legend(self):
        import trimesh

        from udgeom import UDGeom

        geom = UDGeom(stl=trimesh.creation.box(extents=(2.0, 2.0, 2.0)))
        fig = geom.vis.plot_independent_surfaces(show=False, backend="plotly")
        self.assertIn("Independent Surfaces", fig.layout.title.text)
        self.assertTrue(fig.layout.showlegend)
        self.assertTrue(any(getattr(t, "mode", "") == "markers+text" for t in fig.data))

    @staticmethod
    def _vis_with_stub_sim(tmp_path):
        """A UDVis whose sim exposes the minimal data the overlay plots need."""
        import numpy as _np
        import trimesh
        from types import SimpleNamespace
        from pathlib import Path

        from udgeom import UDGeom

        geom = UDGeom(stl=trimesh.creation.box(extents=(2.0, 2.0, 2.0)))
        nf = geom.n_faces
        grid = _np.linspace(0.0, 2.0, 6)
        (Path(tmp_path) / "fluid_boundary_c.txt").write_text("1\n")
        sim = SimpleNamespace(
            _lfgeom=True, geom=geom, expnr="001", path=Path(tmp_path),
            xt=grid, yt=grid, zt=grid, xm=grid, ym=grid, zm=grid,
            dx=0.4, dy=0.4, dzt=_np.full(6, 0.4),
            ffacets="facets", ffactypes="factypes",
            veg={"points": _np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3]])},
            Sc=_np.zeros((6, 6, 6), dtype=bool),
            facs={"typeid": _np.array([1] * (nf // 2) + [2] * (nf - nf // 2))},
            factypes={"id": _np.array([1, 2]), "name": _np.array(["wall", "roof"])},
            load_scalar_sources=lambda: {
                "point": {1: _np.array([[1, 1, 1, 1, 0.5]], dtype=float)},
                "line": {1: _np.array([[0, 0, 0, 2, 2, 2, 1, 0.5]], dtype=float)},
            },
            _load_sparse_file=lambda *a, **k: _np.array([[0, 0, 0], [1, 1, 1]]),
        )
        sim.Sc[1, 1, 1] = True
        sim.Sc[2, 2, 2] = True
        return UDVis(sim, backend="pyvista")

    @unittest.skipUnless(
        _PYVISTA_RENDER_OK,
        "pyvista backend + PYVISTA_OFF_SCREEN required for GL rendering tests",
    )
    def test_pyvista_backend_builds_all_overlay_plots(self):
        import tempfile

        import pyvista as pv

        pv.OFF_SCREEN = True
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        vis = self._vis_with_stub_sim(tmp.name)
        # These default to the instance backend ("pyvista"); no per-call backend.
        for label, plotter in [
            ("plot_veg", vis.plot_veg(show=False)),
            ("plot_veg_outline", vis.plot_veg_outline(show=False)),
            ("plot_scalar_source", vis.plot_scalar_source(show=False)),
            ("plot_solid", vis.plot_solid("c", show=False)),
            ("plot_fluid_boundary", vis.plot_fluid_boundary("c", show=False)),
            ("plot_fac_type", vis.plot_fac_type(show=False)),
        ]:
            self.assertIsInstance(plotter, pv.Plotter, label)
            self.assertGreater(len(plotter.renderer.actors), 0, label)
            plotter.close()

class TestOptionalBackendDependencyErrors(unittest.TestCase):
    """A missing backend dependency must fail with an actionable message naming
    the feature and the install command (not a bare ImportError)."""

    def _render_with_missing(self, backend, missing_modules):
        from udvis.scene import Scene, render_scene
        from exceptions import DependencyError

        blocked = {name: None for name in missing_modules}
        with patch.dict(sys.modules, blocked):
            # DependencyError so callers can catch it typed; also an ImportError.
            with self.assertRaises(DependencyError) as ctx:
                render_scene(Scene(), backend=backend, show=False)
        return str(ctx.exception)

    def test_plotly_backend_missing_points_at_requirements_file(self):
        msg = self._render_with_missing("plotly", ["plotly", "plotly.graph_objects"])
        self.assertIn("plotly", msg.lower())
        self.assertIn("requirements-plotly.txt", msg)

    def test_pyvista_backend_missing_names_install_command(self):
        msg = self._render_with_missing("pyvista", ["pyvista"])
        self.assertIn("pyvista", msg.lower())
        self.assertIn("pip install", msg)


class TestNiceTicks(unittest.TestCase):
    """Tick locator used by the hand-drawn PyVista axes (issue #337)."""

    def test_round_steps_for_typical_domains(self):
        from udvis.scene import _nice_ticks

        values, decimals = _nice_ticks(0.0, 340.0)
        self.assertEqual(values, [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0])
        self.assertEqual(decimals, 0)

        values, decimals = _nice_ticks(0.0, 240.0)
        self.assertEqual(values, [0.0, 50.0, 100.0, 150.0, 200.0])
        self.assertEqual(decimals, 0)

    def test_z_range_gets_intermediate_ticks(self):
        from udvis.scene import _nice_ticks

        values, decimals = _nice_ticks(0.0, 104.0)
        self.assertEqual(values, [0.0, 20.0, 40.0, 60.0, 80.0, 100.0])
        self.assertEqual(decimals, 0)

    def test_at_most_max_ticks_values(self):
        from udvis.scene import _nice_ticks

        for vmax in (1.0, 7.3, 55.0, 104.0, 240.0, 333.0, 1024.0):
            values, _ = _nice_ticks(0.0, vmax, max_ticks=8)
            self.assertLessEqual(len(values), 8)
            self.assertGreaterEqual(len(values), 4)
            self.assertTrue(all(0.0 <= v <= vmax for v in values))

    def test_fractional_steps_report_decimals(self):
        from udvis.scene import _nice_ticks

        values, decimals = _nice_ticks(0.0, 1.0)
        self.assertEqual(values, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        self.assertEqual(decimals, 1)

    def test_offset_range_aligns_to_step_multiples(self):
        from udvis.scene import _nice_ticks

        values, _ = _nice_ticks(-104.0, 104.0)
        self.assertIn(0.0, values)
        self.assertEqual(values[0], -100.0)
        self.assertEqual(values[-1], 100.0)

    def test_degenerate_range_returns_single_tick(self):
        from udvis.scene import _nice_ticks

        values, decimals = _nice_ticks(5.0, 5.0)
        self.assertEqual(values, [5.0])
        self.assertEqual(decimals, 0)


class TestPlotlySceneAxes(unittest.TestCase):
    """Plotly 3-D axes must draw visible lines, clear of the data (issue #337)."""

    @unittest.skipUnless(
        importlib.util.find_spec("plotly"), "plotly backend is optional and not installed"
    )
    def test_axis_lines_shown_on_bounds_with_label_standoff(self):
        from udvis.scene import MeshPrimitive, Scene, render_scene

        verts = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 20.0, 5.0],
                          [0.0, 20.0, 5.0]])
        faces = np.array([[0, 1, 2], [0, 2, 3]])
        scene = Scene(meshes=[MeshPrimitive(verts, faces)],
                      bounds=(verts.min(axis=0), verts.max(axis=0)))
        fig = render_scene(scene, backend="plotly", show=False)

        for name, hi in (("xaxis", 10.0), ("yaxis", 20.0), ("zaxis", 5.0)):
            axis = getattr(fig.layout.scene, name)
            self.assertTrue(axis.showline, name)
            # axes anchored exactly on the domain bounds (no padding); the
            # long outside ticks provide the label standoff instead
            self.assertEqual(axis.range[0], 0.0, name)
            self.assertEqual(axis.range[1], hi, name)
            self.assertEqual(axis.ticks, "outside", name)
            self.assertGreaterEqual(axis.ticklen, 8, name)


class TestVegVoxelMesh(unittest.TestCase):
    """plot_veg draws each vegetation point as a box filling its grid cell."""

    def test_single_cell_spans_exact_cell_bounds(self):
        from udvis.udbase_vis import _veg_voxel_mesh

        xm = np.array([0.0, 2.0, 4.0])
        ym = np.array([0.0, 3.0])
        zm = np.array([0.0, 1.0, 2.5])
        dzt = np.array([1.0, 1.5, 2.0])
        verts, faces = _veg_voxel_mesh([[1, 0, 2]], xm, ym, zm, dx=2.0, dy=3.0, dzt=dzt)
        self.assertEqual(verts.shape, (8, 3))
        self.assertEqual(faces.shape, (12, 3))
        np.testing.assert_allclose(verts.min(axis=0), [2.0, 0.0, 2.5])
        np.testing.assert_allclose(verts.max(axis=0), [4.0, 3.0, 4.5])  # z1 = zm[2] + dzt[2]

    def test_multiple_cells_index_into_own_vertex_block(self):
        from udvis.udbase_vis import _veg_voxel_mesh

        grid = np.arange(5, dtype=float)
        verts, faces = _veg_voxel_mesh(
            [[0, 0, 0], [3, 2, 1]], grid, grid, grid, dx=1.0, dy=1.0,
            dzt=np.ones(5))
        self.assertEqual(verts.shape, (16, 3))
        self.assertEqual(faces.shape, (24, 3))
        # second box's triangles reference only its own 8 vertices
        self.assertTrue((faces[12:] >= 8).all() and (faces[12:] < 16).all())

    def test_boxes_are_closed_watertight_triangulations(self):
        import trimesh

        from udvis.udbase_vis import _veg_voxel_mesh

        grid = np.arange(3, dtype=float)
        verts, faces = _veg_voxel_mesh([[1, 1, 1]], grid, grid, grid,
                                       dx=1.0, dy=1.0, dzt=np.ones(3))
        box = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
        self.assertTrue(box.is_watertight)
        self.assertGreater(box.volume, 0)  # outward winding

    def test_empty_points_give_empty_mesh(self):
        from udvis.udbase_vis import _veg_voxel_mesh

        verts, faces = _veg_voxel_mesh(
            np.empty((0, 3)), np.array([0.0]), np.array([0.0]), np.array([0.0]),
            dx=1.0, dy=1.0, dzt=np.array([1.0]))
        self.assertEqual(len(verts), 0)
        self.assertEqual(len(faces), 0)


class TestZTickValues(unittest.TestCase):
    """Vertical-axis ticks skip the base label that collides with the y axis."""

    def test_base_tick_dropped_at_exact_zero(self):
        from udvis.scene import _z_tick_values

        values, _ = _z_tick_values(0.0, 104.0)
        self.assertEqual(values, [20.0, 40.0, 60.0, 80.0, 100.0])

    def test_base_tick_dropped_despite_float_noise_below_zero(self):
        # STL meshes commonly have zmin at -1e-9 rather than exactly 0; the
        # base tick must still be excluded (regression for issue #337).
        from udvis.scene import _z_tick_values

        values, _ = _z_tick_values(-1.28e-9, 104.42)
        self.assertEqual(values, [20.0, 40.0, 60.0, 80.0, 100.0])

    def test_first_tick_kept_for_elevated_ranges(self):
        from udvis.scene import _z_tick_values

        values, _ = _z_tick_values(50.0, 104.0)
        self.assertEqual(values[0], 60.0)

    def test_degenerate_range_yields_no_ticks(self):
        from udvis.scene import _z_tick_values

        values, _ = _z_tick_values(0.0, 0.0)
        self.assertEqual(values, [])


if __name__ == "__main__":
    unittest.main()
