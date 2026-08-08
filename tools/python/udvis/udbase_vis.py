"""Visualization facade for uDALES postprocessing.

Provides the :class:`UDVis` class attached to each :class:`udbase.UDBase`
instance as ``sim.vis``, offering geometry, field, and statistics
plotting methods backed by matplotlib (2-D plots) and plotly or pyvista
(3-D scene plots).
"""
from __future__ import annotations

import logging
import sys
import warnings
from typing import Any, Dict, List, Optional, Union

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as mplPolygon
import numpy as np

from .scene import (
    BUILDING_RGB,
    GROUND_RGB,
    ColorBar,
    GlyphSet,
    LineSet,
    MeshPrimitive,
    PointSet,
    DEFAULT_BACKEND,
    Scene,
    normalize_backend,
    render_scene,
)

logger = logging.getLogger(__name__)


class UDVis:
    """
    Visualization facade for :class:`udbase.UDBase`.

    The owning `UDBase` instance remains the single source of truth for case,
    grid, geometry, and preprocessing state. `UDVis` provides plotting methods
    that operate on that state through `sim.vis.*`, while legacy `sim.plot_*`
    wrappers can forward here for compatibility.

    Rendering backends are intentionally treated as an internal concern of this
    facade. User code should target `UDVis` methods rather than backend-specific
    helper functions so the rendering implementation can evolve without forcing
    API churn in calling code.
    """

    def __init__(self, sim: Any, backend: str = DEFAULT_BACKEND):
        self.sim = None if hasattr(sim, "stl") else sim
        self.geom = sim if hasattr(sim, "stl") else getattr(sim, "geom", None)
        # Default rendering backend for the 3-D scene plots; set once here so it
        # need not be passed to every call. Individual calls may still override
        # it with a ``backend=`` argument.
        self.backend = backend

    @property
    def backend(self) -> str:
        """Default rendering backend for 3-D plots (``"plotly"`` or ``"pyvista"``)."""
        return self._backend

    @backend.setter
    def backend(self, value: str) -> None:
        self._backend = normalize_backend(value)

    def _resolve_backend(self, backend: Optional[str]) -> str:
        return backend if backend is not None else self.backend

    @staticmethod
    def _set_equal_axes_matplotlib(ax, vertices: np.ndarray) -> None:
        """Apply data bounds directly so plots start at the actual ground level."""
        mins = vertices.min(axis=0)
        maxs = vertices.max(axis=0)
        spans = np.maximum(maxs - mins, 1e-9)

        ax.set_xlim(mins[0], maxs[0])
        ax.set_ylim(mins[1], maxs[1])
        ax.set_zlim(mins[2], maxs[2])
        ax.set_box_aspect(spans)

    @staticmethod
    def _collect_mesh_edges(faces: np.ndarray) -> List[tuple]:
        """Return unique mesh edges from triangular face connectivity."""
        edges = set()
        for tri in np.asarray(faces, dtype=int):
            edges.add(tuple(sorted((tri[0], tri[1]))))
            edges.add(tuple(sorted((tri[1], tri[2]))))
            edges.add(tuple(sorted((tri[2], tri[0]))))
        return sorted(edges)

    @staticmethod
    def _missing_plot_data(message: str):
        """Warn that required plot data is unavailable and return ``None``.

        Uses the warnings machinery (not a bare stderr print) so notebooks/scripts
        can capture or silence the diagnostic; the ``None`` return lets plot
        methods bail out gracefully when their inputs are missing.
        """
        frame = sys._getframe(1)
        registry = frame.f_globals.setdefault("__warningregistry__", {})
        original_showwarning = warnings.showwarning

        def show_without_source_line(
            message, category, filename, lineno, file=None, line=None
        ):
            return original_showwarning(message, category, filename, lineno, file=file, line="")

        warnings.showwarning = show_without_source_line
        try:
            warnings.warn_explicit(
                message,
                UserWarning,
                frame.f_code.co_filename,
                frame.f_lineno,
                module=frame.f_globals.get("__name__"),
                registry=registry,
            )
        finally:
            warnings.showwarning = original_showwarning
        return None

    def show_geometry(
        self,
        color_buildings: bool = True,
        plot_quiver: bool = False,
        normal_scale: float = 0.2,
        show_edges: bool = True,
        show_ground: bool = True,
        show: bool = True,
        backend: Optional[str] = None,
    ):
        """
        Visualize a geometry mesh.

        Parameters
        ----------
        color_buildings : bool, default=True
            If True, color buildings and ground differently.
        plot_quiver : bool, default=False
            If True, overlay face-normal vectors.
        normal_scale : float, default=0.2
            Scale factor for normal vectors.
        show_edges : bool, default=True
            If True, draw mesh edges / outlines.
        show_ground : bool, default=True
            If True, include ground faces in the plot.
        show : bool, default=True
            If True, display immediately and return ``None``; otherwise return
            the backend figure/plotter object.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to the backend set on this UDVis
            instance (``sim.vis.backend``, itself set from the constructor).

        Returns
        -------
        plotly.graph_objects.Figure or pyvista.Plotter or None
        """
        if self.geom is None or getattr(self.geom, "stl", None) is None:
            raise ValueError("No geometry loaded. Cannot visualize.")

        stl = self.geom.stl
        vertices = np.asarray(stl.vertices, dtype=float)
        faces = np.asarray(stl.faces, dtype=int)
        face_centers = np.asarray(stl.triangles_center, dtype=float)
        face_normals = np.asarray(stl.face_normals, dtype=float)
        is_building = face_centers[:, 2] > 0

        lighting = "flat" if (color_buildings and show_ground) else None
        meshes = []
        if color_buildings:
            if show_ground and np.any(~is_building):
                meshes.append(MeshPrimitive(vertices, faces[~is_building],
                                            solid_color=GROUND_RGB, name="ground", lighting=lighting))
            if np.any(is_building):
                meshes.append(MeshPrimitive(vertices, faces[is_building],
                                            solid_color=BUILDING_RGB, name="buildings", lighting=lighting))
        else:
            selected_faces = faces if show_ground else faces[is_building]
            meshes.append(MeshPrimitive(vertices, selected_faces, solid_color=GROUND_RGB, name="geometry"))

        scene = Scene(
            meshes=meshes,
            title=f"Geometry: {len(faces)} facets",
            bounds=(vertices.min(axis=0), vertices.max(axis=0)),
        )

        if show_edges:
            edge_faces = faces if show_ground else faces[is_building]
            segments = np.asarray(self._collect_mesh_edges(edge_faces), dtype=int)
            if len(segments):
                scene.lines.append(LineSet(vertices, segments, color="black", width=2, lift_ground=True))

        if plot_quiver:
            scene.glyphs.append(GlyphSet(points=face_centers, vectors=face_normals,
                                         scale=normal_scale, color="red", name="normals"))

        return render_scene(scene, backend=self._resolve_backend(backend), show=show)

    def show_geometry_outline(
        self,
        angle_threshold: float = 45.0,
        show_ground: bool = True,
        color_buildings: bool = False,
        show: bool = True,
        backend: Optional[str] = None,
    ):
        """
        Plot a geometry mesh with detected outline edges highlighted.

        Parameters
        ----------
        angle_threshold : float, default=45.0
            Angle threshold used to detect outline edges.
        show_ground : bool, default=True
            If True, include ground faces in the visualization.
        color_buildings : bool, default=False
            If True, colour buildings blue and ground grey. If False (default),
            use two-tone grey shading (lighter for horizontal faces) that makes
            the outline edges stand out.
        show : bool, default=True
            If True, display immediately and return ``None``; otherwise return
            the backend figure/plotter object.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to the backend set on this UDVis
            instance (``sim.vis.backend``, itself set from the constructor).

        Returns
        -------
        plotly.graph_objects.Figure or pyvista.Plotter or None
        """
        if self.geom is None or getattr(self.geom, "stl", None) is None:
            raise ValueError("No geometry loaded. Cannot visualize.")

        outline_edges = self.geom._calculate_outline_edges(angle_threshold)
        if len(outline_edges) == 0:
            warnings.warn("No outline edges found.")
            return None

        stl = self.geom.stl
        vertices = np.asarray(stl.vertices, dtype=float)
        faces = np.asarray(stl.faces, dtype=int)
        face_centers = np.asarray(stl.triangles_center, dtype=float)
        face_normals = np.asarray(stl.face_normals, dtype=float)
        is_building = face_centers[:, 2] > 0

        meshes = []
        if color_buildings:
            if show_ground and np.any(~is_building):
                meshes.append(MeshPrimitive(vertices, faces[~is_building],
                                            solid_color=GROUND_RGB, name="ground", lighting="flat"))
            if np.any(is_building):
                meshes.append(MeshPrimitive(vertices, faces[is_building],
                                            solid_color=BUILDING_RGB, name="buildings", lighting="flat"))
        else:
            selected_faces = faces if show_ground else faces[is_building]
            selected_normals = face_normals if show_ground else face_normals[is_building]
            horizontal = np.abs(selected_normals[:, 2]) >= np.cos(np.deg2rad(15.0))
            fc = np.tile(np.array([150.0, 150.0, 150.0, 255.0]), (len(selected_faces), 1))
            fc[horizontal] = np.array([217.0, 217.0, 217.0, 255.0])
            meshes.append(MeshPrimitive(vertices, selected_faces, face_colors=fc,
                                        name="geometry", lighting="flat"))

        segments = np.asarray(outline_edges, dtype=int)
        scene = Scene(
            meshes=meshes,
            lines=[LineSet(vertices, segments, color="black", width=2, lift_ground=True)],
            title=f"Geometry Outline ({len(outline_edges)} edges)",
            bounds=(vertices.min(axis=0), vertices.max(axis=0)),
        )
        return render_scene(scene, backend=self._resolve_backend(backend), show=show)

    def _base_overlay_scene(self, geom, title: Optional[str] = None) -> Scene:
        """Scene with the geometry as a grey base mesh plus its outline edges.

        Shared by the point/line overlay plots (vegetation, scalar sources,
        solid points, fluid boundary) so both backends render the same base.
        """
        stl = geom.stl
        vertices = np.asarray(stl.vertices, dtype=float)
        faces = np.asarray(stl.faces, dtype=int)
        scene = Scene(
            meshes=[MeshPrimitive(vertices, faces,
                                  solid_color=(220 / 255, 220 / 255, 220 / 255), name="geometry")],
            title=title,
            bounds=(vertices.min(axis=0), vertices.max(axis=0)),
        )
        segments = np.asarray(self._collect_mesh_edges(faces), dtype=int)
        if len(segments):
            scene.lines.append(LineSet(vertices, segments, color="black", width=1.5, lift_ground=True))
        return scene

    def plot_veg(self, veg: Optional[Dict[str, Any]] = None, show: bool = False,
                 backend: Optional[str] = None):
        """Plot vegetation points on top of the geometry.

        Parameters
        ----------
        veg : dict, optional
            Vegetation data; loaded from the case when omitted.
        show : bool, default=False
            If True, display immediately and return None. If False (default,
            unlike sibling methods), return the figure/plotter.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to this UDVis instance's backend.

        Returns
        -------
        plotly.graph_objects.Figure or pyvista.Plotter or None
        """
        if not self.sim._lfgeom or self.sim.geom is None:
            return self._missing_plot_data("Geometry data not found for plot_veg")
        if veg is None:
            if not hasattr(self.sim, "veg") or self.sim.veg is None:
                try:
                    self.sim.load_veg(cache=True)
                except (OSError, ValueError):
                    return self._missing_plot_data(f"veg.inp.{self.sim.expnr} data not found")
            if not hasattr(self.sim, "veg") or self.sim.veg is None:
                return self._missing_plot_data(f"veg.inp.{self.sim.expnr} data not found")
            veg = self.sim.veg
        points = np.asarray(veg.get("points", []))
        if points.size == 0:
            return self._missing_plot_data(f"veg.inp.{self.sim.expnr} data not found or contains no points")

        max_points = 50000
        if len(points) > max_points:
            rng = np.random.default_rng(0)
            points = points[rng.choice(len(points), size=max_points, replace=False)]
            logger.info("plot_veg: showing %d of %d points", max_points, len(veg['points']))

        xs = self.sim.xt[points[:, 0].astype(int)]
        ys = self.sim.yt[points[:, 1].astype(int)]
        zs = self.sim.zt[points[:, 2].astype(int)]

        scene = self._base_overlay_scene(
            self.sim.geom, title=f"Geometry with Vegetation ({len(points)} points)")
        scene.points.append(PointSet(
            np.column_stack([xs, ys, zs]),
            color="rgb(34,139,34)", size=2, opacity=0.2, name="vegetation"))
        return render_scene(scene, backend=self._resolve_backend(backend), show=show)

    @staticmethod
    def _scalar_source_color(scalar_index: int) -> str:
        palette = [
            "rgb(0,0,139)",
            "rgb(0,100,0)",
            "rgb(139,0,0)",
            "rgb(128,0,128)",
            "rgb(184,134,11)",
            "rgb(0,139,139)",
            "rgb(199,21,133)",
            "rgb(85,85,85)",
        ]
        return palette[(scalar_index - 1) % len(palette)]

    def plot_scalar_source(
        self,
        scalar_sources: Optional[Dict[str, Dict[int, np.ndarray]]] = None,
        scalar_index: Optional[int] = None,
        show: bool = False,
        backend: Optional[str] = None,
    ):
        """Plot scalar point and line sources on top of the geometry.

        Parameters
        ----------
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to this UDVis instance's backend.
        """
        if not self.sim._lfgeom or self.sim.geom is None:
            raise ValueError("Geometry (STL) file required for plot_scalar_source()")

        if scalar_sources is None:
            scalar_sources = self.sim.load_scalar_sources()
        if scalar_sources is None:
            return self._missing_plot_data("Scalar source data not found")

        point_sources = scalar_sources.get("point", {})
        line_sources = scalar_sources.get("line", {})
        if scalar_index is not None:
            scalar_index = int(scalar_index)
            point_sources = {scalar_index: point_sources[scalar_index]} if scalar_index in point_sources else {}
            line_sources = {scalar_index: line_sources[scalar_index]} if scalar_index in line_sources else {}

        source_indices = sorted(set(point_sources) | set(line_sources))
        if not source_indices:
            return self._missing_plot_data("Scalar source data not found")

        scene = self._base_overlay_scene(self.sim.geom)
        n_point = 0
        n_line = 0
        for ii in source_indices:
            color = self._scalar_source_color(int(ii))
            points = np.asarray(point_sources.get(ii, np.empty((0, 5))), dtype=float)
            if points.size:
                points = np.atleast_2d(points)
                n_point += len(points)
                scene.points.append(PointSet(
                    points[:, :3], color=color, size=5, opacity=0.85,
                    name=f"scalar {ii} point source"))

            lines = np.asarray(line_sources.get(ii, np.empty((0, 8))), dtype=float)
            if lines.size:
                lines = np.atleast_2d(lines)
                n_line += len(lines)
                endpoints = np.empty((2 * len(lines), 3), dtype=float)
                endpoints[0::2] = lines[:, 0:3]
                endpoints[1::2] = lines[:, 3:6]
                segments = np.column_stack([np.arange(0, 2 * len(lines), 2),
                                            np.arange(1, 2 * len(lines), 2)])
                scene.lines.append(LineSet(
                    endpoints, segments, color=color, width=6,
                    name=f"scalar {ii} line source"))

        scene.title = f"Geometry with Scalar Sources ({n_point} points, {n_line} lines)"
        return render_scene(scene, backend=self._resolve_backend(backend), show=show)

    def plot_solid(
        self,
        grid_type: str = "c",
        show: bool = False,
        max_points: int = 100_000,
        backend: Optional[str] = None,
    ):
        """Plot IBM solid points for one grid type on top of the geometry.

        Parameters
        ----------
        grid_type : {'u', 'v', 'w', 'c'}, default='c'
            Which staggered grid's solid mask to visualise.
        show : bool, default=False
            Display immediately and return None; otherwise return the figure.
        max_points : int, default=100_000
            Cap on rendered points; a random subsample is drawn when exceeded.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to this UDVis instance's backend.

        Returns
        -------
        plotly.graph_objects.Figure or pyvista.Plotter or None
        """
        grid_type = grid_type.lower()
        if grid_type not in ("u", "v", "w", "c"):
            raise ValueError(f"grid_type must be one of 'u', 'v', 'w', 'c'; got {grid_type!r}")

        if not self.sim._lfgeom or self.sim.geom is None:
            return self._missing_plot_data("Geometry data not found for plot_solid")

        mask = getattr(self.sim, f"S{grid_type}", None)
        if mask is None:
            return self._missing_plot_data(f"solid_{grid_type}.txt data not found")

        # Map 3-D boolean mask → index triples
        ii, jj, kk = np.where(mask)
        n_total = len(ii)
        if n_total == 0:
            return self._missing_plot_data(f"solid_{grid_type}.txt contains no solid points")

        if n_total > max_points:
            rng = np.random.default_rng(0)
            sel = rng.choice(n_total, size=max_points, replace=False)
            ii, jj, kk = ii[sel], jj[sel], kk[sel]
            logger.info("plot_solid: showing %d of %d %s-solid points", max_points, n_total, grid_type)

        # Map indices to physical coordinates using the appropriate staggered grid
        x_arr = {"u": self.sim.xm, "v": self.sim.xt, "w": self.sim.xt, "c": self.sim.xt}[grid_type]
        y_arr = {"u": self.sim.yt, "v": self.sim.ym, "w": self.sim.yt, "c": self.sim.yt}[grid_type]
        z_arr = {"u": self.sim.zt, "v": self.sim.zt, "w": self.sim.zm, "c": self.sim.zt}[grid_type]

        xs = x_arr[ii]
        ys = y_arr[jj]
        zs = z_arr[kk]

        scene = self._base_overlay_scene(
            self.sim.geom, title=f"Solid points — {grid_type}-grid ({len(xs)} of {n_total} shown)")
        scene.points.append(PointSet(
            np.column_stack([xs, ys, zs]), color="rgb(0,0,139)", size=4, opacity=0.3,
            name=f"solid ({grid_type})"))
        return render_scene(scene, backend=self._resolve_backend(backend), show=show)

    def plot_fluid_boundary(
        self,
        grid_type: str = "c",
        show: bool = False,
        max_points: int = 100_000,
        backend: Optional[str] = None,
    ):
        """Plot fluid-boundary points for one grid type on top of the geometry.

        Fluid-boundary points are the fluid cells adjacent to IBM solid surfaces,
        read directly from ``fluid_boundary_<grid>.txt``.

        Parameters
        ----------
        grid_type : {'u', 'v', 'w', 'c'}, default='c'
            Which staggered grid's fluid-boundary points to visualise.
        show : bool, default=False
            Display immediately and return None; otherwise return the figure.
        max_points : int, default=100_000
            Cap on rendered points; a random subsample is drawn when exceeded.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to this UDVis instance's backend.

        Returns
        -------
        plotly.graph_objects.Figure or pyvista.Plotter or None
        """
        grid_type = grid_type.lower()
        if grid_type not in ("u", "v", "w", "c"):
            raise ValueError(f"grid_type must be one of 'u', 'v', 'w', 'c'; got {grid_type!r}")

        if not self.sim._lfgeom or self.sim.geom is None:
            return self._missing_plot_data("Geometry data not found for plot_fluid_boundary")

        # Read fluid_boundary_<grid>.txt directly (1-based indices → 0-based)
        fb_file = self.sim.path / f"fluid_boundary_{grid_type}.txt"
        if not fb_file.exists():
            return self._missing_plot_data(f"fluid_boundary_{grid_type}.txt data not found")

        try:
            locs = self.sim._load_sparse_file(
                fb_file,
                skiprows=1,
                dtype=int,
                min_cols=3,
                zero_based_cols=[0, 1, 2],
            )
        except (OSError, ValueError) as exc:
            return self._missing_plot_data(f"fluid_boundary_{grid_type}.txt could not be loaded: {exc}")

        if locs.size == 0:
            return self._missing_plot_data(f"fluid_boundary_{grid_type}.txt contains no points")

        n_total = len(locs)
        if n_total > max_points:
            rng = np.random.default_rng(0)
            sel = rng.choice(n_total, size=max_points, replace=False)
            locs = locs[sel]
            logger.info("plot_fluid_boundary: showing %d of %d %s-boundary points", max_points, n_total, grid_type)

        ii, jj, kk = locs[:, 0], locs[:, 1], locs[:, 2]

        # Map to physical coordinates on the appropriate staggered grid
        x_arr = {"u": self.sim.xm, "v": self.sim.xt, "w": self.sim.xt, "c": self.sim.xt}[grid_type]
        y_arr = {"u": self.sim.yt, "v": self.sim.ym, "w": self.sim.yt, "c": self.sim.yt}[grid_type]
        z_arr = {"u": self.sim.zt, "v": self.sim.zt, "w": self.sim.zm, "c": self.sim.zt}[grid_type]

        xs = x_arr[ii]
        ys = y_arr[jj]
        zs = z_arr[kk]

        scene = self._base_overlay_scene(
            self.sim.geom,
            title=f"Fluid-boundary points — {grid_type}-grid ({len(xs)} of {n_total} shown)")
        scene.points.append(PointSet(
            np.column_stack([xs, ys, zs]), color="rgb(139,0,0)", size=4, opacity=0.3,
            name=f"fluid boundary ({grid_type})"))
        return render_scene(scene, backend=self._resolve_backend(backend), show=show)

    def _outline_segments(self, geom, building_ids=None, angle_threshold: float = 45.0) -> np.ndarray:
        """Return outline edge segments (S, 2) for ``geom``, optionally filtered
        to edges whose faces belong to ``building_ids``."""
        outline_edges = geom._calculate_outline_edges(angle_threshold=angle_threshold)
        if len(outline_edges) == 0:
            return np.empty((0, 2), dtype=int)
        if building_ids is not None:
            face_to_building = geom.get_face_to_building_map()
            building_ids = np.asarray(building_ids)
            faces = np.asarray(geom.stl.faces, dtype=int)
            kept = []
            for v0, v1 in outline_edges:
                has = np.where(
                    ((faces[:, 0] == v0) | (faces[:, 1] == v0) | (faces[:, 2] == v0))
                    & ((faces[:, 0] == v1) | (faces[:, 1] == v1) | (faces[:, 2] == v1))
                )[0]
                if len(has) > 0 and np.any(np.isin(face_to_building[has], building_ids)):
                    kept.append((v0, v1))
            outline_edges = kept
        return np.asarray(outline_edges, dtype=int).reshape(-1, 2)

    def plot_fac(self, var: np.ndarray, building_ids: Optional[np.ndarray] = None,
                 show: bool = True, backend: Optional[str] = None):
        """Plot facet data as a 3D surface.

        Parameters
        ----------
        var : ndarray
            One value per facet (length must equal the number of faces).
        building_ids : ndarray, optional
            Restrict the plot to these building ids.
        show : bool, default=True
            If True, display immediately and return ``None``; otherwise only
            build and return the figure/plotter.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to the backend set on this UDVis
            instance (``sim.vis.backend``, itself set from the constructor).

        Returns
        -------
        plotly.graph_objects.Figure or pyvista.Plotter or None
        """
        geom = self.geom if self.sim is None else self.sim.geom
        if geom is None or getattr(geom, "stl", None) is None:
            raise ValueError("This method requires a geometry (STL) file.")
        var = np.asarray(var, dtype=float)
        if len(var) != geom.n_faces:
            raise ValueError(
                f"Variable length ({len(var)}) must match number of facets ({geom.n_faces})"
            )

        stl = geom.stl
        vertices = np.asarray(stl.vertices, dtype=float)
        faces = np.asarray(stl.faces, dtype=int)

        selected_faces, selected_var = faces, var
        if building_ids is not None:
            face_to_building = geom.get_face_to_building_map()
            face_mask = np.isin(face_to_building, np.asarray(building_ids))
            if np.any(face_mask):
                selected_faces = faces[face_mask]
                selected_var = var[face_mask]
            else:
                warnings.warn("No valid faces found for the specified building IDs")
                building_ids = None

        valid = ~np.isnan(selected_var)
        clim = (
            (float(np.nanmin(selected_var[valid])), float(np.nanmax(selected_var[valid])))
            if np.any(valid) else None
        )

        scene = Scene(
            meshes=[MeshPrimitive(vertices, selected_faces, scalars=selected_var, cmap="viridis", clim=clim)],
            colorbar=ColorBar(),
            bounds=(vertices.min(axis=0), vertices.max(axis=0)),
        )
        segments = self._outline_segments(geom, building_ids)
        if len(segments):
            scene.lines.append(LineSet(vertices, segments, color="black", width=2, lift_ground=True))

        return render_scene(scene, backend=self._resolve_backend(backend), show=show)

    def plot_independent_surfaces(self, show: bool = True, return_result: bool = False,
                                  backend: Optional[str] = None):
        """
        Color independent face-connected surfaces by surface id.

        Parameters
        ----------
        show : bool, default=True
            If True, display the figure immediately.
        return_result : bool, default=False
            If True, also return the independent-surface partition result.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to this UDVis instance's backend.

        Returns
        -------
        fig or (fig, result)
            The figure/plotter when ``show=False``; ``None`` when ``show=True``
            (displayed instead). When ``return_result=True``, the
            independent-surface partition is returned alongside, i.e.
            ``(fig, result)`` or ``(None, result)``.
        """
        if self.geom is None or getattr(self.geom, "stl", None) is None:
            raise ValueError("No geometry loaded. Cannot visualize.")

        result = self.geom.calculate_independent_surfaces()
        stl = self.geom.stl
        vertices = np.asarray(stl.vertices, dtype=float)
        faces = np.asarray(stl.faces, dtype=int)
        surface_ids = np.asarray(result["face_surface_ids"], dtype=float)
        n_surfaces = max(int(result["n_surfaces"]), 1)

        scene = Scene(
            meshes=[MeshPrimitive(vertices, faces, scalars=surface_ids,
                                  cmap="viridis", clim=(1, n_surfaces))],
            bounds=(vertices.min(axis=0), vertices.max(axis=0)),
            title=f"Independent Surfaces ({result['n_surfaces']})",
            legend_title="Surface IDs",
        )
        segments = self._outline_segments(self.geom)
        if len(segments):
            scene.lines.append(LineSet(vertices, segments, color="black", width=1.5, lift_ground=True))

        # One labelled marker per surface (id at the surface centroid), plus a
        # legend entry carrying the surface's zmin and face count.
        norm = plt.Normalize(vmin=1, vmax=n_surfaces)
        cmap = plt.get_cmap("viridis")
        centers = np.asarray(stl.triangles_center, dtype=float)
        legend = []
        for surface in result["surfaces"]:
            face_ids = np.asarray(surface["face_ids"], dtype=int)
            if len(face_ids) == 0:
                continue
            c = cmap(norm(surface["surface_id"]))
            color = f"rgb({int(round(255 * c[0]))},{int(round(255 * c[1]))},{int(round(255 * c[2]))})"
            centroid = centers[face_ids].mean(axis=0)
            zmin = float(np.min(centers[face_ids, 2]))
            scene.points.append(PointSet(
                centroid.reshape(1, 3), color=color, size=6,
                labels=[str(surface["surface_id"])], name=f"surface-{surface['surface_id']}"))
            legend.append((
                f"Surface {surface['surface_id']} | zmin={zmin:.2f} m | {surface['n_faces']} faces",
                color,
            ))
        scene.legend = legend

        rendered = render_scene(scene, backend=self._resolve_backend(backend), show=show)
        if return_result:
            return rendered, result
        return rendered

    def plot_fac_type(
        self,
        building_ids: Optional[np.ndarray] = None,
        show_outlines: bool = True,
        angle_threshold: float = 45.0,
        show: bool = True,
        backend: Optional[str] = None,
    ):
        """Plot the different surface types in the geometry, coloured by type.

        Parameters
        ----------
        show_outlines : bool, default=True
            Overlay the detected outline edges.
        angle_threshold : float, default=45.0
            Angle threshold for outline-edge detection.
        show : bool, default=True
            Display immediately and return None; otherwise return the
            figure/plotter.
        backend : {"plotly", "pyvista"}, optional
            Rendering backend; defaults to this UDVis instance's backend.

        Returns
        -------
        plotly.graph_objects.Figure or pyvista.Plotter or None
        """
        if self.sim.geom is None:
            raise ValueError(
                "This method requires a geometry (STL) file. Ensure stl_file is specified in namoptions."
            )
        if not hasattr(self.sim, "facs") or self.sim.facs is None:
            raise ValueError(
                f"This method requires facet data. Ensure {self.sim.ffacets}.{self.sim.expnr} exists."
            )
        if not hasattr(self.sim, "factypes") or self.sim.factypes is None:
            raise ValueError(
                f"This method requires facet type data. Ensure {self.sim.ffactypes}.{self.sim.expnr} exists."
            )

        facids = self.sim.facs["typeid"]
        typeids = self.sim.factypes["id"]
        names = self.sim.factypes["name"]
        unique_ids = np.unique(facids)
        default_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

        stl = self.sim.geom.stl
        vertices = np.asarray(stl.vertices, dtype=float)
        faces = np.asarray(stl.faces, dtype=int)

        # One RGBA per face, coloured by its type; plus a (name, colour) legend.
        face_colors = np.full((len(faces), 4), 255.0)
        legend = []
        for idx, type_id in enumerate(unique_ids):
            type_mask = facids == type_id
            hex_color = default_colors[idx % len(default_colors)].lstrip("#")
            rgb = tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
            face_colors[type_mask, :3] = rgb
            face_colors[type_mask, 3] = 230
            name_idx = np.where(typeids == type_id)[0]
            label = str(names[name_idx[0]]) if len(name_idx) > 0 else f"Type {type_id}"
            legend.append((label, f"rgb({rgb[0]},{rgb[1]},{rgb[2]})"))

        scene = Scene(
            meshes=[MeshPrimitive(vertices, faces, face_colors=face_colors, name="facet types")],
            bounds=(vertices.min(axis=0), vertices.max(axis=0)),
            title="Surface Types",
            legend=legend,
            legend_title="Surface types",
        )
        if show_outlines:
            segments = self._outline_segments(self.sim.geom, angle_threshold=angle_threshold)
            if len(segments):
                scene.lines.append(LineSet(vertices, segments, color="black", width=2, lift_ground=True))

        return render_scene(scene, backend=self._resolve_backend(backend), show=show)

    def plot_building_ids(self, show: bool = True):
        """Plot building IDs from above (x,y view) with distinct colors."""
        if self.sim.geom is None:
            raise ValueError(
                "This method requires a geometry (STL) file. Ensure stl_file is specified in namoptions."
            )

        outlines = self.sim.geom.calculate_outline2d()
        if not outlines:
            raise ValueError("No buildings found in geometry")

        num_buildings = len(outlines)
        labels = [str(i) for i in range(1, num_buildings + 1)]
        color_values = np.random.permutation(num_buildings)

        fig, ax = self.plot_2dmap(color_values, labels, show=show)
        pc = ax.collections[0]
        fig.colorbar(pc, ax=ax, cmap="hsv")
        ax.set_title(f"Building Layout with IDs (Total: {num_buildings})")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_aspect("equal")
        if show:
            plt.show()
        return fig, ax

    def plot_2dmap(self, val: Union[float, np.ndarray], labels: Optional[Union[str, list]] = None, show: bool = True):
        """Plot a 2D map of buildings colored by a value per building."""
        if self.sim.geom is None or not hasattr(self.sim.geom, "stl") or self.sim.geom.stl is None:
            raise ValueError("Geometry data not available. Cannot compute outlines.")

        outlines = self.sim.geom.calculate_outline2d()
        if not outlines:
            raise ValueError("No building outlines found in geometry")

        num_buildings = len(outlines)
        if np.isscalar(val):
            values = np.full(num_buildings, float(val))
        else:
            val_array = np.asarray(val)
            if len(val_array) != num_buildings:
                raise ValueError(
                    f"Length of val ({len(val_array)}) must match number of buildings ({num_buildings})"
                )
            values = val_array.astype(float)

        if labels is not None:
            if isinstance(labels, str):
                label_array = [labels] * num_buildings
            else:
                label_array = list(labels)
                if len(label_array) != num_buildings:
                    raise ValueError(
                        f"Number of labels ({len(label_array)}) must match number of buildings ({num_buildings})"
                    )
        else:
            label_array = None

        fig, ax = plt.subplots(figsize=(10, 8))
        patches = []
        colors = []

        for idx, outline in enumerate(outlines):
            if np.isnan(values[idx]):
                continue

            polygon = outline.get("polygon", None)
            centroid = outline.get("centroid", None)
            if polygon is None or len(polygon) == 0:
                continue

            xy = polygon[:, :2]
            patch = mplPolygon(xy, closed=True)
            patches.append(patch)
            colors.append(values[idx])

            if label_array is not None and centroid is not None and not np.any(np.isnan(centroid[:2])):
                ax.text(
                    centroid[0],
                    centroid[1],
                    label_array[idx],
                    ha="center",
                    va="center",
                    fontsize=10,
                    fontweight="bold",
                    color="black",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.7),
                )

        if patches:
            pc = PatchCollection(patches, cmap="viridis", edgecolor="black", linewidth=0.5)
            pc.set_array(np.array(colors))
            ax.add_collection(pc)
            plt.colorbar(pc, ax=ax)

        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_xlim(0, self.sim.xlen)
        ax.set_ylim(0, self.sim.ylen)
        ax.grid(True, alpha=0.3)

        if show:
            plt.show()
        return fig, ax

    def plot_profiles(self, save: bool = True, show: bool = False):
        """Plot initial condition profiles read from prof.inp.<expnr>.

        Reads the file directly from the case directory so it reflects the
        latest written state, independent of any preprocessing arrays.

        Columns in prof.inp: z, thl, qt, u, v, tke

        Parameters
        ----------
        save : bool, default=True
            Save ``profiles.<expnr>.pdf`` to the case directory.
        show : bool, default=False
            Call ``plt.show()`` after plotting.

        Returns
        -------
        matplotlib.figure.Figure
        """
        try:
            pr = self.sim.load_prof()
        except (OSError, ValueError):
            return self._missing_plot_data(f"prof.inp.{self.sim.expnr} data not found")
        if pr is None or np.asarray(pr).size == 0:
            return self._missing_plot_data(f"prof.inp.{self.sim.expnr} data not found")
        pr = np.atleast_2d(np.asarray(pr))
        if pr.shape[1] < 6:
            return self._missing_plot_data(f"prof.inp.{self.sim.expnr} data not found")

        zt = pr[:, 0]

        fig, axes = plt.subplots(1, 4, sharey=True, figsize=(12, 5))

        axes[0].plot(pr[:, 1], zt)
        axes[0].set_title("Temperature")
        axes[0].set_xlabel("thl [K]")
        axes[0].set_ylabel("z [m]")

        axes[1].plot(pr[:, 2], zt)
        axes[1].set_title("Specific humidity")
        axes[1].set_xlabel("qt [kg/kg]")

        axes[2].plot(pr[:, 3], zt, label="u")
        axes[2].plot(pr[:, 4], zt, "r--", label="v")
        axes[2].set_title("Velocity")
        axes[2].set_xlabel("[m/s]")
        axes[2].legend()

        axes[3].plot(pr[:, 5], zt)
        axes[3].set_title("TKE")
        axes[3].set_xlabel("e [m\u00b2/s\u00b2]")

        fig.tight_layout()

        if save:
            from pathlib import Path
            out = Path(self.sim.path) / f"profiles.{self.sim.expnr}.pdf"
            fig.savefig(out)

        if show:
            plt.show()

        return fig

    @staticmethod
    def _dz_from_zt(zt: np.ndarray) -> np.ndarray:
        zt = np.asarray(zt, dtype=float).reshape(-1)
        if zt.size == 0 or not np.all(np.isfinite(zt)):
            raise ValueError("zt data not found")

        z_faces = np.empty(zt.size + 1, dtype=float)
        z_faces[0] = 0.0
        for k, center in enumerate(zt):
            z_faces[k + 1] = 2.0 * center - z_faces[k]

        zm = z_faces[:-1]
        ztop = z_faces[-1]
        return np.diff(np.append(zm, ztop))

    def plot_dz_variation(self, save: bool = True, show: bool = False):
        """Plot vertical grid-spacing variation from prof.inp.<expnr>.

        The first column in ``prof.inp`` is the cell-center grid ``zt``.
        Cell-face locations are reconstructed from those centers using the
        same relation as ``GridSection.generate_zgrid`` and ``dzt`` is then
        calculated as ``diff(append(zm, ztop))``.

        Parameters
        ----------
        save : bool, default=True
            Save ``dz_variation.<expnr>.pdf`` to the case directory.
        show : bool, default=False
            Call ``plt.show()`` after plotting.

        Returns
        -------
        matplotlib.figure.Figure
        """
        try:
            pr = self.sim.load_prof()
        except (OSError, ValueError):
            return self._missing_plot_data(f"prof.inp.{self.sim.expnr} data not found")
        if pr is None or np.asarray(pr).size == 0:
            return self._missing_plot_data(f"prof.inp.{self.sim.expnr} data not found")

        pr = np.atleast_2d(np.asarray(pr))
        if pr.shape[1] < 1:
            return self._missing_plot_data(f"prof.inp.{self.sim.expnr} data not found")

        try:
            dz = self._dz_from_zt(pr[:, 0])
        except ValueError:
            return self._missing_plot_data(f"prof.inp.{self.sim.expnr} data not found")

        k = np.arange(1, len(dz) + 1)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(k, dz)
        ax.set_title("dz variation")
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$dz$ [m]")
        ax.axis("tight")
        fig.tight_layout()

        if save:
            from pathlib import Path
            out = Path(self.sim.path) / f"dz_variation.{self.sim.expnr}.pdf"
            fig.savefig(out)

        if show:
            plt.show()

        return fig

    def plot_lscale(self, save: bool = True, show: bool = False):
        """Plot large-scale forcing profiles read from lscale.inp.<expnr>.

        Reads the file directly from the case directory so it reflects the
        latest written state, independent of any preprocessing arrays.

        Columns in lscale.inp:
          z, uq, vq, pqx, pqy, wfls, dqtdxls, dqtdyls, dqtdtls, dthlrad

        Parameters
        ----------
        save : bool, default=True
            Save ``lscale.<expnr>.pdf`` to the case directory.
        show : bool, default=False
            Call ``plt.show()`` after plotting.

        Returns
        -------
        matplotlib.figure.Figure
        """
        try:
            ls = self.sim.load_lscale()
        except (OSError, ValueError):
            return self._missing_plot_data(f"lscale.inp.{self.sim.expnr} data not found")
        if ls is None or np.asarray(ls).size == 0:
            return self._missing_plot_data(f"lscale.inp.{self.sim.expnr} data not found")
        ls = np.atleast_2d(np.asarray(ls))
        if ls.shape[1] < 10:
            return self._missing_plot_data(f"lscale.inp.{self.sim.expnr} data not found")

        zt = ls[:, 0]

        fig, axes = plt.subplots(1, 6, sharey=True, figsize=(18, 5))

        axes[0].plot(ls[:, 1], zt, label="uq")
        axes[0].plot(ls[:, 2], zt, "r--", label="vq")
        axes[0].set_title("Geostrophic velocity")
        axes[0].set_xlabel("[m/s]")
        axes[0].set_ylabel("z [m]")
        axes[0].legend()

        axes[1].plot(ls[:, 3], zt, label="pqx")
        axes[1].plot(ls[:, 4], zt, "r--", label="pqy")
        axes[1].set_title("Pressure gradient")
        axes[1].set_xlabel("[m/s\u00b2]")
        axes[1].legend()

        axes[2].plot(ls[:, 5], zt)
        axes[2].set_title("Subsidence")
        axes[2].set_xlabel("wfls [m/s]")

        axes[3].plot(ls[:, 6], zt, label="dqtdxls")
        axes[3].plot(ls[:, 7], zt, "r--", label="dqtdyls")
        axes[3].set_title("Moisture advection")
        axes[3].set_xlabel("[kg/kg/m]")
        axes[3].legend()

        axes[4].plot(ls[:, 8], zt)
        axes[4].set_title("Moisture tendency")
        axes[4].set_xlabel("dqtdtls [kg/kg/s]")

        axes[5].plot(ls[:, 9], zt)
        axes[5].set_title("Radiative forcing")
        axes[5].set_xlabel("dthlrad [K/s]")

        fig.tight_layout()

        if save:
            from pathlib import Path
            out = Path(self.sim.path) / f"lscale.{self.sim.expnr}.pdf"
            fig.savefig(out)

        if show:
            plt.show()

        return fig
