"""Backend-neutral 3-D scene description and renderers for uDALES plotting.

A :class:`Scene` is a plain description of *what* to draw (meshes, lines,
points, glyphs, a title, an optional colour bar and axis labels). Rendering is
delegated to a backend selected by name via :func:`render_scene`:

* ``"plotly"`` — interactive HTML figure (default; notebook-friendly),
* ``"pyvista"`` — PyVista/VTK window or off-screen ``Plotter``.

Plot methods build one ``Scene`` and hand it to whichever backend the caller
asked for, so a new backend is a new renderer here rather than a parallel set of
``*_pyvista`` methods on every plot function.

Colour handling for a mesh is one of, in priority order: ``scalars`` (mapped
through ``cmap``/``clim``), ``face_colors`` (explicit per-face RGB/RGBA), or
``solid_color``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

RGB = Tuple[float, float, float]

# Shared colour constants (0..1) used by the uDALES geometry views.
GROUND_RGB: RGB = (217 / 255, 217 / 255, 217 / 255)
BUILDING_RGB: RGB = (186 / 255, 212 / 255, 245 / 255)


@dataclass
class MeshPrimitive:
    """A triangulated surface with one of three colouring modes."""

    vertices: np.ndarray
    faces: np.ndarray
    solid_color: Optional[RGB] = None
    face_colors: Optional[np.ndarray] = None  # (M,3|4), 0..1 or 0..255
    scalars: Optional[np.ndarray] = None       # (M,), per face
    cmap: str = "viridis"
    clim: Optional[Tuple[float, float]] = None
    nan_color: RGB = (1.0, 1.0, 1.0)
    show_edges: bool = False
    edge_color: str = "black"
    edge_width: float = 0.5
    opacity: float = 1.0
    name: Optional[str] = None
    # Optional Plotly lighting preset: "flat" (matte, evenly lit) or None.
    lighting: Optional[str] = None

    def resolved_face_colors(self) -> np.ndarray:
        """Return an (M, 4) float RGBA array in 0..1 for every face.

        Resolves whichever colouring mode is set into explicit per-face colours.
        This is the Plotly path; PyVista maps scalars through VTK instead, so the
        two backends produce visually equivalent (not bit-identical) colours.
        """
        import matplotlib.pyplot as plt

        n = len(np.asarray(self.faces))
        if n == 0:
            return np.empty((0, 4), dtype=float)
        if self.scalars is not None:
            scal = np.asarray(self.scalars, dtype=float)
            valid = ~np.isnan(scal)
            out = np.tile(np.array([*self.nan_color, 1.0], dtype=float), (n, 1))
            if np.any(valid):
                if self.clim is not None:
                    vmin, vmax = self.clim
                else:
                    vmin, vmax = float(np.nanmin(scal[valid])), float(np.nanmax(scal[valid]))
                norm = plt.Normalize(vmin=vmin, vmax=vmax)
                out[valid] = plt.get_cmap(self.cmap)(norm(scal[valid]))
            return out
        if self.face_colors is not None:
            fc = np.asarray(self.face_colors, dtype=float)
            if fc.size and fc.max() > 1.0:
                fc = fc / 255.0
            if fc.shape[1] == 3:
                fc = np.hstack([fc, np.ones((len(fc), 1))])
            return fc
        rgb = self.solid_color if self.solid_color is not None else GROUND_RGB
        return np.tile(np.array([*rgb, 1.0], dtype=float), (n, 1))


@dataclass
class LineSet:
    """Line segments defined by index pairs into ``vertices``."""

    vertices: np.ndarray
    segments: np.ndarray  # (S, 2) int
    color: str = "black"
    width: float = 2.0
    name: Optional[str] = None
    # Plotly cosmetic: lift edges that sit on the ground plane (z≈0) slightly so
    # they are not hidden by the ground mesh.
    lift_ground: bool = False


@dataclass
class PointSet:
    points: np.ndarray
    color: str = "blue"
    size: float = 4.0
    opacity: float = 1.0
    name: Optional[str] = None
    # Optional per-point text labels (e.g. surface ids), drawn beside the marker.
    labels: Optional[List[str]] = None


@dataclass
class GlyphSet:
    """Oriented vectors (e.g. face normals), drawn as arrows/lines."""

    points: np.ndarray
    vectors: np.ndarray
    scale: float = 0.2
    color: str = "red"
    name: Optional[str] = None


@dataclass
class ColorBar:
    title: str = ""
    fmt: str = "%.1f"


@dataclass
class Scene:
    meshes: List[MeshPrimitive] = field(default_factory=list)
    lines: List[LineSet] = field(default_factory=list)
    points: List[PointSet] = field(default_factory=list)
    glyphs: List[GlyphSet] = field(default_factory=list)
    title: Optional[str] = None
    colorbar: Optional[ColorBar] = None
    axis_labels: Tuple[str, str, str] = ("x (m)", "y (m)", "z (m)")
    show_axes: bool = True
    bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None  # (mins, maxs)
    # Optional legend as (label, colour-spec) entries, rendered by both backends.
    legend: Optional[List[Tuple[str, str]]] = None
    legend_title: Optional[str] = None

    def compute_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        if self.bounds is not None:
            return self.bounds
        chunks = []
        for m in self.meshes:
            chunks.append(np.asarray(m.vertices, dtype=float))
        for p in self.points:
            chunks.append(np.asarray(p.points, dtype=float))
        for ln in self.lines:
            chunks.append(np.asarray(ln.vertices, dtype=float))
        for g in self.glyphs:
            chunks.append(np.asarray(g.points, dtype=float))
        chunks = [c for c in chunks if len(c)]
        if not chunks:
            return np.zeros(3), np.ones(3)
        allpts = np.vstack(chunks)
        return allpts.min(axis=0), allpts.max(axis=0)


# --------------------------------------------------------------------------
# Dispatch
# --------------------------------------------------------------------------

BACKENDS = ("plotly", "pyvista")
_BACKENDS = BACKENDS  # internal alias retained for existing references

# Single default 3-D rendering backend for the whole package. The constructors
# (UDGeom, UDBase, UDVis) default their ``backend`` argument to this value, so
# changing this one line switches the default everywhere. Individual objects
# (``obj.backend = ...``) and calls (``backend=``) still override it.
#
# PyVista is the default: its server/remote-rendering backend scales to the large
# facet counts uDALES produces (geometry stays server-side, only images stream),
# where Plotly's in-browser Mesh3d struggles. Plotly remains available as an
# optional lightweight backend (``backend="plotly"``).
DEFAULT_BACKEND = "pyvista"


def normalize_backend(backend: Optional[str]) -> str:
    """Return the lower-cased backend name, validating it is supported.

    ``None`` maps to :data:`DEFAULT_BACKEND`. Raises ``ValueError`` for any name
    outside :data:`BACKENDS`, so a bad choice fails at assignment time rather
    than deep inside a renderer.
    """
    name = (backend or DEFAULT_BACKEND).lower()
    if name not in BACKENDS:
        raise ValueError(f"backend must be one of {BACKENDS}; got {backend!r}")
    return name


def render_scene(scene: Scene, backend: str = DEFAULT_BACKEND, show: bool = True):
    """Render ``scene`` with the requested backend.

    Parameters
    ----------
    scene : Scene
        Backend-neutral description of the drawing.
    backend : {"plotly", "pyvista"}, default :data:`DEFAULT_BACKEND`
        Rendering backend.
    show : bool, default True
        Display immediately and return ``None``; otherwise return the
        backend-native object (a Plotly ``Figure`` or a PyVista ``Plotter``).

    Returns
    -------
    plotly.graph_objects.Figure or pyvista.Plotter or None
    """
    backend = normalize_backend(backend)
    if backend == "pyvista":
        return _render_pyvista(scene, show=show)
    return _render_plotly(scene, show=show)


# --------------------------------------------------------------------------
# Plotly backend
# --------------------------------------------------------------------------


_PLOTLY_COLORSCALE = {"viridis": "Viridis", "greys": "Greys", "greys_r": "Greys_r"}


def _render_plotly(scene: Scene, show: bool = True):
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise ImportError("plotly is required for the plotly backend. Install with: pip install plotly") from exc

    # NB: do not mutate pio.renderers.default here — that is process-global and
    # would override a renderer the user configured for scripts/CLI. fig.show()
    # uses Plotly's own configured default (auto-detects notebooks).
    fig = go.Figure()

    for mesh in scene.meshes:
        faces = np.asarray(mesh.faces, dtype=int)
        if len(faces) == 0:
            continue
        verts = np.asarray(mesh.vertices, dtype=float)
        rgba = mesh.resolved_face_colors()
        opacity = float(mesh.opacity)
        trace = go.Mesh3d(
            x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
            i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
            facecolor=[f"rgb({int(round(255 * c[0]))},{int(round(255 * c[1]))},{int(round(255 * c[2]))})" for c in rgba],
            opacity=opacity,
            flatshading=True,
            name=mesh.name,
        )
        if mesh.lighting == "flat":
            trace.update(
                flatshading=False,
                lighting=dict(ambient=1.0, diffuse=0.0, specular=0.0, roughness=1.0, fresnel=0.0),
                lightposition=dict(x=0, y=0, z=1),
            )
        fig.add_trace(trace)

    for ln in scene.lines:
        verts = np.asarray(ln.vertices, dtype=float)
        segs = np.asarray(ln.segments, dtype=int)
        ex, ey, ez = [], [], []
        for a, b in segs:
            p0, p1 = verts[a], verts[b]
            z0, z1 = p0[2], p1[2]
            if ln.lift_ground:
                z0 = p0[2] + 0.1 if abs(p0[2]) < 0.01 else p0[2]
                z1 = p1[2] + 0.1 if abs(p1[2]) < 0.01 else p1[2]
            ex.extend([p0[0], p1[0], None])
            ey.extend([p0[1], p1[1], None])
            ez.extend([z0, z1, None])
        fig.add_trace(go.Scatter3d(
            x=ex, y=ey, z=ez, mode="lines",
            line=dict(color=ln.color, width=ln.width),
            name=ln.name, showlegend=False, hoverinfo="skip",
        ))

    for ps in scene.points:
        pts = np.asarray(ps.points, dtype=float)
        if pts.size == 0:
            continue
        trace = go.Scatter3d(
            x=pts[:, 0], y=pts[:, 1], z=pts[:, 2],
            mode="markers+text" if ps.labels else "markers",
            marker=dict(size=ps.size, color=ps.color, opacity=ps.opacity),
            name=ps.name, showlegend=False,
        )
        if ps.labels:
            trace.update(text=list(ps.labels), textposition="top center",
                         textfont=dict(size=12, color=ps.color))
        fig.add_trace(trace)

    # Legend entries (label, colour), drawn as invisible marker traces so they
    # appear in the Plotly legend; PyVista draws them via add_legend.
    if scene.legend:
        for label, color in scene.legend:
            fig.add_trace(go.Scatter3d(
                x=[None], y=[None], z=[None], mode="markers",
                marker=dict(size=8, color=color),
                name=label, showlegend=True, hoverinfo="skip",
            ))

    # Colour bar for scalar meshes (attached via an invisible marker trace so it
    # coexists with the per-face colouring). Rendered by both backends now.
    if scene.colorbar is not None:
        scalar_meshes = [m for m in scene.meshes
                         if m.scalars is not None and len(np.asarray(m.faces))]
        if scalar_meshes:
            m0 = scalar_meshes[0]
            sc = np.asarray(m0.scalars, dtype=float)
            valid = ~np.isnan(sc)
            if m0.clim is not None:
                cmin, cmax = m0.clim
            elif np.any(valid):
                cmin, cmax = float(np.nanmin(sc[valid])), float(np.nanmax(sc[valid]))
            else:
                cmin, cmax = 0.0, 1.0
            fig.add_trace(go.Scatter3d(
                x=[float(scene.compute_bounds()[0][0])],
                y=[float(scene.compute_bounds()[0][1])],
                z=[float(scene.compute_bounds()[0][2])],
                mode="markers",
                marker=dict(
                    size=0.1, opacity=0.0,
                    color=[cmin], cmin=cmin, cmax=cmax, showscale=True,
                    colorscale=_PLOTLY_COLORSCALE.get(m0.cmap.lower(), "Viridis"),
                    colorbar=dict(title=scene.colorbar.title,
                                  tickformat=scene.colorbar.fmt.lstrip("%"),
                                  len=0.7),
                ),
                hoverinfo="skip", showlegend=False,
            ))

    for g in scene.glyphs:
        base = np.asarray(g.points, dtype=float)
        tip = base + g.scale * np.asarray(g.vectors, dtype=float)
        qx, qy, qz = [], [], []
        for p0, p1 in zip(base, tip):
            qx.extend([p0[0], p1[0], None])
            qy.extend([p0[1], p1[1], None])
            qz.extend([p0[2], p1[2], None])
        fig.add_trace(go.Scatter3d(
            x=qx, y=qy, z=qz, mode="lines",
            line=dict(color=g.color, width=3),
            name=g.name, showlegend=False, hoverinfo="skip",
        ))

    mins, maxs = scene.compute_bounds()
    az, el, dist = np.deg2rad(225.0), np.deg2rad(20.0), 1.75
    lx, ly, lz = scene.axis_labels
    # visible=False hides the axis line, ticks, labels and title together.
    _axis_extra = {} if scene.show_axes else dict(visible=False)
    _legend = dict(title=scene.legend_title, itemsizing="constant") if scene.legend else {}
    fig.update_layout(
        title=scene.title,
        showlegend=bool(scene.legend),
        legend=_legend,
        margin=dict(l=0, r=0, b=0, t=40, pad=0),
        scene=dict(
            aspectmode="data",
            xaxis=dict(title=lx, range=[float(mins[0]), float(maxs[0])], showgrid=False, showbackground=False, **_axis_extra),
            yaxis=dict(title=ly, range=[float(mins[1]), float(maxs[1])], showgrid=False, showbackground=False, **_axis_extra),
            zaxis=dict(title=lz, range=[float(mins[2]), float(maxs[2])], showgrid=False, showbackground=False, **_axis_extra),
            camera=dict(
                projection=dict(type="orthographic"),
                eye=dict(
                    x=float(dist * np.cos(el) * np.cos(az)),
                    y=float(dist * np.cos(el) * np.sin(az)),
                    z=float(dist * np.sin(el)),
                ),
            ),
        ),
    )
    if show:
        fig.show()
        return None
    return fig


# --------------------------------------------------------------------------
# PyVista backend
# --------------------------------------------------------------------------


def _faces_to_pyvista(faces: np.ndarray) -> np.ndarray:
    faces = np.asarray(faces, dtype=np.int64)
    out = np.empty((len(faces), 4), dtype=np.int64)
    out[:, 0] = 3
    out[:, 1:] = faces
    return out.ravel()


def _pyvista_color(c):
    """Normalise a colour spec for PyVista.

    Plotly-style ``"rgb(r,g,b)"`` strings (used throughout the plot code) are
    converted to a 0..1 tuple; hex, named colours and tuples pass through, since
    PyVista accepts those directly.
    """
    if isinstance(c, str) and c.startswith("rgb(") and c.endswith(")"):
        r, g, b = (float(v) for v in c[4:-1].split(","))
        return (r / 255.0, g / 255.0, b / 255.0)
    return c


def _offset_surface_behind_lines(actor) -> None:
    """Apply a VTK polygon offset so a surface mesh renders *behind* coincident
    lines. Outline segments sit exactly on the facet faces; without this they
    lose the depth test and building edges vanish (the Plotly backend renders
    them fine, so this keeps the two backends consistent)."""
    mapper = getattr(actor, "mapper", None)
    if mapper is None:
        return
    try:
        mapper.SetResolveCoincidentTopologyToPolygonOffset()
        mapper.SetRelativeCoincidentTopologyPolygonOffsetParameters(2.0, 2.0)
    except AttributeError:
        pass


def _render_pyvista(scene: Scene, show: bool = True):
    try:
        import pyvista as pv
    except ImportError as exc:
        raise ImportError("pyvista is required for the pyvista backend. Install with: pip install pyvista") from exc

    plotter = pv.Plotter(image_scale=4)
    plotter.set_background("white")
    plotter.enable_parallel_projection()

    for mesh in scene.meshes:
        faces = np.asarray(mesh.faces, dtype=np.int64)
        if len(faces) == 0:
            continue
        verts = np.asarray(mesh.vertices, dtype=float)
        poly = pv.PolyData(verts, _faces_to_pyvista(faces))
        common = dict(show_edges=mesh.show_edges, name=mesh.name or None)
        if mesh.show_edges:
            common.update(edge_color=mesh.edge_color, line_width=mesh.edge_width)
        if mesh.scalars is not None:
            poly.cell_data["scalars"] = np.asarray(mesh.scalars, dtype=float)
            clim = mesh.clim
            if clim is None:
                sc = np.asarray(mesh.scalars, dtype=float)
                valid = ~np.isnan(sc)
                clim = (float(np.nanmin(sc[valid])), float(np.nanmax(sc[valid]))) if np.any(valid) else (0.0, 1.0)
            actor = plotter.add_mesh(poly, scalars="scalars", cmap=mesh.cmap, clim=clim,
                                     nan_color=mesh.nan_color, show_scalar_bar=False, **common)
        elif mesh.face_colors is not None:
            fc = np.asarray(mesh.face_colors, dtype=float)
            if fc.size and fc.max() > 1.0:
                fc = fc / 255.0
            poly.cell_data["rgb"] = fc[:, :3]
            actor = plotter.add_mesh(poly, scalars="rgb", rgb=True, show_scalar_bar=False, **common)
        else:
            actor = plotter.add_mesh(poly, color=_pyvista_color(mesh.solid_color or GROUND_RGB), **common)
        # Push surface polygons back in the depth buffer so the outline lines
        # (which lie exactly on the faces) win the depth test and stay visible,
        # matching the Plotly backend. Without this, building edges are occluded.
        _offset_surface_behind_lines(actor)

    for ln in scene.lines:
        segs = np.asarray(ln.segments, dtype=np.int64)
        if len(segs) == 0:
            continue
        verts = np.asarray(ln.vertices, dtype=float)
        cells = np.empty((len(segs), 3), dtype=np.int64)
        cells[:, 0] = 2
        cells[:, 1:] = segs
        line_poly = pv.PolyData()
        line_poly.points = verts
        line_poly.lines = cells.ravel()
        plotter.add_mesh(line_poly, color=_pyvista_color(ln.color), line_width=ln.width, name=ln.name or None)

    for ps in scene.points:
        pts = np.asarray(ps.points, dtype=float)
        if pts.size == 0:
            continue
        plotter.add_mesh(pv.PolyData(pts),
                         color=_pyvista_color(ps.color), point_size=ps.size, opacity=ps.opacity,
                         render_points_as_spheres=True, name=ps.name or None)
        if ps.labels:
            plotter.add_point_labels(
                pts, list(ps.labels), font_size=12, text_color="black",
                font_family="times", shape=None, fill_shape=False,
                show_points=False, always_visible=True,
                name=(ps.name or "labels") + "-labels",
            )

    for g in scene.glyphs:
        gpts = np.asarray(g.points, dtype=float)
        if gpts.size == 0:
            continue
        arrows = pv.PolyData(gpts)
        arrows["vectors"] = np.asarray(g.vectors, dtype=float)
        glyphs = arrows.glyph(orient="vectors", scale=False, factor=g.scale)
        plotter.add_mesh(glyphs, color=_pyvista_color(g.color), name=g.name or None)

    if scene.title:
        plotter.add_text(scene.title, position="upper_edge", font_size=10, color="black")

    if scene.legend:
        plotter.add_legend(
            [[label, _pyvista_color(color)] for label, color in scene.legend],
            bcolor="white", face=None,
        )

    if scene.colorbar is not None:
        plotter.add_scalar_bar(
            title=scene.colorbar.title, color="black", vertical=True,
            position_x=0.92, position_y=0.15, width=0.04, height=0.7,
            title_font_size=16, label_font_size=14, font_family="times",
            fmt=scene.colorbar.fmt,
        )

    if scene.show_axes:
        mins, maxs = scene.compute_bounds()
        _draw_pyvista_axes(plotter, mins, maxs, scene.axis_labels)

    plotter.view_isometric()
    plotter.camera.azimuth = 180
    plotter.camera.elevation = -10

    if show:
        plotter.show()
        return None
    return plotter


def _draw_pyvista_axes(plotter, mins, maxs, labels=("x (m)", "y (m)", "z (m)")) -> None:
    """Draw x/y/z axes (lines, ticks, labels), working around a VTK StaticTriad
    z-axis rendering bug by drawing plain line meshes plus point labels."""
    import pyvista as pv

    xmin_f, ymin_f, zmin_f = float(mins[0]), float(mins[1]), float(mins[2])
    xmax_f, ymax_f, zmax_f = float(maxs[0]), float(maxs[1]), float(maxs[2])
    span_xy = float(np.asarray(maxs - mins)[:2].max())

    def _draw_axis(p0, p1, ticks, label, tick_dir, font_size=14):
        plotter.add_mesh(pv.Line(p0, p1), color="black", line_width=3)
        tick_len = 0.02 * span_xy
        label_offset = 4.0 * tick_len
        pts, lbls = [], []
        for val, pos in ticks:
            t_end = [pos[j] + tick_dir[j] * tick_len for j in range(3)]
            plotter.add_mesh(pv.Line(pos, t_end), color="black", line_width=2)
            pts.append([pos[j] + tick_dir[j] * label_offset for j in range(3)])
            lbls.append(f"{val:.0f}")
        if pts:
            plotter.add_point_labels(
                np.array(pts), lbls, point_size=0, render_points_as_spheres=False,
                font_size=font_size, text_color="black", font_family="times",
                shape=None, fill_shape=False, show_points=False, margin=3, always_visible=True,
            )
        mid = [0.5 * (p0[j] + p1[j]) for j in range(3)]
        title_pos = [mid[j] + tick_dir[j] * tick_len * 8 for j in range(3)]
        plotter.add_point_labels(
            np.array([title_pos]), [label], point_size=0, render_points_as_spheres=False,
            show_points=False, font_size=font_size + 2, text_color="black", font_family="times",
            shape=None, fill_shape=False, bold=True, margin=3, always_visible=True,
        )

    lx, ly, lz = labels
    x_ticks = [(v, (v, ymin_f, zmin_f)) for v in np.linspace(xmin_f, xmax_f, 9)]
    _draw_axis((xmin_f, ymin_f, zmin_f), (xmax_f, ymin_f, zmin_f), x_ticks, lx, (0, -1, 0))
    y_ticks = [(v, (xmin_f, v, zmin_f)) for v in np.linspace(ymin_f, ymax_f, 9)]
    _draw_axis((xmin_f, ymin_f, zmin_f), (xmin_f, ymax_f, zmin_f), y_ticks, ly, (-1, 0, 0))
    z_ticks = [(v, (xmin_f, ymax_f, v)) for v in [zmin_f, zmax_f]]
    _draw_axis((xmin_f, ymax_f, zmin_f), (xmin_f, ymax_f, zmax_f), z_ticks, lz, (-1, 0, 0))
