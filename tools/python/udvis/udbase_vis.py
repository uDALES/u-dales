from __future__ import annotations

import sys
from typing import Any, Dict, List, Optional, Union

import numpy as np


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

    def __init__(self, sim: Any):
        self.sim = sim

    def plot_veg(self, veg: Optional[Dict[str, Any]] = None, show: bool = False):
        """Plot vegetation points on top of the geometry."""
        if not self.sim._lfgeom or self.sim.geom is None:
            raise ValueError("Geometry (STL) file required for plot_veg()")
        if veg is None:
            if not hasattr(self.sim, "veg") or self.sim.veg is None:
                self.sim.load_veg(cache=True)
            if not hasattr(self.sim, "veg") or self.sim.veg is None:
                raise ValueError("veg.inp file required for plot_veg()")
            veg = self.sim.veg
        points = np.asarray(veg.get("points", []))
        if points.size == 0:
            raise ValueError("veg.inp contains no vegetation points")

        max_points = 50000
        if len(points) > max_points:
            rng = np.random.default_rng(0)
            points = points[rng.choice(len(points), size=max_points, replace=False)]
            print(f"plot_veg: showing {max_points} of {len(veg['points'])} points")

        xs = self.sim.xt[points[:, 0].astype(int)]
        ys = self.sim.yt[points[:, 1].astype(int)]
        zs = self.sim.zt[points[:, 2].astype(int)]

        try:
            import plotly.graph_objects as go
        except ImportError as exc:
            raise ImportError("plotly is required for plot_veg. Install with: pip install plotly") from exc

        try:
            import trimesh
        except ImportError as exc:
            raise ImportError("trimesh is required for plot_veg. Install with: pip install trimesh") from exc

        base_mesh = self.sim.geom.stl.copy()
        base_color = np.array([220, 220, 220, 255], dtype=np.uint8)
        base_mesh.visual.face_colors = np.tile(base_color, (len(base_mesh.faces), 1))

        faces = self.sim.geom.stl.faces
        edges = set()
        for tri in faces:
            e0 = tuple(sorted((tri[0], tri[1])))
            e1 = tuple(sorted((tri[1], tri[2])))
            e2 = tuple(sorted((tri[2], tri[0])))
            edges.update([e0, e1, e2])
        outline_edges = list(edges)

        fig = self._render_scene(
            [base_mesh],
            show_outlines=True,
            custom_edges=outline_edges,
            show=False,
        )
        if fig is None:
            return None
        veg_trace = go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="markers",
            marker=dict(size=2, color="rgb(34,139,34)", opacity=0.2),
            name="vegetation",
        )

        fig.add_trace(veg_trace)
        fig.update_layout(title=f"Geometry with Vegetation ({len(points)} points)")
        if show:
            fig.show()
        return fig

    def plot_trees(self, show: bool = False):
        """Backward-compatible alias for plot_veg."""
        return self.plot_veg(show=show)

    def plot_fac(self, var: np.ndarray, building_ids: Optional[np.ndarray] = None, show: bool = True):
        """Plot facet data as a 3D surface."""
        if self.sim.geom is None:
            raise ValueError("This method requires a geometry (STL) file.")
        if len(var) != self.sim.geom.n_faces:
            raise ValueError(
                f"Variable length ({len(var)}) must match number of facets ({self.sim.geom.n_faces})"
            )

        mesh = self._create_colored_mesh(var, building_ids)
        fig = self._render_scene(mesh, building_ids=building_ids, show=show)
        if fig is not None:
            fig.update_layout(scene=dict(aspectmode="data"))

        self._add_building_outlines_to_scene(building_ids)
        return fig

    def _create_colored_mesh(self, var: np.ndarray, building_ids: Optional[np.ndarray] = None):
        """Create a colored trimesh object from facet data."""
        try:
            import trimesh
        except ImportError as exc:
            raise ImportError("trimesh is required. Install with: pip install trimesh") from exc

        import matplotlib.cm as cm
        import matplotlib.pyplot as plt

        vertices = self.sim.geom.stl.vertices
        faces = self.sim.geom.stl.faces

        face_mask = None
        if building_ids is not None:
            face_to_building = self.sim.geom.get_face_to_building_map()
            building_ids = np.asarray(building_ids)
            face_mask = np.isin(face_to_building, building_ids)

            if np.any(face_mask):
                selected_faces = faces[face_mask]
                selected_var = var[face_mask]
            else:
                print("=" * 67, file=sys.stderr)
                print("WARNING: No valid faces found for the specified building IDs", file=sys.stderr)
                print("=" * 67, file=sys.stderr)
                selected_faces = faces
                selected_var = var
                face_mask = None
        else:
            selected_faces = faces
            selected_var = var

        mesh = trimesh.Trimesh(vertices=vertices, faces=selected_faces, process=False)

        valid_mask = ~np.isnan(selected_var)
        if np.any(valid_mask):
            vmin = np.nanmin(selected_var[valid_mask])
            vmax = np.nanmax(selected_var[valid_mask])
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
            cmap = cm.get_cmap("viridis")

            face_colors = np.ones((len(selected_faces), 4))
            face_colors[valid_mask] = cmap(norm(selected_var[valid_mask]))
            mesh.visual.face_colors = face_colors

        return mesh

    def _render_scene(
        self,
        mesh,
        show_outlines: bool = True,
        angle_threshold: float = 45.0,
        building_ids: Optional[np.ndarray] = None,
        custom_edges: Optional[List[tuple]] = None,
        show: bool = True,
    ):
        """Render the mesh scene using trimesh/plotly."""
        try:
            import trimesh
        except ImportError as exc:
            raise ImportError("trimesh is required. Install with: pip install trimesh") from exc

        meshes = mesh if isinstance(mesh, (list, tuple)) else [mesh]

        scene = trimesh.Scene()
        for current_mesh in meshes:
            scene.add_geometry(current_mesh)

        outline_edges = []
        if show_outlines:
            if custom_edges is not None:
                outline_edges = custom_edges
            else:
                outline_edges = self.sim.geom._calculate_outline_edges(angle_threshold=angle_threshold)

            if building_ids is not None and len(outline_edges) > 0:
                face_to_building = self.sim.geom.get_face_to_building_map()
                building_ids = np.asarray(building_ids)
                faces = self.sim.geom.stl.faces

                filtered_edges = []
                for edge in outline_edges:
                    v0, v1 = edge
                    faces_with_edge = np.where(
                        ((faces[:, 0] == v0) | (faces[:, 1] == v0) | (faces[:, 2] == v0))
                        & ((faces[:, 0] == v1) | (faces[:, 1] == v1) | (faces[:, 2] == v1))
                    )[0]

                    if len(faces_with_edge) > 0 and np.any(
                        np.isin(face_to_building[faces_with_edge], building_ids)
                    ):
                        filtered_edges.append(edge)

                outline_edges = filtered_edges

            if len(outline_edges) > 0:
                vertices = self.sim.geom.stl.vertices
                entities = [trimesh.path.entities.Line([edge[0], edge[1]]) for edge in outline_edges]
                entity_colors = np.tile([0, 0, 0, 255], (len(entities), 1))
                path = trimesh.path.Path3D(entities=entities, vertices=vertices, colors=entity_colors)
                scene.add_geometry(path)

        try:
            from IPython.display import display  # noqa: F401

            in_notebook = True
        except ImportError:
            in_notebook = False

        if in_notebook:
            return self._render_plotly(meshes, outline_edges, show=show)
        if show:
            self._render_trimesh(scene, len(outline_edges))
        return None

    def _render_plotly(self, meshes, outline_edges, show: bool = True):
        """Render using plotly for notebook display."""
        try:
            import plotly.graph_objects as go
            import plotly.io as pio

            pio.renderers.default = "notebook"
            traces = []
            for current_mesh in meshes:
                vertices = current_mesh.vertices
                faces = current_mesh.faces
                colors = current_mesh.visual.face_colors

                opacity = 1.0
                if colors.shape[1] == 4:
                    opacity = np.clip(np.mean(colors[:, 3]) / 255.0, 0.0, 1.0)

                traces.append(
                    go.Mesh3d(
                        x=vertices[:, 0],
                        y=vertices[:, 1],
                        z=vertices[:, 2],
                        i=faces[:, 0],
                        j=faces[:, 1],
                        k=faces[:, 2],
                        facecolor=[f"rgb({c[0]},{c[1]},{c[2]})" for c in colors[:, :3]],
                        opacity=opacity,
                        flatshading=True,
                    )
                )

            fig = go.Figure(data=traces)

            if len(outline_edges) > 0:
                edge_x, edge_y, edge_z = [], [], []
                z_offset = 0.1
                base_vertices = meshes[0].vertices
                for edge in outline_edges:
                    p0 = base_vertices[edge[0]]
                    p1 = base_vertices[edge[1]]
                    z0 = p0[2] + z_offset if abs(p0[2]) < 0.01 else p0[2]
                    z1 = p1[2] + z_offset if abs(p1[2]) < 0.01 else p1[2]
                    edge_x.extend([p0[0], p1[0], None])
                    edge_y.extend([p0[1], p1[1], None])
                    edge_z.extend([z0, z1, None])

                fig.add_trace(
                    go.Scatter3d(
                        x=edge_x,
                        y=edge_y,
                        z=edge_z,
                        mode="lines",
                        line=dict(color="black", width=2),
                        showlegend=False,
                        hoverinfo="skip",
                    )
                )

            fig.update_layout(
                scene=dict(
                    aspectmode="data",
                    xaxis_title="x (m)",
                    yaxis_title="y (m)",
                    zaxis_title="z (m)",
                    xaxis=dict(showgrid=False, showbackground=False),
                    yaxis=dict(showgrid=False, showbackground=False),
                    zaxis=dict(showgrid=False, showbackground=False),
                    camera=dict(
                        projection=dict(type="orthographic"),
                        eye=dict(x=-1.25, y=-1.25, z=1.25),
                    ),
                ),
                showlegend=False,
            )

            if show:
                fig.show()
            return fig
        except ImportError:
            print("Plotly not available. Install with: pip install plotly")
            return None

    def _render_trimesh(self, scene, num_outline_edges, show: bool = True):
        """Render using trimesh viewer."""
        if not show:
            return
        try:
            scene.show()
        except Exception as exc:
            print(f"Could not open trimesh viewer: {exc}")
            print("Install pyglet or pyrender: pip install pyglet")

    def _add_building_outlines_to_scene(self, building_ids=None):
        """Placeholder for MATLAB-compatible call structure."""
        pass

    def _add_building_outlines(self, ax, building_ids=None, angle_threshold: float = 45.0):
        """Add building outline edges to a 3D matplotlib plot."""
        if self.sim.geom is None or not hasattr(self.sim.geom, "stl") or self.sim.geom.stl is None:
            return

        outline_edges = self.sim.geom._calculate_outline_edges(angle_threshold=angle_threshold)
        if len(outline_edges) == 0:
            return

        vertices = self.sim.geom.stl.vertices
        z_offset = 0.1
        for edge in outline_edges:
            p0 = vertices[edge[0]]
            p1 = vertices[edge[1]]
            z0 = p0[2] + z_offset if abs(p0[2]) < 0.01 else p0[2]
            z1 = p1[2] + z_offset if abs(p1[2]) < 0.01 else p1[2]

            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [z0, z1], "k-", linewidth=2, alpha=1.0, zorder=10)

    def plot_fac_type(
        self,
        building_ids: Optional[np.ndarray] = None,
        show_outlines: bool = True,
        angle_threshold: float = 45.0,
        show: bool = True,
    ):
        """Plot the different surface types in the geometry."""
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

        try:
            import trimesh
        except ImportError as exc:
            raise ImportError("trimesh is required for this visualization. Install with: pip install trimesh") from exc

        facids = self.sim.facs["typeid"]
        typeids = self.sim.factypes["id"]
        names = self.sim.factypes["name"]
        unique_ids = np.unique(facids)

        import matplotlib.pyplot as plt

        prop_cycle = plt.rcParams["axes.prop_cycle"]
        default_colors = prop_cycle.by_key()["color"]

        vertices = self.sim.geom.stl.vertices
        faces = self.sim.geom.stl.faces

        face_colors = np.ones((len(faces), 4)) * 255

        type_labels = []
        type_colors = []
        for idx, type_id in enumerate(unique_ids):
            type_mask = facids == type_id
            hex_color = default_colors[idx % len(default_colors)].lstrip("#")
            rgb = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
            face_colors[type_mask, :3] = rgb
            face_colors[type_mask, 3] = 230

            name_idx = np.where(typeids == type_id)[0]
            label = names[name_idx[0]] if len(name_idx) > 0 else f"Type {type_id}"
            type_labels.append(label)
            type_colors.append(rgb)

        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
        mesh.visual.face_colors = face_colors

        try:
            from IPython.display import display  # noqa: F401

            in_notebook = True
        except ImportError:
            in_notebook = False

        if in_notebook:
            print(f"Rendering {len(mesh.faces)} faces for notebook display...")

            try:
                import plotly.graph_objects as go
                import plotly.io as pio

                pio.renderers.default = "notebook"
                fig = go.Figure()

                for idx, type_id in enumerate(unique_ids):
                    type_mask = facids == type_id
                    type_face_indices = np.where(type_mask)[0]
                    if len(type_face_indices) == 0:
                        continue

                    rgb = type_colors[idx]
                    color_str = f"rgb({rgb[0]},{rgb[1]},{rgb[2]})"
                    type_faces = faces[type_face_indices]

                    fig.add_trace(
                        go.Mesh3d(
                            x=vertices[:, 0],
                            y=vertices[:, 1],
                            z=vertices[:, 2],
                            i=type_faces[:, 0],
                            j=type_faces[:, 1],
                            k=type_faces[:, 2],
                            color=color_str,
                            opacity=1.0,
                            flatshading=True,
                            name=type_labels[idx],
                            showlegend=True,
                        )
                    )

                if show_outlines:
                    outline_edges = self.sim.geom._calculate_outline_edges(angle_threshold=angle_threshold)
                    if len(outline_edges) > 0:
                        print(f"Added {len(outline_edges)} outline edges")
                        edge_x, edge_y, edge_z = [], [], []
                        z_offset = 0.1
                        for edge in outline_edges:
                            p0 = vertices[edge[0]]
                            p1 = vertices[edge[1]]
                            z0 = p0[2] + z_offset if abs(p0[2]) < 0.01 else p0[2]
                            z1 = p1[2] + z_offset if abs(p1[2]) < 0.01 else p1[2]
                            edge_x.extend([p0[0], p1[0], None])
                            edge_y.extend([p0[1], p1[1], None])
                            edge_z.extend([z0, z1, None])

                        fig.add_trace(
                            go.Scatter3d(
                                x=edge_x,
                                y=edge_y,
                                z=edge_z,
                                mode="lines",
                                line=dict(color="black", width=2),
                                name="Outlines",
                                showlegend=False,
                                hoverinfo="skip",
                            )
                        )

                fig.update_layout(
                    scene=dict(
                        aspectmode="data",
                        xaxis_title="x (m)",
                        yaxis_title="y (m)",
                        zaxis_title="z (m)",
                        xaxis=dict(showgrid=False, showbackground=False),
                        yaxis=dict(showgrid=False, showbackground=False),
                        zaxis=dict(showgrid=False, showbackground=False),
                        camera=dict(
                            projection=dict(type="orthographic"),
                            eye=dict(x=-1.25, y=-1.25, z=1.25),
                        ),
                    ),
                    title="Surface Types",
                    showlegend=True,
                )

                if show:
                    fig.show()
                return fig

            except ImportError:
                print("Plotly not available. Falling back to static rendering.")
                print("Install plotly for interactive 3D: pip install plotly")
                scene = trimesh.Scene(mesh)
                if show_outlines:
                    outline_edges = self.sim.geom._calculate_outline_edges(angle_threshold=angle_threshold)
                    if len(outline_edges) > 0:
                        entities = [trimesh.path.entities.Line([edge[0], edge[1]]) for edge in outline_edges]
                        entity_colors = np.tile([0, 0, 0, 255], (len(entities), 1))
                        path = trimesh.path.Path3D(entities=entities, vertices=vertices, colors=entity_colors)
                        scene.add_geometry(path)
                if show:
                    try:
                        scene.show()
                    except Exception:
                        print("Could not display. Try installing: pip install plotly or pyglet")
        else:
            print(f"Opening trimesh viewer with {len(mesh.faces)} faces...")
            scene = trimesh.Scene(mesh)
            if show_outlines:
                outline_edges = self.sim.geom._calculate_outline_edges(angle_threshold=angle_threshold)
                if len(outline_edges) > 0:
                    print(f"Added {len(outline_edges)} outline edges")
                    entities = [trimesh.path.entities.Line([edge[0], edge[1]]) for edge in outline_edges]
                    entity_colors = np.tile([0, 0, 0, 255], (len(entities), 1))
                    path = trimesh.path.Path3D(entities=entities, vertices=vertices, colors=entity_colors)
                    scene.add_geometry(path)
            if show:
                try:
                    scene.show()
                except Exception as exc:
                    print(f"Could not open trimesh viewer: {exc}")
                    print("You may need to install pyglet or pyrender: pip install pyglet")

    def plot_building_ids(self, show: bool = True):
        """Plot building IDs from above (x,y view) with distinct colors."""
        if self.sim.geom is None:
            raise ValueError(
                "This method requires a geometry (STL) file. Ensure stl_file is specified in namoptions."
            )

        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError("matplotlib is required for visualization. Install with: pip install matplotlib") from exc

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

        try:
            import matplotlib.pyplot as plt
            from matplotlib.collections import PatchCollection
            from matplotlib.patches import Polygon as mplPolygon
        except ImportError as exc:
            raise ImportError("matplotlib is required for visualization. Install with: pip install matplotlib") from exc

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
