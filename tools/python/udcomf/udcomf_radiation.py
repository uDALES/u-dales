"""Pedestrian receptor geometry, visibility, and irradiance.

The radiation model uses isotropic sky radiance and Lambertian facet exitance.
Angular ray integration resolves partial occlusion without mistaking
facet-to-facet View3D factors for receptor view factors.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import operator
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from udprep.solar import nsun_from_angles
from .heights import configured_receptor_heights, height_tag, validate_heights

if TYPE_CHECKING:
    from udbase import UDBase


_PLANE_NAMES = (
    "upface", "downface", "northface", "southface", "eastface", "westface",
)
_SKY = -1
_GROUND = -2


@dataclass(frozen=True)
class ReceptorGrid:
    """Cell-centre receptors on the native horizontal grid.

    ``valid`` excludes scalar cells marked solid by IBM. Flat indices use
    NumPy C-order on the (x, y) grid, and can be processed in small batches.
    """

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    valid: np.ndarray

    def valid_indices(self) -> np.ndarray:
        return np.flatnonzero(self.valid)

    def points(self, flat_indices: np.ndarray) -> np.ndarray:
        indices = np.asarray(flat_indices)
        if indices.ndim != 1 or not np.issubdtype(indices.dtype, np.integer):
            raise ValueError("flat_indices must be a one-dimensional integer array")
        if np.any(indices < 0) or np.any(indices >= self.valid.size):
            raise IndexError("Receptor index is outside the horizontal grid")
        if not np.all(self.valid.ravel()[indices]):
            raise ValueError("Receptor indices include solid or invalid cells")
        i, j = np.unravel_index(indices, self.valid.shape)
        return np.column_stack((self.x[i], self.y[j], self.z[i, j]))


@dataclass(frozen=True)
class ShortwaveState:
    """Archived atmospheric forcing and reflected facet exitance at one time."""

    time: float
    dni: float
    dsky: float
    zenith: float
    azimuth_local: float
    facet_exitance: np.ndarray


@dataclass(frozen=True)
class LongwaveState:
    """Archived sky forcing and emitted facet exitance at one EB output time."""

    time: float
    sky_irradiance: float
    facet_exitance: np.ndarray


@dataclass(frozen=True)
class ShortwaveRayMap:
    """Static first-hit sources for one receptor and angular quadrature."""

    point: tuple[float, float, float]
    n_mu: int
    n_azimuth: int
    sources: np.ndarray


def facet_shortwave_exitance(netsw: np.ndarray, albedo: np.ndarray) -> np.ndarray:
    """Recover Lambertian reflected exitance from saved absorbed shortwave.

    For an opaque facet, absorbed=(1-albedo)*incident and reflected=
    albedo*incident. This is undefined at albedo=1, so those cases must have
    a separately saved outgoing flux instead of inventing one here.
    """
    netsw = np.asarray(netsw, dtype=float)
    albedo = np.asarray(albedo, dtype=float)
    if netsw.ndim != 1 or netsw.shape != albedo.shape:
        raise ValueError("netsw and albedo must be matching one-dimensional facet arrays")
    if not np.isfinite(netsw).all() or np.any(netsw < 0):
        raise ValueError("Facet net absorbed shortwave must be finite and nonnegative")
    if not np.isfinite(albedo).all() or np.any((albedo < 0) | (albedo >= 1)):
        raise ValueError("Facet albedo must be finite and in [0, 1) to recover reflection")
    return albedo * netsw / (1.0 - albedo)


@lru_cache(maxsize=8)
def _sphere_quadrature(n_mu: int, n_azimuth: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre elevation and midpoint azimuth on both hemispheres."""
    if n_mu < 2 or n_azimuth < 8:
        raise ValueError("Use at least two elevation and eight azimuth nodes")
    nodes, node_weights = np.polynomial.legendre.leggauss(n_mu)
    mu = 0.5 * (nodes + 1.0)
    mu_weights = 0.5 * node_weights
    azimuth = 2.0 * np.pi * (np.arange(n_azimuth) + 0.5) / n_azimuth
    directions = []
    weights = []
    for sign in (1.0, -1.0):
        vertical = np.repeat(sign * mu, n_azimuth)
        phi = np.tile(azimuth, n_mu)
        horizontal = np.sqrt(1.0 - vertical**2)
        directions.append(np.column_stack((horizontal * np.cos(phi), horizontal * np.sin(phi), vertical)))
        weights.append(np.repeat(mu_weights, n_azimuth) * (2.0 * np.pi / n_azimuth))
    return np.vstack(directions), np.concatenate(weights)


class UDComfRadiation:
    """Use a loaded UDBase case for pedestrian geometry and visibility.

    Ray tracing is deliberately lazy: no receptor-by-facet matrix is built.
    Facet queries return centroid visibility for selected facet IDs only;
    this is insufficient by itself for finite-area view factors.
    """

    def __init__(self, sim: UDBase):
        self.sim = sim
        self._ray_mesh = None

    def receptor_grid(self, receptor_height: float | None = None) -> ReceptorGrid:
        """Locate fluid-cell receptors at a height above model ground.

        uDALES currently uses a flat model ground at z=0. No variable-terrain
        elevation is inferred from building facets or an unverified data source.
        """
        sim = self.sim
        height = float(sim.receptor_height if receptor_height is None else receptor_height)
        if not np.isfinite(height) or height <= 0:
            raise ValueError("receptor_height must be finite and positive")
        x = np.asarray(sim.xt, dtype=float)
        y = np.asarray(sim.yt, dtype=float)
        edges = np.asarray(sim.zm, dtype=float)
        top = float(sim.zsize)
        solid = getattr(sim, "Sc", None)
        if solid is None:
            raise ValueError("solid_c.txt is required to exclude building receptors")
        if x.ndim != 1 or y.ndim != 1 or edges.ndim != 1 or edges.size == 0:
            raise ValueError("UDBase grid coordinates must be one-dimensional")
        if not (np.isfinite(x).all() and np.isfinite(y).all() and np.isfinite(edges).all()
                and np.isfinite(top) and np.all(np.diff(edges) > 0) and top > edges[-1]):
            raise ValueError("UDBase grid coordinates are not finite and ordered")
        if solid.shape != (x.size, y.size, edges.size):
            raise ValueError("Scalar solid mask shape does not match the UDBase grid")
        if height >= top:
            raise ValueError("receptor_height must be below the model top")

        k = int(np.searchsorted(edges, height, side="right") - 1)
        if k < 0 or k >= edges.size:
            raise ValueError("receptor_height is outside the scalar-cell grid")
        valid = ~np.asarray(solid[:, :, k], dtype=bool)
        z = np.full((x.size, y.size), height, dtype=float)
        return ReceptorGrid(x=x, y=y, z=z, valid=valid)

    def plane_normals(self) -> dict[str, np.ndarray]:
        """Return receiving-plane outward normals in local model coordinates.

        Geographic azimuth is clockwise from true north. The same
        ``azimuth - xazimuth`` convention as udprep.solar is used here.
        """
        rotation = float(self.sim.xazimuth)
        if not np.isfinite(rotation):
            raise ValueError("xazimuth must be finite")
        normals = {
            "upface": np.array([0.0, 0.0, 1.0]),
            "downface": np.array([0.0, 0.0, -1.0]),
        }
        for name, azimuth in (
            ("northface", 0.0), ("eastface", 90.0),
            ("southface", 180.0), ("westface", 270.0),
        ):
            normal = nsun_from_angles(90.0, azimuth - rotation)
            normal[2] = 0.0
            normals[name] = normal
        return normals

    @staticmethod
    def _point(point: np.ndarray) -> np.ndarray:
        point = np.asarray(point, dtype=float)
        if point.shape != (3,) or not np.isfinite(point).all() or point[2] <= 0:
            raise ValueError("Receptor point must be a finite (x, y, z>0) vector")
        return point

    def _mesh(self):
        if self._ray_mesh is not None:
            return self._ray_mesh
        geom = getattr(self.sim, "geom", None)
        stl = getattr(geom, "stl", None)
        if stl is None:
            if getattr(self.sim, "stl_file", None):
                raise ValueError("Load the case geometry in UDBase before visibility queries")
            return None

        import pyvista as pv

        vertices = np.asarray(stl.vertices, dtype=float)
        faces = np.asarray(stl.faces, dtype=np.int64)
        if faces.ndim != 2 or faces.shape[1] != 3 or vertices.ndim != 2 or vertices.shape[1] != 3:
            raise ValueError("Receptor visibility requires a triangular STL mesh")
        if not np.isfinite(vertices).all():
            raise ValueError("STL contains non-finite vertices")
        facet_types = getattr(self.sim, "facs", {}).get("typeid")
        if facet_types is not None and len(facet_types) != len(faces):
            raise ValueError("STL face count does not match facets.inp ordering")
        facet_normals = getattr(self.sim, "facs", {}).get("normals")
        if facet_normals is not None:
            facet_normals = np.asarray(facet_normals, dtype=float)
            mesh_normals = np.asarray(stl.face_normals, dtype=float)
            if facet_normals.shape != mesh_normals.shape:
                raise ValueError("STL normals do not match facets.inp ordering or orientation")
            lengths = np.linalg.norm(facet_normals, axis=1)
            if (not np.isfinite(lengths).all() or np.any(lengths == 0)
                    or not np.isfinite(mesh_normals).all()
                    or np.any(np.einsum("ij,ij->i", facet_normals / lengths[:, None], mesh_normals) < 0.99)):
                raise ValueError("STL normals do not match facets.inp ordering or orientation")
        if faces.size == 0:
            return None
        vtk_faces = np.empty((len(faces), 4), dtype=np.int64)
        vtk_faces[:, 0] = 3
        vtk_faces[:, 1:] = faces
        self._ray_mesh = pv.PolyData(vertices, vtk_faces.ravel())
        return self._ray_mesh

    def _clear_segment(
        self, start: np.ndarray, end: np.ndarray, *, sky: bool,
        target_facet: int | None = None,
    ) -> bool:
        direction = end - start
        length = float(np.linalg.norm(direction))
        if not np.isfinite(length) or length <= 0:
            raise ValueError("Ray endpoint must differ from its receptor")
        unit = direction / length
        tolerance = min(1.0e-5, 1.0e-7 * length)
        if unit[2] < 0 and -start[2] / unit[2] < length - tolerance:
            return False  # Model ground is opaque even if omitted from the STL.
        mesh = self._mesh()
        if mesh is None:
            return True
        origin = start + min(1.0e-7, length * 1.0e-5) * unit
        trace_end = end
        if target_facet is not None:
            # VTK can omit an intersection exactly at the segment endpoint.
            extension = max(
                1.0e-6,
                16.0 * np.finfo(np.float32).eps * max(1.0, np.max(np.abs(end))),
            )
            trace_end = end + extension * unit
        hits, hit_facets = mesh.ray_trace(origin, trace_end, first_point=True)
        if len(hits) == 0:
            return target_facet is None
        if sky:
            return False
        if target_facet is not None:
            return int(hit_facets[0]) == target_facet
        hit_distance = float(np.linalg.norm(np.asarray(hits).reshape(3) - start))
        return hit_distance >= length - tolerance

    def sky_visibility(self, point: np.ndarray, directions: np.ndarray) -> np.ndarray:
        """Test unobstructed directions to the upper sky for one receptor."""
        point = self._point(point)
        directions = np.asarray(directions, dtype=float)
        if directions.ndim != 2 or directions.shape[1] != 3 or not np.isfinite(directions).all():
            raise ValueError("directions must be finite (n, 3) vectors")
        lengths = np.linalg.norm(directions, axis=1)
        if np.any(lengths == 0):
            raise ValueError("Sky directions must be nonzero")
        result = np.zeros(len(directions), dtype=bool)
        if not np.any(directions[:, 2] > 0):
            return result
        mesh = self._mesh()
        if mesh is None:
            return directions[:, 2] > 0
        bounds = np.asarray(mesh.bounds, dtype=float)
        far_corner = np.maximum(np.abs(point - bounds[0]), np.abs(point - bounds[1]))
        distance = float(np.linalg.norm(far_corner) + 1.0)
        for i, direction in enumerate(directions):
            if direction[2] > 0:
                end = point + distance * direction / lengths[i]
                result[i] = self._clear_segment(point, end, sky=True)
        return result

    def facet_visibility(self, point: np.ndarray, facet_indices: np.ndarray) -> np.ndarray:
        """Test centroid visibility of selected front-facing facets from one receptor.

        This is a diagnostic convenience. Use ``facet_sample_visibility``
        with an integration rule for finite-area radiative exchange.
        """
        samples = np.array([[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]])
        return self.facet_sample_visibility(point, facet_indices, samples)[:, 0]

    def facet_sample_visibility(
        self, point: np.ndarray, facet_indices: np.ndarray,
        barycentric_samples: np.ndarray,
    ) -> np.ndarray:
        """Test visibility at specified triangle samples, shape (facet, sample).

        Each sample is a nonnegative barycentric triplet summing to one. The
        caller chooses the quadrature rule; this routine applies no area or
        projected-solid-angle weights.
        """
        point = self._point(point)
        indices = np.asarray(facet_indices)
        if indices.ndim != 1 or not np.issubdtype(indices.dtype, np.integer):
            raise ValueError("facet_indices must be a one-dimensional integer array")
        samples = np.asarray(barycentric_samples, dtype=float)
        if (samples.ndim != 2 or samples.shape[1] != 3 or samples.shape[0] == 0
                or not np.isfinite(samples).all() or np.any(samples < 0)
                or not np.allclose(samples.sum(axis=1), 1.0, rtol=0, atol=1e-12)):
            raise ValueError("barycentric_samples must be nonnegative (n, 3) rows summing to one")
        mesh = self._mesh()
        if mesh is None:
            if indices.size:
                raise ValueError("No facets are loaded for visibility queries")
            return np.zeros((0, len(samples)), dtype=bool)
        if np.any(indices < 0) or np.any(indices >= mesh.n_cells):
            raise IndexError("Facet index is outside the STL mesh")
        stl = self.sim.geom.stl
        normals = np.asarray(stl.face_normals, dtype=float)
        result = np.zeros((len(indices), len(samples)), dtype=bool)
        for i, facet in enumerate(indices):
            triangle = np.asarray(stl.triangles[facet], dtype=float)
            for j, weights in enumerate(samples):
                target = weights @ triangle
                if np.dot(normals[facet], point - target) > 0:
                    result[i, j] = self._clear_segment(
                        point, target, sky=False, target_facet=int(facet)
                    )
        return result

    def direct_solar_visibility(
        self, point: np.ndarray, zenith: float, azimuth_local: float
    ) -> bool:
        """Test the sun ray independently of diffuse-sky visibility."""
        point = self._point(point)
        if not np.isfinite(zenith) or not np.isfinite(azimuth_local):
            raise ValueError("Solar angles must be finite")
        if not 0 <= zenith < 90:
            return False
        direction = nsun_from_angles(float(zenith), float(azimuth_local))
        return bool(self.sky_visibility(point, direction[None, :])[0])

    def load_shortwave_state(self, time_index: int) -> ShortwaveState:
        """Read one preprocessing timestamp through UDBase's existing loaders."""
        index = operator.index(time_index)
        times = np.asarray(self.sim.load_shortwave_forcing("time"), dtype=float)
        if times.ndim != 1 or not 0 <= index < len(times):
            raise IndexError("Shortwave forcing time_index is outside the archive")
        facet_data = self.sim.load_timedepsw(time_index=index)
        time = float(times[index])
        if not np.isclose(float(facet_data["time"]), time, rtol=0.0, atol=0.005):
            raise ValueError("Shortwave forcing and timedepsw timestamps do not match")

        forcing = {}
        for name in ("dni", "dsky", "solar_zenith", "solar_azimuth_local"):
            values = np.asarray(self.sim.load_shortwave_forcing(name), dtype=float)
            if values.shape != times.shape:
                raise ValueError(f"Shortwave forcing {name} has an invalid time axis")
            forcing[name] = float(values[index])
        exitance = facet_shortwave_exitance(
            np.asarray(facet_data["netsw"], dtype=float),
            self.sim.assign_prop_to_fac("al"),
        )
        return ShortwaveState(
            time=time,
            dni=forcing["dni"],
            dsky=forcing["dsky"],
            zenith=forcing["solar_zenith"],
            azimuth_local=forcing["solar_azimuth_local"],
            facet_exitance=exitance,
        )

    def load_longwave_state(self, time_index: int) -> LongwaveState:
        """Read one EB record and interpolate saved sky LW to its timestamp.

        uDALES writes LWout=emissivity*sigma*T_surface**4 after updating facet
        temperature. LWin is emissivity-weighted absorption, not a source
        radiance for pedestrian receptors.
        """
        index = operator.index(time_index)
        times = np.asarray(self.sim.load_fac_eb("t"), dtype=float)
        if times.ndim != 1 or not 0 <= index < len(times):
            raise IndexError("Facet EB time_index is outside the archive")
        if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
            raise ValueError("Facet EB timestamps must be finite and increasing")
        time = float(times[index])
        forcing = self.sim.load_timedeplw()
        forcing_times = np.asarray(forcing["time"], dtype=float)
        sky_flux = np.asarray(forcing["LWsky"], dtype=float)
        if (forcing_times.ndim != 1 or sky_flux.shape != forcing_times.shape
                or not np.isfinite(forcing_times).all()
                or not np.isfinite(sky_flux).all()
                or np.any(np.diff(forcing_times) <= 0) or np.any(sky_flux < 0)):
            raise ValueError("Sky longwave forcing must be finite, nonnegative, and ordered")
        if not forcing_times.size or time < forcing_times[0] or time > forcing_times[-1]:
            raise ValueError("Facet EB timestamp is outside saved sky longwave forcing")
        exitance = np.asarray(self.sim.load_fac_eb("LWout", time_index=index), dtype=float)
        if exitance.ndim != 1 or not np.isfinite(exitance).all() or np.any(exitance < 0):
            raise ValueError("Facet emitted longwave exitance must be finite and nonnegative")
        return LongwaveState(
            time=time,
            sky_irradiance=float(np.interp(time, forcing_times, sky_flux)),
            facet_exitance=exitance,
        )

    def _directional_source(self, point: np.ndarray, direction: np.ndarray) -> int:
        """Return the first facet ID, sky, or opaque unmeshed model ground."""
        mesh = self._mesh()
        if mesh is None:
            return _SKY if direction[2] > 0 else _GROUND
        bounds = np.asarray(mesh.bounds, dtype=float)
        far_corner = np.maximum(np.abs(point - bounds[0]), np.abs(point - bounds[1]))
        distance = float(np.linalg.norm(far_corner) + 1.0)
        origin = point + 1.0e-7 * direction
        end = point + distance * direction
        hits, hit_facets = mesh.ray_trace(origin, end, first_point=True)
        if len(hits):
            if direction[2] < 0:
                ground_distance = -point[2] / direction[2]
                hit_distance = float(np.linalg.norm(np.asarray(hits).reshape(3) - point))
                precision = 16.0 * np.finfo(np.float32).eps * max(
                    1.0, np.max(np.abs(point)), ground_distance
                )
                if hit_distance > ground_distance + precision:
                    return _GROUND
            facet = int(hit_facets[0])
            normal = np.asarray(self.sim.geom.stl.face_normals[facet], dtype=float)
            return facet if np.dot(normal, direction) < 0 else _GROUND
        return _SKY if direction[2] > 0 else _GROUND

    def _plane_geometry(
        self, n_mu: int, n_azimuth: int
    ) -> tuple[np.ndarray, np.ndarray]:
        directions, weights = _sphere_quadrature(
            operator.index(n_mu), operator.index(n_azimuth)
        )
        normals = self.plane_normals()
        normal_matrix = np.stack([normals[name] for name in _PLANE_NAMES])
        projected_weights = (
            np.maximum(normal_matrix @ directions.T, 0.0) * weights[None, :] / np.pi
        )
        # Every receiving plane integrates a constant radiance to pi*L.
        # Enforce that exact zeroth moment after finite angular quadrature.
        projected_weights /= projected_weights.sum(axis=1, keepdims=True)
        return directions, projected_weights

    def trace_shortwave_rays(
        self, point: np.ndarray, *, n_mu: int = 8, n_azimuth: int = 32
    ) -> ShortwaveRayMap:
        """Trace the time-invariant sky/facet directions once for reuse."""
        point = self._point(point)
        n_mu = operator.index(n_mu)
        n_azimuth = operator.index(n_azimuth)
        directions, _ = _sphere_quadrature(n_mu, n_azimuth)
        sources = np.array(
            [self._directional_source(point, direction) for direction in directions],
            dtype=np.int32,
        )
        sources.setflags(write=False)
        return ShortwaveRayMap(tuple(float(value) for value in point), n_mu, n_azimuth, sources)

    def _shortwave_at_receptor(
        self, point: np.ndarray, state: ShortwaveState,
        directions: np.ndarray, projected_weights: np.ndarray,
        sources: np.ndarray | None = None,
    ) -> dict[str, float]:
        point = self._point(point)
        if getattr(self.sim, "ltrees", False):
            raise NotImplementedError(
                "Pedestrian shortwave does not yet account for vegetation attenuation"
            )
        if not all(np.isfinite(value) for value in (
            state.dni, state.dsky, state.zenith, state.azimuth_local
        )) or state.dni < 0 or state.dsky < 0:
            raise ValueError("Atmospheric shortwave forcing must be finite and nonnegative")
        if not 0 <= state.zenith <= 180:
            raise ValueError("Solar zenith must be between 0 and 180 degrees")
        exitance = np.asarray(state.facet_exitance, dtype=float)
        if exitance.ndim != 1 or not np.isfinite(exitance).all() or np.any(exitance < 0):
            raise ValueError("Facet shortwave exitance must be finite and nonnegative")
        mesh = self._mesh()
        n_facets = mesh.n_cells if mesh is not None else 0
        if len(exitance) != n_facets:
            raise ValueError("Facet shortwave exitance count does not match the STL mesh")

        if sources is None:
            source_ids = np.array(
                [self._directional_source(point, direction) for direction in directions],
                dtype=np.int32,
            )
        else:
            source_ids = np.asarray(sources)
            if (source_ids.shape != (len(directions),)
                    or not np.issubdtype(source_ids.dtype, np.integer)):
                raise ValueError("Shortwave ray map has an invalid source array")
        if np.any(source_ids < _GROUND) or np.any(source_ids >= n_facets):
            raise ValueError("Shortwave ray map references an invalid facet")
        ray_exitance = np.zeros(len(source_ids), dtype=float)
        ray_exitance[source_ids == _SKY] = state.dsky
        facet_rays = source_ids >= 0
        ray_exitance[facet_rays] = exitance[source_ids[facet_rays]]
        irradiances = projected_weights @ ray_exitance
        result = {
            f"sw_nondirect_{name}": float(value)
            for name, value in zip(_PLANE_NAMES, irradiances)
        }
        result["sw_direct_normal"] = (
            float(state.dni) if state.dni > 0 and self.direct_solar_visibility(
                point, state.zenith, state.azimuth_local
            ) else 0.0
        )
        return result

    def shortwave_at_receptor(
        self, point: np.ndarray, state: ShortwaveState,
        *, n_mu: int = 8, n_azimuth: int = 32,
        ray_map: ShortwaveRayMap | None = None,
    ) -> dict[str, float]:
        """Return seven irradiances [W m-2] at one fluid receptor.

        ``dsky/pi`` is isotropic upper-sky radiance. Each facet is Lambertian
        with radiance ``facet_exitance/pi``. The six non-direct planes contain
        visible sky diffuse plus all visible facet reflection; the direct
        beam appears only in ``sw_direct_normal``.
        """
        point = self._point(point)
        if ray_map is not None and (
            ray_map.point != tuple(float(value) for value in point)
            or ray_map.n_mu != n_mu or ray_map.n_azimuth != n_azimuth
        ):
            raise ValueError("Shortwave ray map does not match receptor or quadrature")
        directions, projected_weights = self._plane_geometry(n_mu, n_azimuth)
        return self._shortwave_at_receptor(
            point, state, directions, projected_weights,
            sources=ray_map.sources if ray_map is not None else None,
        )

    def shortwave_plane(
        self, state: ShortwaveState, *, flat_indices: np.ndarray | None = None,
        n_mu: int = 8, n_azimuth: int = 32,
        receptor_height: float | None = None,
    ) -> dict[str, np.ndarray]:
        """Calculate native (x, y) fields, leaving invalid/unselected cells NaN.

        ``flat_indices`` permits small batches and future checkpointing. A full
        Paris plane would require a very large number of ray queries and has
        not been performance-qualified by this implementation.
        """
        grid = self.receptor_grid(receptor_height)
        indices = grid.valid_indices() if flat_indices is None else np.asarray(flat_indices)
        points = grid.points(indices)
        directions, projected_weights = self._plane_geometry(n_mu, n_azimuth)
        names = ("sw_direct_normal",) + tuple(f"sw_nondirect_{name}" for name in _PLANE_NAMES)
        fields = {name: np.full(grid.valid.shape, np.nan, dtype=float) for name in names}
        for flat_index, point in zip(indices, points):
            values = self._shortwave_at_receptor(point, state, directions, projected_weights)
            i, j = np.unravel_index(int(flat_index), grid.valid.shape)
            for name, value in values.items():
                fields[name][i, j] = value
        return fields

    def longwave_at_receptor(
        self, point: np.ndarray, state: LongwaveState,
        *, n_mu: int = 8, n_azimuth: int = 32,
        ray_map: ShortwaveRayMap | None = None,
    ) -> dict[str, float]:
        """Return six longwave irradiances [W m-2] at one fluid receptor.

        Sky radiance is LWsky/pi. Facet radiance is archived LWout/pi, which
        already includes surface temperature and emissivity. uDALES does not
        include reflected longwave in its facet exchange, so none is added.
        """
        point = self._point(point)
        if getattr(self.sim, "ltrees", False):
            raise NotImplementedError(
                "Pedestrian longwave does not yet account for vegetation emission or attenuation"
            )
        if ray_map is not None and (
            ray_map.point != tuple(float(value) for value in point)
            or ray_map.n_mu != n_mu or ray_map.n_azimuth != n_azimuth
        ):
            raise ValueError("Ray map does not match receptor or quadrature")
        sky_irradiance = float(state.sky_irradiance)
        exitance = np.asarray(state.facet_exitance, dtype=float)
        if not np.isfinite(sky_irradiance) or sky_irradiance < 0:
            raise ValueError("Sky longwave irradiance must be finite and nonnegative")
        if exitance.ndim != 1 or not np.isfinite(exitance).all() or np.any(exitance < 0):
            raise ValueError("Facet emitted longwave exitance must be finite and nonnegative")
        mesh = self._mesh()
        n_facets = mesh.n_cells if mesh is not None else 0
        if exitance.size != n_facets:
            raise ValueError("Facet longwave exitance count does not match the STL mesh")

        directions, projected_weights = self._plane_geometry(n_mu, n_azimuth)
        if ray_map is None:
            source_ids = np.array(
                [self._directional_source(point, direction) for direction in directions],
                dtype=np.int32,
            )
        else:
            source_ids = np.asarray(ray_map.sources)
            if (source_ids.shape != (len(directions),)
                    or not np.issubdtype(source_ids.dtype, np.integer)):
                raise ValueError("Ray map has an invalid source array")
        if np.any(source_ids < _GROUND) or np.any(source_ids >= n_facets):
            raise ValueError("Ray map references an invalid facet")
        if np.any(source_ids == _GROUND):
            raise ValueError(
                "Longwave ray reaches unmeshed ground or a back-facing facet; "
                "total irradiance cannot be determined"
            )
        ray_exitance = np.zeros(len(source_ids), dtype=float)
        ray_exitance[source_ids == _SKY] = sky_irradiance
        facet_rays = source_ids >= 0
        ray_exitance[facet_rays] = exitance[source_ids[facet_rays]]
        irradiances = projected_weights @ ray_exitance
        return {f"lw_{name}": float(value) for name, value in zip(_PLANE_NAMES, irradiances)}

    def longwave_plane(
        self, state: LongwaveState, *, flat_indices: np.ndarray | None = None,
        n_mu: int = 8, n_azimuth: int = 32,
        receptor_height: float | None = None,
    ) -> dict[str, np.ndarray]:
        """Calculate native (x, y) longwave fields for selected fluid receptors."""
        grid = self.receptor_grid(receptor_height)
        indices = grid.valid_indices() if flat_indices is None else np.asarray(flat_indices)
        points = grid.points(indices)
        names = tuple(f"lw_{name}" for name in _PLANE_NAMES)
        fields = {name: np.full(grid.valid.shape, np.nan, dtype=float) for name in names}
        for flat_index, point in zip(indices, points):
            values = self.longwave_at_receptor(
                point, state, n_mu=n_mu, n_azimuth=n_azimuth
            )
            i, j = np.unravel_index(int(flat_index), grid.valid.shape)
            for name, value in values.items():
                fields[name][i, j] = value
        return fields

    def _radiation_archive(self, kind: str, times: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
        """Load archived facet/sky series once through UDBase for batch processing."""
        if kind == "shortwave":
            archive = self.sim.load_timedepsw()
            if not np.allclose(archive["time"], times, rtol=0, atol=0.005):
                raise ValueError("timedepsw and shortwave forcing timestamps differ")
            exitance = np.asarray(archive["netsw"], dtype=float)
            albedo = np.asarray(self.sim.assign_prop_to_fac("al"), dtype=float)
            if (exitance.shape != (len(albedo), len(times))
                    or not np.isfinite(exitance).all() or np.any(exitance < 0)
                    or not np.isfinite(albedo).all()
                    or np.any((albedo < 0) | (albedo >= 1))):
                raise ValueError("Invalid facet shortwave flux or albedo")
            mesh = self._mesh()
            if exitance.shape[0] != (mesh.n_cells if mesh is not None else 0):
                raise ValueError("Facet shortwave archive does not match the STL mesh")
            exitance *= (albedo / (1.0 - albedo))[:, None]
            forcing = {
                name: np.asarray(self.sim.load_shortwave_forcing(name), dtype=float)
                for name in ("dni", "dsky", "solar_zenith", "solar_azimuth_local")
            }
            if (any(values.shape != times.shape or not np.isfinite(values).all()
                    for values in forcing.values())
                    or np.any(forcing["dni"] < 0) or np.any(forcing["dsky"] < 0)
                    or np.any((forcing["solar_zenith"] < 0) | (forcing["solar_zenith"] > 180))):
                raise ValueError("Invalid archived shortwave forcing")
            return exitance, forcing["dsky"], forcing

        exitance = np.asarray(self.sim.load_fac_eb("LWout"), dtype=float)
        if (exitance.ndim != 2 or exitance.shape[1] != len(times)
                or not np.isfinite(exitance).all() or np.any(exitance < 0)):
            raise ValueError("Invalid archived facet longwave exitance")
        mesh = self._mesh()
        if exitance.shape[0] != (mesh.n_cells if mesh is not None else 0):
            raise ValueError("Facet longwave archive does not match the STL mesh")
        forcing = self.sim.load_timedeplw()
        lw_times = np.asarray(forcing["time"], dtype=float)
        lw_sky = np.asarray(forcing["LWsky"], dtype=float)
        if (lw_times.ndim != 1 or lw_times.size < 2 or lw_sky.shape != lw_times.shape
                or not np.isfinite(lw_times).all() or np.any(np.diff(lw_times) <= 0)
                or not np.isfinite(lw_sky).all() or np.any(lw_sky < 0)
                or times[0] < lw_times[0] or times[-1] > lw_times[-1]):
            raise ValueError("Facet times are outside valid sky-longwave forcing")
        return exitance, np.interp(times, lw_times, lw_sky), {}

    def _source_series(
        self, point: np.ndarray, kind: str, used: np.ndarray,
        exitance: np.ndarray, sky: np.ndarray, forcing: dict,
        projected_weights: np.ndarray, n_mu: int, n_azimuth: int,
    ) -> np.ndarray:
        ray_map = self.trace_shortwave_rays(point, n_mu=n_mu, n_azimuth=n_azimuth)
        sources = ray_map.sources
        if np.any(sources == _GROUND) and kind == "longwave":
            raise ValueError("Longwave ray reaches unmeshed ground or a back-facing facet")
        if np.any(sources >= exitance.shape[0]):
            raise ValueError("Ray map references a facet outside the archived source")
        rays = np.zeros((len(sources), len(used)), dtype=float)
        sky_rays = sources == _SKY
        facet_rays = sources >= 0
        rays[sky_rays] = sky[used]
        rays[facet_rays] = exitance[sources[facet_rays]][:, used]
        nondirect = projected_weights @ rays
        if kind == "longwave":
            return nondirect

        direct = np.zeros(len(used), dtype=float)
        daylight = np.flatnonzero(
            (forcing["dni"][used] > 0) & (forcing["solar_zenith"][used] < 90)
        )
        if daylight.size:
            sun_dirs = np.stack([
                nsun_from_angles(
                    float(forcing["solar_zenith"][used[i]]),
                    float(forcing["solar_azimuth_local"][used[i]]),
                )
                for i in daylight
            ])
            direct[daylight] = (
                forcing["dni"][used[daylight]]
                * self.sky_visibility(point, sun_dirs)
            )
        return np.vstack((direct, nondirect))

    def write_hourly(
        self, kind: str, target_times: np.ndarray, *,
        receptor_height: float | None = None,
        flat_indices: np.ndarray | None = None,
        output_path: Path | None = None,
        checkpoint_dir: Path | None = None,
        tile_size: int = 256,
        n_mu: int = 8,
        n_azimuth: int = 32,
        max_gap: float | None = None,
        resume: bool = True,
        overwrite: bool = False,
    ) -> Path:
        """Write preceding-15-minute hourly irradiance planes.

        Source samples are integrated under piecewise-linear interpolation,
        not claimed to be exact sub-cadence physical fluxes. Missing coverage
        and restart gaps produce NaN windows. Checkpoints are atomic per hour
        and receptor tile; resume requires identical input file signatures.
        """
        from .checkpoints import (
            RadiationCheckpoints, consolidate_radiation, file_signature,
            hourly_windows, selection_hash,
        )

        if kind not in ("shortwave", "longwave"):
            raise ValueError("kind must be 'shortwave' or 'longwave'")
        if getattr(self.sim, "ltrees", False):
            raise NotImplementedError("Pedestrian radiation does not include vegetation")
        tile_size = operator.index(tile_size)
        n_mu = operator.index(n_mu)
        n_azimuth = operator.index(n_azimuth)
        if tile_size <= 0:
            raise ValueError("tile_size must be positive")
        _, projected_weights = self._plane_geometry(n_mu, n_azimuth)
        grid = self.receptor_grid(receptor_height)
        indices = grid.valid_indices() if flat_indices is None else np.asarray(flat_indices)
        grid.points(indices)
        if len(indices) == 0 or len(np.unique(indices)) != len(indices):
            raise ValueError("Select at least one unique fluid receptor")
        indices = np.sort(indices)
        times = np.asarray(
            self.sim.load_shortwave_forcing("time") if kind == "shortwave"
            else self.sim.load_fac_eb("t"), dtype=float,
        )
        windows, gap = hourly_windows(times, target_times, max_gap=max_gap)
        names = (
            ("sw_direct_normal",) + tuple(f"sw_nondirect_{name}" for name in _PLANE_NAMES)
            if kind == "shortwave" else tuple(f"lw_{name}" for name in _PLANE_NAMES)
        )
        case = Path(self.sim.path)
        expnr = self.sim.expnr
        output = (
            Path(output_path) if output_path is not None
            else case / (
                f"pedestrian_{kind}.{expnr}.nc" if receptor_height is None
                else f"pedestrian_{kind}.{height_tag(receptor_height)}.{expnr}.nc"
            )
        )
        if output.exists() and not overwrite:
            raise FileExistsError(f"Radiation output exists: {output}")
        checkpoint_dir = (
            Path(checkpoint_dir) if checkpoint_dir is not None
            else (
                case / "udcomf_radiation.checkpoints" / kind
                if receptor_height is None
                else case / "udcomf_radiation.checkpoints" / kind / height_tag(receptor_height)
            )
        )
        source_files = [
            case / f"namoptions.{expnr}", case / f"facets.inp.{expnr}",
            case / f"factypes.inp.{expnr}", case / "solid_c.txt",
        ]
        profile = case / f"prof.inp.{expnr}"
        if profile.is_file():
            source_files.append(profile)
        stl_file = getattr(self.sim, "stl_file", None)
        if stl_file:
            source_files.append(case / stl_file)
        source_files.extend(
            [case / f"shortwave_forcing.{expnr}.nc", case / f"timedepsw.inp.{expnr}"]
            if kind == "shortwave"
            else [case / f"facEB.{expnr}.nc", case / f"timedeplw.inp.{expnr}"]
        )
        manifest = {
            "version": 2,
            "kind": kind,
            "case": str(case.resolve()),
            "receptor_height_m": float(grid.z.flat[0]),
            "grid_shape": list(grid.valid.shape),
            "selection_sha256": selection_hash(indices),
            "selected_count": len(indices),
            "tile_size": tile_size,
            "n_mu": n_mu,
            "n_azimuth": n_azimuth,
            "target_times_s": [window.end for window in windows],
            "window_duration_s": 900.0,
            "max_source_gap_s": gap,
            "source_files": [file_signature(path) for path in source_files],
        }
        checkpoints = RadiationCheckpoints(checkpoint_dir, manifest, resume=resume)
        missing = []
        for tile_index, offset in enumerate(range(0, len(indices), tile_size)):
            length = min(tile_size, len(indices) - offset)
            missing.append([
                time_index for time_index in range(len(windows))
                if checkpoints.read(time_index, tile_index, (length, len(names))) is None
            ])
        if any(missing):
            complete_windows = [
                windows[time_index] for tile in missing for time_index in tile
                if windows[time_index].complete
            ]
            used = (
                np.unique(np.concatenate([window.indices for window in complete_windows]))
                if complete_windows else np.empty(0, dtype=int)
            )
            if used.size:
                exitance, sky, forcing = self._radiation_archive(kind, times)
            for tile_index, offset in enumerate(range(0, len(indices), tile_size)):
                if not missing[tile_index]:
                    continue
                points = grid.points(indices[offset:offset + tile_size])
                values = {
                    time_index: np.full((len(points), len(names)), np.nan, dtype=np.float32)
                    for time_index in missing[tile_index]
                }
                if used.size:
                    for row, point in enumerate(points):
                        series = self._source_series(
                            point, kind, used, exitance, sky, forcing,
                            projected_weights, n_mu, n_azimuth,
                        )
                        for time_index, tile_values in values.items():
                            window = windows[time_index]
                            if window.complete:
                                positions = np.searchsorted(used, window.indices)
                                tile_values[row] = series[:, positions] @ window.weights
                for time_index, tile_values in values.items():
                    checkpoints.write(time_index, tile_index, tile_values)
        return consolidate_radiation(
            output, checkpoints, grid, indices, names, windows, tile_size,
            overwrite=overwrite,
        )

    def write_hourly_heights(
        self, target_times: np.ndarray, *,
        heights: tuple[float, ...] | None = None,
        kinds: tuple[str, ...] = ("shortwave", "longwave"),
        flat_indices: np.ndarray | None = None,
        output_dir: Path | None = None,
        checkpoint_root: Path | None = None,
        tile_size: int = 256,
        n_mu: int = 8,
        n_azimuth: int = 32,
        max_gap: float | None = None,
        resume: bool = True,
        overwrite: bool = False,
    ) -> dict[float, dict[str, Path]]:
        """Write one native-grid radiation pair per height, processing sequentially."""
        selected = (
            configured_receptor_heights(self.sim.path) if heights is None
            else validate_heights(heights)
        )
        if not kinds or len(set(kinds)) != len(kinds) or any(
            kind not in ("shortwave", "longwave") for kind in kinds
        ):
            raise ValueError("kinds must contain unique shortwave/longwave names")
        case = Path(self.sim.path)
        directory = case if output_dir is None else Path(output_dir)
        checkpoints = (
            case / "udcomf_radiation.checkpoints" if checkpoint_root is None
            else Path(checkpoint_root)
        )
        if not overwrite:
            for height in selected:
                for kind in kinds:
                    path = directory / f"pedestrian_{kind}.{height_tag(height)}.{self.sim.expnr}.nc"
                    if path.exists():
                        raise FileExistsError(f"Radiation output exists: {path}")
        if output_dir is not None:
            directory.mkdir(parents=True, exist_ok=True)
        results: dict[float, dict[str, Path]] = {}
        for height in selected:
            results[height] = {}
            for kind in kinds:
                tag = height_tag(height)
                results[height][kind] = self.write_hourly(
                    kind, target_times, receptor_height=height,
                    flat_indices=flat_indices,
                    output_path=directory / f"pedestrian_{kind}.{tag}.{self.sim.expnr}.nc",
                    checkpoint_dir=checkpoints / kind / tag,
                    tile_size=tile_size, n_mu=n_mu, n_azimuth=n_azimuth,
                    max_gap=max_gap, resume=resume, overwrite=overwrite,
                )
        return results
