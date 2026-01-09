from __future__ import annotations

"""
Numba-accelerated direct shortwave radiation on facets with vegetation attenuation.
"""

from dataclasses import dataclass
from typing import Dict, Tuple

import math
import numpy as np

try:
    import numba as nb
except ImportError:  # pragma: no cover - runtime dependent
    nb = None


@dataclass
class VegData:
    points: np.ndarray  # (n, 3) 1-based (i, j, k)
    lad: np.ndarray     # (n,)
    dec: np.ndarray     # (n,)


def _read_sparse_points(sim_dir: str, expnr: str) -> np.ndarray:
    """Load sparse vegetation point indices (1-based i,j,k) from veg.inp.<expnr>."""
    path = f"{sim_dir}/veg.inp.{expnr}"
    points = []
    with open(path, "r", encoding="ascii") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                continue
            points.append([int(parts[0]), int(parts[1]), int(parts[2])])
    if not points:
        raise ValueError(f"No vegetation points found in {path}")
    return np.asarray(points, dtype=int)


def _read_sparse_params(sim_dir: str, expnr: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load sparse vegetation parameters (LAD and DEC) from veg_params.inp.<expnr>."""
    path = f"{sim_dir}/veg_params.inp.{expnr}"
    lad = []
    dec = []
    with open(path, "r", encoding="ascii") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 7:
                continue
            lad.append(float(parts[1]))
            dec.append(float(parts[4]))
    if not lad:
        raise ValueError(f"No vegetation params found in {path}")
    return np.asarray(lad, dtype=float), np.asarray(dec, dtype=float)


def _load_veg_data(sim_dir: str, expnr: str) -> VegData:
    """Load vegetation points and parameters, validating length consistency."""
    points = _read_sparse_points(sim_dir, expnr)
    lad, dec = _read_sparse_params(sim_dir, expnr)
    if len(points) != len(lad):
        raise ValueError(
            f"veg.inp.{expnr} has {len(points)} points but veg_params has {len(lad)} rows"
        )
    return VegData(points=points, lad=lad, dec=dec)


def _build_veg_fields(sim, veg: VegData, ktot: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Expand sparse vegetation into dense 3D grids up to ktot.

    Returns
    -------
    lad_3d, dec_3d : float arrays (itot, jtot, ktot)
        Vegetation parameters on the grid.
    veg_index : int array (itot, jtot, ktot)
        Maps cell to sparse veg index; -1 means no vegetation.
    """
    lad_3d = np.zeros((sim.itot, sim.jtot, ktot), dtype=float)
    dec_3d = np.zeros((sim.itot, sim.jtot, ktot), dtype=float)
    veg_index = -np.ones((sim.itot, sim.jtot, ktot), dtype=int)

    for idx, (i, j, k) in enumerate(veg.points):
        ii = i - 1
        jj = j - 1
        kk = k - 1
        lad_3d[ii, jj, kk] = veg.lad[idx]
        dec_3d[ii, jj, kk] = veg.dec[idx]
        veg_index[ii, jj, kk] = idx

    return lad_3d, dec_3d, veg_index


def _require_numba() -> None:
    if nb is None:
        raise ImportError("numba is required for directshortwave_numba")


if nb is not None:

    @nb.njit(cache=True)
    def _point_to_cell_numba(
        x: float,
        y: float,
        z: float,
        dx: float,
        dy: float,
        z_edges: np.ndarray,
        z_max: float,
        itot: int,
        jtot: int,
        ktot: int,
        allow_outside_xy: bool,
    ) -> Tuple[int, int, int]:
        if z < 0.0 or z >= z_max:
            return -1, -1, -1
        if not allow_outside_xy:
            if x < 0.0 or x >= dx * itot:
                return -1, -1, -1
            if y < 0.0 or y >= dy * jtot:
                return -1, -1, -1
            i = int(min(itot - 1, max(0, math.floor(x / dx))))
            j = int(min(jtot - 1, max(0, math.floor(y / dy))))
        else:
            i = int(math.floor(x / dx))
            j = int(math.floor(y / dy))
        k = int(min(ktot - 1, max(0, np.searchsorted(z_edges, z, side="right") - 1)))
        return i, j, k

    @nb.njit(cache=True)
    def _ray_box_intersection_numba(
        origin: np.ndarray,
        direction: np.ndarray,
        bounds_min: np.ndarray,
        bounds_max: np.ndarray,
    ) -> Tuple[float, float]:
        t_min = -math.inf
        t_max = math.inf
        for a in range(3):
            if abs(direction[a]) < 1.0e-12:
                if origin[a] < bounds_min[a] or origin[a] > bounds_max[a]:
                    return math.inf, -math.inf
                continue
            inv_d = 1.0 / direction[a]
            t0 = (bounds_min[a] - origin[a]) * inv_d
            t1 = (bounds_max[a] - origin[a]) * inv_d
            if t0 > t1:
                t0, t1 = t1, t0
            t_min = max(t_min, t0)
            t_max = min(t_max, t1)
        return t_min, t_max

    @nb.njit(cache=True)
    def _trace_ray_numba(
        origin: np.ndarray,
        direction: np.ndarray,
        dx: float,
        dy: float,
        z_edges: np.ndarray,
        z_max: float,
        lad_3d: np.ndarray,
        dec_3d: np.ndarray,
        veg_index: np.ndarray,
        solid: np.ndarray,
        has_solid: bool,
        energy_in: np.ndarray,
        solid_hit_energy: np.ndarray,
        veg_absorb: np.ndarray,
        ray_area: float,
        itot: int,
        jtot: int,
        ktot: int,
        dz: np.ndarray,
        irradiance: float,
        periodic_xy: bool,
        max_ray_length: float,
        enable_hit_count: bool,
        hit_count: np.ndarray,
        allow_outside_xy: bool,
    ) -> float:
        x, y, z = origin

        i, j, k = _point_to_cell_numba(
            x,
            y,
            z,
            dx,
            dy,
            z_edges,
            z_max,
            itot,
            jtot,
            ktot,
            allow_outside_xy,
        )
        if i < 0:
            return 0.0

        dir_x, dir_y, dir_z = direction
        step_x = 0 if dir_x == 0.0 else (1 if dir_x > 0.0 else -1)
        step_y = 0 if dir_y == 0.0 else (1 if dir_y > 0.0 else -1)
        step_z = 0 if dir_z == 0.0 else (1 if dir_z > 0.0 else -1)

        if dir_x == 0.0:
            t_max_x = math.inf
            t_delta_x = math.inf
        else:
            if dir_x > 0.0:
                x_next = (i + 1) * dx
                t_max_x = (x_next - x) / dir_x
            else:
                x_prev = i * dx
                t_max_x = (x - x_prev) / (-dir_x)
            t_delta_x = dx / abs(dir_x)

        if dir_y == 0.0:
            t_max_y = math.inf
            t_delta_y = math.inf
        else:
            if dir_y > 0.0:
                y_next = (j + 1) * dy
                t_max_y = (y_next - y) / dir_y
            else:
                y_prev = j * dy
                t_max_y = (y - y_prev) / (-dir_y)
            t_delta_y = dy / abs(dir_y)

        if dir_z == 0.0:
            t_max_z = math.inf
            t_delta_z = math.inf
        else:
            if dir_z > 0.0:
                z_next = z_edges[k + 1]
                t_max_z = (z_next - z) / dir_z
            else:
                z_prev = z_edges[k]
                t_max_z = (z - z_prev) / (-dir_z)
            t_delta_z = dz[k] / abs(dir_z)

        t = 0.0
        r_in = 1.0
        last_i = -1
        last_j = -1
        last_k = -1

        while 0 <= k < ktot and (periodic_xy or allow_outside_xy or (0 <= i < itot and 0 <= j < jtot)):
            inside = 0 <= i < itot and 0 <= j < jtot
            if periodic_xy:
                ii = i % itot
                jj = j % jtot
                if has_solid and solid[ii, jj, k]:
                    if inside and last_i >= 0:
                        solid_hit_energy[last_i, last_j, last_k] += r_in * ray_area * irradiance
                        return 0.0
                    return r_in * ray_area * irradiance
            else:
                if inside:
                    ii = i
                    jj = j
                    if has_solid and solid[ii, jj, k]:
                        if last_i >= 0:
                            solid_hit_energy[last_i, last_j, last_k] += r_in * ray_area * irradiance
                            return 0.0
                        return r_in * ray_area * irradiance
                else:
                    ii = i
                    jj = j

            if enable_hit_count:
                if inside:
                    hit_count[ii, jj, k] += 1
            if inside:
                energy_in[ii, jj, k] += r_in * irradiance * ray_area

            t_next = min(t_max_x, t_max_y, t_max_z)
            if t_next > max_ray_length:
                ds = max_ray_length - t
            else:
                ds = t_next - t
            if ds < 0.0:
                ds = 0.0

            if inside:
                lad = lad_3d[ii, jj, k]
                dec = dec_3d[ii, jj, k]
                if lad > 0.0 and dec > 0.0:
                    tau = lad * dec * ds
                    r_out = r_in * math.exp(-tau)
                    absorbed = (r_in - r_out) * ray_area
                    vidx = veg_index[ii, jj, k]
                    if vidx >= 0:
                        cell_vol = dx * dy * dz[k]
                        veg_absorb[vidx] += absorbed / cell_vol
                    r_in = r_out

            if t_next > max_ray_length:
                return r_in * ray_area * irradiance

            prev_i, prev_j, prev_k = last_i, last_j, last_k
            if inside:
                prev_i, prev_j, prev_k = ii, jj, k
            if t_max_x == t_next:
                i += step_x
                t_max_x += t_delta_x
            if t_max_y == t_next:
                j += step_y
                t_max_y += t_delta_y
            if t_max_z == t_next:
                k += step_z
                if 0 <= k < ktot:
                    t_delta_z = dz[k] / abs(dir_z) if dir_z != 0.0 else math.inf
                    t_max_z = t_next + t_delta_z

            last_i, last_j, last_k = prev_i, prev_j, prev_k
            t = t_next

        if dir_z < 0.0 and k < 0 and last_i >= 0:
            solid_hit_energy[last_i, last_j, last_k] += r_in * ray_area * irradiance
            return 0.0
        return r_in * ray_area * irradiance


def directshortwave(
    sim,
    nsun: np.ndarray,
    irradiance: float,
    ray_scale: float = 1.0,
    periodic_xy: bool = False,
    max_ray_length: float | None = None,
    halo_pad: float = 0.125,
    ray_jitter: float = 0.0,
    ray_jitter_seed: int | None = None,
    return_hit_count: bool = False,
    extend_bounds: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    """
    Numba-accelerated direct shortwave on facets and absorbed radiation in vegetation cells.

    """
    if nb is None:
        raise ImportError("numba is required for directshortwave")

    if sim.geom is None or sim.geom.stl is None:
        raise ValueError("Geometry not loaded; UDBase must be created with load_geometry=True")

    nsun = np.asarray(nsun, dtype=float)
    norm = np.linalg.norm(nsun)
    if norm <= 0.0:
        raise ValueError("nsun must be non-zero")
    nsun_unit = nsun / norm
    direction = -nsun_unit
    if periodic_xy and abs(direction[2]) < 1.0e-12:
        raise ValueError("periodic_xy requires a non-zero vertical sun component")
    if extend_bounds and abs(direction[2]) < 1.0e-12:
        raise ValueError("extend_bounds requires a non-zero vertical sun component")

    ltree = getattr(sim, "ltree", 0)
    if ltree:
        veg = _load_veg_data(str(sim.path), sim.expnr)
    else:
        veg = VegData(
            points=np.empty((0, 3), dtype=int),
            lad=np.empty((0,), dtype=float),
            dec=np.empty((0,), dtype=float),
        )
    solid_full = getattr(sim, "Sc", None)
    kmax_solid = -1
    if solid_full is not None:
        solid_any = np.any(solid_full, axis=(0, 1))
        if np.any(solid_any):
            kmax_solid = int(np.max(np.where(solid_any)[0]))
    kmax_veg = int(np.max(veg.points[:, 2] - 1)) if veg.points.size else -1
    kmax = max(0, min(sim.ktot - 1, max(kmax_solid, kmax_veg)))
    ktot = min(sim.ktot, kmax + 2)

    lad_3d, dec_3d, veg_index = _build_veg_fields(sim, veg, ktot)

    z_edges_full = np.concatenate([sim.zm, [sim.zsize]])
    z_edges = z_edges_full[: ktot + 1]
    z_max = z_edges[-1]
    dz = sim.dzt[:ktot]
    energy_in = np.zeros((sim.itot, sim.jtot, ktot), dtype=float)
    solid_hit_energy = np.zeros((sim.itot, sim.jtot, ktot), dtype=float)
    veg_absorb = np.zeros(len(veg.points), dtype=float)
    solid = solid_full[:, :, :ktot] if solid_full is not None else np.zeros((1, 1, 1), dtype=bool)
    has_solid = solid_full is not None

    eps = 1.0e-6

    if extend_bounds:
        pad_x = z_max * abs(direction[0] / direction[2])
        pad_y = z_max * abs(direction[1] / direction[2])
        bounds_min = np.array([0.0, 0.0, 0.0], dtype=float)
        bounds_max = np.array([sim.xlen, sim.ylen, z_max], dtype=float)
        if direction[0] > 0.0:
            bounds_min[0] -= pad_x
        elif direction[0] < 0.0:
            bounds_max[0] += pad_x
        if direction[1] > 0.0:
            bounds_min[1] -= pad_y
        elif direction[1] < 0.0:
            bounds_max[1] += pad_y
    else:
        bounds_min = np.array([0.0, 0.0, 0.0], dtype=float)
        bounds_max = np.array([sim.xlen, sim.ylen, z_max], dtype=float)

    corners = np.array(
        [
            [bounds_min[0], bounds_min[1], bounds_min[2]],
            [bounds_max[0], bounds_min[1], bounds_min[2]],
            [bounds_min[0], bounds_max[1], bounds_min[2]],
            [bounds_max[0], bounds_max[1], bounds_min[2]],
            [bounds_min[0], bounds_min[1], bounds_max[2]],
            [bounds_max[0], bounds_min[1], bounds_max[2]],
            [bounds_min[0], bounds_max[1], bounds_max[2]],
            [bounds_max[0], bounds_max[1], bounds_max[2]],
        ],
        dtype=float,
    )
    dots = corners @ direction
    p0 = corners[np.argmin(dots)] - direction * (eps * 10.0)

    up = np.array([0.0, 0.0, 1.0], dtype=float)
    if abs(np.dot(up, nsun_unit)) > 0.95:
        up = np.array([0.0, 1.0, 0.0], dtype=float)
    u1 = np.cross(up, nsun_unit)
    u1 /= np.linalg.norm(u1)
    u2 = np.cross(nsun_unit, u1)

    proj = (corners - p0) @ np.vstack([u1, u2]).T
    umin, vmin = proj.min(axis=0)
    umax, vmax = proj.max(axis=0)
    du = umax - umin
    dv = vmax - vmin
    if halo_pad < 0.0:
        raise ValueError("halo_pad must be >= 0")
    umin -= halo_pad * du
    umax += halo_pad * du
    vmin -= halo_pad * dv
    vmax += halo_pad * dv

    if ray_scale <= 0.0:
        raise ValueError("ray_scale must be > 0")
    if max_ray_length is None:
        max_ray_length = 10.0 * max(sim.xlen, sim.ylen)
    if max_ray_length <= 0.0:
        raise ValueError("max_ray_length must be > 0")
    if ray_jitter < 0.0:
        raise ValueError("ray_jitter must be >= 0")

    step = (min(sim.dx, sim.dy, float(np.min(sim.dzt)))) / ray_scale
    u_vals = np.arange(umin, umax + step, step)
    v_vals = np.arange(vmin, vmax + step, step)
    ray_area = step * step
    rng = np.random.default_rng(ray_jitter_seed) if ray_jitter > 0.0 else None

    hit_count = (
        np.zeros((sim.itot, sim.jtot, ktot), dtype=np.int32)
        if return_hit_count
        else np.zeros((1, 1, 1), dtype=np.int32)
    )
    enable_hit_count = return_hit_count

    bud = {
        "in": 0.0,
        "veg": 0.0,
        "sol": 0.0,
        "out": 0.0,
        "fac": 0.0,
    }
    bud["rays"] = int(len(u_vals) * len(v_vals))

    for u in u_vals:
        for v in v_vals:
            if rng is not None:
                du_jit = (rng.random() - 0.5) * step * ray_jitter
                dv_jit = (rng.random() - 0.5) * step * ray_jitter
                u_use = u + du_jit
                v_use = v + dv_jit
            else:
                u_use = u
                v_use = v
            origin = p0 + u1 * u_use + u2 * v_use
            t0, t1 = _ray_box_intersection_numba(origin, direction, bounds_min, bounds_max)
            if t1 < t0:
                continue
            if t1 < 0.0:
                continue
            bud["in"] += irradiance * ray_area
            entry = max(t0, 0.0)
            start = origin + direction * (entry + eps)
            bud["out"] += _trace_ray_numba(
                start,
                direction,
                sim.dx,
                sim.dy,
                z_edges,
                z_max,
                lad_3d,
                dec_3d,
                veg_index,
                solid,
                has_solid,
                energy_in,
                solid_hit_energy,
                veg_absorb,
                ray_area,
                sim.itot,
                sim.jtot,
                ktot,
                dz,
                irradiance,
                periodic_xy,
                max_ray_length,
                enable_hit_count,
                hit_count,
                extend_bounds,
            )

    mesh = sim.geom.stl
    face_normals = mesh.face_normals
    nfaces = len(face_normals)

    sdir = np.zeros(nfaces, dtype=float)
    cos_inc_all = np.dot(face_normals, nsun_unit)
    cos_inc_all = np.where(cos_inc_all > 0.0, cos_inc_all, 0.0)

    if not hasattr(sim, "facsec") or sim.facsec is None or "c" not in sim.facsec:
        raise ValueError("Facet sections not available; sim.facsec['c'] is required.")
    facsec = sim.facsec["c"]
    facids = facsec["facid"].astype(int)
    areas = facsec["area"]
    locs = facsec["locs"].astype(int)
    i_idx = locs[:, 0]
    j_idx = locs[:, 1]
    k_idx = locs[:, 2]

    cell_area = np.zeros((sim.itot, sim.jtot, ktot), dtype=float)
    for idx in range(len(facids)):
        fid = facids[idx]
        cos_inc = cos_inc_all[fid]
        if cos_inc <= 0.0:
            continue
        i = i_idx[idx]
        j = j_idx[idx]
        k = k_idx[idx]
        cell_area[i, j, k] += areas[idx]

    sdir_accum = np.zeros(nfaces, dtype=float)
    area_accum = np.zeros(nfaces, dtype=float)
    for idx in range(len(facids)):
        fid = facids[idx]
        cos_inc = cos_inc_all[fid]
        if cos_inc <= 0.0:
            continue
        i = i_idx[idx]
        j = j_idx[idx]
        k = k_idx[idx]
        if cell_area[i, j, k] <= 0.0:
            continue
        cell_energy = solid_hit_energy[i, j, k]
        if cell_energy <= 0.0:
            continue
        sect_energy = cell_energy * (areas[idx] / cell_area[i, j, k])
        sdir_accum[fid] += sect_energy
        area_accum[fid] += areas[idx]

    mask = area_accum > 0.0
    areas = sim.facs["area"] if hasattr(sim, "facs") and "area" in sim.facs else mesh.area_faces
    sdir[mask] = sdir_accum[mask] / areas[mask]
    # WARNING: Temporary cap to avoid unrealistically large sdir on tiny facets when
    # an entire cell's energy lands on one face. Proper fix should distribute energy
    # by sub-facet geometry or ray/face intersections so flux per facet is bounded.
    sdir = np.minimum(sdir, irradiance * cos_inc_all)

    bud["fac"] = np.sum(sdir * areas)
    bud["sol"] = float(np.sum(solid_hit_energy))
    if veg_absorb.size:
        k_idx = veg.points[:, 2].astype(int) - 1
        cell_vol = sim.dx * sim.dy * dz[k_idx]
        bud["veg"] = float(np.sum(veg_absorb * cell_vol) * irradiance)
    if extend_bounds:
        unmapped = bud["sol"] - bud["fac"]
        if unmapped > 0.0:
            bud["out"] += unmapped
            bud["sol"] = bud["fac"]
    if return_hit_count:
        bud["hit_count"] = hit_count

    return sdir, veg_absorb, bud
