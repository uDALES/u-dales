from __future__ import annotations

"""
Direct shortwave solvers (facsec, moller, scanline) with shared helpers.

This module consolidates the implementations that previously lived in
directshortwave_facsec.py and directshortwave_moller.py so it can be used
standalone.
"""

from dataclasses import dataclass
from typing import Any, Dict, Tuple

import math
import numpy as np
import warnings

try:
    import numba as nb
except ImportError:  # pragma: no cover - runtime dependent
    nb = None


@dataclass
class VegData:
    points: np.ndarray  # (n, 3) 0-based (i, j, k)
    lad: np.ndarray     # (n,)
    dec: np.ndarray     # (n,)


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
        lad_3d[i, j, k] = veg.lad[idx]
        dec_3d[i, j, k] = veg.dec[idx]
        veg_index[i, j, k] = idx

    return lad_3d, dec_3d, veg_index


def _veg_from_data(veg_data: Dict[str, Any] | None) -> VegData:
    if veg_data is None:
        return VegData(
            points=np.empty((0, 3), dtype=int),
            lad=np.empty((0,), dtype=float),
            dec=np.empty((0,), dtype=float),
        )
    points = np.asarray(veg_data.get("points", []), dtype=int)
    if points.ndim == 1 and points.size:
        points = points.reshape(1, -1)
    params = veg_data.get("params", {})
    lad = np.asarray(params.get("lad", []), dtype=float)
    dec = np.asarray(params.get("dec", []), dtype=float)
    if points.size == 0:
        raise ValueError("veg_data contains no vegetation points")
    if len(points) != len(lad) or len(points) != len(dec):
        raise ValueError("veg_data points and params lengths do not match")
    return VegData(points=points, lad=lad, dec=dec)


def _compute_ktot_and_z_edges(
    grid,
    veg_points: np.ndarray | None = None,
) -> Tuple[int, np.ndarray, float, np.ndarray]:
    solid_full = getattr(grid, "Sc", None)
    kmax_solid = -1
    if solid_full is not None:
        solid_any = np.any(solid_full, axis=(0, 1))
        if np.any(solid_any):
            kmax_solid = int(np.max(np.where(solid_any)[0]))
    kmax_veg = int(np.max(veg_points[:, 2])) if veg_points is not None and veg_points.size else -1
    kmax = max(0, min(grid.ktot - 1, max(kmax_solid, kmax_veg)))
    ktot = min(grid.ktot, kmax + 2)
    z_base = getattr(grid, "zf", None)
    if z_base is None or len(z_base) == 0:
        z_base = grid.zm
    z_edges_full = np.concatenate([z_base, [grid.zsize]])
    z_edges = z_edges_full[: ktot + 1]
    z_max = float(z_edges[-1])
    dz = grid.dzt[:ktot]
    return ktot, z_edges, z_max, dz


def _build_cell_facet_lookup(
    triangles: np.ndarray,
    dx: float,
    dy: float,
    z_edges: np.ndarray,
    itot: int,
    jtot: int,
    ktot: int,
    face_mask: np.ndarray | None = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a cell->facet CSR lookup based on triangle AABB overlap."""
    # Build a per-cell list of candidate facets by overlapping triangle AABBs with grid cells.
    cell_idx_list = []
    facet_idx_list = []
    eps = 1.0e-9 * max(dx, dy)
    for fid, tri in enumerate(triangles):
        if face_mask is not None and not face_mask[fid]:
            continue
        vmin = tri.min(axis=0) - eps
        vmax = tri.max(axis=0) + eps
        ix0 = int(math.floor(vmin[0] / dx))
        ix1 = int(math.floor(vmax[0] / dx))
        iy0 = int(math.floor(vmin[1] / dy))
        iy1 = int(math.floor(vmax[1] / dy))
        k0 = int(np.searchsorted(z_edges, vmin[2], side="right") - 1)
        k1 = int(np.searchsorted(z_edges, vmax[2], side="right") - 1)
        if ix1 < 0 or iy1 < 0 or k1 < 0:
            continue
        if ix0 >= itot or iy0 >= jtot or k0 >= ktot:
            continue
        ix0 = max(ix0, 0)
        iy0 = max(iy0, 0)
        k0 = max(k0, 0)
        ix1 = min(ix1, itot - 1)
        iy1 = min(iy1, jtot - 1)
        k1 = min(k1, ktot - 1)
        for i in range(ix0, ix1 + 1):
            for j in range(iy0, iy1 + 1):
                for k in range(k0, k1 + 1):
                    cell_idx_list.append(i + itot * (j + jtot * k))
                    facet_idx_list.append(fid)
    if not cell_idx_list:
        n_cells = itot * jtot * ktot
        cell_offsets = np.zeros(n_cells + 1, dtype=np.int64)
        cell_facets = np.zeros(0, dtype=np.int32)
        cell_has_facets = np.zeros((itot, jtot, ktot), dtype=bool)
        return cell_offsets, cell_facets, cell_has_facets
    cell_idx = np.asarray(cell_idx_list, dtype=np.int64)
    facet_idx = np.asarray(facet_idx_list, dtype=np.int32)
    # Sort and de-duplicate (cell, facet) pairs, then build a CSR lookup.
    order = np.lexsort((facet_idx, cell_idx))
    cell_idx = cell_idx[order]
    facet_idx = facet_idx[order]
    uniq = np.ones(len(cell_idx), dtype=bool)
    uniq[1:] = (cell_idx[1:] != cell_idx[:-1]) | (facet_idx[1:] != facet_idx[:-1])
    cell_idx = cell_idx[uniq]
    facet_idx = facet_idx[uniq]
    n_cells = itot * jtot * ktot
    counts = np.bincount(cell_idx, minlength=n_cells).astype(np.int64)
    cell_offsets = np.zeros(n_cells + 1, dtype=np.int64)
    cell_offsets[1:] = np.cumsum(counts, dtype=np.int64)
    cell_facets = facet_idx.astype(np.int32)
    cell_has_facets = counts.reshape((itot, jtot, ktot), order="F") > 0
    return cell_offsets, cell_facets, cell_has_facets


def _require_numba() -> None:
    if nb is None:
        raise ImportError("numba is required for directshortwave")


if nb is not None:

    @nb.njit(cache=True)
    def _ray_triangle_intersect_numba(
        ox: float,
        oy: float,
        oz: float,
        dx: float,
        dy: float,
        dz: float,
        v0: np.ndarray,
        v1: np.ndarray,
        v2: np.ndarray,
        t_max: float,
    ) -> Tuple[bool, float]:
        eps = 1.0e-12
        e1x = v1[0] - v0[0]
        e1y = v1[1] - v0[1]
        e1z = v1[2] - v0[2]
        e2x = v2[0] - v0[0]
        e2y = v2[1] - v0[1]
        e2z = v2[2] - v0[2]
        pvx = dy * e2z - dz * e2y
        pvy = dz * e2x - dx * e2z
        pvz = dx * e2y - dy * e2x
        det = e1x * pvx + e1y * pvy + e1z * pvz
        if -eps < det < eps:
            return False, 0.0
        inv_det = 1.0 / det
        tvx = ox - v0[0]
        tvy = oy - v0[1]
        tvz = oz - v0[2]
        u = (tvx * pvx + tvy * pvy + tvz * pvz) * inv_det
        if u < 0.0 or u > 1.0:
            return False, 0.0
        qvx = tvy * e1z - tvz * e1y
        qvy = tvz * e1x - tvx * e1z
        qvz = tvx * e1y - tvy * e1x
        v = (dx * qvx + dy * qvy + dz * qvz) * inv_det
        if v < 0.0 or u + v > 1.0:
            return False, 0.0
        t = (e2x * qvx + e2y * qvy + e2z * qvz) * inv_det
        if t < -1.0e-6 or t > t_max + 1.0e-6:
            return False, 0.0
        return True, t

    @nb.njit(cache=True)
    def _ray_cell_intersect_numba(
        ox: float,
        oy: float,
        oz: float,
        dx: float,
        dy: float,
        dz: float,
        triangles: np.ndarray,
        cell_offsets: np.ndarray,
        cell_facets: np.ndarray,
        cell_idx: int,
        t_max: float,
        debug_facid: int,
        debug_min_t: np.ndarray,
        debug_count: np.ndarray,
        debug_test_count: np.ndarray,
    ) -> Tuple[float, int]:
        start = cell_offsets[cell_idx]
        end = cell_offsets[cell_idx + 1]
        if start == end:
            return -1.0, -1
        best = t_max + 1.0
        best_fid = -1
        for idx in range(start, end):
            fid = cell_facets[idx]
            if debug_facid >= 0 and fid == debug_facid:
                debug_test_count[0] += 1
            v0 = triangles[fid, 0]
            v1 = triangles[fid, 1]
            v2 = triangles[fid, 2]
            hit, t = _ray_triangle_intersect_numba(ox, oy, oz, dx, dy, dz, v0, v1, v2, t_max)
            if hit and t < best:
                best = t
                best_fid = fid
            if debug_facid >= 0 and fid == debug_facid and hit:
                if t < debug_min_t[0]:
                    debug_min_t[0] = t
                debug_count[0] += 1
        if best <= t_max:
            return best, best_fid
        return -1.0, -1

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
    def _trace_ray_facsec_numba(
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

    @nb.njit(cache=True)
    def _trace_ray_moller_numba(
        origin: np.ndarray,
        direction: np.ndarray,
        dx: float,
        dy: float,
        z_edges: np.ndarray,
        z_max: float,
        lad_3d: np.ndarray,
        dec_3d: np.ndarray,
        veg_index: np.ndarray,
        cell_has_facets: np.ndarray,
        cell_offsets: np.ndarray,
        cell_facets: np.ndarray,
        triangles: np.ndarray,
        facet_hit_energy: np.ndarray,
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
        debug_facid: int,
        debug_min_t: np.ndarray,
        debug_count: np.ndarray,
        debug_test_count: np.ndarray,
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
            base_inside = 0 <= i < itot and 0 <= j < jtot
            geom_inside = base_inside or periodic_xy
            if periodic_xy:
                ii = i % itot
                jj = j % jtot
            else:
                ii = i
                jj = j
            if enable_hit_count and base_inside:
                hit_count[ii, jj, k] += 1
            if base_inside:
                energy_in[ii, jj, k] += r_in * irradiance * ray_area

            t_next = min(t_max_x, t_max_y, t_max_z)
            t_limit = t_next
            if t_limit > max_ray_length:
                t_limit = max_ray_length
            ds = t_limit - t
            if ds < 0.0:
                ds = 0.0
            if dir_z < 0.0 and k == 0:
                ds += 1.0e-6

            hit_facet = False
            if geom_inside and cell_has_facets[ii, jj, k] and ds > 0.0:
                cell_idx = ii + itot * (jj + jtot * k)
                ox = x + dir_x * t
                oy = y + dir_y * t
                oz = z + dir_z * t
                if periodic_xy:
                    xlen = dx * itot
                    ylen = dy * jtot
                    ox = ox - math.floor(ox / xlen) * xlen
                    oy = oy - math.floor(oy / ylen) * ylen
                t_hit, hit_fid = _ray_cell_intersect_numba(
                    ox,
                    oy,
                    oz,
                    dir_x,
                    dir_y,
                    dir_z,
                    triangles,
                    cell_offsets,
                    cell_facets,
                    cell_idx,
                    ds,
                    debug_facid,
                    debug_min_t,
                    debug_count,
                    debug_test_count,
                )
                if t_hit >= 0.0:
                    ds = max(0.0, t_hit)
                    if hit_fid >= 0:
                        if base_inside:
                            facet_hit_energy[hit_fid] += r_in * ray_area * irradiance
                    hit_facet = True

            if base_inside:
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

            if hit_facet:
                if base_inside:
                    solid_hit_energy[ii, jj, k] += r_in * ray_area * irradiance
                    return 0.0
                return r_in * ray_area * irradiance

            if t_next > max_ray_length:
                return r_in * ray_area * irradiance

            prev_i, prev_j, prev_k = last_i, last_j, last_k
            if base_inside:
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

    @nb.njit(cache=True)
    def _trace_ray_segments_numba(
        origin: np.ndarray,
        direction: np.ndarray,
        dx: float,
        dy: float,
        z_edges: np.ndarray,
        z_max: float,
        itot: int,
        jtot: int,
        ktot: int,
        dz: np.ndarray,
        max_ray_length: float,
        allow_outside_xy: bool,
        max_steps: int,
    ) -> Tuple[np.ndarray, int]:
        buf = np.zeros((max_steps, 10), dtype=np.float64)
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
            return buf, 0
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
        idx = 0
        while 0 <= k < ktot and (allow_outside_xy or (0 <= i < itot and 0 <= j < jtot)):
            if idx >= max_steps:
                break
            inside = 1.0 if (0 <= i < itot and 0 <= j < jtot) else 0.0
            t_next = min(t_max_x, t_max_y, t_max_z)
            buf[idx, 0] = float(i)
            buf[idx, 1] = float(j)
            buf[idx, 2] = float(k)
            buf[idx, 3] = t
            buf[idx, 4] = t_next
            buf[idx, 5] = t_max_x
            buf[idx, 6] = t_max_y
            buf[idx, 7] = t_max_z
            buf[idx, 8] = t_next - t
            buf[idx, 9] = inside
            idx += 1
            if t_next > max_ray_length:
                break
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
            t = t_next
        return buf, idx

    @nb.njit(cache=True)
    def _cast_rays_numba(
        u_vals: np.ndarray,
        v_vals: np.ndarray,
        jitter_u: np.ndarray,
        jitter_v: np.ndarray,
        use_jitter: bool,
        p0: np.ndarray,
        u1: np.ndarray,
        u2: np.ndarray,
        direction: np.ndarray,
        bounds_min: np.ndarray,
        bounds_max: np.ndarray,
        dx: float,
        dy: float,
        z_edges: np.ndarray,
        z_max: float,
        lad_3d: np.ndarray,
        dec_3d: np.ndarray,
        veg_index: np.ndarray,
        cell_has_facets: np.ndarray,
        cell_offsets: np.ndarray,
        cell_facets: np.ndarray,
        triangles: np.ndarray,
        facet_hit_energy: np.ndarray,
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
        debug_facid: int,
        debug_min_t: np.ndarray,
        debug_count: np.ndarray,
        debug_test_count: np.ndarray,
        bud_in: np.ndarray,
        bud_out: np.ndarray,
    ) -> None:
        nu = u_vals.shape[0]
        nv = v_vals.shape[0]
        idx = 0
        for iu in range(nu):
            u = u_vals[iu]
            for iv in range(nv):
                v = v_vals[iv]
                if use_jitter:
                    u_use = u + jitter_u[idx]
                    v_use = v + jitter_v[idx]
                else:
                    u_use = u
                    v_use = v
                idx += 1
                origin = p0 + u1 * u_use + u2 * v_use
                t0, t1 = _ray_box_intersection_numba(origin, direction, bounds_min, bounds_max)
                if t1 < t0:
                    continue
                if t1 < 0.0:
                    continue
                bud_in[0] += irradiance * ray_area
                entry = t0 if t0 > 0.0 else 0.0
                start = origin + direction * (entry + 1.0e-6)
                bud_out[0] += _trace_ray_moller_numba(
                    start,
                    direction,
                    dx,
                    dy,
                    z_edges,
                    z_max,
                    lad_3d,
                    dec_3d,
                    veg_index,
                    cell_has_facets,
                    cell_offsets,
                    cell_facets,
                    triangles,
                    facet_hit_energy,
                    energy_in,
                    solid_hit_energy,
                    veg_absorb,
                    ray_area,
                    itot,
                    jtot,
                    ktot,
                    dz,
                    irradiance,
                    periodic_xy,
                    max_ray_length,
                    enable_hit_count,
                    hit_count,
                    allow_outside_xy,
                    debug_facid,
                    debug_min_t,
                    debug_count,
                    debug_test_count,
                )


class DirectShortwaveSolver:
    """
    Direct shortwave solver with cached geometry/vegetation and per-call sun setup.

    Models
    ------
    method : {"moller", "facsec", "scanline"}
        - "moller": DDA ray casting with Moller-Trumbore triangle hits
          (most accurate, most expensive; supports vegetation + periodic_xy).
        - "facsec": DDA ray casting with solid mask + facet-section reconstruction
          (accurate, faster; supports vegetation + periodic_xy).
        - "scanline": f2py scanline rasterization on the surface mesh
          (fastest; no vegetation support, no periodicity).

    Constructor options
    -------------------
    sim : UDBase-like
        Provides grid geometry (dx, dy, dzt, itot, jtot, ktot, zf/zm, xlen, ylen, zsize)
        and optional fields (Sc, facsec, facs["area"]).
    surface_mesh : trimesh-like, optional
        Triangle mesh with vertices/faces/face_normals. Defaults to sim.geom.stl.
    ray_density : float, default 4.0
        Ray sampling density on the ray plane. Spacing = min(dx, dy, dz) / ray_density.
    ray_jitter : float, default 1.0
        Random jitter factor in [0, 1] of the ray spacing.
    veg_data : dict, optional
        Vegetation data from UDBase.load_veg(zero_based=True).

    Compute options
    ---------------
    compute(nsun, irradiance, periodic_xy=False, resolution=None)
        - nsun: unit or non-unit sun vector.
        - irradiance: direct normal irradiance [W/m^2].
        - periodic_xy: wrap rays in x/y to mimic infinite tiling.
        - resolution: scanline pixel size override (scanline only).
    """

    def __init__(
        self,
        sim,
        method: str,
        *,
        surface_mesh=None,
        ray_density: float = 4.0,
        ray_jitter: float = 1.0,
        veg_data: Dict[str, Any] | None = None,
    ) -> None:
        self.sim = sim
        self.method = method.strip().lower()
        self.ray_density = float(ray_density)
        self.ray_jitter = float(ray_jitter)
        self.veg = _veg_from_data(veg_data)

        if self.ray_density <= 0.0:
            raise ValueError("ray_density must be > 0")
        if self.ray_jitter < 0.0:
            raise ValueError("ray_jitter must be >= 0")
        if surface_mesh is None:
            if sim is None or getattr(sim, "geom", None) is None or sim.geom.stl is None:
                raise ValueError("surface_mesh is required (sim.geom.stl not available)")
            surface_mesh = sim.geom.stl

        if self.method in ("moller", "facsec"):
            if nb is None:
                raise ImportError("numba is required for direct shortwave (moller/facsec)")
        elif self.method == "scanline":
            try:
                from udprep.directshortwave_f2py import directshortwave_f2py_mod as _dsmod
            except ImportError as exc:
                raise RuntimeError(
                    "directshortwave_f2py module not available; "
                    "build it with tools/python/fortran/build_f2py.ps1"
                ) from exc
            self._dsmod = _dsmod
            if self.veg.points.size:
                warnings.warn(
                    "Scanline (f2py) direct shortwave does not include vegetation; "
                    "results ignore tree attenuation.",
                    RuntimeWarning,
                )
                self.veg = VegData(
                    points=np.empty((0, 3), dtype=int),
                    lad=np.empty((0,), dtype=float),
                    dec=np.empty((0,), dtype=float),
                )
        else:
            raise ValueError(f"Unknown direct shortwave method: {method}")

        self.ktot, self.z_edges, self.z_max, self.dz = _compute_ktot_and_z_edges(
            sim, self.veg.points if self.veg.points.size else None
        )
        self.lad_3d, self.dec_3d, self.veg_index = _build_veg_fields(
            sim, self.veg, self.ktot
        )

        self.mesh = surface_mesh
        mesh = surface_mesh
        self.face_normals = mesh.face_normals
        self.nfaces = len(self.face_normals)
        self.face_areas = mesh.area_faces
        if hasattr(sim, "facs") and "area" in sim.facs and sim.facs["area"].shape == self.face_areas.shape:
            self.face_areas = sim.facs["area"]

        if self.method == "moller":
            if hasattr(mesh, "triangles"):
                self.triangles = np.asarray(mesh.triangles, dtype=float)
            else:
                self.triangles = np.asarray(mesh.vertices[mesh.faces], dtype=float)
            (
                self.cell_offsets,
                self.cell_facets,
                self.cell_has_facets,
            ) = _build_cell_facet_lookup(
                self.triangles,
                sim.dx,
                sim.dy,
                self.z_edges,
                sim.itot,
                sim.jtot,
                self.ktot,
            )

        if self.method == "facsec":
            solid_full = getattr(sim, "Sc", None)
            self.has_solid = solid_full is not None
            self.solid = solid_full[:, :, : self.ktot] if self.has_solid else np.zeros((1, 1, 1), dtype=bool)

            if not hasattr(sim, "facsec") or sim.facsec is None or "c" not in sim.facsec:
                raise ValueError("Facet sections not available; sim.facsec['c'] is required.")
            facsec = sim.facsec["c"]
            self.facsec_facids = facsec["facid"].astype(int)
            self.facsec_areas = facsec["area"]
            self.facsec_locs = facsec["locs"].astype(int)

    def compute(
        self,
        nsun: np.ndarray,
        irradiance: float,
        *,
        periodic_xy: bool = False,
        resolution: float | None = None,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        nsun = np.asarray(nsun, dtype=float)
        norm = np.linalg.norm(nsun)
        if norm <= 0.0:
            raise ValueError("nsun must be non-zero")
        nsun_unit = nsun / norm
        direction = -nsun_unit
        if abs(direction[2]) < 1.0e-12:
            raise ValueError("direct shortwave requires a non-zero vertical sun component")

        if self.method == "scanline":
            return self._compute_scanline(nsun_unit, irradiance, resolution=resolution)
        if self.method == "moller":
            return self._compute_moller(nsun_unit, irradiance, direction, periodic_xy)
        if self.method == "facsec":
            return self._compute_facsec(nsun_unit, irradiance, direction, periodic_xy)
        raise ValueError(f"Unknown direct shortwave method: {self.method}")

    def _compute_bounds(self, direction: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        # Expand the domain bounds along the sun direction so rays start outside the scene.
        pad_x = self.z_max * abs(direction[0] / direction[2])
        pad_y = self.z_max * abs(direction[1] / direction[2])
        bounds_min = np.array([0.0, 0.0, 0.0], dtype=float)
        bounds_max = np.array([self.sim.xlen, self.sim.ylen, self.z_max], dtype=float)
        if direction[0] > 0.0:
            bounds_min[0] -= pad_x
        elif direction[0] < 0.0:
            bounds_max[0] += pad_x
        if direction[1] > 0.0:
            bounds_min[1] -= pad_y
        elif direction[1] < 0.0:
            bounds_max[1] += pad_y
        return bounds_min, bounds_max

    def _compute_ray_plane(self, direction: np.ndarray, nsun_unit: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float, float]:
        # Build a sampling plane orthogonal to the sun direction that fully covers the scene.
        bounds_min, bounds_max = self._compute_bounds(direction)
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
        eps = 1.0e-6
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
        return p0, u1, u2, float(umin), float(umax), float(vmin), float(vmax)

    def _build_ray_grid(self, umin: float, umax: float, vmin: float, vmax: float) -> Tuple[np.ndarray, np.ndarray, float]:
        # Uniform ray grid on the sampling plane; density is controlled by ray_density.
        step = (min(self.sim.dx, self.sim.dy, float(np.min(self.sim.dzt)))) / self.ray_density
        u_vals = np.arange(umin, umax + step, step)
        v_vals = np.arange(vmin, vmax + step, step)
        return u_vals, v_vals, step

    def _compute_moller(
        self,
        nsun_unit: np.ndarray,
        irradiance: float,
        direction: np.ndarray,
        periodic_xy: bool,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        # Ray-plane definition + grid on that plane.
        p0, u1, u2, umin, umax, vmin, vmax = self._compute_ray_plane(direction, nsun_unit)
        u_vals, v_vals, step = self._build_ray_grid(umin, umax, vmin, vmax)
        ray_area = step * step
        use_jitter = self.ray_jitter > 0.0

        bud = {"in": 0.0, "veg": 0.0, "sol": 0.0, "out": 0.0, "fac": 0.0}
        bud["rays"] = int(len(u_vals) * len(v_vals))

        if use_jitter:
            rng = np.random.default_rng(0)
            jitter_u = (rng.random(bud["rays"]) - 0.5) * step * self.ray_jitter
            jitter_v = (rng.random(bud["rays"]) - 0.5) * step * self.ray_jitter
        else:
            jitter_u = np.zeros(bud["rays"], dtype=float)
            jitter_v = np.zeros(bud["rays"], dtype=float)

        energy_in = np.zeros((self.sim.itot, self.sim.jtot, self.ktot), dtype=float)
        solid_hit_energy = np.zeros((self.sim.itot, self.sim.jtot, self.ktot), dtype=float)
        veg_absorb = np.zeros(len(self.veg.points), dtype=float)
        facet_hit_energy = np.zeros(self.nfaces, dtype=float)

        bud_in = np.zeros(1, dtype=float)
        bud_out = np.zeros(1, dtype=float)

        # Cast all rays with triangle intersection tests; accumulate energy in arrays.
        _cast_rays_numba(
            u_vals,
            v_vals,
            jitter_u,
            jitter_v,
            use_jitter,
            p0,
            u1,
            u2,
            direction,
            *self._compute_bounds(direction),
            self.sim.dx,
            self.sim.dy,
            self.z_edges,
            self.z_max,
            self.lad_3d,
            self.dec_3d,
            self.veg_index,
            self.cell_has_facets,
            self.cell_offsets,
            self.cell_facets,
            self.triangles,
            facet_hit_energy,
            energy_in,
            solid_hit_energy,
            veg_absorb,
            ray_area,
            self.sim.itot,
            self.sim.jtot,
            self.ktot,
            self.dz,
            irradiance,
            periodic_xy,
            10.0 * max(self.sim.xlen, self.sim.ylen),
            False,
            np.zeros((1, 1, 1), dtype=np.int32),
            True,
            -1,
            np.array([np.inf], dtype=float),
            np.array([0], dtype=np.int64),
            np.array([0], dtype=np.int64),
            bud_in,
            bud_out,
        )

        bud["in"] = float(bud_in[0])
        bud["out"] = float(bud_out[0])

        cos_inc_all = np.dot(self.face_normals, nsun_unit)
        cos_inc_all = np.where(cos_inc_all > 0.0, cos_inc_all, 0.0)

        # Convert accumulated facet energy to flux (W/m^2) with a cosine cap.
        sdir = np.zeros(self.nfaces, dtype=float)
        mask = self.face_areas > 0.0
        sdir[mask] = facet_hit_energy[mask] / self.face_areas[mask]
        sdir = np.minimum(sdir, irradiance * cos_inc_all)

        bud["fac"] = float(np.sum(sdir * self.face_areas))
        bud["sol"] = float(np.sum(solid_hit_energy))
        if veg_absorb.size:
            k_idx = self.veg.points[:, 2].astype(int)
            cell_vol = self.sim.dx * self.sim.dy * self.dz[k_idx]
            bud["veg"] = float(np.sum(veg_absorb * cell_vol) * irradiance)
        # Any solid energy not mapped to facets is treated as escaping the domain.
        unmapped = bud["sol"] - bud["fac"]
        if unmapped > 0.0:
            bud["out"] += unmapped
            bud["sol"] = bud["fac"]
        return sdir, veg_absorb, bud

    def _compute_facsec(
        self,
        nsun_unit: np.ndarray,
        irradiance: float,
        direction: np.ndarray,
        periodic_xy: bool,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        # Ray-plane definition + grid on that plane.
        p0, u1, u2, umin, umax, vmin, vmax = self._compute_ray_plane(direction, nsun_unit)
        u_vals, v_vals, step = self._build_ray_grid(umin, umax, vmin, vmax)
        ray_area = step * step
        rng = np.random.default_rng(0) if self.ray_jitter > 0.0 else None

        energy_in = np.zeros((self.sim.itot, self.sim.jtot, self.ktot), dtype=float)
        solid_hit_energy = np.zeros((self.sim.itot, self.sim.jtot, self.ktot), dtype=float)
        veg_absorb = np.zeros(len(self.veg.points), dtype=float)

        bud = {"in": 0.0, "veg": 0.0, "sol": 0.0, "out": 0.0, "fac": 0.0}
        bud["rays"] = int(len(u_vals) * len(v_vals))

        # Trace rays through the solid mask; energy is accumulated per grid cell.
        for u in u_vals:
            for v in v_vals:
                if rng is not None:
                    du_jit = (rng.random() - 0.5) * step * self.ray_jitter
                    dv_jit = (rng.random() - 0.5) * step * self.ray_jitter
                    u_use = u + du_jit
                    v_use = v + dv_jit
                else:
                    u_use = u
                    v_use = v
                origin = p0 + u1 * u_use + u2 * v_use
                t0, t1 = _ray_box_intersection_numba(origin, direction, *self._compute_bounds(direction))
                if t1 < t0:
                    continue
                if t1 < 0.0:
                    continue
                bud["in"] += irradiance * ray_area
                entry = max(t0, 0.0)
                start = origin + direction * (entry + 1.0e-6)
                bud["out"] += _trace_ray_facsec_numba(
                    start,
                    direction,
                    self.sim.dx,
                    self.sim.dy,
                    self.z_edges,
                    self.z_max,
                    self.lad_3d,
                    self.dec_3d,
                    self.veg_index,
                    self.solid,
                    self.has_solid,
                    energy_in,
                    solid_hit_energy,
                    veg_absorb,
                    ray_area,
                    self.sim.itot,
                    self.sim.jtot,
                    self.ktot,
                    self.dz,
                    irradiance,
                    periodic_xy,
                    10.0 * max(self.sim.xlen, self.sim.ylen),
                    False,
                    np.zeros((1, 1, 1), dtype=np.int32),
                    True,
                )

        cos_inc_all = np.dot(self.face_normals, nsun_unit)
        cos_inc_all = np.where(cos_inc_all > 0.0, cos_inc_all, 0.0)

        facids = self.facsec_facids
        areas = self.facsec_areas
        locs = self.facsec_locs
        i_idx = locs[:, 0]
        j_idx = locs[:, 1]
        k_idx = locs[:, 2]

        # Precompute projected face area per cell for facsec energy redistribution.
        cell_proj_area = np.zeros((self.sim.itot, self.sim.jtot, self.ktot), dtype=float)
        # Redistribute cell energy to facets proportionally by projected area.
        for idx in range(len(facids)):
            fid = facids[idx]
            cos_inc = cos_inc_all[fid]
            if cos_inc <= 0.0:
                continue
            i = i_idx[idx]
            j = j_idx[idx]
            k = k_idx[idx]
            cell_proj_area[i, j, k] += areas[idx] * cos_inc

        sdir_accum = np.zeros(self.nfaces, dtype=float)
        area_accum = np.zeros(self.nfaces, dtype=float)
        for idx in range(len(facids)):
            fid = facids[idx]
            cos_inc = cos_inc_all[fid]
            if cos_inc <= 0.0:
                continue
            i = i_idx[idx]
            j = j_idx[idx]
            k = k_idx[idx]
            if cell_proj_area[i, j, k] <= 0.0:
                continue
            cell_energy = solid_hit_energy[i, j, k]
            if cell_energy <= 0.0:
                continue
            proj_area = areas[idx] * cos_inc
            sect_energy = cell_energy * (proj_area / cell_proj_area[i, j, k])
            sdir_accum[fid] += sect_energy
            area_accum[fid] += areas[idx]

        sdir = np.zeros(self.nfaces, dtype=float)
        mask = area_accum > 0.0
        sdir[mask] = sdir_accum[mask] / self.face_areas[mask]
        sdir = np.minimum(sdir, irradiance * cos_inc_all)

        bud["fac"] = float(np.sum(sdir * self.face_areas))
        bud["sol"] = float(np.sum(solid_hit_energy))
        if veg_absorb.size:
            k_idx = self.veg.points[:, 2].astype(int)
            cell_vol = self.sim.dx * self.sim.dy * self.dz[k_idx]
            bud["veg"] = float(np.sum(veg_absorb * cell_vol) * irradiance)
        unmapped = bud["sol"] - bud["fac"]
        if unmapped > 0.0:
            bud["out"] += unmapped
            bud["sol"] = bud["fac"]
        return sdir, veg_absorb, bud

    def _compute_scanline(
        self,
        nsun_unit: np.ndarray,
        irradiance: float,
        *,
        resolution: float | None = None,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        mesh = self.mesh
        vertices = np.asfortranarray(mesh.vertices, dtype=float)
        faces = np.asfortranarray(mesh.faces, dtype=np.int32) + 1  # 1-based
        incenter = np.asfortranarray(mesh.triangles_center, dtype=float)
        face_normal = np.asfortranarray(mesh.face_normals, dtype=float)
        if resolution is None:
            cell_min = min(self.sim.dx, self.sim.dy, float(np.min(self.sim.dzt)))
            resolution = 0.25 * cell_min / self.ray_density
        if resolution <= 0.0:
            raise ValueError("resolution must be > 0")

        sdir = self._dsmod.calculate_direct_shortwave(
            faces,
            incenter,
            face_normal,
            vertices,
            nsun_unit,
            float(irradiance),
            float(resolution),
        )
        sdir = np.asarray(sdir, dtype=float)
        bud = {
            "fac": float(np.sum(sdir * self.face_areas)),
            "veg": 0.0,
        }
        return sdir, np.zeros(0, dtype=float), bud

    def trace_ray_segments(
        self,
        origin: np.ndarray,
        direction: np.ndarray,
        max_steps: int = 2000,
        allow_outside_xy: bool = True,
    ) -> np.ndarray:
        """Trace a single ray through the DDA and return segment diagnostics."""
        if nb is None:
            raise ImportError("numba is required for trace_ray_segments")
        ktot, z_edges, z_max, dz = _compute_ktot_and_z_edges(self.sim, None)
        max_ray_length = 10.0 * max(self.sim.xlen, self.sim.ylen)
        buf, count = _trace_ray_segments_numba(
            np.asarray(origin, dtype=float),
            np.asarray(direction, dtype=float),
            self.sim.dx,
            self.sim.dy,
            z_edges,
            z_max,
            self.sim.itot,
            self.sim.jtot,
            ktot,
            dz,
            max_ray_length,
            allow_outside_xy,
            max_steps,
        )
        return buf[:count]
