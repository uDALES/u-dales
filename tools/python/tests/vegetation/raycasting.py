from __future__ import annotations
import time

start = time.perf_counter()

import sys
from pathlib import Path

sys.path.insert(0, "/rds/general/user/mvr/home/udales/u-dales/tools/python/")
sys.path.insert(0, r"C:\Users\mvr\OneDrive - Imperial College London\codes\uDALES\u-dales\tools\python")

sim_id = 525
base_path = Path(r"C:\Users\mvr\OneDrive - Imperial College London\codes\uDALES\u-dales\tests\tests_tree_input")

from udbase import UDBase  # noqa: E402
from udprep import convert_block_to_sparse, directshortwave as directshortwave_py  # noqa: E402
from udprep.directshortwave_numba import directshortwave as directshortwave_nb  # noqa: E402

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

elapsed = time.perf_counter() - start
print(f"loading libraries runtime: {elapsed:.3f} s")

#% ----------------------------------------------------------------

start = time.perf_counter()

# Instantiate UDBase with geometry for direct shortwave
sim = UDBase(expnr=sim_id, path=base_path, load_geometry=True)
if sim.trees is None:
    raise RuntimeError(f"trees.inp.{sim_id} not loaded from {base_path}")

elapsed = time.perf_counter() - start
print(f"UDbase startup runtime: {elapsed:.3f} s")

azimuth_deg = 20.0
elevation_deg = 15.0
az = np.deg2rad(azimuth_deg)
el = np.deg2rad(elevation_deg)
nsun = [np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)]
irradiance = 800.0
extend_bounds = True

def _ray_debug_info(sim, nsun_vec, halo_pad=0.125):
    ns = np.asarray(nsun_vec, dtype=float)
    ns /= np.linalg.norm(ns)
    direction = -ns
    z_max = sim.zsize
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
    eps = 1.0e-6
    dots = corners @ direction
    p0 = corners[np.argmin(dots)] - direction * (eps * 10.0)
    up = np.array([0.0, 0.0, 1.0], dtype=float)
    if abs(np.dot(up, ns)) > 0.95:
        up = np.array([0.0, 1.0, 0.0], dtype=float)
    u1 = np.cross(up, ns)
    u1 /= np.linalg.norm(u1)
    u2 = np.cross(ns, u1)
    proj = (corners - p0) @ np.vstack([u1, u2]).T
    umin, vmin = proj.min(axis=0)
    umax, vmax = proj.max(axis=0)
    du = umax - umin
    dv = vmax - vmin
    umin -= halo_pad * du
    umax += halo_pad * du
    vmin -= halo_pad * dv
    vmax += halo_pad * dv
    return {
        "proj_corners": proj,
        "umin": umin,
        "umax": umax,
        "vmin": vmin,
        "vmax": vmax,
        "p0": p0,
        "u1": u1,
        "u2": u2,
        "bounds_min": bounds_min,
        "bounds_max": bounds_max,
    }
start = time.perf_counter()
ray_factor = 6.0
#sdir, veg_absorb, bud = directshortwave_py(
#    sim,
#    nsun=nsun,
#    irradiance=irradiance,
#    ray_scale=ray_factor,
#    ray_jitter=1.0,
#    return_hit_count=True,
#    extend_bounds=extend_bounds,
#)
elapsed_py = time.perf_counter() - start
print(f"Direct shortwave runtime (python): {elapsed_py:.3f} s")

start = time.perf_counter()
sdir, veg_absorb, bud = directshortwave_nb(
    sim,
    nsun=nsun,
    irradiance=irradiance,
    ray_scale=ray_factor,
    ray_jitter=1.0,
    return_hit_count=True,
    extend_bounds=extend_bounds,
)
elapsed_nb = time.perf_counter() - start
print(f"Direct shortwave runtime (numba, first call): {elapsed_nb:.3f} s")

start = time.perf_counter()
sdir, veg_absorb, bud = directshortwave_nb(
    sim,
    nsun=nsun,
    irradiance=irradiance,
    ray_scale=ray_factor,
    ray_jitter=1.0,
    return_hit_count=True,
    extend_bounds=extend_bounds,
)

#% ----------------------------------------------------------------

#facsec_c = sim.facsec["c"]
#facids = facsec_c["facid"].astype(int)
#areas = facsec_c["area"]
#if hasattr(sim, "facs") and sim.facs is not None and "area" in sim.facs:
#    facet_areas = sim.facs["area"]
#else:
#    facet_areas = sim.geom.stl.area_faces
#area_accum = np.zeros_like(facet_areas, dtype=float)
#np.add.at(area_accum, facids, areas)
#area_ratio = area_accum / facet_areas
#print(
#    "Facet area check (sections / facet): "
#    f"min={np.min(area_ratio):.6f}, max={np.max(area_ratio):.6f}, "
#    f"mean={np.mean(area_ratio):.6f}"
#)
#mesh_areas = sim.geom.stl.area_faces
#if hasattr(sim, "facs") and sim.facs is not None and "area" in sim.facs:
#    facs_areas = sim.facs["area"]
#    mesh_ratio = mesh_areas / facs_areas
#    print(
#        "Facet area check (mesh / facs): "
#        f"min={np.min(mesh_ratio):.6f}, max={np.max(mesh_ratio):.6f}, "
#        f"mean={np.mean(mesh_ratio):.6f}"
#    )

elapsed = time.perf_counter() - start
print(f"Sdir facets: {sdir.shape}  veg_absorb: {veg_absorb.shape}")
print(f"Direct shortwave runtime (numba, warm): {elapsed:.3f} s")
if "rays" in bud:
    print(f"Rays cast: {bud['rays']}")
ray_debug = _ray_debug_info(sim, nsun)

#% ----------------------------------------------------------------

print("Conservation check (kW):")
kw = 1.0e-3
incoming_kw = bud["in"] * kw
exit_kw = bud["out"] * kw
absorb_kw = bud["veg"] * kw
facets_kw = bud["fac"] * kw
balance_kW = incoming_kw - exit_kw - absorb_kw - facets_kw
print(f"  incoming   : {incoming_kw:12.2f} kW")
print(f"  exiting    : {exit_kw:12.2f} kW")
print(f"  vegetation : {absorb_kw:12.2f} kW")
print(f"  facets     : {facets_kw:12.2f} kW")
print(f"  balance    : {balance_kW:12.2f} kW")
print(f"  solid_hit  : {bud['sol'] * kw:12.2f} kW (should match facets)")
if "hit_count" in bud:
    hit_count = bud["hit_count"]

#% ----------------------------------------------------------------
start = time.perf_counter()

veg_path = base_path / f"veg.inp.{sim_id}"
veg_pts = np.loadtxt(veg_path, skiprows=1)
if veg_pts.ndim == 1:
    veg_pts = veg_pts.reshape(1, -1)
i = veg_pts[:, 0].astype(int) - 1
j = veg_pts[:, 1].astype(int) - 1
k = veg_pts[:, 2].astype(int) - 1
if i.size > 0:
    i_min = int(np.min(i))
    i_max = int(np.max(i))
    j_min = int(np.min(j))
    j_max = int(np.max(j))
    k_min = int(np.min(k))
    k_max = int(np.max(k))
    veg_bbox = {
        "min": np.array([i_min * sim.dx, j_min * sim.dy, sim.zm[k_min]]),
        "max": np.array(
            [
                (i_max + 1) * sim.dx,
                (j_max + 1) * sim.dy,
                sim.zm[k_max + 1] if (k_max + 1) < len(sim.zm) else sim.zsize,
            ]
        ),
    }
else:
    veg_bbox = None

j_mid = sim.jtot // 2
mask_mid = j == j_mid
veg_slice = np.zeros((sim.itot, sim.ktot), dtype=float)
for idx in np.where(mask_mid)[0]:
    veg_slice[i[idx], k[idx]] += veg_absorb[idx]

fig1, ax1 = plt.subplots(figsize=(9, 6))
xg, zg = np.meshgrid(sim.xt, sim.zt, indexing="ij")
cs = ax1.contourf(xg, zg, veg_slice, levels=30, cmap="viridis")
fig1.colorbar(cs, ax=ax1, label="Veg absorption (W/m3)")
ax1.set_xlabel("x")
ax1.set_ylabel("z")
ax1.set_title(f"Vegetation absorption at y-index {j_mid}")
ax1.set_aspect("equal", adjustable="box")
ray_dir = -np.asarray(nsun, dtype=float)
ray_dir /= np.linalg.norm(ray_dir)
arrow_len = 0.2 * min(sim.xlen, sim.zsize)
ax1.quiver(
    sim.xlen * 0.1,
    sim.zsize * 0.1,
    ray_dir[0],
    ray_dir[2],
    angles="xy",
    scale_units="xy",
    scale=1.0 / arrow_len,
    color="white",
    width=0.005,
)

#% ----------------------------------------------------------------

#% ----------------------------------------------------------------
start = time.perf_counter()

mesh = sim.geom.stl
faces = mesh.faces
vertices = mesh.vertices
tris = vertices[faces]

fig2 = plt.figure(figsize=(9, 6))
ax2 = fig2.add_subplot(111, projection="3d")
poly = Poly3DCollection(tris, linewidths=0.0, alpha=1.0)
poly.set_array(sdir)
poly.set_cmap("inferno")
ax2.add_collection3d(poly)
fig2.colorbar(poly, ax=ax2, label="Facet Sdir (W/m2)")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("z")
ax2.set_title("Facet direct radiation (top view)")
ax2.set_box_aspect([sim.xlen, sim.ylen, sim.zsize])
ax2.view_init(elev=90, azim=-90)
arrow_len_xy = 0.2 * min(sim.xlen, sim.ylen)
ax2.quiver(
    sim.xlen * 0.1,
    sim.ylen * 0.1,
    sim.zsize * 0.05,
    ray_dir[0],
    ray_dir[1],
    0.0,
    length=arrow_len_xy,
    normalize=True,
    color="red",
    linewidth=2.0,
)

veg_x = sim.xt[i]
veg_y = sim.yt[j]
veg_z = sim.zt[k]
ax2.scatter(
    veg_x,
    veg_y,
    np.full_like(veg_x, sim.zsize * 0.02),
    s=2,
    c="white",
    alpha=0.6,
)
if veg_bbox is not None:
    vmin = veg_bbox["min"]
    vmax = veg_bbox["max"]
    vcorners = np.array(
        [
            [vmin[0], vmin[1], vmin[2]],
            [vmax[0], vmin[1], vmin[2]],
            [vmin[0], vmax[1], vmin[2]],
            [vmax[0], vmax[1], vmin[2]],
            [vmin[0], vmin[1], vmax[2]],
            [vmax[0], vmin[1], vmax[2]],
            [vmin[0], vmax[1], vmax[2]],
            [vmax[0], vmax[1], vmax[2]],
        ],
        dtype=float,
    )
    edges = [
        (0, 1), (1, 3), (3, 2), (2, 0),
        (4, 5), (5, 7), (7, 6), (6, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ]
    for e0, e1 in edges:
        ax2.plot(
            [vcorners[e0, 0], vcorners[e1, 0]],
            [vcorners[e0, 1], vcorners[e1, 1]],
            [vcorners[e0, 2], vcorners[e1, 2]],
            color="white",
            linewidth=1.0,
        )

#% ----------------------------------------------------------------
if "hit_count" in bud:
    hit_xy = np.sum(hit_count, axis=2)
    fig3, ax3 = plt.subplots(figsize=(9, 6))
    xg_xy, yg_xy = np.meshgrid(sim.xt, sim.yt, indexing="ij")
    hm = ax3.pcolormesh(
        xg_xy,
        yg_xy,
        hit_xy,
        shading="auto",
        cmap="viridis",
    )
    fig3.colorbar(hm, ax=ax3, label="Ray hit count (sum over z)")
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_title("Ray hit count (plan view)")
    ax3.set_aspect("equal", adjustable="box")

#% ----------------------------------------------------------------
if ray_debug is not None:
    proj = ray_debug["proj_corners"]
    p0 = ray_debug["p0"]
    u1 = ray_debug["u1"]
    u2 = ray_debug["u2"]
    umin = ray_debug["umin"]
    umax = ray_debug["umax"]
    vmin = ray_debug["vmin"]
    vmax = ray_debug["vmax"]
    rect_u = [umin, umax, umax, umin, umin]
    rect_v = [vmin, vmin, vmax, vmax, vmin]
    fig4, ax4 = plt.subplots(figsize=(9, 6))
    edges = [
        (0, 1), (1, 3), (3, 2), (2, 0),
        (4, 5), (5, 7), (7, 6), (6, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ]
    for e0, e1 in edges:
        ax4.plot(
            [proj[e0, 0], proj[e1, 0]],
            [proj[e0, 1], proj[e1, 1]],
            "k-",
        )
    ax4.plot([], [], "k-", label="BBox wireframe")
    ax4.plot(rect_u, rect_v, "r--", label="Ray plane bounds")
    ax4.scatter(proj[:, 0], proj[:, 1], s=20, c="k")
    if veg_bbox is not None:
        vmin_xyz = veg_bbox["min"]
        vmax_xyz = veg_bbox["max"]
        vcorners = np.array(
            [
                [vmin_xyz[0], vmin_xyz[1], vmin_xyz[2]],
                [vmax_xyz[0], vmin_xyz[1], vmin_xyz[2]],
                [vmin_xyz[0], vmax_xyz[1], vmin_xyz[2]],
                [vmax_xyz[0], vmax_xyz[1], vmin_xyz[2]],
                [vmin_xyz[0], vmin_xyz[1], vmax_xyz[2]],
                [vmax_xyz[0], vmin_xyz[1], vmax_xyz[2]],
                [vmin_xyz[0], vmax_xyz[1], vmax_xyz[2]],
                [vmax_xyz[0], vmax_xyz[1], vmax_xyz[2]],
            ],
            dtype=float,
        )
        vproj = (vcorners - p0) @ np.vstack([u1, u2]).T
        edges = [
            (0, 1), (1, 3), (3, 2), (2, 0),
            (4, 5), (5, 7), (7, 6), (6, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ]
        for e0, e1 in edges:
            ax4.plot(
                [vproj[e0, 0], vproj[e1, 0]],
                [vproj[e0, 1], vproj[e1, 1]],
                color="cyan",
                linewidth=1.0,
            )
        ax4.plot([], [], color="cyan", label="Veg bbox")
    ax4.set_xlabel("u")
    ax4.set_ylabel("v")
    ax4.set_title("Bounding box projection onto ray plane")
    ax4.set_aspect("equal", adjustable="box")
    ax4.legend()

#% ----------------------------------------------------------------
if ray_debug is not None:
    fig5 = plt.figure(figsize=(9, 6))
    ax5 = fig5.add_subplot(111, projection="3d")
    poly_sun = Poly3DCollection(tris, linewidths=0.0, alpha=1.0)
    poly_sun.set_array(sdir)
    poly_sun.set_cmap("inferno")
    ax5.add_collection3d(poly_sun)
    fig5.colorbar(poly_sun, ax=ax5, label="Facet Sdir (W/m2)")
    ax5.set_xlabel("x")
    ax5.set_ylabel("y")
    ax5.set_zlabel("z")
    ax5.set_title("Facet radiation (view from sun)")
    ax5.set_box_aspect([sim.xlen, sim.ylen, sim.zsize])
    ray_dir_3d = -np.asarray(nsun, dtype=float)
    ray_dir_3d /= np.linalg.norm(ray_dir_3d)
    view_az = np.rad2deg(np.arctan2(ray_dir_3d[1], ray_dir_3d[0]))
    view_el = np.rad2deg(np.arcsin(ray_dir_3d[2]))
    ax5.view_init(elev=view_el, azim=view_az)
    bounds_min = ray_debug["bounds_min"]
    bounds_max = ray_debug["bounds_max"]
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
    edges = [
        (0, 1), (1, 3), (3, 2), (2, 0),
        (4, 5), (5, 7), (7, 6), (6, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ]
    for e0, e1 in edges:
        ax5.plot(
            [corners[e0, 0], corners[e1, 0]],
            [corners[e0, 1], corners[e1, 1]],
            [corners[e0, 2], corners[e1, 2]],
            color="black",
            linewidth=1.0,
        )
    if veg_bbox is not None:
        vmin = veg_bbox["min"]
        vmax = veg_bbox["max"]
        vcorners = np.array(
            [
                [vmin[0], vmin[1], vmin[2]],
                [vmax[0], vmin[1], vmin[2]],
                [vmin[0], vmax[1], vmin[2]],
                [vmax[0], vmax[1], vmin[2]],
                [vmin[0], vmin[1], vmax[2]],
                [vmax[0], vmin[1], vmax[2]],
                [vmin[0], vmax[1], vmax[2]],
                [vmax[0], vmax[1], vmax[2]],
            ],
            dtype=float,
        )
        for e0, e1 in edges:
            ax5.plot(
                [vcorners[e0, 0], vcorners[e1, 0]],
                [vcorners[e0, 1], vcorners[e1, 1]],
                [vcorners[e0, 2], vcorners[e1, 2]],
                color="cyan",
                linewidth=1.0,
            )

plt.tight_layout()
plt.show()

elapsed = time.perf_counter() - start
print(f"plotting runtime: {elapsed:.3f} s")
