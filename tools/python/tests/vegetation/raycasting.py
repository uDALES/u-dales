from __future__ import annotations
import time

start = time.perf_counter()


import sys
from pathlib import Path

sys.path.insert(0, "/rds/general/user/mvr/home/udales/u-dales/tools/python/")

from udbase import UDBase  # noqa: E402
from udprep import convert_block_to_sparse, directshortwave  # noqa: E402

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

elapsed = time.perf_counter() - start
print(f"loading libraries runtime: {elapsed:.3f} s")

#% ----------------------------------------------------------------

start = time.perf_counter()

sim_id = 525
base_path = Path("/rds/general/user/mvr/home/udales/u-dales/tests/tests_tree_input")


# Instantiate UDBase with geometry for direct shortwave
sim = UDBase(expnr=sim_id, path=base_path, load_geometry=True)
if sim.trees is None:
    raise RuntimeError(f"trees.inp.{sim_id} not loaded from {base_path}")

elapsed = time.perf_counter() - start
print(f"UDbase startup runtime: {elapsed:.3f} s")

azimuth_deg = 20.0
elevation_deg = 10.0
az = np.deg2rad(azimuth_deg)
el = np.deg2rad(elevation_deg)
nsun = [np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)]
irradiance = 800.0
start = time.perf_counter()
sdir, veg_absorb, bud = directshortwave(
    sim,
    nsun=nsun,
    irradiance=irradiance,
)
facsec_c = sim.facsec["c"]
facids = facsec_c["facid"].astype(int)
areas = facsec_c["area"]
if hasattr(sim, "facs") and sim.facs is not None and "area" in sim.facs:
    facet_areas = sim.facs["area"]
else:
    facet_areas = sim.geom.stl.area_faces
area_accum = np.zeros_like(facet_areas, dtype=float)
np.add.at(area_accum, facids, areas)
area_ratio = area_accum / facet_areas
print(
    "Facet area check (sections / facet): "
    f"min={np.min(area_ratio):.6f}, max={np.max(area_ratio):.6f}, "
    f"mean={np.mean(area_ratio):.6f}"
)
mesh_areas = sim.geom.stl.area_faces
if hasattr(sim, "facs") and sim.facs is not None and "area" in sim.facs:
    facs_areas = sim.facs["area"]
    mesh_ratio = mesh_areas / facs_areas
    print(
        "Facet area check (mesh / facs): "
        f"min={np.min(mesh_ratio):.6f}, max={np.max(mesh_ratio):.6f}, "
        f"mean={np.mean(mesh_ratio):.6f}"
    )
elapsed = time.perf_counter() - start
print(f"Sdir facets: {sdir.shape}  veg_absorb: {veg_absorb.shape}")
print(f"Direct shortwave runtime: {elapsed:.3f} s")

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
print(f"  solid_hit (should match facets) : {bud['sol'] * kw:12.2f} kW")

#% ----------------------------------------------------------------
start = time.perf_counter()

veg_path = base_path / f"veg.inp.{sim_id}"
veg_pts = np.loadtxt(veg_path, skiprows=1)
if veg_pts.ndim == 1:
    veg_pts = veg_pts.reshape(1, -1)
i = veg_pts[:, 0].astype(int) - 1
j = veg_pts[:, 1].astype(int) - 1
k = veg_pts[:, 2].astype(int) - 1

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
ax2.scatter(
    veg_x,
    veg_y,
    np.full_like(veg_x, sim.zsize * 0.02),
    s=2,
    c="white",
    alpha=0.6,
)

plt.tight_layout()
plt.show()

elapsed = time.perf_counter() - start
print(f"plotting runtime: {elapsed:.3f} s")
