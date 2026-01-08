from __future__ import annotations
import time
import os
from pathlib import Path

start = time.perf_counter()

import sys

# Add the uDALES python path
udbase_path = Path(f"C:/Users/mvr/OneDrive - Imperial College London/codes/uDALES/u-dales").resolve()
tools_path = (udbase_path / "tools" / "python").resolve()
if tools_path not in sys.path:
    sys.path.insert(0, str(tools_path))


expnr = '105'
expdir = (udbase_path.parents[0] / "experiments" / expnr).resolve

from udbase import UDBase  # noqa: E402
from udprep import convert_block_to_sparse
from udprep.directshortwave import directshortwave as directshortwave_nb  # noqa: E402

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

elapsed = time.perf_counter() - start
print(f"loading libraries runtime: {elapsed:.3f} s")

#% ----------------------------------------------------------------

start = time.perf_counter()

# Instantiate UDBase with geometry for direct shortwave
sim = UDBase(expnr, expdir)

elapsed = time.perf_counter() - start
print(f"UDbase startup runtime: {elapsed:.3f} s")

# convert_block_to_sparse(sim)

azimuth_deg = 20.0
elevation_deg = 30.0
az = np.deg2rad(azimuth_deg)
el = np.deg2rad(elevation_deg)
nsun = [np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)]
irradiance = 800.0
extend_bounds = True

ray_factor = 4.0

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

print(f"Sdir facets: {sdir.shape}  veg_absorb: {veg_absorb.shape}")
print(f"Direct shortwave runtime (numba, warm): {elapsed:.3f} s")
if "rays" in bud:
    print(f"Rays cast: {bud['rays']}")

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
ray_dir = -np.asarray(nsun, dtype=float)
ray_dir /= np.linalg.norm(ray_dir)
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
plt.show()
