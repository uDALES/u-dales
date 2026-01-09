from __future__ import annotations
import time
from pathlib import Path
import sys

# Add the uDALES python path
udbase_path = Path(f"C:/Users/mvr/OneDrive - Imperial College London/codes/uDALES/u-dales").resolve()
tools_path = (udbase_path / "tools" / "python").resolve()
if tools_path not in sys.path:
    sys.path.insert(0, str(tools_path))
udprep_path = (udbase_path / "tools" / "python" / "udprep").resolve()
if str(udprep_path) not in sys.path:
    sys.path.append(str(udprep_path))

expnr = '064'
expdir = (udbase_path.parents[0] / "experiments" / expnr).resolve()

start = time.perf_counter()
from udbase import UDBase  # noqa: E402
from udprep.directshortwave_facsec import directshortwave as directshortwave_nb  # noqa: E402
import udprep.directshortwave_f2py as ds  # noqa: E402

import numpy as np

import plotly.graph_objects as go

elapsed = time.perf_counter() - start
print(f"loading libraries runtime: {elapsed:.3f} s")
# %% ----------------------------------------------------------------

start = time.perf_counter()

# Instantiate UDBase with geometry for direct shortwave
sim = UDBase(expnr, expdir)

elapsed = time.perf_counter() - start
print(f"UDbase startup runtime: {elapsed:.3f} s")

azimuth_deg = 20.0
elevation_deg = 15.0
az = np.deg2rad(azimuth_deg)
el = np.deg2rad(elevation_deg)
nsun = [np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)]
irradiance = 800.0
extend_bounds = True

ray_factor = 6.0
ray_jitter = 1.0

start = time.perf_counter()
sdir, veg_absorb, bud = directshortwave_nb(
    sim,
    nsun=nsun,
    irradiance=irradiance,
    ray_scale=ray_factor,
    ray_jitter=ray_jitter,
    ray_jitter_seed=0,
    return_hit_count=True,
    extend_bounds=extend_bounds,
    periodic_xy=False
)
elapsed_nb = time.perf_counter() - start
print(f"Direct shortwave runtime (numba): {elapsed_nb:.3f} s")

print(f"Sdir facets: {sdir.shape}  veg_absorb: {veg_absorb.shape}")
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

# f2py comparison (no helper functions)
sdir_f2py = None
elapsed_f2py = None

mesh = sim.geom.stl
vertices = np.asarray(mesh.vertices, dtype=float, order="F")
faces = np.asarray(mesh.faces, dtype=np.int32) + 1  # Fortran expects 1-based indices
incenter = np.asarray(mesh.triangles_center, dtype=float, order="F")
face_normal = np.asarray(mesh.face_normals, dtype=float, order="F")
nsun_f = np.asarray(nsun, dtype=float)
resolution = 0.25*min(sim.dx, sim.dy, min(sim.dzt)) / ray_factor
calc = ds.directshortwave_f2py_mod.calculate_direct_shortwave
faces_f = np.asfortranarray(faces, dtype=np.int32)
incenter_f = np.asfortranarray(incenter, dtype=float)
face_normal_f = np.asfortranarray(face_normal, dtype=float)
start_f2py = time.perf_counter()
sdir_f2py = calc(
    faces_f,
    incenter_f,
    face_normal_f,
    vertices,
    nsun_f,
    float(irradiance),
    float(resolution),
)
sdir_f2py = np.asarray(sdir_f2py, dtype=float)
elapsed_f2py = time.perf_counter() - start_f2py
print(f"Direct shortwave runtime (f2py): {elapsed_f2py:.3f} s")
print(f"Sdir facets (f2py): {sdir_f2py.shape}")

# %% ----------------------------------------------------------------

numba_title = f"Facet Sdir (numba, {elapsed_nb:.3f} s)"
fig_numba = sim.plot_fac(sdir, show=False)
vmin = float(np.nanmin(sdir))
vmax = float(np.nanmax(sdir))
fig_numba.add_trace(
    go.Scatter3d(
        x=[0.0],
        y=[0.0],
        z=[0.0],
        mode="markers",
        marker=dict(
            size=1,
            opacity=0.0,
            color=[vmin],
            colorscale="Viridis",
            cmin=vmin,
            cmax=vmax,
            showscale=True,
            colorbar=dict(title="Facet value"),
        ),
        showlegend=False,
        hoverinfo="skip",
    )
)
fig_numba.update_layout(title=numba_title)
fig_numba.show(renderer="browser")

f2py_title = f"Facet Sdir (f2py, {elapsed_f2py:.3f} s)"
fig_f2py = sim.plot_fac(sdir_f2py, show=False)
vmin_f2py = float(np.nanmin(sdir_f2py))
vmax_f2py = float(np.nanmax(sdir_f2py))
fig_f2py.add_trace(
    go.Scatter3d(
        x=[0.0],
        y=[0.0],
        z=[0.0],
        mode="markers",
        marker=dict(
            size=1,
            opacity=0.0,
            color=[vmin_f2py],
            colorscale="Viridis",
            cmin=vmin_f2py,
            cmax=vmax_f2py,
            showscale=True,
            colorbar=dict(title="Facet value"),
        ),
        showlegend=False,
        hoverinfo="skip",
    )
)
fig_f2py.update_layout(title=f2py_title)
fig_f2py.show(renderer="browser")

diff = sdir_f2py - sdir
fig_hist = go.Figure(data=go.Histogram(x=diff, nbinsx=60))
fig_hist.update_layout(
    title="Histogram of Sdir differences (f2py - numba)",
    xaxis_title="Difference (W/m2)",
    yaxis_title="Count",
)
fig_hist.show(renderer="browser")
