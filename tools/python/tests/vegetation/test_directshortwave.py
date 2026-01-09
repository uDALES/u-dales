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
udprep_path = (udbase_path / "tools" / "python" / "udprep").resolve()

expnr = '064'
expdir = (udbase_path.parents[0] / "experiments" / expnr).resolve

from udbase import UDBase  # noqa: E402
from udprep import convert_block_to_sparse
from udprep.directshortwave_nur1 import directshortwave as directshortwave_nb  # noqa: E402

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.graph_objects as go
import plotly.io as pio

elapsed = time.perf_counter() - start
print(f"loading libraries runtime: {elapsed:.3f} s")

def run_directshortwave_f2py(sim, nsun, irradiance, resolution):
    """
    Run the Fortran directShortwave via f2py.

    Build once (from repo root) with:
      python -m numpy.f2py -c -m directshortwave_f2py tools/python/fortran/directShortwave_f2py.f90
    """
    try:
        import directshortwave_f2py as ds
    except ImportError:
        ds = _load_f2py_from_udprep()

    mesh = sim.geom.stl
    vertices = np.asarray(mesh.vertices, dtype=float, order="F")
    faces = np.asarray(mesh.faces, dtype=np.int32) + 1  # Fortran expects 1-based indices
    use_incenter = False
    use_face_normal_from_faces = False
    if use_incenter:
        incenter = triangle_incenter(vertices, faces - 1)
    else:
        incenter = np.asarray(mesh.triangles_center, dtype=float, order="F")
    if use_face_normal_from_faces:
        face_normal = triangle_normals(vertices, faces - 1)
    else:
        face_normal = np.asarray(mesh.face_normals, dtype=float, order="F")
    nsun = np.asarray(nsun, dtype=float)

    calc = getattr(ds, "calculate_direct_shortwave", None) or getattr(ds, "calculateDirectShortwave", None)
    if calc is None and hasattr(ds, "directshortwave_f2py_mod"):
        calc = getattr(ds.directshortwave_f2py_mod, "calculate_direct_shortwave", None)
    if calc is None:
        raise AttributeError(
            "f2py module does not expose calculate_direct_shortwave (check directshortwave_f2py_mod)."
        )

    doc = (calc.__doc__ or "").replace(" ", "").lower()
    if "connectivitylist(nfaces,3)" in doc:
        faces_f = np.asfortranarray(faces, dtype=np.int32)
    elif "connectivitylist(3,nfaces)" in doc:
        faces_f = np.asfortranarray(faces.T, dtype=np.int32)
    else:
        faces_f = np.asfortranarray(faces, dtype=np.int32)

    n_faces = faces.shape[0]
    n_vertices = vertices.shape[0]
    include_sizes = "nfaces" in doc and "nvertices" in doc

    kwargs = {
        "connectivitylist": faces_f,
        "incenter": np.asfortranarray(incenter, dtype=float),
        "facenormal": face_normal,
        "vertices": vertices,
        "nsun": nsun,
        "irradiance": float(irradiance),
        "resolution": float(resolution),
    }
    if include_sizes:
        kwargs["nfaces"] = int(n_faces)
        kwargs["nvertices"] = int(n_vertices)

    try:
        sdir = calc(**kwargs)
    except TypeError:
        args = [faces_f, np.asfortranarray(incenter, dtype=float), face_normal, vertices]
        if include_sizes:
            args.extend([n_faces, n_vertices])
        args.extend([nsun, float(irradiance), float(resolution)])
        sdir = calc(*args)

    return np.asarray(sdir, dtype=float)


def triangle_incenter(vertices, faces):
    tri = vertices[faces]
    a = np.linalg.norm(tri[:, 1] - tri[:, 2], axis=1)
    b = np.linalg.norm(tri[:, 0] - tri[:, 2], axis=1)
    c = np.linalg.norm(tri[:, 0] - tri[:, 1], axis=1)
    w = a + b + c
    w = np.where(w == 0.0, 1.0, w)
    return (
        a[:, None] * tri[:, 0]
        + b[:, None] * tri[:, 1]
        + c[:, None] * tri[:, 2]
    ) / w[:, None]


def triangle_normals(vertices, faces):
    tri = vertices[faces]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm = np.linalg.norm(n, axis=1)
    norm = np.where(norm == 0.0, 1.0, norm)
    return (n / norm[:, None]).astype(float, copy=False, order="F")


def _load_f2py_from_udprep():
    candidates = list(udprep_path.glob("directshortwave_f2py*.pyd"))
    if not candidates:
        raise ImportError(
            "directshortwave_f2py not found. Build it with f2py first."
        )
    module_path = candidates[0]
    import importlib.util
    spec = importlib.util.spec_from_file_location("directshortwave_f2py", module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
# %% ----------------------------------------------------------------

start = time.perf_counter()

# Instantiate UDBase with geometry for direct shortwave
sim = UDBase(expnr, expdir)

elapsed = time.perf_counter() - start
print(f"UDbase startup runtime: {elapsed:.3f} s")

# convert_block_to_sparse(sim)

azimuth_deg = 20.0
elevation_deg = 15.0
az = np.deg2rad(azimuth_deg)
el = np.deg2rad(elevation_deg)
nsun = [np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)]
nsun_unit = np.asarray(nsun, dtype=float)
nsun_unit /= np.linalg.norm(nsun_unit)
irradiance = 800.0
extend_bounds = True

ray_factor = 6.0
ray_jitter = 0.0

# print(ds.__file__)
# print([name for name in dir(ds) if "shortwave" in name.lower() or "mod" in name.lower()])

use_f2py = True
sdir_f2py = None
elapsed_f2py = None
if use_f2py:
    resolution = min(sim.dx, sim.dy) / ray_factor
    start_f2py = time.perf_counter()
    sdir_f2py = run_directshortwave_f2py(sim, nsun, irradiance, resolution)
    elapsed_f2py = time.perf_counter() - start_f2py
    print(f"Direct shortwave runtime (f2py): {elapsed_f2py:.3f} s")
    print(f"Sdir facets (f2py): {sdir_f2py.shape}")

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
    periodic_xy=True
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

save_path = Path(__file__).parent / "directshortwave_numba_output.npz"

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

if sdir_f2py is not None:
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
