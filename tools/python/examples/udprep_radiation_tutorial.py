from __future__ import annotations

"""
User-facing walkthrough for direct shortwave preprocessing.

This tutorial mirrors the style of the shortwave test scripts, but stays focused
on what a user needs to run. It demonstrates:
- Loading a simulation with UDPrep/UDBase.
- Comparing methods in a no-vegetation case.
- Demonstrating periodic boundary support for facsec + moller.
- Running a vegetation case.

Notes
-----
- All geometry/vegetation I/O happens inside UDBase.
- The direct shortwave solvers do not read files directly.
- The scanline method requires a compiled f2py module.
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
import warnings

# -----------------------------------------------------------------------------
# Resolve uDALES python tools
# -----------------------------------------------------------------------------
udales_path = Path(__file__).resolve().parents[3]
PYTOOLS = udales_path / "tools" / "python"
if str(PYTOOLS) not in sys.path:
    sys.path.insert(0, str(PYTOOLS))

from udprep import UDPrep  # noqa: E402
import plotly.graph_objects as go  # noqa: E402

# -----------------------------------------------------------------------------
# Method overview
# -----------------------------------------------------------------------------
# Three direct shortwave methods are available via prep.radiation.calc_direct_sw:
# 1) "scanline"  : f2py Fortran scanline rasterization (fastest, no vegetation)
# 2) "facsec"    : ray-casting with solid mask + facet-section reconstruction
#                  (fast, supports vegetation, supports periodic_xy)
# 3) "moller"    : ray-casting with Moller-Trumbore triangle hits
#                  (most accurate, supports vegetation, supports periodic_xy)


# -----------------------------------------------------------------------------
# Sun vector and irradiance
# -----------------------------------------------------------------------------
azimuth_deg = 15.0
elevation_deg = 10.0

az = np.deg2rad(azimuth_deg)
el = np.deg2rad(elevation_deg)
nsun = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], dtype=float)
irradiance = 800.0

# Shared solver options
# ray_density controls the ray spacing on the sampling plane:
#   spacing = min(dx, dy, dz) / ray_density
# So ray_density=1 uses one ray per smallest grid spacing. The default value is 4, 
# but this can be too coarse for some geometries/facet sizes, and results may be 
# noisy or inaccurate. Always test sensitivity by increasing ray_density (e.g., 2, 4, 8)
# until results stabilize.
ray_density = 3.0 # reducing for faster demo runs; increase for production

plot_az_case1 = 30.0
plot_az_case2 = plot_az_case1
plot_az_case3 =  210.0

def print_budget(label: str, bud: dict[str, float]) -> None:
    """Pretty-print an energy budget summary."""
    kw = 1.0e-3
    incoming_kw = bud.get("in", 0.0) * kw
    exit_kw = bud.get("out", 0.0) * kw
    absorb_kw = bud.get("veg", 0.0) * kw
    facets_kw = bud.get("fac", 0.0) * kw
    balance_kw = incoming_kw - exit_kw - absorb_kw - facets_kw
    print(f"{label} conservation check (kW):")
    print(f"  incoming   : {incoming_kw:12.2f} kW")
    print(f"  exiting    : {exit_kw:12.2f} kW")
    print(f"  vegetation : {absorb_kw:12.2f} kW")
    print(f"  facets     : {facets_kw:12.2f} kW")
    print(f"  balance    : {balance_kw:12.2f} kW")

def plot_shortwave(
    sim,
    sdir: np.ndarray,
    title: str,
    view_azimuth_deg: float,
    veg_data: dict[str, np.ndarray] | None = None,
) -> None:
    """Plot facet field with an explicit colorbar range and fixed view azimuth."""
    fig = sim.plot_fac(sdir, show=False)
    vmin = float(np.nanmin(sdir))
    vmax = float(np.nanmax(sdir))
    fig.add_trace(
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
                colorbar=dict(title="[W/m²]"),
            ),
            showlegend=False,
            hoverinfo="skip",
        )
    )
    if veg_data is not None:
        points = np.asarray(veg_data.get("points", []), dtype=int)
        if points.ndim == 1 and points.size:
            points = points.reshape(1, -1)
        if points.size:
            i = points[:, 0]
            j = points[:, 1]
            k = points[:, 2]
            fig.add_trace(
                go.Scatter3d(
                    x=sim.xt[i],
                    y=sim.yt[j],
                    z=sim.zt[k],
                    mode="markers",
                    marker=dict(size=2, color="rgb(34,139,34)", opacity=0.7),
                    name="vegetation",
                    showlegend=True,
                )
            )
    az = np.deg2rad(view_azimuth_deg)
    el = np.deg2rad(20.0)
    dist = 1.75
    fig.update_layout(
        title=title,
        scene_camera=dict(
            eye=dict(
                x=float(dist * np.cos(el) * np.cos(az)),
                y=float(dist * np.cos(el) * np.sin(az)),
                z=float(dist * np.sin(el)),
            )
        ),
    )
    fig.show(renderer="browser")

# -----------------------------------------------------------------------------
# Case 1: No vegetation, compare facsec vs scanline (f2py)
# -----------------------------------------------------------------------------
# Experiments typically live in the repo root (one level above u-dales).
expnr = "065"
expdir = udales_path.parents[0] / "experiments" / expnr

prep = UDPrep(expnr, expdir, load_geometry=True)
sim = prep.sim

print("Compute direct shortwave (facsec, no vegetation)...")
sdir_facsec, _, bud_facsec = prep.radiation.calc_direct_sw(
    nsun=nsun,
    irradiance=irradiance,
    method="facsec",
    ray_density=ray_density,
    return_hit_count=False,
    periodic_xy=False,
)
print_budget("facsec", bud_facsec)
plot_shortwave(sim, sdir_facsec, "Direct shortwave (facsec, no vegetation)", plot_az_case1)

print("Compute direct shortwave (scanline, no vegetation)...")
sdir_scanline, _, bud_scanline = prep.radiation.calc_direct_sw(
    nsun=nsun,
    irradiance=irradiance,
    method="scanline",
    ray_density=ray_density,
)
print_budget("scanline", bud_scanline)
plot_shortwave(sim, sdir_scanline, "Direct shortwave (scanline, no vegetation)", plot_az_case1)
diff = sdir_scanline - sdir_facsec
print(
    "Scanline - facsec difference stats: "
    f"mean={float(np.mean(diff)):.3f} "
    f"max={float(np.max(diff)):.3f} "
    f"min={float(np.min(diff)):.3f}"
)

# -----------------------------------------------------------------------------
# Case 2: Periodic boundaries (facsec)
# -----------------------------------------------------------------------------
# periodic_xy=True wraps rays across x/y boundaries to mimic an infinite tiling of
# the domain, so shading includes upwind buildings beyond the edges.
print("Compute direct shortwave (facsec, periodic)...")
sdir_facsec_per, veg_absorb_facsec_per, bud_facsec_per = prep.radiation.calc_direct_sw(
    nsun=nsun,
    irradiance=irradiance,
    method="facsec",
    ray_density=ray_density,
    return_hit_count=False,
    periodic_xy=True,
)
print_budget("facsec periodic", bud_facsec_per)
plot_shortwave(sim, sdir_facsec_per, "Direct shortwave (facsec, periodic)", plot_az_case2)

# -----------------------------------------------------------------------------
# Case 3: Vegetation case (trees enabled)
# -----------------------------------------------------------------------------
expnr_veg = "525"
expdir_veg = udales_path / "tests" / "tests_tree_input"
prep_veg = UDPrep(expnr_veg, expdir_veg, load_geometry=True)
sim_veg = prep_veg.sim

print("Compute direct shortwave (facsec, vegetation)...")
sdir_facsec_veg, veg_absorb_facsec_veg, bud_facsec_veg = prep_veg.radiation.calc_direct_sw(
    nsun=nsun,
    irradiance=irradiance,
    method="facsec",
    ray_density=ray_density,
    periodic_xy=False,
    return_hit_count=False,
)
print_budget("facsec vegetation", bud_facsec_veg)
veg_data = sim_veg.load_veg(zero_based=True, cache=True)
plot_shortwave(
    sim_veg,
    sdir_facsec_veg,
    "Direct shortwave (facsec, vegetation)",
    plot_az_case3,
    veg_data=veg_data,
)

# Vegetation absorption slice (y = mid-plane)
j_mid = sim_veg.jtot // 2
veg_slice = np.zeros((sim_veg.itot, sim_veg.ktot), dtype=float)
veg_points = np.asarray(veg_data.get("points", []), dtype=int)
if veg_points.ndim == 1 and veg_points.size:
    veg_points = veg_points.reshape(1, -1)
if veg_points.size:
    if veg_absorb_facsec_veg.size != len(veg_points):
        warnings.warn(
            "veg_absorb length does not match vegetation points; "
            "skipping slice plot.",
            RuntimeWarning,
        )
    else:
        veg_3d = np.zeros((sim_veg.itot, sim_veg.jtot, sim_veg.ktot), dtype=float)
        i = veg_points[:, 0]
        j = veg_points[:, 1]
        k = veg_points[:, 2]
        veg_3d[i, j, k] += veg_absorb_facsec_veg
        veg_slice = veg_3d[:, j_mid, :]

fig_slice, ax_slice = plt.subplots(figsize=(9, 6))
xg, zg = np.meshgrid(sim_veg.xt, sim_veg.zt, indexing="ij")
cs = ax_slice.contourf(xg, zg, veg_slice, levels=30, cmap="viridis")
fig_slice.colorbar(cs, ax=ax_slice, label="Veg absorption (W/m3)")
ax_slice.set_xlabel("x")
ax_slice.set_ylabel("z")
ax_slice.set_title(f"Vegetation absorption at y-index {j_mid}")
ax_slice.set_aspect("equal", adjustable="box")
plt.show()