from __future__ import annotations

"""
User-facing walkthrough for direct shortwave preprocessing.

This tutorial mirrors the style of the shortwave test scripts, but stays focused
on what a user needs to run. It demonstrates:
- Loading a simulation with UDPrep/UDBase.
- Comparing configurations in a no-vegetation case.
- Demonstrating periodic boundary support.
- Running a vegetation case.
- Using calc_short_wave to run direct shortwave, view factors, and reflections.

Notes
-----
- All geometry/vegetation I/O happens inside UDBase.
- The direct shortwave solvers do not read files directly.
- The scanline method requires a compiled f2py module.
"""

from pathlib import Path
import sys
import time

import numpy as np

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
# calc_short_wave runs three pieces:
# - Direct shortwave S_dir on facets (defaults to the facsec method).
# - View factors + sky view factors using View3D.
# - Reflections to get net shortwave K_star (includes diffuse sky shortwave via Dsky).
#
# Direct shortwave methods are available via prep.radiation.calc_short_wave:
# 1) "scanline"  : f2py Fortran scanline rasterization (fastest, no vegetation)
# 2) "facsec"    : ray-casting with solid mask + facet-section reconstruction
#                  (default; supports vegetation, supports periodic_xy)
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
ray_density = 6.0  # reducing for faster demo runs; increase for production

plot_az_case1 = 30.0
plot_az_case2 = plot_az_case1
plot_az_case3 =  210.0

def print_field_stats(label: str, field: np.ndarray) -> None:
    """Print basic stats for a facet field."""
    print(
        f"{label} stats: "
        f"mean={float(np.mean(field)):.3f} "
        f"max={float(np.max(field)):.3f} "
        f"min={float(np.min(field)):.3f}"
    )

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
# Case 1: No vegetation, facsec baseline
# We report/plot both S_dir (direct) and K_star (net incl. reflections).
# -----------------------------------------------------------------------------
# Experiments typically live in the repo root (one level above u-dales).
expnr = "065"
expdir = udales_path.parents[0] / "experiments" / expnr

prep = UDPrep(expnr, expdir, load_geometry=True)
sim = prep.sim
prep.run_all()
exit()
print("Compute shortwave (facsec, no vegetation)...")
S_dir_facsec, K_star_facsec, S_veg_facsec = prep.radiation.calc_short_wave(
    nsun=nsun,
    irradiance=irradiance,
    ray_density=ray_density,
    periodic_xy=False,
    force=True,
)
print_field_stats("facsec S_dir", S_dir_facsec)
print_field_stats("facsec K_star", K_star_facsec)
plot_shortwave(sim, S_dir_facsec, "Direct shortwave (facsec, no vegetation)", plot_az_case1)
plot_shortwave(
    sim,
    K_star_facsec,
    f"Net shortwave (facsec, reflections, Dsky={prep.radiation.Dsky:.1f} W/m^2)",
    plot_az_case1,
)

# -----------------------------------------------------------------------------
# Case 2: Periodic boundaries (facsec)
# -----------------------------------------------------------------------------
# periodic_xy=True wraps rays across x/y boundaries to mimic an infinite tiling of
# the domain, so shading includes upwind buildings beyond the edges.
print("Compute shortwave (facsec, periodic)...")
S_dir_facsec_per, K_star_facsec_per, S_veg_facsec_per = prep.radiation.calc_short_wave(
    nsun=nsun,
    irradiance=irradiance,
    ray_density=ray_density,
    periodic_xy=True,
    force=True,
)
print_field_stats("facsec periodic S_dir", S_dir_facsec_per)
print_field_stats("facsec periodic K_star", K_star_facsec_per)
plot_shortwave(sim, S_dir_facsec_per, "Direct shortwave (facsec, periodic)", plot_az_case2)
plot_shortwave(
    sim,
    K_star_facsec_per,
    f"Net shortwave (facsec, periodic, reflections, Dsky={prep.radiation.Dsky:.1f} W/m^2)",
    plot_az_case2,
)

print("Timing cached calc_short_wave calls (should be fast)...")
t0 = time.perf_counter()
for _ in range(10):
    prep.radiation.calc_short_wave(
        nsun=nsun,
        irradiance=irradiance,
        ray_density=ray_density,
        periodic_xy=True,
        force=False,
    )
dt = time.perf_counter() - t0
print(f"10 cached calls: {dt:.3f} s (avg {dt / 10.0:.3f} s per call)")

# -----------------------------------------------------------------------------
# Case 3: Vegetation case (trees enabled)
# -----------------------------------------------------------------------------
expnr_veg = "525"
expdir_veg = udales_path / "tests" / "tests_tree_input"
prep_veg = UDPrep(expnr_veg, expdir_veg, load_geometry=True)
sim_veg = prep_veg.sim

print("Compute shortwave (facsec, vegetation)...")
S_dir_facsec_veg, K_star_facsec_veg, S_veg_facsec_veg = prep_veg.radiation.calc_short_wave(
    nsun=nsun,
    irradiance=irradiance,
    ray_density=ray_density,
    periodic_xy=False,
    force=True,
)
print_field_stats("facsec vegetation S_dir", S_dir_facsec_veg)
print_field_stats("facsec vegetation K_star", K_star_facsec_veg)
veg_data = sim_veg.load_veg(zero_based=True, cache=True)
plot_shortwave(
    sim_veg,
    S_dir_facsec_veg,
    "Direct shortwave (facsec, vegetation)",
    plot_az_case3,
    veg_data=veg_data,
)
plot_shortwave(
    sim_veg,
    K_star_facsec_veg,
    f"Net shortwave (facsec, vegetation, reflections, Dsky={prep_veg.radiation.Dsky:.1f} W/m^2)",
    plot_az_case3,
    veg_data=veg_data,
)
