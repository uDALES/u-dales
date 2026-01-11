from __future__ import annotations

import numpy as np
from pathlib import Path
import sys

# Add the uDALES python path
udbase_path = Path("C:/Users/mvr/OneDrive - Imperial College London/codes/uDALES/u-dales").resolve()
tools_path = (udbase_path / "tools" / "python").resolve()
if str(tools_path) not in sys.path:
    sys.path.insert(0, str(tools_path))

from udbase import UDBase  # noqa: E402
from udprep.directshortwave import DirectShortwaveSolver  # noqa: E402

import plotly.graph_objects as go

expnr = "064"
expdir = (udbase_path.parents[0] / "experiments" / expnr).resolve()

sim = UDBase(expnr, expdir)
veg_data = sim.load_veg(zero_based=True, cache=True) if getattr(sim, "ltrees", False) else None

azimuth_deg = 20.0
elevation_deg = 15.0
az = np.deg2rad(azimuth_deg)
el = np.deg2rad(elevation_deg)
nsun = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], dtype=float)
nsun_unit = nsun / np.linalg.norm(nsun)

results = {}
for periodic_xy in [False, True]:
    solver = DirectShortwaveSolver(
        sim,
        method="moller",
        ray_density=6.0,
        ray_jitter=0.0,
        veg_data=veg_data,
    )
    sdir, veg_absorb, bud = solver.compute(nsun=nsun, irradiance=800.0, periodic_xy=periodic_xy)
    results[periodic_xy] = (sdir, veg_absorb, bud)

    kw = 1.0e-3
    incoming_kw = bud["in"] * kw
    exit_kw = bud["out"] * kw
    absorb_kw = bud["veg"] * kw
    facets_kw = bud["fac"] * kw
    balance_kw = incoming_kw - exit_kw - absorb_kw - facets_kw
    label = f"periodic_xy={periodic_xy}"
    print(f"{label} conservation check (kW):")
    print(f"  incoming   : {incoming_kw:12.2f} kW")
    print(f"  exiting    : {exit_kw:12.2f} kW")
    print(f"  vegetation : {absorb_kw:12.2f} kW")
    print(f"  facets     : {facets_kw:12.2f} kW")
    print(f"  balance    : {balance_kw:12.2f} kW")
    print(f"  solid_hit  : {bud['sol'] * kw:12.2f} kW (should match facets)")

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
                colorbar=dict(title="Facet value"),
            ),
            showlegend=False,
            hoverinfo="skip",
        )
    )
    fig.update_layout(title=f"Facet Sdir (periodic_xy={periodic_xy})")
    fig.show(renderer="browser")

sdir = results[True][0]
sdir_np = results[False][0]
common = (sdir > 0.0) & (sdir_np > 0.0)
if np.any(common):
    diff = np.abs(sdir[common] - sdir_np[common])
    base = np.maximum(sdir_np[common], 1.0e-12)
    rel = diff / base
    print(f"unoccluded comparison: count={diff.size} mean_rel={float(np.mean(rel)):.4f} max_rel={float(np.max(rel)):.4f}")
else:
    print("unoccluded comparison: no common positive faces to compare")

print("OK")
