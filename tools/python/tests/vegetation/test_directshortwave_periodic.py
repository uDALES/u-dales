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
from udprep.directshortwave import directshortwave  # noqa: E402

import plotly.graph_objects as go


def test_periodic_near_not_darker():
    expnr = "064"
    expdir = (udbase_path.parents[0] / "experiments" / expnr).resolve()

    sim = UDBase(expnr, expdir)

    azimuth_deg = 20.0
    elevation_deg = 15.0
    az = np.deg2rad(azimuth_deg)
    el = np.deg2rad(elevation_deg)
    nsun = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], dtype=float)
    nsun_unit = nsun / np.linalg.norm(nsun)

    sdir, veg_absorb, bud = directshortwave(
        sim,
        nsun=nsun,
        irradiance=800.0,
        ray_scale=6.0,
        ray_jitter=0.0,
        ray_jitter_seed=0,
        return_hit_count=False,
        extend_bounds=True,
        periodic_xy=True,
    )

    sdir_np, veg_absorb_np, bud_np = directshortwave(
        sim,
        nsun=nsun,
        irradiance=800.0,
        ray_scale=6.0,
        ray_jitter=0.0,
        ray_jitter_seed=0,
        return_hit_count=False,
        extend_bounds=True,
        periodic_xy=False,
    )

    def _print_budget(label: str, bud_data: dict) -> None:
        kw = 1.0e-3
        incoming_kw = bud_data["in"] * kw
        exit_kw = bud_data["out"] * kw
        absorb_kw = bud_data["veg"] * kw
        facets_kw = bud_data["fac"] * kw
        balance_kw = incoming_kw - exit_kw - absorb_kw - facets_kw
        print(f"{label} conservation check (kW):")
        print(f"  incoming   : {incoming_kw:12.2f} kW")
        print(f"  exiting    : {exit_kw:12.2f} kW")
        print(f"  vegetation : {absorb_kw:12.2f} kW")
        print(f"  facets     : {facets_kw:12.2f} kW")
        print(f"  balance    : {balance_kw:12.2f} kW")
        print(f"  solid_hit  : {bud_data['sol'] * kw:12.2f} kW (should match facets)")

    _print_budget("periodic_xy=True", bud)
    _print_budget("periodic_xy=False", bud_np)

    common = (sdir > 0.0) & (sdir_np > 0.0)
    if np.any(common):
        diff = np.abs(sdir[common] - sdir_np[common])
        base = np.maximum(sdir_np[common], 1.0e-12)
        rel = diff / base
        print(f"unoccluded comparison: count={diff.size} mean_rel={float(np.mean(rel)):.4f} max_rel={float(np.max(rel)):.4f}")
    else:
        print("unoccluded comparison: no common positive faces to compare")

    def _plot_sdir(title: str, sdir_data: np.ndarray) -> None:
        fig = sim.plot_fac(sdir_data, show=False)
        vmin = float(np.nanmin(sdir_data))
        vmax = float(np.nanmax(sdir_data))
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
        fig.update_layout(title=title)
        fig.show(renderer="browser")

    _plot_sdir("Facet Sdir (periodic_xy=True)", sdir)
    _plot_sdir("Facet Sdir (periodic_xy=False)", sdir_np)


if __name__ == "__main__":
    test_periodic_near_not_darker()
    print("OK")
