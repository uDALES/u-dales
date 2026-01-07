# Direct shortwave computation (tools/python/udprep/directshortwave.py)

This document describes how `directshortwave.py` computes direct shortwave
irradiance on facets and absorbed radiation in vegetation cells.

## Inputs and outputs

Inputs
- `sim`: UDBase object with grid and geometry.
- `nsun`: sun direction vector (any non-zero length).
- `irradiance`: direct normal irradiance (W/m2).
- `ray_scale`: scale factor for ray spacing (>= 1 means coarser rays, < 1 means finer rays).

Outputs
- `Sdir`: direct shortwave on facets, shape `(n_facets,)`, W/m2.
- `veg_absorb`: absorbed radiation per vegetation point, aligned with
  `veg.inp.<expnr>` order, W/m3.

## Vegetation data and sparse mapping

1) Vegetation points are loaded from `veg.inp.<expnr>` and stored as 1-based
   indices (i, j, k).
2) Vegetation parameters are loaded from `veg_params.inp.<expnr>`; only
   `lad` (leaf area density) and `dec` (extinction coefficient) are used.
3) The code builds three 3D arrays that map the sparse vegetation into the
   grid:
   - `lad_3d[itot, jtot, ktot]`
   - `dec_3d[itot, jtot, ktot]`
   - `veg_index[itot, jtot, ktot]` (index back into sparse veg list)

These arrays are truncated in the vertical to reduce memory (see next section).

## Vertical truncation (memory reduction)

The vertical size of 3D arrays is reduced to the highest occupied cell in
buildings or vegetation, plus one safety cell. The logic is:

- `kmax_solid` is taken from `sim.Sc` (solid mask on the cell-centered grid).
  If any solid cells are present, the maximum occupied k index is used.
- `kmax_veg` is taken from the sparse vegetation points (max of k-1).
- `kmax = max(kmax_solid, kmax_veg)`
- `ktot = min(sim.ktot, kmax + 2)`  (adds one safety cell)

All 3D arrays and vertical coordinates are sliced to this `ktot`.

## Grid and bounds setup

- `nsun` is normalized; the ray direction is `direction = -nsun_unit`.
- The vertical coordinate edges are `z_edges = [sim.zm, sim.zsize]`, then
  truncated to `ktot + 1`.
- The domain bounds for ray intersection use `z_max = z_edges[-1]` instead of
  `sim.zsize`.

## Ray generation

1) A ray origin plane is constructed just outside the domain by projecting the
   8 domain corners onto a plane orthogonal to `nsun`.
2) Two orthonormal basis vectors (`u1`, `u2`) span this plane.
3) The projected domain extents in this plane define a rectangle.
4) Rays are launched on a regular grid across this rectangle.
   - `step = min(dx, dy, min(dzt)) / ray_scale`
   - Each ray represents area `ray_area = step * step`

## Ray tracing through the grid

For each ray:

1) Intersect the ray with the domain AABB. If it does not enter, skip.
2) Start the ray just inside the entry point.
3) Traverse the grid with a 3D DDA stepping algorithm:
   - Track `t_max_*` and `t_delta_*` for x, y, z crossings.
   - Advance to the next cell boundary each step.
4) Stop if the ray leaves the domain or hits a solid cell.
5) For each cell crossed:
   - `trans_in[i,j,k]` is set on first visit to the current ray transmittance.
   - If `lad` and `dec` are positive, apply Beer-Lambert attenuation:
     - `tau = lad * dec * ds`
     - `r_out = r_in * exp(-tau)`
     - Absorbed energy is `(r_in - r_out) * ray_area`
   - Accumulate absorbed energy into `veg_absorb` (per vegetation point),
     converted to W/m3 by dividing by cell volume.

`trans_in` is initialized to -1 and stores the incoming transmittance for cells
hit by any ray. Unvisited cells remain negative.

## Facet direct irradiance (Sdir)

The direct irradiance for each facet is computed after the ray pass:

1) Compute `cos_inc = max(0, dot(face_normal, nsun_unit))` for all facets.
2) If facet sections (`sim.facsec['c']`) are available, each section cell
   contributes area-weighted irradiance using the cell transmittance.
3) Otherwise, each facet uses its center point to query the cell transmittance.
4) Final formula (per facet):

   `Sdir = irradiance * cos_inc * trans`

where `trans` is the incoming transmittance stored in `trans_in` for the
containing cell.

## Notes

- Rays terminate when they intersect solid cells; no attenuation is applied
  above the highest solid/vegetation cell due to truncation.
- The safety padding (+1 cell) reduces the risk of truncating rays that graze
  the highest objects.
