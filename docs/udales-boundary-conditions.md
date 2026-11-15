# Boundary conditions

The boundary conditions for uDALES are specified under the `&BC` header in `namoptions.###`.

## Momentum

### Top

Determined by `BCtopm`. Possible values:

- 1: free-slip, i.e. zero flux.
- 2: no-slip, i.e. zero velocity.
- 3: variable vertical velocity (necessary with inflow-outflow lateral boundary conditions).

### Bottom

NB: only relevant if the bottom of the domain is not covered by floor facets. Determined by `BCbotm`. Possible values:

- 1: free slip, i.e. zero flux.
- 2: flux given by wall function involving temperature.
- 3: flux given by neutral wall function.

### Walls

Determined by `iwallmom`. Possible values:

- 1: free-slip, i.e. zero flux.
- 2: flux given by wall function involving temperature.
- 3: flux given by neutral wall function.
- 4: Monin-Obukhov surface layer on the facets whose normal is within `mo_maxslope` (default 30 degrees) of the vertical, so on floors and gentle slopes; the other facets (walls, steep slopes) use the neutral wall function (3).

`iwallmom = 4` is the DALES dev 075ac4e surface layer for prescribed fluxes (`isurf = 4`, or `isurf = 3` with a prescribed friction velocity `ustin > 0` in `&WALLS`), with the Obukhov length computed locally (DALES `lmostlocal = .true.`), per facet section. In each section it takes the air state at the section's reference point (as the other wall functions), the facet roughness lengths `z0`, `z0h` (factypes columns 3 and 4) and the prescribed facet fluxes: `bctfz` for heat, and `bcqfz` for moisture on facets whose type has `lGR = 1` (factypes column 2; in this context the flag only means that the facet exchanges moisture). The drag coefficients use the previous Obukhov length, the new length follows from the bulk Richardson number with the surface `thl` and `qt` diagnosed from the fluxes, `u* = sqrt(Cm) max(|U|, 0.1 m/s)` (at least 0.01 m/s), and the stress is `u*^2` along the tangential wind. A Monin-Obukhov facet takes `bctfz`/`bcqfz` whatever its slope.

Requirements: `libm`, `ltempeq`, `iwalltemp = 1` (prescribed heat flux) and, with `lmoist`, `iwallmoist = 1`; `z0, z0h > 0` on every Monin-Obukhov facet. Because `iwalltemp = 1` prescribes the flux, the facet temperature does not follow it, so the steep facets use the neutral wall function, as `iwallmom = 2` already does under `iwalltemp = 1`. A start-up line reports how many facets use each. `bldump` (`lbldump`) writes the area means over the Monin-Obukhov sections of `u*`, the Obukhov length (as `1/mean(1/L)`), and the diagnosed surface `thl`, `qt`.

First level above horizontal Monin-Obukhov facets: in the fluid cells directly above a Monin-Obukhov facet whose normal is +z (to within 1e-6), `iwallmom = 4` also applies the DALES first-level treatment of the one-equation SGS model (`loneeqn`; DALES `sgs_surface_fix = .false.`). There the vertical shears in the TKE shear production are the surface-layer gradients `du/dz = u* phim(z1/L) / (kappa z1) u/|U|` (and the same for `v`), and `dthvdz` is the surface-layer gradient built from `dthl/dz = -w'thl' phih(z1/L) / (u* kappa z1)` and `dqt/dz` (fluxes upward-positive), which also enters the buoyancy source and the stable length scale. `dthvdz` in these cells is the surface-layer gradient for every SGS model, so it also enters the Vreman `lbuoycorr` correction and the `dt` field dump (on master it is 0 at `kb`). The Monin-Obukhov state of these cells is computed per c-grid facet section from the cell-centre `thl`, `qt` and horizontal wind (as DALES does per column) with the section distance as `z1`, before the subgrid step of each substep (the DALES order) and, at a cold start, once before the start-up thermodynamics. It is a separate state from that of the u, v, w sections that set the stress; both are identical in a horizontally homogeneous state. A start-up line counts the cells. Cells next to sloped Monin-Obukhov facets and to every other wall keep the general treatment: no TKE source at `kb`, `dthvdz(kb) = 0`, and the interior form above `kb`.

Caveats: the per-section state (Obukhov length, surface `thl`, `qt`) is not in the restart files; after a warm start it is cold-started again (surface values from the first level, `L = -0.1`) and re-converges within a few substeps. As with the other wall functions, a section is skipped while its wind is exactly zero. A flat floor of facets for a case without buildings is written by `tools/python/udprep/flat_floor.py`.

### x

Determined by `BCxm`. Possible values:

- 1: periodic
- 2: inflow-outflow, fixed profile
- 3: inflow-outflow, inflow given by time-varying profile from precursor simulation

### y

Determined by `BCym`. Possible values:

- 1: periodic
- 2: inflow-outflow, fixed profile

## Temperature

Global rule: all prescribed boundary fluxes (`wtsurf`, `wqsurf`, `wttop`, `wqtop`, facet `bctf*`/`bcqf*`) are positive in -z, i.e. directed downward (into the surface at the bottom), so a surface heating the air has a negative `wtsurf`.

Only required when `ltempeq = .true.`

### Top

Determined by `BCtopT`. Possible values:

- 1: constant flux given by `wttop`.
- 2: constant temperature given by `thl_top`.

### Bottom

Determined by `BCbotT`. Possible values:

- 1: constant flux given by `wtsurf`.
- 2: flux given by wall function. Temperature of ghost cells below floor facets given by `thls`.

### Walls

Determined by `iwalltemp`. Possible values:

- 1: constant flux given by `bctfxm`/`bctfxp` for facets with normal in -/+ x-direction, `bctfym`/`bctfyp` for facets with normal in -/+ y-direction, and `bctfz` for facets with normal in +z direction.
- 2: flux given by wall function.

## Moisture

Determined by `iwallmoist`. Possible values:

- 1: constant flux given by `bcqfxm`/`bcqfxp` for facets with normal in -/+ x-direction, `bcqfym`/`bctfyp` for facets with normal in -/+ y-direction, and `bcqfz` for facets with normal in +z direction.
- 2: flux given by wall function.

## Scalars

Only relevant when `nsv > 0`.

### Top

Determined by `BCtops`. Possible values:

- 1: constant flux given by `wsvtop`.
- 2: constant value given by `sv_top`.

### Bottom

NB: only relevant if the bottom of the domain is not covered by floor facets. Determined by `BCbots`. Possible values:

- 1: zero flux.

### x

Determined by `BCxs`. Possible values:

- 1: periodic
- 2: inflow-outflow, fixed profile
- 3: inflow-outflow, inflow given by time-varying profile from precursor simulation
- 4: inflow-outflow, fixed profile imposed only in a narrow strip around the domain's central y-plane; used for demonstration/testing purposes rather than production use.

### y

Determined by `BCys`. Possible values:

- 1: periodic
- 2: inflow-outflow, fixed profile
