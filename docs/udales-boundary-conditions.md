# Boundary conditions

The boundary conditions for uDALES are specified under the `&BC` header in `namoptions.###`.

## Momentum

### Top

Determined by `BCtopm`. Possible values:

- 1: free-slip, i.e. zero flux.
- 2: no-slip, i.e. zero velocity.
- 3: variable vertical velocity (necessary with inflow-outflow lateral boundary conditions).

### Bottom

The bottom of the domain is always covered by ground facets, treated as immersed-boundary walls. The wall function is selected by `iwallmom` (see [Walls](#walls) below), with per-facet roughness given in `factypes.inp`.

The legacy flat-surface bottom BC keys (`BCbotm`, `BCbotT`, `BCbotq`, `BCbots`, `wtsurf`, `wqsurf`, `thls`, `qts`, `z0`, `z0h`, `wsvsurfdum`, `lbottom`) have been removed. The bottom boundary condition is now always given by ground facets; old namoptions files that still set these keys fail to parse (`stop 1` at the namelist read).

### Walls

Determined by `iwallmom`. Possible values:

- 1: free-slip, i.e. zero flux.
- 2: flux given by wall function involving temperature.
- 3: flux given by neutral wall function.

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

Only required when `ltempeq = .true.`

### Top

Determined by `BCtopT`. Possible values:

- 1: constant flux given by `wttop`.
- 2: constant temperature given by `thl_top`.

### Bottom

The bottom of the domain is always covered by ground facets, treated as immersed-boundary walls. The wall function is selected by `iwalltemp` (see [Walls](#walls) below), with per-facet initial temperature given in `Tfacinit.inp`.

The legacy flat-surface bottom BC keys (`BCbotm`, `BCbotT`, `BCbotq`, `BCbots`, `wtsurf`, `wqsurf`, `thls`, `qts`, `z0`, `z0h`, `wsvsurfdum`, `lbottom`) have been removed. The bottom boundary condition is now always given by ground facets; old namoptions files that still set these keys fail to parse (`stop 1` at the namelist read).

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
