# Boundary conditions

The boundary conditions for uDALES are specified under the `&BC` header in `namoptions.inp`.

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

TBC