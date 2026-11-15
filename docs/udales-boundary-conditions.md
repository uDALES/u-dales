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

- 1: prescribed flux. A facet whose type is listed in the optional `factypeflux.inp` gets the heat flux of that type [K m/s] (see the facet types in the pre-processing documentation). Otherwise a facet with normal exactly in the -/+ x-direction gets `bctfxm`/`bctfxp`, in the -/+ y-direction `bctfym`/`bctfyp`, and in the +z direction `bctfz`. Any other facet (sloped, or facing down) with no `factypeflux.inp` entry gets zero flux; a start-up line reports how many.
- 2: flux given by wall function.

Values other than 1 and 2 are rejected at start-up.

## Moisture

Determined by `iwallmoist`. Possible values:

- 1: prescribed flux, on green (`lGR`) facets only. A facet whose type is listed in `factypeflux.inp` gets the moisture flux of that type [kg/kg m/s]. Otherwise a facet with normal exactly in the -/+ x-direction gets `bcqfxm`/`bcqfxp`, in the -/+ y-direction `bcqfym`/`bcqfyp`, and in the +z direction `bcqfz`. Any other facet with no `factypeflux.inp` entry gets zero flux; a start-up line reports how many.
- 2: flux given by wall function. This needs the heat transfer coefficient of the temperature wall function, so it requires `ltempeq = .true.` and `iwalltemp = 2`; the run stops at start-up otherwise. Where the heat transfer coefficient times the tangential velocity is zero, the moisture flux is zero.

Heat and moisture fluxes are independent: with `iwalltemp = 2` and `iwallmoist = 1` the moisture flux comes only from `bcqf*` / `factypeflux.inp`.

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
