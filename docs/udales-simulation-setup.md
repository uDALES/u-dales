# Boundary Conditions

The boundary conditions for uDALES are specified under the `&BC` header in `namoptions.inp`.

## Momentum
### Top
Determined by `BCtopm`. Possible values:
- 1: free-slip, i.e. zero flux.
- 2: no-slip, i.e. zero velocity.
- 3: determined by inflow conditions.

### Bottom (below floor facets)
These are de facto useless, as the bottom of the domain is covered by floor facets. Determined by `BCbotm`. Possible values:
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
- 2/3/4: inflow/outflow (write more on this)

### y
Only possible value is 1: periodic.

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
- 1: constant flux given by `bctfxm` and `bctfxp` for x-walls, `bctfym` and `bctfyp` for y-walls, and `bctfz` for z-walls. (write more on this - what are x/y/z-walls)
- 2: flux given by wall function involving temperature.

## Moisture


## Scalars
