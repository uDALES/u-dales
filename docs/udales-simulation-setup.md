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

# Driver simulations

The options for running precursor driver simulations and driven simulations are specified under the `&DRIVER` header in namoptions. This provides two model capabilities:
- 1) To run a precursor simulation where y-z planes at a specified index in the x-direction (`iplane`) are written every `dtdriver` seconds to output files (named \*driver\*)
2) To run a driven simulation where the inlet y-z plane (at `i=ib-1` and `i=ib`) is determined by reading the \*driver\* files of a precursor simulation that has already been run.

## Running precursor simulations

Precursor simulations are indicated by first setting `idriver = 1` in namoptions.inp.xxx. The set-up of the precursor simulation is then completely to the users discretion. Other important variables to specify:
- `tdriverstart` specifies the number of seconds at which the \*driver\* files will start being written. If the precursor simulation is not a warm start then it is recommended to allow a sufficient amount of time to allow the flow field to develop to the desired state before starting to write to the \*driver\* files. 
- `dtdriver` specifies the timestep at which the y-z planes are written to the \*driver\* files. It is important to set this to be small (ideally `dtdriver = dt`) in order to reduce the requirement to interpolate these fields in the resulting driven simulation. However, depending on the case, this can be made larger - with the advantage being a reduction in size of the \*driver\* files that are produced.
- `iplane` is the index in the x-direction that you want to save the instantaneous y-z planes. For many cases the expected value would be `iplane = ie` so that the outlet of the precursor simulation is saved.
- `driverstore` is the number of timesteps that the user wants to write to the \*driver\* files. The total simulation time should therefore be equal or greater to `tdriverstart + dtdriver\*driverstore` seconds to ensure that the writing process completes. The driven simulation based off this precursor will be limited to running `(driverstore-1)*dtdriver` seconds.

Outputs:
- A file for each prognostic variable (e.g. `u0`, `v0`, `w0`, `thl0` (if `ltempeq = .true.`) etc.) for each processor.The file names follow $var"driver\_"$nproc"."$expnr where $var indicates the variable (NOTE: h - potential temperature under this convention). These files will hold the corresponding instantaneous y-z planes over the specified time period and are to be used to drive a driven simulation.

## Running deriven simulations

It is necessary to first have run a simulation following the above instructions. Driven simulations are initiated by setting `idriver = 2`. The following are guidelines for setting up the driven simulation:
- The driven simulation must have the same `jtot`, `ysize`, `kmax` and `zgrid.inp.xxx` as its corresponding precursor simulation.
- It is not necessary to apply a forcing to the driven simulation due to the enforced inlet-outlet boundary conditions.
- `BCxm = 5` is the current index to enforce a driven simulation. **we could automate this at a later date if suitable**
- The \*driver\* files from the precursor simulation must be copied from its output directory to the experiments directory of the driven simulation. For example, if the precursor is 001 and the driven simulation is 002 and you are working from the top uDALES directory: `cp outputs/001/*driver* experiments/002/`.
- 'driverstore' must equal the number of timesteps saved in the \*driver\* files (equal to `driverstore` in the precursor simulation).
- `driverjobnr` must equal the job number of the corresponding precursor simulation. Following the above example: `driverjobnr = 001`.
