Page Under Development...

# Simulation setup

This document contains some notes on setting up uDALES simulations. A complete list of all parameters and switches is documented under [Configuration](./udales-namoptions-overview.md).

## Table of contents

- [Boundary Conditions](#Boundary-Conditions)
- [Buildings and blocks](#Buildings-and-blocks)
- [Driver simulations](#Driver-simulations)

## Boundary Conditions

The boundary conditions for uDALES are specified under the `&BC` header in `namoptions.inp`.

### Momentum

#### Top

Determined by `BCtopm`. Possible values:

- 1: free-slip, i.e. zero flux.
- 2: no-slip, i.e. zero velocity.
- 3: variable vertical velocity (necessary with inflow-outflow lateral boundary conditions).

#### Bottom

NB: only relevant if the bottom of the domain is not covered by floor facets. Determined by `BCbotm`. Possible values:

- 1: free slip, i.e. zero flux.
- 2: flux given by wall function involving temperature.
- 3: flux given by neutral wall function.

#### Walls

Determined by `iwallmom`. Possible values:

- 1: free-slip, i.e. zero flux.
- 2: flux given by wall function involving temperature.
- 3: flux given by neutral wall function.

#### x

Determined by `BCxm`. Possible values:

- 1: periodic
- 2: inflow-outflow, fixed profile
- 3: inflow-outflow, inflow given by time-varying profile from precursor simulation

#### y

Determined by `BCym`. Possible values:

- 1: periodic
- 2: inflow-outflow, fixed profile

### Temperature

Only required when `ltempeq = .true.`

#### Top

Determined by `BCtopT`. Possible values:

- 1: constant flux given by `wttop`.
- 2: constant temperature given by `thl_top`.

##### Bottom

Determined by `BCbotT`. Possible values:

- 1: constant flux given by `wtsurf`.
- 2: flux given by wall function. Temperature of ghost cells below floor facets given by `thls`.

#### Walls

Determined by `iwalltemp`. Possible values:

- 1: constant flux given by `bctfxm`/`bctfxp` for facets with normal in -/+ x-direction, `bctfym`/`bctfyp` for facets with normal in -/+ y-direction, and `bctfz` for facets with normal in +z direction.
- 2: flux given by wall function.

### Moisture

Determined by `iwallmoist`. Possible values:

- 1: constant flux given by `bcqfxm`/`bcqfxp` for facets with normal in -/+ x-direction, `bcqfym`/`bctfyp` for facets with normal in -/+ y-direction, and `bcqfz` for facets with normal in +z direction.
- 2: flux given by wall function.

### Scalars

TBC

## Facets

uDALES uses the immersed boundary method (IBM) to capture buildings in the fluid domain. This geometry is prescribed by an STL file, and this is used to generate the necessary input files for the IBM are generated using uDALES' pre-processing routines (see [pre-processing](./udales-pre-processing.md)). 

The `facets.inp.xxx` file indicates the the type (corresponding to a row of `factypes.inp.xxx` and therefore indicating the roughness length, heat capacity etc. of that surface) and its surface normal. Finally `Tfacinit.inp.xxx` indicates the initial temperature of each facet where applicable.

## Driver simulations

The options for running precursor and driven simulations are specified under the `&DRIVER` header in namoptions.xxx. This provides two model functionalities:

1) To run a precursor simulation where instantaneous y-z planes at a specified index in the x-direction (`iplane`) are written every `dtdriver` seconds to output files (named \*driver\*).
2) To run a driven simulation where the inlet y-z plane (at `i=ib-1` and `i=ib`) is determined by reading the \*driver\* files of a precursor simulation that has already been run.

### Running precursor simulations

Precursor simulations are indicated by first setting `idriver = 1` in namoptions.xxx. The set-up of the precursor simulation is then up to the users discretion. The other variables to specify are:

- `tdriverstart` - specifies the number of seconds after which the \*driver\* files will start being written. If the precursor simulation is not a warm start then it is recommended to allow a sufficient amount of time for the flow field to develop to the desired state before starting to write to the \*driver\* files.
- `dtdriver` specifies the timestep in seconds at which the y-z planes are written to the \*driver\* files. It is important to set this to be small (ideally `dtdriver = dt`) in order to reduce the requirement to interpolate these fields in the resulting driven simulation. However, depending on the case, this can be made larger with the advantage being a reduction in size of the \*driver\* files that are produced.
- `iplane` is the index in the x-direction that you want to save the instantaneous y-z planes. For many cases the expected value is `iplane = ie` so that the outlet of the precursor simulation is saved.
- `driverstore` is the number of timesteps that the user wants to write to the \*driver\* files. The total simulation time should therefore be equal or greater to `tdriverstart + (driverstore-1)*dtdriver` seconds to ensure that the writing process completes. The driven simulation based off this precursor will be limited to a maximum run time of `(driverstore-1)*dtdriver` seconds (unless multiple precursors are run using warm starts).

Outputs:

- A file for each prognostic variable (e.g. `u0`, `v0`, `w0`, `thl0` (if `ltempeq = .true.`) etc.) for each processor.The file names follow '$var"driver\_"$nproc"."$expnr' where $var indicates the variable (NOTE: h is potential temperature under this convention). These files will hold the corresponding instantaneous y-z planes over the specified time period and are to be used to drive a driven simulation.

### Running driven simulations

It is necessary to first have run a simulation following the above instructions. Driven simulations are initiated by setting `idriver = 2`. The following are guidelines for setting up the driven simulation:

- The driven simulation must have the same `jtot`, `ysize`, `kmax` and `zgrid.inp.xxx` as its corresponding precursor simulation.
- The driven simulation must use the same number of cores as the precursor simulation.
- It is not necessary to apply a forcing to the driven simulation due to the enforced inlet-outlet boundary conditions.
- `BCxm = 3` is the current index to enforce the required inlet-outlet boundary conditions. Boundary conditions in the x-direction for other prognostic variables will be overwritten by this and therefore do not need to be set.
- The \*driver\* files from the precursor simulation must be copied from its output directory to the experiments directory of the driven simulation. For example, if the precursor is 001 and the driven simulation is 002 and you are working from the top uDALES directory: `cp outputs/001/*driver* experiments/002/`.
- `driverstore` must be equal to or less than the number of timesteps saved in the \*driver\* files (equivalent value of `driverstore` in the precursor simulation).
- `driverjobnr` must equal the job number of the corresponding precursor simulation. Following the above example: `driverjobnr = 001`.
- `runtime` must equal `(driverstore-1)*dtdriver` seconds or less where `dtdriver` is from the precursor simulation. If this time limit is exceeded the simulation will be stopped as there will be no more data available to determine the inlet.
- NOTE: It is important that buildings are not positioned too close to the downwind edge of the domain. The vortices and wakes downwind in the nearfield of the buildings can cause errors alongside the convective outflow boundary condition. The necessary distance is case specific. If this distance is too small a typical error will be that high velocities are found at the domain edge and these lead to slow simulations and dt tending to nought.
