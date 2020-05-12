# Example simulations

uDALES can simulate a very large variety of urban case studies.
The main variables are:

- domain (size, equidistant or stretched grid, morphology)
- forcing in x and y (pressure gradient, volume flow rate, outflow rate, free stream velocity, driven, coriolis)
- lateral momentum BCs (periodic, driver)
- lateral scalar BCs (periodic, inflow-outflow, driver)
- scalar bottom BC/IBM (zero flux, constant flux, iso, energy balance)
- scalar top BC (zero flux, constant flux, iso)
- passive scalar sources (point, line, network)
- outputs (fielddump, averaging times, statsdumps)

Other variables:

- initial conditions (initial profiles)
- large scale forcings (subsidence, volumetric sources)
- nudging
- chemistry
- trees
- purifiers
- energy balance specifics - wall types, green roofs, radiation etc.

The following simulations are examples in which we combined several different features from above. We will explain the setup and parameters for each of them. A comprehensive description of all parameters can be found in the [Namoptions overview](./udales-namoptions-overview.md) page.

Please note that setting up these parameters also requires running the pre-processing routines, as some of the setups need to be in additional input files - see [Pre-processing](./udales-pre-processing.md).

## Neutral simulations

### 001

#### No buildings, bottom surface roughness

This simulation does not contain any buildings. Therefore we need to set:

```fortran
&INPS
lflat        = .true.
&WALLS
nfcts        = 0
```

A rough surface at the bottom is simulated by setting the following parameters:

```fortran
&BC
wtsurf       = 0.
wqsurf       = 0.
thls         = 288.
z0           = 0.01
z0h          = 0.000067
```

where `z0` is the roughness length for momentum. Note that even though this is a neutral simulation, all of these parameters need to be specified to reasonable values.
The simulation uses periodic lateral boundary conditions by default.

#### Constant pressure gradient

A constant pressure gradient in x and initial wind speed of u = 2 m/s are set by

```fortran
&INPS
dpdx         = 0.0001
u0           = 2.
```

#### Outputs

The simulation output contains the 3D instantaneous fields of u, v, and w, and 1D instantaneous and time-averaged vertical profiles. We specify this by setting:

```fortran
&OUTPUT
lfielddump   = .true.
fieldvars    = 'u0,v0,w0'
lxydump      = .true.
lxytdump     = .true.
```

### 002

This simulation has a similar setup to `001`.

#### Aligned cuboid buildings

Additionally, it contains aligned cube-shaped buildings, which can be automatically added by the pre-processing by setting:

```fortran
&INPS
lcube        = .true.
blockheight  = 16
blockwidth   = 16
canyonwidth  = 16
```

The corresponding number of blocks (including blocks on the floor) and block facets are:

```fortran
&WALLS
nblocks      = 17
nfcts        = 33
```

## Non-neutral simulations

### 101

#### Infinite canyon buildings

This simulation has "infinite canyons" along the y-axis as buildings. We set this by using

```fortran
&INPS
lcanyons     = .true.
```

#### Outputs (2)

Because there is no change in building geometry along y, it makes sense to look at the y-averaged 2D statistics:

```fortran
&OUTPUT
lydump       = .true.
lytdump      = .true.
```

#### Volume flow rate forcing

The simulation is driven by a fixed volume-flow rate forcing, which prescribes the domain-average velocity of u = 1.5 m/s:

```fortran
&PHYSICS
luvolflowr   = .true.
uflowrate    = 1.5
```

#### Temperature

The simulation also considers changes in temperature. We therefore need to set

```fortran
&PHYSICS
lbuoyancy    = .true.
ltempeq      = .true.
```

set an initial temperature profile

```fortran
&INPS
thl0         = 290.
```

and specify the advection scheme for temperature:

```fortran
&DYNAMICS
iadv_thl     = 2
```

#### Isothermal boundary conditions for temperature

A fixed temperature at the roof top (`thls`) and the top of the domain (`tthl_top`) with no-slip boundary conditions are set by

```fortran
&BC
thls         = 295.
thl_top      = 285.
BCtopT       = 2
```

The temperature at the building walls is determined by wall functions:

```fortran
&WALLS
iwalltemp    = 2
```

#### Passive scalars line source

A passive scalar line source is set by

```fortran
&SCALARS
nsv          = 1
lscasrcl     = .true.
SS           = 1.
sigS         = 0.5
```

and by specifying the scalar advection scheme:

```fortran
&DYNAMICS
iadv_sv      = 7
```

#### Outputs (3)

The temperature and scalar concentration are also outputs of the instantaneous fields:

```fortran
&OUTPUT
fieldvars    = 'u0,v0,w0,th,s1'
```

### 102

The basic setup of this simulation is similar to `101`.

#### Staggered cuboid buildings

The simulation has staggered cuboid buildings, which can be set up using

```fortran
&INPS
lstaggered   = .true.
```

#### Volume flow rate forcing (2)

The simulation is driven by a fixed volume-flow rate forcing for u and v:

```fortran
&PHYSICS
luvolflowr   = .true.
uflowrate    = 1.5
lvvolflowr   = .true.
vflowrate    = 0.3
```

#### Constant thermal flux boundary conditions for temperature

The temperature is determined by a constant thermal flux from the roads, building roofs and the top of the domain:

```fortran
&BC
wtsurf       = -0.01
wttop        = -0.01
bctfz        = -0.01
```

#### Passive scalars point source

The simulation contains a scalar field with two scalar point sources:

```fortran
&SCALARS
xS           = 4.
yS           = 8.
zS           = 3.
```

#### Scalar inflow-outflow boundary conditions

The scalar concentration is not determined by periodic boundary conditions like momentum and temperature, but leaves the domain at the outflow plane:

```fortran
&BC
BCxs         = 2
```

#### Warmstart

The simulation is continued from a previous simulation with a similar setup. The restart files (`initd` and `inits` for scalars) containing all relevant field data is saved in the example directory.

```fortran
&RUN
lwarmstart   = .true.
startfile    = 'initd00003199_xxx.102'
```

For a warmstart containing scalar concentration fields you need to additionally set

```fortran
&SCALARS
lreadscal    = .true.
```

## Simulation using energy balance

### 201

#### Pre-defined buildings

The buildings are read in by a file containing the block geometries:

```fortran
&INPS
lblocksfile  = .true.
blocksfile   = buildings.201
```

#### Grid-stretching in z

The vertical grid resolution is stretched by setting up the following parameters:

```fortran
&INPS
zsize        = 160
lzstretch    = .true.
stretchconst = 0.01
lstretchexp  = .true.
hlin         = 40
dzlin        = 1
```

`zsize` is an initial target value for the final domain height and `hlin` determines how many non-stretched grid cells there are at the lower end of the domain. Make sure to always use non-stretched grid cells whereever buildings are present.

#### Outputs (4)

3D time averaged output fields:

```fortran
&OUTPUT
ltdump       = .true.
```

#### Energybalance

The simulations solves the surface energybalance.

```fortran
&ENERGYBALANCE
lEB          = .true.
lconstW      = .true.
dtEB         = 2.
bldT         = 301.0
wsoil        = 314.0
wgrmax       = 451.0
wwilt        = 172.0
wfc          = 314.0
skyLW        = 201.0
GRLAI        = 2.1
rsmin        = 200.0
```

```fortran
&WALLS
iwalltemp    = 2
```

#### Moisture

The energybalance contains a latent heat flux and therefore we need to solve the full wet thermodynamics including moisture:

```fortran
&PHYSICS
lmoist       = .true.
```

The advection scheme for moisture is specified by:

```fortran
&DYNAMICS
iadv_qt      = 2
```

and domain top- and bottom-boundary values are chosen:

```fortran
&BC
qt_top       = 0.0
qts          = 0.0
```

#### Coriolis forcing and nudging

```fortran
&PHYSICS
lcoriol      = .true.
lnudge       = .true.
tnudge       = 10800.
nnudge       = 64
```

## Driver simulation

### 501

This simulation is similar to `101` but does not contain any scalars. It is used as the precursor simulation for simulation `502`.

#### Driver/precursor simulation

To save the output of a simulation into files that can be read by another simulation as driver-inputs, we set:

```fortran
&DRIVER
idriver      = 1
tdriverstart = 5.
dtdriver     = 1.
driverstore  = 6
iplane       = 128
```

### 502

#### Buildings from LIDAR image

The buildings of this simulation are generated from a grey-scale image of the buildings and their heights.

```fortran
&INPS
llidar       = .true.
sourcename   = DAPPLE7.png
dxinp        = 1
dyinp        = 1
dzinp        = 1
centeri      = 400
centerj      = 400
maxh         = 25
pad          = 3
smallarea    = 150
```

#### Driven simulation

This simulation is forced by the data from the stored outlet plane of simulation `501`:

```fortran
&DRIVER
idriver      = 2
driverjobnr  = 501
driverstore  = 6
```

All boundary conditions (momentum, temperature) are therefore set to inflow-outflow:

```fortran
&BC
BCxm         = 5
```

## Simulation with trees

### 801

This simulation is can be found on the branch `tomgrylls/trees-driver-patch`.

## Simulation with air purifiers

### 901

This simulation is can be found on the branch `tomgrylls/trees-driver-patch`.
