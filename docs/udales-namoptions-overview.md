# Input parameters

Below we specify the input parameters in the `namoptions` file of your experiment. 

<!--
This list 
This list refers to the original code-base [DALES](https://github.com/dalesteam/dales). The latest version of the namoptions overview of DALES is documented [here](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf).
-->
## Namelist RUN

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| iexpnr | 000 | Three digit integer. | Experiment case number | - |
| lwarmstart | .false. | .true. or .false. | If .true. simulation reads in restart file to initialize fields.  | - |
| startfile | '' | 'initd00001234_xxx_xxx.000'| Name of restart file, the last three digits should match iexpnr. | - |
| runtime | 300 | > 0 | Simulation time. | [s] |
| trestart | 10000. | 0 < trestart < runtime | Time at which restart files are written, trestart > runtime will prevent restart files being written.  | [s] |
| dtmax | 20 | > 0 | Maximum allowed numerical integration timestep. | [s] |
| ladaptive | .false. | .true. or .false. | Switch for adaptive time-stepping, .true. recommended. | - |
| courant | 1.1 | 1 <= courant <=2 | Courant number, default sets it to 1.5 or 1.1 (if Kappa or upwind scheme is used). | - |
| lrandomize | .true. | .true. or .false. | Switch that determines whether initial field is randomised.| - | 
| irandom | 43 | `INTEGER` > 0 | Seed for random number generation. | - |
| randu | 0.01 |`REAL` > 0 | Amplitude of velocity field randomisation. | [m/s] |
<!---
| randthl | 0.0 | `REAL` > 0 | Amplitude of temperature field randomisation. | [K] |
| randqt | 0.0 | `REAL` > 0 | Amplitude of moisture field randomisation. | [kg/kg] |
--->
| libm | .true. | .true. or .false. | Switch that determines whether the Immersed Boundary Method is turned on. | - |
| lles | .true. | .true. or .false. | Switch that determines whether the subgrid model is turned on or constant ekm and ekh are used (DNS). | - |
| nprocx | -  | `INTEGER` > 0 | Number of pencils in the x-direction (see 2decomp documentation [https://2decomp-fft.github.io/]), must be a divisor of itot. | - |
| nprocy | -  | `INTEGER` > 0 | Number of pencils in the y-direction (see 2decomp documentation [https://2decomp-fft.github.io/]), must be a divisor of jtot and ktot. | - |



## Namelist DOMAIN

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| itot | 96 | `INTEGER` > 0  | Number of points in x-direction. | - |
| jtot | 96 | `INTEGER` > 0  | Number of points in y-direction. | - |
| ktot | 96 | `INTEGER` > 0  | Number of points in z-direction. | - |
| xlen | -1 | `REAL` > 0 | Domain size in x-direction.| - |
| ylen | -1 | `REAL` > 0 | Domain size in x-direction.| - |


## Namelist PHYSICS

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| ps | -1 | `REAL` > 0 | Air pressure at surface, recommend using standard pressure.| [Pa] |
| igrw_damp | 0 | 0, 1. | Switch to enable gravity wave damping.| - |
| ltempeq | .false. | .true. or .false. | Switch for solving temperature equation. | - |
| lbuoyancy | .false. | .true. or .false. | Switch for buoyancy force in temperature equation. | - |
| lmoist | .false. | .true. or .false. | Switch for solving moisture equation. | - |
| lcoriol | .false. | .true. or .false. | Switch for inluding the Coriolis force. | - |
| luoutflowr | .false. | .true. or .false. | Switch that determines whether u-velocity is corrected to get a fixed outflow rate. *Only functional when x-direction is not parellelised.* | |
| lvoutflowr | .false. | .true. or .false. | Switch that determines whether u-velocity is corrected to get a fixed outflow rate. *Only functional when y-direction is not parellelised.* | |
| luvolflowr | .false. | .true. or .false. | Switch that determines whether u-velocity is corrected to get a fixed volume flow rate. | - |
| lvvolflowr | .false. | .true. or .false. | Switch that determines whether u-velocity is corrected to get a fixed volume flow rate. | - |
| luflowr | .false. | .true. or .false. | Switch that determines whether u-velocity is corrected to get a fixed flow velocity. | - |
| lvflowr | .false. | .true. or .false. | Switch that determines whether v-velocity is corrected to get a fixed flow velocity. | - |
| uflowrate | 1. | `REAL` | U-velocity flow rate for out- or volume-flow forcing. | [m/s] |
| vflowrate | 1. | `REAL` | V-velocity flow rate out- or volume-flow forcing. | [m/s] |
| lprofforc | .false. | .true. or .false. | Switch for nudging flow to a profile (forcing). | - |
| lnudge | .false. | .true. or .false. | Switch for nudging flow to profiles (in `prof.inp.xxx`). | - |
| nnudge | 0 | `INTEGER` | Number of points from bottom to nudge. | - |
| tnudge | 60. | `REAL` | Time scale for nudging | - |
| ltimedepsurf | .false. | .true. or .false. | Switch for time-dependent surface heat flux (`bctfz` etc - see BC section). | - |
| ntimedepsurf | 0 | `REAL` | Number of time-dependent surface heat fluxes in file `timedepsurf.inp.xxx`. | - |
| ltimedepnudge | .false. | .true. or .false. | Switch for time-dependent profiles. When `lnudge = .true.`, then this switch causes the nudging profile to vary in time. If using inflow-outflow boundary conditions with inflow given by profile (`BCxm/BCxT/BCxq = 2`), then the inflow profile varies in time. | - |
| ntimedepnudge | 0 | `REAL` | Number of time-dependent nudging profiles in file `timedepnudge.inp.xxx`. | - |
| ltimedepsw | .false. | .true. or .false. | Time-dependent shortwave radiation on facets | - |
| ntimedepsw | 0 | `REAL` | Number of time-dependent shortwave radiative fluxes in file `timedepsw.inp.xxx`. | - |
| ltimedeplw | .false. | .true. or .false. | Time-dependent longwave radiation on facets. | - |
| ntimedeplw | 0 | `REAL` | Number of time-dependent shortwave radiative fluxes in file `timedeplw.inp.xxx`. | - |

<!---
| ifixuinf | 0 | 1, 2 | Choice for free stream forcing. (0 = nothing) | |
| lvinf | .false. | .true., .false. | Use Vinf instead of Uinf for the fixed velocity at infinity | |
| tscale | | | Timescale: domain height*Uinf/utau\*\*2 | |
| dpdx | 0. | | Constant pressure gradient forcing in x. | |
--->


## Namelist DYNAMICS

Possible advection schemes:

1 = 1st order upwind scheme

2 = 2nd order central difference scheme

7 = Kappa (flux limited) scheme. This scheme can only be applied to passive scalars.

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lqlnr | .false. | .true. or .false. | Logical for calculation of liquid water concentration. | - |
| ipoiss | 0 | 0 | Poisson solver. 0 = Fast Fourier Transform. | - |
| iadv_mom | 2 | 2 | Advection scheme for momentum. | - |
| iadv_tke | 2 | 2 | Advection scheme for TKE. Only used if `loneeqn = True`. | - |
| iadv_thl | 2 | 2 | Advection scheme for temperature. | - |
| iadv_qt | 2 | 2 | Advection scheme for moisture. | - |
| iadv_sv | 7 | 1, 2, 7 | Advection scheme for scalars. | - |

## Namelist BC

Switches for boundary conditions: momentum (m), temperature (T), humidity (q) and scalars (s).

Lateral BCs (BCx, BCy): 1 = periodic, 2 in/outflow conditions, inflow given by profile (usually constant), 3: in/outflow conditions, inflow given by precursor simulation.

BCs at the top (BCtop): 1 = freeslip, 2 = noslip, 3 = should be used with inflow/outflow conditions.

BCs at the bottom (BCbot; only effective if not covered with ground facets): 1 = flux, 2 = wall function, 3 = neutral wall function.

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| BCxm | 1 | 1,2,3 | Domain boundary condition for momentum in x. | - |
| BCxT | 1 | 1,2,3 | Domain boundary condition for temperature in x. | - |
| BCxq | 1 | 1,2,3 | Domain boundary condition for humidity in x. | - |
| BCxs | 1 | 1,2,3 | Domain boundary condition for scalars in x. | - |
| BCym | 1 | 1,2 | Domain boundary condition for momentum in y. | - |
| BCyT | 1 | 1,2 | Domain boundary condition for temperature in y. | - |
| BCyq | 1 | 1,2 | Domain boundary condition for humidity in y. | - |
| BCys | 1 | 1,2 | Domain boundary condition for scalars in y. | - |
| BCtopm | 1 | 1,2,3 | Boundary condition for momentum at domain top. | - |
| BCtopT | 1 | 1,2 | Boundary condition for temperature at domain top. | - |
| BCtopq | 1 | 1,2 | Boundary condition for humidity at domain top. | - |
| BCtops | 1 | 1,2 | Boundary condition for scalars at domain top. | - |
| bctfxm | 0 | `REAL` | Temperature flux on facets with surface normal in -x direction. |  [Km/s] |
| bctfxp | 0 | `REAL` | Temperature flux on facets with surface normal in +x direction. | [Km/s] |
| bctfym | 0 | `REAL` | Temperature flux on facets with surface normal in -y direction. | [Km/s] |
| bctfyp | 0 | `REAL` | Temperature flux on facets with surface normal in +y direction. | [Km/s] |
| bctfz | 0 | `REAL` | Temperature flux on facets with surface normal in +z direction. | [Km/s] |
| bcqfxm | 0 | `REAL` | Moisture flux on facets with surface normal in -x direction. | [m/s] |
| bcqfxp | 0 | `REAL` | Moisture flux on facets with surface normal in +x direction. | [m/s] |
| bcqfym | 0 | `REAL` | Moisture flux on facets with surface normal in -y direction. | [m/s] |
| bcqfyp | 0 | `REAL` | Moisture flux on facets with surface normal in +y direction. | [m/s] |
| bcqfz | 0 | `REAL` | Moisture flux on facets with surface normal in +z direction. | [m/s] |
| thl_top | -1. | `REAL` >= 0 | Temperature at the top boundary. | [K] |
| qt_top | -1. | `REAL` >= 0| Humidity at the top boundary. | [kg/kg] |
| wttop | 0. | `REAL` | Temperature flux at the top boundary. | [Km/s] |
| BCbotm | 2 | 1,2,3 | Boundary condition for momentum at domain bottom (if `lbottom = .true.`). | - |
| BCbotT | 1 | 1,2 | Boundary condition for temperature at domain bottom (if `lbottom = .true.`). | - |
| BCbotq | 1 | 1 | Boundary condition for humidity at domain bottom (if `lbottom = .true.`). | - |
| BCbots | 1 | 1 | Boundary condition for scalars at domain bottom (if `lbottom = .true.`). | - |
| wtsurf | -1. |`REAL` | Temperature flux at domain bottom (if `lbottom = .true.`). | [Km/s] |
| wqsurf | -1. | `REAL`| Moisture flux at domain bottom (if `lbottom = .true.`).  | [m/s] |
| thls | -1. |  `REAL`| Temperature at domain bottom (if `lbottom = .true.`).  | [K] |
| qts | -1. | `REAL` | Moisture at domain bottom (if `lbottom = .true.`). Used in modthermodynamics to get a BC for the moisture profile. | [kg/kg] |
| z0 | -1. |  `REAL`| Momentum roughness length of the domain bottom (if `lbottom = .true.`).  | [m] |
| z0h | -1. |  `REAL`| Heat roughness length of the domain bottom (if `lbottom = .true.`).| [m] |

<!---
| wsvtopdum | 0 | | Scalar boundary conditions top. | - |
| wsvsurfdum | | | Scalar flux at domain bottom (if `lbottom = .true.`). | - |
--->
## Namelist NAMSUBGRID

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lvreman | .false. | .true. or .false. | Switch for Vreman (2004) sub-grid scheme. | - |
<!-- | c_vreman | 0.07 | `REAL` | Model constant for Vreman scheme. | - | --->
<!---
| lbuoycorr | .false. | .true. or .false. | Switch for buoyancy correlation in the Vreman scheme. | - |
| loneeqn | .false. | .true. or .false. | Switch for one-equation sub-grid scheme. | - |
| ldelta | .false. | .true. or .false. | Switch for diminished sfs in stable flow. | - |
| lmason | .false. | .true. or .false.| Switch for decreased length scale near the surface | - |
--->
| cf | 2.5 | `REAL` > 0 | Filter constant. | - |
| cn | 0.76 | `REAL` > 0 | Subfilter scale parameter. | - |
| Rigc | 0.25 | `REAL` > 0 | Critical Richardson number. | - |
<!---
| Prandtl | 0.333 |`REAL` > 0 | Prandtl number. | - |
| lsmagorinsky | .false . | .true. or .false. | Switch for Smagorinsky subgrid scheme. | - |
| cs | -1 | > 0 | Smagorinsky constant. | - |
| nmason | 2 | > 0 | Exponent in Mason correction function.| - |
--->

<!---
## Namelist NAMCHECKSIM

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| tcheck | 0 | > 0 |  Time interval between checks of velocity divergence and Courant numbers.| [s] |
--->

## Namelist WALLS

1 = fixed flux

2 = flux determined by wall function involving temperature

3 = flux determined by neutral wall function (set automatically if `ltempeq = .false.`)

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| nblocks | 0 | `INTEGER` | Number of blocks specified in `blocks.inp`. | - |
| nfcts | -1 | `INTEGER` | Number of facets specified in `facets.inp`. | - |
| iwallmom | 2 | 1, 2, 3 (1 means zero flux) | Building wall momentum flux. | - |
| iwalltemp | 1 | 1, 2 |  Building wall temperature flux. | - |
| iwallmoist | 1 | 1, 2 |  Building wall moisture flux. | - |
| iwallscal | 1 | 1, 2 | Building wall scalar flux | - |
| lbottom | .false. | .true., .false. | Switch for using wall function as bottom BC. *Used only if no ground facets.* | - |
| nsolpts_u | 0 | `INTEGER` | Number of solid points on u-grid. | - |
| nsolpts_v | 0 | `INTEGER` | Number of solid points on v-grid. | - |
| nsolpts_w | 0 | `INTEGER` | Number of solid points on w-grid. | - |
| nsolpts_c | 0 | `INTEGER` | Number of solid points on c-grid. | - |
| nbndpts_u | 0 | `INTEGER` | Number of fluid boundary points on u-grid. | - |
| nbndpts_v | 0 | `INTEGER` | Number of fluid boundary points on v-grid. | - |
| nbndpts_w | 0 | `INTEGER` | Number of fluid boundary points on w-grid. | - |
| nbndpts_c | 0 | `INTEGER` | Number of fluid boundary points on c-grid. | - |
| nfctsecs_u | 0 | `INTEGER` | Number of facet sections on u-grid. | - |
| nfctsecs_v | 0 | `INTEGER` | Number of facet sections on v-grid. | - |
| nfctsecs_w | 0 | `INTEGER` | Number of facet sections on w-grid. | - |
| nfctsecs_c | 0 | `INTEGER` | Number of facet sections on c-grid. | - |
| lnorec | .false. | .true. or .false. | Switch for not using reconstruction. | - |


## Namelist ENERGYBALANCE

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lEB | .false. | .true. or .false. | Switch for using the facet energy balance. | - |
| lwriteEBfiles | .false. | .true., .false. | Switch for writing facet temperatures and energy budget to file. | - |
| lconstW | .false. | .true. or .false. | Switch whether soil moisture is assumed as constant in time (.true.) or the evaporated water is from the soil (.false.). | - |
| dtEB | 10. | `REAL` | Time interval between calculations of facet energy balance. | s |
| bldT | 0. | `REAL` | Internal temperature of the buildings. | [K] |
| flrT | 0. | `REAL` | Internal temperature of the ground. | [K] |
| wsoil | 0. | `REAL` | Water content of soil. | [kg/m3] |
| wgrmax | 450. | `REAL` | Maximum water content. | [kg/m3] |
| wwilt | 171. | `REAL` | Water content at wilting point. | [kg/m3] |
| wfc | 313. | `REAL` | Water content at field capacity. | [kg/m3] |
| skyLW | 0. | `REAL` | Long-wave radiation from the sky. | [W/m2] |
| GRLAI | 2. | `REAL` | Leaf area index of a green roof. | [m2/m2] |
| rsmin | 110. | `REAL` | Minimum resistance of soil/plant. | [s/m] |
| nfaclyrs | 3 | `INTEGER` | Number of layers making up each facet (nwalllayers in uDALES v1). | - |
| lvfsparse | .false. | .true. or .false. | Switch for view factors in sparse (text) format. | - |
| nnz | 0 | `INTEGER` | Number of non-zero view factors (only used with sparse view factor format. | - |
| lperiodicEBcorr | .false. | .true. or .false. | Switch for preventing over-heating and moisture saturation in periodic simualtions. | - |
| sinkbase | 0 | `INTEGER` > 0  | k index above which the periodicEBcorr sink is applied (should be above height of tallest building). | - |
| fraction | 0 | `REAL` > 0  | Ratio of domain height to uncapped boundary layer height. | - |

## Namelist SCALARS

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| nsv | 0 | `INTEGER` > 0 | Number of passive scalars. | - |
| lreadscal | .false. | .true., .false. | Switch for reading scalar pollutant field (warm start). *Deprecated, will be removed in the future.* | - |
| lscasrcr | .false. | .true., .false. |  Switch for 2-D network of point sources at lowest level as defined in scals.inp.xxx.  | - |
| lscasrcl | .false. | .true., .false. |  Switch for passive scalar line source when using canyon geometry.  | - |
| lscasrc | .false. | .true., .false. |  Switch for passive scalar point source defined by xS,yS,zS,SS,sigS. | - |


## Namelist DRIVER

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| idriver | 0 | 0, 1, 2 | Options for running precursor driver simulations where \*driver\* files will be written (`= 1`) and reading a completed driver simulation as the inlet to a simulation (`= 2`). Default (`= 0`) will do neither. | - |
| tdriverstart | 0. | `REAL` | Time at which \*driver\* files start being written. In use for `idriver = 1`. | s |
| dtdriver | 0. | `REAL` | Timestep at which \*driver\* file planes are written. In use for `idriver = 1`. | s |
| iplane | - | `INTEGER` | Index of the position on the x-axis of the plane that will be written to \*driver\* files. In use for `idriver = 1`. | |
| driverstore | 0. | `INTEGER` | Number of timesteps (`idriver = 1`) to be written to \*driver\* files or (`idriver = 2`) contained in \*driver\* files to be read. | - |
| driverjobnr | - | - | Job number of the \*driver\* files to be read. These files should be copied into the experiments folder of the driven simulation. In use for `idriver = 2`. | - |
| lsdriver | .false. | .true., .false. |  Switch for reading scalar driver files. In use for `idriver = 2`. | - |

## Namelist OUTPUT

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lfielddump | .false. | .true. or .false. | Switch for instantaneous field output. | - |
| tfielddump | 10000. | | Output time for fields. | [s] |
| fieldvars | '' | Any of the given labels, several are seperated by a comma: u0,v0,w0,th,ql,qt,p0,s1,s2,s3,s4,s5 | Variable names of fields. | - |
| tsample | 5. | | Sample time for statistics. | [s] |
| tstatsdump | 10000. | | Output time for statistics. | [s] |
| ltdump | .false. | .true. or .false. | Switch to output time-averaged statistics. | - |
| lydump | .false. | .true. or .false. | Switch to output y-averaged statistics. *Only functional if x-direction is not parallelised.* | - |
| lytdump | .false. | .true. or .false. | Switch to output y- and time- averaged statistics. *Only functional if x-direction is not parallelised.* | - |
| lxydump | .false. | .true. or .false. | Switch to output x- and y- averaged statistics. | - |
| lxytdump | .false. | .true. or .false. | Switch to output x-, y- and time-averaged statistics. | - |
| lslicedump | .false. | .true. or .false. | Switch to output slices in the xy-plane. | - |
<!---
| ltkedump | .false. | .true. or .false. | *Not supported in the current version.* | - |
--->

## Namelist INPS

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| zsize          | -       | -               | Size of domain in z direction.                                                            | [m]     |
| lzstretch      | false   | true or false   | Switch for stretched z grid.                                                              | -           |
| lstretchexp    | false   | true or false   | Switch for z grid stretched using exp function.                                           | -           |
| lstretchtanh   | false   | true or false   | Switch for z grid stretched using tanh function.                                          | -           |
| lstretch2tanh  | false   | true or false   | Switch for z grid stretched using 2tanh function.                                         | -           |
| stretchconst   | 0.01    | -               | Stretch constant.                                                                         | -           |
| u0             | 0       | -               | Initial u-velocity. Also applied as geostrophic term where applicable.                    | [m/s]         |
| v0             | 0       | -               | Initial v-velocity. Also applied as geostrophic term where applicable.                    | [m/s]         |
| dpdx           | 0       | -               | Pressure gradient in x direction.                                                         | [Pa/m]        |
| dpdy           | 0       | -               | Pressure gradient in y direction.                                                         | [Pa/m]        |
| thl0           | 288     | -               | Temperature at z = 0.                                                                     | -           |
| qt0            | 0       | -               | Specific humidity at z = 0.                                                               | -           |
| lapse          | 0       | -               | Lapse rate.                                                                               | [K/m]         |
| w_s            | 0       | -               | Subsidence.                                                                               | -           |
| R              | 0       | -               | Radiative forcing.                                                                        | [W/m^2]       |
| xS | 0 | `REAL` | Position of scalar point source in x. | [m] |
| yS | 0 | `REAL` | Position of scalar point source in y. | [m] |
| zS | 0 | `REAL` | Position of scalar point source in z. | [m] |
| SSp | 0. | `REAL` | Strength of scalar point source. | [g/s] |
| sigSp | 0. | `REAL` | Standard deviation of scalar point source. | [m] |
| xSb | 0 | `REAL` | Start of scalar line source in x. | [m] |
| ySb | 0 | `REAL` | Start of scalar line source in y. | [m] |
| zSb | 0 | `REAL` | Start of scalar line source in z. | [m] |
| xSe | 0 | `REAL` | End of scalar line source in x. | [m] |
| ySe | 0 | `REAL` | End of scalar line source in y. | [m] |
| zSe | 0 | `REAL` | End of scalar line source in z. | [m] |
| SSl | 0. | `REAL` | Strength of scalar line source. | [g/ms] |
| sigSl | 0. | `REAL` | Standard deviation of scalar line source. | [m] |
| NOb            | 0       | -               | Initial concentration of NO.                                                              | -           |
| NO2b           | 0       | -               | Initial concentration of NO2.                                                             | -           |
| O3b            | 0       | -               | Initial concentration of O3.                                                              | -           |
| stl_file | - | `STRING` | Name of STL file defining the geometry. | - |
| read_types | false | true or false | Switch for reading facet types from file. Default: false (all facets are set to type 1). | - |
| types_path | - | `STRING` | Name of types file. | - |
| facT | 288 | `REAL` > 0 | If `iwallmom = 2` or  `iwalltemp = 2` then this sets the facet temperature, or if `lEB = .true.` it is the initial facet temperature | [K] |
| ifacsec | 1 | 1 or 2 | Option for facet section calculation. 1: Fortran (default, fast), 2: MATLAB (useful for debugging | - |
| ishortwave   | 1        | 1  or 2 | Option for shortwave radiation calculation,  1  uses Fortran and is faster, 2 uses MATLAB and useful for debugging. | -  |
| isolar       | 1        | 1 , 2, 3| Option for solar radiation, 1 uses custom values, 2 uses lattitude and lonigtude, 3 uses weather file.             | -  |
| view3d_out   | 0        | 0 , 1 , 2                    | Output format for View3D, 0 is text, 2 is binary, 2 is sparse.                   | - |
| maxD         | Inf      |   `REAL` > 0   | Maximum distance to check view factors, otherwise they are zero.    | - |
| xazimuth     | 90       |  `REAL`    | The azimuthal angle of the x-axis (with respect to North).     | [degrees]  |
| solarazimuth | 135       |  `REAL`            | Solar azimuth, used if isolar = 1.                                   | [degrees] |
| solarzenith  | 28.4066   |  `REAL`               | Solar zenith, used if isolar = 1 .                                  | [degrees] |
| I            | 800       |  `REAL` > 0               | Shortwave direct normal irradiance (DNI), used if isolar = 1.                  | W/m^2   |
| Dsky         | 418.8041  |  `REAL` > 0               | Diffuse sky irradiance, used if isolar = 1.                          | W/m^2   |
| year      | -       |   `INTEGER`   | Year.                                               | -     |
| month     | -       | 1 <= `INTEGER` <= 12         | Month (where 6 corresponds to June), if isoloar = 2 or 3.                | -     |
| day       | -       |  1 <= `INTEGER` <= 31       | Day, if isloar = 2 or 3.                                                | -     |
| hour      | 6       | 0 <= `INTEGER` <= 23        | Hour (0 for midnight, 23 for 11pm), if isloar = 2 or 3.                 | [hours]     |
| minute    | 0       | 0 <= `INTEGER` <= 59       | Minute, if isloar = 2 or 3.                                             | [mins]     |
| second    | 0       | 0 <= `INTEGER` <= 59         | Second, if isloar = 2 or 3.                                             | [s]     |
| longitude | -0.13   |   `REAL`            | Longitude, if isloar = 2.                                          | [degrees]     |
| latitude  | 51.5    |  `REAL`            | Latitude, if isloar = 2.                                           | [degrees]     |
| timezone  | 0       | -               | Timezone, if isloar = 2.                                           | -     |
| elevation | 0       |  `REAL`         | Elevation, if isloar = 2.                                          | -     |
| weatherfname | - | `STRING` | File containing weather data, if isolar = 3. | - |


## Namelist CHEMISTRY
<!--
*This section will be updated with the next version.*
-->

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lchem | .false. | .true., .false. | Switch for basic chemistry. | - |
| k1 | 0. | | Rate constant (O3 + NO -> NO2 + 02 ). Chemistry model parameter. | |
| JNO2 | 0. | | NO2 photolysis rate. Chemistry model parameter. | |

<!---
## Namelist INLET

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| Uinf | 0. | `REAL` | Fixed velocity at domain top (x-direction). | m/s |
| Vinf | 0. | `REAL` | Fixed velocity at domain top (y-direction). | m/s |
-->
