# Namoptions overview

This list refers to the original code-base [DALES](https://github.com/dalesteam/dales). The latest version of the namoptions overview of DALES is documented [here](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf).

## Namelist DOMAIN

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| itot | 64 | | Number of points in x-direction (imax in uDALES v1). | |
| jtot | 64 | | Number of points in y-direction. | |
| ktot | 96 | | Number of points in z-direction (kmax in uDALES v1).| |
| xlen | -1 | | Domain size in x-direction (xsize in uDALES v1).| |
| ylen | -1 | | Domain size in x-direction (ysize in uDALES v1).| |

## Namelist DYNAMICS

Possible advection schemes:

1 = 1st order upwind scheme

2 = 2nd order central difference scheme

7 = Kappa (flux limited) scheme. This scheme designed for quantities that should never become negative.

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lqlnr | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default switched from .true. | |
| iadv_mom | 2 | 2 | Advection scheme for momentum. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| iadv_tke | -1 | 2 | Advection scheme for TKE. Only used if `loneeqn = True`. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| iadv_thl | -1 | 2, 7 | Advection scheme for temperature. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| iadv_qt | -1 | 2 | Advection scheme for moisture. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| iadv_sv | -1 | 1, 2, 7 | Advection scheme for scalars. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| ipoiss | 0 | 0 | Poisson solver. 0 = Fast Fourier Transform. | - |

## Namelist PHYSICS

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| ps | -1 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lcoriol | .false. | .true., .false. | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default switched to .false. | - |
| igrw_damp | 2 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lmoist | .false. | .true., .false. | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default switched to .false. | - |
| ltempeq | .false. | .true., .false. | Switch for solving temperature equation. | - |
| lbuoyancy | .false. | .true., .false. | Switch for buoyancy force in temperature equation. | - |
| lprofforc | .false. | .true., .false. | Switch for nudging flow to a profile (forcing). | - |
| luoutflowr | .false. | .true., .false. | Switch that determines whether u-velocity is corrected to get a fixed outflow rate. *Only functional when x-direction is not parellelised.* | |
| lvoutflowr | .false. | .true., .false. | Switch that determines whether u-velocity is corrected to get a fixed outflow rate. *Only functional when y-direction is not parellelised.* | |
| luvolflowr | .false. | .true., .false. | Switch that determines whether u-velocity is corrected to get a fixed volume flow rate | |
| lvvolflowr | .false. | .true., .false. | Switch that determines whether u-velocity is corrected to get a fixed volume flow rate | |
| uflowrate | 1. | `REAL` | U-velocity flow rate for out- or volume-flow forcing. | m/s |
| vflowrate | 1. | `REAL` | V-velocity flow rate out- or volume-flow forcing. | m/s |
| ifixuinf | 0 | 1, 2 | Choice for free stream forcing. (0 = nothing) | |
| lvinf | .false. | .true., .false. | Use Vinf instead of Uinf for the fixed velocity at infinity | |
| tscale | | | Timescale: domain height*Uinf/utau\*\*2 | |
| lnudge | .false. | .true., .false. | Switch for nudging flow to profiles (in `prof.inp.xxx`) | |
| nnudge | 0 | `INTEGER` | Number of points from bottom to nudge. | |
| tnudge | 60. | `REAL` | Time scale for nudging | |
| dpdx | 0. | | Constant pressure gradient forcing in x. | |
| ltimedepsurf | .false. | .true., .false. | Switch for time-dependent surface heat flux (`bctfz` etc - see BC section) | |
| ntimedepsurf | 0 | `REAL` | Number of time-dependent surface heat fluxes in file `timedepsurf.inp.xxx` | |
| ltimedepnudge | .false. | .true., .false. | Switch for time-dependent profiles. When `lnudge = .true.`, then this switch causes the nudging profile to vary in time. If using inflow-outflow boundary conditions with inflow given by profile (`BCxm/BCxT/BCxq = 2`), then the inflow profile varies in time. | |
| ntimedepnudge | 0 | `REAL` | Number of time-dependent nudging profiles in file `timedepnudge.inp.xxx` | |
| ltimedepsw | .false. | .true., .false. | Time-dependent shortwave radiation on facets | |
| ntimedepsw | 0 | `REAL` | Number of time-dependent shortwave radiative fluxes in file `timedepsw.inp.xxx` | |
| ltimedeplw | .false. | .true., .false. | Time-dependent longwave radiation on facets | |
| ntimedeplw | 0 | `REAL` | Number of time-dependent shortwave radiative fluxes in file `timedeplw.inp.xxx` | |

## Namelist RUN

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| iexpnr | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| runtime | 300 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| dtmax | 20 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lwarmstart | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| startfile | '' | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| trestart | 10000. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| irandom | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| krand | | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Setting no value will return kmax. | |
| randu | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from 0.5 | |
| randthl | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from 0.1 | |
| randqt | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from 1e-5. | |
| ladaptive | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| courant | -1 | | Default sets it to 1.5 or 1.1 (if Kappa or upwind scheme is used). These are different values than in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| author | '' | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lles | .true. | .true., .false. | Switch that determines whether the subgrid model is turned on or constant ekm and ekh are used (DNS) | - |
| libm | .true. | | Switch that determines whether the Immersed Boundary Method is turned on. | |
| lrandomize | .true. | | Switch that determines whether initial field is randomised. *Currently not independent of domain decomposition.* | 

## Namelist OUTPUT

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lfielddump | .false. | .true., .false. | Switch for instantaneous field output. | - |
| tfielddump | 10000. | | Output time for fields. | s |
| fieldvars | '' | Any of the given labels, several are seperated by a comma: u0,v0,w0,th,ql,qt,p0,s1,s2,s3,s4,s5 | Variable names of fields. | - |
| tsample | 5. | | Sample time for statistics. | s |
| tstatsdump | 10000. | | Output time for statistics. | s |
| ltdump | .false. | .true., .false. | Switch to output time-averaged statistics. | - |
| lydump | .false. | .true., .false. | Switch to output y-averaged statistics. *Only functional if x-direction is not parallelised.* | - |
| lytdump | .false. | .true., .false. | Switch to output y- and time- averaged statistics. *Only functional if x-direction is not parallelised.* | - |
| lxydump | .false. | .true., .false. | Switch to output x- and y- averaged statistics. | - |
| lxytdump | .false. | .true., .false. | Switch to output x-, y- and time-averaged statistics. | - |
| lslicedump | .false. | .true., .false. | Switch to output slices in the xy-plane. | - |
| ltkedump | .false. | .true., .false. | *Not supported in the current version.* | - |

## Namelist NAMSUBGRID

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lvreman | .false. | .true., .false. | Switch for Vreman (2004) sub-grid scheme. | - |
| c_vreman | 0.07 | | Model constant for Vreman scheme. | |
| lbuoycorr | .false. | .true., .false. | Switch for buoyancy correlation in the Vreman scheme. | - |
| loneeqn | .false. | .true., .false. | Switch for one-equation sub-grid scheme. | - |
| ldelta | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lmason | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| cf | 2.5 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| cn | 0.76 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| Rigc | 0.25 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| Prandtl | 0.333 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lsmagorinsky | .false . | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| cs | -1 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| nmason | 2 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |

## Namelist NAMCHECKSIM

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| tcheck | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |

## Namelist BC

Switches for boundary conditions: momentum (m), temperature (T), humidity (q) and scalars (s).

Lateral BCs (BCx, BCy): 1 = periodic, 2 in/outflow conditions, inflow given by profile (usually constant), 3: in/outflow conditions, inflow given by precursor simulation.

BCs at the top (BCtop): 1 = freeslip, 2 = noslip, 3 = should be used with inflow/outflow conditions.

BCs at the bottom (BCbot; only effective if not covered with ground facets): 1 = flux, 2 = wall function, 3 = neutral wall function.

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| BCxm | 1 | 1,2,3 | Domain boundary condition for momentum in x. | |
| BCxT | 1 | 1,2,3 | Domain boundary condition for temperature in x. | |
| BCxq | 1 | 1,2,3 | Domain boundary condition for humidity in x. | |
| BCxs | 1 | 1,2,3 | Domain boundary condition for scalars in x. | |
| BCym | 1 | 1,2 | Domain boundary condition for momentum in y. | |
| BCyT | 1 | 1,2 | Domain boundary condition for temperature in y. | |
| BCyq | 1 | 1,2 | Domain boundary condition for humidity in y. | |
| BCys | 1 | 1,2 | Domain boundary condition for scalars in y. | |
| BCtopm | 1 | 1,2,3 | Boundary condition for momentum at domain top. | |
| BCtopT | 1 | 1,2 | Boundary condition for temperature at domain top. | |
| BCtopq | 1 | 1,2 | Boundary condition for humidity at domain top. | |
| BCtops | 1 | 1,2 | Boundary condition for scalars at domain top. | |
| bctfxm | 0 | `REAL` | Temperature flux on facets with surface normal in -x direction. | |
| bctfxp | 0 | `REAL` | Temperature flux on facets with surface normal in +x direction. | |
| bctfym | 0 | `REAL` | Temperature flux on facets with surface normal in -y direction. | |
| bctfyp | 0 | `REAL` | Temperature flux on facets with surface normal in +y direction. | |
| bctfz | 0 | `REAL` | Temperature flux on facets with surface normal in +z direction. | |
| bcqfxm | 0 | `REAL` | Moisture flux on facets with surface normal in -x direction. | |
| bcqfxp | 0 | `REAL` | Moisture flux on facets with surface normal in +x direction. | |
| bcqfym | 0 | `REAL` | Moisture flux on facets with surface normal in -y direction. | |
| bcqfyp | 0 | `REAL` | Moisture flux on facets with surface normal in +y direction. | |
| bcqfz | 0 | `REAL` | Moisture flux on facets with surface normal in +z direction. | |
| thl_top | -1. | | Temperature at the top boundary. | |
| qt_top | -1. | | Humidity at the top boundary. | |
| wttop | 0. | | Temperature flux at the top boundary. | |
| wsvtopdum | | | Scalar boundary conditions top. | |
| BCbotm | 2 | 1,2,3 | Boundary condition for momentum at domain bottom (if `lbottom = .true.`). | |
| BCbotT | 1 | 1,2 | Boundary condition for temperature at domain bottom (if `lbottom = .true.`). | |
| BCbotq | 1 | 1 | Boundary condition for humidity at domain bottom (if `lbottom = .true.`). | |
| BCbots | 1 | 1 | Boundary condition for scalars at domain bottom (if `lbottom = .true.`). | |
| wtsurf | -1. | | Temperature flux at domain bottom (if `lbottom = .true.`). *Currently need to be set to reasonable values for subroutine bottom.* | |
| wqsurf | -1. | | Moisture flux at domain bottom (if `lbottom = .true.`). *Currently need to be set to reasonable values for subroutine bottom.* | |
| wsvsurfdum | | | Scalar flux at domain bottom (if `lbottom = .true.`). | |
| thls | -1. | | Temperature at domain bottom (if `lbottom = .true.`). *Currently need to be set to reasonable values for subroutine bottom.* | |
| qts | -1. | | Moisture at domain bottom (if `lbottom = .true.`). Used in modthermodynamics to get a BC for the moisture profile. | |
| z0 | -1. | | Momentum roughness length of the domain bottom (if `lbottom = .true.`). *Currently need to be set to reasonable values for subroutine bottom.* | |
| z0h | -1. | | Heat roughness length of the domain bottom (if `lbottom = .true.`). *Currently need to be set to reasonable values for subroutine bottom.* | |

## Namelist ENERGYBALANCE

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lEB | .false. | .true., .false. | Switch for using the facet energy balance. | - |
| lwriteEBfiles | .false. | .true., .false. | Switch for writing facet temperatures and energy budget to file. | - |
| lconstW | .false. | .true., .false. | Switch whether soil moisture is assumed as constant in time (.true.) or the evaporated water is from the soil (.false.). | - |
| dtEB | 10. | `REAL` | Time interval between calculations of facet energy balance. | s |
| bldT | 0. | `REAL` | Internal temperature of the buildings. | K |
| flrT | 0. | `REAL` | Internal temperature of the ground. | K |
| wsoil | 0. | `REAL` | Water content of soil. | kg/m3 |
| wgrmax | 450. | `REAL` | Maximum water content. | kg/m3 |
| wwilt | 171. | `REAL` | Water content at wilting point. | kg/m3 |
| wfc | 313. | `REAL` | Water content at field capacity. | kg/m3 |
| skyLW | 0. | `REAL` | Long-wave radiation from the sky. | |
| GRLAI | 2. | `REAL` | Leaf area index of a green roof. | |
| rsmin | 110. | `REAL` | Minimum resistance of soil/plant. | |
| nfaclyrs | 3 | `INTEGER` | Number of layers making up each facet (nwalllayers in uDALES v1). | |
| lvfsparse | .false. | .true., .false. | Switch for view factors in sparse (text) format. | |
| nnz | 0 | `INTEGER` | Number of non-zero view factors (only used with sparse view factor format | |

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
| lnorec | .false. | .true., .false. | Switch for not using reconstruction. | - |

## Namelist SCALARS

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| nsv | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lreadscal | .false. | .true., .false. | Switch for reading scalar pollutant field (warm start). *Deprecated, will be removed in the future.* | - |
| lscasrcr | .false. | .true., .false. |  Switch for 2-D network of point sources at lowest level as defined in scals.inp.xxx.  | - |
| lscasrcl | .false. | .true., .false. |  Switch for passive scalar line source when using canyon geometry.  | - |
| lscasrc | .false. | .true., .false. |  Switch for passive scalar point source defined by xS,yS,zS,SS,sigS. | - |
| xS | 0 | `REAL` | Position of scalar source in x. | m |
| yS | 0 | `REAL` | Position of scalar source in y. | m |
| zS | 0 | `REAL` | Position of scalar source in z. | m |
| SS | 0. | `REAL` | Strength of scalar source. | g/ms |
| sigS | 0. | `REAL` | Standard deviation of scalar source. | m |

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

## Namelist CHEMISTRY

*This section will be updated with the next version.*

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lchem | .false. | .true., .false. | Switch for basic chemistry. | - |
| k1 | 0. | | Rate constant (O3 + NO -> NO2 + 02 ). Chemistry model parameter. | |
| JNO2 | 0. | | NO2 photolysis rate. Chemistry model parameter. | |

## Namelist INLET

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| Uinf | 0. | `REAL` | Fixed velocity at domain top (x-direction). | m/s |
| Vinf | 0. | `REAL` | Fixed velocity at domain top (y-direction). | m/s |
