# Namoptions overview

This list refers to the original code-base [DALES](https://github.com/dalesteam/dales). The latest version of the namoptions overview of DALES is documented [here](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf).

## Namelist DOMAIN

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| imax | 64 | | See `itot` in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| jtot | 64 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| kmax | 96 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| xsize | -1 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| ysize | -1 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| xlat | 52. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| xlon | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| xday | 1. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| xtime | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| ksp | -1 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). (Setting to -1 calculates default value) | |

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
| ipoiss | 1 | 0, 1 | Poisson solver. 0 = Fast Fourier Transformation, 1 = Cyclic reduction scheme. *Default will change to 0 in the future.* | - |

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
| luoutflowr | .false. | | switch that determines whether u-velocity is corrected to get a fixed outflow rate | |
| lvoutflowr | .false. | | switch that determines whether u-velocity is corrected to get a fixed outflow rate | |
| luvolflowr | .false. | | switch that determines whether u-velocity is corrected to get a fixed volume flow rate | |
| lvvolflowr | .false. | | switch that determines whether u-velocity is corrected to get a fixed volume flow rate | |
| uflowrate | 1. | `REAL` | U-velocity flow rate for out- or volume-flow forcing. | m/s |
| vflowrate | 1. | `REAL` | V-velocity flow rate out- or volume-flow forcing. | m/s |
| ifixuinf | 0 | 1, 2 | Choice for free stream forcing. (0 = nothing) | |
| lvinf | .false. | | use Vinf instead of Uinf for the fixed velocity at infinity | |
| tscale | | | timescale: domain height*Uinf/utau\*\*2 | |
| lnudge | .false. | | switch for applying nudging at the top of the domain | |
| tnudge | 50. | | time scale for nudging | |
| nnudge | 10 | | | |
| dpdx | 0. | | Constant pressure gradient forcing in x. | |

## Namelist RUN

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| iexpnr | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| runtime | 300 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| dtmax | 20 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lwarmstart | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lper2inout | .false. | .true., .false. | Switch that determines type of restart: .true. means switching from periodic to in/outflow: inlet profile is read from `prof.inp`. *Potentially deprecated. May be removed in the future.* | |
| startfile | '' | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| **lstratstart** | .false. | .true., .false. | *Description missing* | |
| trestart | 10000. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| irandom | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| krand | | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Setting no value will return kmax. | |
| randu | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from 0.5 | |
| randthl | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from 0.1 | |
| randqt | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from 1e-5. | |
| ladaptive | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| courant | -1 | | Default sets it to 1.5 or 1.1 (if Kappa or upwind scheme is used). These are different values than in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| **diffnr** | 0.25 | | Diffusion number? Used to determine adaptive time step. | |
| author | '' | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lles | .true. | .true., .false. | Switch that determines whether the subgrid model is turned on or constant ekm and ekh are used (DNS) | - |
| libm | .true. | | Switch that determines whether the Immersed Boundary Method is turned on. *Deprecated. Will be removed in the future.* | |
| lreadmean | .false. | | Switch that determines whether mean variables should be read from means#myid#.#expnr# *Potentially deprecated. May be removed in the future.* | |
| lwalldist | .false. | | Switch that determines whether the wall distances should be computed for the subgrid models. *Potentially deprecated. May be removed in the future.* | |

## Namelist OUTPUT

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lfielddump | .false. | .true., .false. | Switch for instantaneous field output. | - |
| tfielddump | 10000. | | Output time for fields. | s |
| fieldvars | '' | Any of the given labels, several are seperated by a comma: u0,v0,w0,th,ql,qt,p0,s1,s2,s3,s4,s5 | Variable names of fields. | - |
| tsample | 5. | | Sample time for statistics. | s |
| tstatsdump | 10000. | | Output time for statistics. | s |
| ltdump | .false. | .true., .false. | Switch to output time-averaged statistics. | - |
| lydump | .false. | .true., .false. | Switch to output y-averaged statistics. | - |
| lytdump | .false. | .true., .false. | Switch to output y- and time- averaged statistics. | - |
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

Lateral BCs (BCx, BCy): 1 = periodic, > 1 special in/outflow conditions

BCs at the top (BCtop): 1 = freeslip, 2 = noslip, 3 = determined by inflow conditions

BCs at the bottom (BCbot; only effective if not covered with road facets): 1 = flux, 2 = wall function

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| BCxm | 1 | | Domain boundary condition for momentum in x. | |
| BCxT | 1 | | Domain boundary condition for temperature in x. | |
| BCxq | 1 | | Domain boundary condition for humidity in x. | |
| BCxs | 1 | | Domain boundary condition for scalars in x. | |
| BCym | 1 | | Domain boundary condition for momentum in y. | |
| BCyT | 1 | | Domain boundary condition for temperature in y. | |
| BCyq | 1 | | Domain boundary condition for humidity in y. | |
| BCys | 1 | | Domain boundary condition for scalars in y. | |
| BCtopm | 1 | | Boundary condition for momentum at domain top. | |
| BCtopT | 1 | | Boundary condition for temperature at domain top. | |
| BCtopq | 1 | | Boundary condition for humidity at domain top. | |
| BCtops | 1 | | Boundary condition for scalars at domain top. | |
| BCbotm | 2 | | Boundary condition for momentum at domain bottom. | |
| BCbotT | 1 | | Boundary condition for temperature at domain bottom. | |
| BCbotq | 1 | | Boundary condition for humidity at domain bottom. | |
| BCbots | 1 | | Boundary condition for scalars at domain bottom. | |
| bctfxm | 0. | | Bounary Condition Temperature Flux X-minus-wall. | |
| bctfxp | 0. | | Bounary Condition Temperature Flux X-plus-wall. | |
| bctfym | 0. | | Bounary Condition Temperature Flux Y-minus-wall. | |
| bctfyp | 0. | | Bounary Condition Temperature Flux y-plus-wall. | |
| bctfz | 0. | | Bounary Condition Temperature Flux z top-wall. | |
| thl_top | -1. | | Temperature at the top boundary. | |
| qt_top | -1. | | Humidity at the top boundary. | |
| wttop | 0. | | Temperature flux at the top boundary. | |
| qts | -1. | | Used in modthermodynamics to get a BC for the moisture profile. | |
| wsvsurfdum | | | Scalar boundary conditions bottom. | |
| wsvtopdum | | | Scalar boundary conditions top. | |
| wtsurf | -1. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). *Currently need to be set to reasonable values for subroutine bottom.* | |
| wqsurf | -1. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). *Currently need to be set to reasonable values for subroutine bottom.* | |
| thls | -1. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). *Currently need to be set to reasonable values for subroutine bottom.* | |
| z0 | -1. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). *Currently need to be set to reasonable values for subroutine bottom.* | |
| z0h | -1. | | Facet roughness length for heat. *Currently need to be set to reasonable values for subroutine bottom.* | |

## Namelist ENERGYBALANCE

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lEB | .false. | .true., .false. | Switch for using the facet energy balance. | - |
| lwriteEBfiles | .false. | .true., .false. | Switch for writing facet temperatures and energy budget to file. | - |
| lconstW | .false. | .true., .false. | Switch whether soil moisture is assumed as constant in time (.true.) or the evaporated water is from the soil (.false.). | - |
| dtEB | 10. | `REAL` | Time interval between calculations of facet energy balance. | s |
| bldT | 0. | `REAL` | Internal temperature of the buildings, currently also ground temperature at a depth equal to floor facet thickness. | K |
| wsoil | 0. | `REAL` | Water content of soil. | kg/m3 |
| wgrmax | 450. | `REAL` | Maximum water content. | kg/m3 |
| wwilt | 171. | `REAL` | Water content at wilting point. | kg/m3 |
| wfc | 313. | `REAL` | Water content at field capacity. | kg/m3 |
| skyLW | 0. | `REAL` | Long-wave radiation from the sky. | |
| GRLAI | 2. | `REAL` | Leaf area index of a green roof. | |
| rsmin | 110. | `REAL` | Minimum resistance of soil/plant. | |

## Namelist WALLS

1 = fixed flux

2 = flux determined by wall function involving temperature

3 = flux determined by neutral wall function (set automatically if `ltempeq = .false.`)

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| nblocks | 0 | `INTEGER` | Number of blocks specified in `blocks.inp`. | - |
| nfcts | -1 | `INTEGER` | Number of facets specified in `facets.inp`. | - |
| iwallmom | 2 | 2, 3 (1 currently not implemented) | Building wall momentum flux. | - |
| iwalltemp | 1 | 1, 2 |  Building wall temperature flux. | - |
| iwallmoist | 1 | 1, 2 |  Building wall moisture flux. | - |
| iwallscal | 1 | 1, 2 | Building wall scalar flux | - |

## Namelist SCALARS

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| nsv | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lreadscal | .false. | .true., .false. | Switch for reading scalar pollutant field (warm start). *Deprecated, will be removed in the future.* | - |
| **lscasrc** | .false. | .true., .false. |  *Description missing* | |
| **lscasrcl** | .false. | .true., .false. |  *Description missing* | |
| xS | 0 | | | |
| yS | 0 | | | |
| zS | 0 | | | |
| SS | 0. | | | |
| sigS | 0. | | | |

## Namelist DRIVER

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| idriver | 0 | 0, 1, 2 | Options for running precursor driver simulations where \*driver\* files will be written (`= 1`) and reading a completed driver simulation as the inlet to a simulation (`= 2`). Default (`= 0`) will do neither. | - |
| tdriverstart | 0. | `REAL` | Time at which \*driver\* files start being written. In use for `idriver = 1`. | s |
| dtdriver | 0. | `REAL` | Timestep at which \*driver\* file planes are written. In use for `idriver = 1`. | s |
| iplane | - | `INTEGER` | Index of the position on the x-axis of the plane that will be written to \*driver\* files. In use for `idriver = 1`. | |
| driverstore | 0. | `INTEGER` | Number of timesteps (`idriver = 1`) to be written to \*driver\* files or (`idriver = 2`) contained in \*driver\* files to be read. | - |
| driverjobnr | - | - | Job number of the \*driver\* files to be read. These files should be copied into the experiments folder of the driven simulation. In use for `idriver = 2`. | - |

## Namelist CHEMISTRY

*This section will be updated with the next version.*

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lchem | .false. | .true., .false. | Switch for basic chemistry. | - |
| k1 | 0. | | Rate constant (O3 + NO -> NO2 + 02 ). Chemistry model parameter. | |
| JNO2 | 0. | | NO2 photolysis rate. Chemistry model parameter. | |

## Namelist INLET

*This section will be updated with the next version.*

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| Uinf | 0. | | | |
| Vinf | 0. | | | |
| inletav | 0. | | | |
| lstoreplane | .false. | | | |
| linletRA | .false. | | | |
| lfixinlet | .false. | | | |
| lfixutauin | .false. | | | |
| lreadminl | .false. | | | |
| di | 0.09 | | | |
| dti | | | | |
| lwallfunc | .true. | | Switch that determines whether wall functions are used to compute the wall-shear stress. *Deprecated, only in use in modinlet. Will be removed in the future.* | |
