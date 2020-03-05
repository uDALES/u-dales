This list refers to the original code-base [DALES](https://github.com/dalesteam/dales). The latest version of the namoptions overview of DALES is documented [here](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf).

# Namelist DOMAIN

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
| nblocks | 0 | `INTEGER` | Number of blocks specified in `blocks.inp`. | - |
| **ntrees** | 0 | `INTEGER` | Number of trees specified in `trees.inp`. (*Not supported in current version.*) | - |
| **npurif** | 0 | `INTEGER` | Number of air purifiers specified in `purifs.inp`. (*Not supported in current version.*) | - |


# Namelist DYNAMICS

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
| **ipoiss** | 1 | 0, 1 | Poisson solver. 0 = Fast Fourier Transformation, 1 = Cyclic reduction scheme. *Should we change default to 0?* | - |


# Namelist PHYSICS

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| ps | -1 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lcoriol | .false. | .true., .false. | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default switched to .false. | - |
| igrw_damp | 2 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lmoist | .false. | .true., .false. | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default switched to .false. | - |
| ltempeq | .false. | .true., .false. | Switch for solving temperature equation. | - |
| lbuoyancy | .false. | .true., .false. | Switch for buoyancy force in temperature equation. | - |
| uflowrate | 1. | `REAL` | U-velocity flow rate for out- or volume-flow forcing. | m/s |
| vflowrate | 1. | `REAL` | V-velocity flow rate out- or volume-flow forcing. | m/s |
| lprofforc | .false. | .true., .false. | Switch for nudging flow to a profile (forcing). | - |
| sun | 0. | | Tree model parameter. (*Not supported in current version.*) | |
| Bowen | 0. | | Tree model parameter. (*Not supported in current version.*) | |
| cd | 0. | | Tree model parameter. (*Not supported in current version.*) | |
| decay | 0. | | Tree model parameter. (*Not supported in current version.*) | |
| ud | 0. | | Tree model parameter. (*Not supported in current version.*) | |
| lchem | .false. | | Switch for basic chemistry. (*Not supported in current version.*) | |
| k1 | 0. | | Rate constant (O3 + NO -> NO2 + 02 ). Chemistry model parameter. (*Not supported in current version.*) | |
| JNO2 | 0. | | NO2 photolysis rate. Chemistry model parameter. (*Not supported in current version.*) | |
| Qpu | 0. | | Purifiers flow rate. (*Not supported in current version.*) | |
| epu | 0. | | Purifiers efficiency. (*Not supported in current version.*) | |
| **z0** | 0.1 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from -1. | |
| **z0h** | 0.1 | | | |


# Namelist RUN

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| iexpnr | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| runtime | 300 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| dtmax | 20 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lwarmstart | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| startfile | '' | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| **trestart** | | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| irandom | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| krand | *kmax* | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Setting now value will return kmax. | |
| **randu** | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from 0.5 | |
| **randthl** | 0. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from 0.1 | |
| **randqt** | | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default changed from 1e-5. | |
| nsv | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| ladaptive | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| courant | -1 | | Default sets it to 1.5 or 1.1 (if Kappa or upwind scheme is used). These are different values than in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| author | '' | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |
| lreadscal | .false. | .true., .false. | Switch for reading scalar pollutant field (warm start) *Is this working?* | - |
| lstratstart | .false. | .true., .false. | *Description missing* | |
| lscasrc | .false. | .true., .false. |  *Description missing* | |
| lscasrcl | .false. | .true., .false. |  *Description missing* | |
| lper2inout | .false. | .true., .false. | Switch that determines type of restart: .true. means switching from periodic to in/outflow: inlet profile is read from `prof.inp`.  *Is this working?* | |
| lles | .true. | .true., .false. | Switch that determines wether the subgrid model is turned on or constant ekm and ekh are used (DNS) | - |
| diffnr | 0.25 | | | |
| lnudge | .false. | | switch for applying nudging at the top of the domain | |
| tnudge | 50. | | time scale for nudging | |
| nnudge | 10 | | | |
| lwallfunc | .true. | | switch that determines wether wall functions are used to compute the wall-shear stress | |
| lreadmean | .false. | | switch that determines wether mean variables should be read from means#myid#.#expnr# | |
| startmean | | | | |
| lwalldist | .false. | | switch that determines wether the wall distances should be computed | |
| dpdx | 0. | | *Does this still work?* | |
| libm | .true. | | switch that determines wether the Immersed Boundary Method is turned on | |
| ltrees | .false. | | | |
| lpurif | .false. | | | |
| luoutflowr | .false. | | switch that determines wether u-velocity is corrected to get a fixed outflow rate | |
| lvoutflowr | .false. | | switch that determines wether u-velocity is corrected to get a fixed outflow rate | |
| luvolflowr | .false. | | switch that determines wether u-velocity is corrected to get a fixed volume flow rate | |
| lvvolflowr | .false. | | switch that determines wether u-velocity is corrected to get a fixed volume flow rate | |
| ifixuinf | 0 | | | |
| lvinf | .false. | | use Vinf instead of Uinf for the fixed velocity at infinity | |
| tscale | | | timescale: domain height*Uinf/utau**2 | |
| lfielddump | .true. | | Switch for instantaneous field output. | |
| tfielddump | | | Output time for fields. | |
| fieldvars | '' | Any of the given labels, several are seperated by a comma: u0,v0,w0,th,ql,qt,p0,s1,s2,s3,s4,s5 | Variable names of fields. | |
| tsample | | | Sample time for statistics. | |
| tstatsdump | | Output time for statistics. | | |
| ltdump | .false. | | Switch to output time-averaged statistics. | |
| lydump | .false. | | Switch to output y-averaged statistics. | |
| lytdump | .false. | | Switch to output y- and time- averaged statistics. | |
| lxydump | .false. | | Switch to output x- and y- averaged statistics. | |
| lxytdump | .false. | | Switch to output x-, y- and time-averaged statistics. | |
| lslicedump | .false. | | Switch to output slices in the xy-plane. | |
| ltkedump | .false. | | *Not supported in the current version.* | |


# Namelist NAMSUBGRID

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


# Namelist NAMCHECKSIM

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| tcheck | 0 | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | |


# Namelist BC

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| | | | | |


# Namelist ENERGYBALANCE

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lEB | .false. | .true., .false. | Switch for using the facet energy balance. | - |
| **lconstW** | .false. | .true., .false. | *Description missing.* | - |
| dtEB | 10. | `REAL` | Time interval between calculations of facet energy balance. | s |
| nfcts | -1 | `INTEGER` | Number of facets specified in `facets.inp`. | - |
| bldT | 0. | `REAL` | Internal temperature of the buildings, currently also ground temperature at a depth equal to floor facet thickness. | K |
| wsoil | 0. | `REAL` | Water content of soil. | kg/m3 |
| wgrmax | 450. | `REAL` | Maximum water content. | kg/m3 |
| wwilt | 171. | `REAL` | Water content at wilting point. | kg/m3 |
| wfc | 313. | `REAL` | Water content at field capacity. | kg/m3 |
| skyLW | 0. | `REAL` | Long-wave radiation from the sky. | |
| GRLAI | 2. | `REAL` | Leaf area index of a green roof. | |
| rsmin | 110. | `REAL` | Minimum resistance of soil/plant. | |


# Namelist WALLS

1 = fixed flux

2 = flux determined by wall function involving temperature

3 = flux determined by neutral wall function (no temperature)

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| **iwallmom** | 2 | 2, 3 (1 currently not implemented) | Building wall momentum flux. *SHould we change default to 3?* | - |
| iwalltemp | 1 | 1, 2 |  Building wall temperature flux. | - |
| iwallmoist | 1 | 1, 2 |  Building wall moisture flux. | - |


# Namelist INLET

*Work in Progress. This section currently only lists default values as defined in the code, their functionality needs to be checked!*

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| Uinf | 0. | | | |
| Vinf | 0. | | | |
| inletav | 0. | | | |
| xS | 0 | | | |
| yS | 0 | | | |
| zS | 0 | | | |
| SS | 0. | | | |
| sigS | 0. | | | |
| lstoreplane | .false. | | | |
| iplane | | | | |
| linletRA | .false. | | | |
| lfixinlet | .false. | | | |
| lfixutauin | .false. | | | |
| lreadminl | .false. | | | |
| di | 0.09 | | | |
| dti | | | | |
