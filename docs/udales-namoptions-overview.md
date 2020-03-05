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
| nblocks | 0 | | Number of blocks specified in `blocks.inp`. | - |
| **ntrees** | 0 | | Number of trees specified in `trees.inp`. (*Not supported in current version.*) | - |
| **npurif** | 0 | | Number of air purifiers specified in `purifs.inp`. (*Not supported in current version.*) | - |


# Namelist DYNAMICS

Possible advection schemes: 1 = 1st order upwind scheme, 2 = 2nd order central difference scheme, 7 = Kappa (flux limited) scheme. This scheme designed for quantities that should never become negative.

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| lqlnr | .false. | | See [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). Default switched from .true. | |
| iadv_mom | 2 | 2 | Advection scheme for momentum. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| iadv_tke | -1 | 2 | Advection scheme for TKE. Only used if `loneeqn = True`. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| iadv_thl | -1 | 2, 7 | Advection scheme for temperature. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| iadv_qt | -1 | 2 | Advection scheme for moisture. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| iadv_sv | -1 | 1, 2, 7 | Advection scheme for scalars. Also in [DALES](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf). | - |
| **ipoiss** | 1 | 0, 1 | Poisson solver. 0 = Fast Fourier Transformation, 1 = Cyclic reduction scheme. | - |


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
| **numol** | | | | |
| **prandtlmol** | | | | |

# Namelist RUN

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| | | | | |


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
| | | | | |


# Namelist WALLS

| Name | Default | Possible values | Description | Unit |
| ---- | ------- | --------------- | ----------- | ---- |
| | | | | |


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
