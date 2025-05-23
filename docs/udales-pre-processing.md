# uDALES preprocessing

In order to perform a simulation with uDALES, the input data needs to be processed, which creates a number of [input files](#Input-files) that will be read by the uDALES simulation.
This can be done from the command line using the shell script `write_inputs.sh`, which is a wrapper around the MATLAB script `write_inputs.m`. For more info about the functions see the [Developer's guide](#developers-guide). The script requires several variables in order to be execute without errors. Three main files required inside the experiment case directory are
1. an appropriately set namoptions.001 (assumming 001 is the case directory name) file, 
2. an STL file of the  building geometry (except for few special cases)
3. the config.sh file.

Below is an example setup for copying and pasting. You need to specify these parameters in a `config.sh` file within the experiment directory, which is then read by the scripts.

``` sh
# We assume you are running the following commands from your
# top-level project directory.

export DA_TOOLSDIR=$(pwd)/u-dales/tools # Directory of the scripts
export DA_EXPDIR=$(pwd)/experiments #  The top-level directory of the simulation setups
```

Before running the preprocessing, one must build the View3D submodule. This is a one time task and should be done as soon as you clone u-dales from GitHub.
``` sh
# We assume you are running the following commands from the u-dales directory.

# To build on local/common ubuntu or mac systems
./tools/build_preprocessing.sh common

# To build on ICL HPC
./tools/build_preprocessing.sh icl
```

Then, to start the pre-processing, run:

For local ubuntu or mac
``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: write_inputs.sh exp_id
./u-dales/tools/write_inputs.sh experiments/001
```

For ICL HPC
``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: write_inputs.sh esperiments/exp_id run_node_type

# To run preprocessing on HPC log in node (not recomended)
./u-dales/tools/write_inputs.sh experiments/001 l

# To run preprocessing on HPC compute node (recomended)
./u-dales/tools/write_inputs.sh experiments/001 c
```

In above example commands, replace 001 with the number of your simulation.

## Input files

uDALES requires a number of input files, all suffixed by the experiment number, which is omitted in the following documentation. The `namoptions.inp` file contains a list of parameters for uDALES and the pre-processing routines, and the pre-processing is intended to run using solely this file (with some exceptions). The input files are:

- `prof.inp`: initial profiles of flow variables (described in DALES documentation).
- `lscale.inp`: large-scale forcings (described in DALES documentation).
- `facets.inp`: column of facet types and three columns (x,y,z) of surface normal components (nfcts x 4).
- `Tfacinit.inp`: list of initial facet temperatures (nfcts x 1).
- `factypes.inp`: a description of the properties of different types of facet (described below). This is copied to the experiment directory automatically. Each row describes a wall type, with the first 6 columns being: wall type id, a boolean for whether it is a 'green' facet or not, momentum roughness length, heat roughness length, albedo, and emissivity. Assuming that each facet is composed of 3 layers, the next 3 columns give the thickness of each, the next 3 give the volumetric heat capacity, the next 3 give the heat conductivity, and the final 4 give the thermal diffusivity at each interface between the layers.
- `scalar.inp`: if using scalars, the initial scalar profiles.
- `timedepnudge.inp`: if using time-dependent nudging (described in DALES documentation).
- `timedepsurf.inp`: if using time-dependent surface fluxes (described in DALES documentation).

With the IBM, the following additional files are required:
- `solid_u.txt`: Indices of solid points on u-grid (nsolpts_u x 3).
- `solid_v.txt`: Indices of solid points on v-grid (nsolpts_v x 3).
- `solid_w.txt`: Indices of solid points on w-grid (nsolpts_w x 3).
- `solid_c.txt`: Indices of solid points on c-grid (nsolpts_c x 3).
- `fluid_boundary_u.txt`: Indices of fluid boundary points on u-grid (nbndpts_u x 3).
- `fluid_boundary_v.txt`: Indices of fluid boundary points on v-grid (nbndpts_v x 3).
- `fluid_boundary_w.txt`: Indices of fluid boundary points on w-grid (nbndpts_w x 3).
- `fluid_boundary_c.txt`: Indices of fluid boundary points on c-grid (nbndpts_c x 3).
- `facet_sections_u.txt`: Indices of facet sections on u-grid (nfctsecs_u x 3).
- `facet_sections_v.txt`: Indices of facet sections on v-grid (nfctsecs_v x 3).
- `facet_sections_w.txt`: Indices of facet sections on w-grid (nfctsecs_w x 3).
- `facet_sections_c.txt`: Indices of facet sections on c-grid (nfctsecs_c x 3).

With the SEB, the following additional files are required:
- `facetarea.inp`: list of areas of facets (nfcts x 1).
- `vf.nc.inp` or `vfsparse.inp`: view factors stored in netcdf or sparse format respectively.
- `svf.inp`: list of sky view factors of facets (nfcts x 1).
- `netsw.inp`: list of net shortwave radiation on facets (nfcts x 1).
If using time-dependent radiation on facets:
- `timedepsw.inp`: shortwave - row of times (1 x ntimedepsw) followed by array of values (nfcts x ntimedepsw).
- `timedeplw.inp`: longwave - columns for times and values (ntimedeplw x 2). This currently must be generated by a separate user-defined script.

## Developer's guide

The `u-dales/tools/preprocessing.m` matlab class contains the functionality for preprocessing. The constructor reads the parameters in `namoptions` and stores them as member variables, and defines default variables for those not specified. These are then used in the member functions. In these member functions, additional data structures are also stored as member variables, including those used repeatedly and those eventually written to files, so that one can easily view and manipulate them using the matlab IDE.

The `u-dales/tools/write_inputs.m` matlab script calls member functions of `preprocessing.m` in order to write the basic input files (those not relating to the IBM or SEB), followed by routines located in the `IBM` and `SEB` directories within the uDALES tools directory. It is intended to be as short and readable as possible, with the goal being that a developer can edit for a particular purpose. It will work simply as a normal script using the matlab IDE, but when doing this, ensure that `DA_EXPDIR = <top level directory>/experiments/` and `DA_TOOLSDIR = <top level directory>/u-dales/tools/` are defined.

The `u-dales/tools/write_inputs.sh` shell script acts as a wrapper around `write_inputs.m`. Before running the matlab script, it will run the shell script `config.sh` located in the experiment directory, which defines environmental variables `DA_EXPDIR` and `DA_TOOLSDIR`. After running the script, it will also write the correct number of facets to `namoptions`. It is intended to be run from the top level project directory.

Some parameters used by uDALES are used in the pre-processing. They are the following:

### `&RUN`
- `runtime`: Length of simulation period (in seconds). This is used with time-varying solar position (see below).

### `&DOMAIN`

- `itot`: number of cells in x-direction. Default: 64.
- `jtot`: number of cells in y-direction. Default: 64.
- `ktot`: number of cells in z-direction. Default: 96.
- `xlen`: domain size in x-direction (metres).
- `ylen`: domain size in y-direction (metres).

### `&WALLS`

This section describes the parameters used by the IBM.
- `iwallmom`: Momentum flux boundary condition - 1:zero flux, 2: non-neutral wall function, 3: neutral wall function.
- `iwalltemp`: Temperature flux boundary condition - 1: constant flux, 2: wall function.
If either of these are true, then `Tfacinit.inp.xxx` is written using the value of `facT`.

The following parameters are not used, but instead generated by the pre-processing routines and written to the file `info.txt` by the pre-processing script. They must be added to namoptions manually.
- `nfcts`: number of facets. If using `write_inputs.sh` (see below), this will write its value into namoptions. Equal to the number of (non-header) lines in `facets.inp.xxx`.
- `nsolpts_u`: number of solid points on the u-grid. Equal to the number of (non-header) lines in `solid_u.txt`.
- `nsolpts_v`: number of solid points on the v-grid. Equal to the number of (non-header) lines in `solid_v.txt`.
- `nsolpts_w`: number of solid points on the w-grid. Equal to the number of (non-header) lines in `solid_w.txt`.
- `nsolpts_c`: number of solid points on the c-grid. Equal to the number of (non-header) lines in `solid_c.txt`.
- `nsndpts_u`: number of fluid boundary points on the u-grid. Equal to the number of (non-header) lines in `fluid_boundary_u.txt`.
- `nsndpts_v`: number of fluid boundary points on the v-grid. Equal to the number of (non-header) lines in `fluid_boundary_v.txt`.
- `nsndpts_w`: number of fluid boundary points on the w-grid. Equal to the number of (non-header) lines in `fluid_boundary_w.txt`.
- `nsndpts_c`: number of fluid boundary points on the c-grid. Equal to the number of (non-header) lines in `fluid_boundary_c.txt`.
- `nfctsecs_u`: number of facet sections on the u-grid. Equal to the number of (non-header) lines in `facet_sections_u.txt`.
- `nfctsecs_v`: number of facet sections on the v-grid. Equal to the number of (non-header) lines in `facet_sections_v.txt`.
- `nfctsecs_w`: number of facet sections on the w-grid. Equal to the number of (non-header) lines in `facet_sections_w.txt`.
- `nfctsecs_c`: number of facet sections on the c-grid. Equal to the number of (non-header) lines in `facet_sections_c.txt`.

### `&PHYSICS`

- `luoutflowr`: switch that determines whether u-velocity is corrected to get a fixed outflow rate Default: false.
- `lvoutflowr`: switch that determines whether v-velocity is corrected to get a fixed outflow rate. Default: false.
- `luvolflowr`: switch that determines whether u-velocity is corrected to get a fixed volume flow rate. Default: false.
- `lvvolflowr`: switch that determines whether v-velocity is corrected to get a fixed volume flow rate. Default: false.
- `lcoriol`: switch for coriolis force. Default: false.
- `lprofforc`: switch for nudging flow to a profile. Default: false.

Note only one forcing should be specified, i.e. one of `luoutflowr`/`lvoutflowr`,`luvolflowr`/`lvvolflowr`, `lprofforc`, or `lcoriol`.

### `&ENERGYBALANCE`

- `lEB`: switch for energy balance. Default: false.
- `lvfsparse`: switch for view factors in sparse (text) format.
- `nnz`: number of non-zero view factors when using sparse format - this needs to be written to namoptions after pre-processing.
- `dtEB`: surface energy balance timestep.

### `&CHEMISTRY`

- `lchem`: switch for chemistry.

### `&SCALARS`

- `nsv`: number of scalar variables. Default: 0. Note that `nsv > 0` is not yet supported in the pre-processing.

### `&INPS`

The parameters under the `&INPS` header are used only in the pre-processing.

- `zsize`: size of domain in z direction (metres).
- `lzstretch`: switch for stretched z grid. Default: false.
- `lstretchexp`: switch for z grid stretched using exp function. Default: false.
- `lstretchtanh`: switch for z grid stretched using tanh function. Default: false.
- `lstretch2tanh`: switch for z grid stretched using 2tanh function. Default: false.
- `stretchconst`: stretch constant. Default: 0.01.
- `u0`: initial u-velocity (m/s). Also applied as geostrophic term where applicable. Default: 0.
- `v0`: initial v-velocity (m/s). Also applied as geostrophic term where applicable. Default: 0.
- `dpdx`: pressure gradient in x direction (Pa/m). Default: 0.
- `dpdy`: pressure gradient in y direction (Pa/m). Default: 0.
- `thl0`: temperature at z = 0. Default: 288.
- `qt0`: specific humidity at z = 0. Default: 0.
- `lapse`: lapse rate (K/m). Default: 0.
- `w_s`: subsidence. Default: 0.
- `R`: radiative forcing (W/m^2). Default: 0.
- `NOb`: initial concentration of NO. Default: 0.
- `NO2b`: initial concentration of NO2b. Default: 0.
- `O3b`: initial concentration of O3b. Default: 0.

The following parameters are related to the immersed boundary method.
- `stl_file`: Name of STL file defining the geometry.
- `read_types`: Switch for reading facet types from file. Default: false (all facets are set to type 1).
- `types_path`: Name of types file.
- `facT`: if `iwallmom = 2` or  `iwalltemp = 2` then (constant) facet temperature, or if `lEB = .true.` then initial facet temperature. Default: 288.
- `isolid_bound`: Option for classification of solid/fluid points, including boundary points. 1: In-house Fortran routine (default, fast), 2: equivalent MATLAB routine (useful for debugging), 3: inpolyhedron (MATLAB; provided for when option 1 & 2 are not producing expected results): https://www.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume.
- `ifacsec`: Option for facet section calculation. 1: Fortran (default, fast), 2: MATLAB (useful for debugging).

If using the energy balance, the following parameters can also be specified.
- `ishortwave`: Option for shortwave radiation calculation. 1: Fortran (default, fast), 2:  MATLAB (useful for debugging).
- `isolar`: Option for solar radiation (see below for futher detail). 1: custom (default), 2: from latitude & longitude, 3: from weather file.
- `view3d_out`: Output format for View3D: 0: text, 1: binary, 2: sparse (text). Default: 0.
- `maxD`: Maximum distance to check view factors, otherwise they are zero. Default: Inf.
- `xazimuth`: the azimuthal angle of the x-axis (with respect to North). Default : 90 degrees, i.e. East.

If `isolar = 1`, then the solar radiation is determined by:
- `solarazimuth`: solar azimuth (degrees). Default: 135. (solaz in uDALES v1).
- `solarzenith`: solar zenith (degrees). Default: 28.4066. (Z in uDALES v1).
- `I`: direct normal irradiance (DNI) (W/m^2). Default: 800.
- `Dsky`: diffuse sky irradiance (W/m^2). Default: 418.8041. (Dsk in uDALES v1).

If `isolar = 2`, then the solar position is calculated according to the NOAA solar position algorithm, and `I` and `Dsky` are determined by the ASHRAE clear-sky model, and the spacetime location is specified by:
- `year`, e.g. 2023.
- `month`, where e.g. 6 corresponds to June.
- `day`, e.g. 21.
- `hour`, where e.g. 0 corresponds to midnight and 23 corresponds to 11pm. Default: 6.
- `minute`. Default: 0.
- `second`. Default: 0.
- `longitude`. Default: -0.13.
- `latitude`. Default: 51.5.
- `timezone`. Default: 0.
- `elevation`. Default: 0.

If `isolar = 3`, then the parameters are specified in a file containing weather data for a given year.
- `weatherfname`: file name.
- `month`
- `day`
- `hour`. Default: 0.

The solar parameters can also be varied in time using the `ltimedepsw` switch. This occurs on a timescale `dtSP`, which is equal to `dtEB` by default.
