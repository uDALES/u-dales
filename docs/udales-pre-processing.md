# uDALES preprocessing

In order to perform a simulation with uDALES, the input data needs to be processed, which creates a number of [input files](#input-files) that will be read by the uDALES simulation.
This can be done from the command line using the shell script `write_inputs.sh`, which is a wrapper around either the MATLAB script `write_inputs.m` or the Python script `write_inputs.py`. For more info about the functions see the [Developer's guide](#developers-guide). The main files required inside the experiment case directory are

1. an appropriately set namoptions.001 (assumming 001 is the case directory name) file,
2. an STL file of the  building geometry (except for few special cases).

The wrapper also reads `config.sh` from the example directory when it is present. Use it for overrides such as:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

export DA_TOOLSDIR=$(pwd)/u-dales/tools # Optional; defaults to the script directory
export DA_EXPDIR=$(pwd)/experiments # Optional; defaults to the parent directory of the case

# Optional overrides for Imperial HPC compute-node preprocessing
export PREPROC_WALLTIME="24:00:00" # Defaults to 24:00:00
export PREPROC_MEM="128gb" # Defaults to 128gb
```

`tools/write_inputs.sh` scans `namoptions.###` for `nompthreads` and uses it for the preprocessing CPU request. If `nompthreads` is omitted, the preprocessing default is `8`. If `nompthreads` appears more than once anywhere in the file, the wrapper stops and asks for a single value. The wrapper exports the derived value internally as `PREPROC_NCPU` for PBS submission and for the default View3D thread count, unless `VIEW3D_NUM_THREADS` is explicitly overridden.

When running preprocessing on an Imperial HPC compute node with `write_inputs.sh ... c`, the wrapper uses `PREPROC_WALLTIME="24:00:00"` and `PREPROC_MEM="128gb"` unless these are set in `config.sh` or the calling environment. These size the preprocessing PBS job only; they are separate from the solver job variables `WALLTIME` and `MEM` used by `hpc_execute.sh`. `PREPROC_MEM` must be written as a number followed by lowercase `gb`, such as `128gb`; a unitless value such as `128` is rejected before submitting the PBS job. Unless `VIEW3D_MAX_DENSE_MATRIX_GIB` is explicitly set, the default View3D configuration derives the dense-matrix guard from `PREPROC_MEM`: requests above `16gb` leave 16 GiB for overhead, while smaller requests use the requested GiB value. For example, `PREPROC_MEM="128gb"` gives `VIEW3D_MAX_DENSE_MATRIX_GIB=112`.

Before running the Python preprocessing route, build the Python virtual environment. The setup script creates the environment, installs all dependencies, and builds the preprocessing tools (View3D and f2py extension modules). For the MATLAB route, the same setup is still a convenient way to build View3D and the bundled Python helpers used by some preprocessing paths, such as vegetation conversion. For more details on virtual environment setup see [here](./../tools/python/README_VENV.md):

```bash
# For a local machine
bash tools/python/setup_venv.sh common

# For the Imperial HPC machine
bash tools/python/setup_venv.sh icl
```
This is a one-time setup task and should be repeated after pulling changes that affect Python dependencies or preprocessing binaries.

<!---
``` sh
# We assume you are running the following commands from the u-dales directory.

# To build on local/common ubuntu or mac systems
./tools/build_preprocessing.sh common preprocessing_tools

# To build on ICL HPC
./tools/build_preprocessing.sh icl preprocessing_tools
```
--->

Then, to start the pre-processing, run:

For local ubuntu or mac

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: write_inputs.sh <-m|-p> <PATH_TO_CASE> [start]
./u-dales/tools/write_inputs.sh -m experiments/001

# Or run the Python preprocessing route
./u-dales/tools/write_inputs.sh -p experiments/001
```

For ICL HPC

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: write_inputs.sh <-m|-p> experiments/exp_id run_node_type

# To run preprocessing on HPC log in node (not recomended)
./u-dales/tools/write_inputs.sh -m experiments/001 l

# To run preprocessing on HPC compute node (recomended)
# For MATLAB
./u-dales/tools/write_inputs.sh -m experiments/001 c
# For Python
./u-dales/tools/write_inputs.sh -p experiments/001 c
```

In above example commands, replace 001 with the number of your example.

Note: After editing or creating a `namoptions.###` file (use an existing template where possible), run `write_inputs.sh` to generate or update geometry-related input files such as `facets.inp` and `factypes.inp`. If your setup requires generating geometry from STL files during pre-processing, set the `gen_geom` parameter in the `namoptions` file so the pre-processor will build geometry from the provided STL.

## Input files

uDALES requires a number of input files, all suffixed by the experiment number, which is omitted in the following documentation. The `namoptions.###` file contains a list of parameters for uDALES and the pre-processing routines, and the pre-processing is intended to run using solely this file (with some exceptions). The input files are:

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

The `u-dales/tools/write_inputs.sh` shell script acts as a wrapper around either `write_inputs.m` or `write_inputs.py`. Before running the selected route, it will run the shell script `config.sh` located in the experiment directory when that file is present. It defaults `DA_TOOLSDIR` to the directory containing `write_inputs.sh` and `DA_EXPDIR` to the parent directory of the case directory. It scans `namoptions.###` for `nompthreads`, requires at most one occurrence, and exports the derived preprocessing CPU count internally as `PREPROC_NCPU`. The first argument must be either `-m` for MATLAB or `-p` for Python. It is intended to be run from the top level project directory.

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

The following parameters are generated by the pre-processing routines. The current MATLAB and Python preprocessing routes write the generated values back to `namoptions.###` when they are generated, so check the updated namelist after preprocessing completes.

- `nfcts`: number of facets. If using `write_inputs.sh` (see below), this will write its value into namoptions. Equal to the number of (non-header) lines in `facets.inp.xxx`.
- `nsolpts_u`: number of solid points on the u-grid. Equal to the number of (non-header) lines in `solid_u.txt`.
- `nsolpts_v`: number of solid points on the v-grid. Equal to the number of (non-header) lines in `solid_v.txt`.
- `nsolpts_w`: number of solid points on the w-grid. Equal to the number of (non-header) lines in `solid_w.txt`.
- `nsolpts_c`: number of solid points on the c-grid. Equal to the number of (non-header) lines in `solid_c.txt`.
- `nbndpts_u`: number of fluid boundary points on the u-grid. Equal to the number of (non-header) lines in `fluid_boundary_u.txt`.
- `nbndpts_v`: number of fluid boundary points on the v-grid. Equal to the number of (non-header) lines in `fluid_boundary_v.txt`.
- `nbndpts_w`: number of fluid boundary points on the w-grid. Equal to the number of (non-header) lines in `fluid_boundary_w.txt`.
- `nbndpts_c`: number of fluid boundary points on the c-grid. Equal to the number of (non-header) lines in `fluid_boundary_c.txt`.
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
- `nnz`: number of non-zero view factors when using sparse format. This is written to `namoptions.###` by preprocessing when sparse view factors are generated.
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
- `nompthreads`: number of OpenMP threads used by preprocessing routines that support OpenMP, including IBM preprocessing. Default: 8. `write_inputs.sh` also uses this value for the preprocessing CPU request and requires that it appears at most once in `namoptions.###`.
- `u0`: initial u-velocity (m/s). Also applied as geostrophic term where applicable. Default: 0.
- `v0`: initial v-velocity (m/s). Also applied as geostrophic term where applicable. Default: 0.
- `dpdx`: pressure gradient in x direction (Pa/m). Default: 0.
- `dpdy`: pressure gradient in y direction (Pa/m). Default: 0.
- `thl0`: temperature at z = 0. Default: 288.
- `qt0`: specific humidity at z = 0. Default: 0.
- `lapse`: lapse rate (K/m). Default: 0.
- `w_s`: subsidence. Default: 0.
- `dqtdxls`: large-scale advective tendency of specific humidity in x direction (kg/kg/s). Default: 0.
- `dqtdyls`: large-scale advective tendency of specific humidity in y direction (kg/kg/s). Default: 0.
- `dqtdtls`: large-scale advective tendency of specific humidity in time (kg/kg/s). Default: 0.
- `tke`: initial turbulence kinetic energy (m^2/s^2). Default: 0.
- `R`: radiative forcing (W/m^2). Default: 0.
- `NOb`: initial concentration of NO. Default: 0.
- `NO2b`: initial concentration of NO2b. Default: 0.
- `O3b`: initial concentration of O3b. Default: 0.

The following parameters are related to the immersed boundary method.

- `stl_file`: Name of STL file defining the geometry.
- `read_types`: Switch for reading facet types from file. Default: false (all facets are set to type 1).
- `types_path`: Name of types file.
- `facT`: if `iwallmom = 2` or  `iwalltemp = 2` then (constant) facet temperature, or if `lEB = .true.` then initial facet temperature. Default: 288.
- `isolid_bound`: Option for classification of solid/fluid points, including boundary points. 1: In-house Fortran routine (default, fast), 2: equivalent MATLAB routine (useful for debugging), 3: inpolyhedron (MATLAB; provided for when option 1 & 2 are not producing expected results): <https://www.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume>.
- `ifacsec`: Option for facet section calculation. 1: Fortran (default, fast), 2: MATLAB (useful for debugging).

If using the energy balance, the following parameters can also be specified.

- `ishortwave`: Option for shortwave radiation calculation. The MATLAB preprocessing default is `1`; the Python preprocessing default is `3`.
- `isolar`: Option for solar radiation (see below for futher detail). 1: custom (default), 2: from latitude & longitude, 3: from weather file.
- `view3d_out`: Output format for View3D: 0: text, 1: binary, 2: sparse (text). The MATLAB preprocessing default is `0`; the Python preprocessing default is `2`.
- `maxD`: Maximum distance to check view factors, otherwise they are zero. Default: Inf.
- `xazimuth`: the azimuthal angle of the x-axis (with respect to North). Default : 90 degrees, i.e. East.

`ishortwave` values are interpreted as follows during preprocessing:

| Value | MATLAB preprocessing | Python preprocessing |
| ----- | -------------------- | -------------------- |
| `1` | Standalone Fortran scanline rasterization (no vegetation). | f2py scanline rasterization wrapper (same scanline algorithm, no vegetation). |
| `2` | MATLAB scanline debug implementation (slow; no vegetation). | Unsupported; raises an error because this implementation is MATLAB-only. |
| `3` | Unsupported; raises an error because this is Python-only. | `facsec` method, using ray casting with solid mask and facet-section reconstruction. |
| `4` | Unsupported; raises an error because this is Python-only. | `moller` method, using Moller-Trumbore triangle hits. |

The explicit Python API backend `method="scanline_legacy"` remains available as a reference/parity implementation, but it is not a normal documented namelist choice.

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
