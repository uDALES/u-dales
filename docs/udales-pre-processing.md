# Pre-processing

## Input files

uDALES requires a number of input files, all suffixed by the experiment number, which is taken to be `009` for illustratitive purposes. The `namoptions.inp.009` file contains a list of parameters for uDALES and the pre-processing routines, and the pre-processing is intended to run using solely this file (with some exceptions). The input files are:

- `xgrid.inp.009`: x grid values.
- `zgrid.inp.009`: z grid values.
- `prof.inp.009`: initial profiles of e.g. u,v,w,T.
- `lscale.inp.009`: large-scale forcings.
- `blocks.inp.009`: block position and corresponding facets.
- `facets.inp.009`: facet orientation and corresponding block.
- `Tfacinit.inp.009`: initial facet temperatures.
- `walltypes.inp.009`: a description of the properties of different types of facet (this is copied to the experiment directory automatically).

If using scalars, `scalar.inp.009` is also required to specifiy the initial scalar profiles.
When using the energy balance, the following additional files are required:

- `facetarea.inp.009`: areas of facets.
- `vf.nc.inp.009`: view factors of facets.
- `scf.inp.009`: sky view factor of facets.
- `netsw.inp.009`: net shortwave radiation on facets.

## Pre-processing parameters
Some parameters used by uDALES are also used in the pre-processing. They are the following:

### `&RUN`

- `ltrees`: switch for trees. Default: false. MAYBE SHOULDN'T BE HERE FOR NOW
- `lpurif`: switch for purifiers. Default: false. MAYBE SHOULDN'T BE HERE FOR NOW
- `nsv`: number of scalar variables. Default: 0. MAYBE SHOULDN'T BE HERE FOR NOW
- `luoutflowr`: switch that determines whether u-velocity is corrected to get a fixed outflow rate Default: false.
- `lvoutflowr`: switch that determines whether v-velocity is corrected to get a fixed outflow rate. Default: false.
- `luvolflowr`: switch that determines whether u-velocity is corrected to get a fixed volume flow rate. Default: false.
- `lvvolflowr`: switch that determines whether v-velocity is corrected to get a fixed volume flow rate. Default: false.

Note `luoutflowr` should not be used with `luvolflowr`, and similarly with `lvoutflowr` and `lvvolflowr`.


### `&DOMAIN`

- `imax`: number of cells in x-direction. Default: 64.
- `jtot`: number of cells in y-direction. Default: 64.
- `kmax`: number of cells in z-direction. Default: 96.
- `xsize`: domain size in x direction (metres).
- `ysize`: domain size in y direction (metres).
- `nblocks`: MAYBE SHOULDN'T BE HERE
- `ntrees`: number of trees. Default: 0. MAYBE SHOULDN'T BE HERE FOR NOW
- `npurif`: number of purifiers. Default: 0. MAYBE SHOULDN'T BE HERE FOR NOW

### `&ENERGYBALANCE`

- `lEB`: switch for energy balance. Default: false.
- `nfcts`: MAYBE SHOULDN'T BE HERE

### `&PHYSICS`

- `lcoriol`: switch for coriolis force. Default: false.
- `lprofforc`: switch for nudging flow to a profile. Default: false.
- `lchem`: switch for chemistry.
 
### `&INPS`
The parameters under the `&INPS` header are used only in the pre-processing.

- `zsize`: size of domain in z direction (metres).
- `lzstretch`: switch for stretched z grid. Default: false.
- `stretchconst`: stretch constant. Default: 0.01.
- `u0`: initial u-velocity (m/s). Also applied as geostrophic term where applicable. Default: 0.
- `v0`: initial v-velocity (m/s). Also applied as geostrophic term where applicable. Default: 0.
- `dpdx`: pressure gradient in x direction (Pa/m). Default: 0.
- `dpdy`: pressure gradient in y direction (Pa/m). Default: 0.
- `thl0`: temperature (at z=0?). Default: 288. could be `thls`?
- `qt0`: specific humidity (at z=0?). Default: 0. could be `qts`?
- `lapse`: lapse rate (K/m). Default: 0.
- `w_s`: subsidence (?units?). Default: 0.
- `R`: radiative forcing (?units?). Default: 0.

If using the energy balance, the following parameters can also be specified.

- `solaz`: solar azimuth (degrees). Default: 135.
- `Z`: solar zenith (degrees). Default: 28.4066.
- `I`: direct solar radiation (W/m^2). Default: 184.8775.
- `Dsk`: diffuse incoming radiation (W/m^2). Default: 418.8041.

The following parameters relate to generating `blocks.inp`. Only one of the following three methods should be used.

#### Generate blocks from a given file

- `lblocksfile`: switch for generating blocks from a given file.
- `blocksfile`: name of blocks file (must be specified if `lblocksfile` is used).

#### Generate blocks from LIDAR data

- `llidar`: switch for generating blocks from LIDAR data. Default: false.
- `sourcename`: name of image file.
- `dxinp`: resolution of image in x direction (metres per pixel). Default: 1.
- `dyinp`: resolution of image in y direction (metres per pixel). Default: 1.
- `centeri`: center of area of interest in image (horizontal pixel). Default: 0.
- `centerj`: center of area of interest in image (vertical pixel). Default: 0.
- `maxh`: maximum height of buildings in image (metres). Default: 0.
- `pad`: padding (metres). Default: 5.
- `smallarea`: objects smaller than this area (metres) will be deleted. Default: round(150 / (dx * dy)), where dx = xsize / imax and dy = ysize / jtot.

#### Generate simple block geometry

- `lflat`: switch for no blocks. Default: false.
- `lcube`: switch for linear cubes. Default: false.
- `lcastro`: switch for staggered cubes. Default: false.
- `lcanyons`: switch for inifite canyons. Default: false.
- `blockheight`: height of blocks (metres?). Required when using `lcube`, `lcastro`, or `lcanyons`. Default: 16.
- `blockwidth`: width of blocks (metres?). Required when using `lcube`, `lcastro`, or `lcanyons`. Default: 16.
- `canyonwidth`: width of canyon (metres?). Reuired when using `lcanyons`. Default: 16.


## Instructions

The `da_inp` bash script in `u-dales/tools/pre` acts as a wrapper around the matlab pre-processing routines.

```sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: da_inp.sh exp_id
./u-dales/tools/pre/da_inp 009
```

This will write the necessary input files according to the parameters in `namoptions.inp.009`, as well as updating the number of blocks and facets. 
