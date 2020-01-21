# Pre-processing

## Input files

uDALES requires a number of input files, all suffixed by the experiment number, which is taken to be `009` for illustratitive purposes. The `namoptions.inp.009` file contains a list of parameters for uDALES and the the pre-processing routines, and the pre-processing routine is intended to run using solely this file. The input files are:

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
- `netsw.inp.009`: net shortwave on facets.

## Pre-processing parameters
Some parameters used by uDALES are also used in the pre-processing. They are the following:

### `&RUN`

- `ltrees`: switch for trees.
- `lpurif`: switch for purifiers.
- `luoutflowr`: switch for constant outflow rate in x direction.
- `lvoutflor`: switch for constant outflow rate in y direction.
- `luvolflowr`: switch for constant volume flow rate in x direction.
- `lvvolflor`: switch for constant volume flow rate in y direction.
- `ldp`: switch for constant pressure gradient.

### Not sure
`thl0` could be `thls` from `&BC`?
`qt0` could be `qts` from `&BC`?


### `&DOMAIN`

- `imax`
- `jtot`
- `kmax`
- `xsize`
- `ysize`
- `nblocks` MAYBE SHOULDN'T BE HERE
- `ntrees`
- `npurif`

### `&BC`

### `&ENERGYBALANCE`

- `lEB`
- `nfcts` MAYBE SHOULDN'T BE HERE

### `&PHYSICS`
- `lcoriol`
- `lprofforc`
 
Some parameters are used only in the pre-processing.
### `&INPS`

## Instructions

The `da_inp` bash script in `u-dales/tools/pre` acts as a wrapper around the matlab pre-processing routines.

```sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: da_inp.sh exp_id
./u-dales/tools/pre/da_inp 009
```

This will write the necessary input files according to the parameters in `namoptions.inp.009`, as well as updating the number of blocks and facets. 
