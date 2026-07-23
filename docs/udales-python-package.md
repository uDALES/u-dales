# The uDALES Python package

!!! note "Status: under active testing"
    The uDALES Python package is under active testing and its interfaces may
    change. It is intended to replace the MATLAB toolchain, which is
    deprecated and will be retired from uDALES v3.0 onwards.

The package lives in `tools/python` and provides a set of modules for
pre-processing uDALES cases and for loading, analysing, and visualising
simulation output, as a Python-native alternative to the MATLAB tools.

## What it provides

| Module | Purpose |
| --- | --- |
| `udprep` | Pre-processing: builds a case's input files (grid, IBM, vegetation, radiation, SEB) from `namoptions` and STL geometry via the `UDPrep` class, driving the compiled f2py/Fortran extensions and View3D. |
| `udbase` | Loading and analysing simulation output: the `UDBase` class reads a case's NetCDF fields, statistics, and facet data, mirroring the MATLAB `udbase` class. |
| `udgeom` | Geometry handling: the `UDGeom` package for STL mesh checking/repair, building separation, outline calculation, and View3D view-factor computation. |
| `udvis` | 3-D visualisation: `UDVis` gives a stable plotting entry point (`sim.vis`) on top of `UDBase`, keeping the rendering backend an internal detail. |
| `udstats` | Statistics: stateless array-in/array-out helpers (e.g. time averaging, variance/covariance) used internally by `UDBase`. |
| Supporting modules | `udconfig`, `udgrid`, `udnetcdf`, `udfacet`, and `exceptions` provide namelist parsing, grid-coordinate maths, NetCDF loading, and facet-data helpers used internally by `udbase` and `udprep`. |

## Installation

Set up the Python virtual environment from the repository root using
`tools/python/setup_venv.sh`, passing `common` for a local machine or `icl`
for the Imperial College HPC cluster. This creates the virtual environment
at `tools/python/.venv`, installs the runtime and build dependencies, and
builds the compiled preprocessing extensions (View3D and the
`directshortwave`/`ibm_preproc` f2py modules):

```bash
# For a local machine
bash tools/python/setup_venv.sh common

# For the Imperial HPC machine
bash tools/python/setup_venv.sh icl
```

See [README_VENV.md](https://github.com/uDALES/u-dales/blob/master/tools/python/README_VENV.md) for the full reference,
including the Windows PowerShell setup script, manual rebuilds, and
registering the environment as a Jupyter kernel.

## Using it for pre-processing

The `write_inputs.sh` wrapper script can drive either the MATLAB or the
Python preprocessing route; pass `-p` to use Python (`-m` selects MATLAB).
The Python route requires the virtual environment above to be set up first:

```sh
# We assume you are running the following commands from your
# top-level project directory.

./u-dales/tools/write_inputs.sh -p experiments/001
```

This invokes `tools/write_inputs.py`, which drives the `UDPrep` stack to
generate the same set of input files as the MATLAB route (grid, IBM,
vegetation, and, if the energy balance is enabled, radiation/SEB files). See
[Pre-processing](udales-pre-processing.md) for the full set of input files
and namelist options.

## Tutorials

The documented tutorials on this site remain MATLAB-based for now. Python
tutorial notebooks and example scripts live in `tools/python/examples/` in
the repository:

- `udbase_tutorial.ipynb` — loading and analysing simulation output with `UDBase`.
- `fields_tutorial.ipynb` — working with flow field output.
- `facets_tutorial.ipynb` — working with facet data.
- `geometry_tutorial.ipynb` — building and manipulating geometry with `udgeom`.
- `geometry_QA_tutorial.ipynb` — geometry quality-assurance checks.
- `plot_prep_tutorial.ipynb` — visualising preprocessing output.
- `udprep_tutorial.py` — driving the `UDPrep` pre-processing workflow.
- `udprep_radiation_tutorial.py` — radiation/shortwave pre-processing.
- `udprep_vegetation_tutorial.py` — vegetation pre-processing.

## Feedback

The Python package is new and still evolving — please report issues or
unexpected behaviour on the [uDALES GitHub issue tracker](https://github.com/uDALES/u-dales/issues).
