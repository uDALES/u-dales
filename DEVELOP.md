# Development notes

## Set up

Install all required packages for uDALES described in the [prerequisites section](https://udales.github.io/u-dales/udales-installation/#prerequisites), plus optionally [Graphviz](https://graphviz.org/) for generating graphs in the code viewer. E.g. installing all the required packages using Ubuntu's APT:

```sh
sudo apt update && sudo apt install -y gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev graphviz
```

Then, to set up the development environment for testing and generating the docs, download the [latest version of Miniconda](https://docs.conda.io/en/latest/miniconda.html) and install the required dependencies with:

```sh
conda env create -f environment.yml
```

Then activate with `conda activate udales`.

For the Python tooling (pre-processing and the Python package tests), create the project virtual environment with:

```sh
./tools/python/setup_venv.sh common
```

This creates the environment at `tools/python/.venv`, installs all dependencies, and builds the compiled preprocessing extensions (View3D and the f2py modules). See [tools/python/README_VENV.md](https://github.com/uDALES/u-dales/blob/master/tools/python/README_VENV.md) for details.

## Building for development

Building the model is described in the [installation guide](https://udales.github.io/u-dales/udales-installation/); the same instructions apply for development. In addition, developers will usually want a `Debug` build, which enables runtime checks and floating-point exception trapping:

```sh
mkdir -p build/debug
pushd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug ../..
make
```

For quick testing without the wrapper scripts, the solver can be run directly from a directory containing all input files:

``` sh
mpiexec -n <NCPU> <BUILD> <NAMOPTIONS>
```

## Testing

Tests are dispatched from a manifest with `tests/run_tests.py` (suites are defined in `tests/test_suites.yml`). Please refer to the [test docs](https://github.com/uDALES/u-dales/blob/master/tests/README.md) for the test layout and execution contracts.

## Documentation

The user documentation is built with [MkDocs](https://www.mkdocs.org/) (Material theme) and the software docs with [FORD](https://github.com/Fortran-FOSS-Programmers/ford); both are installed by the conda environment above. To build:

```sh
mkdocs build --site-dir build/html
ford docs/udales-docs-software.md
```

For live preview while editing, use `mkdocs serve`.

### Examples input plots

To create domain plots of the examples, run the following from your command line (requires MATLAB):

```sh
 matlab -nosplash -nodesktop -r "cd('tools/examples'); plot_blocks('<CASE_NUMBER>'); quit"
```

where `<CASE_NUMBER>` is e.g. `201`. Plots are then saved in their respective example folders.

### Examples outputs and plots

Run the following script to run and generate outputs for all example cases:

```sh
./tools/examples/run_examples.sh
```

Then, to create a sample plot for case `102` run the following from your command line (requires MATLAB):

```sh
matlab -nosplash -nodesktop -r "cd('tools/examples'); plot_fielddump_slice('102','u','y',32,1); quit"
```

## Versioning

This project uses [semantic versioning](https://semver.org/).
