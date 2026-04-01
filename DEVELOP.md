# Development notes

## Set up

Install all required packages for uDALES described in the [prerequisites section](https://udales.github.io/u-dales/udales-getting-started/#prerequisites-when-not-using-singularity), plus optionally [Graphviz](https://graphviz.org/) for generating graphs in the code viewer. E.g. installing all the required packages using Ubuntu's APT:

```sh
sudo apt update && sudo apt install -y gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev graphviz
```

Then, to set up the development environment for testing and generating the docs, download the [latest version of Miniconda](https://docs.conda.io/en/latest/miniconda.html) and install the required dependencies with:

```sh
conda env create -f environment.yml
```

Then activate with `conda activate udales`.

If a shared virtual environment already exists at `../.venv/`, you can
also use that directly for Python tooling and tests, for example:

```sh
../.venv/bin/python -m unittest tools/python/tests/test_namelist.py
```

The Python radiation tooling also requires the compiled View3D executable and
the `directshortwave_f2py` wrapper. The simplest setup path is:

```sh
PYTHON_BIN=/opt/pbs/python/bin/python3 ./tools/python/setup_venv.sh
```

Use another interpreter only if it provides the matching Python development
headers needed by f2py.

New Linux/WSL setups now default to `../.venv/`, matching the Windows
setup script. If you already have a repo-local `.venv/`, the setup script will
keep using it until you migrate manually.

## Installation

To install uDALES on Linux, macOS, and WSL, use the following commands from the command prompt:

```sh
mkdir -p build/release
pushd build/release
cmake ../..
make
```

To know more about build options, please see [build/default options](https://udales.github.io/u-dales/udales-getting-started/#build-defaultsoptions).

## Running

A uDALES simulation needs to be executed from a directory containing all required input files. Examples of experiments and required inputs are in the `examples` directory. To run a uDALES simulation you need to specify the number of cpus `<NCPU>`, the path to the build file `<BUILD>` and the simulation configuration file `<NAMOPTIONS>` and execute the simulation with the following command:

``` sh
mpiexec -n <NCPU> <BUILD> <NAMOPTIONS>
```

## Testing

Please refer to [Test docs](https://github.com/uDALES/u-dales/blob/master/tests/README.md).

## Documentation

```sh
mkdocs build --site-dir build/html
ford docs/udales-docs-software.md
```

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
