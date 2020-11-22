# Development notes

## Set up

Install all packages described in the [prerequisites section](./docs/udales-getting-started.md#prerequisites) plus [Graphviz](https://graphviz.org/) for generating graphs in the code viewer. E.g. installing all the required packages using Ubuntu's APT:

```sh
sudo apt update && sudo apt install -y gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev graphviz
```

Then, to set up the development environment for testing and doc generation, download the [latest version of Miniconda](https://docs.conda.io/en/latest/miniconda.html) and installed the required dependency with:

```
conda env create -f environment.yml
```

Then activate with `conda activate udales`.


## Installation

To install uDALES on Linux, macOS, and WSL, use the following commands from the command prompt:

```sh
mkdir -p build/release
pushd build/release
cmake ../..
make
```

To know more about build options, please see [build/default options](https://udales.github.io/u-dales/0YiO263pFxExSdkMvWfId3qkVUSF4dREFnwM1jQD9y1KvzeAVAWzGykQemUrkJCM/html/udales-getting-started/#build-defaultsoptions).


## Running

A uDALES simulation needs to be executed from a directory containing all required input files. Examples of experiments and required inputs are in the `examples` directory. To run a uDALES simulation you need to specify the number of cpus `<NCPU>`, the path to the build file `<BUILD>` and the simulation configuration file `<NAMOPTIONS>` and execute the simulation with the following command:

``` sh
mpiexec -n <NCPU> <BUILD> <NAMOPTIONS>
```


## Testing

Please refer to [Test docs](tests/README.md).


## Documentation

```
mkdocs build --site-dir build/html
ford docs/udales-docs-software.md
```

## Versioning

This project uses [semantic versioning](https://semver.org/).