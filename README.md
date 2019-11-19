# uDALES [![Build Status](https://travis-ci.com/uDALES/u-dales.svg?token=3tqUbxqJuLtozjxqDymC&branch=master)](https://travis-ci.com/uDALES/u-dales)

This is the development repository for the uDALES (urban Dutch Atmospheric Large Eddy Simulation) model. It contains the uDALES software source-code, tests, documentation, and examples with tools to pre- and post-process inputs/outputs.

**If you are new to uDALES, please follow our [getting started guide](https://udales.github.io/u-dales/0YiO263pFxExSdkMvWfId3qkVUSF4dREFnwM1jQD9y1KvzeAVAWzGykQemUrkJCM/html/udales-getting-started).**

## Installation

If you are new to uDALES or need to set up your own experiments, please see our [getting started guide](https://udales.github.io/u-dales/0YiO263pFxExSdkMvWfId3qkVUSF4dREFnwM1jQD9y1KvzeAVAWzGykQemUrkJCM/html/udales-getting-started). If you are looking to develop the code instead, you can install uDALES on Linux, macOS, and WSL with the following commands:

``` sh
mkdir -p build/release
pushd build/release
cmake ../..
make
```

To know more about build options, see [build/default options](https://udales.github.io/u-dales/0YiO263pFxExSdkMvWfId3qkVUSF4dREFnwM1jQD9y1KvzeAVAWzGykQemUrkJCM/html/udales-getting-started/#build-defaultsoptions).


## Usage

A uDALES simulation needs to be executed from a directory containing all required input files. Examples of experiments and required inputs are in the `examples` directory. To run a uDALES simulation you need to specify the number of cpus `<NCPU>`, the path to the build file `<BUILD>` and the simulation configuration file `<NAMOPTIONS>` and execute the simulation with the following command:

``` sh
mpiexec -n <NCPU> <BUILD> <NAMOPTIONS>
```


## Documentation

Currently work in progress: see [this link](https://uDALES.github.io/u-dales/0YiO263pFxExSdkMvWfId3qkVUSF4dREFnwM1jQD9y1KvzeAVAWzGykQemUrkJCM/html/index.html).


## Testing

TODO: add info about current testing and CI.


## Contributing

TODO:

## Versioning

TODO: need create first release as soon as cleanup is complete.

## Copyright and License

FIXME: need authorship and year

This project is licensed under the GNU General Public License - see the [LICENSE.md](LICENSE.md) file for details