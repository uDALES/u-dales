# Development notes

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


## Versioning

This project uses [semantic versioning](https://semver.org/).
