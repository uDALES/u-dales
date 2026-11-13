# Developing uDALES

uDALES is developed openly on [GitHub](https://github.com/uDALES/u-dales). Bug reports, feature requests, and pull requests are all welcome — see [how to contribute](udales-contributing.md) for the process.

Development starts from the `master` branch, which is the bleeding-edge version of the model; users running simulations should instead use the [latest release](https://github.com/uDALES/u-dales/releases), as described in the [installation guide](udales-installation.md).

This section collects the developer-facing documentation:

- [Code architecture](udales-architecture.md) — how the solver is organised: repository layout, module responsibilities, and the exact sequence of calls in a timestep. Start here before modifying the solver.
- [Extending uDALES](udales-extending.md) — recipes for common tasks, each traced through a working example in the code: adding a namelist option, an output variable, a statistic, an example case, or a test.
- [Testing and CI](udales-testing.md) — the test suites, how to run them locally, and what continuous integration checks on every pull request.
- [Development notes](udales-development-notes.md) — setting up a development environment, building in debug mode, and building this documentation.
- [How to contribute](udales-contributing.md) — reporting bugs, requesting features, and opening pull requests.
- [Parallelisation](udales-2decomp.md) — how the 2-D pencil domain decomposition with 2DECOMP&FFT works, and what it means for code that loops over the grid.
- [Using agents](udales-agents.md) — conventions for working on this repository with AI coding agents.
- [Fortran API reference](udales-software-docs.md) — API documentation of the Fortran source, generated with [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

The model source lives in `src/`, the pre/post-processing tools in `tools/`, example cases in `examples/`, and the tests in `tests/` (see the [test docs](https://github.com/uDALES/u-dales/blob/master/tests/README.md) for the layout). uDALES uses [semantic versioning](https://semver.org/).
