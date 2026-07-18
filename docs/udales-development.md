# Developing uDALES

uDALES is developed openly on [GitHub](https://github.com/uDALES/u-dales). Bug reports, feature requests, and pull requests are all welcome — see [how to contribute](udales-contributing.md) for the process.

This section collects the developer-facing documentation:

- [Development notes](udales-development-notes.md) — setting up a development environment, building in debug mode, running the tests, and building this documentation.
- [How to contribute](udales-contributing.md) — reporting bugs, requesting features, and opening pull requests.
- [Parallelisation](udales-2decomp.md) — how the 2-D pencil domain decomposition with 2DECOMP&FFT works, and what it means for code that loops over the grid.
- [Using agents](udales-agents.md) — conventions for working on this repository with AI coding agents.
- [Fortran API reference](udales-software-docs.md) — API documentation of the Fortran source, generated with [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

The model source lives in `src/`, the pre/post-processing tools in `tools/`, example cases in `examples/`, and the tests in `tests/` (see the [test docs](https://github.com/uDALES/u-dales/blob/master/tests/README.md) for the layout). uDALES uses [semantic versioning](https://semver.org/).
