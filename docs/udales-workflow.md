# Workflow

The practical workflow is:

1. Build the model so that `u-dales/bin/u-dales` is available.
2. Set `UD_TOPDIR` to the repository root and prepend `$UD_TOPDIR/bin` to `PATH`.
3. Prepare a case under `examples/` or your own experiment directory.
4. Run the case directly with `u-dales namoptions.xxx` or through the wrapper scripts documented in [Running uDALES](./udales-simulation-setup.md).
5. Merge and inspect outputs using the post-processing tools.
6. Run the relevant checks under [`tests/README.md`](../tests/README.md) before merging changes.

For the input interface specifically:

- runtime input is namelist-based
- schema files under `docs/schemas/` are for tooling and editor support
- the input round-trip and source/schema checks live under `tests/integration/input`
