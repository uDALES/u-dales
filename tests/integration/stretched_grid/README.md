# Stretched Vertical Grid Test

This case preprocesses and runs the committed energy-balance case `064` on an
exponentially stretched vertical grid. No other committed case sets
`lzstretch`, so this is the end-to-end test of the stretched-grid code paths.

Status:

- classification: `supported`
- role: preprocessing and solver integration test
- expected in the supported dispatcher path
- shared fixture owner: `tests/cases/064/`

- Entry point: `test_stretched_grid.py`
- Grid: `lstretchexp` with `dzlin = 0.5` up to `hlin = 8`, stretched to 64 m
- Default MPI size: `4` (`nprocx = nprocy = 2`)
- Shared case fixture: `tests/cases/064/`

Prerequisites:

- the Python preprocessing environment with the IBM f2py module and View3D
  (`tools/build_preprocessing.sh <system> preprocessing_tools`)
- for the solver part: MPI launcher/runtime available in `PATH` and a built
  `u-dales` executable (`UDALES_BUILD`, default `build/release/u-dales`); that
  part is skipped when either is missing

The test copies the case into a temporary directory, switches the grid options
in the namelist, runs `tools/write_inputs.py` and checks that

- `prof.inp` describes the requested grid (linear near-wall part, monotone
  stretching, domain top at `zsize`);
- the facet sections on the scalar grid sum to the area of every facet;
- the facet sections on the `u` and `v` grids equal those of the committed
  uniform-grid preprocessing, facet by facet;
- the solver runs for a few time steps with the IBM wall functions, the energy
  balance and the periodic energy-balance correction active, stays divergence
  free to round-off and writes finite facet heat fluxes.

Run it directly with:

```bash
python tests/integration/stretched_grid/test_stretched_grid.py
```
