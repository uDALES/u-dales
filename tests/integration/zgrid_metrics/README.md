# Vertical Grid Metrics Test

This case validates the vertical grid metrics the solver derives from the cell
centres in `prof.inp`, on a uniform and on a stretched grid.

Status:

- classification: `supported`
- role: solver-facing integration test
- expected in the supported dispatcher path
- shared fixture owner: `tests/cases/100/`

- Entry point: `run_test.sh`
- Run mode: `1006`
- MPI sizes: `1` (uniform, stretched, expected failure) and `4` (stretched, 2x2)
- Shared case fixture: `tests/cases/100/`, with `prof_stretched.inp.100`
  replacing the profile for the stretched modes (dzlin = 0.5 m up to
  hlin = 32 m, exponential stretching to 100 m, spacing ratio 3.8)

Prerequisites:

- MPI launcher/runtime available in `PATH`
- built `u-dales` executable available at the path expected by `run_test.sh`
- project module stack loaded where required by the local environment

After `initglobal` the test asserts the identities the rest of the code relies
on: `zh(kb) = 0`, `dzf(k) = zh(k+1) - zh(k)`, each centre `zf(k)` is the midpoint
of its faces (which makes the plain half-to-full-level averages exact on any
grid), `dzh(k) = zf(k) - zf(k-1)` equals the mean of the adjacent `dzf`, the
ghost-level conventions, the inverse and squared metric arrays, the kappa-scheme
copies and the subgrid length scale `delta`.

On a stretched grid it also requires a Poisson configuration that solves the
vertical direction by Gaussian elimination (`ipoiss = 0` with `BCzp = 1`, or
`ipoiss = 3`). The `stretched_badpois` mode combines the stretched profile with
`BCzp = 2` and must fail: `initglobal` refuses the combination before the test
runs, and the test checks the same rule.

Run it directly with:

```bash
cd tests/integration/zgrid_metrics
./run_test.sh
```
