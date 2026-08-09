# CUDA Extended-Halo Initializer

This is a GPU integration check. It is not part of the CPU-only GitHub Actions
gate because it requires an NVHPC CUDA build and a CUDA-capable device. It is
owned and automatically enabled by the `gpu-smoke` and `gpu-mpi` selections.

The test implementation lives in the dedicated Fortran test module
`src/tests_cuda.f90`. The Python GPU runner opts into it by setting
`UDALES_RUN_CUDA_SELFTEST=1`; ordinary Debug GPU simulations do not execute
test-only code. The test fills an array having the scalar/Kappa halo widths
(`ihc`, `jhc`, `khc`) with a non-zero sentinel, launches the production
`initfield` kernel, copies the complete allocation back to the host, and
terminates the run if any element was not reset. A successful rank writes:

```text
CUDA extended-halo initfield self-test passed. rank=0
```

For a standalone Debug GPU run, request the test and verify its log with:

```bash
export UDALES_RUN_CUDA_SELFTEST=1
mpiexec -n 1 /path/to/u-dales namoptions.103 > output.log 2>&1
python tests/integration/cuda_extended_halo_initializer/check_debug_log.py \
  --expected-ranks 1 output.log
```

The standalone log checker remains useful for existing manual simulations. The
automated parity runner sets the environment flag itself and requires exactly
one rank-qualified success marker per MPI rank when `--require-debug-selftest`
is active.
