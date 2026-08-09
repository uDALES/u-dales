# CUDA Extended-Halo Initializer

This is a local GPU integration check. It is not part of the CPU-only GitHub
Actions gate because it requires an NVHPC CUDA build and a CUDA-capable device.

Every Debug GPU executable runs an initializer self-test from `initCUDA`. The
test fills an array having the scalar/Kappa halo widths (`ihc`, `jhc`, `khc`)
with a non-zero sentinel, launches `initfield`, copies the complete allocation
back to the host, and terminates the run if any element was not reset. A
successful run writes:

```text
CUDA extended-halo initfield self-test passed.
```

After running any Debug GPU case, verify that the self-test executed with:

```bash
python tests/integration/cuda_extended_halo_initializer/check_debug_log.py /path/to/output.log
```

This check is intended for local/manual GPU validation.
