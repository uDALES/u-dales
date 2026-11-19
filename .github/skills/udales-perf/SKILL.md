---
name: udales-perf
description: Commission uDALES on a new or changed machine from a performance standpoint - verify compiler settings, probe MPI viability per node class, run a standard benchmark with a noise-robust timing protocol, measure scaling, and record a comparable performance profile.
---

# Performance Commissioning Skill

This skill establishes trustworthy performance numbers for uDALES on a machine:
single-node speed, parallel scaling, and the compiler/placement configuration
that produced them. It exists because timing measurements on shared clusters
are dominated by avoidable systematic errors (wrong node class, broken MPI
launchers, co-tenant contention, CPU-vendor flag traps) that waste far more
time than the measurements themselves.

It is distinct from correctness/regression testing (see `tests/`) and from
day-to-day building/running (see the `udales-exec` skill).

## When to use

- "How fast is uDALES on this machine?" / "Is this machine worth using?"
- After a compiler/toolchain/OS upgrade on a known machine.
- Before and after performance-related code changes, when the effect must be
  quantified (see issue #330 for a worked campaign).
- When choosing core counts / node counts for a production run.

## Phase 0 — Identify the environment

Use the `udales-detect` skill and `udales-exec`'s `references/clusters.md` to
identify the machine and its known build stack. If this machine already has an
entry in `references/machines.md` (this skill), start from that profile.

## Phase 1 — Compiler-settings audit

1. Build release and debug with the repo wrapper (`tools/build_executable.sh`).
2. Verify the flags that were ACTUALLY used — read
   `build/<type>/CMakeFiles/u-dales.dir/flags.make`, not CMakeCache (the cache
   shows defaults, not the effective directory-scoped flags).
3. Check, for Release:
   - NO floating-point-exception trapping (`-ffpe-trap`, `-fpe0`, `-K trap`) —
     it belongs in Debug only; in Release it blocks vectorization (#330).
   - Architecture targeting matches the COMPUTE nodes, not the login node:
     - GNU: `-march=native` only if build and compute nodes share an
       architecture; otherwise an explicit `-march=<arch>`.
     - Intel on AMD hardware: use `-march=core-avx2` (or similar). NEVER
       `-xHost`/`-xCORE-AVX2`: the `-x` variants embed a GenuineIntel CPU
       check and refuse to run, or underdispatch, on AMD.
4. Record compiler + MPI + library versions (module list) for the report.

## Phase 2 — MPI viability probe (do this BEFORE any benchmark jobs)

Clusters often have several node classes (different racks, OS images, queue
routing by job size). An MPI stack that works on one class can fail on
another — e.g. Intel MPI 2021.x hydra cannot bootstrap on newer OS images:
`error setting up the bootstrap proxies` even though ssh works.

For EVERY node class the campaign will touch (each queue/size class you plan
to submit to), submit `scripts/probe_mpi.sh` as a batch job with that class's
resource request and confirm `mpiexec hostname` works there. A 2-minute probe
prevents entire job campaigns dying at `mpiexec` (observed: 15+ jobs lost to
an unprobed node class).

## Phase 3 — Benchmark cases and calibration

Two tiers of the SAME case (urban IBM geometry — Xie–Castro staggered cubes —
the production subgrid model and Poisson solver), differing only in size.
**The case definitions live in the repo under `simulations/`**
(`benchmark-standard-900`, `benchmark-large-901`; see `simulations/README.md`
for usage and the rules for comparable numbers):

- **Standard tier** (every machine): 256×256×128 ≈ 8.4M cells, ~1.1k facets,
  64 ranks on one node. Fast to preprocess and run; the cross-machine
  single-node reference.
- **Large tier** (machines with O(10+) nodes): O(1024³) ≈ 1e9 cells and
  ~1e6 facets — the same cube geometry tiled to a proportionally larger
  domain at production resolution, with the STL triangulation refined to
  reach the facet-count target. Reference implementation:
  `scripts/make_large_case.py <unit.stl> <out.stl>` (defaults: 8×8 tiles at
  160 m pitch, 2 subdivision levels → 1,128,448 facets over 1280×1280 m),
  paired with itot=jtot=2048, ktot=256 (same dx=dy=0.625 m, dz=0.9375 m as
  the standard tier; zsize=240 m), nprocx=nprocy=32, and otherwise the
  standard tier's namoptions. This is the tier that actually stresses
  multi-node behaviour: Poisson transposes, halo exchange at scale, the
  O(nfcts) IBM machinery, memory per node, and preprocessing at production
  scale. Practical notes:
  - fields alone need ~250+ GB (1e9 cells × ~30 arrays × 8 B) — size the
    minimum node count from memory before anything else
  - preprocessing (facet sections, solid/fluid masks) is a serious job at
    this size: run it as a compute-node batch job and archive its outputs —
    generate once, reuse for every benchmark on that machine
  - calibrate `runtime` the same way (one measurement 5-10 min); expect to
    calibrate at the intended production node count, not on one node

Common rules for both tiers:

- every dump/statistics switch off; `trestart` >> `runtime` (no restart files)
- FIXED timestep (`ladaptive = .false.`) so every run does identical work
- `runtime` sized by a short calibration run so one measurement takes 5-10 min
- metric: the `TOTAL CPU time` printed at exit — it brackets exactly the time
  loop (excludes startup), so s/step = TOTAL CPU time / nsteps
- the case in `references/machines.md` entries names grid, decomposition and
  step count so numbers are comparable across machines

## Phase 4 — Timing protocol (the important part)

Shared-cluster timing noise is one-sided (contention only slows) and can be
~8% even for serial reps on one node (CPU-slot placement, external
co-tenants). Rules, in order of preference:

1. **Exclusive nodes** if the scheduler allows them AND MPI works there
   (Phase 2!). Whole-node requests also give defined NUMA behaviour.
2. Otherwise: **paired round-robin chain** — submit all binaries/configs
   under comparison as ONE serial dependency chain (job N+1 depends on job N),
   alternating configurations each round, several rounds. Each round samples
   all configurations under near-identical conditions.
3. Compare only WITHIN rounds, and only cells that ran on the same node
   (record `hostname` and CPU model per run); discard mixed-node rounds.
   Report min-over-clean-rounds plus the spread.
4. Every job snapshots its binary and inputs at submit time and records the
   git commit — a timing row must map to an exact commit (no dirty trees).
5. Also time the debug build (fewer steps; s/step is the metric either way):
   it catches order-of-magnitude debug regressions and documents the
   release/debug ratio (~12x on CX3).

`scripts/bench.sh`, `scripts/verify.sh` and `scripts/compare_restarts.py` are
tested templates (PBS; adapt the directives for other schedulers). bench.sh
implements the provenance guards, snapshotting, self-harvesting to CSV and
dependency chaining (`BENCH_DEPEND`); verify.sh + compare_restarts.py give
restart-file-based verification that a code change did not change physics
(expect agreement at ~1e-15 rel Linf, not bitwise — default Intel `-fp-model
fast` reassociates fused loops).

## Phase 5 — Scaling

Strong-scaling ladders at fixed grid, same protocol as Phase 4 (chained,
node-recorded):

- Standard tier: quarter node → half node → full node → 2, 4 nodes.
- Large tier (large machines): from the minimum node count that fits memory,
  doubling upward as far as the machine sensibly allows.

Report s/step and parallel efficiency vs the smallest configuration of each
ladder. Note where the Poisson transposes start to dominate — on the large
tier this is the number that determines the machine's useful production size.

## Phase 6 — Report and self-update rule

Append to `references/machines.md`: machine, date, toolchain/flags, node
class + placement rules that worked, benchmark s/step (with grid/ranks/steps),
scaling table, and every gotcha discovered. Do not overwrite existing entries;
append dated updates. This file is the cross-machine comparison record.
