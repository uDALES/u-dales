# Machine Performance Profiles

Per-machine performance-commissioning records (udales-perf skill, Phase 6).
Append dated entries; do not overwrite. Cluster build/run basics live in
`../../udales-exec/references/clusters.md` — this file holds what matters for
TIMING: node classes and placement, MPI viability per class, benchmark
numbers and their exact configuration, and measurement gotchas.

## Template

```
### <machine> — <date>

Toolchain/flags: <modules; effective Release flags from flags.make>
Benchmark case: <grid, decomposition, physics, nsteps, metric>
Node classes probed: <class -> MPI OK/BROKEN>
Placement protocol: <exclusive | paired round-robin chain + node matching>
Results: <s/step release/debug; scaling table>
Gotchas: <...>
```

## ICL CX3 — 2026-07-18

Toolchain/flags: `tools/prod` + intel/2021a stack (see clusters.md); Release
`-r8 -heap-arrays 10 -O3 -march=core-avx2` (after #330 Tier 0; `-fpe0` is
Debug-only). `-march=core-avx2` chosen over `-xHost` because compute nodes are
AMD Rome (EPYC 7742) and Intel `-x` flags embed a GenuineIntel check.

Benchmark case: neutral urban IBM LES (Xie–Castro staggered cubes),
256×256×128, nprocx=8 × nprocy=8 = 64 ranks, Vreman SGS, FFT2D Poisson,
fixed dt=0.05 s, runtime=50 s → 1001 steps ≈ 7 min; all output off,
trestart >> runtime. Metric: `TOTAL CPU time` / nsteps (brackets the time
loop only). Debug variant: 101 steps.

Node classes probed (jobs route by size: ≤16 cores small, ≤64 medium,
≤128 large, ≤256 capability; suffix 24/72 = walltime band):
- 64-core medium class (racks 2–6 era, EPYC 7742 128-core nodes): MPI OK.
- **CORRECTED DIAGNOSIS (2026-07-18):** hydra bootstrap failures are NOT an
  OS/toolchain incompatibility. Affected nodes (observed: cx3-12-24,
  cx3-13-13, cx3-14-14, cx3-15-17; healthy: cx3-15-0, racks 2–6) firewall
  TCP connections TO THEMSELVES — the proxy times out connecting
  node→same-node ("check for firewalls!" in the .ER file). Fix:
  `export I_MPI_HYDRA_IFACE=lo` before mpiexec (single-node runs; ranks use
  shm anyway). PROBE-CONFIRMED on cx3-14-14 (2026-07-19): plain mpiexec
  fails, with IFACE=lo it works; 36/36 production chain jobs succeeded with
  the fix. This also reopens whole-node/large-class requests → true
  exclusivity may be attainable; retest. Report the node list to RCS.
  (Side note: newer toolchains live under `module load tools/dev` but
  intel/2023a did not co-load cleanly as of 2026-07.)
- `place=excl` in the medium class: never scheduled within 1.5 h. Whole-node
  requests (128 cores) DO schedule; with the I_MPI_HYDRA_IFACE=lo fix they
  are the route to true exclusivity — retest and update here.

Placement protocol: paired round-robin dependency chain in the medium class
(all binaries alternate each round; 6 release rounds + 2 debug rounds), then
within-round, same-node comparison, min over clean rounds. Observed noise
without pairing: ~8% between serial same-node reps (CPU-slot placement /
external co-tenants); pairing reduces round-to-round spread to <1%.
64-core jobs pair up two-per-node — never compare timings from jobs that
overlapped on a node.

Tier: standard (256×256×128) measured. Large tier: geometry generated
(1,128,448 facets, 1280×1280 m; `~/udales/experiments/901`, 2048×2048×256)
and preprocessing launched 2026-07-18 as a 16-core/800 GB/24 h batch job —
record its runtime/memory here when it finishes. RUNNING the case needs
multi-node MPI, i.e. the large-class/toolchain problem above solved first.

Results (issue #330 campaign, pre-rebase reference, 2026-07-18): base
0.414–0.429 s/step depending on node; #330 Tier-0/1 batches cumulatively
−23% (flags −1%, loop order −1.4%, tstep fusion −6.4%, Poisson −11.3%,
IBM −23%). Release/debug ratio ≈ 12.6×. See the #330 PR for final rebased
numbers.

Gotchas:
- Scheduler feeds a whole queued campaign to one node class; probe first.
- qstat's job table drops finished jobs silently — watchers must distinguish
  "all jobs finished" from "qstat shows nothing"; verify expected output
  files/rows exist before declaring completion.
- PBS `-o/-e` default to the submission cwd; point them at a log dir
  (a stray spool file in the repo makes the working tree dirty).
- Benchmark jobs must snapshot binary + inputs at submission: a rebuild
  while jobs are queued otherwise races them (observed: job exec'd a
  half-deleted binary).

## ARCHER2 — 2026-07-21

Toolchain/flags: PrgEnv-cray/8.4.0, cce/16.0.1, craype/2.7.23,
craype-x86-rome, cray-mpich/8.1.27, cray-libsci/23.09.1.1,
cray-fftw/3.3.10.5, cray-netcdf/4.9.0.7 (crayclang/14.0). Built via
`tools/build_executable.sh archer <release|debug>` (`FC=ftn`). Effective
Release flags (from `build/release/CMakeFiles/u-dales.dir/flags.make`):
`-s real64 -N 1023 -O3 -em -J.` — no FPE trapping, matches the #330
Cray branch already on `performance` (no source changes needed here).
Architecture targeting comes from the `craype-x86-rome` module, not an
explicit flag: `ps` on a running job shows `ftn_driver.exe -hcpu=x86-rome
-hnetwork=ofi -hdynamic -D__CRAY_X86_ROME ...`. Effective Debug flags:
`-s real64 -N 1023 -G0 -O0 -ei -R bcdps -m 0 -K trap=divz,inv,ovf -em -J.`
(FPE trapping + `-ei` uninitialized-value poisoning, Debug-only — correct).

Benchmark case: same frozen case as CX3 (`simulations/benchmark-standard-900`,
neutral urban IBM LES, Xie–Castro staggered cubes, Vreman SGS, FFT2D
Poisson, 256×256×128, nprocx=8×nprocy=8=64 ranks, dt=0.05 s, runtime=50 s →
1001 steps). Metric: `TOTAL CPU time` / nsteps. Debug variant: runtime=5 s
(→101 steps planned; see Gotchas — it never reached step 1).

Node classes probed: `standard` partition, `short`/`standard` QoS, 1 node
(128-core AMD EPYC 7742 "Rome" ×2 = 64 used, `--tasks-per-node=64
--hint=nomultithread`). MPI probe (`probe_mpi_archer2.sh`, job 14447996 on
nid006080): plain `srun` and `srun --mpi=pmi2` both OK — no bootstrap issues
(unlike CX3's hydra/firewall saga). ARCHER2 bills and allocates whole nodes
per Slurm job, so single-node jobs here are inherently exclusive; no
paired-comparison protocol was needed.

Placement protocol: exclusive whole-node allocation (Slurm default on this
system) — reps simply chained serially via `--dependency=afterany` for
convenience, not for co-tenancy avoidance. Harness:
`.github/skills/udales-perf/scripts/bench_archer2.sh <tag> [nreps]
[release|debug]`, a Slurm port of `bench.sh` (git-clean + binary-staleness
guards, per-rep snapshot of binary+inputs, CSV harvest, chained `sbatch`).

Results (release, commit `edd3c0d169`, 3 reps, tag `perfbranch`):

| rep | node       | nsteps | total CPU time (s) | s/step  |
|-----|------------|--------|---------------------|---------|
| 1   | nid001126  | 1001   | 360.46               | 0.36010 |
| 2   | nid001126  | 1001   | 359.03               | 0.35868 |
| 3   | nid004837  | 1001   | 355.02               | 0.35466 |

min 0.35466 s/step, spread ~1.5% across reps/nodes — consistent with
exclusive whole-node allocation (no pairing needed to see this level of
repeatability). This is the current `performance`-branch baseline on
ARCHER2; no #330-style optimization campaign has been run here yet, so
there is no before/after delta to report (only a single toolchain/config
point).

Debug/release ratio: **not measured** — the Debug build crashes with
SIGFPE before completing even step 1 (before its first stdout line),
100% reproducible across two different nodes/jobs. Root-caused via core
dump + `gdb -batch -ex bt`:
```
#0 closure () at src/modsubgrid.f90:357   -- ekh(:,:,:) = ekm(:,:,:)*prandtli
#1 subgrid () at src/modsubgrid.f90:139
#2 udales () at src/program.f90:147       -- first `call subgrid` in the main loop
```
In the `lvreman` branch of `closure()`, `ekm`/`ekh` are only assigned over
the interior loop range (`ib:ie, jb:je, kb:ke`); the very next line does a
*whole-array* multiply `ekh(:,:,:) = ekm(:,:,:)*prandtli` that also touches
the halo, which is still uninitialized at this point (`call closurebc`,
which sets the K-closure boundary/halo values, only runs at the end of
`closure()`, after this line). Confirmed independent of `-ei`: rebuilding
Debug with `-ei` dropped (isolated `build/debug_noei`, CMakeLists.txt edit
reverted immediately after configure) still crashes identically, so this is
a real `-K trap=inv,ovf` catching genuine garbage in the halo, not an
`-ei`-only false positive. Release silently computes over the same garbage
and it is presumably overwritten by the following `closurebc` call, so the
result is likely harmless in practice — but it is a real uninitialized-read
and should be fixed (either restrict the multiply to `ib:ie,jb:je,kb:ke`, or
move it after `closurebc`, or zero-init `ekm`/`ekh` on allocation). Worth a
follow-up issue; not fixed as part of this commissioning pass. Not hit on
CX3/Intel for the same case/flags in the prior campaign — likely just a
difference in what garbage bit pattern happened to be sitting in the
freshly-allocated halo memory on each machine, not a toolchain-specific bug.

Gotchas:
- Slurm `#SBATCH` `-o/-e` default to the submission cwd exactly like PBS —
  a probe script without them dropped `slurm-<jobid>.out` into the repo
  root (dirty tree); always set `-o/-e` explicitly.
- No git identity was configured for this user on this machine; the
  provenance guard requires a clean tree + commit sha, so a repo-local (not
  global, not pushed) `user.name`/`user.email` was set to allow committing
  the new scripts.
- Harness bug (fixed in `bench_archer2.sh`): the per-rep cleanup was
  `cd; rm -rf $rundir` — the bare `cd` (to `$HOME`) can silently fail under
  `set -e` in the batch environment, aborting the script before the `rm`
  runs. Fixed to just `rm -rf $rundir` (no `cd` needed, cwd is not `$rundir`
  by then in a way that matters). Observed: 3 release run dirs (~38 MB
  each incl. `facet_sections_*`/`fluid_boundary_*` inputs) were left behind
  in `outputs/bench900/` despite successful CSV harvest.
- Core dumps from `srun`-launched FPE crashes land as a plain `core` file
  in the run's cwd (`ulimit -c` is `unlimited` by default here); `gdb -batch
  -ex bt -ex "info registers" ./binary core` gives a usable backtrace
  straight away — useful for root-causing Debug-only trap crashes without
  needing to reproduce interactively.
- ARCHER2 login node compiles are slow for at least one file: `-O3` on
  `src/modstatsdump.f90` (2172 lines) took ~565 s (single-core `optcg`
  backend at ~100% CPU) — budget for this, it is not a hang.

