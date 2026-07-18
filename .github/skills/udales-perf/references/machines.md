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
  shm anyway). This also reopens whole-node/large-class requests → true
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
