# uDALES simulations

Curated, ready-to-use setups for real uDALES simulations — reference
configurations that should always be reachable from the repo, as opposed to
the minimal teaching cases in `examples/`. Categories:

- **Benchmarks** — frozen cases for measuring performance across machines
  and code versions.
- **Reference simulations** — production configurations of general value,
  e.g. the planned St Paul's simulation that will provide the data for the
  tutorials.

| Directory | Category | Purpose |
|-----------|----------|---------|
| `benchmark-standard-900` | benchmark | standard tier: 256×256×128 (8.4M cells), ~1.1k facets, 64 ranks (8×8) |
| `benchmark-large-901` | benchmark | large tier: 2048×2048×256 (1.07e9 cells), ~1.13M facets, ≥1024 ranks (32×32) |
| *(planned)* `stpauls-NNN` | reference | St Paul's production simulation; source of the tutorial datasets |

## Conventions

- Directory names are descriptive with a trailing three-digit case number:
  `<name>-NNN`. The tooling derives `iexpnr` from the LAST THREE characters
  of the directory name, so the numeric suffix is required and must match
  `iexpnr` in the namoptions file.
- Each case directory contains the namoptions file, the geometry source
  (commit small STLs directly; large geometry ships as a deterministic
  generator script plus its exact invocation), and any case-specific notes.
- Preprocessed IBM geometry (`solid_*`, `facet_sections_*`, …) IS committed
  when its total size permits (tens of MB of text), so the case runs without
  MATLAB preprocessing. Where it is too large for git (~GB, e.g. the large
  benchmark tier), it stays out: regenerate per machine or fetch an archived
  copy, and record where the archive lives in the case notes.
- To use a case, copy its directory into your experiments directory (the
  copy may keep the full name or just the number — only the trailing digits
  matter), add a `config.sh` (see `docs/udales-simulation-setup.md`), then
  run. Cases that fit ship WITH their preprocessed IBM geometry
  (`benchmark-standard-900` does), so no MATLAB preprocessing is needed;
  cases whose geometry is too large for git (`benchmark-large-901`, ~GB
  scale) require running the preprocessing once per machine
  (`tools/write_inputs.sh <case> c`) or fetching an archived copy (see the
  case notes).
- Treat committed cases as frozen: changes invalidate cross-version and
  cross-machine comparisons. New configurations get new case numbers and a
  row in the table above.

## Benchmarks

Both benchmark tiers are the same physics: neutral urban LES over the
Xie & Castro (2008) staggered-cube array, immersed boundaries with
Werner–Wengle/Uno wall functions, Vreman subgrid model, FFT-based Poisson
solver, fixed timestep (dt = 0.05 s), and **no output of any kind**
(`trestart` >> `runtime`; all dump/statistics switches off). Grid spacing is
identical in both tiers (dx = dy = 0.625 m, dz = 0.9375 m); the large tier
is the standard geometry tiled 8×8 (1280×1280×240 m domain) with the STL
triangulation subdivided twice to reach ~1e6 facets:

```bash
python .github/skills/udales-perf/scripts/make_large_case.py \
    simulations/benchmark-standard-900/xie_castro_2008.stl \
    $DA_EXPDIR/901/xie_castro_tiled.stl
```

The timing metric is the `TOTAL CPU time` printed at exit, which brackets
exactly the time-integration loop, so s/step = TOTAL CPU time / nsteps.
Calibrate `runtime` so one measurement takes 5–10 minutes on your
machine/rank count.

The measurement procedure (compiler-flag audit, MPI viability probes, the
noise-robust timing protocol, scaling ladders) lives in the `udales-perf`
skill: `.github/skills/udales-perf/SKILL.md`, with per-machine results
recorded in its `references/machines.md`. Issue #330 is a worked campaign
using these cases.

Rules for comparable numbers:

- Never enable output or adaptive timestepping in a timing run.
- Report: commit, compiler + flags (from `flags.make`), MPI, node/CPU type,
  ranks × nodes, nsteps, and s/step (min over clean repetitions + spread).

Sizing the large tier: memory for prognostic + work fields is ~250+ GB
total (1e9 cells × ~30 arrays × 8 B) — choose node count from memory first.
The decomposition must satisfy: nprocx divides itot, nprocy divides jtot
and ktot; nprocx = nprocy = 32 is the reference decomposition. Its
preprocessing is a heavyweight batch job (hours, hundreds of GB) — generate
once per machine and archive the outputs.
