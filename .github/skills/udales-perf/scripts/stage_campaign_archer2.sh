#!/usr/bin/env bash
# HANDOFF SCRIPT — not yet run. Written 2026-07-21 for a follow-up agent/session
# to execute the ARCHER2 staged #330 replication campaign. See
# .github/skills/udales-perf/references/machines.md ("ARCHER2 — 2026-07-21")
# for what has already been established (single-point release/debug timing,
# MPI viability, the debug-build FPE-crash finding at src/modsubgrid.f90:357)
# and what is still missing: the actual B0-B4 staged before/after numbers that
# mirror the #330 PR's own campaign (currently only run on ICL CX3).
#
# WHY THIS EXISTS
# ----------------
# A prior pass on this machine benchmarked only the tip of the `performance`
# branch (one data point). That does not reproduce the PR's methodology: the
# PR measured SIX specific commits in sequence (base, then five cumulative
# optimization batches labelled B0-B4) to show the incremental effect of each
# change. This script builds and benchmarks all six commits on ARCHER2 so the
# same table can be produced here, comparable to the existing CX3 entry.
#
# THE SIX STAGE COMMITS (confirmed present in local history; `git log
# --oneline -1 <sha>` on each returns a real commit — verified 2026-07-21):
#
#   tag          commit      subject
#   base(master) 6f5c4494    (pre-#330 tip; = fb90a0cb^)
#   B0  flags    fb90a0cb    perf: restrict FPE trapping to Debug, add arch targeting to Release (#330 Tier 0)
#   B1  looporder8a0bb7c2    perf: i-innermost loop order in advection, tderive, forces, Smagorinsky (#330 Tier 1)
#   B2  tstep    ae25ba0e    perf: fuse tendency zeroing and substep-3 copies into integration nests (#330 Tier 1)
#   B3  poisson  9ac67d41    perf: persistent Poisson work arrays, drop dead solver copies (#330 Tier 1)
#   B4  ibm      ec47eaeb    perf: accumulate IBM stress/flux deltas directly, drop snapshot/diff sweeps (#330 Tier 1)
#
# (The PR text lists B4's pre-rebase SHA as d92ab46d*; on this repo's history
# that content landed as ec47eaeb — verified: it is the tip of the perf
# commits on `origin/performance`.) CX3 reference numbers for the same six
# points (s/step, standard-tier case): base 0.4289, B0 0.4219 (-1.6%),
# B1 0.4218, B2 0.3986 (-7.1%), B3 0.3773 (-12.0%), B4 0.3269 (-23.8%,
# 1.312x). ARCHER2's job is to reproduce this table, NOT to match the exact
# percentages (different CPU/compiler).
#
# APPROACH: git worktree per commit (NOT `git checkout` in the main tree —
# that would require serializing all 6 builds through one working copy and
# risks leaving it on a detached HEAD). Each worktree gets its own
# build/release via the existing tools/build_executable.sh, then is
# benchmarked with the existing scripts/bench_archer2.sh using its
# BENCH_BIN/BENCH_COMMIT override (that override exists exactly for this
# use case — see bench_archer2.sh's header).
#
# WHAT THIS SCRIPT DELIBERATELY DOES NOT DO
# ------------------------------------------
# - Does not touch the main /work/.../u-dales checkout's branch/HEAD.
# - Does not attempt Debug-build timing for any stage: the tip of this branch
#   (B4/HEAD) was already confirmed to crash with a genuine SIGFPE before
#   step 1 for the Debug build of this exact benchmark case (uninitialized
#   halo read at src/modsubgrid.f90:357, independent of the -ei flag — see
#   machines.md gotchas). That bug likely predates B0-B4 (it is in the SGS
#   closure, untouched by any of these commits per their subjects above), so
#   assume all six Debug builds crash identically; do not waste queue time
#   confirming this per-stage unless investigating that bug specifically.
# - Does not run scaling (Phase 5) — single node, 64 ranks, same as the
#   existing ARCHER2 entry.
# - Does not append the final table to machines.md — that is the last step
#   below, left for the agent running this script (numbers are not known
#   until the jobs complete).
#
# HOW TO RUN
# ----------
#   cd /work/n02/n02/dip_assure/udupdate/u-dales
#   git fetch origin   # per AGENTS.md preflight
#   bash .github/skills/udales-perf/scripts/stage_campaign_archer2.sh
#
# Then poll `squeue -u $USER` / inspect
# /work/n02/n02/dip_assure/udupdate/experiments/900/results.csv (columns:
# tag,rep,jobid,commit,node,cpu_model,nsteps,total_cpu_time_s,s_per_step) for
# rows tagged base/B0/B1/B2/B3/B4. Compute min-over-reps per tag, then delta
# vs base, then append a table to machines.md's ARCHER2 section (do not
# overwrite the existing entry — extend it, matching the CX3 table format).
# Finally: `git worktree remove` (or `git worktree remove --force` if a
# build dir is still present) each entry printed by `git worktree list` under
# WORKTREE_ROOT, to reclaim disk.

set -euo pipefail

REPO=/work/n02/n02/dip_assure/udupdate/u-dales
WORKTREE_ROOT=/work/n02/n02/dip_assure/udupdate/perf-worktrees
NREPS=${NREPS:-3}

cd "$REPO"
if [ -n "$(git status --porcelain)" ]; then
    echo "ERROR: $REPO has uncommitted changes. This script only reads history" >&2
    echo "(via worktrees) but a dirty main tree usually means unfinished work -" >&2
    echo "resolve that first so nothing here is confused with it." >&2
    exit 1
fi

mkdir -p "$WORKTREE_ROOT"

# tag:commit pairs, in the order the PR measured them.
STAGES=(
    "base:6f5c4494"
    "B0:fb90a0cb"
    "B1:8a0bb7c2"
    "B2:ae25ba0e"
    "B3:9ac67d41"
    "B4:ec47eaeb"
)

LASTJOB=""
for entry in "${STAGES[@]}"; do
    tag=${entry%%:*}
    commit=${entry#*:}
    wtdir="$WORKTREE_ROOT/$tag"

    echo "=== stage $tag ($commit) ==="

    if [ -d "$wtdir" ]; then
        echo "worktree $wtdir already exists, reusing (remove it first for a clean rebuild)"
    else
        git worktree add --detach "$wtdir" "$commit"
    fi

    (
        cd "$wtdir"
        # View3D/f2py preprocessing extras are NOT needed to build/run the
        # solver binary benchmarked here, so the submodule is left alone.
        ./tools/build_executable.sh archer release
    )

    bin="$wtdir/build/release/u-dales"
    [ -x "$bin" ] || { echo "ERROR: build for $tag did not produce $bin" >&2; exit 1; }

    full_commit=$(cd "$wtdir" && git rev-parse HEAD)

    if [ -n "$LASTJOB" ]; then
        export BENCH_DEPEND=$LASTJOB
    else
        unset BENCH_DEPEND || true
    fi
    export BENCH_BIN="$bin"
    export BENCH_COMMIT="$full_commit"

    out=$(bash "$REPO/.github/skills/udales-perf/scripts/bench_archer2.sh" "$tag" "$NREPS" release)
    echo "$out"
    LASTJOB=$(echo "$out" | awk '/^LASTJOB/ {print $2}')
done

echo
echo "All stages submitted (chained). Watch with:"
echo "  squeue -u \$USER"
echo "Results land in:"
echo "  $REPO/experiments/900/results.csv"
echo
echo "When done, reclaim disk with, for each stage:"
echo "  git -C $REPO worktree remove $WORKTREE_ROOT/<tag>"
